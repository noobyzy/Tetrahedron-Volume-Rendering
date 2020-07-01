
#include <iostream>
#include <algorithm>
#include <tetra.h>



#include"renderer.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define WIDTH 1024
#define HEIGHT 1024

#define MIN4(a,b,c,d) (((a)<(b)?(a):(b))<((c)<(d)?(c):(d))?((a)<(b)?(a):(b)):((c)<(d)?(c):(d)))
#define MAX4(a,b,c,d) (((a)>(b)?(a):(b))>((c)>(d)?(c):(d))?((a)>(b)?(a):(b)):((c)>(d)?(c):(d)))
#define SIGN(x) (((x) > 0)?1:-1)

void renderSingle(){
	/*
	 * 1. Volume Setting
	 */
	
	Volume vol("data/test2.bin");
	/*
	 * 2. Camera Setting
	 */
	Eigen::Vector3f cameraPosition= vol.bbox.getCenter()-2.5*Eigen::Vector3f(0,0,vol.size_physics.z());
	Eigen::Vector3f cameraLookAt= vol.bbox.getCenter();
	Eigen::Vector3f cameraUp(0, 1, 0);
	float verticalFov = 45;
	Eigen::Vector2i filmRes(1024, 1024);

	Camera camera(cameraPosition, cameraLookAt, cameraUp, verticalFov, filmRes);
	PointLight light(cameraPosition+Eigen::Vector3f(0,5,0),Eigen::Vector3f(1,1,1));



	/*
	 * 3. Renderer Setting
	 */

	myClassifier classifer;
	TrilinearInterpolator interpolator;

	Renderer renderer;
	renderer.setCamera(&camera);
	renderer.setVolume(&vol);
	renderer.setClassifier((Classifier *)&classifer);
	renderer.addLight((Light *)&light);
	renderer.setInterpolator((Interpolator *)&interpolator);
	
	/*
	 * 3. Render
	 */
	renderer.renderFrontToBack("./foward_pre.png");
}

void multipleRender(){
/*
	 * 1. Volume Setting
	 */
	
	Volume vol("data/test2.bin");
	/*
	 * 2. Camera Setting
	 */
	Eigen::Vector3f cameraPosition= vol.bbox.getCenter()-2.5*Eigen::Vector3f(0,0,vol.size_physics.z());
	Eigen::Vector3f cameraLookAt= vol.bbox.getCenter();
	Eigen::Vector3f cameraUp(0, 1, 0);
	float verticalFov = 45;
	Eigen::Vector2i filmRes(1024, 1024);

	Camera camera(cameraPosition, cameraLookAt, cameraUp, verticalFov, filmRes);
	PointLight light(cameraPosition+Eigen::Vector3f(0,5,0),Eigen::Vector3f(1,1,1));



	/*
	 * 3. Renderer Setting
	 */

	myClassifier classifer;
	TrilinearInterpolator interpolator;

	Renderer renderer;
	renderer.setCamera(&camera);
	renderer.setVolume(&vol);
	renderer.setClassifier((Classifier *)&classifer);
	renderer.addLight((Light *)&light);
	renderer.setInterpolator((Interpolator *)&interpolator);
	for(int frame = 50;frame<200;frame++){
		  char filein[256];
    	sprintf(filein,"data/density/density_frame_out%05d.bin",frame);
		vol=Volume (filein);
		renderer.setVolume(&vol);
		char fileout[256];
    	sprintf(fileout,"result/frame%05d.png",frame);
		renderer.renderFrontToBack(fileout);
	}
}

std::vector<Eigen::Vector2f> ComputeScreenSpaceProjections(std::vector<MyVertex> Vertices, Camera camera){
	std::vector<Eigen::Vector2f> SSC;
	for(int i=0; i<Vertices.size(); ++i){
		Ray ray(camera.m_Pos, Vertices[i].coordinate - camera.m_Pos);
		Eigen::Vector3f p = (1.0f / ray.m_Dir.dot(camera.m_Forward))*ray.m_Dir - camera.m_Forward;
		SSC.push_back(Eigen::Vector2f(p.dot(camera.m_Right), p.dot(camera.m_Up)));
	}
	return SSC;
}

/* test the side of the point */
int line_side(Eigen::Vector2f va, Eigen::Vector2f vb, Eigen::Vector2f p){
	float x0 = va.x();
	float y0 = va.y();
	float x1 = vb.x();
	float y1 = vb.y();

	float result = p.y() - ((y0 - y1)/(x0 - x1)) * p.x() - y1 + ((y0 - y1)/(x0 - x1)) * x1;
	if(result > 0){
		return 1;
	} else if(result == 0){
		return 0;
	} else {
		return -1;
	}
}

/* to justify whether a point is outside the triangle */
bool is_inside_triangle(Eigen::Vector2f va, Eigen::Vector2f vb, Eigen::Vector2f vc, Eigen::Vector2f p){
	bool line_c_test, line_b_test, line_a_test;
	line_c_test = (line_side(va, vb, p) == line_side(va, vb, vc));
	line_b_test = (line_side(va, vc, p) == line_side(va, vc, vb));
	line_a_test = (line_side(vc, vb, p) == line_side(vc, vb, va));
	return (line_c_test && line_b_test && line_a_test);
}

float cross_product(Eigen::Vector2f v1, Eigen::Vector2f v2){
	return (v1.x()*v2.y() - v2.x()*v1.y());
}

/* get the projection of each tetrahedron, get the piexl projected by the tetrahedron, build a intersection list for each pixel */
void ExtractIntersectionRecords(std::vector<Tetrahedron> tetra_list, std::vector<MyVertex> vertex_list, std::vector<Eigen::Vector2f> SSC, std::vector<std::vector<std::vector<int>>>& PerPixelIntersectionList){
	// iterate each tetrahedron
	for(int i = 0; i < tetra_list.size(); i++){
		
		/* the four projected points on screen */
		Eigen::Vector2f v1_proj = SSC[tetra_list[i].v1_idx];
		Eigen::Vector2f v2_proj = SSC[tetra_list[i].v2_idx];
		Eigen::Vector2f v3_proj = SSC[tetra_list[i].v3_idx];
		Eigen::Vector2f v4_proj = SSC[tetra_list[i].v4_idx];

		/* lower bound and upper bound of x and y */
		int xlb = std::max((int)std::floor(MIN4(v1_proj.x(), v2_proj.x(), v3_proj.x(), v4_proj.x())), 0);
		int ylb = std::max((int)std::floor(MIN4(v1_proj.y(), v2_proj.y(), v3_proj.y(), v4_proj.y())), 0);

		int xub = std::min((int)std::ceil(MAX4(v1_proj.x(), v2_proj.x(), v3_proj.x(), v4_proj.x())), 1023);
		int yub = std::min((int)std::ceil(MAX4(v1_proj.y(), v2_proj.y(), v3_proj.y(), v4_proj.y())), 1023);


		if(is_inside_triangle(v1_proj, v2_proj, v3_proj, v4_proj) || is_inside_triangle(v1_proj, v2_proj, v4_proj, v3_proj) || is_inside_triangle(v1_proj, v3_proj, v4_proj, v2_proj) || is_inside_triangle(v2_proj, v3_proj, v4_proj, v1_proj)){
			/* at least one vertex is inside the triangle consist of other 3 vertices*/
			// the projected shape is a triangle

			// determine the 3 vertices of triangle
			Eigen::Vector2f tri_v1, tri_v2, tri_v3;
			if(is_inside_triangle(v1_proj, v2_proj, v3_proj, v4_proj) == true){
				tri_v1 = v1_proj;
				tri_v2 = v2_proj;
				tri_v3 = v3_proj;
			} else if(is_inside_triangle(v1_proj, v2_proj, v4_proj, v3_proj) == true){
				tri_v1 = v1_proj;
				tri_v2 = v2_proj;
				tri_v3 = v4_proj;
			} else if(is_inside_triangle(v1_proj, v3_proj, v4_proj, v2_proj) == true){
				tri_v1 = v1_proj;
				tri_v2 = v3_proj;
				tri_v3 = v4_proj;
			} else {
				tri_v1 = v2_proj;
				tri_v2 = v3_proj;
				tri_v3 = v4_proj;
			}

			// iterate through the square (xlb, xub, ylb, yub), check whether the pixel inside the triangle
			for(int pixel_i = xlb; pixel_i <= xub; pixel_i++){
				for(int pixel_j = ylb; pixel_j <= yub; pixel_j++){
					Eigen::Vector2f temp_pixel = Eigen::Vector2f((float)pixel_i, (float)pixel_j);
					if(is_inside_triangle(tri_v1, tri_v2, tri_v3, temp_pixel) == true){
						PerPixelIntersectionList[pixel_i][pixel_j].push_back(i);
					}
				}
			}


		} else {
			// the projected shape is a quardrilateral

			// determine the 4 points CounterClockWise
			float Min_x = MIN4(v1_proj.x(), v2_proj.x(), v3_proj.x(), v4_proj.x());
			Eigen::Vector2f leftmost_point;
			Eigen::Vector2f remain_1, remain_2, remain_3;
			if(Min_x == v1_proj.x()){
				leftmost_point = v1_proj;
				remain_1 = v2_proj;
				remain_2 = v3_proj;
				remain_3 = v4_proj;
			}else if(Min_x == v2_proj.x()){
				leftmost_point = v2_proj;
				remain_1 = v1_proj;
				remain_2 = v3_proj;
				remain_3 = v4_proj;
			}else if(Min_x == v3_proj.x()){
				leftmost_point = v3_proj;
				remain_1 = v1_proj;
				remain_2 = v2_proj;
				remain_3 = v4_proj;
			}else{
				leftmost_point = v4_proj;
				remain_1 = v1_proj;
				remain_2 = v2_proj;
				remain_3 = v3_proj;
			}

			//find the rightmost point
			Eigen::Vector2f rightmost_point, diff_1, diff_2;
			if(line_side(leftmost_point, remain_1, remain_2) != line_side(leftmost_point, remain_1, remain_3)){
				rightmost_point = remain_1;
				diff_1 = remain_2;
				diff_2 = remain_3;
			} else if(line_side(leftmost_point, remain_2, remain_1) != line_side(leftmost_point, remain_2, remain_3)){
				rightmost_point = remain_2;
				diff_1 = remain_1;
				diff_2 = remain_3;
			} else{
				rightmost_point = remain_3;
				diff_1 = remain_1;
				diff_2 = remain_2;
			}

			Eigen::Vector2f P, Q, l1, r1;
			P = rightmost_point - leftmost_point;
			Q = diff_1 - leftmost_point;
			//P is on the clockwise direction of Q
			if((P.x()*Q.y() - Q.x()*P.y()) > 0){
				r1 = diff_2;
				l1 = diff_1;
			} else {
				r1 = diff_1;
				l1 = diff_2;
			}
			
			//now that the counterclockwise order of 4 vertices is leftmost, l1, rightmost, r1
			Eigen::Vector2f A, B, C, D;
			A = leftmost_point;
			B = l1;
			C = rightmost_point;
			D = r1;

			Eigen::Vector2f AB = B - A;
			Eigen::Vector2f BC = C - B;
			Eigen::Vector2f CD = D - C;
			Eigen::Vector2f DA = A - D;

			for(int pixel_i = xlb; pixel_i <= xub; pixel_i++){
				for(int pixel_j = ylb; pixel_j <= yub; pixel_j++){
					Eigen::Vector2f temp_pixel = Eigen::Vector2f((float)pixel_i, (float)pixel_j);
					Eigen::Vector2f AP = temp_pixel - A;
					Eigen::Vector2f BP = temp_pixel - B;
					Eigen::Vector2f CP = temp_pixel - C;
					Eigen::Vector2f DP = temp_pixel - D;
					if(SIGN(cross_product(AB, AP)) * SIGN(cross_product(BC, BP)) * SIGN(cross_product(CD, CP)) * SIGN(cross_product(DA, DP)) == 1){
						PerPixelIntersectionList[pixel_i][pixel_j].push_back(i);
					}
				}
			}
		}
	}
}


void ComposeIntersectionEffects(Intersection_effect *list) {
	Eigen::Vector3f c_color;
	float c_opacity;
	for (int i = 0; i < list.size(); i++) {
		Eigen::Vector3f r_color = *(list + i).color;
		float r_float = *(list + i).opacity;
		c_color += (1 - c_opacity) * r_color;
		c_opacity += (1 - c_opacity) * r_opacity;
		if (c_opacity >= 0.9999f) {
			break;
		}
	}
	return c_color;
}


int main()
{
	/*
	 * 1. Volume Setting
	 */
	
	Volume vol("data/test2.bin");
	/*
	 * 2. Camera Setting
	 */
	Eigen::Vector3f cameraPosition= vol.bbox.getCenter()-2.5*Eigen::Vector3f(0,0,vol.size_physics.z());
	Eigen::Vector3f cameraLookAt= vol.bbox.getCenter();
	Eigen::Vector3f cameraUp(0, 1, 0);
	float verticalFov = 45;
	Eigen::Vector2i filmRes(1024, 1024);

	Camera camera(cameraPosition, cameraLookAt, cameraUp, verticalFov, filmRes);
	PointLight light(cameraPosition+Eigen::Vector3f(0,5,0),Eigen::Vector3f(1,1,1));


	/*******************************************************************************************************/
	/*
	x is the pixel in width
	y is the pixel in heigtht
	z is a vector of int, storing the index of tetrahedron projected above the pixel
	*/
	std::vector<std::vector<std::vector<int>>> PerPixelIntersectionList;
	//initialize the 3d vector
	PerPixelIntersectionList.resize(WIDTH);
	for(int i = 0; i<WIDTH; i++){
		PerPixelIntersectionList[i].resize(HEIGHT);
	}
	
	std::vector<Eigen::Vector2f> SSC = ComputeScreenSpaceProjections(vol.raw_data, camera);
	
	ExtractIntersectionRecords(tetra_list, vertex_list, SSC, PerPixelIntersectionList);

	return 0;
}