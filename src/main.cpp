﻿#include <algorithm>
#include "volume.h"
#include "camera.hpp"
#include "light.hpp"
#include <iostream>
#include <omp.h>
#include "include/colormap/colormap.h"
#include "utils.hpp"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define TETRA_ID 492
#define WIDTH 1024
#define HEIGHT 1024

#define DISTCONST 0.2f
#define NUM_Samples 3

#define MIN4(a,b,c,d) (((a)<(b)?(a):(b))<((c)<(d)?(c):(d))?((a)<(b)?(a):(b)):((c)<(d)?(c):(d)))
#define MAX4(a,b,c,d) (((a)>(b)?(a):(b))>((c)>(d)?(c):(d))?((a)>(b)?(a):(b)):((c)>(d)?(c):(d)))
#define SIGN(x) (((x) > 0)?1:-1)

using namespace colormap;

int RayTetraIntersection(Eigen::Vector3f & ip0, Eigen::Vector3f & ip1,
						Tetrahedron tetra,
						Camera *camera, Eigen::Vector2i pixel_idx,
						std::vector<MyVertex> *Allvertices);

float InterpolateScalar(Tetrahedron t, Eigen::Vector3f p,
						std::vector<MyVertex> *Allvertices);

int intersect_triangle(Eigen::Vector3f & ip,
						Eigen::Vector3f A, Eigen::Vector3f B, Eigen::Vector3f C,
						Ray ray);

/* project vertices in 3D to screen space in 2d */
void ComputeScreenSpaceProjections(std::vector<Eigen::Vector2f> & SSC,
									std::vector<MyVertex> *Vertices, Camera *camera){

	for(int i=0; i<Vertices->size(); ++i){
		Ray ray(camera->m_Pos, Vertices->at(i).coordinate - camera->m_Pos);
		if(ray.m_Dir.dot(camera->m_Forward) <= 0){
			// vertex is behind the camera
			SSC.push_back(Eigen::Vector2f(-1.0f, -1.0f)); // means that the vertex has an invalid SSC
			continue;
		}
		Eigen::Vector3f p = (1.0f / ray.m_Dir.dot(camera->m_Forward))*ray.m_Dir - camera->m_Forward;

		float y_idx = p.dot(camera->m_Up)/camera->m_Up.squaredNorm() * camera->m_Film.m_Res.y()/2.0f + camera->m_Film.m_Res.y()/2.0f;

		SSC.push_back(Eigen::Vector2f(p.dot(camera->m_Right)/camera->m_Right.squaredNorm() * camera->m_Film.m_Res.x()/2.0f + camera->m_Film.m_Res.x()/2.0f,
									  y_idx));
	}
}

/* get the projection of each tetrahedron, get the piexl projected by the tetrahedron, build a intersection list for each pixel */
void ExtractIntersectionRecords(std::vector<Tetrahedron>* tetra_list, std::vector<Eigen::Vector2f>* SSC, std::vector<std::vector<std::vector<int>>>& PerPixelIntersectionList){
	int xmin = 1023;
	int ymin = 1023;
	int xmax = 0;
	int ymax = 0;
	// iterate each tetrahedron
	for(int i = 0; i < tetra_list->size(); i++){
		
		/* the four projected points on screen */
		Eigen::Vector2f v1_proj = SSC->at(tetra_list->at(i).v1_idx);
		Eigen::Vector2f v2_proj = SSC->at(tetra_list->at(i).v2_idx);
		Eigen::Vector2f v3_proj = SSC->at(tetra_list->at(i).v3_idx);
		Eigen::Vector2f v4_proj = SSC->at(tetra_list->at(i).v4_idx);

		if(v1_proj.x() < 0.0f || v1_proj.y() < 0.0f || v1_proj.x() >1023.0f || v1_proj.y() >1023.0f ||
		   v2_proj.x() < 0.0f || v2_proj.y() < 0.0f || v2_proj.x() >1023.0f || v2_proj.y() >1023.0f ||
		   v3_proj.x() < 0.0f || v3_proj.y() < 0.0f || v3_proj.x() >1023.0f || v3_proj.y() >1023.0f ||
		   v4_proj.x() < 0.0f || v4_proj.y() < 0.0f || v4_proj.x() >1023.0f || v4_proj.y() >1023.0f){
			   continue;
		   }
		/* lower bound and upper bound of x and y */
		int xlb = std::max((int)std::ceil(MIN4(v1_proj.x(), v2_proj.x(), v3_proj.x(), v4_proj.x())), 0);
		int ylb = std::max((int)std::ceil(MIN4(v1_proj.y(), v2_proj.y(), v3_proj.y(), v4_proj.y())), 0);

		int xub = std::min((int)std::floor(MAX4(v1_proj.x(), v2_proj.x(), v3_proj.x(), v4_proj.x())), 1023);
		int yub = std::min((int)std::floor(MAX4(v1_proj.y(), v2_proj.y(), v3_proj.y(), v4_proj.y())), 1023);

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
			
			//now that the clockwise order of 4 vertices is leftmost, r1, rightmost, l1
			Eigen::Vector2f A, B, C, D;
			A = leftmost_point;
			B = r1;
			C = rightmost_point;
			D = l1;

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
					if(cross_product(AB, AP)>0 && cross_product(BC, BP)>0 && cross_product(CD, CP)>0 && cross_product(DA, DP)>0){
						PerPixelIntersectionList[pixel_i][pixel_j].push_back(i);
					}
				}
			}
		}
	}

}

void CalculateIntersectionEffect(std::vector<Intersection_effect> & effectlist_for_this_pixel /* filled in this function */, 
								Camera *camera, 
								Eigen::Vector2i pixel_idx,
								std::vector<int> *Intersectionlist_for_this_pixel,
								std::vector<Tetrahedron> *Alltetra,
								std::vector<MyVertex> *Allvertices,
								int NumOfSamples = 3)
{
	
	for(int tetra_iter=0; tetra_iter < Intersectionlist_for_this_pixel->size(); ++tetra_iter){
		
		Intersection_effect record;
		Eigen::Vector3f ip0, ip1;
		
		int exist = RayTetraIntersection(ip0, ip1, Alltetra->at(Intersectionlist_for_this_pixel->at(tetra_iter)), camera, pixel_idx, Allvertices);
		if(exist==0) continue;
		
		record.dist = (ip0 - camera->m_Pos).norm();
		Eigen::Vector3f d = (ip1-ip0) / NumOfSamples;
		
		Eigen::Vector3f _color(0.0f,0.0f,0.0f); float _opacity = 0.0f;
		
		record.color = _color; 
		record.opacity = _opacity;
		
		for(int i=0; i<NumOfSamples; ++i){
			Eigen::Vector3f ip = ip0 + d*i;
			float s = InterpolateScalar(Alltetra->at(Intersectionlist_for_this_pixel->at(tetra_iter)), ip, Allvertices); // interpolated density of the sampled point
			
			transform::HotMetal transform;
			IDL::Prism idl;
			MATLAB::Hsv matlab;
			Color mc = idl.getColor(s);
			_color = Eigen::Vector3f((float)mc.r, (float)mc.g, (float)mc.b);
			_color = s * _color;
			record.color += DISTCONST * d.norm() * (1-record.opacity) * _color;
			record.opacity += DISTCONST * d.norm() * (1-record.opacity) * (float)mc.a; // opacity_src
		}
		effectlist_for_this_pixel.push_back(record);
	}
}

int RayTetraIntersection(Eigen::Vector3f & ip0, Eigen::Vector3f & ip1,
						Tetrahedron tetra,
						Camera *camera, Eigen::Vector2i pixel_idx,
						std::vector<MyVertex> *Allvertices)
{
	Ray ray = camera->generateRay((float)pixel_idx.x(), (float)pixel_idx.y());
	std::vector<Eigen::Vector3f> reg, reg2, finalreg;
	int a[4];
	Eigen::Vector3f p1,p2,p3,p4, regp;
	reg.resize(4);

	a[0] = intersect_triangle(reg[0], (Allvertices->at(tetra.v1_idx)).coordinate,
									  (Allvertices->at(tetra.v2_idx)).coordinate,
									  (Allvertices->at(tetra.v3_idx)).coordinate, ray);
	a[1] = intersect_triangle(reg[1], (Allvertices->at(tetra.v1_idx)).coordinate,
									  (Allvertices->at(tetra.v2_idx)).coordinate,
									  (Allvertices->at(tetra.v4_idx)).coordinate, ray);
	a[2] = intersect_triangle(reg[2], (Allvertices->at(tetra.v1_idx)).coordinate,
									  (Allvertices->at(tetra.v3_idx)).coordinate,
									  (Allvertices->at(tetra.v4_idx)).coordinate, ray);
	a[3] = intersect_triangle(reg[3], (Allvertices->at(tetra.v2_idx)).coordinate,
									  (Allvertices->at(tetra.v3_idx)).coordinate,
									  (Allvertices->at(tetra.v4_idx)).coordinate, ray);
	for(int i=0; i<4; ++i){
		if(a[i]){
			reg2.push_back(reg[i]);
		}
	}
	if(reg2.size()<=1) return 0;
	for(int k=0; k<reg2.size(); ++k){
		if(finalreg.empty()){
			finalreg.push_back(reg2[k]);
		}else{
			if(reg2[k]!=finalreg[0]){
				finalreg.push_back(reg2[k]);
			}
		}
	}

	float dist2eye[2];
	for(int i=0; i<finalreg.size(); ++i){
		dist2eye[i] = (finalreg[i]-ray.m_Ori).norm();
	}
	if(dist2eye[0] < dist2eye[1]){
		ip0 = finalreg[0]; ip1 = finalreg[1];
	}else{
		ip0 = finalreg[1]; ip1 = finalreg[0];
	}
	return 1;

}

float InterpolateScalar(Tetrahedron t, Eigen::Vector3f p,
						std::vector<MyVertex> *Allvertices)
{
	std::vector<Eigen::Vector3f> M;
	Eigen::Vector3f N, R, Scalar;
	float s;
	{
		M.push_back((Allvertices->at(t.v2_idx)).coordinate - (Allvertices->at(t.v1_idx)).coordinate);
		M.push_back((Allvertices->at(t.v3_idx)).coordinate - (Allvertices->at(t.v1_idx)).coordinate);
		M.push_back((Allvertices->at(t.v4_idx)).coordinate - (Allvertices->at(t.v1_idx)).coordinate);
	}
	N = p - (Allvertices->at(t.v1_idx)).coordinate;
	{
		Scalar.x() = (Allvertices->at(t.v2_idx)).density - (Allvertices->at(t.v1_idx)).density;
		Scalar.y() = (Allvertices->at(t.v3_idx)).density - (Allvertices->at(t.v1_idx)).density;
		Scalar.z() = (Allvertices->at(t.v4_idx)).density - (Allvertices->at(t.v1_idx)).density;
	}

	Eigen::Matrix3f _M;
	_M << M[0].x(), M[1].x(), M[2].x(),
		  M[0].y(), M[1].y(), M[2].y(),
		  M[0].z(), M[1].z(), M[2].z();
	R = _M.inverse() * Scalar;
	s = R.dot(N) + (Allvertices->at(t.v1_idx)).density;
	return s;
}

Eigen::Vector3f ComposeIntersectionEffects(std::vector<Intersection_effect> *list) {
	Eigen::Vector3f c_color = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
	float c_opacity = 0.0f;
	for (int i = 0; i < list->size(); i++) {
		Eigen::Vector3f r_color = list->at(i).color;
		float r_opacity = list->at(i).opacity;
		c_color += (1 - c_opacity) * r_color;
		c_opacity += (1 - c_opacity) * r_opacity;
		if (c_opacity >= 0.9999f) {
			break;
		}
	}
	return c_color;
}

bool SortFunc(const Intersection_effect& ef_a, const Intersection_effect& ef_b){
	return ef_a.dist < ef_b.dist;
}


int main()
{	
	std::cout << "Start Render" <<std::endl;
	std::vector<Tetrahedron> tetra_list;
	std::vector<MyVertex> vertex_list;

	/*
	 * 1. Volume Setting
	 */
	Volume vol(tetra_list, vertex_list, "data/post.txt");

	/*
	 * 2. Camera Setting
	 */
	
	Eigen::Vector3f cameraPosition= vol.bbox.getCenter()-2.5*Eigen::Vector3f(0.3,0.7,vol.size_physics.z());
	Eigen::Vector3f cameraLookAt= vol.bbox.getCenter();
	Eigen::Vector3f cameraUp(0, 1, 0);
	float verticalFov = 45;
	Eigen::Vector2i filmRes(WIDTH, HEIGHT);

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
	
	std::vector<Eigen::Vector2f> SSC;
	ComputeScreenSpaceProjections(SSC, &vertex_list, &camera);
	ExtractIntersectionRecords(&tetra_list, &SSC, PerPixelIntersectionList);
	
	#pragma omp parallel for
	for (int i = 0; i < camera.m_Film.m_Res.x(); i++){
		for(int j = 0; j < camera.m_Film.m_Res.y(); j++){
			std::vector<int> list = PerPixelIntersectionList[i][j];
			if(i%100==0 && j%512==0 ) std::cout<< "i, j, list size: " << i << " " << j << " " << list.size() << std::endl;
			if(list.size() == 0){
				camera.setPixel(i, j, Eigen::Vector3f(0.0f, 0.0f, 0.0f));
			} else {
			
				std::vector<Intersection_effect> IntEffectList;
			
				CalculateIntersectionEffect(IntEffectList, &camera, Eigen::Vector2i(i, j), &list, &tetra_list, &vertex_list, NUM_Samples);
				
				std::sort(IntEffectList.begin(), IntEffectList.end(), SortFunc);
				
				Eigen::Vector3f color_c = ComposeIntersectionEffects(&IntEffectList);
				
				camera.setPixel(i, j, color_c);
			}

		}
	}
	std::cout << "End Render\n";
	camera.m_Film.write("./tetra_vol.png");
	return 0;
}