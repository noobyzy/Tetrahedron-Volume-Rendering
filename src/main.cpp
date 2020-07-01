
#include <iostream>
#include <algorithm>
#include <tetra.h>



#include"renderer.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define WIDTH 1024
#define HEIGHT 1024
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

void ExtractIntersectionRecords(std::vector<Tetrahedron> tetra_list, std::vector<MyVertex> vertex_list, std::vector<Eigen::Vector2f> SSC, std::vector<std::vector<std::vector<int>>>& PerPixelIntersectionList){

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
		RayTetraIntersection(ip0, ip1, Alltetra->at(tetra_iter), camera, pixel_idx, Allvertices);
		record.dist = (ip0 - camera->m_Pos).norm();
		Eigen::Vector3f d = (ip1-ip0) / NumOfSamples;
		Eigen::Vector3f _color(0.0,0.0,0.0); float _opacity = 0.0;
		record.color = _color; record.opacity = _opacity;
		for(int i=0; i<NumOfSamples; ++i){
			Eigen::Vector3f ip = ip0 + d*i;
			float s = InterpolateScalar(Alltetra->at(tetra_iter), ip, Allvertices);
			tinycolormap::Color tinycolor = tinycolormap::GetColor(s, tinycolormap::ColormapType::Jet);
			_color.x() = tinycolor.r(); _color.y() = tinycolor.g(); _color.z() = tinycolor.b();
			
		}

	}
}

void RayTetraIntersection(Eigen::Vector3f & ip0, Eigen::Vector3f & ip1,
						Tetrahedron tetra,
						Camera *camera, Eigen::Vector2i pixel_idx,
						std::vector<MyVertex> *Allvertices)
{
	Ray ray = camera->generateRay((float)pixel_idx.x(), (float)pixel_idx.y());
	std::vector<Eigen::Vector3f> reg, reg2;
	int a[4];
	Eigen::Vector3f p1,p2,p3,p4, regp;
	reg.push_back(p1); reg.push_back(p2); reg.push_back(p3); reg.push_back(p4);
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
	float dist2eye[2];
	for(int i=0; i<2; ++i){
		dist2eye[i] = (reg2[i]-ray.m_Ori).norm()
	}
	if(dist2eye[0] < dist2eye[1]){
		ip0 = reg2[0]; ip1 = reg2[1];
	}else{
		ip0 = reg2[1]; ip1 = reg2[0];
	}
	return;

}

int intersect_triangle(Eigen::Vector3f & ip,
						Eigen::Vector3f A, Eigen::Vector3f B, Eigen::Vector3f C,
						Ray ray)
{
	Eigen::Vector3f E1 = B-A;
	Eigen::Vector3f E2 = C-A;
	Eigen::Vector3f N = E1.cross(E2);
	float det = -ray.m_Dir.dot(N);
	float invdet = 1.0/det;
	Eigen::Vector3f A0 = ray.m_Ori - A;
	Eigen::Vector3f DAO = A0.cross(ray.m_Dir);
	float u = E2.dot(DAO) * invdet;
	float v = -E1.dot(DAO) * invdet;
	float t = A0.dot(N) * invdet;
	if(det>=ray.m_fMin && t>= 0.0 && u>=0.0 && v>=0.0 && (u+v)<=0.0){
		ip = A + u*E1 + v*E2;
		return 1;
	}
	return 0;
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


	std::vector<std::vector<std::vector<int>>> PerPixelIntersectionList;
	
	std::vector<Eigen::Vector2f> SSC = ComputeScreenSpaceProjections(vol.raw_data, camera);

	


	return 0;
}