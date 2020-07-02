﻿#include <algorithm>
#include "volume.h"
#include "opticsData.h"
#include "compositor.h"
#include "classifier.h"
#include "camera.hpp"
#include "light.hpp"
#include "classifier.h"
#include "tetra.h"
#include <iostream>
#include <omp.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define WIDTH 1024
#define HEIGHT 1024

#define MIN4(a,b,c,d) (((a)<(b)?(a):(b))<((c)<(d)?(c):(d))?((a)<(b)?(a):(b)):((c)<(d)?(c):(d)))
#define MAX4(a,b,c,d) (((a)>(b)?(a):(b))>((c)>(d)?(c):(d))?((a)>(b)?(a):(b)):((c)>(d)?(c):(d)))
#define SIGN(x) (((x) > 0)?1:-1)

void RayTetraIntersection(Eigen::Vector3f & ip0, Eigen::Vector3f & ip1,
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
	float xmin = 12345.0f;
	float ymin = 12345.0f;
	float zmin = 12345.0f;
	float xmax = 0.0f;
	float ymax = 0.0f;
	float zmax = 0.0f;


	for(int i=0; i<Vertices->size(); ++i){
		if(Vertices->at(i).coordinate.x() < xmin) xmin = Vertices->at(i).coordinate.x();
		if(Vertices->at(i).coordinate.y() < ymin) ymin = Vertices->at(i).coordinate.y();
		if(Vertices->at(i).coordinate.z() < zmin) zmin = Vertices->at(i).coordinate.z();
		if(Vertices->at(i).coordinate.x() > xmax) xmax = Vertices->at(i).coordinate.x();
		if(Vertices->at(i).coordinate.y() > ymax) ymax = Vertices->at(i).coordinate.y();
		if(Vertices->at(i).coordinate.z() > zmax) zmax = Vertices->at(i).coordinate.z();



		Ray ray(camera->m_Pos, Vertices->at(i).coordinate - camera->m_Pos);
		Eigen::Vector3f p = (1.0f / ray.m_Dir.dot(camera->m_Forward))*ray.m_Dir - camera->m_Forward;
		SSC.push_back(Eigen::Vector2f(p.dot(camera->m_Right)/camera->m_Right.squaredNorm() * camera->m_Film.m_Res.x()/2.0f,
									  p.dot(camera->m_Up)/camera->m_Up.squaredNorm() * camera->m_Film.m_Res.y()/2.0f));
	}
	std::cout<<"xmin="<<xmin<<std::endl;
	std::cout<<"ymin="<<ymin<<std::endl;
	std::cout<<"zmin="<<zmin<<std::endl;
	std::cout<<"xmax="<<xmax<<std::endl;
	std::cout<<"ymax="<<ymax<<std::endl;
	std::cout<<"zmax="<<zmax<<std::endl;



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

		/* lower bound and upper bound of x and y */
		int xlb = std::max((int)std::floor(MIN4(v1_proj.x(), v2_proj.x(), v3_proj.x(), v4_proj.x())), 0);
		int ylb = std::max((int)std::floor(MIN4(v1_proj.y(), v2_proj.y(), v3_proj.y(), v4_proj.y())), 0);

		int xub = std::min((int)std::ceil(MAX4(v1_proj.x(), v2_proj.x(), v3_proj.x(), v4_proj.x())), 1023);
		int yub = std::min((int)std::ceil(MAX4(v1_proj.y(), v2_proj.y(), v3_proj.y(), v4_proj.y())), 1023);

		if(xlb<xmin) xmin = xlb;
		if(ylb<ymin) ymin = ylb;
		if(xub>xmax) xmax = xub;
		if(yub>ymax) ymax = yub;
		

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
	std::cout<<"xmin="<<xmin<<std::endl;
	std::cout<<"ymin="<<ymin<<std::endl;
	std::cout<<"xmax="<<xmax<<std::endl;
	std::cout<<"ymax="<<ymax<<std::endl;

}

void CalculateIntersectionEffect(std::vector<Intersection_effect> & effectlist_for_this_pixel /* filled in this function */, 
								Camera *camera, 
								Eigen::Vector2i pixel_idx,
								std::vector<int> *Intersectionlist_for_this_pixel,
								std::vector<Tetrahedron> *Alltetra,
								std::vector<MyVertex> *Allvertices,
								int NumOfSamples = 3)
{
	float DISTCONST = 0.5;
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
			float s = InterpolateScalar(Alltetra->at(tetra_iter), ip, Allvertices); // interpolated density of the sampled point
			tinycolormap::Color tinycolor = tinycolormap::GetColor(s, tinycolormap::ColormapType::Jet);
			_color.x() = tinycolor.r(); _color.y() = tinycolor.g(); _color.z() = tinycolor.b();
			record.color += DISTCONST * d.norm() * (1-record.opacity) * _color;
			record.opacity += DISTCONST * d.norm() * (1-record.opacity) * (1-exp(-s*d.norm())); // opacity_src
		}
		effectlist_for_this_pixel.push_back(record);
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
		dist2eye[i] = (reg2[i]-ray.m_Ori).norm();
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
	if(det>=ray.m_fMin && t>= 0.0 && u>=0.0 && v>=0.0 && (u+v)<=1.0){
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
	std::vector<Tetrahedron> tetra_list;
	std::vector<MyVertex> vertex_list;

	/*
	 * 1. Volume Setting
	 */
	
	
	Volume vol("data/test2.bin", tetra_list, vertex_list);
	/*
	 * 2. Camera Setting
	 */
	
	Eigen::Vector3f cameraPosition= vol.bbox.getCenter()-2.5*Eigen::Vector3f(0,0,vol.size_physics.z());
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
	SSC.reserve(vertex_list.size());
	ComputeScreenSpaceProjections(SSC, &vertex_list, &camera);
	
	ExtractIntersectionRecords(&tetra_list, &SSC, PerPixelIntersectionList);
	

	for (int i = 0; i < camera.m_Film.m_Res.x(); i++){
		for(int j = 0; j < camera.m_Film.m_Res.y(); j++){
			std::vector<int> list = PerPixelIntersectionList[i][j];
			
			std::vector<Intersection_effect> IntEffectList;
			CalculateIntersectionEffect(IntEffectList, &camera, Eigen::Vector2i(i, j), &list, &tetra_list, &vertex_list, 3);
			std::sort(IntEffectList.begin(), IntEffectList.end(), SortFunc);
			Eigen::Vector3f color_c = ComposeIntersectionEffects(&IntEffectList);
			camera.setPixel(i, j, color_c);
		}
	}
	camera.m_Film.write("./tetra_vol.png");
	return 0;
}