#include "tetra.h"

float cross_product(Eigen::Vector2f v1, Eigen::Vector2f v2){
	return (v1.x()*v2.y() - v2.x()*v1.y());
}

/* test the side of the point */
int line_side(Eigen::Vector2f va, Eigen::Vector2f vb, Eigen::Vector2f p){

	Eigen::Vector2f PA = va - p;
	Eigen::Vector2f AB = va - vb;
	float temp = cross_product(PA, AB);
	if(temp > 0){
		return 1;
	}else if(temp == 0){
		return 0;
	}else{
		return -1;
	}
	return 1;
}

/* to justify whether a point is outside the triangle */
bool is_inside_triangle(Eigen::Vector2f va, Eigen::Vector2f vb, Eigen::Vector2f vc, Eigen::Vector2f p){
	bool line_c_test, line_b_test, line_a_test;
	line_c_test = (line_side(va, vb, p) == line_side(va, vb, vc));
	line_b_test = (line_side(va, vc, p) == line_side(va, vc, vb));
	line_a_test = (line_side(vc, vb, p) == line_side(vc, vb, va));
	return (line_c_test && line_b_test && line_a_test);
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
	if(t>= 0.0 && u>=0.0 && v>=0.0 && (u+v)<=1.0){
		ip = A + u*E1 + v*E2;
		return 1;
	}
	return 0;
}