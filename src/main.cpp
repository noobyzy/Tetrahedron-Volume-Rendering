
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