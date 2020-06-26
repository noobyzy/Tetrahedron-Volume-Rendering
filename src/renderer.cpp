#include "renderer.h"
void Renderer::setCamera(Camera * cam){this->camera=cam;};
void Renderer::addLight(Light * light){this->lights.push_back(light);};
void Renderer::setVolume(Volume * vol){this->volume=vol;};
void Renderer::setClassifier(Classifier * classifier){this->classifier=classifier;};
void Renderer::setInterpolator(Interpolator * interpolator){this->interpolator=interpolator;};
void Renderer::renderFrontToBack(std::string imgPath){
#pragma omp parallel for
    for (int i = 0; i < camera->m_Film.m_Res.x() * camera->m_Film.m_Res.y(); i++) {
        int dy = i / camera->m_Film.m_Res.x();
        int dx = i - dy * camera->m_Film.m_Res.x();
        Eigen::Vector3f color(0,0,0);Eigen::Vector3f alpha(0,0,0);
        Ray ray = camera->generateRay((float)dx , (float)dy );
        float t_start; float t_end;
        //Intersection calculate
        if(this->volume->getRayStartEnd(ray,t_start,t_end)){
			
			//Here your render code 
            float step = 0.01f;
            Compositor comp;
            for (float t = t_start; t <= t_end; t += step) {
                Eigen::Vector3f samplePos = ray.m_Ori + t * ray.m_Dir;
                volumeData sampleData;
                sampleData = interpolator->interpolate(samplePos, volume->getVoxel(samplePos));
                opticsData sampleOpt = classifier->transfer(sampleData, step, camera,
                    lights, ray.m_Dir);
                comp.compositeFrontToBack(color, alpha, sampleOpt.getColor(), sampleOpt.getOpacity());
            }
        }
        camera->setPixel(dx, dy, color );
    }
    camera->m_Film.write(imgPath);
};
void Renderer::renderBackToFront(std::string imgPath){
#pragma omp parallel for
    for (int i = 0; i < camera->m_Film.m_Res.x() * camera->m_Film.m_Res.y(); i++) {
        int dy = i / camera->m_Film.m_Res.x();
        int dx = i - dy * camera->m_Film.m_Res.x();
        Eigen::Vector3f color(0,0,0);
        Ray ray = camera->generateRay((float)dx , (float)dy );
        float t_start; float t_end;
        //Intersection calculate
        if(this->volume->getRayStartEnd(ray,t_start,t_end)){
			
			//Here your render code 
        }
        camera->setPixel(dx, dy, color );
    }
    camera->m_Film.write(imgPath);
};