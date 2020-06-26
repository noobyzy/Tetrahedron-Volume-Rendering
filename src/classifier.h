#pragma once
#include "volume.h"
#include "camera.hpp"
#include <iostream>
#include <cmath>
//https://github.com/yuki-koyama/tinycolormap
#include "tinycolormap.hpp"

class Classifier{
    public:
    virtual opticsData transfer(volumeData v_data,float dt,Camera* cam ,std::vector<Light *> lights,Eigen::Vector3f raydir, float grad_max_norm=1)=0;
};

class myClassifier:public Classifier{
    public:
    opticsData transfer(volumeData v_data,float dt,Camera* cam,std::vector<Light *> lights, Eigen::Vector3f raydir, float grad_max_norm=1){
        opticsData optics;
		//Write your own classifier, make use of input args(may not all)
		// should set transparency and color of optics
        tinycolormap::Color color = tinycolormap::GetColor(v_data.density, tinycolormap::ColormapType::Jet);
        // Phong lighting model
        float shininess = 16.0f;
        float lightPower = 15.0f;
        Eigen::Vector3f ambientColor(0.5,0.5,0.5);
        Eigen::Vector3f diffuseColor(0.6,0.6,0.6);
        Eigen::Vector3f specularColor(1.0, 1.0, 1.0);
        Eigen::Vector3f ambient, diffuse, specular;
        Eigen::Vector3f lightDir = (-(v_data.position) + lights[0]->m_Pos).normalized();
        Eigen::Vector3f normal = v_data.gradient.normalized();
        
        ambient = ambientColor;

        diffuse = std::max(normal.dot(lightDir), 0.0f) * diffuseColor;

        Eigen::Vector3f reflectDir = (2 * lightDir.dot(normal) * normal - lightDir).normalized();
        float temp = reflectDir.dot(raydir);
        float spec = powf(temp, shininess) / 3;
        specular = spec * specularColor;

        optics.color = Eigen::Vector3f(color.r(), color.g(), color.b()).cwiseProduct(
            ambient + diffuse + specular) * dt * v_data.density;

        optics.transparency = exp(-v_data.density * dt) * Eigen::Vector3f::Ones();

        return optics;
    }
};