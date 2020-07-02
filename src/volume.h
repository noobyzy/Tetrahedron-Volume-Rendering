#pragma once
#include <vector>
#include <string>
#include "Eigen/Dense"
#include <iostream>
#include <fstream>
#include "aabb.hpp"
#include "volumeData.h"
#include "voxel.h"
#include "ray.hpp"
#include "tetra.h"
class Volume{
    
public:
    std::vector<MyVertex> raw_data;
    std::vector<Tetrahedron> tetra_data;
    Eigen::Vector3i size;Eigen::Vector3f size_physics;
    AABB bbox ;
    double dx;
    float grad_maxnorm;
    Volume();
    ~Volume();
    Volume(std::string volumefile);
    bool getRayStartEnd(Ray& ray,float& t_start,float & t_end);
    int indexToData(Eigen::Vector3i index);
    void getVoxel(int x_idx, int y_idx, int z_idx);
};