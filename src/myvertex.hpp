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
class MyVertex{
    public:
    Eigen::Vector3f coordinate;
    float density;
    MyVertex(float x, float y, float z, float den):coordinate(Eigen::Vector3f(x,y,z)), density(den){/* default constructor*/}
}