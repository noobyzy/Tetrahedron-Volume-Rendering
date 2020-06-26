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
};

class Tetrahedron{
public:
    int v1_idx;
    int v2_idx;
    int v3_idx;
    int v4_idx;
    Tetrahedron(int v1, int v2, int v3, int v4):v1_idx(v1), v2_idx(v2), v3_idx(v3), v4_idx(v4){/* default constructor */}
};

class Intersection_record{
    public:
    int pixel_width_idx;
    int pixel_height_idx;
    int tetra_idx;
    Intersection_record(int w, int h, int t):pixel_width_idx(w), pixel_height_idx(h), tetra_idx(t){ /* DC */}
};

class Intersection_effect{
    float dist; // the eye distance to the first tetrahedron
    Eigen::Vector3f color;
    float opacity;
};



