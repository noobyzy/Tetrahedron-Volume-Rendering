#include "volume.h"
#include <iostream>
#include <fstream>
// using namespace std;
template <typename T>
void SwapEnd(T &var)
{
    char *varArray = reinterpret_cast<char *>(&var);
    for (long i = 0; i < static_cast<long>(sizeof(var) / 2); i++)
        std::swap(varArray[sizeof(var) - 1 - i], varArray[i]);
}

Volume::Volume(){};
Volume::~Volume(){};

Volume::Volume(std::vector<Tetrahedron> & tetra_data, std::vector<MyVertex> & vertex_list)
{
    std::fstream volumefile;

    volumefile.open("data/data.txt", std::ios::in);

    int trash;
    std::string str;
    int vertex_num, tetra_num, v1_buffer, v2_buffer, v3_buffer, v4_buffer;
    float x_buffer, y_buffer, z_buffer, density_buffer, density_max=-1,
          xmin=400, xmax=-400, ymin=400, ymax=-400, zmin=400, zmax=-400;
    std::vector<float> density_list;

    volumefile >> vertex_num;
    std::cout << "vertex_num: " << vertex_num << std::endl;
    density_list.reserve(vertex_num);
    volumefile >> tetra_num;
    std::cout << "tetra_num: " << tetra_num << std::endl;
    for (int i = 0; i < 3; i++) {
        volumefile >> trash;
        std::cout << trash << std::endl;
    }
    for(int i=0; i<vertex_num; ++i){
        volumefile >> trash; // vertex index
        volumefile >> x_buffer;
        volumefile >> y_buffer;
        volumefile >> z_buffer;
        if(xmin > x_buffer) xmin = x_buffer;
        if(xmax < x_buffer) xmax = x_buffer;
        if(ymin > y_buffer) ymin = y_buffer;
        if(ymax < y_buffer) ymax = y_buffer;
        if(zmin > z_buffer) zmin = z_buffer;
        if(zmax < z_buffer) zmax = z_buffer;
        vertex_list.push_back(MyVertex(x_buffer, y_buffer, z_buffer, 0.0));
    }
    for(int i=0; i<tetra_num; ++i){
        volumefile >> trash; // tetra index
        volumefile >> trash; // 1
        volumefile >> str;   // tet
        volumefile >> v1_buffer;
        volumefile >> v2_buffer;
        volumefile >> v3_buffer;
        volumefile >> v4_buffer;
        tetra_data.push_back(Tetrahedron(v1_buffer-1, v2_buffer-1, v3_buffer-1, v4_buffer-1));
    }
    volumefile >> trash; //1
    volumefile >> trash; //1
    volumefile >> str; // density
    for(int i=0; i<vertex_num; ++i){
        volumefile >> trash; // vertex index
        volumefile >> density_buffer;
        if(density_max < density_buffer) density_max = density_buffer;
        vertex_list[i].coordinate = Eigen::Vector3f(vertex_list[i].coordinate.x()-xmin,
                                                    vertex_list[i].coordinate.y()-ymin,
                                                    vertex_list[i].coordinate.z()-zmin);    
        density_list.push_back(density_buffer);
    }
    for(int i=0; i<vertex_num; ++i){
        vertex_list[i].density = density_list[i] / density_max;
    }
    this->bbox = AABB(0,0,0, xmax-xmin, ymax-ymin, zmax-zmin);
    this->size_physics = Eigen::Vector3f(xmax-xmin, ymax-ymin, zmax-zmin);
    // for (int i = 0; i < 3; i++) {
    //     volumefile >> ch;
    //     std::cout << ch;
    // }
    //     volumefile >> n;
    //     std::cout << n << std::endl;


    volumefile.close();
    return;

    // if (fp)
    // {    
    //     fread(this->size.data(),sizeof(int),3,fp);
    //     fread(this->size_physics.data(),sizeof(float),3,fp);
    //     this->bbox = AABB(0, 0, 0, size_physics.x(), size_physics.y(), size_physics.z());
    //     this->dx = size_physics.x() / (this->size.x()-1);
    //     int count = this->size.x() * this->size.y() * this->size.z();
    //     // this->raw_data.clear();
    //     float den = 0;
    //     float min_den = 9999;
    //     float max_den = -9999;
    //     for (int curind = 0; curind < count; curind++)
    //     {
    //         int z_idx = curind / (this->size.y() * this->size.x());
    //         int y_idx = (curind - z_idx * this->size.x() * this->size.y()) / this->size.x();
    //         int x_idx = curind - z_idx * this->size.x() * this->size.y() - y_idx * this->size.x();
    //         fread(&den, sizeof(float), 1, fp);
    //         min_den = std::min(min_den,den);
    //         max_den = std::max(max_den,den);
    //         vertex_list.push_back(MyVertex(this->dx * x_idx, this->dx * y_idx, this->dx * z_idx, den));
    //         if(x_idx<this->size.x() - 1 && y_idx<this->size.y() - 1 && z_idx<this->size.z() - 1){
    //             getVoxel(x_idx, y_idx, z_idx, tetra_data);
    //         }
    //     }
    //     //printf("density range [%.3f ,%.3f]\n",min_den,max_den);
    // }
    // fclose(fp);
};


bool Volume::getRayStartEnd(Ray &ray, float &t_start, float &t_end)
{
    return bbox.rayIntersection(ray, t_start, t_end);
};

int Volume::indexToData(Eigen::Vector3i index)
{
    return index.z() * this->size.x() * this->size.y() + index.y() * this->size.x() + index.x();
}
void Volume::getVoxel(int x_idx, int y_idx, int z_idx, std::vector<Tetrahedron> & tetra_data)
{
    Eigen::Vector3i base(x_idx, y_idx, z_idx);
    int v1 = indexToData(Eigen::Vector3i(x_idx, y_idx, z_idx));

    // tetrehedraon 1
    int v2 = indexToData(base + Eigen::Vector3i(0, 0, 1));
    int v3 = indexToData(base + Eigen::Vector3i(0, 1, 1));
    int v4 = indexToData(base + Eigen::Vector3i(1, 0, 1));
    Tetrahedron tetra1(v1, v2, v3, v4);
    // pushback
    tetra_data.push_back(tetra1);

    // tetrehedraon 2
    v2 = indexToData(base + Eigen::Vector3i(0, 1, 1));
    v3 = indexToData(base + Eigen::Vector3i(0, 1, 0));
    v4 = indexToData(base + Eigen::Vector3i(1, 1, 0));
    Tetrahedron tetra2(v1, v2, v3, v4);
    // pushback
    tetra_data.push_back(tetra2);

    // tetrehedraon 3
    v2 = indexToData(base + Eigen::Vector3i(1, 0, 0));
    v3 = indexToData(base + Eigen::Vector3i(1, 1, 0));
    v4 = indexToData(base + Eigen::Vector3i(1, 0, 1));
    Tetrahedron tetra3(v1, v2, v3, v4);
    // pushback
    tetra_data.push_back(tetra3);

    // tetrehedraon 4
    v2 = indexToData(base + Eigen::Vector3i(1, 1, 0));
    v3 = indexToData(base + Eigen::Vector3i(0, 1, 1));
    v4 = indexToData(base + Eigen::Vector3i(1, 0, 1));
    Tetrahedron tetra4(v1, v2, v3, v4);
    // pushback
    tetra_data.push_back(tetra4);

    // tetrehedraon 5
    v1 = indexToData(base + Eigen::Vector3i(1, 1, 1));
    v2 = indexToData(base + Eigen::Vector3i(1, 0, 1));
    v3 = indexToData(base + Eigen::Vector3i(1, 1, 0));
    v4 = indexToData(base + Eigen::Vector3i(0, 1, 1));
    Tetrahedron tetra5(v1, v2, v3, v4);
    // pushback
    tetra_data.push_back(tetra5);

    return;
};

