#include "volume.h"
template <typename T>
void SwapEnd(T &var)
{
    char *varArray = reinterpret_cast<char *>(&var);
    for (long i = 0; i < static_cast<long>(sizeof(var) / 2); i++)
        std::swap(varArray[sizeof(var) - 1 - i], varArray[i]);
}

Volume::Volume(){};
Volume::~Volume(){};

Volume::Volume(std::string volumefile, std::vector<Tetrahedron> & tetra_data, std::vector<MyVertex> & vertex_list)
{
    FILE *fp = fopen(volumefile.c_str(), "rb");
    if (fp)
    {   
        int vertex_amount, tetra_amount; 
        int waste;
        fread(&vertex_amount, sizeof(int), 1, fp);
        fread(&tetra_amount, sizeof(int), 1, fp);
        
        fread(&waste, sizeof(int), 1, fp);
        fread(&waste, sizeof(int), 1, fp);
        fread(&waste, sizeof(int), 1, fp);

        for(int i = 0; i < vertex_amount; i++){
            fread(&waste, sizeof(int), 1, fp);
            float posx, posy, posz;
            fread(&posx, sizeof(float), 1, fp);
            fread(&posy, sizeof(float), 1, fp);
            fread(&posz, sizeof(float), 1, fp);
            vertex_list.push_back(MyVertex(posx, posy, posz, 0.0f));
        }

        for(int i = 0; i < tetra_amount; i++){
            fread(&waste, sizeof(int), 1, fp);
            fread(&waste, sizeof(int), 1, fp);
            char waste_tet;
            fread(&waste_tet, sizeof(char), 1, fp);
            fread(&waste_tet, sizeof(char), 1, fp);
            fread(&waste_tet, sizeof(char), 1, fp);

            int idA, idB, idC, idD;
            fread(&idA, sizeof(float), 1, fp);
            fread(&idB, sizeof(float), 1, fp);
            fread(&idC, sizeof(float), 1, fp);
            fread(&idD, sizeof(float), 1, fp);

            idA -= 1;
            idB -= 1;
            idC -= 1;
            idD -= 1;
            tetra_data.push_back(Tetrahedron(idA, idB, idC, idD));
        }

        for(int i = 0; i < vertex_amount; i++){
            fread(&waste, sizeof(int), 1, fp);
            float density;
            fread(&density, sizeof(float), 1, fp);
            vertex_list[i].density = density;
        } 
    }
    fclose(fp);
    /*FILE *fp = fopen(volumefile.c_str(), "rb");
    if (fp)
    {    
        fread(this->size.data(),sizeof(int),3,fp);
        fread(this->size_physics.data(),sizeof(float),3,fp);
        this->bbox = AABB(0, 0, 0, size_physics.x(), size_physics.y(), size_physics.z());
        this->dx = size_physics.x() / (this->size.x()-1);
        int count = this->size.x() * this->size.y() * this->size.z();
        // this->raw_data.clear();
        float den = 0;
        float min_den = 9999;
        float max_den = -9999;
        for (int curind = 0; curind < count; curind++)
        {
            int z_idx = curind / (this->size.y() * this->size.x());
            int y_idx = (curind - z_idx * this->size.x() * this->size.y()) / this->size.x();
            int x_idx = curind - z_idx * this->size.x() * this->size.y() - y_idx * this->size.x();
            fread(&den, sizeof(float), 1, fp);
            min_den = std::min(min_den,den);
            max_den = std::max(max_den,den);
            vertex_list.push_back(MyVertex(this->dx * x_idx, this->dx * y_idx, this->dx * z_idx, den));
            if(x_idx<this->size.x() - 1 && y_idx<this->size.y() - 1 && z_idx<this->size.z() - 1){
                getVoxel(x_idx, y_idx, z_idx, tetra_data);
            }
        }
        //printf("density range [%.3f ,%.3f]\n",min_den,max_den);
    }
    fclose(fp);*/
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

