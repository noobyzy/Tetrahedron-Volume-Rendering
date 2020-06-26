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

Volume::Volume(std::string volumefile)
{
    FILE *fp = fopen(volumefile.c_str(), "rb");
    if (fp)
    {    
            fread(this->size.data(),sizeof(int),3,fp);
            fread(this->size_physics.data(),sizeof(float),3,fp);
    this->bbox = AABB(0, 0, 0, size_physics.x(), size_physics.y(), size_physics.z());
    this->dx = size_physics.x() / (this->size.x()-1);
    int count = this->size.x() * this->size.y() * this->size.z();
        this->raw_data.clear();
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
            this->raw_data.push_back(MyVertex(this->dx * x_idx, this->dx * y_idx, this->dx * z_idx, den));
        }
        printf("density range [%.3f ,%.3f]\n",min_den,max_den);
    }
    fclose(fp);
};

bool Volume::getRayStartEnd(Ray &ray, float &t_start, float &t_end)
{
    return bbox.rayIntersection(ray, t_start, t_end);
};

MyVertex& Volume::indexToData(Eigen::Vector3i index)
{
    return raw_data[index.z() * size.y() * size.x() + index.y() * size.x() + index.x()];
};