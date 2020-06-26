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
    this->computeGradients();
};
bool Volume::getRayStartEnd(Ray &ray, float &t_start, float &t_end)
{
    return bbox.rayIntersection(ray, t_start, t_end);
};

volumeData &Volume::indexToData(Eigen::Vector3i index)
{
    return this->raw_data.at(index.z() * size.y() * size.x() + index.y() * size.x() + index.x());
}
voxel Volume::getVoxel(Eigen::Vector3f position)
{
    position = position.cwiseMin(this->size_physics);
    Eigen::Vector3f index = (position - bbox.lb);
    index = index / this->dx;
    index=index.cwiseMax(Eigen::Vector3f::Zero());
    index=index.cwiseMin((this->size-Eigen::Vector3i(2,2,2)).cast<float>());

    Eigen::Vector3i lb(std::floor(index.x()), std::floor(index.y()), std::floor(index.z()));
    volumeData &c000 = indexToData(lb);

    volumeData &c100 = indexToData(lb + Eigen::Vector3i(1, 0, 0));
    volumeData &c010 = indexToData(lb + Eigen::Vector3i(0, 1, 0));
    volumeData &c110 = indexToData(lb + Eigen::Vector3i(1, 1, 0));
    volumeData &c001 = indexToData(lb + Eigen::Vector3i(0, 0, 1));
    volumeData &c101 = indexToData(lb + Eigen::Vector3i(1, 0, 1));
    volumeData &c011 = indexToData(lb + Eigen::Vector3i(0, 1, 1));
    volumeData &c111 = indexToData(lb + Eigen::Vector3i(1, 1, 1));
    return voxel(c000, c100, c010, c110, c001, c101, c011, c111);
};
void Volume::computeGradients()
{
    this->grad_maxnorm = std::numeric_limits<float>::min();
    int count=this->raw_data.size();
    for (int i = 0; i < count; i++)
    {
		int z_idx = i / (this->size.y() * this->size.x());
		int y_idx = (i - z_idx * this->size.x() * this->size.y()) / this->size.x();
		int x_idx = i - z_idx * this->size.x() * this->size.y() - y_idx * this->size.x();
			

        float x_grad = 0, y_grad = 0, z_grad = 0;
        // handle boudary conditions, set boundary's graident as its neighbours' gradient
        if (z_idx == 0) { z_idx = 1; }
        if (z_idx == this->size.z() - 1) { z_idx -= 1; }
        if (y_idx == 0) { y_idx = 1; }
        if (y_idx == this->size.y() - 1) { y_idx -= 1; }
        if (x_idx == 0) { x_idx = 1; }
        if (x_idx == this->size.x() - 1) { x_idx -= 1; }

        // Compute your gradient here
        x_grad = (indexToData(Eigen::Vector3i(x_idx + 1, y_idx, z_idx)).density - indexToData(Eigen::Vector3i(x_idx - 1, y_idx, z_idx)).density) /
            (indexToData(Eigen::Vector3i(x_idx + 1, y_idx, z_idx)).position.x() - indexToData(Eigen::Vector3i(x_idx - 1, y_idx, z_idx)).position.x())/2;
        y_grad = (indexToData(Eigen::Vector3i(x_idx, y_idx + 1, z_idx)).density - indexToData(Eigen::Vector3i(x_idx, y_idx - 1, z_idx)).density) /
            (indexToData(Eigen::Vector3i(x_idx, y_idx + 1, z_idx)).position.y() - indexToData(Eigen::Vector3i(x_idx, y_idx - 1, z_idx)).position.y())/2;
        z_grad = (indexToData(Eigen::Vector3i(x_idx, y_idx, z_idx + 1)).density - indexToData(Eigen::Vector3i(x_idx, y_idx, z_idx - 1)).density) /
            (indexToData(Eigen::Vector3i(x_idx, y_idx, z_idx + 1)).position.z() - indexToData(Eigen::Vector3i(x_idx, y_idx, z_idx - 1)).position.z())/2;

		Eigen::Vector3f grad = Eigen::Vector3f(x_grad, y_grad, z_grad);
		
		this->grad_maxnorm = std::max(this->grad_maxnorm,grad.norm());
        this->raw_data[i].gradient =grad;
    }
    printf("max_gradient norm %f \n",grad_maxnorm);
};
