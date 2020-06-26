#pragma once
#include "volumeData.h"
#include "Eigen/Dense"
//ref https://en.wikipedia.org/wiki/Trilinear_interpolation
class Interpolator
{
public:
    virtual volumeData interpolate(Eigen::Vector3f &p, const voxel &voxel) = 0;
};
class NearestNeighbourInterpolator : public Interpolator
{
    inline volumeData interpolate(Eigen::Vector3f &p, const voxel &voxel)
    {
        //Here your NNInterpolator code
        volumeData resultdata;
        resultdata.position = p;
        float xmid = (voxel.c000.position.x() + voxel.c100.position.x()) / 2.0;
        float ymid = (voxel.c000.position.y() + voxel.c010.position.y()) / 2.0;
        float zmid = (voxel.c000.position.z() + voxel.c001.position.z()) / 2.0;

        if (p.x() >= xmid && p.y() >= ymid && p.z() >= zmid) {
            resultdata.density = voxel.c111.density;
            resultdata.gradient = voxel.c111.gradient;
            return resultdata;
        }
        if (p.x() <= xmid && p.y() >= ymid && p.z() >= zmid) {
            resultdata.density = voxel.c011.density;
            resultdata.gradient = voxel.c011.gradient;
            return resultdata;
        }
        if (p.x() >= xmid && p.y() <= ymid && p.z() >= zmid) {
            resultdata.density = voxel.c101.density;
            resultdata.gradient = voxel.c101.gradient;
            return resultdata;
        }
        if (p.x() >= xmid && p.y() >= ymid && p.z() <= zmid) {
            resultdata.density = voxel.c110.density;
            resultdata.gradient = voxel.c110.gradient;
            return resultdata;
        }
        if (p.x() <= xmid && p.y() <= ymid && p.z() >= zmid) {
            resultdata.density = voxel.c001.density;
            resultdata.gradient = voxel.c001.gradient;
            return resultdata;
        }
        if (p.x() >= xmid && p.y() <= ymid && p.z() <= zmid) {
            resultdata.density = voxel.c100.density;
            resultdata.gradient = voxel.c100.gradient;
            return resultdata;
        }
        if (p.x() <= xmid && p.y() >= ymid && p.z() <= zmid) {
            resultdata.density = voxel.c010.density;
            resultdata.gradient = voxel.c010.gradient;
            return resultdata;
        }
        if (p.x() <= xmid && p.y() <= ymid && p.z() <= zmid) {
            resultdata.density = voxel.c000.density;
            resultdata.gradient = voxel.c000.gradient;
            return resultdata;
        }
    };
};
class TrilinearInterpolator : public Interpolator
{

    inline volumeData interpolate(Eigen::Vector3f &p, const voxel &voxel)
    {
        
        //Here your TrilinearInterpolator code
        volumeData resultdata;
        float xd = (p.x() - voxel.c000.position.x()) / (voxel.c100.position.x() - voxel.c000.position.x());
        float yd = (p.y() - voxel.c000.position.y()) / (voxel.c010.position.y() - voxel.c000.position.y());
        float zd = (p.z() - voxel.c000.position.z()) / (voxel.c001.position.z() - voxel.c000.position.z());

        float d00 = voxel.c000.density * (1 - xd) + voxel.c100.density * xd;
        float d01 = voxel.c001.density * (1 - xd) + voxel.c101.density * xd;
        float d10 = voxel.c010.density * (1 - xd) + voxel.c110.density * xd;
        float d11 = voxel.c011.density * (1 - xd) + voxel.c111.density * xd;
        Eigen::Vector3f g00 = voxel.c000.gradient * (1 - xd) + voxel.c100.gradient * xd;
        Eigen::Vector3f g01 = voxel.c001.gradient * (1 - xd) + voxel.c101.gradient * xd;
        Eigen::Vector3f g10 = voxel.c010.gradient * (1 - xd) + voxel.c110.gradient * xd;
        Eigen::Vector3f g11 = voxel.c011.gradient * (1 - xd) + voxel.c111.gradient * xd;

        float d0 = d00 * (1 - yd) + d10 * yd;
        float d1 = d01 * (1 - yd) + d11 * yd;
        Eigen::Vector3f g0 = g00 * (1 - yd) + g10 * yd;
        Eigen::Vector3f g1 = g01 * (1 - yd) + g11 * yd;

        resultdata.position = p;
        resultdata.density = d0 * (1 - zd) + d1 * zd;
        resultdata.gradient = g0 * (1 - zd) + g1 * zd;
        return resultdata;

    
    };

};