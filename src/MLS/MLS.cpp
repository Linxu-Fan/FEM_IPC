#include "mesh.h"

Eigen::Vector4d MLSPoints::computeP(const Eigen::Vector3d& x)
{
    Eigen::Vector4d P;
    P << 1.0, x(0), x(1), x(2);
    return P;
}

double MLSPoints::computeWeight(const Eigen::Vector3d& x_s, const Eigen::Vector3d& x_i, double d_i)
{
    double r = (x_s - x_i).norm();
    return exp(-(r / d_i) * (r / d_i));
}

Eigen::Vector3d MLSPoints::computeWeightDerivative(const Eigen::Vector3d& x_s, const Eigen::Vector3d& x_i, double d_i)
{
    Eigen::Vector3d diff = x_s - x_i;
    double r = diff.norm();
    double w = computeWeight(x_s, x_i, d_i);
    double w_prime = (-2.0 * w) / (d_i * d_i);
    return w_prime * (diff / r);
}

void MLSPoints::init_MLS(Eigen::Vector3d& pos_Rest_, double volume_, std::vector<int>& index_node_, std::vector<Eigen::Vector3d>& pos_node_Rest, std::string kernel, double kernelSize)
{
    /*pos_Rest = pos_Rest_;
    pos = pos_Rest_;
	volume = volume_;
	index_node = index_node_;
    centroid_ref = ComputeWeightedCentroid(pos_node_Rest);
	cal_MLS_weights(pos_node_Rest, kernel, kernelSize);*/

}

void MLSPoints::update_MLS(std::vector<Eigen::Vector3d>& pos_node)
{
    //pos = cal_MLS_point_pos(pos_node);
    //F = ComputeDeformationGradient(pos_node);
}



