#ifndef MLS_H
#define MLS_H

#include "utils.h"



class MLSPoints
{
public:
	double volume = 0; // volume of the MLS point
	Eigen::Vector3d pos = { 0,0,0 }; // position of the MLS point
	Eigen::Vector3d pos_Rest = { 0,0,0 }; // position of the MLS point in the rest configuration
	Eigen::Matrix3d F = Eigen::Matrix3d::Identity(); // deformation gradient of the MLS point
	std::vector<int> index_node; // index of the nodes in the MLS point's support domain; Index is the index of the node in the original simulation mesh
	std::vector<double> weight; // weight of each node in the MLS point's support domain
	std::vector<Eigen::Matrix<double, 9, 3>> dFdx; // the partial derivative of deformation gradient wrt the node position of the tetrahedral
	Eigen::Vector3d centroid_ref = {0,0,0}; // weighted centroid of the MLS point in the reference configuration

	// initialize the MLS point
	void init_MLS(Eigen::Vector3d& pos_Rest_, double volume_, std::vector<int>& index_node_, std::vector<Eigen::Vector3d>& pos_node_Rest, std::string kernel, double kernelSize);

	// calculate the weight of each point
	void cal_MLS_weights(std::vector<Eigen::Vector3d>& pos_node_Rest, std::string kernel, double kernelSize);

	Eigen::Vector3d cal_MLS_point_pos(const std::vector<Eigen::Vector3d>& pos_node);

	void update_MLS(std::vector<Eigen::Vector3d>& pos_node); // update the MLS point's information after advection

	double ComputeGaussianWeight(double distance_squared, double kernelSize);

	Eigen::Vector3d ComputeWeightedCentroid(const std::vector<Eigen::Vector3d>& pos_node);

	Eigen::Matrix3d ComputeDeformationGradient(const std::vector<Eigen::Vector3d>& pos_node);

};




#endif


