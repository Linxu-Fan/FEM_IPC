#ifndef MLS_H
#define MLS_H

#include "utils.h"



class MLSPoints
{
public:
	double volume = 0; // volume of the MLS point
	Eigen::Vector3d pos = { 0,0,0 }; // position of the MLS point
	Eigen::Vector3d pos_Rest = { 0,0,0 }; // position of the MLS point in the rest configuration
	std::vector<int> index_node; // index of the nodes in the MLS point's support domain; Index is the index of the node in the original simulation mesh
	Eigen::Matrix3d F = Eigen::Matrix3d::Identity(); // deformation gradient of the MLS point
	std::vector<Eigen::Matrix<double, 9, 3>> dFdx; // the partial derivative of deformation gradient wrt the node position of the tetrahedral


	void MLS_approximation(const std::vector<Vector3d>& pos_node_Rest, 
		const std::vector<Vector3d>& pos_node, double& radius);

	// initialize the MLS point
	void init_MLS(Eigen::Vector3d& pos_Rest_, double volume_, 
		std::vector<int>& index_node_, std::string kernel, double& kernelSize);

	double computeWeight(const Eigen::Vector3d& x_s, const Eigen::Vector3d& x_i, double d_i);

	Eigen::Vector3d computeWeightDerivative(const Eigen::Vector3d& x_s, const Eigen::Vector3d& x_i, double d_i);



	// new MLS approximation
	Eigen::MatrixXd computeBasisGradient(int d);




	// quadratic MLS approximation
	// Function to construct the quadratic shifted basis for a point
	Eigen::VectorXd computeQuadraticShiftedBasis(const Eigen::VectorXd& p_i, const Eigen::VectorXd& p_s);

	// Function to compute the gradient of the quadratic shifted basis at p_s
	Eigen::MatrixXd computeBasisGradientAtPs(int d, int m);


};




#endif


