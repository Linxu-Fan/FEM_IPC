#ifndef ELASTICENERGY_H
#define ELASTICENERGY_H

#include "utils.h"

namespace ElasticEnergy
{
	// compute df/dx: the partial derivative of deformation gradient wrt the node position per tetrahedral
	Eigen::Matrix<double, 9, 12> dF_wrt_dx(Eigen::Matrix3d& DmInv);

	// compute the elastic energy
	double Val(Material& mat, std::string model, Eigen::Matrix3d& F, double dt, double vol);

	Eigen::Matrix3d calPK1(Material& mat, std::string model, Eigen::Matrix3d& F); // compute the first Piola-Kirchhoff stress tensor

	// compute the energy gradient dPHI/dF (PK1 stress). Return a vectorized gradient which is a 9x1 matrix; tetVertInd_BC: boundary condition of each vertex
	// startIndex is the start index of grad_triplet
	void Grad(std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_grad, 
		Material& mat, std::string model, Eigen::Matrix3d& F, double dt, double vol, 
		Eigen::Matrix<double, 9, 12>& dFdx, Eigen::Vector4i& tetVertInd, Eigen::Vector4i& tetVertInd_BC); // timestep, and volume of the element

	Eigen::Matrix<double, 9, 9> calPK1_wrt_F(Material& mat, std::string model, Eigen::Matrix3d& F); // compute the derivative of first Piola-Kirchhoff stress tensor wrt F

	// compute the energy hessian dPHI2/d2F. Return a vectorized gradient which is a 9x9 matrix; tetVertInd_BC: boundary condition of each vertex
	// startIndex is the start index of hessian_triplet
	void Hess(std::vector<Eigen::Triplet<double>>& hessian_triplet, int& startIndex_hess, 
		Material& mat, std::string model, Eigen::Matrix3d& F, double dt, double vol, 
		Eigen::Matrix<double, 9, 12>& dFdx, Eigen::Vector4i& tetVertInd, Eigen::Vector4i& tetVertInd_BC); // timestep, and volume of the element

}

#endif