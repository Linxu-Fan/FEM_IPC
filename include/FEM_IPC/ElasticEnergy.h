#ifndef ELASTICENERGY_H
#define ELASTICENERGY_H

#include "utils.h"

namespace ElasticEnergy
{
	// compute df/dx: the partial derivative of deformation gradient wrt the node position per tetrahedral
	Eigen::Matrix<double, 9, 12> dF_wrt_dx(Eigen::Matrix3d& DmInv);

	// compute the elastic energy
	double Val(Material& mat, std::string model, Eigen::Matrix3d& F, double dt, double vol);

	// compute the energy gradient dPHI/dF (PK1 stress). Return a vectorized gradient which is a 9x1 matrix
	std::vector<std::pair<int, double>> Grad(Material& mat, std::string model, Eigen::Matrix3d& F, double dt, double vol, Eigen::Matrix<double, 9, 12>& dFdx, Eigen::Vector4i& tetVertInd); // timestep, and volume of the element

	// compute the energy hessian dPHI2/d2F. Return a vectorized gradient which is a 9x9 matrix
	std::vector<Eigen::Triplet<double>> Hess(Material& mat, std::string model, Eigen::Matrix3d& F, double dt, double vol, Eigen::Matrix<double, 9, 12>& dFdx, Eigen::Vector4i& tetVertInd); // timestep, and volume of the element

}

#endif