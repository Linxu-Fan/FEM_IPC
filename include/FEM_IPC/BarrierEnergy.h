#ifndef BARRIERENERGY_H
#define BARRIERENERGY_H

#include "utils.h"

namespace BarrierEnergy
{

	// compute the barrier energy
	double Val(Material& mat, std::string model, Eigen::Matrix3d& F, double dt, double vol);

	// compute the energy gradient dPHI/dF (PK1 stress). Return a vectorized gradient which is a 9x1 matrix; tetVertInd_BC: boundary condition of each vertex
	std::vector<std::pair<int, double>> Grad(Material& mat, std::string model, Eigen::Matrix3d& F, double dt, double vol, Eigen::Matrix<double, 9, 12>& dFdx, Eigen::Vector4i& tetVertInd, Eigen::Vector4i& tetVertInd_BC); // timestep, and volume of the element

	// compute the energy hessian dPHI2/d2F. Return a vectorized gradient which is a 9x9 matrix; tetVertInd_BC: boundary condition of each vertex
	std::vector<Eigen::Triplet<double>> Hess(Material& mat, std::string model, Eigen::Matrix3d& F, double dt, double vol, Eigen::Matrix<double, 9, 12>& dFdx, Eigen::Vector4i& tetVertInd, Eigen::Vector4i& tetVertInd_BC); // timestep, and volume of the element

}

#endif