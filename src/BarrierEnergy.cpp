#include "BarrierEnergy.h"



// compute the elastic energy
double BarrierEnergy::Val(Material& mat, std::string model, Eigen::Matrix3d& F, double dt, double vol)
{
	double energy = 0;
	return dt * dt * energy * vol;
}

// compute the energy gradient dPHI/dF (PK1 stress). Return a vectorized gradient which is a 9x1 matrix
std::vector<std::pair<int, double>> BarrierEnergy::Grad(Material& mat, std::string model, Eigen::Matrix3d& F, double dt, double vol, Eigen::Matrix<double, 9, 12>& dFdx, Eigen::Vector4i& tetVertInd, Eigen::Vector4i& tetVertInd_BC)
{
	Vector9d grad_tmp;
	std::vector<std::pair<int, double>> res;
	return  res;
}

// compute the energy hessian dPHI2/d2F. Return a vectorized gradient which is a 9x9 matrix
std::vector<Eigen::Triplet<double>> BarrierEnergy::Hess(Material& mat, std::string model, Eigen::Matrix3d& F, double dt, double vol, Eigen::Matrix<double, 9, 12>& dFdx, Eigen::Vector4i& tetVertInd, Eigen::Vector4i& tetVertInd_BC)
{

	std::vector<Eigen::Triplet<double>> res;
	return res;

}
