#include "InertiaEnergy.h"


// compute the elastic energy
double InertiaEnergy::Val(double nodeMass, double dt, Eigen::Vector3d& xt, Eigen::Vector3d& v, Eigen::Vector3d& x, Eigen::Vector3d& fe)
{
	Eigen::Vector3d x_minus_xt = x - (xt + dt * v + dt * dt / nodeMass * fe);
	return x_minus_xt.dot(x_minus_xt) * nodeMass / 2.0;
}

// compute the energy gradient wrt vertex's position.
std::vector<Eigen::Triplet<double>> InertiaEnergy::Grad(double nodeMass, double dt, Eigen::Vector3d& xt, Eigen::Vector3d& v, Eigen::Vector3d& x, Eigen::Vector3d& fe, int vertIndex)
{
	Eigen::Vector3d x_minus_xt = x - (xt + dt * v + dt * dt / nodeMass * fe);
	Eigen::Vector3d gradVec = nodeMass * Eigen::Matrix3d::Identity() * x_minus_xt;

	std::vector<Eigen::Triplet<double>> res;
	for (int dI = 0; dI < 3; dI++)
	{
		res.emplace_back(vertIndex * 3 + dI, vertIndex * 3 + dI, gradVec[dI]);
	}
	return res;
}

// the hessian is just an Identity matrix
std::vector<Eigen::Triplet<double>> InertiaEnergy::Hess(double nodeMass, int vertIndex)
{
	std::vector<Eigen::Triplet<double>> res;
	for (int dI = 0; dI < 3; dI++)
	{
		res.emplace_back(vertIndex * 3 + dI, vertIndex * 3 + dI, nodeMass);
	}
	return res;
}