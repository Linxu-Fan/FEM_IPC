#include "InertiaEnergy.h"


// compute the elastic energy
double InertiaEnergy::Val(double nodeMass, double dt, Eigen::Vector3d& xt, Eigen::Vector3d& v, Eigen::Vector3d& x, Eigen::Vector3d& fe)
{
	Eigen::Vector3d x_minus_xt = x - (xt + dt * v + dt * dt / nodeMass * fe);
	return x_minus_xt.dot(x_minus_xt) * nodeMass / 2.0;
}

// compute the energy gradient wrt vertex's position.
Eigen::Vector3d InertiaEnergy::Grad(double nodeMass, double dt, Eigen::Vector3d& xt, Eigen::Vector3d& v, Eigen::Vector3d& x, Eigen::Vector3d& fe)
{
	Eigen::Vector3d x_minus_xt = x - (xt + dt * v + dt * dt / nodeMass * fe);
	return nodeMass * Eigen::Matrix3d::Identity() * x_minus_xt;
}

// the hessian is just an Identity matrix
Eigen::Matrix3d InertiaEnergy::Hess()
{
	return Eigen::Matrix3d::Identity();
}