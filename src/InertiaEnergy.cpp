#include "InertiaEnergy.h"


// compute the elastic energy
double InertiaEnergy::Val(double nodeMass, double dt, Eigen::Vector3d& xt, Eigen::Vector3d& v, Eigen::Vector3d& x, Eigen::Vector3d& extForce, FEMParamters& param)
{
	Eigen::Vector3d x_minus_xt = x - (xt + dt * v + dt * dt / nodeMass * (nodeMass * param.gravity + extForce));
	return x_minus_xt.dot(x_minus_xt) * nodeMass / 2.0;
}

// compute the energy gradient wrt vertex's position.
void InertiaEnergy::Grad(std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_grad, double nodeMass, double dt, Eigen::Vector3d& xt, Eigen::Vector3d& v, Eigen::Vector3d& x, Eigen::Vector3d& extForce, int vertIndex, FEMParamters& param, int BC)
{
	Eigen::Vector3d x_minus_xt = x - (xt + dt * v + dt * dt / nodeMass * (nodeMass * param.gravity + extForce));
	Eigen::Vector3d gradVec = nodeMass * x_minus_xt;

	for (int dI = 0; dI < 3; dI++)
	{	
		if (BC != 1)
		{
			std::pair<int, double> pa = std::make_pair(vertIndex * 3 + dI, gradVec[dI]);
			grad_triplet[startIndex_grad + dI] = pa;
		}
		else
		{
			std::pair<int, double> pa = std::make_pair(vertIndex * 3 + dI, 0);
			grad_triplet[startIndex_grad + dI] = pa;
		}

	}

}

// the hessian is just an Identity matrix
void InertiaEnergy::Hess(std::vector<Eigen::Triplet<double>>& hessian_triplet, int& startIndex_hess, double nodeMass, int vertIndex, int BC)
{
	for (int dI = 0; dI < 3; dI++)
	{
		if (BC != 1)
		{
			Eigen::Triplet<double> pa = { vertIndex * 3 + dI, vertIndex * 3 + dI, nodeMass };
			hessian_triplet[startIndex_hess + dI] = pa;
		}
		else
		{
			Eigen::Triplet<double> pa = { vertIndex * 3 + dI, vertIndex * 3 + dI, 0 };
			hessian_triplet[startIndex_hess + dI] = pa;
		}
	}
}