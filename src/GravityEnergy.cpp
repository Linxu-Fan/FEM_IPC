#include "GravityEnergy.h"


// compute the elastic energy
double GravityEnergy::Val(double nodeMass, double dt, Eigen::Vector3d& x, FEMParamters& param)
{	
	return - dt * dt * nodeMass * x.dot(param.gravity);
}

// compute the energy gradient wrt vertex's position.
std::vector<Eigen::Triplet<double>> GravityEnergy::Grad(double nodeMass, double dt, Eigen::Vector3d& x, FEMParamters& param, int vertIndex)
{
	Eigen::Vector3d gravityForceEng = -dt * dt * nodeMass * param.gravity;

	std::vector<Eigen::Triplet<double>> res;
	for (int dI = 0; dI < 3; dI++)
	{
		res.emplace_back(vertIndex * 3 + dI, vertIndex * 3 + dI, gravityForceEng[dI]);
	}
	return res;
}
