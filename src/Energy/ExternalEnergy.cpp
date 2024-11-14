#include "ExternalEnergy.h"


// compute the external energy: including applied force and gravity
double ExternalEnergy::Val(double nodeMass, double dt, Eigen::Vector3d& x, FEMParamters& param, Eigen::Vector3d& extForce)
{	
	return - dt * dt * x.dot(nodeMass * param.gravity + extForce);
}

// compute the energy gradient wrt vertex's position.
void ExternalEnergy::Grad(std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_grad, double nodeMass, double dt, Eigen::Vector3d& x, FEMParamters& param, Eigen::Vector3d& extForce, int vertIndex, int BC)
{
	Eigen::Vector3d gravityForceEng = -dt * dt * (nodeMass * param.gravity + extForce);

	for (int dI = 0; dI < 3; dI++)
	{
		grad_triplet[startIndex_grad + dI] = { vertIndex * 3 + dI, gravityForceEng[dI] };
	}
}
