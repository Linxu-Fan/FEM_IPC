#include "ExternalEnergy.h"


// compute the external energy: including applied force and gravity
double ExternalEnergy::Val(double nodeMass, double dt, Eigen::Vector3d& x, FEMParamters& param, Eigen::Vector3d& extForce)
{	
	return - dt * dt * x.dot(nodeMass * param.gravity + extForce);
}

// compute the energy gradient wrt vertex's position.
std::vector<std::pair<int, double>> ExternalEnergy::Grad(double nodeMass, double dt, Eigen::Vector3d& x, FEMParamters& param, Eigen::Vector3d& extForce, int vertIndex, int BC)
{
	Eigen::Vector3d gravityForceEng = -dt * dt * (nodeMass * param.gravity + extForce);

	std::vector<std::pair<int, double>> res;
	for (int dI = 0; dI < 3; dI++)
	{
		if (BC != 1)
		{
			res.emplace_back(vertIndex * 3 + dI, gravityForceEng[dI]);
		}
		else
		{
			res.emplace_back(vertIndex * 3 + dI, 0);
		}
		
	}
	return res;
}
