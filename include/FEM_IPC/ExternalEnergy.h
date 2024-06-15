#ifndef EXTERNALENERGY_H
#define EXTERNALENERGY_H

#include "utils.h"


namespace ExternalEnergy
{
	// compute the external energy: including applied force and gravity
	double Val(double nodeMass, double dt, Eigen::Vector3d& x, FEMParamters& param, Eigen::Vector3d& extForce);

	// compute the energy gradient wrt vertex's position.
	std::vector<std::pair<int, double>> Grad(double nodeMass, double dt, Eigen::Vector3d& x, FEMParamters& param, Eigen::Vector3d& extForce, int vertIndex);

}



#endif