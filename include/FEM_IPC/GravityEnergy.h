#ifndef GRAVITYENERGY_H
#define GRAVITYENERGY_H

#include "utils.h"
#include "materials.h"


namespace GravityEnergy
{
	// compute the inertia energy
	double Val(double nodeMass, double dt, Eigen::Vector3d& x, FEMParamters& param);

	// compute the energy gradient wrt vertex's position.
	std::vector<Eigen::Triplet<double>> Grad(double nodeMass, double dt, Eigen::Vector3d& x, FEMParamters& param, int vertIndex);

}



#endif