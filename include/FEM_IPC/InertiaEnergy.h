#ifndef INERTIAENERGY_H
#define INERTIAENERGY_H

#include "utils.h"


namespace InertiaEnergy
{
	// compute the inertia energy; vertIndex: the index of the vertex
	double Val(double nodeMass, double dt, Eigen::Vector3d& xt, Eigen::Vector3d& v, Eigen::Vector3d& x, Eigen::Vector3d& fe);

	// compute the energy gradient wrt vertex's position.
	std::vector<std::pair<int, double>> Grad(double nodeMass, double dt, Eigen::Vector3d& xt, Eigen::Vector3d& v, Eigen::Vector3d& x, Eigen::Vector3d& fe, int vertIndex);

	// the hessian is just an Identity matrix
	std::vector<Eigen::Triplet<double>> Hess(double nodeMass, int vertIndex);

}



#endif