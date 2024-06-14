#ifndef INERTIAENERGY_H
#define INERTIAENERGY_H

#include "utils.h"
#include "materials.h"


namespace InertiaEnergy
{
	// compute the inertia energy
	double Val(double nodeMass, double dt, Eigen::Vector3d& xt, Eigen::Vector3d& v, Eigen::Vector3d& x, Eigen::Vector3d& fe);

	// compute the energy gradient wrt vertex's position.
	Eigen::Vector3d Grad(double nodeMass, double dt, Eigen::Vector3d& xt, Eigen::Vector3d& v, Eigen::Vector3d& x, Eigen::Vector3d& fe);

	// the hessian is just an Identity matrix
	Eigen::Matrix3d Hess();

}



#endif