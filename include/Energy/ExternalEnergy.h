#ifndef EXTERNALENERGY_H
#define EXTERNALENERGY_H

#include "utils.h"


namespace ExternalEnergy
{
	// compute the external energy: including applied force and gravity
	double Val(double nodeMass, double dt, Eigen::Vector3d& x, 
		FEMParamters& param, Eigen::Vector3d& extForce);

	// compute the energy gradient wrt vertex's position.
	void Grad(std::vector<std::pair<int, double>>& grad_triplet, 
		int& startIndex_grad, double nodeMass, double dt, 
		Eigen::Vector3d& x, FEMParamters& param, 
		Eigen::Vector3d& extForce, int vertIndex, int BC); // BC: boundary condition of the vertex

}



#endif