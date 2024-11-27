#ifndef INERTIAENERGY_H
#define INERTIAENERGY_H

#include "utils.h"


namespace InertiaEnergy
{
	// compute the inertia energy; vertIndex: the index of the vertex
	double Val(double nodeMass, double dt, Eigen::Vector3d& xt, Eigen::Vector3d& v, 
		Eigen::Vector3d& x, Eigen::Vector3d& extForce, FEMParamters& param, int& vertIndex, 
		std::vector<boundaryCondition>& boundaryCondition_node, int timestep);

	// compute the energy gradient wrt vertex's position.
	void Grad(std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_grad,
		double nodeMass, double dt, Eigen::Vector3d& xt, Eigen::Vector3d& v, Eigen::Vector3d& x,
		Eigen::Vector3d& extForce, int& vertIndex, FEMParamters& param, 
		std::vector<boundaryCondition>& boundaryCondition_node, int timestep);

	// the hessian is just an Identity matrix
	void Hess(std::vector<Eigen::Triplet<double>>& hessian_triplet, int& startIndex_hess, 
		double nodeMass, int& vertIndex, std::vector<boundaryCondition>& boundaryCondition_node, 
		int timestep, FEMParamters& param);

}



#endif