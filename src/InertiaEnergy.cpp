#include "InertiaEnergy.h"


// compute the elastic energy
double InertiaEnergy::Val(double nodeMass, double dt, Eigen::Vector3d& xt, Eigen::Vector3d& v, Eigen::Vector3d& x, 
	Eigen::Vector3d& extForce, FEMParamters& param, int& vertIndex, std::vector<boundaryCondition>& boundaryCondition_node, int timestep)
{
	double energy = 0;
	Eigen::Vector3d x_minus_xt = x - (xt + dt * v + dt * dt / nodeMass * (nodeMass * param.gravity + extForce));
	energy += x_minus_xt.dot(x_minus_xt) * nodeMass / 2.0;


	if (boundaryCondition_node[vertIndex].type == 1  && timestep >= boundaryCondition_node[vertIndex].appliedTime[0] && timestep <= boundaryCondition_node[vertIndex].appliedTime[1])
	{
		energy += param.IPC_B3Stiffness * nodeMass / 2.0 * (x - boundaryCondition_node[vertIndex].location[timestep]).dot(x - boundaryCondition_node[vertIndex].location[timestep]);
	}


	return energy;
}

// compute the energy gradient wrt vertex's position.
void InertiaEnergy::Grad(std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_grad, 
	double nodeMass, double dt, Eigen::Vector3d& xt, Eigen::Vector3d& v, Eigen::Vector3d& x, 
	Eigen::Vector3d& extForce, int& vertIndex, FEMParamters& param, std::vector<boundaryCondition>& boundaryCondition_node, int timestep)
{
	Eigen::Vector3d x_minus_xt = x - (xt + dt * v + dt * dt / nodeMass * (nodeMass * param.gravity + extForce));
	Eigen::Vector3d gradVec = nodeMass * x_minus_xt;

	if (boundaryCondition_node[vertIndex].type == 1 && timestep >= boundaryCondition_node[vertIndex].appliedTime[0] && timestep <= boundaryCondition_node[vertIndex].appliedTime[1])
	{
		gradVec += param.IPC_B3Stiffness * nodeMass * (x - boundaryCondition_node[vertIndex].location[timestep]);
	}

	for (int dI = 0; dI < 3; dI++)
	{	
		grad_triplet[startIndex_grad + dI] = { vertIndex * 3 + dI, gradVec[dI] };
	}

}

// the hessian is just an Identity matrix
void InertiaEnergy::Hess(std::vector<Eigen::Triplet<double>>& hessian_triplet, int& startIndex_hess, 
	double nodeMass, int& vertIndex, std::vector<boundaryCondition>& boundaryCondition_node, int timestep, FEMParamters& param)
{
	double hessVal = nodeMass;
	if (boundaryCondition_node[vertIndex].type == 1 && timestep >= boundaryCondition_node[vertIndex].appliedTime[0] && timestep <= boundaryCondition_node[vertIndex].appliedTime[1])
	{
		hessVal += param.IPC_B3Stiffness * nodeMass;
	}

	for (int dI = 0; dI < 3; dI++)
	{
		hessian_triplet[startIndex_hess + dI] = { vertIndex * 3 + dI, vertIndex * 3 + dI, hessVal };
	}
}