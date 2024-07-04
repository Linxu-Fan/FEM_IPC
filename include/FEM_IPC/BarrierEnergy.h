#ifndef BARRIERENERGY_H
#define BARRIERENERGY_H

#include "distance.h"
#include "mesh.h"

struct BarrierEnergyRes
{
	double Val = 0;
	std::vector<std::pair<int, double>> Grad;
	std::vector<Eigen::Triplet<double>> Hess;
};

namespace BarrierEnergy
{

	// compute the barrier energy
	double Val(bool pointTriangle, double dis2, Mesh& tetMesh, Eigen::Vector4i vtInd, double d_hat, double k_stiff, double dt);

	// compute the energy gradient 
	std::vector<std::pair<int, double>> Grad(bool pointTriangle, int type, double dis2, Mesh& tetMesh, Eigen::Vector4i vtInd, Eigen::Vector4i vtInd_BC, double d_hat, double k_stiff, double dt);

	// compute the energy hessian 
	std::vector<Eigen::Triplet<double>> Hess(bool pointTriangle, int type, double dis2, Mesh& tetMesh, Eigen::Vector4i vtInd, Eigen::Vector4i vtInd_BC, double d_hat, double k_stiff, double dt);

}

#endif