#ifndef BARRIERENERGY_H
#define BARRIERENERGY_H

#include "distance.h"
#include "mesh.h"

struct BarrierEnergyRes
{
	// store the culled constraint set
	bool pointTriangle = true; // true: this is a point-triangle contact; false: edge-edge contact
	Eigen::Vector2i PT_Index = {-99, -99}; // 1st int: point index; 2nd int: triangle index
	Eigen::Vector4i PP_Index = {-99, -99, -99, -99 }; // if pointTriangle is true: 1st int: point index; 2nd 3rd & 4th int: triangle vertices indices; else: four points' indices in two edges

	// the actual energy value
	double Val = 0;
	std::vector<std::pair<int, double>> Grad;
	std::vector<Eigen::Triplet<double>> Hess;
};

namespace BarrierEnergy
{

	// compute the barrier energy
	double Val(bool pointTriangle, double dis2, Mesh& tetMesh, Eigen::Vector4i& vtInd, double d_hat, double k_stiff, double dt);

	//// compute the energy gradient 
	//std::vector<std::pair<int, double>> Grad(bool pointTriangle, int type, double dis2, Mesh& tetMesh, Eigen::Vector4i& vtInd, Eigen::Vector4i& vtInd_BC, double d_hat, double k_stiff, double dt);

	//// compute the energy hessian 
	//std::vector<Eigen::Triplet<double>> Hess(bool pointTriangle, int type, double dis2, Mesh& tetMesh, Eigen::Vector4i& vtInd, Eigen::Vector4i& vtInd_BC, double d_hat, double k_stiff, double dt);

	// compute the energy gradient and hessian of point-triangle contact
	std::pair<std::vector<std::pair<int, double>>, std::vector<Eigen::Triplet<double>>> gradAndHess_PT(int type, double dis2, 
		Mesh& tetMesh, Eigen::Vector4i& vtInd, Eigen::Vector4i& vtInd_BC, double d_hat, double k_stiff, double dt);
	
	// compute the energy gradient and hessian of edge-edge contact
	std::pair<std::vector<std::pair<int, double>>, std::vector<Eigen::Triplet<double>>> gradAndHess_EE(int type, double dis2,
		Mesh& tetMesh, Eigen::Vector4i& vtInd, Eigen::Vector4i& vtInd_BC, double d_hat, double k_stiff, double dt);

	// store gradient and hessian value
	void store_grad_hess(std::vector<std::pair<int, double>>& res_grad, std::vector<Eigen::Triplet<double>>& res_hess, 
		Vector6d& grad_, Matrix6d& hess_,
		std::vector<int>& activePts, std::vector<int>& activePtsBC);

	void store_grad_hess(std::vector<std::pair<int, double>>& res_grad, std::vector<Eigen::Triplet<double>>& res_hess,
		Vector9d& grad_, Matrix9d& hess_,
		std::vector<int>& activePts, std::vector<int>& activePtsBC);

	void store_grad_hess(std::vector<std::pair<int, double>>& res_grad, std::vector<Eigen::Triplet<double>>& res_hess,
		Vector12d& grad_, Matrix12d& hess_,
		std::vector<int>& activePts, std::vector<int>& activePtsBC);

}

#endif