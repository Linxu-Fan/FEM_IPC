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
	Eigen::Vector4i vtInd_BC = { -99, -99, -99, -99 }; // four points' boundary conditions

	// the actual energy value
	double Val = 0;
	std::vector<std::pair<int, double>> Grad;
	std::vector<Eigen::Triplet<double>> Hess;
};

namespace BarrierEnergy
{

	// compute the energy gradient and hessian of point-triangle contact
	void valGradAndHess_PT(BarrierEnergyRes& BEres, int type, double dis2, Mesh& tetMesh, double d_hat, double k_stiff, double dt);
	
	// compute the energy gradient and hessian of edge-edge contact
	void valGradAndHess_EE(BarrierEnergyRes& BEres, int type, double dis2, Mesh& tetMesh, double d_hat, double k_stiff, double dt);

	// store gradient and hessian value for point-triangle case
	void store_grad_hess_PT(BarrierEnergyRes& BEres, Vector6d& grad_, Matrix6d& hess_, Eigen::Vector2i& activePtsLocalInd);
	void store_grad_hess_PT(BarrierEnergyRes& BEres, Vector9d& grad_, Matrix9d& hess_, Eigen::Vector3i& activePtsLocalInd);
	void store_grad_hess_PT(BarrierEnergyRes& BEres, Vector12d& grad_, Matrix12d& hess_, Eigen::Vector4i& activePtsLocalInd);

	// store gradient and hessian value for edge-edge case
	void store_grad_hess_EE(BarrierEnergyRes& BEres, Vector12d& grad_final, Matrix12d& hess_final);

	// project edge-edge's gradient and hessian into full 12x1 vector and 12x12 matrix to faciliate mollifier mulitiplication
	void project_grad_to_full(Eigen::Vector2i& activePtsLocalInd, Vector6d& grad_, Matrix6d& hess_, Vector12d& grad_full, Matrix12d& hess__full);
	void project_grad_to_full(Eigen::Vector3i& activePtsLocalInd, Vector9d& grad_, Matrix9d& hess_, Vector12d& grad_full, Matrix12d& hess__full);


}

#endif