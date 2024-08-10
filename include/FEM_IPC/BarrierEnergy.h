#ifndef BARRIERENERGY_H
#define BARRIERENERGY_H

#include "distance.h"
#include "mesh.h"


namespace BarrierEnergy
{
	// compute the energy of point-triangle contact
	double val_PT(double contactArea, double dis2, double d_hat, double k_stiff, double dt);

	// compute the energy of edge-edge contact
	double val_EE(double contactArea, double dis2, Mesh& tetMesh, Eigen::Vector4i& ptIndices, double d_hat, double k_stiff, double dt);


	// compute the energy gradient and hessian of point-triangle contact
	void gradAndHess_PT(BarrierEnergyRes& GH, Eigen::Vector4i& ptIndices, int type, double dis2, Mesh& tetMesh, double d_hat, double k_stiff, double dt);

	// compute the energy gradient and hessian of edge-edge contact
	void gradAndHess_EE(BarrierEnergyRes& GH, Eigen::Vector4i& ptIndices, int type, double dis2, Mesh& tetMesh, double d_hat, double k_stiff, double dt);


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


namespace Ground
{

	double val(double coor_z, double d_hat, double distributedArea, double k_stiff, double dt);

	void gradAndHess(BarrierEnergyRes& BEres, int index_i, double coor_z, double d_hat, double distributedArea, double k_stiff, double dt);


	void valGradAndHess(BarrierEnergyRes& BEres, int index_i, double dis, double d_hat, double distributedArea, double k_stiff, double dt);

}


#endif