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

	// project edge-edge's gradient and hessian into full 12x1 vector and 12x12 matrix to faciliate mollifier mulitiplication
	void project_grad_to_full(Eigen::Vector2i& activePtsLocalInd, Vector6d& grad_, Matrix6d& hess_, Vector12d& grad_full, Matrix12d& hess__full);
	void project_grad_to_full(Eigen::Vector3i& activePtsLocalInd, Vector9d& grad_, Matrix9d& hess_, Vector12d& grad_full, Matrix12d& hess__full);




	double compute_b(double d, double dHat);

	double compute_g_b(double d, double dHat);

	double compute_H_b(double d, double dHat);


}


namespace Ground
{

	double val(double coor_z, double d_hat, double distributedArea, double k_stiff, double dt);

	void gradAndHess(BarrierEnergyRes& BEres, int index_i, double coor_z, double d_hat, double distributedArea, double k_stiff, double dt);


}


#endif