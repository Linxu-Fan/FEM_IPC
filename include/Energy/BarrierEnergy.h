#ifndef BARRIERENERGY_H
#define BARRIERENERGY_H

#include "distance.h"
#include "mesh.h"


namespace BarrierEnergy
{
	// compute the energy of point-triangle contact
	double val_PT(double& contactArea, double& dis2, FEMParamters& parameters);

	// compute the energy of edge-edge contact
	double val_EE(double& contactArea, double& dis2, Mesh& tetSimMesh,
		Eigen::Vector4i& ptIndices, FEMParamters& parameters);

	// compute the energy gradient and hessian of point-triangle contact
	void gradAndHess_PT(std::vector<Eigen::Triplet<double>>& hessian_triplet, 
		std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_hess, 
		int& startIndex_grad,  
		Eigen::Vector4i& ptIndices, int& type, double& dis2, Mesh& tetSimMesh,
		FEMParamters& parameters, bool ABD = false);

	void cal_PT_PP_gradAndHess(std::vector<Eigen::Triplet<double>>& hessian_triplet,
		std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_hess,
		Eigen::Vector2i& activePtsLocalInd,
		int& startIndex_grad, Mesh& tetSimMesh,
		FEMParamters& parameters, double& contactArea, double& g_bd, double& h_bd, bool ABD = false);

	void cal_PT_PE_gradAndHess(std::vector<Eigen::Triplet<double>>& hessian_triplet,
		std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_hess,
		Eigen::Vector3i& activePtsLocalInd,
		int& startIndex_grad, Mesh& tetSimMesh,
		FEMParamters& parameters, double& contactArea, double& g_bd, double& h_bd, bool ABD = false);

	void cal_PT_PT_gradAndHess(std::vector<Eigen::Triplet<double>>& hessian_triplet,
		std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_hess,		 
		Eigen::Vector4i& activePtsLocalInd,
		int& startIndex_grad, Mesh& tetSimMesh,
		FEMParamters& parameters, double& contactArea, double& g_bd, double& h_bd, bool ABD = false);



	// compute the energy gradient and hessian of edge-edge contact
	void gradAndHess_EE(std::vector<Eigen::Triplet<double>>& hessian_triplet, 
		std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_hess, 
		int& startIndex_grad,
		Eigen::Vector4i& ptIndices, int& type, double& dis2, Mesh& tetSimMesh,
		FEMParamters& parameters, bool ABD = false);

	//void cal_EE_PP_gradAndHess(std::vector<Eigen::Triplet<double>>& hessian_triplet,
	//	std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_hess,
	//	int& startIndex_grad, std::vector<boundaryCondition>& boundaryCondition_node,
	//	Eigen::Vector2i& activePtsLocalInd, int& type, double& dis2, Mesh& tetSimMesh,
	//	FEMParamters& parameters,double& g_bd, double& h_bd, double& contactArea,
	//	Eigen::Vector3d& P1, Eigen::Vector3d& P2, bool ABD = false);



	// project edge-edge's gradient and hessian into full 12x1 vector and 12x12 matrix to faciliate mollifier mulitiplication
	void project_grad_to_full(Eigen::Vector2i& activePtsLocalInd, Vector6d& grad_, 
		Matrix6d& hess_, Vector12d& grad_full, Matrix12d& hess__full);
	void project_grad_to_full(Eigen::Vector3i& activePtsLocalInd, Vector9d& grad_, 
		Matrix9d& hess_, Vector12d& grad_full, Matrix12d& hess__full);


	double compute_b(double& d2, double& dHat2);

	double compute_g_b(double& d2, double& dHat2);

	double compute_H_b(double& d2, double& dHat2);


}


namespace Ground
{

	double val(double& coor_z2, double& contactArea, FEMParamters& parameters);

	void gradAndHess(std::vector<Eigen::Triplet<double>>& hessian_triplet, 
		std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_hess, 
		int& startIndex_grad, 
		int& index_i, double& coor_z, double& contactArea, FEMParamters& parameters, bool ABD = false);


}


#endif