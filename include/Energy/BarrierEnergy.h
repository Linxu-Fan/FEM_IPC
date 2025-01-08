#ifndef BARRIERENERGY_H
#define BARRIERENERGY_H

#include "distance.h"
#include "mesh.h"


template <int size>
void assemble_gradAndHess_ABD(
	std::vector<Eigen::Triplet<double>>& hessian_triplet,
	std::vector<std::pair<int, double>>& grad_triplet,
	int& startIndex_grad,
	int& startIndex_hess,
	std::vector<int>& pt_obj_index,
	Eigen::Matrix<double, size * 3, 1>& grad_,
	Eigen::Matrix<double, size * 3, size * 3>& hess_,
	const std::vector<Eigen::Vector3d>& pos_node_Rest);



namespace BarrierEnergy
{
	// compute the energy of point-triangle contact
	double val_PT(
		double& contactArea, 
		double& dis2, 
		FEMParamters& parameters);

	// compute the energy of edge-edge contact
	double val_EE(
		double& contactArea, 
		double& dis2, 
		const double val_ek,
		FEMParamters& parameters);


	void cal_gradAndHess_PP_ABD(
		Vector6d& grad,
		Matrix6d& hessian,
		const Eigen::Vector3d& pos1,
		const Eigen::Vector3d& pos2,
		std::vector<Eigen::Vector3d>& contactForce_node,
		FEMParamters& parameters,
		double& contactArea,
		double& dis2);


	void cal_gradAndHess_PE_ABD(
		Vector9d& grad,
		Matrix9d& hessian,
		const Eigen::Vector3d& pos1,
		const Eigen::Vector3d& pos2,
		const Eigen::Vector3d& pos3,
		std::vector<Eigen::Vector3d>& contactForce_node,
		FEMParamters& parameters,
		double& contactArea,
		double& dis2);

	void cal_gradAndHess_PT_ABD(
		Vector12d& grad,
		Matrix12d& hessian,
		const Eigen::Vector3d& pos1,
		const Eigen::Vector3d& pos2,
		const Eigen::Vector3d& pos3,
		const Eigen::Vector3d& pos4,
		std::vector<Eigen::Vector3d>& contactForce_node,
		FEMParamters& parameters,
		double& contactArea,
		double& dis2);


	// compute the energy gradient and hessian of point-triangle contact
	void gradAndHess_PT(
		std::vector<Eigen::Triplet<double>>& hessian_triplet,
		std::vector<std::pair<int, double>>& grad_triplet,
		int& startIndex_hess,
		int& startIndex_grad,
		double& dis2,
		const Vector5i& PT_Contact,
		triMesh& triSimMesh,
		std::vector<Eigen::Vector3d>& contactForce_node,
		FEMParamters& parameters);





	// compute the energy gradient and hessian of edge-edge contact
	void gradAndHess_EE(
		std::vector<Eigen::Triplet<double>>& hessian_triplet,
		std::vector<std::pair<int, double>>& grad_triplet,
		int& startIndex_hess,
		int& startIndex_grad,
		const Vector5i& EE_Contact,
		double& dis2,
		triMesh& triSimMesh,
		std::vector<Eigen::Vector3d>& contactForce_node,
		FEMParamters& parameters);


	// project edge-edge's gradient and hessian into full 12x1 vector and 12x12 matrix to faciliate mollifier mulitiplication
	template<int size>
	void project_grad_to_full(std::vector<int>& activePtsLocalInd, Eigen::Matrix<double, size * 3, 1>& grad_,
		Eigen::Matrix<double, size * 3, size * 3>& hess_, Vector12d& grad_full, Matrix12d& hess__full);


	double compute_b(double& d2, double& dHat2);

	double compute_g_b(double& d2, double& dHat2);

	double compute_H_b(double& d2, double& dHat2);


}


namespace Ground
{

	double val(double& coor_z2, double& contactArea, FEMParamters& parameters);

	void gradAndHess(
		std::vector<Eigen::Triplet<double>>& hessian_triplet,
		std::vector<std::pair<int, double>>& grad_triplet,
		int& startIndex_hess,
		int& startIndex_grad,
		const Vector2i PG_Contact,
		triMesh& triSimMesh,
		std::vector<Eigen::Vector3d>& contactForce_node,
		FEMParamters& parameters);


}


#endif