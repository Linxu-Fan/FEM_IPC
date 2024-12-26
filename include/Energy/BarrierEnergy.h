#ifndef BARRIERENERGY_H
#define BARRIERENERGY_H

#include "distance.h"
#include "mesh.h"

template<int size>
void assemble_gradAndHess(
	std::vector<Eigen::Triplet<double>>& hessian_triplet,
	std::vector<std::pair<int, double>>& grad_triplet,
	int& startIndex_hess,
	std::vector<int>& activePts,
	Eigen::Matrix<double, size * 3, 1>& grad_,
	Eigen::Matrix<double, size * 3, size * 3>& hess_,
	int& startIndex_grad);


template<int size>
void assemble_gradAndHess_ABD(
	std::vector<Eigen::Triplet<double>>& hessian_triplet,
	std::vector<std::pair<int, double>>& grad_triplet,
	int& startIndex_hess,
	std::vector<int>& activePts,
	Eigen::Matrix<double, size * 3, 1>& grad_,
	Eigen::Matrix<double, size * 3, size * 3>& hess_,
	int& startIndex_grad,
	const std::vector<Eigen::Vector3d>& pos_node,
	const std::vector<Eigen::Vector3d>& pos_node_Rest,
	const std::vector<Eigen::Vector2i>& index_node);


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
		const std::vector<Eigen::Vector3d>& pos_node, 
		const std::vector<Eigen::Vector3d>& pos_node_Rest,
		Eigen::Vector4i& ptIndices, 
		FEMParamters& parameters);

	// compute the energy gradient and hessian of point-triangle contact
	void gradAndHess_PT(
		double& contactArea,
		std::vector<Eigen::Triplet<double>>& hessian_triplet,
		std::vector<std::pair<int, double>>& grad_triplet,
		int& startIndex_hess,
		int& startIndex_grad,
		Eigen::Vector4i& ptIndices,
		int& type,
		double& dis2,
		const std::vector<Eigen::Vector3d>& pos_node,
		const std::vector<Eigen::Vector3d>& pos_node_Rest,
		const std::vector<Eigen::Vector2i>& index_node,
		std::vector<Eigen::Vector3d>& contactForce_node,
		FEMParamters& parameters,
		bool ABD = false);

	void cal_and_assemble_gradAndHess(
		std::vector<Eigen::Triplet<double>>& hessian_triplet,
		std::vector<std::pair<int, double>>& grad_triplet,
		int& startIndex_hess,
		std::vector<int>& activePtsLocalInd,
		int& startIndex_grad,
		const std::vector<Eigen::Vector3d>& pos_node,
		const std::vector<Eigen::Vector3d>& pos_node_Rest,
		const std::vector<Eigen::Vector2i>& index_node,
		std::vector<Eigen::Vector3d>& contactForce_node,
		FEMParamters& parameters,
		double& contactArea,
		double& g_bd,
		double& h_bd,
		bool ABD);



	// compute the energy gradient and hessian of edge-edge contact
	void gradAndHess_EE(
		double contactArea,
		std::vector<Eigen::Triplet<double>>& hessian_triplet,
		std::vector<std::pair<int, double>>& grad_triplet,
		int& startIndex_hess,
		int& startIndex_grad,
		Eigen::Vector4i& ptIndices,
		int& type,
		double& dis2,
		const std::vector<Eigen::Vector3d>& pos_node,
		const std::vector<Eigen::Vector3d>& pos_node_Rest,
		const std::vector<Eigen::Vector2i>& index_node,
		std::vector<Eigen::Vector3d>& contactForce_node,
		FEMParamters& parameters,
		bool ABD = false);


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
		const std::vector<Eigen::Vector3d>& pos_node,
		const std::vector<Eigen::Vector3d>& pos_node_Rest,
		const std::vector<Eigen::Vector2i>& index_node,
		int& index_i, 
		double& coor_z2, 
		double& contactArea,
		std::vector<Eigen::Vector3d>& contactForce_node,
		FEMParamters& parameters, 
		bool ABD = false);


}


#endif