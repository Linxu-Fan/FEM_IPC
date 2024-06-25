#ifndef BARRIERENERGY_H
#define BARRIERENERGY_H

#include "utils.h"
#include "distance.h"


// the result returned by pointTriangleDistance function
struct pointTriangleDis
{
	int cloestPt = 6; // closest point is A(0), B(1), C(2), AB(3), BC(4), CA(5), inside(6)
	Eigen::Vector3d cloestPtCoor = {0,0,0}; // closestPoint coordinate
	double distance = 0; // the distance between this point and triangle. We use SQUARED distance.
	// size of the following two vectors is THREE
	std::vector<Eigen::Vector3d> distanceGrad; // gradient of the distance wrt each point, including the current point and three points in the triangle
	std::vector<Eigen::Matrix3d> distanceHess; // hessian of the distance wrt each point

	// barrier energy, grad, and hess.
	double Val = 0;
	std::vector<std::pair<int, double>> Grad;
	std::vector<Eigen::Triplet<double>> Hess;
};


// the result returned by edgeEdgeDistance function
struct edgeEdgeDis
{
	bool nearParallel = false;

	double distance = 0; // = (edge1CloestPtCoor - edge2CloestPtCoor).squaredNorm(). We use SQUARED distance.
	// size of the following two vectors is FOUR
	std::vector<Eigen::Vector3d> distanceGrad; // gradient of the distance wrt each point
	std::vector<Eigen::Matrix3d> distanceHess; // hessian of the distance wrt each point

	// barrier energy, grad, and hess.
	double Val = 0;
	std::vector<std::pair<int, double>> Grad;
	std::vector<Eigen::Triplet<double>> Hess;
};


namespace BarrierEnergy
{

	// compute the energy between a point and a triangle. 
	pointTriangleDis pointTriangleBarriereEnergy(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C);

	// compute the energy between two edges in 3D
	edgeEdgeDis edgeEdgeBarriereEnergy(Eigen::Vector4i vtInd, Eigen::Vector3d P1, Eigen::Vector3d P2, Eigen::Vector3d Q1, Eigen::Vector3d Q2, double d_hat);


	// compute the barrier energy
	double Val(double d_ij, double d_hat, double k_stiff, double contactArea);

	// compute the energy gradient 
	std::vector<std::pair<int, double>> Grad(double d_ij, double d_hat, double k_stiff, double contactArea, int vertIndex, int BC);

	// compute the energy hessian 
	std::vector<Eigen::Triplet<double>> Hess(double d_ij, double d_hat, double k_stiff, double contactArea, int vertIndex, int BC);

}

#endif