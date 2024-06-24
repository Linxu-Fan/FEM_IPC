#ifndef BARRIERENERGY_H
#define BARRIERENERGY_H

#include "utils.h"


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
	int edge1CloestPtInd = 1; // closest point in edge1: 0: endpoint P1, 1: point inside P1P2, 2: endpoint P2
	int edge2CloestPtInd = 1; // closest point in edge2: 0: endpoint Q1, 1: point inside Q1Q2, 2: endpoint Q2
	Eigen::Vector3d edge1CloestPtCoor = { 0,0,0 }; // the coordinate of the edge1's colest point
	Eigen::Vector3d edge2CloestPtCoor = { 0,0,0 }; // the coordinate of the edge2's colest point
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

	// compute the distance between a point and a triangle. 
	// int: closest point is A(0), B(1), C(2), AB(3), BC(4), CA(5), inside(6)
	// Eigen::Vector3d: closestPoint coordinate
	// double: distance
	pointTriangleDis pointTriangleDistance(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C);

	// compute the cloest distance between two edges in 3D
	edgeEdgeDis edgeEdgeDistance(Eigen::Vector3d P1, Eigen::Vector3d P2, Eigen::Vector3d Q1, Eigen::Vector3d Q2);


	// compute the barrier energy
	double Val(double d_ij, double d_hat, double k_stiff, double contactArea);

	// compute the energy gradient 
	std::vector<std::pair<int, double>> Grad(double d_ij, double d_hat, double k_stiff, double contactArea, int vertIndex, int BC);

	// compute the energy hessian 
	std::vector<Eigen::Triplet<double>> Hess(double d_ij, double d_hat, double k_stiff, double contactArea, int vertIndex, int BC);

}

#endif