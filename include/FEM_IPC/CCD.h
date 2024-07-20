#ifndef CCD_H
#define CCD_H

#include "BarrierEnergy.h"

// Check if two edges' boundingbox intersect or not
// All points' coordinates are real coordinates not increments
bool edgeEdgeCCDBroadphase(const Eigen::Vector3d& P1, const Eigen::Vector3d& P2, const Eigen::Vector3d& dP1, const Eigen::Vector3d& dP2, const Eigen::Vector3d& Q1, const Eigen::Vector3d& Q2, const Eigen::Vector3d& dQ1, const Eigen::Vector3d& dQ2, double dist_threshold = 1.0E-3);

// check if a point and a triangle's boundingboxes intersect or not
bool pointTriangleCCDBroadphase(const Eigen::Vector3d& P, const Eigen::Vector3d& dP, const Eigen::Vector3d& A, const Eigen::Vector3d& dA, const Eigen::Vector3d& B, const Eigen::Vector3d& dB, const Eigen::Vector3d& C, const Eigen::Vector3d& dC, double dist_threshold = 1.0e-3);

// Calcualte the first time of contact: toc
// Actually calculate the time that first brings the distance to eta * currentDis
// !!!!!!!!!!!!!!Use incremental coordiante
double edgeEdgeCCDNarrowphase(const Eigen::Vector3d& P1, const Eigen::Vector3d& P2, const Eigen::Vector3d& dP1, const Eigen::Vector3d& dP2, const Eigen::Vector3d& Q1, const Eigen::Vector3d& Q2, const Eigen::Vector3d& dQ1, const Eigen::Vector3d& dQ2, double eta = 0.1);

double pointTriangleCCDNarrowphase(const Eigen::Vector3d& P, const Eigen::Vector3d& dP, const Eigen::Vector3d& A, const Eigen::Vector3d& dA, const Eigen::Vector3d& B, const Eigen::Vector3d& dB, const Eigen::Vector3d& C, const Eigen::Vector3d& dC, double eta = 0.1);

// spatial hash data structure
struct spatialHash
{
	std::set<int> vertIndices; // vertex whose boundingbox intersects with this cell
	std::set<int> edgeIndices; // edge whose boundingbox intersects with this cell
	std::set<int> triaIndices; // triangle whose boundingbox intersects with this cell
};

// initial spatial hash



#endif