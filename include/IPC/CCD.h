#ifndef CCD_H
#define CCD_H

#include "BarrierEnergy.h"

// Check if two edges' boundingbox intersect or not
// All points' coordinates are real coordinates not increments
bool edgeEdgeCCDBroadphase(const Eigen::Vector3d& P1, const Eigen::Vector3d& P2, const Eigen::Vector3d& dP1, 
	const Eigen::Vector3d& dP2, const Eigen::Vector3d& Q1, const Eigen::Vector3d& Q2, const Eigen::Vector3d& dQ1, const Eigen::Vector3d& dQ2, double dist_threshold);

// Check if two edges' boundingbox intersect or not (no advection)
bool edgeEdgeCCDBroadphase(const Eigen::Vector3d& P1, const Eigen::Vector3d& P2, const Eigen::Vector3d& Q1, const Eigen::Vector3d& Q2,  double dist_threshold);



// check if a point and a triangle's boundingboxes intersect or not
bool pointTriangleCCDBroadphase(const Eigen::Vector3d& P, const Eigen::Vector3d& dP, const Eigen::Vector3d& A, 
	const Eigen::Vector3d& dA, const Eigen::Vector3d& B, const Eigen::Vector3d& dB, const Eigen::Vector3d& C, const Eigen::Vector3d& dC, double dist_threshold);

// check if a point and a triangle's boundingboxes intersect or not (no advection)
bool pointTriangleCCDBroadphase(const Eigen::Vector3d& P,  const Eigen::Vector3d& A, const Eigen::Vector3d& B,  const Eigen::Vector3d& C, double dist_threshold);


// Calcualte the first time of contact: toc
// Actually calculate the time that first brings the distance to eta * currentDis
// !!!!!!!!!!!!!!Use incremental coordiante
double edgeEdgeCCDNarrowphase(const Eigen::Vector3d& P1, const Eigen::Vector3d& P2, const Eigen::Vector3d& dP1, 
	const Eigen::Vector3d& dP2, const Eigen::Vector3d& Q1, const Eigen::Vector3d& Q2, const Eigen::Vector3d& dQ1, const Eigen::Vector3d& dQ2, double eta = 0.1);

double pointTriangleCCDNarrowphase(const Eigen::Vector3d& P, const Eigen::Vector3d& dP, const Eigen::Vector3d& A, 
	const Eigen::Vector3d& dA, const Eigen::Vector3d& B, const Eigen::Vector3d& dB, const Eigen::Vector3d& C, const Eigen::Vector3d& dC, double eta = 0.1);

// spatial hash data structure
struct spatialHashCellData
{
	Eigen::Vector3i bottomLeftCorner = {0,0,0}; // bottom left corner of this spatial hash cell
	std::set<std::string> hashObjects; // store the names of objects within this hash
	std::set<int> vertIndices; // vertex whose boundingbox intersects with this cell
	std::set<int> edgeIndices; // edge whose boundingbox intersects with this cell
	std::vector<int> edgeIndices_vec; // edge whose boundingbox intersects with this cell but in vector form
	std::set<int> triaIndices; // triangle whose boundingbox intersects with this cell
};

// initialize spatial hash
void initSpatialHash(bool advected, FEMParamters& parameters, Mesh& tetSimMesh, std::vector<Eigen::Vector3d>& direction, double cellSize, std::vector<spatialHashCellData>& spatialHash_vec, std::map<std::string, int>& hashNameIndex, int timestep);

// use spatial hash to calculate the maximum feasible step
double calMaxStep_spatialHash(FEMParamters& parameters, Mesh& tetSimMesh, std::vector<Eigen::Vector3d>& direction, double cellSize, double dist_threshold, double eta);

#endif