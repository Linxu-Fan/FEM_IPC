#ifndef CCD_H
#define CCD_H

#include "BarrierEnergy.h"

// Check if two edges' boundingbox intersect or not
// All points' coordinates are real coordinates not increments
bool edgeEdgeCCDBroadphase(const Eigen::Vector3d& P1, const Eigen::Vector3d& P2,
	const Eigen::Vector3d& dP1, const Eigen::Vector3d& dP2, const Eigen::Vector3d& Q1,
	const Eigen::Vector3d& Q2, const Eigen::Vector3d& dQ1, const Eigen::Vector3d& dQ2,
	double dist_threshold);

// Check if two edges' boundingbox intersect or not (no advection)
bool edgeEdgeCCDBroadphase(const Eigen::Vector3d& P1, const Eigen::Vector3d& P2,
	const Eigen::Vector3d& Q1, const Eigen::Vector3d& Q2, double dist_threshold);



// check if a point and a triangle's boundingboxes intersect or not
bool pointTriangleCCDBroadphase(const Eigen::Vector3d& P,
	const Eigen::Vector3d& dP, const Eigen::Vector3d& A,
	const Eigen::Vector3d& dA, const Eigen::Vector3d& B,
	const Eigen::Vector3d& dB, const Eigen::Vector3d& C,
	const Eigen::Vector3d& dC, double dist_threshold);

// check if a point and a triangle's boundingboxes intersect or not (no advection)
bool pointTriangleCCDBroadphase(const Eigen::Vector3d& P,
	const Eigen::Vector3d& A, const Eigen::Vector3d& B,
	const Eigen::Vector3d& C, double dist_threshold);


// Calcualte the first time of contact: toc
// Actually calculate the time that first brings the distance to eta * currentDis
// !!!!!!!!!!!!!!Use incremental coordiante
double edgeEdgeCCDNarrowphase(const Eigen::Vector3d& P1,
	const Eigen::Vector3d& P2, const Eigen::Vector3d& dP1,
	const Eigen::Vector3d& dP2, const Eigen::Vector3d& Q1,
	const Eigen::Vector3d& Q2, const Eigen::Vector3d& dQ1,
	const Eigen::Vector3d& dQ2, double eta = 0.1);

double pointTriangleCCDNarrowphase(const Eigen::Vector3d& P,
	const Eigen::Vector3d& dP, const Eigen::Vector3d& A,
	const Eigen::Vector3d& dA, const Eigen::Vector3d& B,
	const Eigen::Vector3d& dB, const Eigen::Vector3d& C,
	const Eigen::Vector3d& dC, double eta = 0.1);

// spatial hash data structure
struct spatialHashCellData
{
	Eigen::Vector3i bottomLeftCorner = { 0,0,0 }; // bottom left corner of this spatial hash cell
	std::set<std::string> hashObjects; // store the names of objects within this hash
	std::set<int> vertIndices; // vertex whose boundingbox intersects with this cell
	std::set<int> edgeIndices; // edge whose boundingbox intersects with this cell
	std::vector<int> edgeIndices_vec; // edge whose boundingbox intersects with this cell but in vector form
	std::set<int> triaIndices; // triangle whose boundingbox intersects with this cell
};




// get the minimum and maximum corner of a vector of points
// bool advected: if consider the advection of points. If not, direction is ZERO
std::pair<Eigen::Vector3d, Eigen::Vector3d> getMinMaxCorner_boundaryVertices(
	const bool advected,
	const FEMParamters& parameters,
	objMeshFormat& surfaceInfo,
	const std::vector<Eigen::Vector3d>& direction,
	const std::vector<Eigen::Vector3d>& pos_node);


// get the minimum and maximum corner of a vector of points of each object
// bool advected: if consider the advection of points. If not, direction is ZERO
// we add another layer of cells outside of the original cells
std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> getMinMaxCorner_boundaryVertices_eachObject(
	const bool advected,
	const FEMParamters& parameters,
	objMeshFormat& surfaceInfo,
	const std::vector<Eigen::Vector3d>& direction,
	const std::vector<Eigen::Vector3d>& pos_node,
	const std::vector<std::string>& note_node,
	std::map<std::string, int>& tetMeshIndex);


// find the intersection of objects' hash
std::set<std::string> getTheIntersectedHash(
	const bool advected, 
	FEMParamters& parameters,
	objMeshFormat& surfaceInfo,
	const std::vector<Eigen::Vector3d>& direction,
	const std::vector<Eigen::Vector3d>& pos_node,
	const std::vector<std::string>& note_node,
	std::map<std::string, int>& tetMeshIndex,
	Eigen::Vector3d& minFloorGlobal);





void initSpatialHash(
	const bool advected,
	FEMParamters& parameters,
	objMeshFormat& surfaceInfo,
	const std::vector<Eigen::Vector3d>& direction,
	std::vector<spatialHashCellData>& spatialHash_vec,
	std::map<std::string, int>& hashNameIndex, 
	const std::vector<Eigen::Vector3d>& pos_node,
	const std::vector<std::string>& note_node,
	std::map<std::string, int>& tetMeshIndex, 
	const int timestep);


// use spatial hash to calculate the maximum feasible step
double calMaxStep_spatialHash(
	FEMParamters& parameters,
	objMeshFormat& surfaceInfo,
	const std::vector<Eigen::Vector3d>& direction,
	const std::vector<Eigen::Vector3d>& pos_node,
	const std::vector<std::string>& note_node,
	std::map<std::string, int>& tetMeshIndex,
	const int timestep);




#endif