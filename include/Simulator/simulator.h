#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "utils.h"  
#include "CCD.h" 
#include "ExternalEnergy.h" 
#include "InertiaEnergy.h" 
#include "ElasticEnergy.h" 
#include "BarrierEnergy.h" 
#include "mpmSimulator.h"



// calculate the barrier energy of contact
double compute_Barrier_energy(
	triMesh& triSimMesh,
	FEMParamters& parameters,
	contact_Info& contact_pairs,
	const int timestep);



double calMaxStep(
	FEMParamters& parameters,
	triMesh& triSimMesh,
	contact_Info& contact_pairs,
	const int timestep);

// calculate the maximum feasible step size
double calMaxStepSize(
	FEMParamters& parameters,
	surface_Info& surfaceInfo,
	const std::vector<Eigen::Vector3d>& direction,
	const std::vector<Eigen::Vector3d>& pos_node,
	const std::vector<std::string>& note_node,
	std::map<std::string, int>& tetMeshIndex,
	const int timestep);



/**
 * @brief Find potential contact object pair in the bounding box level using the original bounding-box with transformation
 *
 * @note If in advect_mode, moving_direction is mandatory. Then the object's bounding box covers the potential moving trajectory.
 *			   If not in advect_mode, calculate the actual bounding box with/without moving_direction
 * 
 * @return vector of contact pairs. 1st int: index of the first object; 2nd int: index of the second object
 * 
 */
std::vector<std::pair<int, int>> find_contact_pair_BBX_level(
	const double& dilation, 
	triMesh& triSimMesh, 
	const std::vector<Vector12d>& moving_direction_ABD = {});

/**
 * @brief Find potential contact object pair in the BVH level 
 *
 * @note If in advect_mode, moving_direction is mandatory. Then the object's bounding box covers the potential moving trajectory.
 *			   If not in advect_mode, calculate the actual bounding box with/without moving_direction
 *
 * @return vector of contact pairs. 1st int: index of the first object; 2nd int: index of the second object
 *
 */
std::vector<std::pair<int, int>> find_contact_pair_BVH_level(
	const double& dilation,
	const std::vector<std::pair<int, int>>& BBX_pair,
	triMesh& triSimMesh,
	std::vector<int>& BVH_objects,
	bool advection);


void find_contact_pair_element_level(
	const std::vector<std::pair<int, int>>& BVH_pair,
	triMesh& triSimMesh,
	FEMParamters& parameters,
	contact_Info& contact_pairs);



void find_contact_pair_IPC_level(
	triMesh& triSimMesh,
	FEMParamters& parameters,
	contact_Info& contact_pairs);



// calculate the contact info
void find_contact(
	triMesh& triSimMesh,
	FEMParamters& parameters,
	contact_Info& contact_pairs,
	const std::vector<Vector12d>& moving_direction_ABD = {});




Vector12d constructQVec(Eigen::Vector3d translation, Eigen::Matrix3d deformation);

void implicitFEM_ABD_triMesh(triMesh& triSimMesh, FEMParamters& parameters);

double compute_IP_energy_ABD_triMesh(triMesh& triSimMesh, FEMParamters& parameters, int timestep, contact_Info& contact_pairs);

void solve_linear_system_ABD_triMesh(triMesh& triSimMesh, 
	FEMParamters& parameters, 
	int timestep, 
	std::vector<Vector12d>& movingDir,
	std::map<int, std::vector<Vector6d>>& broken_objects,
	contact_Info& contact_pairs,
	const int& iteration);

void step_forward_ABD_triMesh(FEMParamters& parameters, triMesh& triSimMesh, std::vector<Vector12d>& ABD_direction, double step);

double calculate_maximum_velocity(FEMParamters& parameters, triMesh& triSimMesh, std::vector<Vector12d>& direction_ABD);

/**
 * @brief check if start fracture simulation or not based on the contact force
 *
 * @param int: broken object's index; std::vector<Vector6d>: broken object's contact force
 */
void if_start_fracture_sim(triMesh& triSimMesh, std::map<int, std::vector<Vector6d>>& broken_objects);


void fracture_sim(FEMParamters& parameters, triMesh& triSimMesh, std::map<int, std::vector<Vector6d>>& broken_objects, std::map<int, objMeshFormat>& crackSurface_object);

/**
 * @brief cut objects with generated crack surfaces
 *
 */
void cut_object_with_cracks(FEMParamters& parameters, triMesh& triSimMesh, std::map<int, objMeshFormat>& crackSurface_object);



#endif
