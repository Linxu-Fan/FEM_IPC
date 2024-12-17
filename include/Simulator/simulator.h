#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "utils.h"  
#include "mesh.h" 
#include "CCD.h" 
#include "ExternalEnergy.h" 
#include "InertiaEnergy.h" 
#include "ElasticEnergy.h" 
#include "BarrierEnergy.h" 


//////////////////////////////////////////
// General FEM simulation
//////////////////////////////////////////

// explicit integration. Only used for one object simulation without contact
void explicitFEM(Mesh& tetSimMesh, FEMParamters& parameters);

// implicit integration
void implicitFEM(Mesh& tetSimMesh, FEMParamters& parameters);

// compute external force excluding gravity
Eigen::Vector3d compute_external_force(Mesh& tetSimMesh, int vertInd, int timestep);

// update the MLS points' information after advection
void updateMLS_after_advection(Mesh& tetSimMesh, FEMParamters& parameters);

// compute the incremental potential energy
double compute_IP_energy(Mesh& tetSimMesh, FEMParamters& parameters, int timestep);

// calculate a temporary vector for parallel computing MLS gradient and hessian
std::pair<std::vector<Eigen::Vector2i>, int> cal_temporary_MLS_tet_vector(Mesh& tetSimMesh); // Eigen::Vector2i (first int: index of the tetrahedral; second int: accumulative number of the MLS points in the tetrahedral); int: total number of nodes influenced by MLS points

// compute energy gradient and hessian of the linear system and solve the syetem
std::vector<Eigen::Vector3d> solve_linear_system(Mesh& tetSimMesh, FEMParamters& parameters, int timestep);

// move points' position according to the direction; Note this is a trial movement
void step_forward(FEMParamters& parameters, Mesh& tetSimMesh, 
	std::vector<Eigen::Vector3d>& currentPosition, std::vector<Eigen::Vector3d>& direction, double step);

// calculate the barrier energy of contact
double compute_Barrier_energy(Mesh& tetSimMesh, FEMParamters& parameters, int timestep);

// calculate the contact info: a) contact or not; b) contact type; 3) contact energy, derivative and hessian; 
void calContactInfo(Mesh& tetSimMesh, FEMParamters& parameters, int timestep, 
	std::vector<Vector5i>& PG_PG, std::vector<Vector5i>& PT_PP, 
	std::vector<Vector5i>& PT_PE, std::vector<Vector5i>& PT_PT, std::vector<Vector5i>& EE_EE);

// calculate the maximum feasible step size
double calMaxStepSize(Mesh& tetSimMesh, FEMParamters& parameters, int timestep,  
	std::vector<Eigen::Vector3d>& direction);







//////////////////////////////////////////
// ABD simulation
//////////////////////////////////////////

// implicit integration
void implicitFEM_ABD(Mesh_ABD& tetSimMesh, FEMParamters& parameters);

// construct q vector according translation and deformation
Vector12d constructQVec(Eigen::Vector3d translation, Eigen::Matrix3d deformation); 

// compute the incremental potential energy
double compute_IP_energy_ABD(Mesh_ABD& tetSimMesh, FEMParamters& parameters, int timestep);

// compute energy gradient and hessian of the linear system and solve the syetem
std::vector<Vector12d> solve_linear_system_ABD(Mesh_ABD& tetSimMesh, FEMParamters& parameters, int timestep);

// move points' position according to the direction; Note this is a trial movement
void step_forward_ABD(FEMParamters& parameters, Mesh_ABD& tetSimMesh, std::vector<Eigen::Vector3d>& current_ABD_translation,
	std::vector<Eigen::Matrix3d>& current_ABD_deformation, std::vector<Vector12d>& ABD_direction, std::vector<Eigen::Vector3d>& currentPosition, double step);

// convert the ABD state update to position update direction
void convert_to_position_direction(FEMParamters& parameters, Mesh_ABD& tetSimMesh, std::vector<Vector12d>& direction_ABD, std::vector<Eigen::Vector3d>& position_direction);


#endif
