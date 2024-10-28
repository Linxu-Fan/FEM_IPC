#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "utils.h"  
#include "mesh.h" 
#include "CCD.h" 
#include "ExternalEnergy.h" 
#include "InertiaEnergy.h" 
#include "ElasticEnergy.h" 
#include "BarrierEnergy.h" 

// explicit integration. Only used for one object simulation without contact
void explicitFEM(Mesh& tetSimMesh, FEMParamters& parameters);

// implicit integration
void implicitFEM(Mesh& tetSimMesh, FEMParamters& parameters);

// compute external force excluding gravity
Eigen::Vector3d compute_external_force(Mesh& tetSimMesh, int vertInd, int timestep);

// compute the incremental potential energy
double compute_IP_energy(Mesh& tetSimMesh, FEMParamters& parameters, int timestep);

// compute energy gradient and hessian of the linear system and solve the syetem
std::vector<Eigen::Vector3d> solve_linear_system(Mesh& tetSimMesh, FEMParamters& parameters, int timestep);

// move points' position according to the direction; Note this is a trial movement
void step_forward(FEMParamters& parameters, Mesh& tetSimMesh, std::vector<Eigen::Vector3d>& currentPosition, std::vector<Eigen::Vector3d>& direction, double step);

// calculate the barrier energy of contact
double compute_Barrier_energy(Mesh& tetSimMesh, FEMParamters& parameters, int timestep);

// calculate the contact info: a) contact or not; b) contact type; 3) contact energy, derivative and hessian; 
void calContactInfo(Mesh& tetSimMesh, FEMParamters& parameters, int timestep, std::vector<Vector5i>& PG_PG, std::vector<Vector5i>& PT_PP, std::vector<Vector5i>& PT_PE, std::vector<Vector5i>& PT_PT, std::vector<Vector5i>& EE_EE);

// calculate the maximum feasible step size
double calMaxStepSize(Mesh& tetSimMesh, FEMParamters& parameters, int timestep,  std::vector<Eigen::Vector3d>& direction);

#endif
