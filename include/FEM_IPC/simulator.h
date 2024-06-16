#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "utils.h"  
#include "mesh.h" 

// implicit integration
void implicitFEM(Mesh& tetMesh, FEMParamters& parameters);

// compute external force excluding gravity
Eigen::Vector3d compute_external_force(Mesh& tetMesh, int vertInd, int timestep);

// compute the incremental potential energy
double compute_IP_energy(Mesh& tetMesh, FEMParamters& parameters, int timestep);

// compute energy gradient and hessian of the linear system and solve the syetem
std::vector<Eigen::Vector3d> solve_linear_system(Mesh& tetMesh, FEMParamters& parameters, int timestep);

// do line search to find the optimal step
void lineSearch(Mesh& tetMesh, std::vector<Eigen::Vector3d>& direction, FEMParamters& parameters, double& lastEnergyVal, int timestep);

// move points' position according to the direction; Note this is a trial movement
void step_forward(Mesh& tetMesh, std::vector<Eigen::Vector3d>& currentPosition, std::vector<Eigen::Vector3d>& direction, double step);

#endif
