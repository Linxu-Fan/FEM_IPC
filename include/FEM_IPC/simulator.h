#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "utils.h" 
#include "materials.h" 
#include "mesh.h" 

// implicit integration
void implicitFEM(Mesh& tetMesh, Material& mat, FEMParamters& parameters);

// compute energy gradient and hessian of the linear system and solve the syetem
std::vector<Eigen::Vector3d> solve_linear_system(Mesh& tetMesh, Material& mat, FEMParamters& parameters, int timestep);

// do line search to find the optimal step
void lineSearch(Mesh& tetMesh, Material& mat, std::vector<Eigen::Vector3d>& direction, FEMParamters& parameters, double& lastEnergyVal, int timestep);

// move points' position according to the direction
void step_forward(Mesh& tetMesh, std::vector<Eigen::Vector3d>& currentPosition, std::vector<Eigen::Vector3d>& direction, double step);

#endif
