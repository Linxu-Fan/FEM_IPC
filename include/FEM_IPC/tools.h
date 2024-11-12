#ifndef TOOLS_H
#define TOOLS_H

#include "utils.h"  
#include "mesh.h" 
#include "CCD.h" 
#include "distance.h" 
#include "BarrierEnergy.h" 

Eigen::Vector3d compute_contact_force(objMeshFormat& obj1, objMeshFormat& obj2, Eigen::Vector3d& bbx_min, Eigen::Vector3d& bbx_max, FEMParamters& parameters);



#endif