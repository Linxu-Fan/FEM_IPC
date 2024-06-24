#ifndef DISTANCE_H
#define DISTANCE_H

#include "utils.h"

// calculate the type of point-triangle distance type
// A(0), B(1), C(2), AB(3), BC(4), CA(5), inside(6)
int pointTriangleDisType(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C);

// calcualte the SQUARED point-triangle distance
double pointTriangleDis2(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C);

// calcualte the SQUARED point-point distance
double pointPointDis2(Eigen::Vector3d A, Eigen::Vector3d B);

// calcualte the SQUARED point-edge distance
double pointEdgeDis2(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B);

// calcualte the SQUARED point-triangle distance (only when the point is projected into the interior of the triangle)
double pointTriangleDis2_inside(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C);

#endif