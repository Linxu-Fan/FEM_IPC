#include "simulator.h"

Vector12d constructQVec(Eigen::Vector3d translation, Eigen::Matrix3d deformation) // construct q vector according translation and deformation
{
	Vector12d qb = Vector12d::Zero();
	qb.block(0, 0, 3, 1) = translation;
	qb.block(3, 0, 3, 1) = deformation.col(0);
	qb.block(6, 0, 3, 1) = deformation.col(1);
	qb.block(9, 0, 3, 1) = deformation.col(2);

	return qb;
}
