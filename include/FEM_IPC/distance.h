#ifndef DISTANCE_H
#define DISTANCE_H

#include "utils.h"


struct disGradHess
{
	double val = 0;
	Eigen::VectorXd grad;
	Eigen::MatrixXd hess;

	Vector6d grad_6;
	Vector9d grad_9;
	Vector12d grad_12;
	Matrix6d hess_6;
	Matrix9d hess_9;
	Matrix12d hess_12;

};


// calculate the type of point-triangle distance type
// A(0), B(1), C(2), AB(3), BC(4), CA(5), inside(6)
int pointTriangleDisType(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C);

// calcualte the SQUARED point-triangle distance
double pointTriangleDis2(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C);

// calcualte the SQUARED point-triangle distance given the type
double pointTriangleDis2(int type, Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C);
// calcualte the SQUARED point-triangle distance given the type gradient
double pointTriangleDis2(int type, Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C);
// calcualte the SQUARED point-triangle distance given the type hessian
double pointTriangleDis2(int type, Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C);


// calculate the type of edge-edge distance
// int 0: P1 and Q1; 1: P1 and Q1Q2; 2: P1 and Q2
// int 3: P1P2 and Q1; 4: P1P2 and Q1Q2; 5: P1P2 and Q2
// int 6: P2 and Q1; 7: P2 and Q1Q2; 8: P2 and Q2
int edgeEdgeDisType(Eigen::Vector3d P1, Eigen::Vector3d P2, Eigen::Vector3d Q1, Eigen::Vector3d Q2);

// calculate the distance between two edges
double edgeEdgeDis2(int type, Eigen::Vector3d P1, Eigen::Vector3d P2, Eigen::Vector3d Q1, Eigen::Vector3d Q2);

// calculate the distance between two edges
double edgeEdgeDis2(Eigen::Vector3d P1, Eigen::Vector3d P2, Eigen::Vector3d Q1, Eigen::Vector3d Q2);







// calcualte the SQUARED point-point distance value
double pointPointDis2_val(Eigen::Vector3d A, Eigen::Vector3d B);
// calcualte the SQUARED point-point distance gradient
Vector6d pointPointDis2_grad(Eigen::Vector3d A, Eigen::Vector3d B);
// calcualte the SQUARED point-point distance hessian
Matrix6d pointPointDis2_hess(Eigen::Vector3d A, Eigen::Vector3d B);



// calcualte the SQUARED point-edge distance value
double pointEdgeDis2_val(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B);
// calcualte the SQUARED point-edge distance gradient
Vector9d pointEdgeDis2_grad(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B);
// calcualte the SQUARED point-edge distance hessian
Matrix9d pointEdgeDis2_hess(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B);




// calcualte the SQUARED point-triangle distance (only when the point is projected into the interior of the triangle) value
double pointTriangleInsideDis2_val(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C);
// calcualte the SQUARED point-triangle distance (only when the point is projected into the interior of the triangle) gradient
Vector12d pointTriangleInsideDis2_grad(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C);
// calcualte the SQUARED point-triangle distance (only when the point is projected into the interior of the triangle) hessian
Matrix12d pointTriangleInsideDis2_hess(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C);





// calculate the distance between two edges (only when two cloest points are inside) value
double edgeEdgeInsideDis2_val(Eigen::Vector3d P1, Eigen::Vector3d P2, Eigen::Vector3d Q1, Eigen::Vector3d Q2);
// calculate the distance between two edges (only when two cloest points are inside) gradient
Vector12d edgeEdgeInsideDis2_grad(Eigen::Vector3d P1, Eigen::Vector3d P2, Eigen::Vector3d Q1, Eigen::Vector3d Q2);
// calculate the distance between two edges (only when two cloest points are inside) hessian
Matrix12d edgeEdgeInsideDis2_hess(Eigen::Vector3d P1, Eigen::Vector3d P2, Eigen::Vector3d Q1, Eigen::Vector3d Q2);






/////////////////////////////////////////////////////////////////////////
// The following codes calculate the edge-edge mollifier
/////////////////////////////////////////////////////////////////////////

// calculate the value of eps_x
double cal_EEM_eps_x(Eigen::Vector3d& P1_Rest, Eigen::Vector3d& P2_Rest, Eigen::Vector3d& Q1_Rest, Eigen::Vector3d& Q2_Rest);

// calculate the mollifier c
double cal_EEM_c(Eigen::Vector3d& P1, Eigen::Vector3d& P2, Eigen::Vector3d& Q1, Eigen::Vector3d& Q2);

// calculate the mollifier value
double cal_EEM_ek(double eps_x, double c);

// calculate the mollifier value wrt c
double cal_EEM_ek_wrt_c(double eps_x, double c);

// calculate the second-order derivative of mollifier value wrt c
double cal_EEM_2ek_wrt_c2(double eps_x, double c);

// calculate the mollifier c wrt x
Vector12d cal_EEM_c_wrt_x(Eigen::Vector3d& P1, Eigen::Vector3d& P2, Eigen::Vector3d& Q1, Eigen::Vector3d& Q2);

// calculate the second-order derivative of mollifier c wrt x
Matrix12d cal_EEM_2c_wrt_x2(Eigen::Vector3d& P1, Eigen::Vector3d& P2, Eigen::Vector3d& Q1, Eigen::Vector3d& Q2);






#endif