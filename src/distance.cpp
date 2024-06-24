#include "Distance.h"


// calculate the type of point-triangle distance type
// A(0), B(1), C(2), AB(3), BC(4), CA(5), inside(6)
int pointTriangleDisType(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C)
{
    // Compute vectors
    Eigen::Vector3d AB = B - A;
    Eigen::Vector3d AC = C - A;
    Eigen::Vector3d AP = P - A;

    // Check if P in vertex region outside A
    double d1 = AB.dot(AP);
    double d2 = AC.dot(AP);
    if (d1 <= 0.0 && d2 <= 0.0)
    {
        return 0;
    }

    // Check if P in vertex region outside B
    Eigen::Vector3d BP = P - B;
    double d3 = AB.dot(BP);
    double d4 = AC.dot(BP);
    if (d3 >= 0.0 && d4 <= d3)
    {
        return 1;
    }

    // Check if P in edge region of AB
    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        return 3;
    }

    // Check if P in vertex region outside C
    Eigen::Vector3d CP = P - C;
    double d5 = AB.dot(CP);
    double d6 = AC.dot(CP);
    if (d6 >= 0.0 && d5 <= d6)
    {
        return 2;
    }

    // Check if P in edge region of AC
    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) 
    {
        return 5;
    }

    // Check if P in edge region of BC
    double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) 
    {
        return 4;
    }

    return 6;
}

// calcualte the SQUARED point-triangle distance
double pointTriangleDis2(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C)
{
    int type = pointTriangleDisType( P,  A,  B,  C);
    switch (type)
    {
    case 0:
        return pointPointDis2( P, A);
    case 1:
        return pointPointDis2(P, B);
    case 2:
        return pointPointDis2(P, C);
    case 3:
        return pointEdgeDis2(P, A, B);
    case 4:
        return pointEdgeDis2(P, B, C);
    case 5:
        return pointEdgeDis2(P, C, A);
    case 6:
        return pointTriangleDis2_inside(P, A, B, C);
    default:
        return pointTriangleDis2_inside(P, A, B, C);
    }

}

// calcualte the SQUARED point-point distance
double pointPointDis2(Eigen::Vector3d A, Eigen::Vector3d B)
{
    return (A - B).squaredNorm();
}

// calcualte the SQUARED point-edge distance
double pointEdgeDis2(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B)
{
    double d1 = (B - A).dot(P - A);
    double d3 = (B - A).dot(P - B);
    double v = d1 / (d1 - d3);
    Eigen::Vector3d closestPoint = A + v * (B - A);
    return (P - closestPoint).squaredNorm();   
}

// calcualte the SQUARED point-triangle distance (only when the point is projected into the interior of the triangle)
double pointTriangleDis2_inside(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C)
{
    Eigen::Vector3d norm = ((B - A).cross(C - A)).normalized();
    double dis = std::abs((P - A).dot(norm));
    return dis * dis;
}

