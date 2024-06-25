#include "BarrierEnergy.h"


// compute the distance between a point and a triangle
pointTriangleDis  BarrierEnergy::pointTriangleBarriereEnergy(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C)
{
    pointTriangleDis result;

    // Compute vectors
    Eigen::Vector3d AB = B - A;
    Eigen::Vector3d AC = C - A;
    Eigen::Vector3d AP = P - A;

    // Check if P in vertex region outside A
    double d1 = AB.dot(AP);
    double d2 = AC.dot(AP);
    if (d1 <= 0.0 && d2 <= 0.0) 
    {
        result.cloestPt = 0;
        result.cloestPtCoor = A;
        result.distance = (P - A).squaredNorm();
        return result;
    }

    // Check if P in vertex region outside B
    Eigen::Vector3d BP = P - B;
    double d3 = AB.dot(BP);
    double d4 = AC.dot(BP);
    if (d3 >= 0.0 && d4 <= d3) 
    {
        result.cloestPt = 1;
        result.cloestPtCoor = B;
        result.distance = (P - B).squaredNorm();
        return result;
    }

    // Check if P in edge region of AB
    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) 
    {
        double v = d1 / (d1 - d3);
        Eigen::Vector3d closestPoint = A + v * AB;

        result.cloestPt = 3;
        result.cloestPtCoor = closestPoint;
        result.distance = (P - closestPoint).squaredNorm();
        return result;
    }

    // Check if P in vertex region outside C
    Eigen::Vector3d CP = P - C;
    double d5 = AB.dot(CP);
    double d6 = AC.dot(CP);
    if (d6 >= 0.0 && d5 <= d6) 
    {

        result.cloestPt = 2;
        result.cloestPtCoor = C;
        result.distance = (P - C).squaredNorm();
        return result;
    }

    // Check if P in edge region of AC
    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        double w = d2 / (d2 - d6);
        Eigen::Vector3d closestPoint = A + w * AC;

        result.cloestPt = 5;
        result.cloestPtCoor = closestPoint;
        result.distance = (P - closestPoint).squaredNorm();
        return result;
    }

    // Check if P in edge region of BC
    double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        Eigen::Vector3d closestPoint = B + w * (C - B);

        result.cloestPt = 4;
        result.cloestPtCoor = closestPoint;
        result.distance = (P - closestPoint).squaredNorm();
        return result;
    }

    // P inside face region. Compute Q through its barycentric coordinates (u, v, w)
    double denom = 1.0 / (va + vb + vc);
    double v = vb * denom;
    double w = vc * denom;
    Eigen::Vector3d closestPoint = A + AB * v + AC * w; // P is closest to the triangle's face

    result.cloestPt = 6;
    result.cloestPtCoor = closestPoint;
    result.distance = (P - closestPoint).squaredNorm();
    return result;
  
}

// compute the cloest distance between two edges in 3D
edgeEdgeDis BarrierEnergy::edgeEdgeBarriereEnergy(Eigen::Vector4i vtInd, Eigen::Vector3d P1, Eigen::Vector3d P2, Eigen::Vector3d Q1, Eigen::Vector3d Q2, double d_hat)
{
    edgeEdgeDis result;

    int type = edgeEdgeDisType(P1, P2, Q1, Q2);
    disGradHess dis2_ = edgeEdgeDis2( type, P1, P2, Q1, Q2);
    double dis2 = dis2_.val;

    double d_hat2 = d_hat* d_hat;
    double b2 = -(dis2 - d_hat2) * (dis2 - d_hat2) * std::log(dis2 / (d_hat2));
    double partial_b2_partial_d = -4.0 * (dis2 - d_hat2) * std::sqrt(dis2) * std::log(dis2 / d_hat2) - 2.0 / std::sqrt(dis2) * (dis2 - d_hat2) * (dis2 - d_hat2);
    double partial2_b2_partial2_d = -4.0 * (3.0 * dis2 - d_hat2) * std::log(dis2 / d_hat2) - 8.0 * (dis2 - d_hat2) + 2.0 / dis2 * (dis2 - d_hat2) * (dis2 - d_hat2) - 8.0 * (dis2 - d_hat2);




    return result;
}

// compute the elastic energy
double BarrierEnergy::Val(double d_ij, double d_hat, double k_stiff, double contactArea)
{
	double energy = 0;
	return energy;
}

// compute the energy gradient 
std::vector<std::pair<int, double>> BarrierEnergy::Grad(double d_ij, double d_hat, double k_stiff, double contactArea, int vertIndex, int BC)
{
	Vector9d grad_tmp;
	std::vector<std::pair<int, double>> res;
	return  res;
}

// compute the energy hessian 
std::vector<Eigen::Triplet<double>> BarrierEnergy::Hess(double d_ij, double d_hat, double k_stiff, double contactArea, int vertIndex, int BC)
{

	std::vector<Eigen::Triplet<double>> res;
	return res;

}
