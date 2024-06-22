#include "BarrierEnergy.h"


// compute the distance between a point and a triangle
std::tuple<int, Eigen::Vector3d, double>  BarrierEnergy::pointTriangleDistance(Eigen::Vector3d P, Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C)
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
        std::tuple<int, Eigen::Vector3d, double> res = std::make_tuple(0, A, (P - A).norm());
        return res;
    }

    // Check if P in vertex region outside B
    Eigen::Vector3d BP = P - B;
    double d3 = AB.dot(BP);
    double d4 = AC.dot(BP);
    if (d3 >= 0.0 && d4 <= d3) 
    {
        std::tuple<int, Eigen::Vector3d, double> res = std::make_tuple(1, B, (P - B).norm());
        return res;
    }

    // Check if P in edge region of AB
    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) 
    {
        double v = d1 / (d1 - d3);
        Eigen::Vector3d closestPoint = A + v * AB;
        std::tuple<int, Eigen::Vector3d, double> res = std::make_tuple(3, closestPoint, (P - closestPoint).norm());
        return res;
    }

    // Check if P in vertex region outside C
    Eigen::Vector3d CP = P - C;
    double d5 = AB.dot(CP);
    double d6 = AC.dot(CP);
    if (d6 >= 0.0 && d5 <= d6) 
    {
        std::tuple<int, Eigen::Vector3d, double> res = std::make_tuple(2, C, (P - C).norm());
        return res;
    }

    // Check if P in edge region of AC
    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        double w = d2 / (d2 - d6);
        Eigen::Vector3d closestPoint = A + w * AC;
        std::tuple<int, Eigen::Vector3d, double> res = std::make_tuple(5, closestPoint, (P - closestPoint).norm());
        return res;
    }

    // Check if P in edge region of BC
    double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        Eigen::Vector3d closestPoint = B + w * (C - B);
        std::tuple<int, Eigen::Vector3d, double> res = std::make_tuple(4, closestPoint, (P - closestPoint).norm());
        return res;
    }

    // P inside face region. Compute Q through its barycentric coordinates (u, v, w)
    double denom = 1.0 / (va + vb + vc);
    double v = vb * denom;
    double w = vc * denom;
    Eigen::Vector3d closestPoint = A + AB * v + AC * w; // P is closest to the triangle's face
    std::tuple<int, Eigen::Vector3d, double> res = std::make_tuple(6, closestPoint, (P - closestPoint).norm());
    return res;
  
}


edgeEdgeDis BarrierEnergy::edgeEdgeDistance(Eigen::Vector3d P1, Eigen::Vector3d P2, Eigen::Vector3d Q1, Eigen::Vector3d Q2)
{
    edgeEdgeDis result;

    Eigen::Vector3d closestPointP,closestPointQ;

    // Vector along segment P1P2
    Eigen::Vector3d u = P2 - P1;
    // Vector along segment Q1Q2
    Eigen::Vector3d v = Q2 - Q1;
    // Vector from Q1 to P1
    Eigen::Vector3d w = P1 - Q1;

    // Coefficients for the equations
    double a = u.dot(u); // Squared length of segment P1P2
    double b = u.dot(v); // Dot product of u and v
    double c = v.dot(v); // Squared length of segment Q1Q2
    double d = u.dot(w); // Dot product of u and w
    double e = v.dot(w); // Dot product of v and w
    double D = a * c - b * b; // Determinant of the system, used to check parallelism

    // Numerators and denominators for the parameters s and t
    double sN, tN; // Numerators for s and t
    double sD = D, tD = D; // Denominators for s and t, initially set to D

    // Parameters for the closest points on the segments
    double sc, tc;

    // If the segments are almost parallel (D is very small)
    if (D < 1e-8) {
        sN = 0.0; // Force using the first point on segment P1P2
        sD = 1.0; // Avoid division by zero
        tN = e;   // Use the dot product of v and w
        tD = c;   // Use the squared length of segment Q1Q2

        result.nearParallel = true;
    }
    else {
        // Calculate the parameters of the closest points
        sN = (b * e - c * d);
        tN = (a * e - b * d);

        // Case 1: s < 0
        if (sN < 0.0) {
            sN = 0.0; // Clamp s to 0
            tN = e;   // Recompute tN for t parameter with s clamped to 0
            tD = c;   // Use the squared length of segment Q1Q2
        }
        // Case 2: s > 1
        else if (sN > sD) {
            sN = sD;  // Clamp s to 1 (since sD = D)
            tN = e + b; // Recompute tN for t parameter with s clamped to 1
            tD = c;   // Use the squared length of segment Q1Q2
        }
    }

    // Case 3: t < 0 (with s already clamped)
    if (tN < 0.0) {
        tN = 0.0;   // Clamp t to 0
        if (-d < 0.0) {
            sN = 0.0; // Clamp s to 0
        }
        else if (-d > a) {
            sN = sD;  // Clamp s to 1 (since sD = D)
        }
        else {
            sN = -d;  // Otherwise, set sN to -d
            sD = a;   // And set sD to the squared length of segment P1P2
        }
    }
    // Case 4: t > 1 (with s already clamped)
    else if (tN > tD) {
        tN = tD;    // Clamp t to 1 (since tD = D)
        if ((-d + b) < 0.0) {
            sN = 0; // Clamp s to 0
        }
        else if ((-d + b) > a) {
            sN = sD;  // Clamp s to 1 (since sD = D)
        }
        else {
            sN = (-d + b); // Otherwise, set sN to (-d + b)
            sD = a;        // And set sD to the squared length of segment P1P2
        }
    }

    // Divide to get the final sc and tc
    sc = (std::abs(sN) < 1e-8 ? 0.0 : sN / sD); // If sN is very small, set sc to 0, otherwise compute sc
    tc = (std::abs(tN) < 1e-8 ? 0.0 : tN / tD); // If tN is very small, set tc to 0, otherwise compute tc

    // Calculate the closest points on the segments using the parameters sc and tc
    closestPointP = P1 + sc * u; // Closest point on segment P1P2
    closestPointQ = Q1 + tc * v; // Closest point on segment Q1Q2

    // Handle the 9 cases explicitly
    if (sc == 0.0 && tc == 0.0) {
        // Case 1: Closest points are P1 and Q1
        closestPointP = P1;
        closestPointQ = Q1;

        result.edge1CloestPtInd = 0;
        result.edge2CloestPtInd = 0;
    }
    else if (sc == 0.0 && tc == 1.0) {
        // Case 2: Closest points are P1 and Q2
        closestPointP = P1;
        closestPointQ = Q2;

        result.edge1CloestPtInd = 0;
        result.edge2CloestPtInd = 2;
    }
    else if (sc == 1.0 && tc == 0.0) {
        // Case 3: Closest points are P2 and Q1
        closestPointP = P2;
        closestPointQ = Q1;

        result.edge1CloestPtInd = 2;
        result.edge2CloestPtInd = 0;
    }
    else if (sc == 1.0 && tc == 1.0) {
        // Case 4: Closest points are P2 and Q2
        closestPointP = P2;
        closestPointQ = Q2;

        result.edge1CloestPtInd = 2;
        result.edge2CloestPtInd = 2;
    }
    else if (sc == 0.0) {
        // Case 5: Closest point is P1 and an arbitrary point on Q1Q2
        closestPointP = P1;
        closestPointQ = Q1 + tc * v;

        result.edge1CloestPtInd = 0;
        result.edge2CloestPtInd = 1;
    }
    else if (sc == 1.0) {
        // Case 6: Closest point is P2 and an arbitrary point on Q1Q2
        closestPointP = P2;
        closestPointQ = Q1 + tc * v;

        result.edge1CloestPtInd = 2;
        result.edge2CloestPtInd = 1;
    }
    else if (tc == 0.0) {
        // Case 7: Closest point is an arbitrary point on P1P2 and Q1
        closestPointP = P1 + sc * u;
        closestPointQ = Q1;

        result.edge1CloestPtInd = 1;
        result.edge2CloestPtInd = 0;
    }
    else if (tc == 1.0) {
        // Case 8: Closest point is an arbitrary point on P1P2 and Q2
        closestPointP = P1 + sc * u;
        closestPointQ = Q2;

        result.edge1CloestPtInd = 1;
        result.edge2CloestPtInd = 2;
    }
    else {
        // Case 9: Both closest points are arbitrary points on their respective segments
        closestPointP = P1 + sc * u;
        closestPointQ = Q1 + tc * v;

        result.edge1CloestPtInd = 1;
        result.edge2CloestPtInd = 1;
    }

    // Vector between the closest points
    Eigen::Vector3d dP = closestPointP - closestPointQ;
    result.edge1CloestPtCoor = closestPointP;
    result.edge2CloestPtCoor = closestPointQ;
    result.distance = dP.norm();

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
