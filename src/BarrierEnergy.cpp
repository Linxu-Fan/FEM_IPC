#include "BarrierEnergy.h"


// compute the barrier energy
double BarrierEnergy::Val(bool pointTriangle, double dis2, Mesh& tetMesh, Eigen::Vector4i& vtInd, double d_hat, double k_stiff, double dt)
{
    double result = 0;
    // calculate the point-triangle barrier energy
    if (pointTriangle)
    {
        int pt = vtInd[0];
        double contactArea = tetMesh.boundaryVertices_area[pt];
        // calculate the barrier energy
        double bEnergy = - squaredDouble(dis2 - d_hat * d_hat) * log(dis2 / (d_hat * d_hat));
        result = dt * dt * k_stiff * contactArea * bEnergy;
    }
    else
    {
        int e11 = vtInd[0], e12 = vtInd[1], e21 = vtInd[2], e22 = vtInd[3];
        int emin = std::min(e11, e12), emax = std::max(e11, e12);
        double contactArea = tetMesh.boundaryEdges_area[emin][emax];
        // calculate the barrier energy
        double bEnergy = -squaredDouble(dis2 - d_hat * d_hat) * log(dis2 / (d_hat * d_hat));
        result = dt * dt * k_stiff * contactArea * bEnergy;
    }

	return result;
}

//// compute the energy gradient 
//std::vector<std::pair<int, double>> BarrierEnergy::Grad(bool pointTriangle, int type, double dis2, Mesh& tetMesh, Eigen::Vector4i& vtInd, Eigen::Vector4i& vtInd_BC, double d_hat, double k_stiff, double dt)
//{
//    std::vector<std::pair<int, double>> res;
//
//    // the partial derivative of barrier energy b wrt distance d
//    double d_hat2 = squaredDouble(d_hat);
//    double partial_b_partial_d2is = -2.0 * (dis2 - d_hat2) * log(dis2 / d_hat2) - 1.0 / dis2 * squaredDouble(dis2 - d_hat2);
//    if (pointTriangle)
//    {
//        int pt = vtInd[0], t1 = vtInd[1], t2 = vtInd[2], t3 = vtInd[3];
//        double contactArea = tetMesh.boundaryVertices_area[pt];
//        Eigen::Vector3d P = tetMesh.pos_node[pt], A = tetMesh.pos_node[t1], B = tetMesh.pos_node[t2], C = tetMesh.pos_node[t3];
//        std::vector<int> activePts, activePtsBC;
//
//        switch (type)
//        {
//        case 0:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t1);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[1]);
//
//            Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P, A) * partial_b_partial_d2is;
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//                }
//            }
//        }
//        break;
//        case 1:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t2);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[2]);
//
//            Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P, B) * partial_b_partial_d2is;
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//                }
//            }
//        }       
//        break;
//
//        case 2:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t3);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P, C) * partial_b_partial_d2is;
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//                }
//            }
//        }
//        break;
//
//        case 3:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t1);
//            activePts.push_back(t2);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[1]);
//            activePtsBC.push_back(vtInd_BC[2]);
//
//            Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(P, A, B) * partial_b_partial_d2is;
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//                }
//            }
//        }
//        break;
//
//        case 4:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t2);
//            activePts.push_back(t3);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[2]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(P, B, C) * partial_b_partial_d2is;
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//                }
//            }
//        }
//        break;
//
//        case 5:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t3);
//            activePts.push_back(t1);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[3]);
//            activePtsBC.push_back(vtInd_BC[1]);
//
//            Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(P, C, A) * partial_b_partial_d2is;
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//                }
//            }
//        }
//        break;
//
//        case 6:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t1);
//            activePts.push_back(t2);
//            activePts.push_back(t3);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[1]);
//            activePtsBC.push_back(vtInd_BC[2]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector12d grad_ = dt * dt * k_stiff * contactArea * pointTriangleInsideDis2_grad(P, A, B, C) * partial_b_partial_d2is;
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//                }
//            }
//        }
//        break;
//
//        }
//    
//
//    }
//    else
//    {
//        int P1 = vtInd[0], P2 = vtInd[1], Q1 = vtInd[2], Q2 = vtInd[3];
//        int emin = std::min(P1, P2), emax = std::max(P1, P2);
//        double contactArea = tetMesh.boundaryEdges_area[emin][emax];
//        Eigen::Vector3d P1Coor = tetMesh.pos_node[P1], P2Coor = tetMesh.pos_node[P2], Q1Coor = tetMesh.pos_node[Q1], Q2Coor = tetMesh.pos_node[Q2];
//        std::vector<int> activePts, activePtsBC;
//
//        switch (type)
//        {
//        case 0:
//        {
//            activePts.push_back(P1);
//            activePts.push_back(Q1);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[2]);
//
//            Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P1Coor, Q1Coor) * partial_b_partial_d2is;
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//                }
//            }
//        }
//        break;
//            
//        case 1:
//        {
//            activePts.push_back(P1);
//            activePts.push_back(Q1);
//            activePts.push_back(Q2);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[2]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(P1Coor, Q1Coor, Q2Coor) * partial_b_partial_d2is;
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//                }
//            }
//        }
//        break;
//
//            
//        case 2:
//        {
//            activePts.push_back(P1);
//            activePts.push_back(Q2);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P1Coor, Q2Coor) * partial_b_partial_d2is;
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//                }
//            }
//        }
//        break;
//            
//        case 3:
//        {
//            activePts.push_back(Q1);
//            activePts.push_back(P1);
//            activePts.push_back(P2);
//            activePtsBC.push_back(vtInd_BC[2]);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[1]);
//
//            Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(Q1Coor, P1Coor, P2Coor) * partial_b_partial_d2is;
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//                }
//            }
//        }
//        break;
//
//            
//        case 4:
//        {
//            activePts.push_back(P1);
//            activePts.push_back(P2);
//            activePts.push_back(Q1);
//            activePts.push_back(Q2);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[1]);
//            activePtsBC.push_back(vtInd_BC[2]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//
//            Vector12d grad_ = dt * dt * k_stiff * contactArea * edgeEdgeInsideDis2_grad(P1Coor, P2Coor, Q1Coor, Q2Coor) * partial_b_partial_d2is;
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//                }
//            }
//        }
//        break;
//
//            
//        case 5:
//        {
//            activePts.push_back(Q2);
//            activePts.push_back(P1);
//            activePts.push_back(P2);
//            activePtsBC.push_back(vtInd_BC[3]);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[1]);
//
//            Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(Q2Coor, P1Coor, P2Coor) * partial_b_partial_d2is;
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//                }
//            }
//        }
//        break;
//
//            
//        case 6:
//        {
//            activePts.push_back(P2);
//            activePts.push_back(Q1);
//            activePtsBC.push_back(vtInd_BC[1]);
//            activePtsBC.push_back(vtInd_BC[2]);
//
//            Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P2Coor, Q1Coor) * partial_b_partial_d2is;
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//                }
//            }
//        }
//        break;
//
//            
//        case 7:
//        {
//            activePts.push_back(P2);
//            activePts.push_back(Q1);
//            activePts.push_back(Q2);
//            activePtsBC.push_back(vtInd_BC[1]);
//            activePtsBC.push_back(vtInd_BC[2]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(P2Coor, Q1Coor, Q2Coor) * partial_b_partial_d2is;
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//                }
//            }
//        }
//        break;
//
//            
//        case 8:
//        {
//            activePts.push_back(P2);
//            activePts.push_back(Q2);
//            activePtsBC.push_back(vtInd_BC[1]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P2Coor, Q2Coor) * partial_b_partial_d2is;
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//                }
//            }
//        }
//        break;
//
//           
//        }
//    }
//
//	return  res;
//}
//
//// compute the energy hessian 
//std::vector<Eigen::Triplet<double>> BarrierEnergy::Hess(bool pointTriangle, int type, double dis2, Mesh& tetMesh, Eigen::Vector4i& vtInd, Eigen::Vector4i& vtInd_BC, double d_hat, double k_stiff, double dt)
//{
//
//	std::vector<Eigen::Triplet<double>> res;
//    // the partial derivative of barrier energy b wrt distance d
//    double d = std::sqrt(dis2);
//    double d_hat2 = squaredDouble(d_hat);
//    double partial_b_partial_d2is = -2.0 * (dis2 - d_hat2) * log(dis2 / d_hat2) - 1.0 / dis2 * squaredDouble(dis2 - d_hat2);
//    double partial_b_partial_d = -4.0 * d * (dis2 - d_hat2) * log(dis2 / d_hat2) - 2.0 / d * squaredDouble(dis2 - d_hat2);
//    double partial_2b_partial_d2 = -(4.0 * (3.0 * d * dis2 - d_hat2) * log(dis2 / d_hat2) - 2.0 / dis2 * squaredDouble(dis2 - d_hat2) + 16.0 * (dis2 - d_hat2));
//
//
//    if (pointTriangle)
//    {
//        int pt = vtInd[0], t1 = vtInd[1], t2 = vtInd[2], t3 = vtInd[3];
//        double contactArea = tetMesh.boundaryVertices_area[pt];
//        Eigen::Vector3d P = tetMesh.pos_node[pt], A = tetMesh.pos_node[t1], B = tetMesh.pos_node[t2], C = tetMesh.pos_node[t3];
//        std::vector<int> activePts, activePtsBC;
//
//        switch (type)
//        {
//        case 0:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t1);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[1]);
//
//            Vector6d grad_ = pointPointDis2_grad(P, A);
//            Matrix6d hess_ = pointPointDis2_hess(P, A);
//
//            Matrix6d energyHess = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2 * grad_ * grad_.transpose() + partial_b_partial_d * hess_);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = energyHess(m * 3 + xd, n * 3 + yd);
//                                    res.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//
//        }
//        break;
//
//
//        case 1:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t2);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[2]);
//
//            Vector6d grad_ = pointPointDis2_grad(P, B);
//            Matrix6d hess_ = pointPointDis2_hess(P, B);
//
//            Matrix6d energyHess = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2 * grad_ * grad_.transpose() + partial_b_partial_d * hess_);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = energyHess(m * 3 + xd, n * 3 + yd);
//                                    res.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//
//        }
//        break;
//           
//        case 2:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t3);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector6d grad_ = pointPointDis2_grad(P, C);
//            Matrix6d hess_ = pointPointDis2_hess(P, C);
//
//            Matrix6d energyHess = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2 * grad_ * grad_.transpose() + partial_b_partial_d * hess_);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = energyHess(m * 3 + xd, n * 3 + yd);
//                                    res.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//            
//        case 3:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t1);
//            activePts.push_back(t2);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[1]);
//            activePtsBC.push_back(vtInd_BC[2]);
//
//            Vector9d grad_ = pointEdgeDis2_grad(P, A, B);
//            Matrix9d hess_ = pointEdgeDis2_hess(P, A, B);
//
//            Matrix9d energyHess = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2 * grad_ * grad_.transpose() + partial_b_partial_d * hess_);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = energyHess(m * 3 + xd, n * 3 + yd);
//                                    res.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//            
//
//        case 4:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t2);
//            activePts.push_back(t3);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[2]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector9d grad_ = pointEdgeDis2_grad(P, B, C);
//            Matrix9d hess_ = pointEdgeDis2_hess(P, B, C);
//
//            Matrix9d energyHess = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2 * grad_ * grad_.transpose() + partial_b_partial_d * hess_);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = energyHess(m * 3 + xd, n * 3 + yd);
//                                    res.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//            
//
//        case 5:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t3);
//            activePts.push_back(t1);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[3]);
//            activePtsBC.push_back(vtInd_BC[1]);
//
//            Vector9d grad_ = pointEdgeDis2_grad(P, C, A);
//            Matrix9d hess_ = pointEdgeDis2_hess(P, C, A);
//
//            Matrix9d energyHess = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2 * grad_ * grad_.transpose() + partial_b_partial_d * hess_);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = energyHess(m * 3 + xd, n * 3 + yd);
//                                    res.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//            
//
//        case 6:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t1);
//            activePts.push_back(t2);
//            activePts.push_back(t3);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[1]);
//            activePtsBC.push_back(vtInd_BC[2]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector12d grad_ = pointTriangleInsideDis2_grad(P, A, B, C);
//            Matrix12d hess_ = pointTriangleInsideDis2_hess(P, A, B, C);
//
//            Matrix12d energyHess = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2 * grad_ * grad_.transpose() + partial_b_partial_d * hess_);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = energyHess(m * 3 + xd, n * 3 + yd);
//                                    res.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//            
//
//        }
//
//
//    }
//    else
//    {
//        int P1 = vtInd[0], P2 = vtInd[1], Q1 = vtInd[2], Q2 = vtInd[3];
//        int emin = std::min(P1, P2), emax = std::max(P1, P2);
//        double contactArea = tetMesh.boundaryEdges_area[emin][emax];
//        Eigen::Vector3d P1Coor = tetMesh.pos_node[P1], P2Coor = tetMesh.pos_node[P2], Q1Coor = tetMesh.pos_node[Q1], Q2Coor = tetMesh.pos_node[Q2];
//        std::vector<int> activePts, activePtsBC;
//
//        switch (type)
//        {
//        case 0:
//        {
//            activePts.push_back(P1);
//            activePts.push_back(Q1);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[2]);
//
//            Vector6d grad_ = pointPointDis2_grad(P1Coor, Q1Coor);
//            Matrix6d hess_ = pointPointDis2_hess(P1Coor, Q1Coor);
//
//            Matrix6d energyHess = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2 * grad_ * grad_.transpose() + partial_b_partial_d * hess_);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = energyHess(m * 3 + xd, n * 3 + yd);
//                                    res.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//            
//
//        case 1:
//        {
//            activePts.push_back(P1);
//            activePts.push_back(Q1);
//            activePts.push_back(Q2);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[2]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector9d grad_ = pointEdgeDis2_grad(P1Coor, Q1Coor, Q2Coor);
//            Matrix9d hess_ = pointEdgeDis2_hess(P1Coor, Q1Coor, Q2Coor);
//
//            Matrix9d energyHess = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2 * grad_ * grad_.transpose() + partial_b_partial_d * hess_);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = energyHess(m * 3 + xd, n * 3 + yd);
//                                    res.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//            
//
//
//        case 2:
//        {
//            activePts.push_back(P1);
//            activePts.push_back(Q2);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector6d grad_ = pointPointDis2_grad(P1Coor, Q2Coor);
//            Matrix6d hess_ = pointPointDis2_hess(P1Coor, Q2Coor);
//
//            Matrix6d energyHess = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2 * grad_ * grad_.transpose() + partial_b_partial_d * hess_);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = energyHess(m * 3 + xd, n * 3 + yd);
//                                    res.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//            
//
//        case 3:
//        {
//            activePts.push_back(Q1);
//            activePts.push_back(P1);
//            activePts.push_back(P2);
//            activePtsBC.push_back(vtInd_BC[2]);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[1]);
//
//            Vector9d grad_ = pointEdgeDis2_grad(Q1Coor, P1Coor, P2Coor);
//            Matrix9d hess_ = pointEdgeDis2_hess(Q1Coor, P1Coor, P2Coor);
//
//            Matrix9d energyHess = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2 * grad_ * grad_.transpose() + partial_b_partial_d * hess_);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = energyHess(m * 3 + xd, n * 3 + yd);
//                                    res.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//            
//
//
//        case 4:
//        {
//            activePts.push_back(P1);
//            activePts.push_back(P2);
//            activePts.push_back(Q1);
//            activePts.push_back(Q2);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[1]);
//            activePtsBC.push_back(vtInd_BC[2]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//
//            Vector12d grad_ = edgeEdgeInsideDis2_grad(P1Coor, P2Coor, Q1Coor, Q2Coor);
//            Matrix12d hess_ = edgeEdgeInsideDis2_hess(P1Coor, P2Coor, Q1Coor, Q2Coor);
//
//            Matrix12d energyHess = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2 * grad_ * grad_.transpose() + partial_b_partial_d * hess_);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = energyHess(m * 3 + xd, n * 3 + yd);
//                                    res.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//
//        }
//        break;
//            
//
//        case 5:
//        {
//            activePts.push_back(Q2);
//            activePts.push_back(P1);
//            activePts.push_back(P2);
//            activePtsBC.push_back(vtInd_BC[3]);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[1]);
//
//            Vector9d grad_ = pointEdgeDis2_grad(Q2Coor, P1Coor, P2Coor);
//            Matrix9d hess_ = pointEdgeDis2_hess(Q2Coor, P1Coor, P2Coor);
//
//            Matrix9d energyHess = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2 * grad_ * grad_.transpose() + partial_b_partial_d * hess_);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = energyHess(m * 3 + xd, n * 3 + yd);
//                                    res.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//            
//
//
//        case 6:
//        {
//            activePts.push_back(P2);
//            activePts.push_back(Q1);
//            activePtsBC.push_back(vtInd_BC[1]);
//            activePtsBC.push_back(vtInd_BC[2]);
//
//            Vector6d grad_ = pointPointDis2_grad(P2Coor, Q1Coor);
//            Matrix6d hess_ = pointPointDis2_hess(P2Coor, Q1Coor);
//
//            Matrix6d energyHess = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2 * grad_ * grad_.transpose() + partial_b_partial_d * hess_);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = energyHess(m * 3 + xd, n * 3 + yd);
//                                    res.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//            
//
//
//        case 7:
//        {
//            activePts.push_back(P2);
//            activePts.push_back(Q1);
//            activePts.push_back(Q2);
//            activePtsBC.push_back(vtInd_BC[1]);
//            activePtsBC.push_back(vtInd_BC[2]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector9d grad_ = pointEdgeDis2_grad(P2Coor, Q1Coor, Q2Coor);
//            Matrix9d hess_ = pointEdgeDis2_hess(P2Coor, Q1Coor, Q2Coor);
//
//            Matrix9d energyHess = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2 * grad_ * grad_.transpose() + partial_b_partial_d * hess_);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = energyHess(m * 3 + xd, n * 3 + yd);
//                                    res.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//            
//
//
//        case 8:
//        {
//            activePts.push_back(P2);
//            activePts.push_back(Q2);
//            activePtsBC.push_back(vtInd_BC[1]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector6d grad_ = pointPointDis2_grad(P2Coor, Q2Coor);
//            Matrix6d hess_ = pointPointDis2_hess(P2Coor, Q2Coor);
//
//            Matrix6d energyHess = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2 * grad_ * grad_.transpose() + partial_b_partial_d * hess_);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = energyHess(m * 3 + xd, n * 3 + yd);
//                                    res.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//
//        }
//        break;
//            
//
//
//        }
//    }
//
//    return  res;
//}
//
//// compute the energy gradient and hessian 
//std::pair<std::vector<std::pair<int, double>>, std::vector<Eigen::Triplet<double>>> BarrierEnergy::gradAndHess(bool pointTriangle, int type,
//    double dis2, Mesh& tetMesh, Eigen::Vector4i& vtInd, Eigen::Vector4i& vtInd_BC, double d_hat, double k_stiff, double dt)
//{
//    std::vector<std::pair<int, double>> res_grad;
//    std::vector<Eigen::Triplet<double>> res_hess;
//
//    // the partial derivative of barrier energy b wrt distance d
//    double d_hat2 = squaredDouble(d_hat);
//    double partial_b_partial_d2is = -2.0 * (dis2 - d_hat2) * log(dis2 / d_hat2) - 1.0 / dis2 * squaredDouble(dis2 - d_hat2); // 3
//    double partial_2b_partial_d2is2 = -2.0 * log(dis2 / d_hat2) - 4.0 / dis2 * (dis2 - d_hat2) + squaredDouble(dis2 - d_hat2) / dis2 / dis2; // 1                                                        
//
//
//
//    if (pointTriangle)
//    {
//        int pt = vtInd[0], t1 = vtInd[1], t2 = vtInd[2], t3 = vtInd[3];
//        double contactArea = tetMesh.boundaryVertices_area[pt];
//        Eigen::Vector3d P = tetMesh.pos_node[pt], A = tetMesh.pos_node[t1], B = tetMesh.pos_node[t2], C = tetMesh.pos_node[t3];
//        std::vector<int> activePts, activePtsBC;
//
//        switch (type)
//        {
//        case 0:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t1);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[1]);
//
//            Vector6d partial_d2is_partial_x = pointPointDis2_grad(P, A); // 2
//            Matrix6d partial_2d2is_partial_x2 = pointPointDis2_hess(P, A); // 4
//
//            Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P, A) * partial_b_partial_d2is;
//            Matrix6d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res_grad.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//
//
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = hess_(m * 3 + xd, n * 3 + yd);
//                                    res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//        case 1:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t2);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[2]);
//
//            Vector6d partial_d2is_partial_x = pointPointDis2_grad(P, B); // 2
//            Matrix6d partial_2d2is_partial_x2 = pointPointDis2_hess(P, B); // 4
//
//            Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P, B) * partial_b_partial_d2is;
//            Matrix6d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res_grad.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//
//
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = hess_(m * 3 + xd, n * 3 + yd);
//                                    res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//
//        case 2:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t3);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector6d partial_d2is_partial_x = pointPointDis2_grad(P, C); // 2
//            Matrix6d partial_2d2is_partial_x2 = pointPointDis2_hess(P, C); // 4
//
//            Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P, C) * partial_b_partial_d2is;
//            Matrix6d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res_grad.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//
//
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = hess_(m * 3 + xd, n * 3 + yd);
//                                    res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//
//        case 3:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t1);
//            activePts.push_back(t2);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[1]);
//            activePtsBC.push_back(vtInd_BC[2]);
//
//            Vector9d partial_d2is_partial_x = pointEdgeDis2_grad(P, A, B); // 2
//            Matrix9d partial_2d2is_partial_x2 = pointEdgeDis2_hess(P, A, B); // 4
//
//            Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(P, A, B) * partial_b_partial_d2is;
//            Matrix9d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res_grad.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//
//
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = hess_(m * 3 + xd, n * 3 + yd);
//                                    res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//
//        case 4:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t2);
//            activePts.push_back(t3);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[2]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector9d partial_d2is_partial_x = pointEdgeDis2_grad(P, B, C); // 2
//            Matrix9d partial_2d2is_partial_x2 = pointEdgeDis2_hess(P, B, C); // 4
//
//            Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(P, B, C) * partial_b_partial_d2is;
//            Matrix9d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res_grad.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//
//
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = hess_(m * 3 + xd, n * 3 + yd);
//                                    res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//
//        case 5:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t3);
//            activePts.push_back(t1);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[3]);
//            activePtsBC.push_back(vtInd_BC[1]);
//
//            Vector9d partial_d2is_partial_x = pointEdgeDis2_grad(P, C, A); // 2
//            Matrix9d partial_2d2is_partial_x2 = pointEdgeDis2_hess(P, C, A); // 4
//
//            Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(P, C, A) * partial_b_partial_d2is;
//            Matrix9d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res_grad.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//
//
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = hess_(m * 3 + xd, n * 3 + yd);
//                                    res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//
//        case 6:
//        {
//            activePts.push_back(pt);
//            activePts.push_back(t1);
//            activePts.push_back(t2);
//            activePts.push_back(t3);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[1]);
//            activePtsBC.push_back(vtInd_BC[2]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector12d partial_d2is_partial_x = pointTriangleInsideDis2_grad(P, A, B, C); // 2
//            Matrix12d partial_2d2is_partial_x2 = pointTriangleInsideDis2_hess(P, A, B, C); // 4
//
//            Vector12d grad_ = dt * dt * k_stiff * contactArea * pointTriangleInsideDis2_grad(P, A, B, C) * partial_b_partial_d2is;
//            Matrix12d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res_grad.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//
//
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = hess_(m * 3 + xd, n * 3 + yd);
//                                    res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//
//        }
//
//
//    }
//    else
//    {
//        int P1 = vtInd[0], P2 = vtInd[1], Q1 = vtInd[2], Q2 = vtInd[3];
//        int emin = std::min(P1, P2), emax = std::max(P1, P2);
//        double contactArea = tetMesh.boundaryEdges_area[emin][emax];
//        Eigen::Vector3d P1Coor = tetMesh.pos_node[P1], P2Coor = tetMesh.pos_node[P2], Q1Coor = tetMesh.pos_node[Q1], Q2Coor = tetMesh.pos_node[Q2];
//        std::vector<int> activePts, activePtsBC;
//
//        switch (type)
//        {
//        case 0:
//        {
//            activePts.push_back(P1);
//            activePts.push_back(Q1);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[2]);
//
//            Vector6d partial_d2is_partial_x = pointPointDis2_grad(P1Coor, Q1Coor); // 2
//            Matrix6d partial_2d2is_partial_x2 = pointPointDis2_hess(P1Coor, Q1Coor); // 4
//
//            Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P1Coor, Q1Coor) * partial_b_partial_d2is;
//            Matrix6d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res_grad.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//
//
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = hess_(m * 3 + xd, n * 3 + yd);
//                                    res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//
//        case 1:
//        {
//            activePts.push_back(P1);
//            activePts.push_back(Q1);
//            activePts.push_back(Q2);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[2]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector9d partial_d2is_partial_x = pointEdgeDis2_grad(P1Coor, Q1Coor, Q2Coor); // 2
//            Matrix9d partial_2d2is_partial_x2 = pointEdgeDis2_hess(P1Coor, Q1Coor, Q2Coor); // 4
//
//            Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(P1Coor, Q1Coor, Q2Coor) * partial_b_partial_d2is;
//            Matrix9d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res_grad.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//
//
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = hess_(m * 3 + xd, n * 3 + yd);
//                                    res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//
//
//        case 2:
//        {
//            activePts.push_back(P1);
//            activePts.push_back(Q2);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector6d partial_d2is_partial_x = pointPointDis2_grad(P1Coor, Q2Coor); // 2
//            Matrix6d partial_2d2is_partial_x2 = pointPointDis2_hess(P1Coor, Q2Coor); // 4
//
//            Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P1Coor, Q2Coor) * partial_b_partial_d2is;
//            Matrix6d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res_grad.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//
//
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = hess_(m * 3 + xd, n * 3 + yd);
//                                    res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//
//        case 3:
//        {
//            activePts.push_back(Q1);
//            activePts.push_back(P1);
//            activePts.push_back(P2);
//            activePtsBC.push_back(vtInd_BC[2]);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[1]);
//
//            Vector9d partial_d2is_partial_x = pointEdgeDis2_grad(Q1Coor, P1Coor, P2Coor); // 2
//            Matrix9d partial_2d2is_partial_x2 = pointEdgeDis2_hess(Q1Coor, P1Coor, P2Coor); // 4
//
//            Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(Q1Coor, P1Coor, P2Coor) * partial_b_partial_d2is;
//            Matrix9d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res_grad.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//
//
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = hess_(m * 3 + xd, n * 3 + yd);
//                                    res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//
//
//        case 4:
//        {
//            activePts.push_back(P1);
//            activePts.push_back(P2);
//            activePts.push_back(Q1);
//            activePts.push_back(Q2);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[1]);
//            activePtsBC.push_back(vtInd_BC[2]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector12d partial_d2is_partial_x = edgeEdgeInsideDis2_grad(P1Coor, P2Coor, Q1Coor, Q2Coor); // 2
//            Matrix12d partial_2d2is_partial_x2 = edgeEdgeInsideDis2_hess(P1Coor, P2Coor, Q1Coor, Q2Coor); // 4
//
//
//            Vector12d grad_ = dt * dt * k_stiff * contactArea * edgeEdgeInsideDis2_grad(P1Coor, P2Coor, Q1Coor, Q2Coor) * partial_b_partial_d2is;
//            Matrix12d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res_grad.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//
//
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = hess_(m * 3 + xd, n * 3 + yd);
//                                    res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//
//
//        case 5:
//        {
//            activePts.push_back(Q2);
//            activePts.push_back(P1);
//            activePts.push_back(P2);
//            activePtsBC.push_back(vtInd_BC[3]);
//            activePtsBC.push_back(vtInd_BC[0]);
//            activePtsBC.push_back(vtInd_BC[1]);
//
//            Vector9d partial_d2is_partial_x = pointEdgeDis2_grad(Q2Coor, P1Coor, P2Coor); // 2
//            Matrix9d partial_2d2is_partial_x2 = pointEdgeDis2_hess(Q2Coor, P1Coor, P2Coor); // 4
//
//            Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(Q2Coor, P1Coor, P2Coor) * partial_b_partial_d2is;
//            Matrix9d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res_grad.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//
//
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = hess_(m * 3 + xd, n * 3 + yd);
//                                    res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//
//
//        case 6:
//        {
//            activePts.push_back(P2);
//            activePts.push_back(Q1);
//            activePtsBC.push_back(vtInd_BC[1]);
//            activePtsBC.push_back(vtInd_BC[2]);
//
//            Vector6d partial_d2is_partial_x = pointPointDis2_grad(P2Coor, Q1Coor); // 2
//            Matrix6d partial_2d2is_partial_x2 = pointPointDis2_hess(P2Coor, Q1Coor); // 4
//
//            Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P2Coor, Q1Coor) * partial_b_partial_d2is;
//            Matrix6d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res_grad.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//
//
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = hess_(m * 3 + xd, n * 3 + yd);
//                                    res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//
//
//        case 7:
//        {
//            activePts.push_back(P2);
//            activePts.push_back(Q1);
//            activePts.push_back(Q2);
//            activePtsBC.push_back(vtInd_BC[1]);
//            activePtsBC.push_back(vtInd_BC[2]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector9d partial_d2is_partial_x = pointEdgeDis2_grad(P2Coor, Q1Coor, Q2Coor); // 2
//            Matrix9d partial_2d2is_partial_x2 = pointEdgeDis2_hess(P2Coor, Q1Coor, Q2Coor); // 4
//
//            Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(P2Coor, Q1Coor, Q2Coor) * partial_b_partial_d2is;
//            Matrix9d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res_grad.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//
//
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = hess_(m * 3 + xd, n * 3 + yd);
//                                    res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//
//
//        case 8:
//        {
//            activePts.push_back(P2);
//            activePts.push_back(Q2);
//            activePtsBC.push_back(vtInd_BC[1]);
//            activePtsBC.push_back(vtInd_BC[3]);
//
//            Vector6d partial_d2is_partial_x = pointPointDis2_grad(P2Coor, Q2Coor); // 2
//            Matrix6d partial_2d2is_partial_x2 = pointPointDis2_hess(P2Coor, Q2Coor); // 4
//
//            Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P2Coor, Q2Coor) * partial_b_partial_d2is;
//            Matrix6d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
//            for (int m = 0; m < activePts.size(); m++)
//            {
//                if (activePtsBC[m] != 1)
//                {
//                    int x1_Ind = activePts[m]; // the first vertex index
//                    for (int xd = 0; xd < 3; xd++)
//                    {
//                        double value = grad_(m * 3 + xd, 1);
//                        res_grad.emplace_back(x1_Ind * 3 + xd, value);
//                    }
//
//
//                    for (int n = 0; n < activePts.size(); n++)
//                    {
//                        if (activePtsBC[n] != 1)
//                        {
//                            int x2_Ind = activePts[n]; // the second vertex index	
//                            for (int xd = 0; xd < 3; xd++)
//                            {
//                                for (int yd = 0; yd < 3; yd++)
//                                {
//                                    double value = hess_(m * 3 + xd, n * 3 + yd);
//                                    res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        break;
//
//
//        }
//    }
//
//    return std::make_pair(res_grad, res_hess);
//}




// compute the energy gradient and hessian 
std::pair<std::vector<std::pair<int, double>>, std::vector<Eigen::Triplet<double>>> BarrierEnergy::gradAndHess_PT(int type,
    double dis2, Mesh& tetMesh, Eigen::Vector4i& vtInd, Eigen::Vector4i& vtInd_BC, double d_hat, double k_stiff, double dt)
{
    std::vector<std::pair<int, double>> res_grad;
    std::vector<Eigen::Triplet<double>> res_hess;

    // the partial derivative of barrier energy b wrt distance d
    double d_hat2 = squaredDouble(d_hat);
    double partial_b_partial_d2is = -2.0 * (dis2 - d_hat2) * log(dis2 / d_hat2) - 1.0 / dis2 * squaredDouble(dis2 - d_hat2); // 3
    double partial_2b_partial_d2is2 = -2.0 * log(dis2 / d_hat2) - 4.0 / dis2 * (dis2 - d_hat2) + squaredDouble(dis2 - d_hat2) / dis2 / dis2; // 1                                                        

    int pt = vtInd[0], t1 = vtInd[1], t2 = vtInd[2], t3 = vtInd[3];
    double contactArea = tetMesh.boundaryVertices_area[pt];
    Eigen::Vector3d P = tetMesh.pos_node[pt], A = tetMesh.pos_node[t1], B = tetMesh.pos_node[t2], C = tetMesh.pos_node[t3];
    std::vector<int> activePts, activePtsBC;

    switch (type)
    {
    case 0:
    {
        activePts.push_back(pt);
        activePts.push_back(t1);
        activePtsBC.push_back(vtInd_BC[0]);
        activePtsBC.push_back(vtInd_BC[1]);

        Vector6d partial_d2is_partial_x = pointPointDis2_grad(P, A); // 2
        Matrix6d partial_2d2is_partial_x2 = pointPointDis2_hess(P, A); // 4

        Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P, A) * partial_b_partial_d2is;
        Matrix6d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
        store_grad_hess(res_grad, res_hess,grad_, hess_,activePts, activePtsBC);
    }
    break;
    case 1:
    {
        activePts.push_back(pt);
        activePts.push_back(t2);
        activePtsBC.push_back(vtInd_BC[0]);
        activePtsBC.push_back(vtInd_BC[2]);

        Vector6d partial_d2is_partial_x = pointPointDis2_grad(P, B); // 2
        Matrix6d partial_2d2is_partial_x2 = pointPointDis2_hess(P, B); // 4

        Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P, B) * partial_b_partial_d2is;
        Matrix6d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
        store_grad_hess(res_grad, res_hess, grad_, hess_, activePts, activePtsBC);
    }
    break;

    case 2:
    {
        activePts.push_back(pt);
        activePts.push_back(t3);
        activePtsBC.push_back(vtInd_BC[0]);
        activePtsBC.push_back(vtInd_BC[3]);

        Vector6d partial_d2is_partial_x = pointPointDis2_grad(P, C); // 2
        Matrix6d partial_2d2is_partial_x2 = pointPointDis2_hess(P, C); // 4

        Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P, C) * partial_b_partial_d2is;
        Matrix6d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
        store_grad_hess(res_grad, res_hess, grad_, hess_, activePts, activePtsBC);
    }
    break;

    case 3:
    {
        activePts.push_back(pt);
        activePts.push_back(t1);
        activePts.push_back(t2);
        activePtsBC.push_back(vtInd_BC[0]);
        activePtsBC.push_back(vtInd_BC[1]);
        activePtsBC.push_back(vtInd_BC[2]);

        Vector9d partial_d2is_partial_x = pointEdgeDis2_grad(P, A, B); // 2
        Matrix9d partial_2d2is_partial_x2 = pointEdgeDis2_hess(P, A, B); // 4

        Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(P, A, B) * partial_b_partial_d2is;
        Matrix9d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
        store_grad_hess(res_grad, res_hess, grad_, hess_, activePts, activePtsBC);
    }
    break;

    case 4:
    {
        activePts.push_back(pt);
        activePts.push_back(t2);
        activePts.push_back(t3);
        activePtsBC.push_back(vtInd_BC[0]);
        activePtsBC.push_back(vtInd_BC[2]);
        activePtsBC.push_back(vtInd_BC[3]);

        Vector9d partial_d2is_partial_x = pointEdgeDis2_grad(P, B, C); // 2
        Matrix9d partial_2d2is_partial_x2 = pointEdgeDis2_hess(P, B, C); // 4

        Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(P, B, C) * partial_b_partial_d2is;
        Matrix9d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
        store_grad_hess(res_grad, res_hess, grad_, hess_, activePts, activePtsBC);
    }
    break;

    case 5:
    {
        activePts.push_back(pt);
        activePts.push_back(t3);
        activePts.push_back(t1);
        activePtsBC.push_back(vtInd_BC[0]);
        activePtsBC.push_back(vtInd_BC[3]);
        activePtsBC.push_back(vtInd_BC[1]);

        Vector9d partial_d2is_partial_x = pointEdgeDis2_grad(P, C, A); // 2
        Matrix9d partial_2d2is_partial_x2 = pointEdgeDis2_hess(P, C, A); // 4

        Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(P, C, A) * partial_b_partial_d2is;
        Matrix9d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
        store_grad_hess(res_grad, res_hess, grad_, hess_, activePts, activePtsBC);
    }
    break;

    case 6:
    {
        activePts.push_back(pt);
        activePts.push_back(t1);
        activePts.push_back(t2);
        activePts.push_back(t3);
        activePtsBC.push_back(vtInd_BC[0]);
        activePtsBC.push_back(vtInd_BC[1]);
        activePtsBC.push_back(vtInd_BC[2]);
        activePtsBC.push_back(vtInd_BC[3]);

        Vector12d partial_d2is_partial_x = pointTriangleInsideDis2_grad(P, A, B, C); // 2
        Matrix12d partial_2d2is_partial_x2 = pointTriangleInsideDis2_hess(P, A, B, C); // 4

        Vector12d grad_ = dt * dt * k_stiff * contactArea * pointTriangleInsideDis2_grad(P, A, B, C) * partial_b_partial_d2is;
        Matrix12d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
        store_grad_hess(res_grad, res_hess, grad_, hess_, activePts, activePtsBC);
    }
    break;

    }

    return std::make_pair(res_grad, res_hess);
}

// compute the energy gradient and hessian 
std::pair<std::vector<std::pair<int, double>>, std::vector<Eigen::Triplet<double>>> BarrierEnergy::gradAndHess_EE(int type,
    double dis2, Mesh& tetMesh, Eigen::Vector4i& vtInd, Eigen::Vector4i& vtInd_BC, double d_hat, double k_stiff, double dt)
{
    std::vector<std::pair<int, double>> res_grad;
    std::vector<Eigen::Triplet<double>> res_hess;

    // the partial derivative of barrier energy b wrt distance d
    double d_hat2 = squaredDouble(d_hat);
    double partial_b_partial_d2is = -2.0 * (dis2 - d_hat2) * log(dis2 / d_hat2) - 1.0 / dis2 * squaredDouble(dis2 - d_hat2); // 3
    double partial_2b_partial_d2is2 = -2.0 * log(dis2 / d_hat2) - 4.0 / dis2 * (dis2 - d_hat2) + squaredDouble(dis2 - d_hat2) / dis2 / dis2; // 1                                                        


    int P1 = vtInd[0], P2 = vtInd[1], Q1 = vtInd[2], Q2 = vtInd[3];
    int emin = std::min(P1, P2), emax = std::max(P1, P2);
    double contactArea = tetMesh.boundaryEdges_area[emin][emax];
    Eigen::Vector3d P1Coor = tetMesh.pos_node[P1], P2Coor = tetMesh.pos_node[P2], Q1Coor = tetMesh.pos_node[Q1], Q2Coor = tetMesh.pos_node[Q2];
    std::vector<int> activePts, activePtsBC;

    switch (type)
    {
    case 0:
    {
        activePts.push_back(P1);
        activePts.push_back(Q1);
        activePtsBC.push_back(vtInd_BC[0]);
        activePtsBC.push_back(vtInd_BC[2]);

        Vector6d partial_d2is_partial_x = pointPointDis2_grad(P1Coor, Q1Coor); // 2
        Matrix6d partial_2d2is_partial_x2 = pointPointDis2_hess(P1Coor, Q1Coor); // 4

        Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P1Coor, Q1Coor) * partial_b_partial_d2is;
        Matrix6d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
        store_grad_hess(res_grad, res_hess, grad_, hess_, activePts, activePtsBC);
    }
    break;

    case 1:
    {
        activePts.push_back(P1);
        activePts.push_back(Q1);
        activePts.push_back(Q2);
        activePtsBC.push_back(vtInd_BC[0]);
        activePtsBC.push_back(vtInd_BC[2]);
        activePtsBC.push_back(vtInd_BC[3]);

        Vector9d partial_d2is_partial_x = pointEdgeDis2_grad(P1Coor, Q1Coor, Q2Coor); // 2
        Matrix9d partial_2d2is_partial_x2 = pointEdgeDis2_hess(P1Coor, Q1Coor, Q2Coor); // 4

        Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(P1Coor, Q1Coor, Q2Coor) * partial_b_partial_d2is;
        Matrix9d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
        store_grad_hess(res_grad, res_hess, grad_, hess_, activePts, activePtsBC);
    }
    break;


    case 2:
    {
        activePts.push_back(P1);
        activePts.push_back(Q2);
        activePtsBC.push_back(vtInd_BC[0]);
        activePtsBC.push_back(vtInd_BC[3]);

        Vector6d partial_d2is_partial_x = pointPointDis2_grad(P1Coor, Q2Coor); // 2
        Matrix6d partial_2d2is_partial_x2 = pointPointDis2_hess(P1Coor, Q2Coor); // 4

        Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P1Coor, Q2Coor) * partial_b_partial_d2is;
        Matrix6d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
        store_grad_hess(res_grad, res_hess, grad_, hess_, activePts, activePtsBC);
    }
    break;

    case 3:
    {
        activePts.push_back(Q1);
        activePts.push_back(P1);
        activePts.push_back(P2);
        activePtsBC.push_back(vtInd_BC[2]);
        activePtsBC.push_back(vtInd_BC[0]);
        activePtsBC.push_back(vtInd_BC[1]);

        Vector9d partial_d2is_partial_x = pointEdgeDis2_grad(Q1Coor, P1Coor, P2Coor); // 2
        Matrix9d partial_2d2is_partial_x2 = pointEdgeDis2_hess(Q1Coor, P1Coor, P2Coor); // 4

        Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(Q1Coor, P1Coor, P2Coor) * partial_b_partial_d2is;
        Matrix9d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
        store_grad_hess(res_grad, res_hess, grad_, hess_, activePts, activePtsBC);
    }
    break;


    case 4:
    {
        activePts.push_back(P1);
        activePts.push_back(P2);
        activePts.push_back(Q1);
        activePts.push_back(Q2);
        activePtsBC.push_back(vtInd_BC[0]);
        activePtsBC.push_back(vtInd_BC[1]);
        activePtsBC.push_back(vtInd_BC[2]);
        activePtsBC.push_back(vtInd_BC[3]);

        Vector12d partial_d2is_partial_x = edgeEdgeInsideDis2_grad(P1Coor, P2Coor, Q1Coor, Q2Coor); // 2
        Matrix12d partial_2d2is_partial_x2 = edgeEdgeInsideDis2_hess(P1Coor, P2Coor, Q1Coor, Q2Coor); // 4


        Vector12d grad_ = dt * dt * k_stiff * contactArea * edgeEdgeInsideDis2_grad(P1Coor, P2Coor, Q1Coor, Q2Coor) * partial_b_partial_d2is;
        Matrix12d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
        store_grad_hess(res_grad, res_hess, grad_, hess_, activePts, activePtsBC);
    }
    break;


    case 5:
    {
        activePts.push_back(Q2);
        activePts.push_back(P1);
        activePts.push_back(P2);
        activePtsBC.push_back(vtInd_BC[3]);
        activePtsBC.push_back(vtInd_BC[0]);
        activePtsBC.push_back(vtInd_BC[1]);

        Vector9d partial_d2is_partial_x = pointEdgeDis2_grad(Q2Coor, P1Coor, P2Coor); // 2
        Matrix9d partial_2d2is_partial_x2 = pointEdgeDis2_hess(Q2Coor, P1Coor, P2Coor); // 4

        Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(Q2Coor, P1Coor, P2Coor) * partial_b_partial_d2is;
        Matrix9d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
        store_grad_hess(res_grad, res_hess, grad_, hess_, activePts, activePtsBC);
    }
    break;


    case 6:
    {
        activePts.push_back(P2);
        activePts.push_back(Q1);
        activePtsBC.push_back(vtInd_BC[1]);
        activePtsBC.push_back(vtInd_BC[2]);

        Vector6d partial_d2is_partial_x = pointPointDis2_grad(P2Coor, Q1Coor); // 2
        Matrix6d partial_2d2is_partial_x2 = pointPointDis2_hess(P2Coor, Q1Coor); // 4

        Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P2Coor, Q1Coor) * partial_b_partial_d2is;
        Matrix6d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
        store_grad_hess(res_grad, res_hess, grad_, hess_, activePts, activePtsBC);
    }
    break;


    case 7:
    {
        activePts.push_back(P2);
        activePts.push_back(Q1);
        activePts.push_back(Q2);
        activePtsBC.push_back(vtInd_BC[1]);
        activePtsBC.push_back(vtInd_BC[2]);
        activePtsBC.push_back(vtInd_BC[3]);

        Vector9d partial_d2is_partial_x = pointEdgeDis2_grad(P2Coor, Q1Coor, Q2Coor); // 2
        Matrix9d partial_2d2is_partial_x2 = pointEdgeDis2_hess(P2Coor, Q1Coor, Q2Coor); // 4

        Vector9d grad_ = dt * dt * k_stiff * contactArea * pointEdgeDis2_grad(P2Coor, Q1Coor, Q2Coor) * partial_b_partial_d2is;
        Matrix9d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
        store_grad_hess(res_grad, res_hess, grad_, hess_, activePts, activePtsBC);
    }
    break;


    case 8:
    {
        activePts.push_back(P2);
        activePts.push_back(Q2);
        activePtsBC.push_back(vtInd_BC[1]);
        activePtsBC.push_back(vtInd_BC[3]);

        Vector6d partial_d2is_partial_x = pointPointDis2_grad(P2Coor, Q2Coor); // 2
        Matrix6d partial_2d2is_partial_x2 = pointPointDis2_hess(P2Coor, Q2Coor); // 4

        Vector6d grad_ = dt * dt * k_stiff * contactArea * pointPointDis2_grad(P2Coor, Q2Coor) * partial_b_partial_d2is;
        Matrix6d hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_d2is2 * partial_d2is_partial_x * partial_d2is_partial_x.transpose() + partial_b_partial_d2is * partial_2d2is_partial_x2);
        store_grad_hess(res_grad, res_hess, grad_, hess_, activePts, activePtsBC);
    }
    break;


    }
 
    return std::make_pair(res_grad, res_hess);
}



// store gradient and hessian value
void BarrierEnergy::store_grad_hess(std::vector<std::pair<int, double>>& res_grad, std::vector<Eigen::Triplet<double>>& res_hess,
    Vector6d& grad_, Matrix6d& hess_,
    std::vector<int>& activePts, std::vector<int>& activePtsBC)
{
    for (int m = 0; m < activePts.size(); m++)
    {
        if (activePtsBC[m] != 1)
        {
            int x1_Ind = activePts[m]; // the first vertex index
            for (int xd = 0; xd < 3; xd++)
            {
                double value = grad_(m * 3 + xd, 1);
                res_grad.emplace_back(x1_Ind * 3 + xd, value);
            }


            for (int n = 0; n < activePts.size(); n++)
            {
                if (activePtsBC[n] != 1)
                {
                    int x2_Ind = activePts[n]; // the second vertex index	
                    for (int xd = 0; xd < 3; xd++)
                    {
                        for (int yd = 0; yd < 3; yd++)
                        {
                            double value = hess_(m * 3 + xd, n * 3 + yd);
                            res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
                        }
                    }
                }
            }
        }
    }

}


// store gradient and hessian value
void BarrierEnergy::store_grad_hess(std::vector<std::pair<int, double>>& res_grad, std::vector<Eigen::Triplet<double>>& res_hess,
    Vector9d& grad_, Matrix9d& hess_,
    std::vector<int>& activePts, std::vector<int>& activePtsBC)
{
    for (int m = 0; m < activePts.size(); m++)
    {
        if (activePtsBC[m] != 1)
        {
            int x1_Ind = activePts[m]; // the first vertex index
            for (int xd = 0; xd < 3; xd++)
            {
                double value = grad_(m * 3 + xd, 1);
                res_grad.emplace_back(x1_Ind * 3 + xd, value);
            }


            for (int n = 0; n < activePts.size(); n++)
            {
                if (activePtsBC[n] != 1)
                {
                    int x2_Ind = activePts[n]; // the second vertex index	
                    for (int xd = 0; xd < 3; xd++)
                    {
                        for (int yd = 0; yd < 3; yd++)
                        {
                            double value = hess_(m * 3 + xd, n * 3 + yd);
                            res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
                        }
                    }
                }
            }
        }
    }

}


// store gradient and hessian value
void BarrierEnergy::store_grad_hess(std::vector<std::pair<int, double>>& res_grad, std::vector<Eigen::Triplet<double>>& res_hess,
    Vector12d& grad_, Matrix12d& hess_,
    std::vector<int>& activePts, std::vector<int>& activePtsBC)
{
    for (int m = 0; m < activePts.size(); m++)
    {
        if (activePtsBC[m] != 1)
        {
            int x1_Ind = activePts[m]; // the first vertex index
            for (int xd = 0; xd < 3; xd++)
            {
                double value = grad_(m * 3 + xd, 1);
                res_grad.emplace_back(x1_Ind * 3 + xd, value);
            }


            for (int n = 0; n < activePts.size(); n++)
            {
                if (activePtsBC[n] != 1)
                {
                    int x2_Ind = activePts[n]; // the second vertex index	
                    for (int xd = 0; xd < 3; xd++)
                    {
                        for (int yd = 0; yd < 3; yd++)
                        {
                            double value = hess_(m * 3 + xd, n * 3 + yd);
                            res_hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
                        }
                    }
                }
            }
        }
    }

}
