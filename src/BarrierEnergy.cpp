#include "BarrierEnergy.h"

// compute the energy gradient and hessian 
double BarrierEnergy::val_PT(double contactArea, double dis2, double d_hat, double k_stiff, double dt)
{
    return dt * dt * k_stiff * contactArea * compute_b(dis2, d_hat);
}

// compute the energy gradient and hessian 
double BarrierEnergy::val_EE(double contactArea, double dis2, Mesh& tetMesh, Eigen::Vector4i& ptIndices,  double d_hat, double k_stiff, double dt)
{
    // the partial derivative of barrier energy b wrt distance d
    double val_b = compute_b(dis2, d_hat);

    int P1 = ptIndices[0], P2 = ptIndices[1], Q1 = ptIndices[2], Q2 = ptIndices[3];
    int emin = std::min(P1, P2), emax = std::max(P1, P2);
    Eigen::Vector3d P1Coor = tetMesh.pos_node[P1], P2Coor = tetMesh.pos_node[P2], Q1Coor = tetMesh.pos_node[Q1], Q2Coor = tetMesh.pos_node[Q2];
    Eigen::Vector3d P1Coor_Rest = tetMesh.pos_node_Rest[P1], P2Coor_Rest = tetMesh.pos_node_Rest[P2], Q1Coor_Rest = tetMesh.pos_node_Rest[Q1], Q2Coor_Rest = tetMesh.pos_node_Rest[Q2];

    // the following codes calculate the mollifier-related value
    double eps_x = cal_EEM_eps_x(P1Coor_Rest, P2Coor_Rest, Q1Coor_Rest, Q2Coor_Rest);
    double val_ek = 0;
    DIS::compute_e(P1Coor, P2Coor, Q1Coor, Q2Coor, eps_x, val_ek);
    return dt * dt * k_stiff * contactArea * val_ek * val_b;
}


// compute the energy gradient and hessian 
void BarrierEnergy::gradAndHess_PT(BarrierEnergyRes& GH, Eigen::Vector4i& ptIndices, int type, double dis2, Mesh& tetMesh, double d_hat, double k_stiff, double dt)
{
    GH.PT_Indices.push_back(ptIndices);

    // the partial derivative of barrier energy b wrt distance d
    double g_bd = compute_g_b(dis2, d_hat); // 3
    double h_bd = compute_H_b(dis2, d_hat); // 1                                                        

    int pt = ptIndices[0], t1 = ptIndices[1], t2 = ptIndices[2], t3 = ptIndices[3];
    double contactArea = tetMesh.boundaryVertices_area[pt];

    Eigen::Vector3d P = tetMesh.pos_node[pt], A = tetMesh.pos_node[t1], B = tetMesh.pos_node[t2], C = tetMesh.pos_node[t3];

    switch (type)
    {
    case 0:
    {
        Vector6d g_dx = Vector6d::Zero();
        Matrix6d h_dx = Matrix6d::Zero();
        DIS::g_PP(P, A, g_dx); // 2
        DIS::H_PP(h_dx); // 4

        GH.V6.push_back(dt * dt * k_stiff * contactArea * g_dx * g_bd);

        Matrix6d hessian = dt * dt * k_stiff * contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
        makePD<double, 6>(hessian);
        GH.H6x6.push_back(hessian);
        Eigen::Vector2i activePtsLocalInd = { pt , t1 };
        GH.D2Index.push_back(activePtsLocalInd);

    }
    break;
    case 1:
    {
        Vector6d g_dx = Vector6d::Zero();
        Matrix6d h_dx = Matrix6d::Zero();
        DIS::g_PP(P, B, g_dx); // 2
        DIS::H_PP(h_dx); // 4

        GH.V6.push_back(dt * dt * k_stiff * contactArea * g_dx * g_bd);
        Matrix6d hessian = dt * dt * k_stiff * contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
        makePD<double, 6>(hessian);
        GH.H6x6.push_back(hessian);
        Eigen::Vector2i activePtsLocalInd = { pt , t2 };
        GH.D2Index.push_back(activePtsLocalInd);
    }
    break;

    case 2:
    {
        Vector6d g_dx = Vector6d::Zero();
        Matrix6d h_dx = Matrix6d::Zero();
        DIS::g_PP(P, C, g_dx); // 2
        DIS::H_PP(h_dx); // 4

        GH.V6.push_back(dt * dt * k_stiff * contactArea * g_dx * g_bd);
        Matrix6d hessian = dt * dt * k_stiff * contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
        makePD<double, 6>(hessian);
        GH.H6x6.push_back(hessian);
        Eigen::Vector2i activePtsLocalInd = { pt , t3 };
        GH.D2Index.push_back(activePtsLocalInd);
    }
    break;

    case 3:
    {
        Vector9d g_dx = Vector9d::Zero();
        Matrix9d h_dx = Matrix9d::Zero();
        DIS::g_PE(P, A, B, g_dx); // 2
        DIS::H_PE(P, A, B, h_dx); // 4

        GH.V9.push_back(dt * dt * k_stiff * contactArea * g_dx * g_bd);
        Matrix9d hessian = dt * dt * k_stiff * contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
        makePD<double, 9>(hessian);
        GH.H9x9.push_back(hessian);
        Eigen::Vector3i activePtsLocalInd = { pt , t1, t2 };
        GH.D3Index.push_back(activePtsLocalInd);
    }
    break;

    case 4:
    {
        Vector9d g_dx = Vector9d::Zero();
        Matrix9d h_dx = Matrix9d::Zero();
        DIS::g_PE(P, B, C, g_dx); // 2
        DIS::H_PE(P, B, C, h_dx); // 4

        GH.V9.push_back(dt * dt * k_stiff * contactArea * g_dx * g_bd);
        Matrix9d hessian = dt * dt * k_stiff * contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
        makePD<double, 9>(hessian);
        GH.H9x9.push_back(hessian);
        Eigen::Vector3i activePtsLocalInd = { pt , t2 , t3 };
        GH.D3Index.push_back(activePtsLocalInd);
    }
    break;

    case 5:
    {
        Vector9d g_dx = Vector9d::Zero();
        Matrix9d h_dx = Matrix9d::Zero();
        DIS::g_PE(P, C, A, g_dx); // 2
        DIS::H_PE(P, C, A, h_dx); // 4

        GH.V9.push_back(dt * dt * k_stiff * contactArea * g_dx * g_bd);
        Matrix9d hessian = dt * dt * k_stiff * contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
        makePD<double, 9>(hessian);
        GH.H9x9.push_back(hessian);
        Eigen::Vector3i activePtsLocalInd = { pt , t3 , t1 };
        GH.D3Index.push_back(activePtsLocalInd);
    }
    break;

    case 6:
    {
        Vector12d g_dx = Vector12d::Zero();
        Matrix12d h_dx = Matrix12d::Zero();
        DIS::g_PT(P, A, B, C, g_dx); // 2
        DIS::H_PT(P, A, B, C, h_dx); // 4

        GH.V12.push_back(dt * dt * k_stiff * contactArea * g_dx * g_bd);
        Matrix12d hessian = dt * dt * k_stiff * contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
        makePD<double, 12>(hessian);
        GH.H12x12.push_back(hessian);
        Eigen::Vector4i activePtsLocalInd = { pt , t1 , t2 , t3 };
        GH.D4Index.push_back(activePtsLocalInd);
    }
    break;

    }

}

// compute the energy gradient and hessian 
void BarrierEnergy::gradAndHess_EE(BarrierEnergyRes& GH, Eigen::Vector4i& ptIndices, int type, double dis2, Mesh& tetMesh, double d_hat, double k_stiff, double dt)
{
    GH.EE_Indices.push_back(ptIndices);

    // the partial derivative of barrier energy b wrt distance d
    double g_bd = compute_g_b(dis2, d_hat); // 3
    double h_bd = compute_H_b(dis2, d_hat); // 1     

    int P1 = ptIndices[0], P2 = ptIndices[1], Q1 = ptIndices[2], Q2 = ptIndices[3];
    int emin = std::min(P1, P2), emax = std::max(P1, P2);
    Eigen::Vector3d P1Coor = tetMesh.pos_node[P1], P2Coor = tetMesh.pos_node[P2], Q1Coor = tetMesh.pos_node[Q1], Q2Coor = tetMesh.pos_node[Q2];
    Eigen::Vector3d P1Coor_Rest = tetMesh.pos_node_Rest[P1], P2Coor_Rest = tetMesh.pos_node_Rest[P2], Q1Coor_Rest = tetMesh.pos_node_Rest[Q1], Q2Coor_Rest = tetMesh.pos_node_Rest[Q2];

    double contactArea = tetMesh.boundaryEdges_area[emin][emax];
    double val_b = compute_b(dis2, d_hat);
    Vector12d grad_b = Vector12d::Zero();
    Matrix12d hess_b = Matrix12d::Zero();

    // the following codes calculate the mollifier-related value
    double eps_x = cal_EEM_eps_x(P1Coor_Rest, P2Coor_Rest, Q1Coor_Rest, Q2Coor_Rest);
    double val_ek = 0;
    Vector12d grad_ek = Vector12d::Zero();
    Matrix12d hess_ek = Matrix12d::Zero();
    DIS::compute_e(P1Coor, P2Coor, Q1Coor, Q2Coor, eps_x, val_ek);
    DIS::compute_e_g(P1Coor, P2Coor, Q1Coor, Q2Coor, eps_x, grad_ek);
    DIS::compute_e_g(P1Coor, P2Coor, Q1Coor, Q2Coor, eps_x, grad_ek);
    DIS::compute_e_H(P1Coor, P2Coor, Q1Coor, Q2Coor, eps_x, hess_ek);



    switch (type)
    {
    case 0:
    {
        Vector6d g_dx = Vector6d::Zero(), grad_ = Vector6d::Zero();
        Matrix6d h_dx = Matrix6d::Zero(), hess_ = Matrix6d::Zero();
        DIS::g_PP(P1Coor, Q1Coor, g_dx); // 2
        DIS::H_PP(h_dx); // 4


        grad_ = contactArea * g_dx * g_bd;
        hess_ = contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
        Eigen::Vector2i activePtsLocalInd = { 0, 2 };
        project_grad_to_full(activePtsLocalInd, grad_, hess_, grad_b, hess_b); // since we add edge-edge mollifier, we have to project the hessian to full 12x12 matrix

        // now calcualte the final val, grad and hess considering the mollifier   
        GH.V12.push_back(dt * dt * k_stiff * contactArea * (grad_ek * val_b + val_ek * grad_b));
        Matrix12d hessian = dt * dt * k_stiff * contactArea * (hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b);
        makePD<double, 12>(hessian);
        GH.H12x12.push_back(hessian);
        GH.D4Index.push_back(ptIndices);
    }
    break;

    case 1:
    {
        Vector6d g_dx = Vector6d::Zero(), grad_ = Vector6d::Zero();
        Matrix6d h_dx = Matrix6d::Zero(), hess_ = Matrix6d::Zero();
        DIS::g_PP(P1Coor, Q2Coor, g_dx); // 2
        DIS::H_PP(h_dx); // 4


        grad_ = contactArea * g_dx * g_bd;
        hess_ = contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
        Eigen::Vector2i activePtsLocalInd = { 0, 3 };
        project_grad_to_full(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        GH.V12.push_back(dt * dt * k_stiff * contactArea * (grad_ek * val_b + val_ek * grad_b));
        Matrix12d hessian = dt * dt * k_stiff * contactArea * (hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b);
        makePD<double, 12>(hessian);
        GH.H12x12.push_back(hessian);
        GH.D4Index.push_back(ptIndices);
    }
    break;


    case 2:
    {
        Vector9d g_dx = Vector9d::Zero(), grad_ = Vector9d::Zero();
        Matrix9d h_dx = Matrix9d::Zero(), hess_ = Matrix9d::Zero();
        DIS::g_PE(P1Coor, Q1Coor, Q2Coor, g_dx); // 2
        DIS::H_PE(P1Coor, Q1Coor, Q2Coor, h_dx); // 4

        grad_ = contactArea * g_dx * g_bd;
        hess_ = contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
        Eigen::Vector3i activePtsLocalInd = { 0, 2, 3 };
        project_grad_to_full(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        GH.V12.push_back(dt * dt * k_stiff * contactArea * (grad_ek * val_b + val_ek * grad_b));
        Matrix12d hessian = dt * dt * k_stiff * contactArea * (hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b);
        makePD<double, 12>(hessian);
        GH.H12x12.push_back(hessian);
        GH.D4Index.push_back(ptIndices);
    }
    break;

    case 3:
    {
        Vector6d g_dx = Vector6d::Zero(), grad_ = Vector6d::Zero();
        Matrix6d h_dx = Matrix6d::Zero(), hess_ = Matrix6d::Zero();
        DIS::g_PP(P2Coor, Q1Coor, g_dx); // 2
        DIS::H_PP(h_dx); // 4



        grad_ = contactArea * g_dx * g_bd;
        hess_ = contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
        Eigen::Vector2i activePtsLocalInd = { 1, 2 };
        project_grad_to_full(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        GH.V12.push_back(dt * dt * k_stiff * contactArea * (grad_ek * val_b + val_ek * grad_b));
        Matrix12d hessian = dt * dt * k_stiff * contactArea * (hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b);
        makePD<double, 12>(hessian);
        GH.H12x12.push_back(hessian);
        GH.D4Index.push_back(ptIndices);
    }
    break;


    case 4:
    {
        Vector6d g_dx = Vector6d::Zero(), grad_ = Vector6d::Zero();
        Matrix6d h_dx = Matrix6d::Zero(), hess_ = Matrix6d::Zero();
        DIS::g_PP(P2Coor, Q2Coor, g_dx); // 2
        DIS::H_PP(h_dx); // 4


        grad_ = contactArea * g_dx * g_bd;
        hess_ = contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
        Eigen::Vector2i activePtsLocalInd = { 1, 3 };
        project_grad_to_full(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        GH.V12.push_back(dt* dt* k_stiff* contactArea* (grad_ek* val_b + val_ek * grad_b));
        Matrix12d hessian = dt * dt * k_stiff * contactArea * (hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b);
        makePD<double, 12>(hessian);
        GH.H12x12.push_back(hessian);
        GH.D4Index.push_back(ptIndices);
    }
    break;


    case 5:
    {
        Vector9d g_dx = Vector9d::Zero(), grad_ = Vector9d::Zero();
        Matrix9d h_dx = Matrix9d::Zero(), hess_ = Matrix9d::Zero();
        DIS::g_PE(P2Coor, Q1Coor, Q2Coor, g_dx); // 2
        DIS::H_PE(P2Coor, Q1Coor, Q2Coor, h_dx); // 4


        grad_ = contactArea * g_dx * g_bd;
        hess_ = contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
        Eigen::Vector3i activePtsLocalInd = { 1, 2, 3 };
        project_grad_to_full(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        GH.V12.push_back(dt* dt* k_stiff* contactArea* (grad_ek* val_b + val_ek * grad_b));
        Matrix12d hessian = dt * dt * k_stiff * contactArea * (hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b);
        makePD<double, 12>(hessian);
        GH.H12x12.push_back(hessian);
        GH.D4Index.push_back(ptIndices);
    }
    break;


    case 6:
    {
        Vector9d g_dx = Vector9d::Zero(), grad_ = Vector9d::Zero();
        Matrix9d h_dx = Matrix9d::Zero(), hess_ = Matrix9d::Zero();
        DIS::g_PE(Q1Coor, P1Coor, P2Coor, g_dx); // 2
        DIS::H_PE(Q1Coor, P1Coor, P2Coor, h_dx); // 4



        grad_ = contactArea * g_dx * g_bd;
        hess_ = contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
        Eigen::Vector3i activePtsLocalInd = { 2, 0, 1 };
        project_grad_to_full(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        GH.V12.push_back(dt* dt* k_stiff* contactArea* (grad_ek* val_b + val_ek * grad_b));
        Matrix12d hessian = dt * dt * k_stiff * contactArea * (hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b);
        makePD<double, 12>(hessian);
        GH.H12x12.push_back(hessian);
        GH.D4Index.push_back(ptIndices);
    }
    break;


    case 7:
    {
        Vector9d g_dx = Vector9d::Zero(), grad_ = Vector9d::Zero();
        Matrix9d h_dx = Matrix9d::Zero(), hess_ = Matrix9d::Zero();
        DIS::g_PE(Q2Coor, P1Coor, P2Coor, g_dx); // 2
        DIS::H_PE(Q2Coor, P1Coor, P2Coor, h_dx); // 4


        grad_ = contactArea * g_dx * g_bd;
        hess_ = contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
        Eigen::Vector3i activePtsLocalInd = { 3, 0, 1 };
        project_grad_to_full(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        GH.V12.push_back(dt* dt* k_stiff* contactArea* (grad_ek* val_b + val_ek * grad_b));
        Matrix12d hessian = dt * dt * k_stiff * contactArea * (hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b);
        makePD<double, 12>(hessian);
        GH.H12x12.push_back(hessian);
        GH.D4Index.push_back(ptIndices);
    }
    break;


    case 8:
    {
        Vector12d g_dx = Vector12d::Zero(), grad_b = Vector12d::Zero();
        Matrix12d h_dx = Matrix12d::Zero(), hess_b = Matrix12d::Zero();
        DIS::g_EE(P1Coor, P2Coor, Q1Coor, Q2Coor, g_dx); // 2
        DIS::H_EE(P1Coor, P2Coor, Q1Coor, Q2Coor, h_dx); // 4


        grad_b = contactArea * g_dx * g_bd;
        hess_b = contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);


        // now calcualte the final val, grad and hess considering the mollifier
        GH.V12.push_back(dt* dt* k_stiff* contactArea* (grad_ek* val_b + val_ek * grad_b));
        Matrix12d hessian = dt * dt * k_stiff * contactArea * (hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b);
        makePD<double, 12>(hessian);
        GH.H12x12.push_back(hessian);
        GH.D4Index.push_back(ptIndices);
    }
    break;


    }

}


void BarrierEnergy::project_grad_to_full(Eigen::Vector2i& activePtsLocalInd, Vector6d& grad_, Matrix6d& hess_, Vector12d& grad_full, Matrix12d& hess__full)
{   
    for (int i = 0; i < 2; i++)
    {
        int index_i = activePtsLocalInd[i];

        // project gradient
        for (int s = 0; s < 3; s++)
        {
            grad_full[index_i * 3 + s] = grad_[i * 3 + s];
        }

        // project hessian
        for (int j = 0; j < 2; j++)
        {
            int index_j = activePtsLocalInd[j];

            hess__full.block(index_i * 3, index_j * 3, 3, 3) = hess_.block(i * 3, j * 3, 3, 3);
        }

    }
    
}

void BarrierEnergy::project_grad_to_full(Eigen::Vector3i& activePtsLocalInd, Vector9d& grad_, Matrix9d& hess_, Vector12d& grad_full, Matrix12d& hess__full)
{
    for (int i = 0; i < 3; i++)
    {
        int index_i = activePtsLocalInd[i];

        // project gradient
        for (int s = 0; s < 3; s++)
        {
            grad_full[index_i * 3 + s] = grad_[i * 3 + s];
        }

        // project hessian
        for (int j = 0; j < 3; j++)
        {
            int index_j = activePtsLocalInd[j];

            hess__full.block(index_i * 3, index_j * 3, 3, 3) = hess_.block(i * 3, j * 3, 3, 3);
        }

    }

}



double BarrierEnergy::compute_b(double d, double dHat)
{
    return -(d - dHat) * (d - dHat) * log(d / dHat);
}

double BarrierEnergy::compute_g_b(double d, double dHat)
{
    double t = d - dHat;
    return t * std::log(d / dHat) * -2.0 - (t * t) / d;
}

double BarrierEnergy::compute_H_b(double d, double dHat)
{
    double t = d - dHat;
    return (std::log(d / dHat) * -2.0 - t * 4.0 / d) + 1.0 / (d * d) * (t * t);
}



double Ground::val(double coor_z2, double d_hat, double distributedArea, double k_stiff, double dt)
{
    return dt * dt * k_stiff * distributedArea * BarrierEnergy::compute_b(coor_z2, d_hat);
}


void Ground::gradAndHess(BarrierEnergyRes& BEres, int index_i,  double coor_z2, double d_hat, double distributedArea, double k_stiff, double dt)
{
    double g_bd = BarrierEnergy::compute_g_b(coor_z2, d_hat); // 3
    double h_bd = BarrierEnergy::compute_H_b(coor_z2, d_hat); // 1                                                        
    
    Eigen::Vector3d g_dx = {0, 0, 2.0 * std::sqrt(coor_z2) };
    Eigen::Matrix3d h_dx = Eigen::Matrix3d::Zero();
    h_dx(2, 2) = 2.0;

    BEres.V3.push_back(dt * dt * k_stiff * distributedArea * g_dx * g_bd);

    Matrix3d hessian = dt * dt * k_stiff * distributedArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
    makePD<double, 3>(hessian);
    BEres.H3x3.push_back(hessian);
    BEres.D1Index.push_back(index_i);

}