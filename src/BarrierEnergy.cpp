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
    double partial_b_partial_dis = compute_g_b(dis2, d_hat); // 3
    double partial_2b_partial_dis2 = compute_H_b(dis2, d_hat); // 1                                                        

    int pt = ptIndices[0], t1 = ptIndices[1], t2 = ptIndices[2], t3 = ptIndices[3];
    double contactArea = tetMesh.boundaryVertices_area[pt];

    Eigen::Vector3d P = tetMesh.pos_node[pt], A = tetMesh.pos_node[t1], B = tetMesh.pos_node[t2], C = tetMesh.pos_node[t3];

    switch (type)
    {
    case 0:
    {
        Vector6d partial_dis_partial_x = Vector6d::Zero();
        Matrix6d partial_2dis_partial_x2 = Matrix6d::Zero();
        DIS::g_PP(P, A, partial_dis_partial_x); // 2
        DIS::H_PP(partial_2dis_partial_x2); // 4

        GH.V6.push_back(dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis);

        Matrix6d hessian = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        makePD<double, 6>(hessian);
        GH.H6x6.push_back(hessian);
        Eigen::Vector2i activePtsLocalInd = { pt , t1 };
        GH.D2Index.push_back(activePtsLocalInd);

    }
    break;
    case 1:
    {
        Vector6d partial_dis_partial_x = Vector6d::Zero();
        Matrix6d partial_2dis_partial_x2 = Matrix6d::Zero();
        DIS::g_PP(P, B, partial_dis_partial_x); // 2
        DIS::H_PP(partial_2dis_partial_x2); // 4

        GH.V6.push_back(dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis);
        Matrix6d hessian = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        makePD<double, 6>(hessian);
        GH.H6x6.push_back(hessian);
        Eigen::Vector2i activePtsLocalInd = { pt , t2 };
        GH.D2Index.push_back(activePtsLocalInd);
    }
    break;

    case 2:
    {
        Vector6d partial_dis_partial_x = Vector6d::Zero();
        Matrix6d partial_2dis_partial_x2 = Matrix6d::Zero();
        DIS::g_PP(P, C, partial_dis_partial_x); // 2
        DIS::H_PP(partial_2dis_partial_x2); // 4

        GH.V6.push_back(dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis);
        Matrix6d hessian = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        makePD<double, 6>(hessian);
        GH.H6x6.push_back(hessian);
        Eigen::Vector2i activePtsLocalInd = { pt , t3 };
        GH.D2Index.push_back(activePtsLocalInd);
    }
    break;

    case 3:
    {
        Vector9d partial_dis_partial_x = Vector9d::Zero();
        Matrix9d partial_2dis_partial_x2 = Matrix9d::Zero();
        DIS::g_PE(P, A, B, partial_dis_partial_x); // 2
        DIS::H_PE(P, A, B, partial_2dis_partial_x2); // 4

        GH.V9.push_back(dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis);
        Matrix9d hessian = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        makePD<double, 9>(hessian);
        GH.H9x9.push_back(hessian);
        Eigen::Vector3i activePtsLocalInd = { pt , t1, t2 };
        GH.D3Index.push_back(activePtsLocalInd);
    }
    break;

    case 4:
    {
        Vector9d partial_dis_partial_x = Vector9d::Zero();
        Matrix9d partial_2dis_partial_x2 = Matrix9d::Zero();
        DIS::g_PE(P, B, C, partial_dis_partial_x); // 2
        DIS::H_PE(P, B, C, partial_2dis_partial_x2); // 4

        GH.V9.push_back(dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis);
        Matrix9d hessian = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        makePD<double, 9>(hessian);
        GH.H9x9.push_back(hessian);
        Eigen::Vector3i activePtsLocalInd = { pt , t2 , t3 };
        GH.D3Index.push_back(activePtsLocalInd);
    }
    break;

    case 5:
    {
        Vector9d partial_dis_partial_x = Vector9d::Zero();
        Matrix9d partial_2dis_partial_x2 = Matrix9d::Zero();
        DIS::g_PE(P, C, A, partial_dis_partial_x); // 2
        DIS::H_PE(P, C, A, partial_2dis_partial_x2); // 4

        GH.V9.push_back(dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis);
        Matrix9d hessian = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        makePD<double, 9>(hessian);
        GH.H9x9.push_back(hessian);
        Eigen::Vector3i activePtsLocalInd = { pt , t3 , t1 };
        GH.D3Index.push_back(activePtsLocalInd);
    }
    break;

    case 6:
    {
        Vector12d partial_dis_partial_x = Vector12d::Zero();
        Matrix12d partial_2dis_partial_x2 = Matrix12d::Zero();
        DIS::g_PT(P, A, B, C, partial_dis_partial_x); // 2
        DIS::H_PT(P, A, B, C, partial_2dis_partial_x2); // 4

        GH.V12.push_back(dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis);
        Matrix12d hessian = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
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
    double partial_b_partial_dis = compute_g_b(dis2, d_hat); // 3
    double partial_2b_partial_dis2 = compute_H_b(dis2, d_hat); // 1     

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
        Vector6d partial_dis_partial_x = Vector6d::Zero(), grad_ = Vector6d::Zero();
        Matrix6d partial_2dis_partial_x2 = Matrix6d::Zero(), hess_ = Matrix6d::Zero();
        DIS::g_PP(P1Coor, Q1Coor, partial_dis_partial_x); // 2
        DIS::H_PP(partial_2dis_partial_x2); // 4


        grad_ = contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
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
        Vector6d partial_dis_partial_x = Vector6d::Zero(), grad_ = Vector6d::Zero();
        Matrix6d partial_2dis_partial_x2 = Matrix6d::Zero(), hess_ = Matrix6d::Zero();
        DIS::g_PP(P1Coor, Q2Coor, partial_dis_partial_x); // 2
        DIS::H_PP(partial_2dis_partial_x2); // 4


        grad_ = contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
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
        Vector9d partial_dis_partial_x = Vector9d::Zero(), grad_ = Vector9d::Zero();
        Matrix9d partial_2dis_partial_x2 = Matrix9d::Zero(), hess_ = Matrix9d::Zero();
        DIS::g_PE(P1Coor, Q1Coor, Q2Coor, partial_dis_partial_x); // 2
        DIS::H_PE(P1Coor, Q1Coor, Q2Coor, partial_2dis_partial_x2); // 4

        grad_ = contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
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
        Vector6d partial_dis_partial_x = Vector6d::Zero(), grad_ = Vector6d::Zero();
        Matrix6d partial_2dis_partial_x2 = Matrix6d::Zero(), hess_ = Matrix6d::Zero();
        DIS::g_PP(P2Coor, Q1Coor, partial_dis_partial_x); // 2
        DIS::H_PP(partial_2dis_partial_x2); // 4



        grad_ = contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
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
        Vector6d partial_dis_partial_x = Vector6d::Zero(), grad_ = Vector6d::Zero();
        Matrix6d partial_2dis_partial_x2 = Matrix6d::Zero(), hess_ = Matrix6d::Zero();
        DIS::g_PP(P2Coor, Q2Coor, partial_dis_partial_x); // 2
        DIS::H_PP(partial_2dis_partial_x2); // 4


        grad_ = contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
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
        Vector9d partial_dis_partial_x = Vector9d::Zero(), grad_ = Vector9d::Zero();
        Matrix9d partial_2dis_partial_x2 = Matrix9d::Zero(), hess_ = Matrix9d::Zero();
        DIS::g_PE(P2Coor, Q1Coor, Q2Coor, partial_dis_partial_x); // 2
        DIS::H_PE(P2Coor, Q1Coor, Q2Coor, partial_2dis_partial_x2); // 4


        grad_ = contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
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
        Vector9d partial_dis_partial_x = Vector9d::Zero(), grad_ = Vector9d::Zero();
        Matrix9d partial_2dis_partial_x2 = Matrix9d::Zero(), hess_ = Matrix9d::Zero();
        DIS::g_PE(Q1Coor, P1Coor, P2Coor, partial_dis_partial_x); // 2
        DIS::H_PE(Q1Coor, P1Coor, P2Coor, partial_2dis_partial_x2); // 4



        grad_ = contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
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
        Vector9d partial_dis_partial_x = Vector9d::Zero(), grad_ = Vector9d::Zero();
        Matrix9d partial_2dis_partial_x2 = Matrix9d::Zero(), hess_ = Matrix9d::Zero();
        DIS::g_PE(Q2Coor, P1Coor, P2Coor, partial_dis_partial_x); // 2
        DIS::H_PE(Q2Coor, P1Coor, P2Coor, partial_2dis_partial_x2); // 4


        grad_ = contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
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
        Vector12d partial_dis_partial_x = Vector12d::Zero(), grad_b = Vector12d::Zero();
        Matrix12d partial_2dis_partial_x2 = Matrix12d::Zero(), hess_b = Matrix12d::Zero();
        DIS::g_EE(P1Coor, P2Coor, Q1Coor, Q2Coor, partial_dis_partial_x); // 2
        DIS::H_EE(P1Coor, P2Coor, Q1Coor, Q2Coor, partial_2dis_partial_x2); // 4


        grad_b = contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_b = contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);


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



double Ground::val(double coor_z, double d_hat, double distributedArea, double k_stiff, double dt)
{
    double dis2 = coor_z * coor_z;
    double d_hat2 = squaredDouble(d_hat);
    double bEnergy = -squaredDouble(dis2 - d_hat * d_hat) * log(dis2 / (d_hat * d_hat));
    return dt * dt * k_stiff * distributedArea * bEnergy;
}


void Ground::gradAndHess(BarrierEnergyRes& BEres, int index_i,  double coor_z, double d_hat, double distributedArea, double k_stiff, double dt)
{
    double dis2 = coor_z * coor_z;
    double d_hat2 = squaredDouble(d_hat);

    Eigen::Vector3d grad_b = Eigen::Vector3d::Zero();
    Eigen::Matrix3d hess_b = Eigen::Matrix3d::Zero();

    double partial_b_partial_dis = (-2.0 * (dis2 - d_hat2) * log(dis2 / d_hat2) - 1.0 / dis2 * squaredDouble(dis2 - d_hat2)) * 2.0 * std::sqrt(dis2); // 3
    double partial_2b_partial_dis2 = -4.0 * (3.0 * dis2 - d_hat2) * log(dis2 / d_hat2) - (4.0 / std::sqrt(dis2) + 8.0) * (dis2 - d_hat2) + 2.0 / dis2 * squaredDouble(dis2 - d_hat2); // 1                                                        
    
    Eigen::Vector3d partial_dis_partial_x = {0, 0, 1.0 };
    Eigen::Matrix3d partial_2dis_partial_x2 = Eigen::Matrix3d::Zero();

    BEres.V3.push_back(dt * dt * k_stiff * distributedArea * partial_dis_partial_x * partial_b_partial_dis);
    BEres.H3x3.push_back(dt * dt * k_stiff * distributedArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2));
    BEres.D1Index.push_back(index_i);
}