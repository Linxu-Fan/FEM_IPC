#include "BarrierEnergy.h"

// compute the energy gradient and hessian 
void BarrierEnergy::valGradAndHess_PT(BarrierEnergyRes& BEres, int type, double dis2, Mesh& tetMesh, double d_hat, double k_stiff, double dt)
{

    // the partial derivative of barrier energy b wrt distance d
    double d_hat2 = squaredDouble(d_hat);
    double partial_b_partial_dis = (- 2.0 * (dis2 - d_hat2) * log(dis2 / d_hat2) - 1.0 / dis2 * squaredDouble(dis2 - d_hat2)) * 2.0 * std::sqrt(dis2); // 3
    double partial_2b_partial_dis2 = -4.0 * (3.0 * dis2 - d_hat2) * log(dis2 / d_hat2) - (4.0 / std::sqrt(dis2) + 8.0) * (dis2 - d_hat2) + 2.0 / dis2 * squaredDouble(dis2 - d_hat2); // 1                                                        

    int pt = BEres.PP_Index[0], t1 = BEres.PP_Index[1], t2 = BEres.PP_Index[2], t3 = BEres.PP_Index[3];
    double contactArea = tetMesh.boundaryVertices_area[pt];
    double bEnergy = -squaredDouble(dis2 - d_hat * d_hat) * log(dis2 / (d_hat * d_hat));


    BEres.Val = dt * dt * k_stiff * contactArea * bEnergy;


    Eigen::Vector3d P = tetMesh.pos_node[pt], A = tetMesh.pos_node[t1], B = tetMesh.pos_node[t2], C = tetMesh.pos_node[t3];

    switch (type)
    {
    case 0:
    {
        Vector6d partial_dis_partial_x = Vector6d::Zero(), grad_ = Vector6d::Zero();
        Matrix6d partial_2dis_partial_x2 = Matrix6d::Zero(), hess_ = Matrix6d::Zero();
        DIS::g_PP(P, A, partial_dis_partial_x); // 2
        DIS::H_PP(partial_2dis_partial_x2); // 4

        grad_ = dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        Eigen::Vector2i activePtsLocalInd = { 0 , 1 };
        store_grad_hess_PT(BEres, grad_, hess_, activePtsLocalInd);
    }
    break;
    case 1:
    {
        Vector6d partial_dis_partial_x = Vector6d::Zero(), grad_ = Vector6d::Zero();
        Matrix6d partial_2dis_partial_x2 = Matrix6d::Zero(), hess_ = Matrix6d::Zero();
        DIS::g_PP(P, B, partial_dis_partial_x); // 2
        DIS::H_PP(partial_2dis_partial_x2); // 4

        grad_ = dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        Eigen::Vector2i activePtsLocalInd = { 0 , 2 };
        store_grad_hess_PT(BEres, grad_, hess_, activePtsLocalInd);
    }
    break;

    case 2:
    {
        Vector6d partial_dis_partial_x = Vector6d::Zero(), grad_ = Vector6d::Zero();
        Matrix6d partial_2dis_partial_x2 = Matrix6d::Zero(), hess_ = Matrix6d::Zero();
        DIS::g_PP(P, C, partial_dis_partial_x); // 2
        DIS::H_PP(partial_2dis_partial_x2); // 4

        grad_ = dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        Eigen::Vector2i activePtsLocalInd = { 0 , 3 };
        store_grad_hess_PT(BEres, grad_, hess_, activePtsLocalInd);
    }
    break;

    case 3:
    {
        Vector9d partial_dis_partial_x = Vector9d::Zero(), grad_ = Vector9d::Zero();
        Matrix9d partial_2dis_partial_x2 = Matrix9d::Zero(), hess_ = Matrix9d::Zero();
        DIS::g_PE(P, A, B, partial_dis_partial_x); // 2
        DIS::H_PE(P, A, B, partial_2dis_partial_x2); // 4

        grad_ = dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        Eigen::Vector3i activePtsLocalInd = { 0 , 1, 2 };
        store_grad_hess_PT(BEres, grad_, hess_, activePtsLocalInd);
    }
    break;

    case 4:
    {
        Vector9d partial_dis_partial_x = Vector9d::Zero(), grad_ = Vector9d::Zero();
        Matrix9d partial_2dis_partial_x2 = Matrix9d::Zero(), hess_ = Matrix9d::Zero();
        DIS::g_PE(P, B, C, partial_dis_partial_x); // 2
        DIS::H_PE(P, B, C, partial_2dis_partial_x2); // 4

        grad_ = dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        Eigen::Vector3i activePtsLocalInd = { 0 , 2 , 3 };
        store_grad_hess_PT(BEres, grad_, hess_, activePtsLocalInd);
    }
    break;

    case 5:
    {
        Vector9d partial_dis_partial_x = Vector9d::Zero(), grad_ = Vector9d::Zero();
        Matrix9d partial_2dis_partial_x2 = Matrix9d::Zero(), hess_ = Matrix9d::Zero();
        DIS::g_PE(P, C, A, partial_dis_partial_x); // 2
        DIS::H_PE(P, C, A, partial_2dis_partial_x2); // 4

        grad_ = dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        Eigen::Vector3i activePtsLocalInd = { 0 , 3 , 1 };
        store_grad_hess_PT(BEres, grad_, hess_, activePtsLocalInd);
    }
    break;

    case 6:
    {
        Vector12d partial_dis_partial_x = Vector12d::Zero(), grad_ = Vector12d::Zero();
        Matrix12d partial_2dis_partial_x2 = Matrix12d::Zero(), hess_ = Matrix12d::Zero();
        DIS::g_PT(P, A, B, C, partial_dis_partial_x); // 2
        DIS::H_PT(P, A, B, C, partial_2dis_partial_x2); // 4

        grad_ = dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        Eigen::Vector4i activePtsLocalInd = { 0 , 1 , 2 , 3 };
        store_grad_hess_PT(BEres, grad_, hess_, activePtsLocalInd);
    }
    break;

    }

}

// compute the energy gradient and hessian 
void BarrierEnergy::valGradAndHess_EE(BarrierEnergyRes& BEres, int type, double dis2, Mesh& tetMesh, double d_hat, double k_stiff, double dt)
{

    // the partial derivative of barrier energy b wrt distance d
    double d_hat2 = squaredDouble(d_hat);
    double partial_b_partial_dis = (-2.0 * (dis2 - d_hat2) * log(dis2 / d_hat2) - 1.0 / dis2 * squaredDouble(dis2 - d_hat2)) * 2.0 * std::sqrt(dis2); // 3
    double partial_2b_partial_dis2 = -4.0 * (3.0 * dis2 - d_hat2) * log(dis2 / d_hat2) - (4.0 / std::sqrt(dis2) + 8.0) * (dis2 - d_hat2) + 2.0 / dis2 * squaredDouble(dis2 - d_hat2); // 1                                                        


    int P1 = BEres.PP_Index[0], P2 = BEres.PP_Index[1], Q1 = BEres.PP_Index[2], Q2 = BEres.PP_Index[3];
    int emin = std::min(P1, P2), emax = std::max(P1, P2);
    Eigen::Vector3d P1Coor = tetMesh.pos_node[P1], P2Coor = tetMesh.pos_node[P2], Q1Coor = tetMesh.pos_node[Q1], Q2Coor = tetMesh.pos_node[Q2];
    Eigen::Vector3d P1Coor_Rest = tetMesh.pos_node_Rest[P1], P2Coor_Rest = tetMesh.pos_node_Rest[P2], Q1Coor_Rest = tetMesh.pos_node_Rest[Q1], Q2Coor_Rest = tetMesh.pos_node_Rest[Q2];

    double contactArea = tetMesh.boundaryEdges_area[emin][emax];
    double bEnergy = -squaredDouble(dis2 - d_hat * d_hat) * log(dis2 / (d_hat * d_hat));
    double val_b = dt * dt * k_stiff * contactArea * bEnergy;
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


        grad_ = dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        Eigen::Vector2i activePtsLocalInd = { 0, 2 };
        project_grad_to_full(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        BEres.Val = val_ek * val_b;
        Vector12d grad_final = grad_ek * val_b + val_ek * grad_b;
        Matrix12d hess_final = hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b;
        store_grad_hess_EE(BEres, grad_final, hess_final);
    }
    break;

    case 1:
    {
        Vector6d partial_dis_partial_x = Vector6d::Zero(), grad_ = Vector6d::Zero();
        Matrix6d partial_2dis_partial_x2 = Matrix6d::Zero(), hess_ = Matrix6d::Zero();
        DIS::g_PP(P1Coor, Q2Coor, partial_dis_partial_x); // 2
        DIS::H_PP(partial_2dis_partial_x2); // 4


        grad_ = dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        Eigen::Vector2i activePtsLocalInd = { 0, 3 };
        project_grad_to_full(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        BEres.Val = val_ek * val_b;
        Vector12d grad_final = grad_ek * val_b + val_ek * grad_b;
        Matrix12d hess_final = hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b;
        store_grad_hess_EE(BEres, grad_final, hess_final);
    }
    break;


    case 2:
    {
        Vector9d partial_dis_partial_x = Vector9d::Zero(), grad_ = Vector9d::Zero();
        Matrix9d partial_2dis_partial_x2 = Matrix9d::Zero(), hess_ = Matrix9d::Zero();
        DIS::g_PE(P1Coor, Q1Coor, Q2Coor, partial_dis_partial_x); // 2
        DIS::H_PE(P1Coor, Q1Coor, Q2Coor, partial_2dis_partial_x2); // 4

        grad_ = dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        Eigen::Vector3i activePtsLocalInd = { 0, 2, 3 };
        project_grad_to_full(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        BEres.Val = val_ek * val_b;
        Vector12d grad_final = grad_ek * val_b + val_ek * grad_b;
        Matrix12d hess_final = hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b;
        store_grad_hess_EE(BEres, grad_final, hess_final);
    }
    break;

    case 3:
    {
        Vector6d partial_dis_partial_x = Vector6d::Zero(), grad_ = Vector6d::Zero();
        Matrix6d partial_2dis_partial_x2 = Matrix6d::Zero(), hess_ = Matrix6d::Zero();
        DIS::g_PP(P2Coor, Q1Coor, partial_dis_partial_x); // 2
        DIS::H_PP(partial_2dis_partial_x2); // 4



        grad_ = dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        Eigen::Vector2i activePtsLocalInd = { 1, 2 };
        project_grad_to_full(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        BEres.Val = val_ek * val_b;
        Vector12d grad_final = grad_ek * val_b + val_ek * grad_b;
        Matrix12d hess_final = hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b;
        store_grad_hess_EE(BEres, grad_final, hess_final);
    }
    break;


    case 4:
    {
        Vector6d partial_dis_partial_x = Vector6d::Zero(), grad_ = Vector6d::Zero();
        Matrix6d partial_2dis_partial_x2 = Matrix6d::Zero(), hess_ = Matrix6d::Zero();
        DIS::g_PP(P2Coor, Q2Coor, partial_dis_partial_x); // 2
        DIS::H_PP(partial_2dis_partial_x2); // 4


        grad_ = dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        Eigen::Vector2i activePtsLocalInd = { 1, 3 };
        project_grad_to_full(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        BEres.Val = val_ek * val_b;
        Vector12d grad_final = grad_ek * val_b + val_ek * grad_b;
        Matrix12d hess_final = hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b;
        store_grad_hess_EE(BEres, grad_final, hess_final);
    }
    break;


    case 5:
    {
        Vector9d partial_dis_partial_x = Vector9d::Zero(), grad_ = Vector9d::Zero();
        Matrix9d partial_2dis_partial_x2 = Matrix9d::Zero(), hess_ = Matrix9d::Zero();
        DIS::g_PE(P2Coor, Q1Coor, Q2Coor, partial_dis_partial_x); // 2
        DIS::H_PE(P2Coor, Q1Coor, Q2Coor, partial_2dis_partial_x2); // 4


        grad_ = dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        Eigen::Vector3i activePtsLocalInd = { 1, 2, 3 };
        project_grad_to_full(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        BEres.Val = val_ek * val_b;
        Vector12d grad_final = grad_ek * val_b + val_ek * grad_b;
        Matrix12d hess_final = hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b;
        store_grad_hess_EE(BEres, grad_final, hess_final);
    }
    break;


    case 6:
    {
        Vector9d partial_dis_partial_x = Vector9d::Zero(), grad_ = Vector9d::Zero();
        Matrix9d partial_2dis_partial_x2 = Matrix9d::Zero(), hess_ = Matrix9d::Zero();
        DIS::g_PE(Q1Coor, P1Coor, P2Coor, partial_dis_partial_x); // 2
        DIS::H_PE(Q1Coor, P1Coor, P2Coor, partial_2dis_partial_x2); // 4



        grad_ = dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        Eigen::Vector3i activePtsLocalInd = { 2, 0, 1 };
        project_grad_to_full(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        BEres.Val = val_ek * val_b;
        Vector12d grad_final = grad_ek * val_b + val_ek * grad_b;
        Matrix12d hess_final = hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b;
        store_grad_hess_EE(BEres, grad_final, hess_final);
    }
    break;


    case 7:
    {
        Vector9d partial_dis_partial_x = Vector9d::Zero(), grad_ = Vector9d::Zero();
        Matrix9d partial_2dis_partial_x2 = Matrix9d::Zero(), hess_ = Matrix9d::Zero();
        DIS::g_PE(Q2Coor, P1Coor, P2Coor, partial_dis_partial_x); // 2
        DIS::H_PE(Q2Coor, P1Coor, P2Coor, partial_2dis_partial_x2); // 4


        grad_ = dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        Eigen::Vector3i activePtsLocalInd = { 3, 0, 1 };
        project_grad_to_full(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        BEres.Val = val_ek * val_b;
        Vector12d grad_final = grad_ek * val_b + val_ek * grad_b;
        Matrix12d hess_final = hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b;
        store_grad_hess_EE(BEres, grad_final, hess_final);
    }
    break;


    case 8:
    {
        Vector12d partial_dis_partial_x = Vector12d::Zero(), grad_ = Vector12d::Zero();
        Matrix12d partial_2dis_partial_x2 = Matrix12d::Zero(), hess_ = Matrix12d::Zero();
        DIS::g_EE(P1Coor, P2Coor, Q1Coor, Q2Coor, partial_dis_partial_x); // 2
        DIS::H_EE(P1Coor, P2Coor, Q1Coor, Q2Coor, partial_2dis_partial_x2); // 4


        grad_ = dt * dt * k_stiff * contactArea * partial_dis_partial_x * partial_b_partial_dis;
        hess_ = dt * dt * k_stiff * contactArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);
        grad_b = grad_;
        hess_b = hess_;

        // now calcualte the final val, grad and hess considering the mollifier
        BEres.Val = val_ek * val_b;
        Vector12d grad_final = grad_ek * val_b + val_ek * grad_b;
        Matrix12d hess_final = hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b;
        store_grad_hess_EE(BEres, grad_final, hess_final);
    }
    break;


    }

}


// store gradient and hessian value
void BarrierEnergy::store_grad_hess_PT(BarrierEnergyRes& BEres, Vector6d& grad_, Matrix6d& hess_, Eigen::Vector2i& activePtsLocalInd)
{
    for (int m = 0; m < activePtsLocalInd.size(); m++)
    {
        int index_m = activePtsLocalInd[m];
        if (BEres.vtInd_BC[index_m] != 1)
        {
            int x1_Ind = BEres.PP_Index[index_m]; // the first vertex index
            for (int xd = 0; xd < 3; xd++)
            {
                double value = grad_(m * 3 + xd, 1);
                BEres.Grad.emplace_back(x1_Ind * 3 + xd, value);
            }


            for (int n = 0; n < activePtsLocalInd.size(); n++)
            {
                int index_n = activePtsLocalInd[n];
                if (BEres.vtInd_BC[index_n] != 1)
                {
                    int x2_Ind = BEres.PP_Index[index_n]; // the second vertex index	
                    for (int xd = 0; xd < 3; xd++)
                    {
                        for (int yd = 0; yd < 3; yd++)
                        {
                            double value = hess_(m * 3 + xd, n * 3 + yd);
                            BEres.Hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
                        }
                    }
                }
            }
        }
    }

}

// store gradient and hessian value
void BarrierEnergy::store_grad_hess_PT(BarrierEnergyRes& BEres, Vector9d& grad_, Matrix9d& hess_, Eigen::Vector3i& activePtsLocalInd)
{
    for (int m = 0; m < activePtsLocalInd.size(); m++)
    {
        int index_m = activePtsLocalInd[m];
        if (BEres.vtInd_BC[index_m] != 1)
        {
            int x1_Ind = BEres.PP_Index[index_m]; // the first vertex index
            for (int xd = 0; xd < 3; xd++)
            {
                double value = grad_(m * 3 + xd, 1);
                BEres.Grad.emplace_back(x1_Ind * 3 + xd, value);
            }


            for (int n = 0; n < activePtsLocalInd.size(); n++)
            {
                int index_n = activePtsLocalInd[n];
                if (BEres.vtInd_BC[index_n] != 1)
                {
                    int x2_Ind = BEres.PP_Index[index_n]; // the second vertex index	
                    for (int xd = 0; xd < 3; xd++)
                    {
                        for (int yd = 0; yd < 3; yd++)
                        {
                            double value = hess_(m * 3 + xd, n * 3 + yd);
                            BEres.Hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
                        }
                    }
                }
            }
        }
    }

}

// store gradient and hessian value
void BarrierEnergy::store_grad_hess_PT(BarrierEnergyRes& BEres, Vector12d& grad_, Matrix12d& hess_, Eigen::Vector4i& activePtsLocalInd)
{
    for (int m = 0; m < activePtsLocalInd.size(); m++)
    {
        int index_m = activePtsLocalInd[m];
        if (BEres.vtInd_BC[index_m] != 1)
        {
            int x1_Ind = BEres.PP_Index[index_m]; // the first vertex index
            for (int xd = 0; xd < 3; xd++)
            {
                double value = grad_(m * 3 + xd, 1);
                BEres.Grad.emplace_back(x1_Ind * 3 + xd, value);
            }


            for (int n = 0; n < activePtsLocalInd.size(); n++)
            {
                int index_n = activePtsLocalInd[n];
                if (BEres.vtInd_BC[index_n] != 1)
                {
                    int x2_Ind = BEres.PP_Index[index_n]; // the second vertex index	
                    for (int xd = 0; xd < 3; xd++)
                    {
                        for (int yd = 0; yd < 3; yd++)
                        {
                            double value = hess_(m * 3 + xd, n * 3 + yd);
                            BEres.Hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
                        }
                    }
                }
            }
        }
    }

}

void BarrierEnergy::store_grad_hess_EE(BarrierEnergyRes& BEres, Vector12d& grad_final, Matrix12d& hess_final)
{
    for (int m = 0; m < 4; m++)
    {
        if (BEres.vtInd_BC[m] != 1)
        {
            int x1_Ind = BEres.PP_Index[m]; // the first vertex index
            for (int xd = 0; xd < 3; xd++)
            {
                double value = grad_final(m * 3 + xd, 1);
                BEres.Grad.emplace_back(x1_Ind * 3 + xd, value);
            }

            for (int n = 0; n < 4; n++)
            {
                if (BEres.vtInd_BC[n] != 1)
                {
                    int x2_Ind = BEres.PP_Index[n]; // the second vertex index	
                    for (int xd = 0; xd < 3; xd++)
                    {
                        for (int yd = 0; yd < 3; yd++)
                        {
                            double value = hess_final(m * 3 + xd, n * 3 + yd);
                            BEres.Hess.emplace_back(x1_Ind * 3 + xd, x2_Ind * 3 + yd, value);
                        }
                    }
                }
            }
        }
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



void Ground::valGradAndHess(BarrierEnergyRes& BEres, int index_i,  double coor_z, double d_hat, double distributedArea, double k_stiff, double dt)
{
    double dis2 = coor_z * coor_z;
    double d_hat2 = squaredDouble(d_hat);

    double bEnergy = -squaredDouble(dis2 - d_hat * d_hat) * log(dis2 / (d_hat * d_hat));
    double val_b = dt * dt * k_stiff * distributedArea * bEnergy;

    Eigen::Vector3d grad_b = Eigen::Vector3d::Zero();
    Eigen::Matrix3d hess_b = Eigen::Matrix3d::Zero();


    double partial_b_partial_dis = (-2.0 * (dis2 - d_hat2) * log(dis2 / d_hat2) - 1.0 / dis2 * squaredDouble(dis2 - d_hat2)) * 2.0 * std::sqrt(dis2); // 3
    double partial_2b_partial_dis2 = -4.0 * (3.0 * dis2 - d_hat2) * log(dis2 / d_hat2) - (4.0 / std::sqrt(dis2) + 8.0) * (dis2 - d_hat2) + 2.0 / dis2 * squaredDouble(dis2 - d_hat2); // 1                                                        

      
    Eigen::Vector3d partial_dis_partial_x = {0, 0, 1.0 };
    Eigen::Matrix3d partial_2dis_partial_x2 = Eigen::Matrix3d::Zero();

    Eigen::Vector3d grad_ = dt * dt * k_stiff * distributedArea * partial_dis_partial_x * partial_b_partial_dis;
    Eigen::Matrix3d hess_ = dt * dt * k_stiff * distributedArea * (partial_2b_partial_dis2 * partial_dis_partial_x * partial_dis_partial_x.transpose() + partial_b_partial_dis * partial_2dis_partial_x2);

    BEres.Val = val_b;
    for (int xd = 0; xd < 3; xd++)
    {
        double value = grad_[xd];
        BEres.Grad.emplace_back(index_i * 3 + xd, value);
    }
    for (int xd = 0; xd < 3; xd++)
    {
        for (int yd = 0; yd < 3; yd++)
        {
            double value = hess_(xd, yd);
            BEres.Hess.emplace_back(index_i * 3 + xd, index_i * 3 + yd, value);
        }
    }


}