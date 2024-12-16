#include "BarrierEnergy.h"

// compute the energy gradient and hessian 
double BarrierEnergy::val_PT(double& contactArea, double& dis2, FEMParamters& parameters)
{
    //contactArea = 1.0;
    double d_hat2 = squaredDouble(parameters.IPC_dis);
    return parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * compute_b(dis2, d_hat2);
}

// compute the energy gradient and hessian 
double BarrierEnergy::val_EE(double& contactArea, double& dis2, Mesh& tetSimMesh,
    Eigen::Vector4i& ptIndices, FEMParamters& parameters)
{
    //contactArea = 1.0;
    double d_hat2 = squaredDouble(parameters.IPC_dis);

    // the partial derivative of barrier energy b wrt distance d
    double val_b = compute_b(dis2, d_hat2);

    int P1 = ptIndices[0], P2 = ptIndices[1], Q1 = ptIndices[2], Q2 = ptIndices[3];
    int emin = std::min(P1, P2), emax = std::max(P1, P2);
    Eigen::Vector3d P1Coor = tetSimMesh.pos_node[P1], P2Coor = tetSimMesh.pos_node[P2], 
        Q1Coor = tetSimMesh.pos_node[Q1], Q2Coor = tetSimMesh.pos_node[Q2];
    Eigen::Vector3d P1Coor_Rest = tetSimMesh.pos_node_Rest[P1], P2Coor_Rest = tetSimMesh.pos_node_Rest[P2], 
        Q1Coor_Rest = tetSimMesh.pos_node_Rest[Q1], Q2Coor_Rest = tetSimMesh.pos_node_Rest[Q2];

    // the following codes calculate the mollifier-related value
    double eps_x = DIS::cal_EEM_eps_x(P1Coor_Rest, P2Coor_Rest, Q1Coor_Rest, Q2Coor_Rest);
    double val_ek = 0;
    DIS::compute_e(P1Coor, P2Coor, Q1Coor, Q2Coor, eps_x, val_ek);
    return parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * val_ek * val_b;
}


// compute the energy gradient and hessian 
void BarrierEnergy::gradAndHess_PT(std::vector<Eigen::Triplet<double>>& hessian_triplet, 
    std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_hess, 
    int& startIndex_grad, Eigen::Vector4i& ptIndices, int& type, double& dis2, Mesh& tetSimMesh,
    FEMParamters& parameters, bool ABD)
{
    double d_hat2 = squaredDouble(parameters.IPC_dis);
    // the partial derivative of barrier energy b wrt distance d
    double g_bd = compute_g_b(dis2, d_hat2); // 3
    double h_bd = compute_H_b(dis2, d_hat2); // 1  
    double contactArea = tetSimMesh.boundaryVertices_area[ptIndices[0]];
    //contactArea = 1.0;


    switch (type)
    {
    case 0:
    {
        std::vector<int> activePtsLocalInd = { ptIndices[0] , ptIndices[1] };
        BarrierEnergy::cal_and_assemble_gradAndHess(hessian_triplet, grad_triplet, startIndex_hess,
            activePtsLocalInd,startIndex_grad, tetSimMesh,
            parameters, contactArea, g_bd, h_bd, ABD);
    }
    break;

    case 1:
    {
        std::vector<int> activePtsLocalInd = { ptIndices[0] , ptIndices[2] };
        BarrierEnergy::cal_and_assemble_gradAndHess(hessian_triplet, grad_triplet, startIndex_hess,
            activePtsLocalInd, startIndex_grad, tetSimMesh,
            parameters, contactArea, g_bd, h_bd, ABD);
    }
    break;

    case 2:
    {
        std::vector<int> activePtsLocalInd = { ptIndices[0] , ptIndices[3] };
        BarrierEnergy::cal_and_assemble_gradAndHess(hessian_triplet, grad_triplet, startIndex_hess,
            activePtsLocalInd, startIndex_grad, tetSimMesh,
            parameters, contactArea, g_bd, h_bd, ABD);
    }
    break;

    case 3:
    {
        std::vector<int> activePtsLocalInd = { ptIndices[0] , ptIndices[1], ptIndices[2] };
        BarrierEnergy::cal_and_assemble_gradAndHess(hessian_triplet, grad_triplet, startIndex_hess,
            activePtsLocalInd, startIndex_grad, tetSimMesh,
            parameters, contactArea, g_bd, h_bd, ABD);
    }
    break;

    case 4:
    {
        std::vector<int> activePtsLocalInd = { ptIndices[0] , ptIndices[2] , ptIndices[3] };
        BarrierEnergy::cal_and_assemble_gradAndHess(hessian_triplet, grad_triplet, startIndex_hess,
            activePtsLocalInd, startIndex_grad, tetSimMesh,
            parameters, contactArea, g_bd, h_bd, ABD);
    }
    break;

    case 5:
    {
        std::vector<int> activePtsLocalInd = { ptIndices[0] , ptIndices[3] , ptIndices[1] };
        BarrierEnergy::cal_and_assemble_gradAndHess(hessian_triplet, grad_triplet, startIndex_hess,
            activePtsLocalInd, startIndex_grad, tetSimMesh,
            parameters, contactArea, g_bd, h_bd, ABD);
    }
    break;

    case 6:
    {
        std::vector<int> activePtsLocalInd = { ptIndices[0] , ptIndices[1] , ptIndices[2] , ptIndices[3] };
        BarrierEnergy::cal_and_assemble_gradAndHess(hessian_triplet, grad_triplet, startIndex_hess,
            activePtsLocalInd, startIndex_grad, tetSimMesh,
            parameters, contactArea, g_bd, h_bd, ABD);
    }
    break;

    }


}

void BarrierEnergy::cal_and_assemble_gradAndHess(std::vector<Eigen::Triplet<double>>& hessian_triplet,
    std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_hess,
    std::vector<int>& activePtsLocalInd,
    int& startIndex_grad, Mesh& tetSimMesh,
    FEMParamters& parameters, double& contactArea, double& g_bd, double& h_bd, bool ABD)
{
    if (activePtsLocalInd.size() == 2)
    {
        Vector6d g_dx = Vector6d::Zero();
        Matrix6d h_dx = Matrix6d::Zero();
        DIS::g_PP(tetSimMesh.pos_node[activePtsLocalInd[0]], tetSimMesh.pos_node[activePtsLocalInd[1]], g_dx); // 2
        DIS::H_PP(h_dx); // 4

        Vector6d grad = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * g_dx * g_bd;
        Matrix6d hessian = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
        makePD<double, 6>(hessian);
        assemble_gradAndHess<2>(hessian_triplet,grad_triplet, startIndex_hess, activePtsLocalInd,grad, hessian, startIndex_grad,tetSimMesh, ABD);
   
    }
    else if (activePtsLocalInd.size() == 3)
    {
        Vector9d g_dx = Vector9d::Zero();
        Matrix9d h_dx = Matrix9d::Zero();
        DIS::g_PE(tetSimMesh.pos_node[activePtsLocalInd[0]], tetSimMesh.pos_node[activePtsLocalInd[1]], tetSimMesh.pos_node[activePtsLocalInd[2]], g_dx); // 2
        DIS::H_PE(tetSimMesh.pos_node[activePtsLocalInd[0]], tetSimMesh.pos_node[activePtsLocalInd[1]], tetSimMesh.pos_node[activePtsLocalInd[2]], h_dx); // 4

        Vector9d grad = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * g_dx * g_bd;
        Matrix9d hessian = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea *
            (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
        makePD<double, 9>(hessian);
        assemble_gradAndHess<3>(hessian_triplet, grad_triplet, startIndex_hess, activePtsLocalInd, grad, hessian, startIndex_grad, tetSimMesh, ABD);
    }
    else if (activePtsLocalInd.size() == 4)
    {
        Vector12d g_dx = Vector12d::Zero();
        Matrix12d h_dx = Matrix12d::Zero();
        DIS::g_PT(tetSimMesh.pos_node[activePtsLocalInd[0]], tetSimMesh.pos_node[activePtsLocalInd[1]], tetSimMesh.pos_node[activePtsLocalInd[2]], tetSimMesh.pos_node[activePtsLocalInd[3]], g_dx); // 2
        DIS::H_PT(tetSimMesh.pos_node[activePtsLocalInd[0]], tetSimMesh.pos_node[activePtsLocalInd[1]], tetSimMesh.pos_node[activePtsLocalInd[2]], tetSimMesh.pos_node[activePtsLocalInd[3]], h_dx); // 4

        Vector12d grad = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * g_dx * g_bd;
        Matrix12d hessian = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
        makePD<double, 12>(hessian);
        assemble_gradAndHess<4>(hessian_triplet, grad_triplet, startIndex_hess, activePtsLocalInd, grad, hessian, startIndex_grad, tetSimMesh, ABD);
    }
    else
    {

    }



   
}




template <int size>
void assemble_gradAndHess(std::vector<Eigen::Triplet<double>>& hessian_triplet,
    std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_hess,
    std::vector<int>& activePts, Eigen::Matrix<double, size * 3, 1>& grad_, Eigen::Matrix<double, size * 3, size * 3>& hess_,
    int& startIndex_grad, Mesh& tetSimMesh, bool ABD)
{

    if (!ABD)
    {
        for (int j = 0; j < activePts.size(); j++)
        {
            int pt1 = activePts[j];
            for (int xd = 0; xd < 3; xd++)
            {
                double value = grad_[j * 3 + xd];
                grad_triplet[startIndex_grad + j * 3 + xd] = { pt1 * 3 + xd, value };
            }

            for (int q = 0; q < activePts.size(); q++)
            {
                int pt2 = activePts[q];

                for (int xd = 0; xd < 3; xd++)
                {
                    for (int yd = 0; yd < 3; yd++)
                    {
                        hessian_triplet[startIndex_hess + j * 36 + q * 9 + xd * 3 + yd] = { pt1 * 3 + xd, pt2 * 3 + yd, hess_(j * 3 + xd, q * 3 + yd) };
                    }
                }
            }


        }
    }
    else
    {

        for (int j = 0; j < activePts.size(); j++)
        {
            int pt1 = activePts[j];
            Eigen::Matrix<double, 3, 12> Jx1 = build_Jx_matrix_for_ABD(tetSimMesh.pos_node_Rest[pt1]);
            Vector12d energy_wrt_q_grad = Jx1.transpose() * grad_.block<3, 1>(j * 3, 0);
            int AB_index_1 = tetSimMesh.index_node[pt1][0]; // the ABD object's index

            // assemble gradient
            for (int xd = 0; xd < 12; xd++)
            {
                double value = energy_wrt_q_grad[xd];
                grad_triplet[startIndex_grad + xd] = { AB_index_1 * 12 + xd, value };
            }

            // assemble hessian
            for (int q = 0; q < activePts.size(); q++)
            {
                int pt2 = activePts[q];
                Eigen::Matrix<double, 3, 12> Jx2 = build_Jx_matrix_for_ABD(tetSimMesh.pos_node_Rest[pt2]);
                int AB_index_2 = tetSimMesh.index_node[pt2][0]; // the ABD object's index
                Matrix12d energy_wrt_q_hess = Jx1.transpose() * hess_.block<3, 3>(j * 3, q * 3) * Jx2;


                for (int xd = 0; xd < 12; xd++)
                {
                    for (int yd = 0; yd < 12; yd++)
                    {
                        hessian_triplet[startIndex_hess + j * 144 * 2 + q * 144 + xd * 12 + yd] = { AB_index_1 * 12 + xd, AB_index_2 * 12 + yd, energy_wrt_q_hess(xd, yd) };
                    }
                }
            }
        }

    }


    

}




// compute the energy gradient and hessian 
void BarrierEnergy::gradAndHess_EE(std::vector<Eigen::Triplet<double>>& hessian_triplet, 
    std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_hess, 
    int& startIndex_grad, 
    Eigen::Vector4i& ptIndices, int& type, double& dis2, Mesh& tetSimMesh,
    FEMParamters& parameters, bool ABD)
{
    double d_hat2 = squaredDouble(parameters.IPC_dis);
    // the partial derivative of barrier energy b wrt distance d
    double g_bd = compute_g_b(dis2, d_hat2); // 3
    double h_bd = compute_H_b(dis2, d_hat2); // 1     

    int P1 = ptIndices[0], P2 = ptIndices[1], Q1 = ptIndices[2], Q2 = ptIndices[3];
    int emin = std::min(P1, P2), emax = std::max(P1, P2);
    Eigen::Vector3d P1Coor = tetSimMesh.pos_node[P1], P2Coor = tetSimMesh.pos_node[P2], Q1Coor = tetSimMesh.pos_node[Q1], Q2Coor = tetSimMesh.pos_node[Q2];
    Eigen::Vector3d P1Coor_Rest = tetSimMesh.pos_node_Rest[P1], P2Coor_Rest = tetSimMesh.pos_node_Rest[P2], Q1Coor_Rest = tetSimMesh.pos_node_Rest[Q1], Q2Coor_Rest = tetSimMesh.pos_node_Rest[Q2];

    double contactArea = tetSimMesh.boundaryEdges_area[emin][emax];
    //contactArea = 1.0;
    double val_b = compute_b(dis2, d_hat2);
    Vector12d grad_b = Vector12d::Zero();
    Matrix12d hess_b = Matrix12d::Zero();

    // the following codes calculate the mollifier-related value
    double eps_x = DIS::cal_EEM_eps_x(P1Coor_Rest, P2Coor_Rest, Q1Coor_Rest, Q2Coor_Rest);
    double val_ek = 0;
    Vector12d grad_ek = Vector12d::Zero();
    Matrix12d hess_ek = Matrix12d::Zero();
    DIS::compute_e(P1Coor, P2Coor, Q1Coor, Q2Coor, eps_x, val_ek);
    DIS::compute_e_g(P1Coor, P2Coor, Q1Coor, Q2Coor, eps_x, grad_ek);
    DIS::compute_e_g(P1Coor, P2Coor, Q1Coor, Q2Coor, eps_x, grad_ek);
    DIS::compute_e_H(P1Coor, P2Coor, Q1Coor, Q2Coor, eps_x, hess_ek);

    Vector12d grad;
    Matrix12d hessian;

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
        std::vector<int> activePtsLocalInd = { 0, 2 };
        project_grad_to_full<2>(activePtsLocalInd, grad_, hess_, grad_b, hess_b); // since we add edge-edge mollifier, we have to project the hessian to full 12x12 matrix

        // now calcualte the final val, grad and hess considering the mollifier   
        grad = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea* (grad_ek * val_b + val_ek * grad_b);
        hessian = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea* (hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b);
        makePD<double, 12>(hessian);

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
        std::vector<int> activePtsLocalInd = { 0, 3 };
        project_grad_to_full<2>(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        grad = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (grad_ek * val_b + val_ek * grad_b);
        hessian = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b);
        makePD<double, 12>(hessian);
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
        std::vector<int> activePtsLocalInd = { 0, 2, 3 };
        project_grad_to_full<3>(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        grad = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (grad_ek * val_b + val_ek * grad_b);
        hessian = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b);
        makePD<double, 12>(hessian);

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
        std::vector<int> activePtsLocalInd = { 1, 2 };
        project_grad_to_full<2>(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        grad = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (grad_ek * val_b + val_ek * grad_b);
        hessian = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b);
        makePD<double, 12>(hessian);
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
        std::vector<int> activePtsLocalInd = { 1, 3 };
        project_grad_to_full<2>(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        grad = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (grad_ek * val_b + val_ek * grad_b);
        hessian = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b);
        makePD<double, 12>(hessian);

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
        std::vector<int> activePtsLocalInd = { 1, 2, 3 };
        project_grad_to_full<3>(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        grad = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (grad_ek * val_b + val_ek * grad_b);
        hessian = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b);
        makePD<double, 12>(hessian);

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
        std::vector<int> activePtsLocalInd = { 2, 0, 1 };
        project_grad_to_full<3>(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        grad = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (grad_ek * val_b + val_ek * grad_b);
        hessian = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b);
        makePD<double, 12>(hessian);
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
        std::vector<int> activePtsLocalInd = { 3, 0, 1 };
        project_grad_to_full<3>(activePtsLocalInd, grad_, hess_, grad_b, hess_b);

        // now calcualte the final val, grad and hess considering the mollifier
        grad = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (grad_ek * val_b + val_ek * grad_b);
        hessian = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b);
        makePD<double, 12>(hessian);
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
        grad = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (grad_ek * val_b + val_ek * grad_b);
        hessian = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (hess_ek * val_b + grad_ek * grad_b.transpose() + grad_b * grad_ek.transpose() + val_ek * hess_b);
        makePD<double, 12>(hessian);
    }
    break;


    }

    std::vector<int> activePtsLocalInd = { ptIndices[0] , ptIndices[1] , ptIndices[2] , ptIndices[3] };
    assemble_gradAndHess<4>(hessian_triplet, grad_triplet, startIndex_hess, activePtsLocalInd, grad, hessian, startIndex_grad, tetSimMesh, ABD);
   

}



template<int size>
void BarrierEnergy::project_grad_to_full(std::vector<int>& activePtsLocalInd, Eigen::Matrix<double, size * 3, 1>& grad_,
    Eigen::Matrix<double, size * 3, size * 3>& hess_, Vector12d& grad_full, Matrix12d& hess__full)
{
    for (int i = 0; i < activePtsLocalInd.size(); i++)
    {
        int index_i = activePtsLocalInd[i];

        // project gradient
        for (int s = 0; s < 3; s++)
        {
            grad_full[index_i * 3 + s] = grad_[i * 3 + s];
        }

        // project hessian
        for (int j = 0; j < activePtsLocalInd.size(); j++)
        {
            int index_j = activePtsLocalInd[j];

            hess__full.block(index_i * 3, index_j * 3, 3, 3) = hess_.block(i * 3, j * 3, 3, 3);
        }

    }
}




double BarrierEnergy::compute_b(double& d2, double& dHat2)
{
    return -(d2 - dHat2) * (d2 - dHat2) * log(d2 / dHat2);
}

double BarrierEnergy::compute_g_b(double& d2, double& dHat2)
{
    double t = d2 - dHat2;
    return t * std::log(d2 / dHat2) * -2.0 - (t * t) / d2;
}

double BarrierEnergy::compute_H_b(double& d2, double& dHat2)
{
    double t = d2 - dHat2;
    return (std::log(d2 / dHat2) * -2.0 - t * 4.0 / d2) + 1.0 / (d2 * d2) * (t * t);
}



double Ground::val(double& coor_z2, double& contactArea, FEMParamters& parameters)
{
    double d_hat2 = parameters.IPC_dis * parameters.IPC_dis;
    return parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * BarrierEnergy::compute_b(coor_z2, d_hat2);
}

void Ground::gradAndHess(std::vector<Eigen::Triplet<double>>& hessian_triplet, 
    std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_hess, 
    int& startIndex_grad, Mesh& tetSimMesh,
    int& index_i, double& coor_z2, double& contactArea, FEMParamters& parameters, bool ABD)
{
    double d_hat2 = parameters.IPC_dis * parameters.IPC_dis;
    double g_bd = BarrierEnergy::compute_g_b(coor_z2, d_hat2); // 3
    double h_bd = BarrierEnergy::compute_H_b(coor_z2, d_hat2); // 1                                                        

    Eigen::Vector3d g_dx = { 0, 0, 2.0 * std::sqrt(coor_z2) };
    Eigen::Matrix3d h_dx = Eigen::Matrix3d::Zero();
    h_dx(2, 2) = 2.0;

    Vector3d grad = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * g_dx * g_bd;

    Matrix3d hessian = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);
    makePD<double, 3>(hessian);
    std::vector<int> activePtsLocalInd = { index_i };
    assemble_gradAndHess<1>(hessian_triplet, grad_triplet, startIndex_hess, activePtsLocalInd, grad, hessian, startIndex_grad, tetSimMesh, ABD);



}