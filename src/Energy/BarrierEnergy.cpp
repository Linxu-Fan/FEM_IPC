#include "BarrierEnergy.h"

// compute the energy gradient and hessian 
double BarrierEnergy::val_PT(
    double& contactArea, 
    double& dis2, 
    FEMParamters& parameters)
{
    //contactArea = 1.0;
    double d_hat2 = squaredDouble(parameters.IPC_dis);
    return parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * compute_b(dis2, d_hat2);
}

// compute the energy gradient and hessian 
double BarrierEnergy::val_EE(
    double& contactArea, 
    double& dis2, 
    const double val_ek,
    FEMParamters& parameters)
{
    double d_hat2 = squaredDouble(parameters.IPC_dis);
    // the partial derivative of barrier energy b wrt distance d
    double val_b = compute_b(dis2, d_hat2);

    return parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * val_ek * val_b;
}



void BarrierEnergy::cal_gradAndHess_PP_ABD(
    Vector6d& grad,
    Matrix6d& hessian,
    const Eigen::Vector3d& pos1,
    const Eigen::Vector3d& pos2,
    std::vector<Eigen::Vector3d>& contactForce_node,
    FEMParamters& parameters,
    double& contactArea,
    double& dis2)
{
    double d_hat2 = squaredDouble(parameters.IPC_dis);
    // the partial derivative of barrier energy b wrt distance d
    double g_bd = compute_g_b(dis2, d_hat2); // 3
    double h_bd = compute_H_b(dis2, d_hat2); // 1  

    Vector6d g_dx = Vector6d::Zero();
    Matrix6d h_dx = Matrix6d::Zero();
    DIS::g_PP(pos1, pos2, g_dx); // 2
    DIS::H_PP(h_dx); // 4

    grad = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * g_dx * g_bd;
    hessian = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);

    contactForce_node.push_back(grad.block(0, 0, 3, 1));
    contactForce_node.push_back(grad.block(3, 0, 3, 1));

    makePD<double, 6>(hessian);

}



void BarrierEnergy::cal_gradAndHess_PE_ABD(
    Vector9d& grad,
    Matrix9d& hessian,
    const Eigen::Vector3d& pos1,
    const Eigen::Vector3d& pos2,
    const Eigen::Vector3d& pos3,
    std::vector<Eigen::Vector3d>& contactForce_node,
    FEMParamters& parameters,
    double& contactArea,
    double& dis2)
{
    double d_hat2 = squaredDouble(parameters.IPC_dis);
    // the partial derivative of barrier energy b wrt distance d
    double g_bd = compute_g_b(dis2, d_hat2); // 3
    double h_bd = compute_H_b(dis2, d_hat2); // 1  


    Vector9d g_dx = Vector9d::Zero();
    Matrix9d h_dx = Matrix9d::Zero();
    DIS::g_PE(pos1, pos2, pos3, g_dx); // 2
    DIS::H_PE(pos1, pos2, pos3, h_dx); // 4

    grad = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * g_dx * g_bd;
    hessian = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea *
        (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);

    contactForce_node.push_back(grad.block(0, 0, 3, 1));
    contactForce_node.push_back(grad.block(3, 0, 3, 1));
    contactForce_node.push_back(grad.block(6, 0, 3, 1));

    makePD<double, 9>(hessian);

}


void BarrierEnergy::cal_gradAndHess_PT_ABD(
    Vector12d& grad,
    Matrix12d& hessian,
    const Eigen::Vector3d& pos1,
    const Eigen::Vector3d& pos2,
    const Eigen::Vector3d& pos3,
    const Eigen::Vector3d& pos4,
    std::vector<Eigen::Vector3d>& contactForce_node,
    FEMParamters& parameters,
    double& contactArea,
    double& dis2)
{
    double d_hat2 = squaredDouble(parameters.IPC_dis);
    // the partial derivative of barrier energy b wrt distance d
    double g_bd = compute_g_b(dis2, d_hat2); // 3
    double h_bd = compute_H_b(dis2, d_hat2); // 1  


    Vector12d g_dx = Vector12d::Zero();
    Matrix12d h_dx = Matrix12d::Zero();
    DIS::g_PT(pos1, pos2, pos3, pos4, g_dx); // 2
    DIS::H_PT(pos1, pos2, pos3, pos4, h_dx); // 4

    grad = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * g_dx * g_bd;
    hessian = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);

    contactForce_node.push_back(grad.block(0, 0, 3, 1));
    contactForce_node.push_back(grad.block(3, 0, 3, 1));
    contactForce_node.push_back(grad.block(6, 0, 3, 1));
    contactForce_node.push_back(grad.block(9, 0, 3, 1));

    makePD<double, 12>(hessian);

}



template <int size>
void assemble_gradAndHess_ABD(
    std::vector<Eigen::Triplet<double>>& hessian_triplet,
    std::vector<std::pair<int, double>>& grad_triplet,
    int& startIndex_grad,
    int& startIndex_hess,
    std::vector<int>& pt_obj_index,
    Eigen::Matrix<double, size * 3, 1>& grad_,
    Eigen::Matrix<double, size * 3, size * 3>& hess_,    
    const std::vector<Eigen::Vector3d>& pos_node_Rest)
{

    for (int j = 0; j < pt_obj_index.size(); j++)
    {
        Eigen::Matrix<double, 3, 12> Jx1 = build_Jx_matrix_for_ABD(pos_node_Rest[j]);
        Vector12d energy_wrt_q_grad = Jx1.transpose() * grad_.block(j * 3, 0, 3, 1);
        int AB_index_1 = pt_obj_index[j]; // the ABD object's index

        // assemble gradient
        for (int xd = 0; xd < 12; xd++)
        {
            double value = energy_wrt_q_grad[xd];
            grad_triplet[startIndex_grad + j * 12 + xd] = { AB_index_1 * 12 + xd, value };
        }

        // assemble hessian
        for (int q = 0; q < pt_obj_index.size(); q++)
        {
            Eigen::Matrix<double, 3, 12> Jx2 = build_Jx_matrix_for_ABD(pos_node_Rest[q]);
            Matrix12d energy_wrt_q_hess = Jx1.transpose() * hess_.block(j * 3, q * 3, 3, 3) * Jx2;
            int AB_index_2 = pt_obj_index[q]; // the ABD object's index


            for (int xd = 0; xd < 12; xd++)
            {
                for (int yd = 0; yd < 12; yd++)
                {
                    hessian_triplet[startIndex_hess + j * pt_obj_index.size() * 144 + q * 144 + xd * 12 + yd] = { AB_index_1 * 12 + xd, AB_index_2 * 12 + yd, energy_wrt_q_hess(xd, yd) };
                }
            }
        }
    }




}



void BarrierEnergy::gradAndHess_PT(
    std::vector<Eigen::Triplet<double>>& hessian_triplet,
    std::vector<std::pair<int, double>>& grad_triplet,
    int& startIndex_hess,
    int& startIndex_grad,
    double& dis2,
    const Vector5i& PT_Contact,
    triMesh& triSimMesh,
    std::vector<Eigen::Vector3d>& contactForce_node,
    FEMParamters& parameters)
{
    int obj_1 = PT_Contact[0], vert = PT_Contact[1], obj_2 = PT_Contact[2], tria = PT_Contact[3];
    double contactArea = triSimMesh.allObjects[obj_1].objectSurfaceMesh.boundaryVertices_area[vert];

    Eigen::Vector3d P = triSimMesh.allObjects[obj_1].pos_node_surface[vert];
    Eigen::Vector3i tri = triSimMesh.allObjects[obj_1].objectSurfaceMesh.boundaryTriangles[tria];


    std::vector<int> pt_obj_index;
    std::vector<Eigen::Vector3d> pos_node_Rest;


    switch (PT_Contact[4])
    {
    case 0:
    {
        Vector6d grad = Vector6d::Zero();
        Matrix6d hessian = Matrix6d::Zero();
        cal_gradAndHess_PP_ABD(grad, hessian, P, triSimMesh.allObjects[obj_2].pos_node_surface[tri[0]], contactForce_node, parameters, contactArea, dis2);


        pt_obj_index.push_back(obj_1);
        pt_obj_index.push_back(obj_2);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_1].objectSurfaceMesh.vertices[vert]);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_2].objectSurfaceMesh.vertices[tri[0]]);
        assemble_gradAndHess_ABD<2>(hessian_triplet, grad_triplet, startIndex_grad, startIndex_hess, pt_obj_index, grad, hessian, pos_node_Rest);

    }
    break;

    case 1:
    {
        Vector6d grad = Vector6d::Zero();
        Matrix6d hessian = Matrix6d::Zero();
        cal_gradAndHess_PP_ABD(grad, hessian, P, triSimMesh.allObjects[obj_2].pos_node_surface[tri[1]], contactForce_node, parameters, contactArea, dis2);


        pt_obj_index.push_back(obj_1);
        pt_obj_index.push_back(obj_2);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_1].objectSurfaceMesh.vertices[vert]);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_2].objectSurfaceMesh.vertices[tri[1]]);
        assemble_gradAndHess_ABD<2>(hessian_triplet, grad_triplet, startIndex_grad, startIndex_hess, pt_obj_index, grad, hessian, pos_node_Rest);


    }
    break;

    case 2:
    {
        Vector6d grad = Vector6d::Zero();
        Matrix6d hessian = Matrix6d::Zero();
        cal_gradAndHess_PP_ABD(grad, hessian, P, triSimMesh.allObjects[obj_2].pos_node_surface[tri[2]], contactForce_node, parameters, contactArea, dis2);


        pt_obj_index.push_back(obj_1);
        pt_obj_index.push_back(obj_2);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_1].objectSurfaceMesh.vertices[vert]);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_2].objectSurfaceMesh.vertices[tri[2]]);
        assemble_gradAndHess_ABD<2>(hessian_triplet, grad_triplet, startIndex_grad, startIndex_hess, pt_obj_index, grad, hessian, pos_node_Rest);



    }
    break;

    case 3:
    {
        Vector9d grad = Vector9d::Zero();
        Matrix9d hessian = Matrix9d::Zero();
        cal_gradAndHess_PE_ABD(grad, hessian, P, triSimMesh.allObjects[obj_2].pos_node_surface[tri[0]],
            triSimMesh.allObjects[obj_2].pos_node_surface[tri[1]],
            contactForce_node, parameters, contactArea, dis2);


        pt_obj_index.push_back(obj_1);
        pt_obj_index.push_back(obj_2);
        pt_obj_index.push_back(obj_2);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_1].objectSurfaceMesh.vertices[vert]);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_2].objectSurfaceMesh.vertices[tri[0]]);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_2].objectSurfaceMesh.vertices[tri[1]]);
        assemble_gradAndHess_ABD<3>(hessian_triplet, grad_triplet, startIndex_grad, startIndex_hess, pt_obj_index, grad, hessian, pos_node_Rest);
    }
    break;

    case 4:
    {
        Vector9d grad = Vector9d::Zero();
        Matrix9d hessian = Matrix9d::Zero();
        cal_gradAndHess_PE_ABD(grad, hessian, P, triSimMesh.allObjects[obj_2].pos_node_surface[tri[1]],
            triSimMesh.allObjects[obj_2].pos_node_surface[tri[2]],
            contactForce_node, parameters, contactArea, dis2);


        pt_obj_index.push_back(obj_1);
        pt_obj_index.push_back(obj_2);
        pt_obj_index.push_back(obj_2);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_1].objectSurfaceMesh.vertices[vert]);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_2].objectSurfaceMesh.vertices[tri[1]]);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_2].objectSurfaceMesh.vertices[tri[2]]);
        assemble_gradAndHess_ABD<3>(hessian_triplet, grad_triplet, startIndex_grad, startIndex_hess, pt_obj_index, grad, hessian, pos_node_Rest);
    }
    break;

    case 5:
    {
        Vector9d grad = Vector9d::Zero();
        Matrix9d hessian = Matrix9d::Zero();
        cal_gradAndHess_PE_ABD(grad, hessian, P, triSimMesh.allObjects[obj_2].pos_node_surface[tri[2]],
            triSimMesh.allObjects[obj_2].pos_node_surface[tri[0]],
            contactForce_node, parameters, contactArea, dis2);


        pt_obj_index.push_back(obj_1);
        pt_obj_index.push_back(obj_2);
        pt_obj_index.push_back(obj_2);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_1].objectSurfaceMesh.vertices[vert]);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_2].objectSurfaceMesh.vertices[tri[2]]);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_2].objectSurfaceMesh.vertices[tri[0]]);
        assemble_gradAndHess_ABD<3>(hessian_triplet, grad_triplet, startIndex_grad, startIndex_hess, pt_obj_index, grad, hessian, pos_node_Rest);
    }
    break;

    case 6:
    {
        Vector12d grad = Vector12d::Zero();
        Matrix12d hessian = Matrix12d::Zero();
        cal_gradAndHess_PT_ABD(grad, hessian, P, triSimMesh.allObjects[obj_2].pos_node_surface[tri[0]],
            triSimMesh.allObjects[obj_2].pos_node_surface[tri[1]], triSimMesh.allObjects[obj_2].pos_node_surface[tri[2]],
            contactForce_node, parameters, contactArea, dis2);


        pt_obj_index.push_back(obj_1);
        pt_obj_index.push_back(obj_2);
        pt_obj_index.push_back(obj_2);
        pt_obj_index.push_back(obj_2);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_1].objectSurfaceMesh.vertices[vert]);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_2].objectSurfaceMesh.vertices[tri[0]]);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_2].objectSurfaceMesh.vertices[tri[1]]);
        pos_node_Rest.push_back(triSimMesh.allObjects[obj_2].objectSurfaceMesh.vertices[tri[2]]);
        assemble_gradAndHess_ABD<4>(hessian_triplet, grad_triplet, startIndex_grad, startIndex_hess, pt_obj_index, grad, hessian, pos_node_Rest);
    }
    break;

    }


}




// compute the energy gradient and hessian 
void BarrierEnergy::gradAndHess_EE(
    std::vector<Eigen::Triplet<double>>& hessian_triplet, 
    std::vector<std::pair<int, double>>& grad_triplet, 
    int& startIndex_hess, 
    int& startIndex_grad, 
    const Vector5i& EE_Contact,
    double& dis2, 
    triMesh& triSimMesh,
    std::vector<Eigen::Vector3d>& contactForce_node,
    FEMParamters& parameters)
{
    double d_hat2 = squaredDouble(parameters.IPC_dis);
    // the partial derivative of barrier energy b wrt distance d
    double g_bd = compute_g_b(dis2, d_hat2); // 3
    double h_bd = compute_H_b(dis2, d_hat2); // 1     


    int obj_1 = EE_Contact[0], obj_2 = EE_Contact[2];
    int edge1 = EE_Contact[1], edge2 = EE_Contact[3];
    Eigen::Vector2i edge1_vert = triSimMesh.allObjects[obj_1].objectSurfaceMesh.index_boundaryEdge[edge1];
    Eigen::Vector2i edge2_vert = triSimMesh.allObjects[obj_2].objectSurfaceMesh.index_boundaryEdge[edge2];
    int P1 = edge1_vert[0], P2 = edge1_vert[1], Q1 = edge2_vert[0], Q2 = edge2_vert[1];
    double contactArea = triSimMesh.allObjects[obj_1].objectSurfaceMesh.boundaryEdges_area[P1][P2];


    Eigen::Vector3d P1Coor = triSimMesh.allObjects[obj_1].pos_node_surface[P1], P2Coor = triSimMesh.allObjects[obj_1].pos_node_surface[P2];
    Eigen::Vector3d Q1Coor = triSimMesh.allObjects[obj_2].pos_node_surface[Q1], Q2Coor = triSimMesh.allObjects[obj_2].pos_node_surface[Q2];
    Eigen::Vector3d P1Coor_Rest = triSimMesh.allObjects[obj_1].objectSurfaceMesh.vertices[P1], P2Coor_Rest = triSimMesh.allObjects[obj_1].objectSurfaceMesh.vertices[P2];
    Eigen::Vector3d Q1Coor_Rest = triSimMesh.allObjects[obj_2].objectSurfaceMesh.vertices[Q1], Q2Coor_Rest = triSimMesh.allObjects[obj_2].objectSurfaceMesh.vertices[Q2];

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

    switch (EE_Contact[4])
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

    std::vector<int> pt_obj_index = { obj_1 , obj_1 , obj_2 , obj_2 };

    contactForce_node.push_back(grad.block(0, 0, 3, 1));
    contactForce_node.push_back(grad.block(3, 0, 3, 1));
    contactForce_node.push_back(grad.block(6, 0, 3, 1));
    contactForce_node.push_back(grad.block(9, 0, 3, 1));

    std::vector<Eigen::Vector3d> pos_node_Rest;
    pos_node_Rest.push_back(P1Coor_Rest);
    pos_node_Rest.push_back(P2Coor_Rest);
    pos_node_Rest.push_back(Q1Coor_Rest);
    pos_node_Rest.push_back(Q2Coor_Rest);
    assemble_gradAndHess_ABD<4>(hessian_triplet, grad_triplet, startIndex_grad, startIndex_hess, pt_obj_index, grad, hessian, pos_node_Rest);
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

void Ground::gradAndHess(
    std::vector<Eigen::Triplet<double>>& hessian_triplet, 
    std::vector<std::pair<int, double>>& grad_triplet, 
    int& startIndex_hess, 
    int& startIndex_grad, 
    const Vector2i PG_Contact,
    triMesh& triSimMesh,
    std::vector<Eigen::Vector3d>& contactForce_node,
    FEMParamters& parameters)
{
    int obj_1 = PG_Contact[0];
    int pt_index = PG_Contact[1];
    Eigen::Vector3d P = triSimMesh.allObjects[obj_1].pos_node_surface[pt_index];
    double coor_z2 = P[2] * P[2];

    double contactArea = triSimMesh.allObjects[obj_1].objectSurfaceMesh.boundaryTriangles_area[pt_index];

    double d_hat2 = parameters.IPC_dis * parameters.IPC_dis;
    double g_bd = BarrierEnergy::compute_g_b(coor_z2, d_hat2); // 3
    double h_bd = BarrierEnergy::compute_H_b(coor_z2, d_hat2); // 1                                                        

    Eigen::Vector3d g_dx = { 0, 0, 2.0 * std::sqrt(coor_z2) };
    Eigen::Matrix3d h_dx = Eigen::Matrix3d::Zero();
    h_dx(2, 2) = 2.0;

    Vector3d grad = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * g_dx * g_bd;

    Matrix3d hessian = parameters.dt * parameters.dt * parameters.IPC_kStiffness * contactArea * (h_bd * g_dx * g_dx.transpose() + g_bd * h_dx);

    contactForce_node.push_back(grad);

    makePD<double, 3>(hessian);

    std::vector<int> pt_obj_index = { obj_1 };
    std::vector<Eigen::Vector3d> pos_node_Rest;
    Eigen::Vector3d P1_Rest = triSimMesh.allObjects[obj_1].objectSurfaceMesh.vertices[pt_index];
    pos_node_Rest.push_back(P1_Rest);

    assemble_gradAndHess_ABD<1>(hessian_triplet, grad_triplet, startIndex_grad, startIndex_hess, pt_obj_index, grad, hessian, pos_node_Rest);

}