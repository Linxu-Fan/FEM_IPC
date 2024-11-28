#include "tools.h"  


// This cpp file contains a few useful tools

Eigen::Vector3d compute_contact_force(objMeshFormat& obj1, objMeshFormat& obj2, Eigen::Vector3d& bbx_min, Eigen::Vector3d& bbx_max, FEMParamters& parameters)
{
    std::vector<Eigen::Vector3d> forceObj2(obj2.vertices.size(), Eigen::Vector3d::Zero());
    
	// PT contact force: obj1(P) -> OBJ2(T)
	for (int v = 0; v < obj1.vertices.size(); v++)
	{
		Eigen::Vector3d P = obj1.vertices[v];
		if (insideBoundingBox(P, bbx_min, bbx_max))
		{
			for (int f = 0; f < obj2.faces.size(); f++)
			{
				Eigen::Vector3i tri = obj2.faces[f];
				Eigen::Vector3d A = obj2.vertices[tri[0]];
				Eigen::Vector3d B = obj2.vertices[tri[1]];
				Eigen::Vector3d C = obj2.vertices[tri[2]];
				if (insideBoundingBox(A, bbx_min, bbx_max) && insideBoundingBox(B, bbx_min, bbx_max) && insideBoundingBox(C, bbx_min, bbx_max))
				{
					if (pointTriangleCCDBroadphase(P, A, B, C, parameters.IPC_dis))
					{
						int type = DIS::dType_PT(P, A, B, C);
						double dis2 = 0;
						DIS::computePointTriD(P, A, B, C, dis2);

						if (dis2 <= squaredDouble(parameters.IPC_dis)) // only calculate the energy when the distance is smaller than the threshold
						{
                            double ipc_dis2 = parameters.IPC_dis * parameters.IPC_dis;
                            double g_bd = BarrierEnergy::compute_g_b(dis2, ipc_dis2);
                            switch (type)
                            {
                            case 0:
                            {
                                Vector6d g_dx = Vector6d::Zero();
                                DIS::g_PP(P, A, g_dx); // 2
                                Vector6d grad = g_dx * g_bd;

                                Eigen::Vector3d FA = { grad[3] , grad[4] , grad[5] };
                                forceObj2[tri[0]] += FA;
                            }
                            break;

                            case 1:
                            {
                                Vector6d g_dx = Vector6d::Zero();
                                DIS::g_PP(P, B, g_dx); // 2
                                Vector6d grad = g_dx * g_bd;

                                Eigen::Vector3d FB = { grad[3] , grad[4] , grad[5] };
                                forceObj2[tri[1]] += FB;
                            }
                            break;

                            case 2:
                            {
                                Vector6d g_dx = Vector6d::Zero();
                                DIS::g_PP(P, C, g_dx); // 2
                                Vector6d grad = g_dx * g_bd;

                                Eigen::Vector3d FC = { grad[3] , grad[4] , grad[5] };
                                forceObj2[tri[2]] += FC;
                            }
                            break;

                            case 3:
                            {
                                Vector9d g_dx = Vector9d::Zero();
                                DIS::g_PE(P, A, B, g_dx); // 2
                                Vector9d grad = g_dx * g_bd;

                                Eigen::Vector3d FA = { grad[3] , grad[4] , grad[5] };
                                Eigen::Vector3d FB = { grad[6] , grad[7] , grad[8] };
                                forceObj2[tri[0]] += FA;
                                forceObj2[tri[1]] += FB;
                            }
                            break;

                            case 4:
                            {
                                Vector9d g_dx = Vector9d::Zero();
                                DIS::g_PE(P, B, C, g_dx); // 2
                                Vector9d grad = g_dx * g_bd;

                                Eigen::Vector3d FB = { grad[3] , grad[4] , grad[5] };
                                Eigen::Vector3d FC = { grad[6] , grad[7] , grad[8] };
                                forceObj2[tri[1]] += FB;
                                forceObj2[tri[2]] += FC;
                            }
                            break;

                            case 5:
                            {
                                Vector9d g_dx = Vector9d::Zero();
                                DIS::g_PE(P, C, A, g_dx); // 2
                                Vector9d grad = g_dx * g_bd;

                                Eigen::Vector3d FC = { grad[3] , grad[4] , grad[5] };
                                Eigen::Vector3d FA = { grad[6] , grad[7] , grad[8] };
                                forceObj2[tri[2]] += FC;
                                forceObj2[tri[0]] += FA;

                            }
                            break;

                            case 6:
                            {
                                Vector12d g_dx = Vector12d::Zero();
                                DIS::g_PT(P, A, B, C, g_dx); // 2
                                Vector12d grad = g_dx * g_bd;

                                Eigen::Vector3d FA = { grad[3] , grad[4] , grad[5] };
                                Eigen::Vector3d FB = { grad[6] , grad[7] , grad[8] };
                                Eigen::Vector3d FC = { grad[9] , grad[10] , grad[11] };
                                forceObj2[tri[0]] += FA;
                                forceObj2[tri[1]] += FB;
                                forceObj2[tri[2]] += FC;

                            }
                            break;

                            }

						}
					}

				}
			}
		}

	}

    // PT contact force: obj1(T) -> OBJ2(P)
    for (int v = 0; v < obj2.vertices.size(); v++)
    {
        Eigen::Vector3d P = obj2.vertices[v];
        if (insideBoundingBox(P, bbx_min, bbx_max))
        {
            for (int f = 0; f < obj1.faces.size(); f++)
            {
                Eigen::Vector3i tri = obj1.faces[f];
                Eigen::Vector3d A = obj1.vertices[tri[0]];
                Eigen::Vector3d B = obj1.vertices[tri[1]];
                Eigen::Vector3d C = obj1.vertices[tri[2]];
                if (insideBoundingBox(A, bbx_min, bbx_max) && insideBoundingBox(B, bbx_min, bbx_max) && insideBoundingBox(C, bbx_min, bbx_max))
                {
                    if (pointTriangleCCDBroadphase(P, A, B, C, parameters.IPC_dis))
                    {
                        int type = DIS::dType_PT(P, A, B, C);
                        double dis2 = 0;
                        DIS::computePointTriD(P, A, B, C, dis2);

                        if (dis2 <= squaredDouble(parameters.IPC_dis)) // only calculate the energy when the distance is smaller than the threshold
                        {
                            double ipc_dis2 = parameters.IPC_dis * parameters.IPC_dis;
                            double g_bd = BarrierEnergy::compute_g_b(dis2, ipc_dis2);
                            switch (type)
                            {
                            case 0:
                            {
                                Vector6d g_dx = Vector6d::Zero();
                                DIS::g_PP(P, A, g_dx); // 2
                                Vector6d grad = g_dx * g_bd;

                                Eigen::Vector3d FP = { grad[0] , grad[1] , grad[2] };
                                forceObj2[v] += FP;
                            }
                            break;

                            case 1:
                            {
                                Vector6d g_dx = Vector6d::Zero();
                                DIS::g_PP(P, B, g_dx); // 2
                                Vector6d grad = g_dx * g_bd;

                                Eigen::Vector3d FP = { grad[0] , grad[1] , grad[2] };
                                forceObj2[v] += FP;
                            }
                            break;

                            case 2:
                            {
                                Vector6d g_dx = Vector6d::Zero();
                                DIS::g_PP(P, C, g_dx); // 2
                                Vector6d grad = g_dx * g_bd;

                                Eigen::Vector3d FP = { grad[0] , grad[1] , grad[2] };
                                forceObj2[v] += FP;
                            }
                            break;

                            case 3:
                            {
                                Vector9d g_dx = Vector9d::Zero();
                                DIS::g_PE(P, A, B, g_dx); // 2
                                Vector9d grad = g_dx * g_bd;

                                Eigen::Vector3d FP = { grad[0] , grad[1] , grad[2] };
                                forceObj2[v] += FP;
                            }
                            break;

                            case 4:
                            {
                                Vector9d g_dx = Vector9d::Zero();
                                DIS::g_PE(P, B, C, g_dx); // 2
                                Vector9d grad = g_dx * g_bd;

                                Eigen::Vector3d FP = { grad[0] , grad[1] , grad[2] };
                                forceObj2[v] += FP;
                            }
                            break;

                            case 5:
                            {
                                Vector9d g_dx = Vector9d::Zero();
                                DIS::g_PE(P, C, A, g_dx); // 2
                                Vector9d grad = g_dx * g_bd;

                                Eigen::Vector3d FP = { grad[0] , grad[1] , grad[2] };
                                forceObj2[v] += FP;

                            }
                            break;

                            case 6:
                            {
                                Vector12d g_dx = Vector12d::Zero();
                                DIS::g_PT(P, A, B, C, g_dx); // 2
                                Vector12d grad = g_dx * g_bd;

                                Eigen::Vector3d FP = { grad[0] , grad[1] , grad[2] };
                                forceObj2[v] += FP;

                            }
                            break;

                            }

                        }
                    }

                }
            }
        }

    }

    // EE contact force: obj1(E) <-> OBJ2(E)
    for (int e1 = 0; e1 < obj1.edges.size(); e1++)
    {
        int E1_1 = obj1.edges[e1][0], E1_2 = obj1.edges[e1][1];
        Eigen::Vector3d P1 = obj1.vertices[E1_1];
        Eigen::Vector3d P2 = obj1.vertices[E1_2];
        if (insideBoundingBox(P1, bbx_min, bbx_max) && insideBoundingBox(P2, bbx_min, bbx_max))
        {
            for (int e2 = 0; e2 < obj2.edges.size(); e2++)
            {
                int E2_1 = obj2.edges[e2][0], E2_2 = obj2.edges[e2][1];
                Eigen::Vector3d Q1 = obj2.vertices[E2_1];
                Eigen::Vector3d Q2 = obj2.vertices[E2_2];
                if (insideBoundingBox(Q1, bbx_min, bbx_max) && insideBoundingBox(Q2, bbx_min, bbx_max))
                {
                    if (edgeEdgeCCDBroadphase(P1, P2, Q1, Q2, parameters.IPC_dis))
                    {
                        int type = DIS::dType_EE(P1, P2, Q1, Q2);
                        double dis2 = 0;
                        DIS::computeEdgeEdgeD(P1, P2, Q1, Q2, dis2);

                        if (dis2 <= squaredDouble(parameters.IPC_dis)) // only calculate the energy when the distance is smaller than the threshold
                        {
                            double ipc_dis2 = parameters.IPC_dis * parameters.IPC_dis;
                            double g_bd = BarrierEnergy::compute_g_b(dis2, ipc_dis2);
                            

                            switch (type)
                            {
                            case 0:
                            {
                                Vector6d g_dx = Vector6d::Zero(), grad_ = Vector6d::Zero();
                                DIS::g_PP(P1, Q1, g_dx); // 2
                                grad_ = g_dx * g_bd;


                                Eigen::Vector3d FQ1 = { grad_[3] , grad_[4] , grad_[5] };
                                forceObj2[E2_1] += FQ1;

                            }
                            break;

                            case 1:
                            {
                                Vector6d g_dx = Vector6d::Zero(), grad_ = Vector6d::Zero();
                                DIS::g_PP(P1, Q2, g_dx); // 2
                                grad_ = g_dx * g_bd;

                                Eigen::Vector3d FQ2 = { grad_[3] , grad_[4] , grad_[5] };
                                forceObj2[E2_2] += FQ2;

    
                            }
                            break;


                            case 2:
                            {
                                Vector9d g_dx = Vector9d::Zero(), grad_ = Vector9d::Zero();
                                DIS::g_PE(P1, Q1, Q2, g_dx); // 2
                                grad_ = g_dx * g_bd;

                                Eigen::Vector3d FQ1 = { grad_[3] , grad_[4] , grad_[5] };
                                forceObj2[E2_1] += FQ1;
                                Eigen::Vector3d FQ2 = { grad_[6] , grad_[7] , grad_[8] };
                                forceObj2[E2_2] += FQ2;

                            }
                            break;

                            case 3:
                            {
                                Vector6d g_dx = Vector6d::Zero(), grad_ = Vector6d::Zero();
                                DIS::g_PP(P2, Q1, g_dx); // 2
                                grad_ = g_dx * g_bd;


                                Eigen::Vector3d FQ1 = { grad_[3] , grad_[4] , grad_[5] };
                                forceObj2[E2_1] += FQ1;

                            }
                            break;


                            case 4:
                            {
                                Vector6d g_dx = Vector6d::Zero(), grad_ = Vector6d::Zero();
                                DIS::g_PP(P2, Q2, g_dx); // 2
                                grad_ = g_dx * g_bd;

                                Eigen::Vector3d FQ2 = { grad_[3] , grad_[4] , grad_[5] };
                                forceObj2[E2_2] += FQ2;

                            }
                            break;


                            case 5:
                            {
                                Vector9d g_dx = Vector9d::Zero(), grad_ = Vector9d::Zero();
                                DIS::g_PE(P2, Q1, Q2, g_dx); // 2
                                grad_ = g_dx * g_bd;

                                Eigen::Vector3d FQ1 = { grad_[3] , grad_[4] , grad_[5] };
                                forceObj2[E2_1] += FQ1;
                                Eigen::Vector3d FQ2 = { grad_[6] , grad_[7] , grad_[8] };
                                forceObj2[E2_2] += FQ2;

                            }
                            break;


                            case 6:
                            {
                                Vector9d g_dx = Vector9d::Zero(), grad_ = Vector9d::Zero();
                                DIS::g_PE(Q1, P1, P2, g_dx); // 2
                                grad_ = g_dx * g_bd;

                                Eigen::Vector3d FQ1 = { grad_[0] , grad_[1] , grad_[2] };
                                forceObj2[E2_1] += FQ1;

                            }
                            break;


                            case 7:
                            {
                                Vector9d g_dx = Vector9d::Zero(), grad_ = Vector9d::Zero();
                                DIS::g_PE(Q2, P1, P2, g_dx); // 2
                                grad_ = g_dx * g_bd;

                                Eigen::Vector3d FQ2 = { grad_[0] , grad_[1] , grad_[2] };
                                forceObj2[E2_2] += FQ2;

                            }
                            break;


                            case 8:
                            {
                                Vector12d g_dx = Vector12d::Zero(), grad_ = Vector12d::Zero();
                                DIS::g_EE(P1, P2, Q1, Q2, g_dx); // 2
                                grad_ = g_dx * g_bd;


                                Eigen::Vector3d FQ1 = { grad_[6] , grad_[7] , grad_[8] };
                                forceObj2[E2_1] += FQ1;
                                Eigen::Vector3d FQ2 = { grad_[9] , grad_[10] , grad_[11] };
                                forceObj2[E2_2] += FQ2;

                            }
                            break;


                            }


                        }
                    }
                }
            }
        }
    }

    Eigen::Vector3d totalForce = Eigen::Vector3d::Zero();
    for (int i = 0; i < forceObj2.size(); i++)
    {
        totalForce += forceObj2[i];
    }

    return totalForce;
}



