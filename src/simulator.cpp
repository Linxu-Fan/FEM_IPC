#include "simulator.h"


// implicit integration
void implicitFEM(Mesh& tetMesh, FEMParamters& parameters)
{

	for (int timestep = 0; timestep < parameters.num_timesteps; timestep++)
	{
		std::cout << "timestep = " << timestep << std::endl;

		if (timestep % parameters.outputFrequency == 0)
		{
			tetMesh.exportSurfaceMesh("surfMesh", timestep);
		}

		tetMesh.pos_node_prev = tetMesh.pos_node;
		std::vector<Eigen::Vector3d> currentPosition = tetMesh.pos_node;
		double lastEnergyVal = compute_IP_energy(tetMesh, parameters, timestep);
		std::cout << "	lastEnergyVal = " << lastEnergyVal << std::endl;
		for(int ite = 0; ite < 50; ite++)
		{	
					
			std::cout << "		ite = " << ite << std::endl;
			
			std::vector<Eigen::Vector3d> direction = solve_linear_system(tetMesh, parameters, timestep);
			double dist_to_converge = infiniteNorm(direction);
			if (ite && (dist_to_converge < sqrt(parameters.searchResidual * tetMesh.calBBXDiagSize() * parameters.dt * parameters.dt)))
			{
				break;
			}
			
			std::cout << "		Calculate step size"  << std::endl;
			double step = calMaxStepSize(tetMesh, parameters, timestep, direction);
			std::cout << "		Step forward" << std::endl;
			step_forward(parameters,tetMesh, currentPosition, direction, step);
			double newEnergyVal = compute_IP_energy(tetMesh, parameters, timestep);
			std::cout << std::scientific << std::setprecision(4) << "		step = " << step<<"; newEnergyVal = "<< newEnergyVal << "; dist_to_converge = " << dist_to_converge << "; threshold = " << sqrt(parameters.searchResidual * tetMesh.calBBXDiagSize() * parameters.dt * parameters.dt) << std::endl;
			while (newEnergyVal >= lastEnergyVal && step >= 1.0e-7)
			{
				step /= 2.0;
				std::cout << "			step = " << step << std::endl;
				step_forward(parameters, tetMesh, currentPosition, direction, step);
				newEnergyVal = compute_IP_energy(tetMesh, parameters, timestep);
			}
			currentPosition = tetMesh.pos_node;
			lastEnergyVal = newEnergyVal;

			std::cout << "		lastEnergyVal = " << lastEnergyVal << std::endl;



		}


		// update the velocity of the node
		for (int i = 0; i < tetMesh.pos_node.size(); i++)
		{
			tetMesh.vel_node[i] = (tetMesh.pos_node[i] - tetMesh.pos_node_prev[i]) / parameters.dt;
		}

	}

}

// compute external force excluding gravity
Eigen::Vector3d compute_external_force(Mesh& tetMesh, int vertInd, int timestep)
{
	Eigen::Vector3d extForce = Eigen::Vector3d::Zero();


	if (tetMesh.boundaryCondition_node[vertInd].type == 2)
	{
		extForce = tetMesh.boundaryCondition_node[vertInd].force;
	}

	return extForce;
}

// compute the incremental potential energy
double compute_IP_energy(Mesh& tetMesh, FEMParamters& parameters, int timestep)
{
	double energyVal = 0;
	tetMesh.update_F(parameters.numOfThreads);
	
	// energy contribution per vertex
	std::vector<double> node_ext_ine_energy_vec(tetMesh.pos_node.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int vI = 0; vI < tetMesh.pos_node.size(); vI++)
	{
		double nodeMass = tetMesh.mass_node[vI];
		Eigen::Vector3d xt = tetMesh.pos_node_prev[vI];
		Eigen::Vector3d v = tetMesh.vel_node[vI];
		Eigen::Vector3d x = tetMesh.pos_node[vI];
		Eigen::Vector3d extForce = compute_external_force(tetMesh, vI, timestep);

		// the external energy contribution
		double eng = 0;
		eng += ExternalEnergy::Val(nodeMass, parameters.dt, x, parameters, extForce);
		// the inertia energy contribution
		eng += InertiaEnergy::Val(nodeMass, parameters.dt, xt, v, x, extForce, parameters);
		node_ext_ine_energy_vec[vI] = eng;
	}
	energyVal += std::accumulate(node_ext_ine_energy_vec.begin(), node_ext_ine_energy_vec.end(), 0);


	// energy contribution per element
	std::vector<double> tex_est_energy_vec(tetMesh.tetra_F.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int eI = 0; eI < tetMesh.tetra_F.size(); eI++)
	{
		// the internal elastic energy contribution
		int matInd = tetMesh.materialInd[eI];
		tex_est_energy_vec[eI] = ElasticEnergy::Val(tetMesh.materialMesh[matInd], parameters.model, tetMesh.tetra_F[eI], parameters.dt, tetMesh.tetra_vol[eI]);
	}	
	energyVal += std::accumulate(tex_est_energy_vec.begin(), tex_est_energy_vec.end(), 0);


	// energy contribution from barrier
	energyVal += compute_Barrier_energy(tetMesh, parameters, timestep);



//	// ground barrier
//// !!!!!!!!!!!!!!!!!Cannot use the following code directly because the vector didn't consider the size of ground contact
//	if (parameters.enableGround == true)
//	{
//		std::vector<double> bov_ground_energy_vec(tetMesh.boundaryVertices_vec.size());
//#pragma omp parallel for num_threads(parameters.numOfThreads)
//		for (int ft = 0; ft < tetMesh.boundaryVertices_vec.size(); ft++)
//		{
//			int ptInd = tetMesh.boundaryVertices_vec[ft];
//			Eigen::Vector3d P = tetMesh.pos_node[ptInd];
//			double eng = 0;
//			if (P[2] <= parameters.IPC_dis)
//			{
//				eng = Ground::val(P[2] * P[2], parameters.IPC_dis * parameters.IPC_dis, tetMesh.boundaryVertices_area[ptInd], parameters.IPC_kStiffness, parameters.dt);
//			}
//			bov_ground_energy_vec[ft] = eng;
//		}
//		energyVal = std::accumulate(bov_ground_energy_vec.begin(), bov_ground_energy_vec.end(), energyVal);
//	}

	




	return energyVal;
}


double compute_Barrier_energy(Mesh& tetMesh, FEMParamters& parameters, int timestep)
{
	std::vector<spatialHashCellData> spatialHash_vec;
	initSpatialHash(parameters, tetMesh, parameters.IPC_hashSize, spatialHash_vec);

	std::vector<double> energyValue_perHash(spatialHash_vec.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int y = 0; y < spatialHash_vec.size(); y++)
	{
		double energy = 0.0;
		// calcualte the PT pair
		for (std::set<int>::iterator itP = spatialHash_vec[y].vertIndices.begin(); itP != spatialHash_vec[y].vertIndices.end(); itP++)
		{
			for (std::set<int>::iterator itT = spatialHash_vec[y].triaIndices.begin(); itT != spatialHash_vec[y].triaIndices.end(); itT++)
			{
				int vert = *itP, tri = *itT;
				if (tetMesh.boundaryVertices[vert].find(tri) == tetMesh.boundaryVertices[vert].end()) // this triangle is not incident with the point
				{
					Eigen::Vector3i triVerts = tetMesh.boundaryTriangles[tri];
					Eigen::Vector3d P = tetMesh.pos_node[vert];
					Eigen::Vector3d A = tetMesh.pos_node[triVerts[0]];
					Eigen::Vector3d B = tetMesh.pos_node[triVerts[1]];
					Eigen::Vector3d C = tetMesh.pos_node[triVerts[2]];

					if (pointTriangleCCDBroadphase(P, A, B, C, parameters.IPC_dis))
					{
						int type = DIS::dType_PT(P, A, B, C);
						double dis2 = 0;
						DIS::computePointTriD(P, A, B, C, dis2);

						if (dis2 <= squaredDouble(parameters.IPC_dis)) // only calculate the energy when the distance is smaller than the threshold
						{
							energy += BarrierEnergy::val_PT(tetMesh.boundaryVertices_area[vert], dis2, squaredDouble(parameters.IPC_dis), parameters.IPC_kStiffness, parameters.dt);
						}
					}
				}
			}
		}

		// calcualte the EE pair
		for (std::set<int>::iterator itE1 = spatialHash_vec[y].edgeIndices.begin(); itE1 != spatialHash_vec[y].edgeIndices.end(); itE1++)
		{
			for (std::set<int>::iterator itE2 = spatialHash_vec[y].edgeIndices.begin(); itE2 != spatialHash_vec[y].edgeIndices.end(); itE2++)
			{
				if (*itE1 != *itE2)
				{
					Eigen::Vector2i E1 = tetMesh.index_boundaryEdge[*itE1], E2 = tetMesh.index_boundaryEdge[*itE2];
					int P1I = E1[0], P2I = E1[1], Q1I = E2[0], Q2I = E2[1];
					if (P1I != Q1I && P1I != Q2I && P2I != Q1I && P2I != Q2I) // not duplicated and incident edges
					{
						Eigen::Vector3d P1 = tetMesh.pos_node[P1I];
						Eigen::Vector3d P2 = tetMesh.pos_node[P2I];
						Eigen::Vector3d Q1 = tetMesh.pos_node[Q1I];
						Eigen::Vector3d Q2 = tetMesh.pos_node[Q2I];

						if (edgeEdgeCCDBroadphase(P1, P2, Q1, Q2, parameters.IPC_dis))
						{
							int type = DIS::dType_EE(P1, P2, Q1, Q2);
							double dis2 = 0;
							DIS::computeEdgeEdgeD(P1, P2, Q1, Q2, dis2);

							if (dis2 <= squaredDouble(parameters.IPC_dis)) // only calculate the energy when the distance is smaller than the threshold
							{
								Eigen::Vector4i ptIndices = { P1I , P2I , Q1I , Q2I };
								energy += BarrierEnergy::val_EE(tetMesh.boundaryEdges_area[P1I][P2I], dis2, tetMesh, ptIndices, squaredDouble(parameters.IPC_dis), parameters.IPC_kStiffness, parameters.dt);
							}
						}
					}
				}
			}
		}

		energyValue_perHash[y] = energy;
	}

	return std::accumulate(energyValue_perHash.begin(), energyValue_perHash.end(), 0);
}


// compute energy gradient and hessian of the linear system and solve the syetem
std::vector<Eigen::Vector3d> solve_linear_system(Mesh& tetMesh, FEMParamters& parameters, int timestep)
{
	std::vector<Eigen::Vector3d> movingDir(tetMesh.pos_node.size());

	// update defromation gradient
	tetMesh.update_F(parameters.numOfThreads);

	std::cout << "	Calculate contact!" << std::endl;

	BarrierEnergyRes pTeEBarrVec;
	pTeEBarrVec.grad_triplet_vec.resize(tetMesh.pos_node.size() + tetMesh.tetrahedrals.size());
	pTeEBarrVec.hessian_triplet_vec.resize(tetMesh.pos_node.size() + tetMesh.tetrahedrals.size());

	std::cout << "	Calculate ine_ext!" << std::endl;

	// energy contribution per vertex
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int vI = 0; vI < tetMesh.pos_node.size(); vI++)
	{
		double nodeMass = tetMesh.mass_node[vI];
		Eigen::Vector3d xt = tetMesh.pos_node_prev[vI];
		Eigen::Vector3d v = tetMesh.vel_node[vI];
		Eigen::Vector3d x = tetMesh.pos_node[vI];
		Eigen::Vector3d extForce = compute_external_force(tetMesh, vI, timestep);


		// the inertia energy contribution
		std::vector<std::pair<int, double>> inerEngGrad = InertiaEnergy::Grad(nodeMass, parameters.dt, xt, v, x, extForce, vI, parameters, tetMesh.boundaryCondition_node[vI].type);
		std::vector<Eigen::Triplet<double>> inerEngHess = InertiaEnergy::Hess(nodeMass, vI, tetMesh.boundaryCondition_node[vI].type);
		
		// the external energy contribution
		std::vector<std::pair<int, double>> extEngGrad = ExternalEnergy::Grad(nodeMass, parameters.dt, x, parameters, extForce, vI, tetMesh.boundaryCondition_node[vI].type);
		extEngGrad.insert(extEngGrad.end(), inerEngGrad.begin(), inerEngGrad.end());

		pTeEBarrVec.grad_triplet_vec[vI] = extEngGrad;
		pTeEBarrVec.hessian_triplet_vec[vI] = inerEngHess;

	}


	std::cout << "	Calculate elas!" << std::endl;

	int startIndex = tetMesh.pos_node.size();
	// energy contribution per element
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int eI = 0; eI < tetMesh.tetrahedrals.size(); eI++)
	{
		Eigen::Vector4i tetVertInd_BC;
		for (int hg = 0; hg < 4; hg++)
		{
			int vertInd = tetMesh.tetrahedrals[eI][hg];
			tetVertInd_BC[hg] = tetMesh.boundaryCondition_node[vertInd].type;
		}

		// the internal elastic energy contribution
		int matInd = tetMesh.materialInd[eI];
		Eigen::Matrix<double, 9, 12> dFdx = ElasticEnergy::dF_wrt_dx(tetMesh.tetra_DM_inv[eI]);
		std::vector<std::pair<int, double>> elasEngGrad = ElasticEnergy::Grad(tetMesh.materialMesh[matInd], parameters.model, tetMesh.tetra_F[eI], parameters.dt, tetMesh.tetra_vol[eI], dFdx, tetMesh.tetrahedrals[eI], tetVertInd_BC);
		std::vector<Eigen::Triplet<double>> elasEngHess = ElasticEnergy::Hess(tetMesh.materialMesh[matInd], parameters.model, tetMesh.tetra_F[eI], parameters.dt, tetMesh.tetra_vol[eI], dFdx, tetMesh.tetrahedrals[eI], tetVertInd_BC);
		
		pTeEBarrVec.grad_triplet_vec[startIndex + eI] = elasEngGrad;
		pTeEBarrVec.hessian_triplet_vec[startIndex + eI] = elasEngHess;
	}


	calContactInfo(tetMesh, parameters, timestep, pTeEBarrVec);


	std::cout << "	Vectorize!" << std::endl;

	// hessian is the left-hand side, and grad is the right-hand side
	std::vector<Eigen::Triplet<double>> hessian_triplet;
	std::vector<std::pair<int, double>> grad_triplet;
	for (int s = 0; s < pTeEBarrVec.grad_triplet_vec.size(); s++)
	{		
		grad_triplet.insert(grad_triplet.end(), pTeEBarrVec.grad_triplet_vec[s].begin(), pTeEBarrVec.grad_triplet_vec[s].end());
		hessian_triplet.insert(hessian_triplet.end(), pTeEBarrVec.hessian_triplet_vec[s].begin(), pTeEBarrVec.hessian_triplet_vec[s].end());
	}


	std::cout << "	Solve!" << std::endl;

	// assemable the left-hand side 
	Eigen::SparseMatrix<double> leftHandSide(3 * tetMesh.pos_node.size(), 3 * tetMesh.pos_node.size());
	leftHandSide.setFromTriplets(hessian_triplet.begin(), hessian_triplet.end());

	// assemable the right-hand side 
	Eigen::VectorXd rightHandSide = Eigen::VectorXd(tetMesh.pos_node.size() * 3);
	rightHandSide.setZero();
	for (int i = 0; i < grad_triplet.size(); i++)
	{	
		std::pair<int, double> ele = grad_triplet[i];
		rightHandSide[ele.first] += ele.second;
	}

	
	// apply the fixed boundary condition
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int vI = 0; vI < tetMesh.pos_node.size(); vI++)
	{
		if (tetMesh.boundaryCondition_node[vI].type == 1) // fixed points
		{
			// modify the hessian's diagonal elements
			for (int row = vI * 3; row < vI * 3 + 3; row++)
			{
				leftHandSide.coeffRef(row, row) = 1.0;
			}
		}
	}
	

	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
	solver.compute(leftHandSide);
	Eigen::VectorXd result = solver.solve(-rightHandSide);
	
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int i = 0; i < tetMesh.pos_node.size(); i++)
	{
		movingDir[i] = result.block<3, 1>(3 * i, 0);
	}

	return movingDir;
}

// move points' position according to the direction; Note this is a trial movement
void step_forward(FEMParamters& parameters, Mesh& tetMesh, std::vector<Eigen::Vector3d>& currentPosition, std::vector<Eigen::Vector3d>& direction, double step)
{
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int vI = 0; vI < tetMesh.pos_node.size(); vI++)
	{
		tetMesh.pos_node[vI] = currentPosition[vI] + step * direction[vI];
	}
}

void calContactInfo(Mesh& tetMesh, FEMParamters& parameters, int timestep, BarrierEnergyRes& pTeEBarrVec)
{
	std::vector<spatialHashCellData> spatialHash_vec;
	initSpatialHash(parameters, tetMesh, parameters.IPC_hashSize, spatialHash_vec);

	std::vector<double> energyValue_perHash(spatialHash_vec.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int y = 0; y < spatialHash_vec.size(); y++)
	{
		std::vector<Eigen::Triplet<double>> hessian_triplet;
		std::vector<std::pair<int, double>> grad_triplet;

		// calcualte the PT pair
		for (std::set<int>::iterator itP = spatialHash_vec[y].vertIndices.begin(); itP != spatialHash_vec[y].vertIndices.end(); itP++)
		{
			for (std::set<int>::iterator itT = spatialHash_vec[y].triaIndices.begin(); itT != spatialHash_vec[y].triaIndices.end(); itT++)
			{
				int vert = *itP, tri = *itT;
				if (tetMesh.boundaryVertices[vert].find(tri) == tetMesh.boundaryVertices[vert].end()) // this triangle is not incident with the point
				{
					Eigen::Vector3i triVerts = tetMesh.boundaryTriangles[tri];
					Eigen::Vector3d P = tetMesh.pos_node[vert];
					Eigen::Vector3d A = tetMesh.pos_node[triVerts[0]];
					Eigen::Vector3d B = tetMesh.pos_node[triVerts[1]];
					Eigen::Vector3d C = tetMesh.pos_node[triVerts[2]];

					if (pointTriangleCCDBroadphase(P, A, B, C, parameters.IPC_dis))
					{
						int type = DIS::dType_PT(P, A, B, C);
						double dis2 = 0;
						DIS::computePointTriD(P, A, B, C, dis2);

						if (dis2 <= squaredDouble(parameters.IPC_dis)) // only calculate the energy when the distance is smaller than the threshold
						{
							Eigen::Vector4i ptIndices = { vert , triVerts[0] , triVerts[1] , triVerts[2] };
							BarrierEnergy::gradAndHess_PT(hessian_triplet, grad_triplet, tetMesh.boundaryCondition_node, ptIndices, type, dis2, tetMesh, squaredDouble(parameters.IPC_dis), parameters.IPC_kStiffness, parameters.dt);
						}
					}
				}
			}
		}

		// calcualte the EE pair
		for (std::set<int>::iterator itE1 = spatialHash_vec[y].edgeIndices.begin(); itE1 != spatialHash_vec[y].edgeIndices.end(); itE1++)
		{
			for (std::set<int>::iterator itE2 = spatialHash_vec[y].edgeIndices.begin(); itE2 != spatialHash_vec[y].edgeIndices.end(); itE2++)
			{
				if (*itE1 != *itE2)
				{
					Eigen::Vector2i E1 = tetMesh.index_boundaryEdge[*itE1], E2 = tetMesh.index_boundaryEdge[*itE2];
					int P1I = E1[0], P2I = E1[1], Q1I = E2[0], Q2I = E2[1];
					if (P1I != Q1I && P1I != Q2I && P2I != Q1I && P2I != Q2I) // not duplicated and incident edges
					{
						Eigen::Vector3d P1 = tetMesh.pos_node[P1I];
						Eigen::Vector3d P2 = tetMesh.pos_node[P2I];
						Eigen::Vector3d Q1 = tetMesh.pos_node[Q1I];
						Eigen::Vector3d Q2 = tetMesh.pos_node[Q2I];

						if (edgeEdgeCCDBroadphase(P1, P2, Q1, Q2, parameters.IPC_dis))
						{
							int type = DIS::dType_EE(P1, P2, Q1, Q2);
							double dis2 = 0;
							DIS::computeEdgeEdgeD(P1, P2, Q1, Q2, dis2);

							if (dis2 <= squaredDouble(parameters.IPC_dis)) // only calculate the energy when the distance is smaller than the threshold
							{
								Eigen::Vector4i ptIndices = { P1I , P2I , Q1I , Q2I };
								BarrierEnergy::gradAndHess_EE(hessian_triplet, grad_triplet, tetMesh.boundaryCondition_node, ptIndices, type, dis2, tetMesh, squaredDouble(parameters.IPC_dis), parameters.IPC_kStiffness, parameters.dt);
							}
						}
					}
				}
			}
		}

		energyValue_perHash[y] = energy;
	}





























	// energy contribution from barrier
	// point-triangle barrier
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for(int ft = 0; ft < tetMesh.boundaryVertices_vec.size(); ft++)
	{
		int ptInd = tetMesh.boundaryVertices_vec[ft];
		Eigen::Vector3d P = tetMesh.pos_node[ptInd];

		std::vector<Eigen::Triplet<double>> hessian_triplet;
		std::vector<std::pair<int, double>> grad_triplet;
		for (int tI = 0; tI < tetMesh.boundaryTriangles.size(); tI++)
		{
			if (tetMesh.boundaryVertices[ptInd].find(tI) == tetMesh.boundaryVertices[ptInd].end()) // this triangle is not incident with the point
			{
				Eigen::Vector3i tri = tetMesh.boundaryTriangles[tI];
				Eigen::Vector3d A = tetMesh.pos_node[tri[0]];
				Eigen::Vector3d B = tetMesh.pos_node[tri[1]];
				Eigen::Vector3d C = tetMesh.pos_node[tri[2]];


				int type = DIS::dType_PT(P, A, B, C);
				double dis2 =0;
				DIS::computePointTriD(P, A, B, C, dis2);
				
				if (dis2 <= squaredDouble(parameters.IPC_dis)) // only calculate the energy when the distance is smaller than the threshold
				{
					Eigen::Vector4i ptIndices = { ptInd , tri[0] , tri[1] , tri[2] };
					BarrierEnergy::gradAndHess_PT(hessian_triplet, grad_triplet, tetMesh.boundaryCondition_node, ptIndices, type, dis2, tetMesh, parameters.IPC_dis * parameters.IPC_dis, parameters.IPC_kStiffness, parameters.dt);
				}
				
			}
		}
		pTeEBarrVec.hessian_triplet_vec[ft] = hessian_triplet;
		pTeEBarrVec.grad_triplet_vec[ft] = grad_triplet;
	}


	// edge-edge barrier
	int startIndex = tetMesh.boundaryVertices_vec.size();
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int ft = 0; ft < tetMesh.index_boundaryEdge_vec.size(); ft++)
	{
		int edge1 = tetMesh.index_boundaryEdge_vec[ft];

		std::vector<Eigen::Triplet<double>> hessian_triplet;
		std::vector<std::pair<int, double>> grad_triplet;
		for (int gt = 0; gt < tetMesh.index_boundaryEdge_vec.size(); gt++)
		{
			int edge2 = tetMesh.index_boundaryEdge_vec[gt];
			if (ft != gt)
			{
				int e1p1 = tetMesh.index_boundaryEdge[edge1][0], e1p2 = tetMesh.index_boundaryEdge[edge1][1], e2p1 = tetMesh.index_boundaryEdge[edge2][0], e2p2 = tetMesh.index_boundaryEdge[edge2][1];
				if (e1p1 != e2p1 && e1p1 != e2p2 && e1p2 != e2p1 && e1p2 != e2p2) // not duplicated and incident edges
				{
					Eigen::Vector3d P1 = tetMesh.pos_node[e1p1];
					Eigen::Vector3d P2 = tetMesh.pos_node[e1p2];
					Eigen::Vector3d Q1 = tetMesh.pos_node[e2p1];
					Eigen::Vector3d Q2 = tetMesh.pos_node[e2p2];


					int type = DIS::dType_EE(P1, P2, Q1, Q2);
					double dis2 = 0;
					DIS::computeEdgeEdgeD(P1, P2, Q1, Q2, dis2);

					if (dis2 <= squaredDouble(parameters.IPC_dis)) // only calculate the energy when the distance is smaller than the threshold
					{
						Eigen::Vector4i ptIndices = { e1p1 , e1p2 , e2p1 , e2p2 };
						BarrierEnergy::gradAndHess_EE(hessian_triplet, grad_triplet, tetMesh.boundaryCondition_node, ptIndices, type, dis2, tetMesh, parameters.IPC_dis * parameters.IPC_dis, parameters.IPC_kStiffness, parameters.dt);
					}
				}
			}
		}
		pTeEBarrVec.hessian_triplet_vec[startIndex + ft] = hessian_triplet;
		pTeEBarrVec.grad_triplet_vec[startIndex + ft] = grad_triplet;
	
	}


	//// ground barrier
	//// !!!!!!!!!!!!!!!!!Cannot use the following code directly because the vector didn't consider the size of ground contact
	//if (parameters.enableGround == true)
	//{
	//	std::vector<Eigen::Triplet<double>> hessian_triplet;
	//	std::vector<std::pair<int, double>> grad_triplet;
	//	for (int ft = 0; ft < tetMesh.boundaryVertices_vec.size(); ft++)
	//	{
	//		int ptInd = tetMesh.boundaryVertices_vec[ft];
	//		Eigen::Vector3d P = tetMesh.pos_node[ptInd];
	//		if (P[2] <= parameters.IPC_dis)
	//		{
	//			Ground::gradAndHess(hessian_triplet, grad_triplet, tetMesh.boundaryCondition_node, ptInd, P[2] * P[2], parameters.IPC_dis * parameters.IPC_dis, tetMesh.boundaryVertices_area[ptInd], parameters.IPC_kStiffness, parameters.dt);
	//		}
	//	}
	//	pTeEBarrVec.grad_triplet_vec.push_back(grad_triplet);
	//	pTeEBarrVec.hessian_triplet_vec.push_back(hessian_triplet);
	//}
	
}

// calculate the maximum feasible step size
double calMaxStepSize(Mesh& tetMesh, FEMParamters& parameters, int timestep, std::vector<Eigen::Vector3d>& direction)
{
	////std::cout << "Calculating maximum step!" << std::endl;
	//std::set<int> culledSet; // vertices who are in a contact

	//// Step 1: calculate the culled constraint
	//for (int i  = 0; i < pTeEBarrVec.PT_Indices.size(); i++)
	//{
	//	culledSet.insert(pTeEBarrVec.PT_Indices[i].begin(), pTeEBarrVec.PT_Indices[i].end());
	//}
	//for (int i = 0; i < pTeEBarrVec.EE_Indices.size(); i++)
	//{
	//	culledSet.insert(pTeEBarrVec.EE_Indices[i].begin(), pTeEBarrVec.EE_Indices[i].end());
	//}
	////std::cout << "culledSet.size() = " << culledSet.size() << std::endl;
	//
	//// Step 2: calculate alpha_F
	//double alpha_F = 1.0; // Eq.3 in IPC paper's supplementary document
	//for (int i = 0; i < direction.size(); i++)
	//{
	//	if (culledSet.find(i) == culledSet.end())
	//	{
	//		alpha_F = std::min(alpha_F, parameters.IPC_dis / 2.0 / direction[i].norm());
	//	}
	//}
	////std::cout << "alpha_F = " << alpha_F <<"; pTeEBarrVec.PT_Indices.size() = "<< pTeEBarrVec.PT_Indices.size() << std::endl;

	//// Step 3: calculate alpha_C_hat
	//double alpha_C_hat = 1.0;
	//for (int i = 0; i < pTeEBarrVec.PT_Indices.size(); i++)
	//{
	//	Eigen::Vector4i PP_Index = pTeEBarrVec.PT_Indices[i];
	//	int p0 = PP_Index[0], p1 = PP_Index[1], p2 = PP_Index[2], p3 = PP_Index[3];
	//	//std::cout << "	PT " << i << " ";
	//	bool intersect = pointTriangleCCDBroadphase(tetMesh.pos_node[p0], direction[p0], tetMesh.pos_node[p1], 
	//		direction[p1], tetMesh.pos_node[p2], direction[p2], tetMesh.pos_node[p3], direction[p3], parameters.IPC_dis);
	//	//std::cout << "	intersect = " << intersect << std::endl;
	//	if (intersect)
	//	{
	//		//std::cout << "		tetMesh.pos_node[p0] = "<< tetMesh.pos_node[p0] << std::endl;
	//		//std::cout << "		tetMesh.pos_node[p1] = "<< tetMesh.pos_node[p1] << std::endl;
	//		//std::cout << "		tetMesh.pos_node[p2] = "<< tetMesh.pos_node[p2] << std::endl;
	//		//std::cout << "		tetMesh.pos_node[p3] = "<< tetMesh.pos_node[p3] << std::endl;
	//		//std::cout << "		direction[p0] = "<< direction[p0] << std::endl;
	//		//std::cout << "		direction[p1] = "<< direction[p1] << std::endl;
	//		//std::cout << "		direction[p2] = "<< direction[p2] << std::endl;
	//		//std::cout << "		direction[p3] = "<< direction[p3] << std::endl;
	//		double alpha_tmp = pointTriangleCCDNarrowphase(tetMesh.pos_node[p0], direction[p0], tetMesh.pos_node[p1],
	//			direction[p1], tetMesh.pos_node[p2], direction[p2], tetMesh.pos_node[p3], direction[p3], parameters.IPC_eta);
	//		//std::cout << "		xxx2" << std::endl;
	//		alpha_C_hat = std::min(alpha_C_hat, alpha_tmp);
	//	}
	//	
	//}
	//for (int i = 0; i < pTeEBarrVec.EE_Indices.size(); i++)
	//{
	//	Eigen::Vector4i PP_Index = pTeEBarrVec.EE_Indices[i];
	//	int p0 = PP_Index[0], p1 = PP_Index[1], p2 = PP_Index[2], p3 = PP_Index[3];
	//	//std::cout << "	EE " << i << " ";
	//	bool intersect = edgeEdgeCCDBroadphase(tetMesh.pos_node[p0], direction[p0], tetMesh.pos_node[p1],
	//		direction[p1], tetMesh.pos_node[p2], direction[p2], tetMesh.pos_node[p3], direction[p3], parameters.IPC_dis);
	//	if (intersect)
	//	{
	//		double alpha_tmp = edgeEdgeCCDNarrowphase(tetMesh.pos_node[p0], direction[p0], tetMesh.pos_node[p1],
	//			direction[p1], tetMesh.pos_node[p2], direction[p2], tetMesh.pos_node[p3], direction[p3], parameters.IPC_eta);
	//		alpha_C_hat = std::min(alpha_C_hat, alpha_tmp);
	//	}
	//	
	//}
	////std::cout << "alpha_C_hat = " << alpha_C_hat << std::endl;

	//// Step 4: calculate full CCD if necessary
	////if (alpha_F >= 0.5 * alpha_C_hat)
	//if(0)
	//{
	//	//std::cout << "Partial calculation!" << std::endl;
	//	if (parameters.enableGround == true)
	//	{
	//		double stepGround = 1.0;
	//		for (std::map<int, std::set<int>>::iterator it = tetMesh.boundaryVertices.begin(); it != tetMesh.boundaryVertices.end(); it++)
	//		{
	//			int ptInd = it->first;
	//			if (direction[ptInd][2] < 0)
	//			{
	//				double coor_z = tetMesh.pos_node[ptInd][2];
	//				stepGround = std::min(stepGround, coor_z * (1.0 - parameters.IPC_eta) / std::abs(direction[ptInd][2]));
	//			}
	//			
	//		}
	//		alpha_C_hat = std::min(stepGround, alpha_C_hat);
	//	}
	//	return std::min(alpha_F, alpha_C_hat);
	//}
	//else // use spatial hash to calculate the actual full CCD
	{
		double CCD_step = calMaxStep_spatialHash(parameters, tetMesh, direction, parameters.IPC_hashSize, parameters.IPC_dis, parameters.IPC_eta);
		if (parameters.enableGround == true)
		{
			double stepGround = 1.0;
			for (std::map<int, std::set<int>>::iterator it = tetMesh.boundaryVertices.begin(); it != tetMesh.boundaryVertices.end(); it++)
			{
				int ptInd = it->first;
				if (direction[ptInd][2] < 0)
				{
					double coor_z = tetMesh.pos_node[ptInd][2];
					stepGround = std::min(stepGround, coor_z * (1.0 - parameters.IPC_eta) / std::abs(direction[ptInd][2]));
				}

			}
			CCD_step = std::min(stepGround, CCD_step);
		}
		return CCD_step;
	}

}


