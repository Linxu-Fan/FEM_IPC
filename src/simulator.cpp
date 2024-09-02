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
		for(int ite = 0; ite < 20; ite++)
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
			while (newEnergyVal > lastEnergyVal && step >= 1.0e-5)
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
	energyVal = std::accumulate(node_ext_ine_energy_vec.begin(), node_ext_ine_energy_vec.end(), 0.0);



	// energy contribution per element
	std::vector<double> tex_est_energy_vec(tetMesh.tetra_F.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int eI = 0; eI < tetMesh.tetra_F.size(); eI++)
	{
		// the internal elastic energy contribution
		int matInd = tetMesh.materialInd[eI];
		tex_est_energy_vec[eI] = ElasticEnergy::Val(tetMesh.materialMesh[matInd], parameters.model, tetMesh.tetra_F[eI], parameters.dt, tetMesh.tetra_vol[eI]);
	}	
	energyVal += std::accumulate(tex_est_energy_vec.begin(), tex_est_energy_vec.end(), 0.0);



	// energy contribution from barrier
	energyVal += compute_Barrier_energy(tetMesh, parameters, timestep);



	return energyVal;
}


double compute_Barrier_energy(Mesh& tetMesh, FEMParamters& parameters, int timestep)
{
	std::vector<spatialHashCellData> spatialHash_vec;
	std::map<std::string, int> hashNameIndex;
	std::vector<Eigen::Vector3d> direction(tetMesh.pos_node.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int i = 0; i < tetMesh.pos_node.size(); i++)
	{
		direction[i] = Eigen::Vector3d::Zero();
	}
	initSpatialHash(false, parameters, tetMesh, direction, parameters.IPC_hashSize, spatialHash_vec, hashNameIndex,timestep);

	

	std::vector<double> energyValue_perHash_PT(spatialHash_vec.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int y = 0; y < spatialHash_vec.size(); y++)
	{
		double energy = 0.0;
		// calcualte the PT pair
		for (std::set<int>::iterator itP = spatialHash_vec[y].vertIndices.begin(); itP != spatialHash_vec[y].vertIndices.end(); itP++)
		{
			Eigen::Vector3i bottomLeftCorner = spatialHash_vec[y].bottomLeftCorner;
			for (int xx = bottomLeftCorner[0] - 1; xx <= bottomLeftCorner[0] + 1; xx++)
			{
				for (int yy = bottomLeftCorner[1] - 1; yy <= bottomLeftCorner[1] + 1; yy++)
				{
					for (int zz = bottomLeftCorner[2] - 1; zz <= bottomLeftCorner[2] + 1; zz++)
					{
						Eigen::Vector3i index = { xx , yy , zz };
						std::string ID = calculateID(index);
						if (hashNameIndex.find(ID) != hashNameIndex.end())
						{
							int neigHashIndex = hashNameIndex[ID];
							for (std::set<int>::iterator itT = spatialHash_vec[neigHashIndex].triaIndices.begin(); itT != spatialHash_vec[neigHashIndex].triaIndices.end(); itT++)
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
					}
				}
			}
		}
		energyValue_perHash_PT[y] = energy;
	}
	double energyHash = std::accumulate(energyValue_perHash_PT.begin(), energyValue_perHash_PT.end(), 0.0);


	std::vector<std::map<std::string,double>> energyValue_perHash_EE(spatialHash_vec.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int y = 0; y < spatialHash_vec.size(); y++)
	{
		std::map<std::string, double> energyValue;
		// calcualte the EE pair
		for (std::set<int>::iterator itE1 = spatialHash_vec[y].edgeIndices.begin(); itE1 != spatialHash_vec[y].edgeIndices.end(); itE1++)
		{
			Eigen::Vector3i bottomLeftCorner = spatialHash_vec[y].bottomLeftCorner;
			for (int xx = bottomLeftCorner[0]; xx <= bottomLeftCorner[0] + 1; xx++)
			{
				for (int yy = bottomLeftCorner[1]; yy <= bottomLeftCorner[1] + 1; yy++)
				{
					for (int zz = bottomLeftCorner[2]; zz <= bottomLeftCorner[2] + 1; zz++)
					{
						Eigen::Vector3i index = { xx , yy , zz };
						std::string ID = calculateID(index);
						if (hashNameIndex.find(ID) != hashNameIndex.end())
						{
							int neigHashIndex = hashNameIndex[ID];
							for (std::set<int>::iterator itE2 = spatialHash_vec[neigHashIndex].edgeIndices.begin(); itE2 != spatialHash_vec[neigHashIndex].edgeIndices.end(); itE2++)
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
												double energy = BarrierEnergy::val_EE(tetMesh.boundaryEdges_area[P1I][P2I], dis2, tetMesh, ptIndices, squaredDouble(parameters.IPC_dis), parameters.IPC_kStiffness, parameters.dt);
												std::string EEPair = std::to_string(std::min(*itE1, *itE2)) + "#" + std::to_string(std::max(*itE1, *itE2));
												energyValue[EEPair] = energy;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		energyValue_perHash_EE[y] = energyValue;
	}


	// add the edge-edge energy value
	std::set<std::string> addedEnergy;
	for (int i = 0; i < energyValue_perHash_EE.size(); i++)
	{
		for (std::map<std::string, double>::iterator it = energyValue_perHash_EE[i].begin(); it != energyValue_perHash_EE[i].end(); it++)
		{
			std::string name = it->first;
			if (addedEnergy.find(name) == addedEnergy.end())
			{
				addedEnergy.insert(name);
				energyHash += it->second;
			}
		}
	}


	// ground barrier
	if (parameters.enableGround == true)
	{
		std::vector<double> bov_ground_energy_vec(tetMesh.boundaryVertices_vec.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int ft = 0; ft < tetMesh.boundaryVertices_vec.size(); ft++)
		{
			int ptInd = tetMesh.boundaryVertices_vec[ft];
			Eigen::Vector3d P = tetMesh.pos_node[ptInd];
			double eng = 0;
			if (P[2] <= parameters.IPC_dis)
			{
				eng = Ground::val(P[2] * P[2], parameters.IPC_dis * parameters.IPC_dis, tetMesh.boundaryVertices_area[ptInd], parameters.IPC_kStiffness, parameters.dt);
			}
			bov_ground_energy_vec[ft] = eng;
		}
		double energyGround = std::accumulate(bov_ground_energy_vec.begin(), bov_ground_energy_vec.end(), 0.0);
		energyHash += energyGround;
	}


	return energyHash;
}


// compute energy gradient and hessian of the linear system and solve the syetem
std::vector<Eigen::Vector3d> solve_linear_system(Mesh& tetMesh, FEMParamters& parameters, int timestep)
{
	std::vector<Eigen::Vector3d> movingDir(tetMesh.pos_node.size());

	// update defromation gradient
	tetMesh.update_F(parameters.numOfThreads);



	std::vector<Vector5i> PG_PG, PT_PP, PT_PE, PT_PT, EE_EE;
	calContactInfo(tetMesh, parameters, timestep, PG_PG, PT_PP, PT_PE, PT_PT, EE_EE);
	std::cout << "EE_EE.size() = " << EE_EE.size() << std::endl;


	int gradSize = tetMesh.pos_node.size() * 6 + tetMesh.tetrahedrals.size() * 12 + PG_PG.size() * 3
					+ PT_PP.size() * 6 + PT_PE.size() * 9 + PT_PT.size() * 12 + EE_EE.size() * 12;
	int hessSize = tetMesh.pos_node.size() * 3 + tetMesh.tetrahedrals.size() * 144 + PG_PG.size() * 9
				   + PT_PP.size() * 36 + PT_PE.size() * 81 + PT_PT.size() * 144 + EE_EE.size() * 144;
	std::vector<std::pair<int, double>> grad_triplet(gradSize);
	std::vector<Eigen::Triplet<double>> hessian_triplet(hessSize);
	



	// energy contribution per vertex
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int vI = 0; vI < tetMesh.pos_node.size(); vI++)
	{
		double nodeMass = tetMesh.mass_node[vI];
		Eigen::Vector3d xt = tetMesh.pos_node_prev[vI];
		Eigen::Vector3d v = tetMesh.vel_node[vI];
		Eigen::Vector3d x = tetMesh.pos_node[vI];
		Eigen::Vector3d extForce = compute_external_force(tetMesh, vI, timestep);

		// !!!! The order is very important. "InertiaEnergy" must be the first and "ExternalEnergy" the second

		// the inertia energy contribution
		int actGradIndex = vI * 6;
		InertiaEnergy::Grad(grad_triplet, actGradIndex, nodeMass, parameters.dt, xt, v, x, extForce, vI, parameters, tetMesh.boundaryCondition_node[vI].type);
		int actHessIndex = vI * 3;
		InertiaEnergy::Hess(hessian_triplet, actHessIndex, nodeMass, vI, tetMesh.boundaryCondition_node[vI].type);
		
		// the external energy contribution
		actGradIndex = vI * 6 + 3;
		ExternalEnergy::Grad(grad_triplet, actGradIndex, nodeMass, parameters.dt, x, parameters, extForce, vI, tetMesh.boundaryCondition_node[vI].type);

	}


	int startIndex_grad = tetMesh.pos_node.size() * 6, startIndex_hess = tetMesh.pos_node.size() * 3;
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
		int actualStartIndex_hess = startIndex_hess + eI * 144, actualStartIndex_grad = startIndex_grad + eI * 12;
		ElasticEnergy::Grad(grad_triplet, actualStartIndex_grad, tetMesh.materialMesh[matInd], parameters.model, tetMesh.tetra_F[eI], parameters.dt, tetMesh.tetra_vol[eI], dFdx, tetMesh.tetrahedrals[eI], tetVertInd_BC);
		ElasticEnergy::Hess(hessian_triplet, actualStartIndex_hess, tetMesh.materialMesh[matInd], parameters.model, tetMesh.tetra_F[eI], parameters.dt, tetMesh.tetra_vol[eI], dFdx, tetMesh.tetrahedrals[eI], tetVertInd_BC);
		
		
	}


	if (gradSize > tetMesh.pos_node.size() * 6 + tetMesh.tetrahedrals.size() * 12)
	{
		std::cout << "Export EE pair" << std::endl;
		std::cout << "PT_PP.size() = " << PT_PP.size() << std::endl;
		std::cout << "PT_PE.size() = " << PT_PE.size() << std::endl;
		std::cout << "PT_PT.size() = " << PT_PT.size() << std::endl;
		std::cout << "EE_EE.size() = " << EE_EE.size() << std::endl;

		std::ofstream outfile2("./output/" + std::to_string(timestep) + ".obj", std::ios::trunc);
		for (int i = 0; i < EE_EE.size(); i++)
		{
			Vector5i cont_EE = EE_EE[i];
			int E1 = cont_EE[1], E2 = cont_EE[2];
			int e1p1 = tetMesh.index_boundaryEdge[E1][0], e1p2 = tetMesh.index_boundaryEdge[E1][1], e2p1 = tetMesh.index_boundaryEdge[E2][0], e2p2 = tetMesh.index_boundaryEdge[E2][1];

			Eigen::Vector3d P1 = tetMesh.pos_node[e1p1];
			Eigen::Vector3d P2 = tetMesh.pos_node[e1p2];
			Eigen::Vector3d Q1 = tetMesh.pos_node[e2p1];
			Eigen::Vector3d Q2 = tetMesh.pos_node[e2p2];

			outfile2 << std::scientific << std::setprecision(8) << "v " << P1[0] << " " << P1[1] << " " << P1[2] << std::endl;
			outfile2 << std::scientific << std::setprecision(8) << "v " << P2[0] << " " << P2[1] << " " << P2[2] << std::endl;
			outfile2 << std::scientific << std::setprecision(8) << "v " << Q1[0] << " " << Q1[1] << " " << Q1[2] << std::endl;
			outfile2 << std::scientific << std::setprecision(8) << "v " << Q2[0] << " " << Q2[1] << " " << Q2[2] << std::endl;

		}


		for (int i = 0; i < EE_EE.size(); i++)
		{
			outfile2 << std::scientific << std::setprecision(8) << "l " << std::to_string(i * 2 + 1) << " " << std::to_string(i * 2 + 1) << std::endl;
			outfile2 << std::scientific << std::setprecision(8) << "l " << std::to_string(i * 2 + 3) << " " << std::to_string(i * 2 + 4) << std::endl;

		}
		outfile2.close();

		std::cout <<  std::endl;
		
	}


	// calculate the barrier gradient and hessian
	{
		startIndex_grad += tetMesh.tetrahedrals.size() * 12, startIndex_hess += tetMesh.tetrahedrals.size() * 144;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < PG_PG.size(); i++)
		{
			int ptInd = PG_PG[i][1];
			Eigen::Vector3d P = tetMesh.pos_node[ptInd];
			int actualStartIndex_hess = startIndex_hess + i * 9, actualStartIndex_grad = startIndex_grad + i * 3;
			Ground::gradAndHess(hessian_triplet, grad_triplet, actualStartIndex_hess, actualStartIndex_grad, tetMesh.boundaryCondition_node, ptInd, P[2] * P[2], parameters.IPC_dis * parameters.IPC_dis, tetMesh.boundaryVertices_area[ptInd], parameters.IPC_kStiffness, parameters.dt);
		}

		startIndex_grad += PG_PG.size() * 3, startIndex_hess += PG_PG.size() * 9;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < PT_PP.size(); i++)
		{
			Vector5i cont_PT = PT_PP[i];

			int ptInd = cont_PT[1];
			Eigen::Vector3d P = tetMesh.pos_node[ptInd];
			Eigen::Vector3i tri = tetMesh.boundaryTriangles[cont_PT[2]];
			Eigen::Vector3d A = tetMesh.pos_node[tri[0]];
			Eigen::Vector3d B = tetMesh.pos_node[tri[1]];
			Eigen::Vector3d C = tetMesh.pos_node[tri[2]];						
			int type = cont_PT[3];

			double dis2 =0;
			DIS::computePointTriD(P, A, B, C, dis2);							
			Eigen::Vector4i ptIndices = { ptInd , tri[0] , tri[1] , tri[2] };
			int actualStartIndex_hess = startIndex_hess + i * 36, actualStartIndex_grad = startIndex_grad + i * 6;
			BarrierEnergy::gradAndHess_PT(hessian_triplet, grad_triplet, actualStartIndex_hess, actualStartIndex_grad, tetMesh.boundaryCondition_node, ptIndices, type, dis2, tetMesh, parameters.IPC_dis * parameters.IPC_dis, parameters.IPC_kStiffness, parameters.dt);
		
		}

		startIndex_grad += PT_PP.size() * 6, startIndex_hess += PT_PP.size() * 36;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < PT_PE.size(); i++)
		{
			Vector5i cont_PT = PT_PE[i];

			int ptInd = cont_PT[1];
			Eigen::Vector3d P = tetMesh.pos_node[ptInd];
			Eigen::Vector3i tri = tetMesh.boundaryTriangles[cont_PT[2]];
			Eigen::Vector3d A = tetMesh.pos_node[tri[0]];
			Eigen::Vector3d B = tetMesh.pos_node[tri[1]];
			Eigen::Vector3d C = tetMesh.pos_node[tri[2]];
			int type = cont_PT[3];

			double dis2 = 0;
			DIS::computePointTriD(P, A, B, C, dis2);
			Eigen::Vector4i ptIndices = { ptInd , tri[0] , tri[1] , tri[2] };
			int actualStartIndex_hess = startIndex_hess + i * 81, actualStartIndex_grad = startIndex_grad + i * 9;
			BarrierEnergy::gradAndHess_PT(hessian_triplet, grad_triplet, actualStartIndex_hess, actualStartIndex_grad, tetMesh.boundaryCondition_node, ptIndices, type, dis2, tetMesh, parameters.IPC_dis * parameters.IPC_dis, parameters.IPC_kStiffness, parameters.dt);

		}

		startIndex_grad += PT_PE.size() * 9, startIndex_hess += PT_PE.size() * 81;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < PT_PT.size(); i++)
		{
			Vector5i cont_PT = PT_PT[i];

			int ptInd = cont_PT[1];
			Eigen::Vector3d P = tetMesh.pos_node[ptInd];
			Eigen::Vector3i tri = tetMesh.boundaryTriangles[cont_PT[2]];
			Eigen::Vector3d A = tetMesh.pos_node[tri[0]];
			Eigen::Vector3d B = tetMesh.pos_node[tri[1]];
			Eigen::Vector3d C = tetMesh.pos_node[tri[2]];
			int type = cont_PT[3];

			double dis2 = 0;
			DIS::computePointTriD(P, A, B, C, dis2);
			Eigen::Vector4i ptIndices = { ptInd , tri[0] , tri[1] , tri[2] };
			int actualStartIndex_hess = startIndex_hess + i * 144, actualStartIndex_grad = startIndex_grad + i * 12;
			BarrierEnergy::gradAndHess_PT(hessian_triplet, grad_triplet, actualStartIndex_hess, actualStartIndex_grad, tetMesh.boundaryCondition_node, ptIndices, type, dis2, tetMesh, parameters.IPC_dis * parameters.IPC_dis, parameters.IPC_kStiffness, parameters.dt);

		}

		startIndex_grad += PT_PT.size() * 12, startIndex_hess += PT_PT.size() * 144;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < EE_EE.size(); i++)
		{
			Vector5i cont_EE = EE_EE[i];
			int E1 = cont_EE[1], E2 = cont_EE[2];
			int e1p1 = tetMesh.index_boundaryEdge[E1][0], e1p2 = tetMesh.index_boundaryEdge[E1][1], e2p1 = tetMesh.index_boundaryEdge[E2][0], e2p2 = tetMesh.index_boundaryEdge[E2][1];

			Eigen::Vector3d P1 = tetMesh.pos_node[e1p1];
			Eigen::Vector3d P2 = tetMesh.pos_node[e1p2];
			Eigen::Vector3d Q1 = tetMesh.pos_node[e2p1];
			Eigen::Vector3d Q2 = tetMesh.pos_node[e2p2];
			int type = cont_EE[3];
						
			double dis2 = 0;
			DIS::computeEdgeEdgeD(P1, P2, Q1, Q2, dis2);
			Eigen::Vector4i ptIndices = { e1p1 , e1p2 , e2p1 , e2p2 };
			int actualStartIndex_hess = startIndex_hess + i * 144, actualStartIndex_grad = startIndex_grad + i * 12;
			BarrierEnergy::gradAndHess_EE(hessian_triplet, grad_triplet, actualStartIndex_hess, actualStartIndex_grad, tetMesh.boundaryCondition_node, ptIndices, type, dis2, tetMesh, parameters.IPC_dis * parameters.IPC_dis, parameters.IPC_kStiffness, parameters.dt);

		}

	}


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


void calContactInfo(Mesh& tetMesh, FEMParamters& parameters, int timestep, std::vector<Vector5i> & PG_PG, std::vector<Vector5i>& PT_PP, std::vector<Vector5i>& PT_PE, std::vector<Vector5i>& PT_PT, std::vector<Vector5i>& EE_EE)
{
	std::vector<spatialHashCellData> spatialHash_vec;
	std::map<std::string, int> hashNameIndex;
	std::vector<Eigen::Vector3d> direction(tetMesh.pos_node.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int i = 0; i < tetMesh.pos_node.size(); i++)
	{
		direction[i] = Eigen::Vector3d::Zero();
	}
	initSpatialHash(false, parameters, tetMesh, direction, parameters.IPC_hashSize, spatialHash_vec, hashNameIndex, timestep);


	// Step 1: find contact pairs and types
	std::vector<std::vector<Vector5i>> cont_pair_vec(spatialHash_vec.size()); // 1st int: 0(PT), 1(EE), 2(PG); 2nd int: index of P(E1)(P: for ground contact case); 3rd int: index of T(E2); 4th int: type; 5th int: actual involved elements, i.e. PP(2), PE(3), PT(4) and EE(4)  
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int y = 0; y < spatialHash_vec.size(); y++)
	{
		std::vector<Vector5i> cont_pair;
		for (std::set<int>::iterator itP = spatialHash_vec[y].vertIndices.begin(); itP != spatialHash_vec[y].vertIndices.end(); itP++)
		{
			Eigen::Vector3i bottomLeftCorner = spatialHash_vec[y].bottomLeftCorner;
			for (int xx = bottomLeftCorner[0] - 1; xx <= bottomLeftCorner[0] + 1; xx++)
			{
				for (int yy = bottomLeftCorner[1] - 1; yy <= bottomLeftCorner[1] + 1; yy++)
				{
					for (int zz = bottomLeftCorner[2] - 1; zz <= bottomLeftCorner[2] + 1; zz++)
					{
						Eigen::Vector3i index = { xx , yy , zz };
						std::string ID = calculateID(index);
						if (hashNameIndex.find(ID) != hashNameIndex.end())
						{
							int neigHashIndex = hashNameIndex[ID];
							for (std::set<int>::iterator itT = spatialHash_vec[neigHashIndex].triaIndices.begin(); itT != spatialHash_vec[neigHashIndex].triaIndices.end(); itT++)
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
											if (type <= 2)
											{
												Vector5i ct = { 0 , vert , tri , type, 2 };
												cont_pair.push_back(ct);
											}
											else if (type > 2 && type <= 5)
											{
												Vector5i ct = { 0 , vert , tri , type, 3 };
												cont_pair.push_back(ct);
											}
											else if (type == 6)
											{
												Vector5i ct = { 0 , vert , tri , type, 4 };
												cont_pair.push_back(ct);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	
		for (std::set<int>::iterator itE1 = spatialHash_vec[y].edgeIndices.begin(); itE1 != spatialHash_vec[y].edgeIndices.end(); itE1++)
		{
			Eigen::Vector3i bottomLeftCorner = spatialHash_vec[y].bottomLeftCorner;
			for (int xx = bottomLeftCorner[0]; xx <= bottomLeftCorner[0] + 1; xx++)
			{
				for (int yy = bottomLeftCorner[1]; yy <= bottomLeftCorner[1] + 1; yy++)
				{
					for (int zz = bottomLeftCorner[2]; zz <= bottomLeftCorner[2] + 1; zz++)
					{
						Eigen::Vector3i index = { xx , yy , zz };
						std::string ID = calculateID(index);
						if (hashNameIndex.find(ID) != hashNameIndex.end())
						{
							int neigHashIndex = hashNameIndex[ID];
							for (std::set<int>::iterator itE2 = spatialHash_vec[neigHashIndex].edgeIndices.begin(); itE2 != spatialHash_vec[neigHashIndex].edgeIndices.end(); itE2++)
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
												Vector5i ct = { 1 , *itE1 , *itE2 , type , 4 };
												cont_pair.push_back(ct);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		cont_pair_vec[y] = cont_pair;
	}


	// Step 2: delete empty pairs
	std::set<std::string> storedEEPairs;
	for (int i = 0; i < cont_pair_vec.size(); i++)
	{
		for (int j = 0; j < cont_pair_vec[i].size(); j++)
		{
			if (cont_pair_vec[i][j][0] == 0)
			{
				if (cont_pair_vec[i][j][4] == 2)
				{
					PT_PP.push_back(cont_pair_vec[i][j]);
				}
				else if (cont_pair_vec[i][j][4] == 3)
				{
					PT_PE.push_back(cont_pair_vec[i][j]);
				}
				else
				{
					PT_PT.push_back(cont_pair_vec[i][j]);
				}
			}
			else
			{
				std::string EEPair = std::to_string(std::min(cont_pair_vec[i][j][1], cont_pair_vec[i][j][2])) + "#" + std::to_string(std::max(cont_pair_vec[i][j][1], cont_pair_vec[i][j][2]));
				if (storedEEPairs.find(EEPair) == storedEEPairs.end())
				{
					EE_EE.push_back(cont_pair_vec[i][j]);
					storedEEPairs.insert(EEPair);
				}
			}
		}
	}


	// Step 3: find ground contact pairs if any
	if (parameters.enableGround == true)
	{
		for (int ft = 0; ft < tetMesh.boundaryVertices_vec.size(); ft++)
		{
			int ptInd = tetMesh.boundaryVertices_vec[ft];
			if (tetMesh.pos_node[ptInd][2] <= parameters.IPC_dis)
			{
				Vector5i ct = {2 , ptInd , 0 , 0 , 0};
				PG_PG.push_back(ct);
			}
		}

	}

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


