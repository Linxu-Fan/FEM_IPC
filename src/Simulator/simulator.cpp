#include "simulator.h"

void explicitFEM(Mesh& tetSimMesh, FEMParamters& parameters)
{
	// create a temporary vector for parallelization
	std::vector<std::vector<Eigen::Vector2i>> nodeSharedByElement(tetSimMesh.pos_node.size()); // Eigen::Vector2i 1st int: tetrahedral index; 2nd int: index of this node in this tetrahedral
	// find all elements that share a node
	for (int eleInd = 0; eleInd < tetSimMesh.tetrahedrals.size(); eleInd++)
	{
		for (int j = 0; j < 4; j++)
		{
			int nodeInd = tetSimMesh.tetrahedrals[eleInd][j];
			Eigen::Vector2i indices = { eleInd , j};
			nodeSharedByElement[nodeInd].push_back(indices);
		}
	}
	std::vector<Eigen::Matrix3d> tetrahedralStress(tetSimMesh.tetrahedrals.size(), Eigen::Matrix3d::Zero());

	// Simulation loop
	for (int timestep = 0; timestep < parameters.num_timesteps; ++timestep)
	{
		std::cout << "timestep = " << timestep << std::endl;

		if (timestep % parameters.outputFrequency == 0)
		{
			tetSimMesh.exportSurfaceMesh("surfMesh", timestep);

			{
				std::ofstream outfile9("./output/Stress_" + std::to_string(timestep) + ".txt", std::ios::trunc);
				for (int nf = 0; nf < tetSimMesh.tetrahedrals.size(); nf++)
				{
					int matInd = tetSimMesh.materialInd[nf];
					Material mat = tetSimMesh.materialMesh[matInd];

					// Compute First Piola-Kirchhoff stress tensor P
					Eigen::Matrix3d F = tetSimMesh.tetra_F[nf];
					double J = F.determinant();
					Eigen::Matrix3d FinvT = F.inverse().transpose();
					Eigen::Matrix3d PK1 = mat.mu * (F - FinvT) + mat.lambda * log(J) * FinvT;
					Eigen::Matrix3d cauchyStress = 1.0 / J * PK1 * F.transpose();

					Eigen::EigenSolver<Eigen::MatrixXd> es2(cauchyStress);
					Eigen::Vector3d eigenValues = { es2.eigenvalues()[0].real() ,  es2.eigenvalues()[1].real() ,  es2.eigenvalues()[2].real() };
					Eigen::Matrix3d eigenVectors = es2.eigenvectors().real();
					double maxEigenValue = std::max(std::max(eigenValues[0], eigenValues[1]), eigenValues[2]);

					Eigen::Vector4i tet = tetSimMesh.tetrahedrals[nf];
					Eigen::Vector3d pos = 0.25 * (tetSimMesh.pos_node[tet[0]] + tetSimMesh.pos_node[tet[1]] + tetSimMesh.pos_node[tet[2]] + tetSimMesh.pos_node[tet[3]]);

					outfile9 << std::scientific << std::setprecision(8) << pos[0] << " " << pos[1] << " " << pos[2] << " " << maxEigenValue << " " << tetSimMesh.Dp[nf] << std::endl;
				}				
				outfile9.close();
			}


		}


		std::fill(tetSimMesh.elastForce_node.begin(), tetSimMesh.elastForce_node.end(), Eigen::Vector3d::Zero());


		// Apply gravity
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int nf = 0; nf < tetSimMesh.elastForce_node.size(); nf++)
		{
			Eigen::Vector3d f_node = compute_external_force(tetSimMesh, nf, timestep);
			f_node += tetSimMesh.mass_node[nf] * parameters.gravity;
			tetSimMesh.elastForce_node[nf] += f_node;
		}

		// Compute deformation gradient F
		tetSimMesh.update_F(parameters.numOfThreads);
		updateMLS_after_advection(tetSimMesh, parameters);


		// Compute forces from elements
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int nf = 0; nf < tetSimMesh.tetrahedrals.size(); nf++)
		{
			int matInd = tetSimMesh.materialInd[nf];
			Material mat = tetSimMesh.materialMesh[matInd];


			// Compute First Piola-Kirchhoff stress tensor P
			Eigen::Matrix3d F = tetSimMesh.tetra_F[nf];
			double J = F.determinant();
			Eigen::Matrix3d FinvT = F.inverse().transpose();
			Eigen::Matrix3d PK1 = mat.mu * (F - FinvT) + mat.lambda * log(J) * FinvT;
			Eigen::Matrix3d P = PK1;
			{
				Eigen::Matrix3d cauchyStress = 1.0 / J * PK1 * F.transpose();
				// compute eigenvalue and eigenvector
				Eigen::EigenSolver<Eigen::MatrixXd> es(cauchyStress);
				Eigen::Vector3d eigenValues = { es.eigenvalues()[0].real() ,  es.eigenvalues()[1].real() ,  es.eigenvalues()[2].real() };
				Eigen::Matrix3d eigenVectors;
				eigenVectors << es.eigenvectors().col(0)[0].real(), es.eigenvectors().col(1)[0].real(), es.eigenvectors().col(2)[0].real(),
					es.eigenvectors().col(0)[1].real(), es.eigenvectors().col(1)[1].real(), es.eigenvectors().col(2)[1].real(),
					es.eigenvectors().col(0)[2].real(), es.eigenvectors().col(1)[2].real(), es.eigenvectors().col(2)[2].real();
				double maxEigenValue = std::max(std::max(eigenValues[0], eigenValues[1]), eigenValues[2]);


				if (maxEigenValue > mat.thetaf)
				{
					double tempDp = (1 + mat.Hs) * (1 - mat.thetaf / maxEigenValue);
					if (maxEigenValue > (1 + 1 / mat.Hs) * mat.thetaf)
					{
						tetSimMesh.Dp[nf] = 1.0;
					}
					else
					{
						if (tempDp > tetSimMesh.Dp[nf])
						{
							tetSimMesh.Dp[nf] = tempDp;
						};
					};
				};

				Eigen::Vector3d sigmaPlus = { 0 , 0 , 0 };
				for (int i = 0; i < 3; i++)
				{
					if (eigenValues[i] > 0)
					{
						if (tetSimMesh.Dp[nf] >= 1.0)
						{
							sigmaPlus[i] = 0;
						}
						else
						{
							sigmaPlus[i] = (1 - tetSimMesh.Dp[nf]) * eigenValues[i];
						};
					}
					else
					{
						sigmaPlus[i] = eigenValues[i];
					};
				};

				Eigen::Matrix3d sigma = Eigen::Matrix3d::Zero();
				for (int i = 0; i < 3; i++)
				{
					sigma = sigma + sigmaPlus[i] * eigenVectors.col(i) * (eigenVectors.col(i).transpose());
				};

				cauchyStress = sigma;

				P = J * cauchyStress * F.transpose().inverse();
			}


			// Compute forces on nodes
			Eigen::Matrix3d H = -tetSimMesh.tetra_vol[nf] * P * tetSimMesh.tetra_DM_inv[nf].transpose();
			tetrahedralStress[nf] = H;


			//// Distribute forces to nodes
			//Eigen::Vector3d f0 = H.col(0);
			//Eigen::Vector3d f1 = H.col(1);
			//Eigen::Vector3d f2 = H.col(2);
			//Eigen::Vector3d f3 = -f0 - f1 - f2;


			//Eigen::Vector4i tet = tetSimMesh.tetrahedrals[nf];
			//tetSimMesh.elastForce_node[tet[0]] += f3;
			//tetSimMesh.elastForce_node[tet[1]] += f0;
			//tetSimMesh.elastForce_node[tet[2]] += f1;
			//tetSimMesh.elastForce_node[tet[3]] += f2;

		}

		// Time integration (Explicit Euler)
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int nf = 0; nf < tetSimMesh.elastForce_node.size(); nf++)
		{
			// Calcualte the elastic force contribution
			for (int k = 0; k < nodeSharedByElement[nf].size(); k++)
			{
				Eigen::Vector2i indices = nodeSharedByElement[nf][k];
				int tetIndex = indices[0], nodeIndex = indices[1];
				if (nodeIndex == 0)
				{
					tetSimMesh.elastForce_node[nf] += -(tetrahedralStress[tetIndex].col(0)+ tetrahedralStress[tetIndex].col(1)+ tetrahedralStress[tetIndex].col(2));
				}
				else
				{
					tetSimMesh.elastForce_node[nf] += tetrahedralStress[tetIndex].col(nodeIndex - 1);
				}
			}


			if (tetSimMesh.pos_node[nf][1] >= -9.0)
			{
				// Update velocity and position
				Eigen::Vector3d acceleration = tetSimMesh.elastForce_node[nf] / tetSimMesh.mass_node[nf];
				tetSimMesh.vel_node[nf] += acceleration * parameters.dt;
				tetSimMesh.pos_node[nf] += tetSimMesh.vel_node[nf] * parameters.dt;

			}

		}



	}


}

// implicit integration
void implicitFEM(Mesh& tetSimMesh, FEMParamters& parameters)
{

	for (int timestep = 0; timestep < parameters.num_timesteps; timestep++)
	{
		std::cout << "timestep = " << timestep << std::endl;

		if (timestep % parameters.outputFrequency == 0)
		{
			tetSimMesh.exportSurfaceMesh("surfMesh", timestep);

			{
				std::ofstream outfile2("./output/MLS-Pts"+std::to_string(timestep) + ".obj", std::ios::trunc);
				for (std::map<int, std::vector<MLSPoints>>::iterator it = tetSimMesh.MLSPoints_tet_map.begin();
					it != tetSimMesh.MLSPoints_tet_map.end(); it++) // each tetrahedral that is replaced by MLS points
				{
					for (int j = 0; j < it->second.size(); j++)
					{
						Eigen::Vector3d pos = it->second[j].pos;
						outfile2 << "v " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
					}
				}
				outfile2.close();
			}

		}

		tetSimMesh.pos_node_prev = tetSimMesh.pos_node;

		double lastEnergyVal = compute_IP_energy(tetSimMesh, parameters, timestep);
		for (int ite = 0; ite < 15; ite++)
		{
			std::vector<Eigen::Vector3d> currentPosition = tetSimMesh.pos_node;
			//if (timestep % parameters.outputFrequency == 0)
			//{
			//	tetSimMesh.exportSurfaceMesh("surfMesh_ite_"+std::to_string(ite)+"_", timestep);
			//}

			std::cout << "		ite = " << ite << "; lastEnergyVal = " << lastEnergyVal << std::endl;
			std::vector<Eigen::Vector3d> direction = solve_linear_system(tetSimMesh, parameters, timestep);
			double dist_to_converge = infiniteNorm(direction);

			/*if (timestep == 0)
			{
				for (int m = 0; m < direction.size(); m++)
				{
					if (std::abs(direction[m][2] + 7.84e-5) > 1.0e-6)
					{
						std::cout << "			direction[" << m << "] = " << direction[m][2] << std::endl;
					}
				}
			}*/


			
			std::cout << std::scientific << std::setprecision(4) << "			dist_to_converge = " 
				<< dist_to_converge / parameters.dt << "m/s; threshold = " << parameters.searchResidual 
				<< "m/s" << std::endl;
			if (ite && (dist_to_converge / parameters.dt < parameters.searchResidual))
			{
				break;
			}

			std::cout << "			Calculate step size;" << std::endl;
			double step = calMaxStepSize(tetSimMesh, parameters, timestep, direction);
			//step = 1.0;
			std::cout << "			Step forward = " << step << std::endl;
			step_forward(parameters, tetSimMesh, currentPosition, direction, step);
			double newEnergyVal = compute_IP_energy(tetSimMesh, parameters, timestep);
			std::cout << std::scientific << std::setprecision(4) << "			step = " 
				<< step << "; newEnergyVal = " << newEnergyVal << "; dis = " 
				<< tetSimMesh.pos_node[7][2] - tetSimMesh.pos_node[0][2] << std::endl;
			while (newEnergyVal > lastEnergyVal && step >= 1.0e-5)
			{
				step /= 2.0;

				step_forward(parameters, tetSimMesh, currentPosition, direction, step);
				newEnergyVal = compute_IP_energy(tetSimMesh, parameters, timestep);
				std::cout << "				step = " << step << "; newEnergyVal = " 
					<< newEnergyVal << "; dis = " << tetSimMesh.pos_node[7][2] - tetSimMesh.pos_node[0][2]
					<< std::endl;
				if (std::abs(newEnergyVal - lastEnergyVal) / lastEnergyVal < 0.001) // traped in the local mimima
				{
					break;
				}
			}
			lastEnergyVal = newEnergyVal;


			std::cout << "			lastEnergyVal = " << lastEnergyVal << std::endl << std::endl;



		}


		// update the velocity of the node
		for (int i = 0; i < tetSimMesh.pos_node.size(); i++)
		{
			tetSimMesh.vel_node[i] = (tetSimMesh.pos_node[i] - tetSimMesh.pos_node_prev[i]) / parameters.dt;
		}

	}

}

// compute external force excluding gravity
Eigen::Vector3d compute_external_force(Mesh& tetSimMesh, int vertInd, int timestep)
{
	Eigen::Vector3d extForce = Eigen::Vector3d::Zero();

	if (tetSimMesh.boundaryCondition_node[vertInd].type == 2)
	{
		if (timestep >= tetSimMesh.boundaryCondition_node[vertInd].appliedTime[0] 
			&& timestep <= tetSimMesh.boundaryCondition_node[vertInd].appliedTime[1])
		{
			extForce = tetSimMesh.boundaryCondition_node[vertInd].force[timestep];
		}
	}

	return extForce;
}

void updateMLS_after_advection(Mesh& tetSimMesh, FEMParamters& parameters)
{
	for (std::map<int, std::vector<MLSPoints>>::iterator it = tetSimMesh.MLSPoints_tet_map.begin(); 
		it != tetSimMesh.MLSPoints_tet_map.end(); it++) // each tetrahedral that is replaced by MLS points
	{
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int j = 0; j < it->second.size(); j++)
		{
			it->second[j].MLS_approximation(tetSimMesh.pos_node_Rest, tetSimMesh.pos_node, parameters.MLS_radius);
		}
	}
}

// compute the incremental potential energy
double compute_IP_energy(Mesh& tetSimMesh, FEMParamters& parameters, int timestep)
{
	double energyVal = 0;
	tetSimMesh.update_F(parameters.numOfThreads);
	updateMLS_after_advection(tetSimMesh, parameters);

	// energy contribution per vertex
	std::vector<double> node_ext_ine_energy_vec(tetSimMesh.pos_node.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int vI = 0; vI < tetSimMesh.pos_node.size(); vI++)
	{
		double nodeMass = tetSimMesh.mass_node[vI];
		Eigen::Vector3d xt = tetSimMesh.pos_node_prev[vI];
		Eigen::Vector3d v = tetSimMesh.vel_node[vI];
		Eigen::Vector3d x = tetSimMesh.pos_node[vI];
		Eigen::Vector3d extForce = compute_external_force(tetSimMesh, vI, timestep);

		// the external energy contribution
		double eng = 0;
		eng += ExternalEnergy::Val(nodeMass, parameters.dt, x, parameters, extForce);

		// the inertia energy contribution
		eng += InertiaEnergy::Val(nodeMass, parameters.dt, xt, v, x, extForce, parameters, 
			vI, tetSimMesh.boundaryCondition_node, timestep);
		node_ext_ine_energy_vec[vI] = eng;

	}
	energyVal = std::accumulate(node_ext_ine_energy_vec.begin(), node_ext_ine_energy_vec.end(), 0.0);

	//std::cout << "		Vertex energy = " << energyVal << std::endl;


	int numEleDeg = 0, numMLSDeg = 0;

	// energy contribution per element
	std::vector<double> tex_est_energy_vec(tetSimMesh.tetra_F.size());
//#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int eI = 0; eI < tetSimMesh.tetra_F.size(); eI++)
	{
		// the internal elastic energy contribution
		int matInd = tetSimMesh.materialInd[eI];
		if (tetSimMesh.MLSPoints_tet_map.find(eI) == tetSimMesh.MLSPoints_tet_map.end()) // if this tetrahedral is not replaced by MLS points
		{
			tex_est_energy_vec[eI] = ElasticEnergy::Val(tetSimMesh.materialMesh[matInd], 
				parameters.model, tetSimMesh.tetra_F[eI], parameters.dt, tetSimMesh.tetra_vol[eI]);	

			if (tetSimMesh.tetra_F[eI].determinant() <= 0)
			{
				numEleDeg += 1;
			}
			
		}
		else
		{
			double energyMLS_tmp = 0;
			for (int m = 0; m < tetSimMesh.MLSPoints_tet_map[eI].size(); m++)
			{
				energyMLS_tmp += ElasticEnergy::Val(tetSimMesh.materialMesh[matInd], 
					parameters.model, tetSimMesh.MLSPoints_tet_map[eI][m].F, parameters.dt, 
					tetSimMesh.MLSPoints_tet_map[eI][m].volume);

				if (tetSimMesh.MLSPoints_tet_map[eI][m].F.determinant() <= 0)
				{
					numMLSDeg += 1;
				}
			}
			tex_est_energy_vec[eI] = energyMLS_tmp;
		}
		
	}
	energyVal += std::accumulate(tex_est_energy_vec.begin(), tex_est_energy_vec.end(), 0.0);

	//std::cout << "		Element energy = " << std::accumulate(tex_est_energy_vec.begin(), tex_est_energy_vec.end(), 0.0) << std::endl;
	
	if (numEleDeg + numMLSDeg > 0)
	{
		std::cout << "		numEleDeg = " << numEleDeg << "; numMLSDeg = " << numMLSDeg << std::endl;
		std::cout << std::endl;
	}
	

	// energy contribution from barrier
	energyVal += compute_Barrier_energy(tetSimMesh, parameters, timestep);

	//std::cout << "		Barrier energy = " << compute_Barrier_energy(tetSimMesh, parameters, timestep) << std::endl;

	return energyVal;
}

double compute_Barrier_energy(Mesh& tetSimMesh, FEMParamters& parameters, int timestep)
{
	std::vector<Vector5i> PG_PG, PT_PP, PT_PE, PT_PT, EE_EE;
	calContactInfo(tetSimMesh, parameters, timestep, PG_PG, 
		PT_PP, PT_PE, PT_PT, EE_EE);


	double barrierEnergy = 0;
	std::vector<double> energy_PT_PP(PT_PP.size(),0), energy_PT_PE(PT_PE.size(),0), 
		energy_PT_PT(PT_PT.size(),0), energy_EE_EE(EE_EE.size(),0), energy_PG_PG(PG_PG.size(),0);
	{
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < PG_PG.size(); i++)
		{
			int ptInd = PG_PG[i][1];
			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
			double z2 = P[2] * P[2];
			energy_PG_PG[i] = Ground::val(z2,tetSimMesh.boundaryVertices_area[ptInd], parameters);
		}
		barrierEnergy += std::accumulate(energy_PG_PG.begin(), energy_PG_PG.end(), 0.0);


#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < PT_PP.size(); i++)
		{
			Vector5i cont_PT = PT_PP[i];

			int ptInd = cont_PT[1];
			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
			Eigen::Vector3i tri = tetSimMesh.boundaryTriangles[cont_PT[2]];
			Eigen::Vector3d A = tetSimMesh.pos_node[tri[0]];
			Eigen::Vector3d B = tetSimMesh.pos_node[tri[1]];
			Eigen::Vector3d C = tetSimMesh.pos_node[tri[2]];

			double dis2 = 0;
			DIS::computePointTriD(P, A, B, C, dis2, cont_PT[3]);
			energy_PT_PP[i] = BarrierEnergy::val_PT(tetSimMesh.boundaryVertices_area[ptInd], dis2, parameters);
		}
		barrierEnergy += std::accumulate(energy_PT_PP.begin(), energy_PT_PP.end(), 0.0);


#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < PT_PE.size(); i++)
		{
			Vector5i cont_PT = PT_PE[i];

			int ptInd = cont_PT[1];
			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
			Eigen::Vector3i tri = tetSimMesh.boundaryTriangles[cont_PT[2]];
			Eigen::Vector3d A = tetSimMesh.pos_node[tri[0]];
			Eigen::Vector3d B = tetSimMesh.pos_node[tri[1]];
			Eigen::Vector3d C = tetSimMesh.pos_node[tri[2]];

			double dis2 = 0;
			DIS::computePointTriD(P, A, B, C, dis2, cont_PT[3]);
			energy_PT_PE[i] = BarrierEnergy::val_PT(tetSimMesh.boundaryVertices_area[ptInd], dis2, parameters);

		}
		barrierEnergy += std::accumulate(energy_PT_PE.begin(), energy_PT_PE.end(), 0.0);


#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < PT_PT.size(); i++)
		{
			Vector5i cont_PT = PT_PT[i];

			int ptInd = cont_PT[1];
			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
			Eigen::Vector3i tri = tetSimMesh.boundaryTriangles[cont_PT[2]];
			Eigen::Vector3d A = tetSimMesh.pos_node[tri[0]];
			Eigen::Vector3d B = tetSimMesh.pos_node[tri[1]];
			Eigen::Vector3d C = tetSimMesh.pos_node[tri[2]];

			double dis2 = 0;
			DIS::computePointTriD(P, A, B, C, dis2, cont_PT[3]);
			energy_PT_PT[i] = BarrierEnergy::val_PT(tetSimMesh.boundaryVertices_area[ptInd], dis2, parameters);

		}
		barrierEnergy += std::accumulate(energy_PT_PT.begin(), energy_PT_PT.end(), 0.0);


#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < EE_EE.size(); i++)
		{
			Vector5i cont_EE = EE_EE[i];
			int E1 = cont_EE[1], E2 = cont_EE[2];
			int e1p1 = tetSimMesh.index_boundaryEdge[E1][0], e1p2 = tetSimMesh.index_boundaryEdge[E1][1], 
				e2p1 = tetSimMesh.index_boundaryEdge[E2][0], e2p2 = tetSimMesh.index_boundaryEdge[E2][1];

			Eigen::Vector3d P1 = tetSimMesh.pos_node[e1p1];
			Eigen::Vector3d P2 = tetSimMesh.pos_node[e1p2];
			Eigen::Vector3d Q1 = tetSimMesh.pos_node[e2p1];
			Eigen::Vector3d Q2 = tetSimMesh.pos_node[e2p2];

			double dis2 = 0;
			DIS::computeEdgeEdgeD(P1, P2, Q1, Q2, dis2, cont_EE[3]);
			Eigen::Vector4i ptIndices = { e1p1 , e1p2 , e2p1 , e2p2 };
			energy_EE_EE[i] = BarrierEnergy::val_EE(tetSimMesh.boundaryEdges_area[e1p1][e1p2], 
				dis2, tetSimMesh, ptIndices, parameters);

		}
		barrierEnergy += std::accumulate(energy_EE_EE.begin(), energy_EE_EE.end(), 0.0);

	}


	return barrierEnergy;
}

std::pair<std::vector<Eigen::Vector2i>, int> cal_temporary_MLS_tet_vector(Mesh& tetSimMesh)
{
	std::vector<Eigen::Vector2i> temp_MLS_tet_vector(tetSimMesh.MLSPoints_tet_map.size(), Eigen::Vector2i::Zero());

	int num_MLS_tet = 0;
	int totalNodes_tets = 0;

	for (std::map<int, std::vector<MLSPoints>>::iterator it = tetSimMesh.MLSPoints_tet_map.begin(); 
		it != tetSimMesh.MLSPoints_tet_map.end(); it++) // each tetrahedral that is replaced by MLS points
	{
		temp_MLS_tet_vector[num_MLS_tet][0] = it->first;
		temp_MLS_tet_vector[num_MLS_tet][1] = totalNodes_tets;
		for (int j = 0; j < it->second.size(); j++)
		{
			totalNodes_tets += it->second[j].index_node.size();
		}
		num_MLS_tet += 1;
	}

	return std::make_pair(temp_MLS_tet_vector, totalNodes_tets);
}

// compute energy gradient and hessian of the linear system and solve the syetem
std::vector<Eigen::Vector3d> solve_linear_system(Mesh& tetSimMesh, FEMParamters& parameters, int timestep)
{
	std::vector<Eigen::Vector3d> movingDir(tetSimMesh.pos_node.size());

	// update defromation gradient
	tetSimMesh.update_F(parameters.numOfThreads);
	updateMLS_after_advection(tetSimMesh, parameters);


	//double startTime1, endTime1;
	//startTime1 = omp_get_wtime();

	std::vector<Vector5i> PG_PG, PT_PP, PT_PE, PT_PT, EE_EE;
	calContactInfo(tetSimMesh, parameters, timestep, PG_PG, PT_PP, PT_PE, PT_PT, EE_EE);

	std::cout << "			PT_PP.size() = " << PT_PP.size();
	std::cout << "; PT_PE.size() = " << PT_PE.size();
	std::cout << "; PT_PT.size() = " << PT_PT.size();
	std::cout << "; EE_EE.size() = " << EE_EE.size() << std::endl;


	//endTime1 = omp_get_wtime();
	//std::cout << "	Cal contact Time is : " << endTime1 - startTime1 << "s" << std::endl;


	// !!!!!!!!!!!!!!!!!!!!!!
	// Some tetrahedrals (a small part ) are replaced by MLS points, so the energy and gradient of these tetrahedrals are not calculated
	// But in order to facilitate parallelization, the energy and gradient of these tetrahedrals are set as 0, but the results are not used
	int gradSize_beforeMLS = tetSimMesh.pos_node.size() * 6 + tetSimMesh.tetrahedrals.size() * 12 + PG_PG.size() * 3
		+ PT_PP.size() * 6 + PT_PE.size() * 9 + PT_PT.size() * 12 + EE_EE.size() * 12;
	int hessSize_beforeMLS = tetSimMesh.pos_node.size() * 3 + tetSimMesh.tetrahedrals.size() * 144 + PG_PG.size() * 9
		+ PT_PP.size() * 36 + PT_PE.size() * 81 + PT_PT.size() * 144 + EE_EE.size() * 144;
	// adding those nodes by MLS points
	std::pair<std::vector<Eigen::Vector2i>, int> temp_MLS_tet_vector_and_numNodes = cal_temporary_MLS_tet_vector(tetSimMesh);
	int gradSize = gradSize_beforeMLS + temp_MLS_tet_vector_and_numNodes.second * 3;
	int hessSize = hessSize_beforeMLS + temp_MLS_tet_vector_and_numNodes.second * 9;
	std::vector<std::pair<int, double>> grad_triplet(gradSize);
	std::vector<Eigen::Triplet<double>> hessian_triplet(hessSize);


	// energy contribution per vertex
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int vI = 0; vI < tetSimMesh.pos_node.size(); vI++)
	{
		double nodeMass = tetSimMesh.mass_node[vI];
		Eigen::Vector3d xt = tetSimMesh.pos_node_prev[vI];
		Eigen::Vector3d v = tetSimMesh.vel_node[vI];
		Eigen::Vector3d x = tetSimMesh.pos_node[vI];
		Eigen::Vector3d extForce = compute_external_force(tetSimMesh, vI, timestep);

		// !!!! The order is very important. "InertiaEnergy" must be the first and "ExternalEnergy" the second

		// the inertia energy contribution
		int actGradIndex = vI * 6;
		InertiaEnergy::Grad(grad_triplet, actGradIndex, nodeMass, parameters.dt, 
			xt, v, x, extForce, vI, parameters, tetSimMesh.boundaryCondition_node, timestep);
		int actHessIndex = vI * 3;
		InertiaEnergy::Hess(hessian_triplet, actHessIndex, nodeMass, vI, 
			tetSimMesh.boundaryCondition_node, timestep, parameters);

		// the external energy contribution
		actGradIndex = vI * 6 + 3;
		ExternalEnergy::Grad(grad_triplet, actGradIndex, nodeMass, parameters.dt, 
			x, parameters, extForce, vI, tetSimMesh.boundaryCondition_node[vI].type);

	}


	int startIndex_grad = tetSimMesh.pos_node.size() * 6, startIndex_hess = tetSimMesh.pos_node.size() * 3;
	// energy contribution per element
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int eI = 0; eI < tetSimMesh.tetrahedrals.size(); eI++)
	{

		if (tetSimMesh.MLSPoints_tet_map.find(eI) == tetSimMesh.MLSPoints_tet_map.end()) // if this tetrahedral is not replaced by MLS points
		{
			Eigen::Vector4i tetVertInd_BC;
			for (int hg = 0; hg < 4; hg++)
			{
				int vertInd = tetSimMesh.tetrahedrals[eI][hg];
				tetVertInd_BC[hg] = tetSimMesh.boundaryCondition_node[vertInd].type;
			}

			// the internal elastic energy contribution
			int matInd = tetSimMesh.materialInd[eI];
			Eigen::Matrix<double, 9, 12> dFdx = ElasticEnergy::dF_wrt_dx(tetSimMesh.tetra_DM_inv[eI]);
			int actualStartIndex_hess = startIndex_hess + eI * 144, 
				actualStartIndex_grad = startIndex_grad + eI * 12;
			ElasticEnergy::Grad(grad_triplet, actualStartIndex_grad, tetSimMesh.materialMesh[matInd], 
				parameters.model, tetSimMesh.tetra_F[eI], parameters.dt, tetSimMesh.tetra_vol[eI], dFdx, 
				tetSimMesh.tetrahedrals[eI], tetVertInd_BC);
			ElasticEnergy::Hess(hessian_triplet, actualStartIndex_hess, tetSimMesh.materialMesh[matInd], 
				parameters.model, tetSimMesh.tetra_F[eI], parameters.dt, tetSimMesh.tetra_vol[eI], dFdx,
				tetSimMesh.tetrahedrals[eI], tetVertInd_BC);

		}
		else // placeholder in order to facilitate parallelization
		{
			Eigen::Triplet<double> pa = { 0, 0, 0.0 };
			for (int empi = 0; empi < 12; empi++)
			{
				grad_triplet[startIndex_grad + empi] = {0, 0.0 };
				for (int empj = 0; empj < 12; empj++)
				{
					hessian_triplet[startIndex_hess + empi * 12 + empj] = pa;
				}
			}
		}


	}


	//double endTime2 = omp_get_wtime();
	//std::cout << "	Element grad_hess Time is : " << endTime2 - endTime1 << "s" << std::endl;



	// calculate the barrier gradient and hessian
	{
		startIndex_grad += tetSimMesh.tetrahedrals.size() * 12, 
			startIndex_hess += tetSimMesh.tetrahedrals.size() * 144;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < PG_PG.size(); i++)
		{
			int ptInd = PG_PG[i][1];
			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
			int actualStartIndex_hess = startIndex_hess + i * 9, 
				actualStartIndex_grad = startIndex_grad + i * 3;
			double z2 = P[2] * P[2];
			Ground::gradAndHess(hessian_triplet, grad_triplet, actualStartIndex_hess, 
				actualStartIndex_grad, ptInd, 
				z2, tetSimMesh.boundaryVertices_area[ptInd],parameters);
		}

		startIndex_grad += PG_PG.size() * 3, startIndex_hess += PG_PG.size() * 9;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < PT_PP.size(); i++)
		{
			Vector5i cont_PT = PT_PP[i];

			int ptInd = cont_PT[1];
			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
			Eigen::Vector3i tri = tetSimMesh.boundaryTriangles[cont_PT[2]];
			Eigen::Vector3d A = tetSimMesh.pos_node[tri[0]];
			Eigen::Vector3d B = tetSimMesh.pos_node[tri[1]];
			Eigen::Vector3d C = tetSimMesh.pos_node[tri[2]];
			int type = cont_PT[3];

			double dis2 = 0;
			DIS::computePointTriD(P, A, B, C, dis2);
			Eigen::Vector4i ptIndices = { ptInd , tri[0] , tri[1] , tri[2] };
			int actualStartIndex_hess = startIndex_hess + i * 36, 
				actualStartIndex_grad = startIndex_grad + i * 6;
			BarrierEnergy::gradAndHess_PT(hessian_triplet, grad_triplet, actualStartIndex_hess, 
				actualStartIndex_grad, ptIndices, type, 
				dis2, tetSimMesh, parameters);

		}

		startIndex_grad += PT_PP.size() * 6, startIndex_hess += PT_PP.size() * 36;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < PT_PE.size(); i++)
		{
			Vector5i cont_PT = PT_PE[i];

			int ptInd = cont_PT[1];
			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
			Eigen::Vector3i tri = tetSimMesh.boundaryTriangles[cont_PT[2]];
			Eigen::Vector3d A = tetSimMesh.pos_node[tri[0]];
			Eigen::Vector3d B = tetSimMesh.pos_node[tri[1]];
			Eigen::Vector3d C = tetSimMesh.pos_node[tri[2]];
			int type = cont_PT[3];

			double dis2 = 0;
			DIS::computePointTriD(P, A, B, C, dis2);
			Eigen::Vector4i ptIndices = { ptInd , tri[0] , tri[1] , tri[2] };
			int actualStartIndex_hess = startIndex_hess + i * 81,
				actualStartIndex_grad = startIndex_grad + i * 9;
			BarrierEnergy::gradAndHess_PT(hessian_triplet, grad_triplet, actualStartIndex_hess, 
				actualStartIndex_grad, ptIndices, type, 
				dis2, tetSimMesh, parameters);

		}

		startIndex_grad += PT_PE.size() * 9, startIndex_hess += PT_PE.size() * 81;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < PT_PT.size(); i++)
		{
			Vector5i cont_PT = PT_PT[i];

			int ptInd = cont_PT[1];
			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
			Eigen::Vector3i tri = tetSimMesh.boundaryTriangles[cont_PT[2]];
			Eigen::Vector3d A = tetSimMesh.pos_node[tri[0]];
			Eigen::Vector3d B = tetSimMesh.pos_node[tri[1]];
			Eigen::Vector3d C = tetSimMesh.pos_node[tri[2]];
			int type = cont_PT[3];

			double dis2 = 0;
			DIS::computePointTriD(P, A, B, C, dis2);
			Eigen::Vector4i ptIndices = { ptInd , tri[0] , tri[1] , tri[2] };
			int actualStartIndex_hess = startIndex_hess + i * 144, 
				actualStartIndex_grad = startIndex_grad + i * 12;
			BarrierEnergy::gradAndHess_PT(hessian_triplet, grad_triplet, actualStartIndex_hess, 
				actualStartIndex_grad, ptIndices, type, 
				dis2, tetSimMesh, parameters);

		}

		startIndex_grad += PT_PT.size() * 12, startIndex_hess += PT_PT.size() * 144;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < EE_EE.size(); i++)
		{
			Vector5i cont_EE = EE_EE[i];
			int E1 = cont_EE[1], E2 = cont_EE[2];
			int e1p1 = tetSimMesh.index_boundaryEdge[E1][0], e1p2 = tetSimMesh.index_boundaryEdge[E1][1], 
				e2p1 = tetSimMesh.index_boundaryEdge[E2][0], e2p2 = tetSimMesh.index_boundaryEdge[E2][1];

			Eigen::Vector3d P1 = tetSimMesh.pos_node[e1p1];
			Eigen::Vector3d P2 = tetSimMesh.pos_node[e1p2];
			Eigen::Vector3d Q1 = tetSimMesh.pos_node[e2p1];
			Eigen::Vector3d Q2 = tetSimMesh.pos_node[e2p2];
			int type = cont_EE[3];

			double dis2 = 0;
			DIS::computeEdgeEdgeD(P1, P2, Q1, Q2, dis2);
			Eigen::Vector4i ptIndices = { e1p1 , e1p2 , e2p1 , e2p2 };
			int actualStartIndex_hess = startIndex_hess + i * 144, 
				actualStartIndex_grad = startIndex_grad + i * 12;
			BarrierEnergy::gradAndHess_EE(hessian_triplet, grad_triplet, actualStartIndex_hess, 
				actualStartIndex_grad, ptIndices, type, 
				dis2, tetSimMesh, parameters);

		}

	}


	// now add energy gradient and hessian contributions from MLS points
	std::vector<Eigen::Vector2i> temp_MLS_tet_vector = temp_MLS_tet_vector_and_numNodes.first;
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int it = 0; it < temp_MLS_tet_vector.size(); it++) // each tetrahedral that is replaced by MLS points
	{
		int tetInd = temp_MLS_tet_vector[it][0];
		int matInd = tetSimMesh.materialInd[tetInd];
		int relativeIndex = temp_MLS_tet_vector[it][1];

		int totalNodes_per_tet = 0;
		for (int j = 0; j < tetSimMesh.MLSPoints_tet_map[tetInd].size(); j++)
		{
			MLSPoints mp = tetSimMesh.MLSPoints_tet_map[tetInd][j];
			
			Eigen::Matrix3d PK1 = ElasticEnergy::calPK1(tetSimMesh.materialMesh[matInd], parameters.model, mp.F);
			Vector9d grad_tmp = flatenMatrix3d(PK1);
			Eigen::Matrix<double, 9, 9> hess_tmp = ElasticEnergy::calPK1_wrt_F(tetSimMesh.materialMesh[matInd], 
				parameters.model, mp.F);

			// calculate gradient & hessian
			for (int k = 0; k < mp.index_node.size(); k++)
			{
				int nodeInd = mp.index_node[k];
				Eigen::Vector3d grad_MLS_node = parameters.dt * parameters.dt * mp.volume * mp.dFdx[k].transpose() * grad_tmp;
				Eigen::Matrix3d hess_MLS_node = parameters.dt * parameters.dt * mp.volume * mp.dFdx[k].transpose() * hess_tmp * mp.dFdx[k];

				for (int xd = 0; xd < 3; xd++)
				{
					int currentNode_grad_Index = gradSize_beforeMLS + relativeIndex * 3 + totalNodes_per_tet * 3;
					grad_triplet[currentNode_grad_Index + xd] = { nodeInd * 3 + xd, grad_MLS_node[xd]};

					int currentNode_hess_Index = hessSize_beforeMLS + relativeIndex * 9 
						+ totalNodes_per_tet * 9 + xd * 3;
					for (int yd = 0; yd < 3; yd++)
					{						
						Eigen::Triplet<double> pa = { nodeInd * 3 + xd, nodeInd * 3 + yd, hess_MLS_node(xd, yd)};
						hessian_triplet[currentNode_hess_Index + yd] = pa;
					}
				}

				totalNodes_per_tet += 1;
			}

		}

	}



	//double endTime3 = omp_get_wtime();
	//std::cout << "	Cal grad_hess Time is : " << endTime3 - endTime2 << "s" << std::endl;

	// assemable the left-hand side 
	Eigen::SparseMatrix<double> leftHandSide(3 * tetSimMesh.pos_node.size(), 3 * tetSimMesh.pos_node.size());
	leftHandSide.setFromTriplets(hessian_triplet.begin(), hessian_triplet.end());

	// assemable the right-hand side 
	Eigen::VectorXd rightHandSide = Eigen::VectorXd(tetSimMesh.pos_node.size() * 3);
	rightHandSide.setZero();
	for (int i = 0; i < grad_triplet.size(); i++)
	{
		std::pair<int, double> ele = grad_triplet[i];
		rightHandSide[ele.first] += ele.second;
	}

	//double endTime4 = omp_get_wtime();
	//std::cout << "	Assemble grad_hess Time is : " << endTime4 - endTime3 << "s" << std::endl;

	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
	solver.compute(leftHandSide);
	Eigen::VectorXd result = solver.solve(-rightHandSide);

	

	//double endTime5 = omp_get_wtime();
	//std::cout << "	Solve grad_hess Time is : " << endTime5 - endTime4 << "s" << std::endl;


#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int i = 0; i < tetSimMesh.pos_node.size(); i++)
	{
		movingDir[i] = result.block<3, 1>(3 * i, 0);
	}



	return movingDir;
}

// move points' position according to the direction; Note this is a trial movement
void step_forward(FEMParamters& parameters, Mesh& tetSimMesh, std::vector<Eigen::Vector3d>& currentPosition, 
	std::vector<Eigen::Vector3d>& direction, double step)
{
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int vI = 0; vI < tetSimMesh.pos_node.size(); vI++)
	{
		tetSimMesh.pos_node[vI] = currentPosition[vI] + step * direction[vI];
	}

	updateMLS_after_advection(tetSimMesh, parameters);
}

void calContactInfo(Mesh& tetSimMesh, FEMParamters& parameters, int timestep, 
	std::vector<Vector5i>& PG_PG, std::vector<Vector5i>& PT_PP, std::vector<Vector5i>& PT_PE,
	std::vector<Vector5i>& PT_PT, std::vector<Vector5i>& EE_EE)
{
	std::vector<spatialHashCellData> spatialHash_vec;
	std::map<std::string, int> hashNameIndex;
	std::vector<Eigen::Vector3d> direction(tetSimMesh.pos_node.size(), Eigen::Vector3d::Zero());
	initSpatialHash(false, parameters, tetSimMesh, direction, parameters.IPC_hashSize, 
		spatialHash_vec, hashNameIndex, timestep);



	// Step 1: find contact pairs and types
	std::vector<std::vector<Vector5i>> cont_pair_vec(spatialHash_vec.size()); // 1st int: 0(PT), 1(EE), 2(PG); 2nd int: index of P(E1)(P: for ground contact case); 3rd int: index of T(E2); 4th int: type; 5th int: actual involved elements, i.e. PP(2), PE(3), PT(4) and EE(4)  
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int y = 0; y < spatialHash_vec.size(); y++)
	{
		std::vector<Vector5i> cont_pair;
		for (std::set<int>::iterator itP = spatialHash_vec[y].vertIndices.begin();
			itP != spatialHash_vec[y].vertIndices.end(); itP++)
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
							for (std::set<int>::iterator itT = spatialHash_vec[neigHashIndex].triaIndices.begin();
								itT != spatialHash_vec[neigHashIndex].triaIndices.end(); itT++)
							{
								int vert = *itP, tri = *itT;
								if (tetSimMesh.note_node[vert] 
									!= tetSimMesh.note_node[tetSimMesh.boundaryTriangles[tri][0]])
								{
									if (tetSimMesh.boundaryVertices[vert].find(tri) 
										== tetSimMesh.boundaryVertices[vert].end()) // this triangle is not incident with the point
									{
										Eigen::Vector3i triVerts = tetSimMesh.boundaryTriangles[tri];
										Eigen::Vector3d P = tetSimMesh.pos_node[vert];
										Eigen::Vector3d A = tetSimMesh.pos_node[triVerts[0]];
										Eigen::Vector3d B = tetSimMesh.pos_node[triVerts[1]];
										Eigen::Vector3d C = tetSimMesh.pos_node[triVerts[2]];

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
		}

		for (std::set<int>::iterator itE1 = spatialHash_vec[y].edgeIndices.begin(); 
			itE1 != spatialHash_vec[y].edgeIndices.end(); itE1++)
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
							for (std::set<int>::iterator itE2 = spatialHash_vec[neigHashIndex].edgeIndices.begin(); 
								itE2 != spatialHash_vec[neigHashIndex].edgeIndices.end(); itE2++)
							{
								if (*itE1 != *itE2)
								{
									Eigen::Vector2i E1 = tetSimMesh.index_boundaryEdge[*itE1], 
										E2 = tetSimMesh.index_boundaryEdge[*itE2];
									int P1I = E1[0], P2I = E1[1], Q1I = E2[0], Q2I = E2[1];
									if (tetSimMesh.note_node[P1I] != tetSimMesh.note_node[Q2I])
									{
										if (P1I != Q1I && P1I != Q2I && P2I != Q1I && P2I != Q2I) // not duplicated and incident edges
										{
											Eigen::Vector3d P1 = tetSimMesh.pos_node[P1I];
											Eigen::Vector3d P2 = tetSimMesh.pos_node[P2I];
											Eigen::Vector3d Q1 = tetSimMesh.pos_node[Q1I];
											Eigen::Vector3d Q2 = tetSimMesh.pos_node[Q2I];

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
				std::string EEPair = std::to_string(std::min(cont_pair_vec[i][j][1], 
					cont_pair_vec[i][j][2])) + "#" + std::to_string(std::max(cont_pair_vec[i][j][1], 
						cont_pair_vec[i][j][2]));
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
		for (int ft = 0; ft < tetSimMesh.boundaryVertices_vec.size(); ft++)
		{
			int ptInd = tetSimMesh.boundaryVertices_vec[ft];
			if (tetSimMesh.pos_node[ptInd][2] <= parameters.IPC_dis)
			{
				Vector5i ct = { 2 , ptInd , 0 , 0 , 0 };
				PG_PG.push_back(ct);
			}
		}

	}

}

// calculate the maximum feasible step size
double calMaxStepSize(Mesh& tetSimMesh, FEMParamters& parameters, int timestep, 
	std::vector<Eigen::Vector3d>& direction)
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
	//	bool intersect = pointTriangleCCDBroadphase(tetSimMesh.pos_node[p0], direction[p0], tetSimMesh.pos_node[p1], 
	//		direction[p1], tetSimMesh.pos_node[p2], direction[p2], tetSimMesh.pos_node[p3], direction[p3], parameters.IPC_dis);
	//	//std::cout << "	intersect = " << intersect << std::endl;
	//	if (intersect)
	//	{
	//		//std::cout << "		tetSimMesh.pos_node[p0] = "<< tetSimMesh.pos_node[p0] << std::endl;
	//		//std::cout << "		tetSimMesh.pos_node[p1] = "<< tetSimMesh.pos_node[p1] << std::endl;
	//		//std::cout << "		tetSimMesh.pos_node[p2] = "<< tetSimMesh.pos_node[p2] << std::endl;
	//		//std::cout << "		tetSimMesh.pos_node[p3] = "<< tetSimMesh.pos_node[p3] << std::endl;
	//		//std::cout << "		direction[p0] = "<< direction[p0] << std::endl;
	//		//std::cout << "		direction[p1] = "<< direction[p1] << std::endl;
	//		//std::cout << "		direction[p2] = "<< direction[p2] << std::endl;
	//		//std::cout << "		direction[p3] = "<< direction[p3] << std::endl;
	//		double alpha_tmp = pointTriangleCCDNarrowphase(tetSimMesh.pos_node[p0], direction[p0], tetSimMesh.pos_node[p1],
	//			direction[p1], tetSimMesh.pos_node[p2], direction[p2], tetSimMesh.pos_node[p3], direction[p3], parameters.IPC_eta);
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
	//	bool intersect = edgeEdgeCCDBroadphase(tetSimMesh.pos_node[p0], direction[p0], tetSimMesh.pos_node[p1],
	//		direction[p1], tetSimMesh.pos_node[p2], direction[p2], tetSimMesh.pos_node[p3], direction[p3], parameters.IPC_dis);
	//	if (intersect)
	//	{
	//		double alpha_tmp = edgeEdgeCCDNarrowphase(tetSimMesh.pos_node[p0], direction[p0], tetSimMesh.pos_node[p1],
	//			direction[p1], tetSimMesh.pos_node[p2], direction[p2], tetSimMesh.pos_node[p3], direction[p3], parameters.IPC_eta);
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
	//		for (std::map<int, std::set<int>>::iterator it = tetSimMesh.boundaryVertices.begin(); it != tetSimMesh.boundaryVertices.end(); it++)
	//		{
	//			int ptInd = it->first;
	//			if (direction[ptInd][2] < 0)
	//			{
	//				double coor_z = tetSimMesh.pos_node[ptInd][2];
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
		double CCD_step = calMaxStep_spatialHash(parameters, tetSimMesh, direction,
			parameters.IPC_hashSize, parameters.IPC_dis, parameters.IPC_eta);
		if (parameters.enableGround == true)
		{
			double stepGround = 1.0;
			for (std::map<int, std::set<int>>::iterator it = tetSimMesh.boundaryVertices.begin(); 
				it != tetSimMesh.boundaryVertices.end(); it++)
			{
				int ptInd = it->first;
				if (direction[ptInd][2] < 0)
				{
					double coor_z = tetSimMesh.pos_node[ptInd][2];
					stepGround = std::min(stepGround, coor_z * (1.0 - parameters.IPC_eta) / std::abs(direction[ptInd][2]));
				}

			}
			CCD_step = std::min(stepGround, CCD_step);
		}
		return CCD_step;
	}

}


