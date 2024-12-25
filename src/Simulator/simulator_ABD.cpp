

#include "simulator.h"


// implicit integration
void implicitFEM_ABD(Mesh_ABD& tetSimMesh, FEMParamters& parameters)
{

	for (int timestep = 0; timestep < parameters.num_timesteps; timestep++)
	{
		std::cout << "timestep = " << timestep << std::endl;

		if (timestep % parameters.outputFrequency == 0)
		{
			tetSimMesh.exportSurfaceMesh("surfMesh", timestep);
		}

		tetSimMesh.pos_node_prev = tetSimMesh.pos_node;
		tetSimMesh.translation_prev_ABD = tetSimMesh.translation_ABD;
		tetSimMesh.deformation_prev_ABD = tetSimMesh.deformation_ABD;

		std::vector<Eigen::Vector3d> position_direction(tetSimMesh.pos_node.size(), Eigen::Vector3d::Zero());

		double lastEnergyVal = compute_IP_energy_ABD(tetSimMesh, parameters, timestep);
		for (int ite = 0; ite < 15; ite++)
		{
			std::vector<Eigen::Vector3d> currentPosition = tetSimMesh.pos_node;
			std::vector<Eigen::Vector3d> current_ABD_translation = tetSimMesh.translation_ABD;
			std::vector<Eigen::Matrix3d> current_ABD_deformation = tetSimMesh.deformation_ABD;


			std::cout << "		ite = " << ite << "; lastEnergyVal = " << lastEnergyVal << std::endl;
			std::vector<Vector12d> direction_ABD = solve_linear_system_ABD(tetSimMesh, parameters, timestep);
			convert_to_position_direction(parameters, tetSimMesh, direction_ABD, position_direction);
			double dist_to_converge = infiniteNorm(position_direction);

			//if (timestep == 182)
			//{
			//	std::cout << "dir = " << position_direction[0] << std::endl;
			//	std::cout << "dir = " << position_direction[0] << std::endl;
			//}


			//std::cout << std::scientific << std::setprecision(4) << "tetSimMesh.massMatrix_ABD[0] = " << std::endl << tetSimMesh.massMatrix_ABD[0] << std::endl;

			std::cout << std::scientific << std::setprecision(4) << "			dist_to_converge = "
				<< dist_to_converge / parameters.dt << "m/s; threshold = " << parameters.searchResidual
				<< "m/s" << std::endl;
			if (ite && (dist_to_converge / parameters.dt < parameters.searchResidual))
			{
				break;
			}

			std::cout << "			Calculate step size;" << std::endl;
			double step = calMaxStepSize(
				parameters,
				tetSimMesh.surfaceInfo,
				position_direction,
				tetSimMesh.pos_node,
				tetSimMesh.note_node,
				tetSimMesh.tetMeshIndex,
				timestep);
			//step = 1.0;
			std::cout << "			Step forward = " << step << std::endl;
			step_forward_ABD(parameters, tetSimMesh, current_ABD_translation, current_ABD_deformation, direction_ABD, currentPosition, step);
			double newEnergyVal = compute_IP_energy_ABD(tetSimMesh, parameters, timestep);
			std::cout << std::scientific << std::setprecision(4) << "			step = "
				<< step << "; newEnergyVal = " << newEnergyVal  << std::endl;
			while (newEnergyVal > lastEnergyVal && step >= 1.0e-5)
			{
				step /= 2.0;

				step_forward_ABD(parameters, tetSimMesh, current_ABD_translation, current_ABD_deformation, direction_ABD, currentPosition, step);
				newEnergyVal = compute_IP_energy_ABD(tetSimMesh, parameters, timestep);
				std::cout << "				step = " << step << "; newEnergyVal = "
					<< newEnergyVal << "; dis = " << tetSimMesh.pos_node[7][2] - tetSimMesh.pos_node[0][2]
					<< std::endl;
				if (std::abs(newEnergyVal - lastEnergyVal) / lastEnergyVal < 0.001) // traped in the local mimima
				{
					break;
				}
			}

			// The energy has been greatly minimized. It is time to stop
			if (newEnergyVal / lastEnergyVal < 0.001)
			{
				break;
			}

			lastEnergyVal = newEnergyVal;


			std::cout << "			lastEnergyVal = " << lastEnergyVal << std::endl << std::endl;



		}


		// update the velocity of the node
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < tetSimMesh.pos_node.size(); i++)
		{
			tetSimMesh.vel_node[i] = (tetSimMesh.pos_node[i] - tetSimMesh.pos_node_prev[i]) / parameters.dt;

			if (i == 0)
			{
				std::cout << "			tetSimMesh.vel_node[0] = " << tetSimMesh.vel_node[0].transpose() << std::endl;
			}

		}


		// update the velocity of the ABD objects
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < tetSimMesh.num_meshes; i++)
		{
			tetSimMesh.translation_vel_ABD[i] = (tetSimMesh.translation_ABD[i] - tetSimMesh.translation_prev_ABD[i]) / parameters.dt;
			tetSimMesh.deformation_vel_ABD[i] = (tetSimMesh.deformation_ABD[i] - tetSimMesh.deformation_prev_ABD[i]) / parameters.dt;
		}

	}

}

Vector12d constructQVec(Eigen::Vector3d translation, Eigen::Matrix3d deformation) // construct q vector according translation and deformation
{
	Vector12d qb = Vector12d::Zero();
	qb.block(0, 0, 3, 1) = translation;
	qb.block(3, 0, 3, 1) = deformation.col(0);
	qb.block(6, 0, 3, 1) = deformation.col(1);
	qb.block(9, 0, 3, 1) = deformation.col(2);

	return qb;
}

// compute the incremental potential energy
double compute_IP_energy_ABD(Mesh_ABD& tetSimMesh, FEMParamters& parameters, int timestep)
{
	double energyVal = 0;	

	// initeria energy contribution per affine body
	std::vector<Vector12d> affine_force(tetSimMesh.num_meshes, Vector12d::Zero());
	for (int i = 0; i < tetSimMesh.pos_node.size(); i++)
	{
		if (tetSimMesh.boundaryCondition_node[i].type == 2)
		{
			Eigen::Matrix<double, 3, 12> Jx = build_Jx_matrix_for_ABD(tetSimMesh.pos_node_Rest[i].transpose());

			int AB_index = tetSimMesh.index_node[i][0];
			Eigen::Vector3d ext_force = compute_external_force(tetSimMesh.boundaryCondition_node, i, timestep);
			Vector12d affine_force_node = Jx.transpose() * ext_force;

			affine_force[AB_index] += affine_force_node;		
		}
	}
	std::vector<double> AB_ine_energy_vec(tetSimMesh.num_meshes,0);
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int AI = 0; AI < AB_ine_energy_vec.size(); AI++)
	{
		// reconstruct q_prev vector
		Vector12d qb_prev = constructQVec(tetSimMesh.translation_prev_ABD[AI], tetSimMesh.deformation_prev_ABD[AI]);
		// reconstruct qb_vel vector
		Vector12d qb_vel = constructQVec(tetSimMesh.translation_vel_ABD[AI], tetSimMesh.deformation_vel_ABD[AI]);	
		// reconstruct q vector
		Vector12d qb = constructQVec(tetSimMesh.translation_ABD[AI], tetSimMesh.deformation_ABD[AI]);

		Vector12d qb_Minus_qbHat = qb - (qb_prev + parameters.dt * qb_vel + parameters.dt * parameters.dt * tetSimMesh.massMatrix_ABD[AI].inverse() * affine_force[AI]);
		AB_ine_energy_vec[AI] = 0.5 * qb_Minus_qbHat.transpose() * tetSimMesh.massMatrix_ABD[AI] * qb_Minus_qbHat;
	}	
	energyVal += std::accumulate(AB_ine_energy_vec.begin(), AB_ine_energy_vec.end(), 0.0);



	// affine energy contribution per affine body
	std::vector<double> AB_aff_energy_vec(tetSimMesh.num_meshes,0);
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int AI = 0; AI < AB_aff_energy_vec.size(); AI++)
	{
		// original implementation
		AB_aff_energy_vec[AI] = parameters.dt * parameters.dt * parameters.ABD_Coeff * tetSimMesh.volume_ABD[AI] * (tetSimMesh.deformation_ABD[AI].transpose() * tetSimMesh.deformation_ABD[AI] - Eigen::Matrix3d::Identity()).squaredNorm();

		//// Linear ARAP implementation
		//double J = tetSimMesh.deformation_ABD[AI].determinant();
		//AB_aff_energy_vec[AI] = parameters.dt * parameters.dt * tetSimMesh.volume_ABD[AI] * parameters.ABD_Coeff / 2.0 * (0.5 * (tetSimMesh.deformation_ABD[AI]
		//						+ tetSimMesh.deformation_ABD[AI].transpose()) - Eigen::Matrix3d::Identity()).squaredNorm();

	}
	energyVal += std::accumulate(AB_aff_energy_vec.begin(), AB_aff_energy_vec.end(), 0.0);



	// energy contribution from barrier
	energyVal += compute_Barrier_energy(
		parameters,
		tetSimMesh.surfaceInfo,
		tetSimMesh.pos_node,
		tetSimMesh.pos_node_Rest,
		tetSimMesh.note_node,
		tetSimMesh.tetMeshIndex,
		timestep);


	return energyVal;
}

// compute energy gradient and hessian of the linear system and solve the syetem
std::vector<Vector12d> solve_linear_system_ABD(Mesh_ABD& tetSimMesh, FEMParamters& parameters, int timestep)
{
	// we only need to update boundary vertices' position to calculate the barrier energy update
	std::vector<Vector12d> movingDir(tetSimMesh.num_meshes);



	// contribution from contact
	//		calculate contacts
	std::vector<Vector5i> PG_PG, PT_PP, PT_PE, PT_PT, EE_EE;
	calContactInfo(parameters,
		tetSimMesh.surfaceInfo,
		tetSimMesh.pos_node,
		tetSimMesh.note_node,
		tetSimMesh.tetMeshIndex,
		timestep, PG_PG, PT_PP,
		PT_PE, PT_PT, EE_EE);


	std::cout << "			PT_PP.size() = " << PT_PP.size();
	std::cout << "; PT_PE.size() = " << PT_PE.size();
	std::cout << "; PT_PT.size() = " << PT_PT.size();
	std::cout << "; EE_EE.size() = " << EE_EE.size() << std::endl;



	int gradSize = tetSimMesh.num_meshes * 12 + tetSimMesh.num_meshes * 9 + PG_PG.size() * 12
		+ PT_PP.size() * 24 + PT_PE.size() * 36 + PT_PT.size() * 48 + EE_EE.size() * 48;
	int hessSize = tetSimMesh.num_meshes * 144 + tetSimMesh.num_meshes * 81 + PG_PG.size() * 144
		+ PT_PP.size() * 144 * 4 + PT_PE.size() * 144 * 9 + PT_PT.size() * 144 * 16 + EE_EE.size() * 144 * 16;
	std::vector<std::pair<int, double>> grad_triplet(gradSize, std::make_pair(0,0.0));
	std::vector<Eigen::Triplet<double>> hessian_triplet(hessSize, Eigen::Triplet<double>(0, 0, 0.0));



	// contribution from inertia term 
	std::vector<Vector12d> affine_force(tetSimMesh.num_meshes, Vector12d::Zero());
	for (int i = 0; i < tetSimMesh.pos_node.size(); i++)
	{
		if (tetSimMesh.boundaryCondition_node[i].type == 2)
		{
			Eigen::Matrix<double, 3, 12> Jx = build_Jx_matrix_for_ABD(tetSimMesh.pos_node_Rest[i].transpose());

			int AB_index = tetSimMesh.index_node[i][0];
			Eigen::Vector3d ext_force = compute_external_force(tetSimMesh.boundaryCondition_node, i, timestep);
			Vector12d affine_force_node = Jx.transpose() * ext_force;

			affine_force[AB_index] += affine_force_node;
		}
	}
	
	
	int startIndex_grad = 0, startIndex_hess = 0;
	std::vector<double> AB_ine_energy_vec(tetSimMesh.num_meshes, 0);
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int AI = 0; AI < AB_ine_energy_vec.size(); AI++)
	{
		// reconstruct q_prev vector
		Vector12d qb_prev = constructQVec(tetSimMesh.translation_prev_ABD[AI], tetSimMesh.deformation_prev_ABD[AI]);
		// reconstruct qb_vel vector
		Vector12d qb_vel = constructQVec(tetSimMesh.translation_vel_ABD[AI], tetSimMesh.deformation_vel_ABD[AI]);
		// reconstruct q vector
		Vector12d qb = constructQVec(tetSimMesh.translation_ABD[AI], tetSimMesh.deformation_ABD[AI]);

		Vector12d qb_Minus_qbHat = qb - (qb_prev + parameters.dt * qb_vel + parameters.dt * parameters.dt * tetSimMesh.massMatrix_ABD[AI].inverse() * affine_force[AI]);

		Vector12d grad = tetSimMesh.massMatrix_ABD[AI] * qb_Minus_qbHat;
		Matrix12d hess = tetSimMesh.massMatrix_ABD[AI];

		for (int i = 0; i < 12; i++)
		{
			grad_triplet[startIndex_grad + AI * 12 + i] = { AI * 12 + i, grad[i] };
			for (int j = 0; j < 12; j++)
			{
				hessian_triplet[startIndex_grad + AI * 144 + i * 12 + j] = { 12 * AI + i, 12 * AI + j, hess(i,j)};
			}
		}
	}




	// contribution from affine energy
	startIndex_grad = tetSimMesh.num_meshes * 12 , startIndex_hess = tetSimMesh.num_meshes * 144;
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int AI = 0; AI < tetSimMesh.num_meshes; AI++)
	{
		Eigen::Matrix3d F = tetSimMesh.deformation_ABD[AI];
		Eigen::Matrix3d I_ = Eigen::Matrix3d::Identity();

		// original implementation
		Eigen::Matrix3d PK1 = parameters.dt * parameters.dt * tetSimMesh.volume_ABD[AI] * parameters.ABD_Coeff * 4.0 * F *
							  (F.transpose() * F - Eigen::Matrix3d::Identity());
		Vector9d grad_E_wrt_F = flatenMatrix3d(PK1);


		Matrix9d matrix_D = Matrix9d::Zero(); // Eq.4.27 in Kim_2022
		for (int i = 0; i < 3; i++)
		{
			Eigen::Vector3d ai = F.col(i);
			for (int j = 0; j < 3; j++)
			{
				Eigen::Vector3d aj = F.col(j);
				matrix_D.block(i * 3, j * 3, 3, 3) = aj * ai.transpose();
			}
		}
		Matrix9d matrix_H_II = 4.0 * (kroneckerProduct_matrices(I_, F * F.transpose()) + kroneckerProduct_matrices(F.transpose() * F, I_) + matrix_D);
		Matrix9d hess_E_wrt_F = parameters.dt * parameters.dt * tetSimMesh.volume_ABD[AI] * parameters.ABD_Coeff * 4.0 * (0.25 * matrix_H_II - Matrix9d::Identity());

		for (int i = 0; i < 9; i++)
		{
			grad_triplet[startIndex_grad + AI * 9 + i] = { AI * 12 + 3 + i, grad_E_wrt_F[i] };
			for (int j = 0; j < 9; j++)
			{
				hessian_triplet[startIndex_hess + AI * 81 + i * 9 + j] = { 12 * AI + 3 + i, 12 * AI + 3 + j, hess_E_wrt_F(i,j) };
			}
		}


		//// Linear ARAP implementation
		//Eigen::Matrix3d PK1 = parameters.dt * parameters.dt * tetSimMesh.volume_ABD[AI] * parameters.ABD_Coeff / 2 * (tetSimMesh.deformation_ABD[AI] 
		//					  + tetSimMesh.deformation_ABD[AI].transpose() - 2 * Eigen::Matrix3d::Identity());;
		//Vector9d grad_E_wrt_F = flatenMatrix3d(PK1);

		////std::cout << "grad_E_wrt_F = " << grad_E_wrt_F.transpose() << std::endl;
		//Eigen::Matrix<double, 9, 9> hess_E_wrt_F = parameters.dt * parameters.dt * tetSimMesh.volume_ABD[AI] * parameters.ABD_Coeff * Eigen::Matrix<double, 9, 9>::Identity();
		//for (int i = 0; i < 9; i++)
		//{
		//	grad_triplet[startIndex_grad + AI * 9 + i] = { AI * 12 + 3 + i, grad_E_wrt_F[i] };
		//	for (int j = 0; j < 9; j++)
		//	{
		//		hessian_triplet[startIndex_hess + AI * 81 + i * 9 + j] = { 12 * AI + 3 + i, 12 * AI + 3 + j, hess_E_wrt_F(i,j) };
		//	}
		//}

	}




	// calculate the barrier gradient and hessian
	{
		startIndex_grad += tetSimMesh.num_meshes * 9,
			startIndex_hess += tetSimMesh.num_meshes * 81;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < PG_PG.size(); i++)
		{
			int ptInd = PG_PG[i][1];
			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
			int actualStartIndex_hess = startIndex_hess + i * 144,
				actualStartIndex_grad = startIndex_grad + i * 12;
			double z2 = P[2] * P[2];
			Ground::gradAndHess(
				hessian_triplet,
				grad_triplet,
				actualStartIndex_hess,
				actualStartIndex_grad,
				tetSimMesh.pos_node,
				tetSimMesh.pos_node_Rest,
				tetSimMesh.index_node,
				ptInd,
				z2,
				tetSimMesh.surfaceInfo.boundaryVertices_area[ptInd],
				parameters, 
				true);
		}

		startIndex_grad += PG_PG.size() * 12, startIndex_hess += PG_PG.size() * 144;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < PT_PP.size(); i++)
		{
			Vector5i cont_PT = PT_PP[i];

			int ptInd = cont_PT[1];
			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
			Eigen::Vector3i tri = tetSimMesh.surfaceInfo.boundaryTriangles[cont_PT[2]];
			Eigen::Vector3d A = tetSimMesh.pos_node[tri[0]];
			Eigen::Vector3d B = tetSimMesh.pos_node[tri[1]];
			Eigen::Vector3d C = tetSimMesh.pos_node[tri[2]];
			int type = cont_PT[3];

			double dis2 = 0;
			DIS::computePointTriD(P, A, B, C, dis2);
			Eigen::Vector4i ptIndices = { ptInd , tri[0] , tri[1] , tri[2] };
			int actualStartIndex_hess = startIndex_hess + i * 144 * 4,
				actualStartIndex_grad = startIndex_grad + i * 24;
			double contactArea = tetSimMesh.surfaceInfo.boundaryVertices_area[ptInd];
			BarrierEnergy::gradAndHess_PT(
				contactArea,
				hessian_triplet,
				grad_triplet,
				actualStartIndex_hess,
				actualStartIndex_grad,
				ptIndices,
				type,
				dis2,
				tetSimMesh.pos_node,
				tetSimMesh.pos_node_Rest,
				tetSimMesh.index_node,
				parameters,
				true);

		}

		startIndex_grad += PT_PP.size() * 24, startIndex_hess += PT_PP.size() * 144 * 4;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < PT_PE.size(); i++)
		{
			Vector5i cont_PT = PT_PE[i];

			int ptInd = cont_PT[1];
			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
			Eigen::Vector3i tri = tetSimMesh.surfaceInfo.boundaryTriangles[cont_PT[2]];
			Eigen::Vector3d A = tetSimMesh.pos_node[tri[0]];
			Eigen::Vector3d B = tetSimMesh.pos_node[tri[1]];
			Eigen::Vector3d C = tetSimMesh.pos_node[tri[2]];
			int type = cont_PT[3];

			double dis2 = 0;
			DIS::computePointTriD(P, A, B, C, dis2);
			Eigen::Vector4i ptIndices = { ptInd , tri[0] , tri[1] , tri[2] };
			int actualStartIndex_hess = startIndex_hess + i * 144 * 9,
				actualStartIndex_grad = startIndex_grad + i * 36;
			double contactArea = tetSimMesh.surfaceInfo.boundaryVertices_area[ptInd];
			BarrierEnergy::gradAndHess_PT(
				contactArea,
				hessian_triplet,
				grad_triplet,
				actualStartIndex_hess,
				actualStartIndex_grad,
				ptIndices,
				type,
				dis2,
				tetSimMesh.pos_node,
				tetSimMesh.pos_node_Rest,
				tetSimMesh.index_node,
				parameters,
				true);

		}

		startIndex_grad += PT_PE.size() * 36, startIndex_hess += PT_PE.size() * 144 * 9;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < PT_PT.size(); i++)
		{
			Vector5i cont_PT = PT_PT[i];

			int ptInd = cont_PT[1];
			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
			Eigen::Vector3i tri = tetSimMesh.surfaceInfo.boundaryTriangles[cont_PT[2]];
			Eigen::Vector3d A = tetSimMesh.pos_node[tri[0]];
			Eigen::Vector3d B = tetSimMesh.pos_node[tri[1]];
			Eigen::Vector3d C = tetSimMesh.pos_node[tri[2]];
			int type = cont_PT[3];

			double dis2 = 0;
			DIS::computePointTriD(P, A, B, C, dis2);
			Eigen::Vector4i ptIndices = { ptInd , tri[0] , tri[1] , tri[2] };
			int actualStartIndex_hess = startIndex_hess + i * 144 * 16,
				actualStartIndex_grad = startIndex_grad + i * 48;
			double contactArea = tetSimMesh.surfaceInfo.boundaryVertices_area[ptInd];
			BarrierEnergy::gradAndHess_PT(
				contactArea,
				hessian_triplet,
				grad_triplet,
				actualStartIndex_hess,
				actualStartIndex_grad,
				ptIndices,
				type,
				dis2,
				tetSimMesh.pos_node,
				tetSimMesh.pos_node_Rest,
				tetSimMesh.index_node,
				parameters,
				true);

		}

		startIndex_grad += PT_PT.size() * 48, startIndex_hess += PT_PT.size() * 144 * 16;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < EE_EE.size(); i++)
		{
			Vector5i cont_EE = EE_EE[i];
			int E1 = cont_EE[1], E2 = cont_EE[2];
			int e1p1 = tetSimMesh.surfaceInfo.index_boundaryEdge[E1][0], e1p2 = tetSimMesh.surfaceInfo.index_boundaryEdge[E1][1],
				e2p1 = tetSimMesh.surfaceInfo.index_boundaryEdge[E2][0], e2p2 = tetSimMesh.surfaceInfo.index_boundaryEdge[E2][1];

			Eigen::Vector3d P1 = tetSimMesh.pos_node[e1p1];
			Eigen::Vector3d P2 = tetSimMesh.pos_node[e1p2];
			Eigen::Vector3d Q1 = tetSimMesh.pos_node[e2p1];
			Eigen::Vector3d Q2 = tetSimMesh.pos_node[e2p2];
			int type = cont_EE[3];

			double dis2 = 0;
			DIS::computeEdgeEdgeD(P1, P2, Q1, Q2, dis2);
			Eigen::Vector4i ptIndices = { e1p1 , e1p2 , e2p1 , e2p2 };
			int actualStartIndex_hess = startIndex_hess + i * 144 * 16,
				actualStartIndex_grad = startIndex_grad + i * 48;

			int emin = std::min(ptIndices[0], ptIndices[1]), emax = std::max(ptIndices[0], ptIndices[1]);
			double contactArea = tetSimMesh.surfaceInfo.boundaryEdges_area[emin][emax];

			BarrierEnergy::gradAndHess_EE(
				contactArea,
				hessian_triplet,
				grad_triplet,
				actualStartIndex_hess,
				actualStartIndex_grad,
				ptIndices,
				type,
				dis2,
				tetSimMesh.pos_node,
				tetSimMesh.pos_node_Rest,
				tetSimMesh.index_node,
				parameters,
				true);

		}

	}



	//double endTime3 = omp_get_wtime();
	//std::cout << "	Cal grad_hess Time is : " << endTime3 - endTime2 << "s" << std::endl;

	// assemable the left-hand side 
	Eigen::SparseMatrix<double> leftHandSide(tetSimMesh.num_meshes * 12, tetSimMesh.num_meshes * 12);
	leftHandSide.setFromTriplets(hessian_triplet.begin(), hessian_triplet.end());

	// assemable the right-hand side 
	Eigen::VectorXd rightHandSide = Eigen::VectorXd(tetSimMesh.num_meshes * 12);
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

	

#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int i = 0; i < tetSimMesh.num_meshes; i++)
	{
		movingDir[i] = result.block<12, 1>(12 * i, 0);
	}



	return movingDir;
}

// move points' position according to the direction; Note this is a trial movement
void step_forward_ABD(FEMParamters& parameters, Mesh_ABD& tetSimMesh, std::vector<Eigen::Vector3d>& current_ABD_translation,
	std::vector<Eigen::Matrix3d>& current_ABD_deformation, std::vector<Vector12d>& ABD_direction, std::vector<Eigen::Vector3d>& currentPosition, double step)
{

	// update the ABD state
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int vI = 0; vI < tetSimMesh.num_meshes; vI++)
	{
		tetSimMesh.translation_ABD[vI] = current_ABD_translation[vI] + step * ABD_direction[vI].block(0, 0, 3, 1);

		Eigen::Matrix3d trial_ABD_deformation = Eigen::Matrix3d::Zero();
		trial_ABD_deformation.col(0) = ABD_direction[vI].block(3, 0, 3, 1);
		trial_ABD_deformation.col(1) = ABD_direction[vI].block(6, 0, 3, 1);
		trial_ABD_deformation.col(2) = ABD_direction[vI].block(9, 0, 3, 1);

		tetSimMesh.deformation_ABD[vI] = current_ABD_deformation[vI] + step * trial_ABD_deformation;
	}

	// update the actual vertex position
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int vI = 0; vI < tetSimMesh.pos_node.size(); vI++)
	{
		Eigen::Matrix<double, 3, 12> Jx = build_Jx_matrix_for_ABD(tetSimMesh.pos_node_Rest[vI].transpose());
		int ABD_index = tetSimMesh.index_node[vI][0];
		tetSimMesh.pos_node[vI] = currentPosition[vI] + step * Jx * ABD_direction[ABD_index];
	}


}

// convert the ABD state update to position update direction
void convert_to_position_direction(FEMParamters& parameters, Mesh_ABD& tetSimMesh, std::vector<Vector12d>& direction_ABD, std::vector<Eigen::Vector3d>& position_direction)
{

#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int vI = 0; vI < tetSimMesh.pos_node.size(); vI++)
	{
		Eigen::Matrix<double, 3, 12> Jx = build_Jx_matrix_for_ABD(tetSimMesh.pos_node_Rest[vI].transpose());
		int ABD_index = tetSimMesh.index_node[vI][0];

		//Eigen::Vector3d current_position = Jx * direction_ABD[ABD_index];
		//position_direction[vI] = current_position - tetSimMesh.pos_node_prev[vI];


		position_direction[vI] = Jx * direction_ABD[ABD_index];
	}
};