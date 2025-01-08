

#include "simulator.h"


// implicit integration
void implicitFEM_ABD_triMesh(triMesh& triSimMesh, FEMParamters& parameters)
{

	for (int timestep = 0; timestep < parameters.num_timesteps; timestep++)
	{
		std::cout << "timestep = " << timestep << std::endl;

		if (timestep % parameters.outputFrequency == 0)
		{
			triSimMesh.exportSurfaceMesh("surfMesh", timestep);
		}

		for (int i = 0; i < triSimMesh.allObjects.size(); i++)
		{
			triSimMesh.allObjects[i].translation_prev_ABD = triSimMesh.allObjects[i].translation_ABD;
			triSimMesh.allObjects[i].deformation_prev_ABD = triSimMesh.allObjects[i].deformation_ABD;
			triSimMesh.allObjects[i].affine_prev = triSimMesh.allObjects[i].affine;


			triSimMesh.allObjects[i].pos_node_interior_prev = triSimMesh.allObjects[i].pos_node_interior;
			triSimMesh.allObjects[i].pos_node_surface_prev = triSimMesh.allObjects[i].pos_node_surface;
		}


		std::map<int, std::vector<Vector6d>> broken_objects;
		contact_Info contact_pairs;
		find_contact(triSimMesh,parameters,contact_pairs);

		double lastEnergyVal = compute_IP_energy_ABD_triMesh(triSimMesh, parameters, timestep, contact_pairs);


		for (int ite = 0; ite < 15; ite++)
		{

			std::cout << "		ite = " << ite << "; lastEnergyVal = " << lastEnergyVal << std::endl;
			std::vector<Vector12d> direction_ABD;
			
			solve_linear_system_ABD_triMesh(triSimMesh, parameters, timestep, direction_ABD, broken_objects, contact_pairs, ite);
			if (broken_objects.size() != 0 && timestep < 3)
			{
				break;
			}
			double dist_to_converge = calculate_maximum_velocity(parameters, triSimMesh, direction_ABD);




			std::cout << std::scientific << std::setprecision(4) << "			dist_to_converge = "
				<< dist_to_converge / parameters.dt << "m/s; threshold = " << parameters.searchResidual
				<< "m/s" << std::endl;
			if (ite && (dist_to_converge / parameters.dt < parameters.searchResidual))
			{
				break;
			}


			double startTime1, endTime1, endTime2;
			startTime1 = omp_get_wtime();

			std::cout << "			Calculate step size;" << std::endl;


			find_contact(triSimMesh, parameters, contact_pairs, direction_ABD);

			double step = calMaxStep(
				parameters,
				triSimMesh,
				contact_pairs,
				timestep);

			endTime1 = omp_get_wtime();


			std::cout << "			Step forward = " << step << std::endl;
			step_forward_ABD_triMesh(parameters, triSimMesh, direction_ABD, step);
			find_contact(triSimMesh, parameters, contact_pairs);
			double newEnergyVal = compute_IP_energy_ABD_triMesh(triSimMesh, parameters, timestep, contact_pairs);
			std::cout << std::scientific << std::setprecision(4) << "			step = "
				<< step << "; newEnergyVal = " << newEnergyVal  << std::endl;


			while (newEnergyVal > lastEnergyVal && step >= 1.0e-5)
			{
				step /= 2.0;

				step_forward_ABD_triMesh(parameters, triSimMesh, direction_ABD, step);
				find_contact(triSimMesh, parameters, contact_pairs);
				newEnergyVal = compute_IP_energy_ABD_triMesh(triSimMesh, parameters, timestep, contact_pairs);
				std::cout << "				step = " << step << "; newEnergyVal = "
					<< newEnergyVal << std::endl;
				if (std::abs(newEnergyVal - lastEnergyVal) / lastEnergyVal < 0.001) // traped in the local mimima
				{
					break;
				}
			}

			// The energy has been greatly minimized. It is time to stop
			if (std::abs(newEnergyVal) / lastEnergyVal < 0.000001)
			{
				break;
			}
			lastEnergyVal = newEnergyVal;


			std::cout << "			lastEnergyVal = " << lastEnergyVal << std::endl << std::endl;



		}


		if (broken_objects.size() != 0 && timestep < 3)
		{
			std::map<int, objMeshFormat> crackSurface_object;
			fracture_sim(parameters, triSimMesh, broken_objects, crackSurface_object);
			cut_object_with_cracks(parameters, triSimMesh, crackSurface_object);
			triSimMesh.updateGlobalSimulationTriMesh_ABD();
			continue;
		}


		// update the velocity of the ABD objects
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < triSimMesh.allObjects.size(); i++)
		{
			triSimMesh.allObjects[i].translation_vel_ABD = (triSimMesh.allObjects[i].translation_ABD - triSimMesh.allObjects[i].translation_prev_ABD) / parameters.dt;
			triSimMesh.allObjects[i].deformation_vel_ABD = (triSimMesh.allObjects[i].deformation_ABD - triSimMesh.allObjects[i].deformation_prev_ABD) / parameters.dt;
		}


	}

}

// compute the incremental potential energy
double compute_IP_energy_ABD_triMesh(triMesh& triSimMesh, FEMParamters& parameters, int timestep, contact_Info& contact_pairs)
{
	double energyVal = 0;	

	// initeria energy contribution per affine body
	std::vector<Vector12d> affine_force(triSimMesh.num_meshes, Vector12d::Zero());
	if (parameters.gravity.norm() != 0)
	{
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int obj = 0; obj < triSimMesh.allObjects.size(); obj++)
		{
			for (int i = 0; i < triSimMesh.allObjects[obj].pos_node_interior.size(); i++) // interior boundary conditions
			{
				Eigen::Matrix<double, 3, 12> Jx = build_Jx_matrix_for_ABD(triSimMesh.allObjects[obj].pos_node_Rest_interior[i].transpose());
				Vector12d affine_force_node = Jx.transpose() * parameters.gravity;

				affine_force[obj] += affine_force_node;
			}
		}

	}
	//for (int i = 0; i < triSimMesh.pos_node_surface.size(); i++) // surface boundary conditions
	//{
	//	if (triSimMesh.boundaryCondition_node_surface[i].type == 2)
	//	{
	//		Eigen::Matrix<double, 3, 12> Jx = build_Jx_matrix_for_ABD(triSimMesh.pos_node_Rest_surface[i].transpose());

	//		int AB_index = triSimMesh.index_node_surface[i][0];
	//		Eigen::Vector3d ext_force = compute_external_force(triSimMesh.boundaryCondition_node_surface, i, timestep);
	//		Vector12d affine_force_node = Jx.transpose() * ext_force;

	//		affine_force[AB_index] += affine_force_node;
	//	}
	//}
	std::vector<double> AB_ine_energy_vec(triSimMesh.num_meshes,0);
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int AI = 0; AI < AB_ine_energy_vec.size(); AI++)
	{
		// reconstruct q_prev vector
		Vector12d qb_prev = constructQVec(triSimMesh.allObjects[AI].translation_prev_ABD, triSimMesh.allObjects[AI].deformation_prev_ABD);
		// reconstruct qb_vel vector
		Vector12d qb_vel = constructQVec(triSimMesh.allObjects[AI].translation_vel_ABD, triSimMesh.allObjects[AI].deformation_vel_ABD);
		// reconstruct q vector
		Vector12d qb = constructQVec(triSimMesh.allObjects[AI].translation_ABD, triSimMesh.allObjects[AI].deformation_ABD);

		Vector12d qb_Minus_qbHat = qb - (qb_prev + parameters.dt * qb_vel + parameters.dt * parameters.dt * triSimMesh.allObjects[AI].massMatrix_ABD.inverse() * affine_force[AI]);
		AB_ine_energy_vec[AI] = 0.5 * qb_Minus_qbHat.transpose() * triSimMesh.allObjects[AI].massMatrix_ABD * qb_Minus_qbHat;


	}	
	energyVal += std::accumulate(AB_ine_energy_vec.begin(), AB_ine_energy_vec.end(), 0.0);


	// affine energy contribution per affine body
	std::vector<double> AB_aff_energy_vec(triSimMesh.num_meshes,0);
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int AI = 0; AI < AB_aff_energy_vec.size(); AI++)
	{
		// original implementation
		AB_aff_energy_vec[AI] = parameters.dt * parameters.dt * parameters.ABD_Coeff * triSimMesh.allObjects[AI].volume
			* (triSimMesh.allObjects[AI].deformation_ABD.transpose() * triSimMesh.allObjects[AI].deformation_ABD - Eigen::Matrix3d::Identity()).squaredNorm();

	}
	energyVal += std::accumulate(AB_aff_energy_vec.begin(), AB_aff_energy_vec.end(), 0.0);


	// energy contribution from barrier
	energyVal += compute_Barrier_energy(
		triSimMesh,
		parameters,
		contact_pairs,
		timestep);


	return energyVal;
}

// compute energy gradient and hessian of the linear system and solve the syetem
void solve_linear_system_ABD_triMesh(
	triMesh& triSimMesh, 
	FEMParamters& parameters, 
	int timestep, 
	std::vector<Vector12d>& movingDir,
	std::map<int, std::vector<Vector6d>>& broken_objects,
	contact_Info& contact_pairs,
	const int& iteration)
{
	// we only need to update boundary vertices' position to calculate the barrier energy update
	movingDir.clear();
	movingDir.resize(triSimMesh.num_meshes, Vector12d::Zero());



	// contribution from contact
	//		calculate contacts




	int gradSize = triSimMesh.num_meshes * 12 + triSimMesh.num_meshes * 9 + contact_pairs.Point_Ground.size() * 12
		+ contact_pairs.Point_Triangle_PP_index.size() * 24 + contact_pairs.Point_Triangle_PE_index.size() * 36 
		+ contact_pairs.Point_Triangle_PT_index.size() * 48 + contact_pairs.Edge_Edge.size() * 48;
	int hessSize = triSimMesh.num_meshes * 144 + triSimMesh.num_meshes * 81 + contact_pairs.Point_Ground.size() * 144
		+ contact_pairs.Point_Triangle_PP_index.size() * 144 * 4 + contact_pairs.Point_Triangle_PE_index.size() * 144 * 9 
		+ contact_pairs.Point_Triangle_PT_index.size() * 144 * 16 + contact_pairs.Edge_Edge.size() * 144 * 16;
	std::vector<std::pair<int, double>> grad_triplet(gradSize, std::make_pair(0,0.0));
	std::vector<Eigen::Triplet<double>> hessian_triplet(hessSize, Eigen::Triplet<double>(0, 0, 0.0));



	// contribution from inertia term 
	std::vector<Vector12d> affine_force(triSimMesh.allObjects.size(), Vector12d::Zero());
	if (parameters.gravity.norm() != 0)
	{
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int obj = 0; obj < triSimMesh.allObjects.size(); obj++)
		{
			for (int i = 0; i < triSimMesh.allObjects[obj].pos_node_interior.size(); i++) // interior boundary conditions
			{
				Eigen::Matrix<double, 3, 12> Jx = build_Jx_matrix_for_ABD(triSimMesh.allObjects[obj].pos_node_Rest_interior[i].transpose());
				Vector12d affine_force_node = Jx.transpose() * parameters.gravity;

				affine_force[obj] += affine_force_node;
			}
		}


	}

	int startIndex_grad = 0, startIndex_hess = 0;
	std::vector<double> AB_ine_energy_vec(triSimMesh.allObjects.size(), 0);
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int AI = 0; AI < AB_ine_energy_vec.size(); AI++)
	{
		// reconstruct q_prev vector
		Vector12d qb_prev = constructQVec(triSimMesh.allObjects[AI].translation_prev_ABD, triSimMesh.allObjects[AI].deformation_prev_ABD);
		// reconstruct qb_vel vector
		Vector12d qb_vel = constructQVec(triSimMesh.allObjects[AI].translation_vel_ABD, triSimMesh.allObjects[AI].deformation_vel_ABD);
		// reconstruct q vector
		Vector12d qb = constructQVec(triSimMesh.allObjects[AI].translation_ABD, triSimMesh.allObjects[AI].deformation_ABD);

		Vector12d qb_Minus_qbHat = qb - (qb_prev + parameters.dt * qb_vel + parameters.dt * parameters.dt * triSimMesh.allObjects[AI].massMatrix_ABD.inverse() * affine_force[AI]);


		Vector12d grad = triSimMesh.allObjects[AI].massMatrix_ABD * qb_Minus_qbHat;
		Matrix12d hess = triSimMesh.allObjects[AI].massMatrix_ABD;

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
	startIndex_grad = triSimMesh.num_meshes * 12 , startIndex_hess = triSimMesh.num_meshes * 144;
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int AI = 0; AI < triSimMesh.allObjects.size(); AI++)
	{
		Eigen::Matrix3d F = triSimMesh.allObjects[AI].deformation_ABD;
		Eigen::Matrix3d I_ = Eigen::Matrix3d::Identity();

		// original implementation
		Eigen::Matrix3d PK1 = parameters.dt * parameters.dt * triSimMesh.allObjects[AI].volume * parameters.ABD_Coeff * 4.0 * F *
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
		Matrix9d hess_E_wrt_F = parameters.dt * parameters.dt * triSimMesh.allObjects[AI].volume * parameters.ABD_Coeff * 4.0 * (0.25 * matrix_H_II - Matrix9d::Identity());

		for (int i = 0; i < 9; i++)
		{
			grad_triplet[startIndex_grad + AI * 9 + i] = { AI * 12 + 3 + i, grad_E_wrt_F[i] };
			for (int j = 0; j < 9; j++)
			{
				hessian_triplet[startIndex_hess + AI * 81 + i * 9 + j] = { 12 * AI + 3 + i, 12 * AI + 3 + j, hess_E_wrt_F(i,j) };
			}
		}


		//// Linear ARAP implementation
		//Eigen::Matrix3d PK1 = parameters.dt * parameters.dt * triSimMesh.volume_ABD[AI] * parameters.ABD_Coeff / 2 * (triSimMesh.deformation_ABD[AI] 
		//					  + triSimMesh.deformation_ABD[AI].transpose() - 2 * Eigen::Matrix3d::Identity());;
		//Vector9d grad_E_wrt_F = flatenMatrix3d(PK1);

		////std::cout << "grad_E_wrt_F = " << grad_E_wrt_F.transpose() << std::endl;
		//Eigen::Matrix<double, 9, 9> hess_E_wrt_F = parameters.dt * parameters.dt * triSimMesh.volume_ABD[AI] * parameters.ABD_Coeff * Eigen::Matrix<double, 9, 9>::Identity();
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
	if (contact_pairs.Point_Ground.size() + contact_pairs.Point_Triangle_PP_index.size() + contact_pairs.Point_Triangle_PE_index.size() 
		+ contact_pairs.Point_Triangle_PT_index.size() + contact_pairs.Edge_Edge.size() != 0)
	{


		startIndex_grad += triSimMesh.num_meshes * 9,
			startIndex_hess += triSimMesh.num_meshes * 81;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < contact_pairs.Point_Ground.size(); i++)
		{
			Vector2i PG_Contact = contact_pairs.Point_Ground[i];

			int actualStartIndex_hess = startIndex_hess + i * 144,
				actualStartIndex_grad = startIndex_grad + i * 12;
			std::vector<Eigen::Vector3d> contactForce_node;


			Ground::gradAndHess(
				hessian_triplet,
				grad_triplet,
				actualStartIndex_hess,
				actualStartIndex_grad,
				PG_Contact,
				triSimMesh,
				contactForce_node,
				parameters);


		}

		startIndex_grad += contact_pairs.Point_Ground.size() * 12, startIndex_hess += contact_pairs.Point_Ground.size() * 144;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < contact_pairs.Point_Triangle_PP_index.size(); i++)
		{
			int ct_index = contact_pairs.Point_Triangle_PP_index[i];
			Vector5i PT_Contact = contact_pairs.Point_Triangle[ct_index];
			double dis2 = contact_pairs.Point_Triangle_Dis[ct_index];
			std::vector<Eigen::Vector3d> contactForce_node;

			int actualStartIndex_hess = startIndex_hess + i * 144 * 4,
				actualStartIndex_grad = startIndex_grad + i * 24;
			BarrierEnergy::gradAndHess_PT(hessian_triplet,
				grad_triplet,
				actualStartIndex_hess,
				actualStartIndex_grad,
				dis2,
				PT_Contact,
				triSimMesh,
				contactForce_node,
				parameters);		
		}

		startIndex_grad += contact_pairs.Point_Triangle_PP_index.size() * 24, startIndex_hess += contact_pairs.Point_Triangle_PP_index.size() * 144 * 4;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < contact_pairs.Point_Triangle_PE_index.size(); i++)
		{
			int ct_index = contact_pairs.Point_Triangle_PE_index[i];
			Vector5i PT_Contact = contact_pairs.Point_Triangle[ct_index];
			double dis2 = contact_pairs.Point_Triangle_Dis[ct_index];
			std::vector<Eigen::Vector3d> contactForce_node;

			int actualStartIndex_hess = startIndex_hess + i * 144 * 9,
				actualStartIndex_grad = startIndex_grad + i * 36;
			BarrierEnergy::gradAndHess_PT(hessian_triplet,
				grad_triplet,
				actualStartIndex_hess,
				actualStartIndex_grad,
				dis2,
				PT_Contact,
				triSimMesh,
				contactForce_node,
				parameters);
		}

		startIndex_grad += contact_pairs.Point_Triangle_PE_index.size() * 36, startIndex_hess += contact_pairs.Point_Triangle_PE_index.size() * 144 * 9;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < contact_pairs.Point_Triangle_PT_index.size(); i++)
		{
			int ct_index = contact_pairs.Point_Triangle_PT_index[i];
			Vector5i PT_Contact = contact_pairs.Point_Triangle[ct_index];
			double dis2 = contact_pairs.Point_Triangle_Dis[ct_index];
			std::vector<Eigen::Vector3d> contactForce_node;

			int actualStartIndex_hess = startIndex_hess + i * 144 * 16,
				actualStartIndex_grad = startIndex_grad + i * 48;
			BarrierEnergy::gradAndHess_PT(hessian_triplet,
				grad_triplet,
				actualStartIndex_hess,
				actualStartIndex_grad,
				dis2,
				PT_Contact,
				triSimMesh,
				contactForce_node,
				parameters);
		}

		startIndex_grad += contact_pairs.Point_Triangle_PT_index.size() * 48, startIndex_hess += contact_pairs.Point_Triangle_PT_index.size() * 144 * 16;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < contact_pairs.Edge_Edge.size(); i++)
		{
			Vector5i EE_Contact = contact_pairs.Edge_Edge[i];
			if (EE_Contact[4] > 0)
			{
				double dis2 = contact_pairs.Edge_Edge_Dis[i];
				std::vector<Eigen::Vector3d> contactForce_node;

				int actualStartIndex_hess = startIndex_hess + i * 144 * 16,
					actualStartIndex_grad = startIndex_grad + i * 48;


				BarrierEnergy::gradAndHess_EE(
					hessian_triplet,
					grad_triplet,
					actualStartIndex_hess,
					actualStartIndex_grad,
					EE_Contact,
					dis2,
					triSimMesh,
					contactForce_node,
					parameters);

			}
				

		}


		// check if any object will be broken
		if (iteration == -1)
		{
			if_start_fracture_sim(triSimMesh, broken_objects);
			if (broken_objects.size() != 0)
			{
				return;
			}
		}

	}





	//double endTime3 = omp_get_wtime();
	//std::cout << "	Cal grad_hess Time is : " << endTime3 - endTime2 << "s" << std::endl;

	// assemable the left-hand side 
	Eigen::SparseMatrix<double> leftHandSide(triSimMesh.num_meshes * 12, triSimMesh.num_meshes * 12);
	leftHandSide.setFromTriplets(hessian_triplet.begin(), hessian_triplet.end());

	// assemable the right-hand side 
	Eigen::VectorXd rightHandSide = Eigen::VectorXd(triSimMesh.num_meshes * 12);
	rightHandSide.setZero();
	for (int i = 0; i < grad_triplet.size(); i++)
	{		
		std::pair<int, double> ele = grad_triplet[i];
		rightHandSide[ele.first] += ele.second;
	}

	//double endTime4 = omp_get_wtime();
	//std::cout << "	Assemble grad_hess Time is : " << endTime4 - endTime3 << "s" << std::endl;

	//std::cout << "leftHandSide = " << leftHandSide << std::endl;
	//std::cout << "rightHandSide = " << rightHandSide << std::endl;

	leftHandSide.makeCompressed();
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::IncompleteLUT<double>> solver;
	//Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
	solver.compute(leftHandSide);
	Eigen::VectorXd result = solver.solve(-rightHandSide);
	solver.setTolerance(1e-12); // 设置容差
	solver.setMaxIterations(10000); // 设置更大的最大迭代次数

	std::cout << "			Newton solver iterations: " << solver.iterations() << "; error: " << solver.error() << std::endl;


	//std::cout << "result = " << result << std::endl;

#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int i = 0; i < triSimMesh.num_meshes; i++)
	{
		movingDir[i] = result.block<12, 1>(12 * i, 0);
	}


}

// move points' position according to the direction; Note this is a trial movement
void step_forward_ABD_triMesh(FEMParamters& parameters, triMesh& triSimMesh, std::vector<Vector12d>& ABD_direction, double step)
{

	// update the ABD state
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int vI = 0; vI < triSimMesh.allObjects.size(); vI++)
	{
		triSimMesh.allObjects[vI].translation_ABD = triSimMesh.allObjects[vI].translation_prev_ABD + step * ABD_direction[vI].block(0, 0, 3, 1);

		Eigen::Matrix3d trial_ABD_deformation = Eigen::Matrix3d::Zero();
		trial_ABD_deformation.col(0) = ABD_direction[vI].block(3, 0, 3, 1);
		trial_ABD_deformation.col(1) = ABD_direction[vI].block(6, 0, 3, 1);
		trial_ABD_deformation.col(2) = ABD_direction[vI].block(9, 0, 3, 1);

		triSimMesh.allObjects[vI].deformation_ABD = triSimMesh.allObjects[vI].deformation_prev_ABD + step * trial_ABD_deformation;

		triSimMesh.allObjects[vI].affine = triSimMesh.allObjects[vI].affine_prev + step * ABD_direction[vI];
	}

	// update the actual vertex position on the surface
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int vI = 0; vI < triSimMesh.allObjects.size(); vI++)
	{
		for (int j = 0; j < triSimMesh.allObjects[vI].pos_node_surface.size(); j++)
		{
			Eigen::Matrix<double, 3, 12> Jx = build_Jx_matrix_for_ABD(triSimMesh.allObjects[vI].objectSurfaceMesh.vertices[j].transpose());
			triSimMesh.allObjects[vI].pos_node_surface[j] = triSimMesh.allObjects[vI].pos_node_surface_prev[j] +step * Jx * ABD_direction[vI];
		}
	}


	// update the actual vertex position in the interior
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int vI = 0; vI < triSimMesh.allObjects.size(); vI++)
	{
		for (int j = 0; j < triSimMesh.allObjects[vI].pos_node_interior.size(); j++)
		{
			Eigen::Matrix<double, 3, 12> Jx = build_Jx_matrix_for_ABD(triSimMesh.allObjects[vI].pos_node_interior[j].transpose());
			triSimMesh.allObjects[vI].pos_node_interior[j] = triSimMesh.allObjects[vI].pos_node_interior_prev[j] + step * Jx * ABD_direction[vI];
		}
	}



}

// convert the ABD state update to position update direction
double calculate_maximum_velocity(FEMParamters& parameters, triMesh& triSimMesh, std::vector<Vector12d>& direction_ABD)
{
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int obj = 0; obj < triSimMesh.allObjects.size(); obj++)
	{
		for (int vI = 0; vI < triSimMesh.allObjects[obj].pos_node_surface.size(); vI++)
		{
			Eigen::Matrix<double, 3, 12> Jx = build_Jx_matrix_for_ABD(triSimMesh.allObjects[obj].objectSurfaceMesh.vertices[vI].transpose());
			triSimMesh.allObjects[obj].pos_node_surface_direction[vI] = Jx * direction_ABD[obj];
		}
	}

	std::vector<double> max_vel(triSimMesh.allObjects.size(),0);
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int obj = 0; obj < triSimMesh.allObjects.size(); obj++)
	{
		max_vel[obj] = infiniteNorm(triSimMesh.allObjects[obj].pos_node_surface_direction);
	}

	double maxValue = *std::max_element(max_vel.begin(), max_vel.end());
	
	return maxValue;

};


void if_start_fracture_sim(triMesh& triSimMesh, std::map<int, std::vector<Vector6d>>& broken_objects)
{
	broken_objects.clear();
	for (int i = 0; i < triSimMesh.allObjects.size(); i++)
	{
		if (triSimMesh.allObjects[i].breakable == true)
		{
			//double max_force = -999999999.0;
			//std::vector<Vector6d> contactForce;
			//for (int j = triSimMesh.allObjects[i].objectSurfaceMeshes_node_start_end[0]; j < triSimMesh.allObjects[i].objectSurfaceMeshes_node_start_end[1]; j++)
			//{
			//	Vector6d contact_position_force;
			//	Eigen::Vector3d position = triSimMesh.pos_node_surface[j];
			//	Eigen::Vector3d force = triSimMesh.contactForce_node_surface[j];
			//	if (force.norm() > 0.0001) // to avoid numerical error
			//	{
			//		contact_position_force.block(0, 0, 3, 1) = position;
			//		contact_position_force.block(3, 0, 3, 1) = force * 100;
			//		contactForce.push_back(contact_position_force);
			//		//std::cout << "position = " << position.transpose() << "; " << "force = " << force.transpose() << std::endl;
			//	}
			//	max_force = std::max(force.norm(), max_force);
			//}
			//std::cout << "			max_force = " << max_force << std::endl;
			//if (max_force > triSimMesh.allObjects[i].objectMaterial.fracture_start_force )
			//{
			//	broken_objects[i] = contactForce;
			//}		
		}
	}

}

void fracture_sim(FEMParamters& parameters, triMesh& triSimMesh, std::map<int, std::vector<Vector6d>>& broken_objects, std::map<int, objMeshFormat>& crackSurface_object)
{
	crackSurface_object.clear(); // int: object's index

	triSimMesh.updateEachObjectSurfaceMesh();
	for (auto it = broken_objects.begin(); it != broken_objects.end(); it++)
	{
		int objectIndex = it->first;
		std::vector<Vector6d> contactForce = it->second;


		int num_particles = 10000;
		std::vector<Eigen::Vector3d> mpm_sim_particles = triSimMesh.allObjects[objectIndex].objectSurfaceMesh.sample_points_inside_mesh(num_particles);
		double per_particle_volume = triSimMesh.allObjects[objectIndex].objectSurfaceMesh.volume / (double)num_particles;
		mpmSimulator::MPMParamters mpm_parameters;
		mpm_parameters.dt = 1.0E-7;
		mpm_parameters.mat_mpm = triSimMesh.allObjects[objectIndex].objectMaterial;
		mpm_parameters.numOfThreads = parameters.numOfThreads;
		mpm_parameters.dx = 2.0 * std::cbrt(per_particle_volume);

		std::cout << "Start fracture simulation in object "<< objectIndex <<": " << triSimMesh.allObjects[objectIndex].objectNote << std::endl;

		std::tuple<bool, objMeshFormat, objMeshFormat, std::vector<objMeshFormat>> crack_res = mpmSimulator::crackSimulation(
			mpm_sim_particles,
			per_particle_volume,
			triSimMesh.allObjects[objectIndex].objectMaterial,
			mpm_parameters,
			contactForce,
			500);

		if (std::get<0>(crack_res))
		{
			// output crack surface and fragments
			objMeshFormat partial_cut = std::get<1>(crack_res);
			objMeshFormat full_cut = std::get<2>(crack_res);
			
			crackSurface_object[objectIndex] = full_cut;
		}

	}
}

void cut_object_with_cracks(FEMParamters& parameters, triMesh& triSimMesh, std::map<int, objMeshFormat>& crackSurface_object)
{
	triSimMesh.updateEachObjectSurfaceMesh();

	for (auto it = crackSurface_object.begin(); it != crackSurface_object.end(); it++)
	{
		int objectIndex = it->first;
		objMeshFormat full_cut = it->second;
		float voxel_size = static_cast<float>(parameters.IPC_dis);
		objMeshFormat full_cut_surf = full_cut.reconstruct_with_vdb(voxel_size);

		objMeshFormat children = triSimMesh.allObjects[objectIndex].objectSurfaceMesh.boolean_difference_with_mesh(full_cut_surf);
		children.outputFile("child", -999, false);
		children.triangulate();

		
		children.sepConnectedComponents();


		std::cout << "Object " << "_" << objectIndex << "_" << " is cut into " << std::to_string(children.componentsSep.size()) << " pieces!" << std::endl;

		ABD_Object parent_object = triSimMesh.allObjects[objectIndex];
		// conver the object into ABD_Object
		for (int cc = 0; cc < children.componentsSep.size(); cc++)
		{
			objMeshFormat childMesh = children.componentsSep[cc];
			childMesh.updateMesh();

			ABD_Object child_object;
			child_object.objectNote = parent_object.objectNote + "_F_" + std::to_string(cc);
			child_object.objectMaterial = parent_object.objectMaterial;
			child_object.breakable = parent_object.breakable;
			child_object.per_point_volume = parent_object.per_point_volume;
			child_object.translation_prev_ABD = parent_object.translation_prev_ABD;
			child_object.translation_vel_ABD = parent_object.translation_vel_ABD;
			child_object.translation_ABD = parent_object.translation_ABD;
			child_object.deformation_prev_ABD = parent_object.deformation_prev_ABD;
			child_object.deformation_vel_ABD = parent_object.deformation_vel_ABD;
			child_object.deformation_ABD = parent_object.deformation_ABD;

			child_object.objectSurfaceMesh = childMesh;
			child_object.volume = childMesh.volume;
			child_object.need_update_rest_position = true;

			// the first child replaces the parent object
			if (cc == 0)
			{
				triSimMesh.allObjects[objectIndex] = child_object;
			}
			else // remaining ones are appended to the vector
			{
				triSimMesh.allObjects.push_back(child_object);
			}
		}

	}

	std::cout << "We now have " << std::to_string(triSimMesh.allObjects.size()) << " objects in the scene!" << std::endl;

}

