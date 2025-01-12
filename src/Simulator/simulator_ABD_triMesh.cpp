

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
			
			//for (int i = 0; i < triSimMesh.allObjects.size(); i++)
			//{
			//	objMeshFormat test;
			//	test.vertices = triSimMesh.allObjects[i].pos_node_surface;
			//	test.faces = triSimMesh.allObjects[i].objectSurfaceMesh.faces;
			//	test.outputFile(triSimMesh.allObjects[i].objectNote, timestep);
			//}

		}

		for (int i = 0; i < triSimMesh.allObjects.size(); i++)
		{
			triSimMesh.allObjects[i].affine_prev = triSimMesh.allObjects[i].affine;
		}


		std::map<int, std::vector<Vector6d>> broken_objects;
		contact_Info contact_pairs;
		find_contact(triSimMesh,parameters,contact_pairs);

		double lastEnergyVal = compute_IP_energy_ABD_triMesh(triSimMesh, parameters, timestep, contact_pairs);


		for (int ite = 0; ite < 15; ite++)
		{


			for (int i = 0; i < triSimMesh.allObjects.size(); i++)
			{
				triSimMesh.allObjects[i].affine_last = triSimMesh.allObjects[i].affine;

				triSimMesh.allObjects[i].pos_node_interior_prev = triSimMesh.allObjects[i].pos_node_interior;
				triSimMesh.allObjects[i].pos_node_surface_prev = triSimMesh.allObjects[i].pos_node_surface;
			}


			std::cout << "		ite = " << ite << "; lastEnergyVal = " << lastEnergyVal << std::endl;
			std::vector<Vector12d> direction_ABD;
			
			solve_linear_system_ABD_triMesh(triSimMesh, parameters, timestep, direction_ABD, broken_objects, contact_pairs, ite);
			if (broken_objects.size() != 0 && timestep == 16)
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


		if (broken_objects.size() != 0 && timestep == 16)
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
			triSimMesh.allObjects[i].affine_vel = (triSimMesh.allObjects[i].affine - triSimMesh.allObjects[i].affine_prev) / parameters.dt;
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
		Vector12d qb_Minus_qbHat = triSimMesh.allObjects[AI].affine - (triSimMesh.allObjects[AI].affine_prev + parameters.dt * triSimMesh.allObjects[AI].affine_vel + parameters.dt * parameters.dt * triSimMesh.allObjects[AI].massMatrix_ABD.inverse() * affine_force[AI]);
		AB_ine_energy_vec[AI] = 0.5 * qb_Minus_qbHat.transpose() * triSimMesh.allObjects[AI].massMatrix_ABD * qb_Minus_qbHat;


	}	
	energyVal += std::accumulate(AB_ine_energy_vec.begin(), AB_ine_energy_vec.end(), 0.0);

	//std::cout << "Energy after ab = " << energyVal << std::endl;

	// affine energy contribution per affine body
	std::vector<double> AB_aff_energy_vec(triSimMesh.num_meshes,0);
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int AI = 0; AI < AB_aff_energy_vec.size(); AI++)
	{
		Eigen::Matrix3d deformation = Eigen::Matrix3d::Zero();
		deformation.col(0) = triSimMesh.allObjects[AI].affine.segment<3>(3).transpose();
		deformation.col(1) = triSimMesh.allObjects[AI].affine.segment<3>(6).transpose();
		deformation.col(2) = triSimMesh.allObjects[AI].affine.segment<3>(9).transpose();

	
		// original implementation
		AB_aff_energy_vec[AI] = parameters.dt * parameters.dt * parameters.ABD_Coeff * triSimMesh.allObjects[AI].objectSurfaceMesh.volume
			* (deformation.transpose() * deformation - Eigen::Matrix3d::Identity()).squaredNorm();

	}
	energyVal += std::accumulate(AB_aff_energy_vec.begin(), AB_aff_energy_vec.end(), 0.0);

	//std::cout << "Energy after aff = " << energyVal << std::endl;

	

	// energy contribution from barrier
	energyVal += compute_Barrier_energy(
		triSimMesh,
		parameters,
		contact_pairs,
		timestep);

	//std::cout << "Energy after contact = " << energyVal << std::endl;

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

		Vector12d qb_Minus_qbHat = triSimMesh.allObjects[AI].affine - (triSimMesh.allObjects[AI].affine_prev + parameters.dt * triSimMesh.allObjects[AI].affine_vel + parameters.dt * parameters.dt * triSimMesh.allObjects[AI].massMatrix_ABD.inverse() * affine_force[AI]);

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
		Eigen::Matrix3d F = Eigen::Matrix3d::Zero();
		F.col(0) = triSimMesh.allObjects[AI].affine.segment<3>(3).transpose();
		F.col(1) = triSimMesh.allObjects[AI].affine.segment<3>(6).transpose();
		F.col(2) = triSimMesh.allObjects[AI].affine.segment<3>(9).transpose();
		Eigen::Matrix3d I_ = Eigen::Matrix3d::Identity();

		// original implementation
		Eigen::Matrix3d PK1 = parameters.dt * parameters.dt * triSimMesh.allObjects[AI].objectSurfaceMesh.volume * parameters.ABD_Coeff * 4.0 * F *
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
		Matrix9d hess_E_wrt_F = parameters.dt * parameters.dt * triSimMesh.allObjects[AI].objectSurfaceMesh.volume * parameters.ABD_Coeff * 4.0 * (0.25 * matrix_H_II - Matrix9d::Identity());

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
		// Store the contact force
		std::vector<std::vector<Eigen::Vector3d>> contact_force_PG(contact_pairs.Point_Ground.size());
		std::vector<std::vector<Eigen::Vector3d>> contact_force_PP(contact_pairs.Point_Triangle_PP_index.size());
		std::vector<std::vector<Eigen::Vector3d>> contact_force_PE(contact_pairs.Point_Triangle_PE_index.size());
		std::vector<std::vector<Eigen::Vector3d>> contact_force_PT(contact_pairs.Point_Triangle_PT_index.size());
		std::vector<std::vector<Eigen::Vector3d>> contact_force_EE(contact_pairs.Edge_Edge.size());


		// compute the gradient and hessian
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

				contact_force_PG[i] = contactForce_node;
			}


			startIndex_grad += contact_pairs.Point_Ground.size() * 12, startIndex_hess += contact_pairs.Point_Ground.size() * 144;
#pragma omp parallel for num_threads(parameters.numOfThreads)
			for (int i = 0; i < contact_pairs.Point_Triangle_PP_index.size(); i++)
			{
				int ct_index = contact_pairs.Point_Triangle_PP_index[i];
				Vector5i PT_Contact = contact_pairs.Point_Triangle[ct_index];
				double dis2 = contact_pairs.Point_Triangle_PP_Dis[i];
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

				contact_force_PP[i] = contactForce_node;

			}


			startIndex_grad += contact_pairs.Point_Triangle_PP_index.size() * 24, startIndex_hess += contact_pairs.Point_Triangle_PP_index.size() * 144 * 4;
#pragma omp parallel for num_threads(parameters.numOfThreads)
			for (int i = 0; i < contact_pairs.Point_Triangle_PE_index.size(); i++)
			{
				int ct_index = contact_pairs.Point_Triangle_PE_index[i];
				Vector5i PT_Contact = contact_pairs.Point_Triangle[ct_index];
				double dis2 = contact_pairs.Point_Triangle_PE_Dis[i];
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

				contact_force_PE[i] = contactForce_node;
			}


			startIndex_grad += contact_pairs.Point_Triangle_PE_index.size() * 36, startIndex_hess += contact_pairs.Point_Triangle_PE_index.size() * 144 * 9;
#pragma omp parallel for num_threads(parameters.numOfThreads)
			for (int i = 0; i < contact_pairs.Point_Triangle_PT_index.size(); i++)
			{
				int ct_index = contact_pairs.Point_Triangle_PT_index[i];
				Vector5i PT_Contact = contact_pairs.Point_Triangle[ct_index];
				double dis2 = contact_pairs.Point_Triangle_PT_Dis[i];;
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

				contact_force_PT[i] = contactForce_node;
			}


			startIndex_grad += contact_pairs.Point_Triangle_PT_index.size() * 48, startIndex_hess += contact_pairs.Point_Triangle_PT_index.size() * 144 * 16;
#pragma omp parallel for num_threads(parameters.numOfThreads)
			for (int i = 0; i < contact_pairs.Edge_Edge.size(); i++)
			{
				Vector5i EE_Contact = contact_pairs.Edge_Edge[i];

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

				contact_force_EE[i] = contactForce_node;

			}

		}
		

		// compute the actual contact force
		// 1st int: object index; 2nd int: surface vertex index; Eigen::Vector3d: contact force
		std::map<int, std::map<int, Eigen::Vector3d>> contact_force_all;
		{
			for (int i = 0; i < contact_pairs.Point_Ground.size(); i++)
			{
				Vector2i PG_Contact = contact_pairs.Point_Ground[i];
				contact_force_all[PG_Contact[0]][PG_Contact[1]] += contact_force_PG[i][0];
			}


			for (int i = 0; i < contact_pairs.Point_Triangle_PP_index.size(); i++)
			{
				int ct_index = contact_pairs.Point_Triangle_PP_index[i];
				Vector5i PP_Contact = contact_pairs.Point_Triangle[ct_index];
				
				int obj1_index = PP_Contact[0], obj2_index = PP_Contact[2];
				int P_index = PP_Contact[1], tri_index = PP_Contact[3];					
				int T0_index = -99;
				Eigen::Vector3i tri = triSimMesh.allObjects[obj2_index].objectSurfaceMesh.faces[tri_index];
				if (PP_Contact[4] == 0)
				{
					T0_index = tri[0];
				}
				else if (PP_Contact[4] == 1)
				{
					T0_index = tri[1];
				}
				else if(PP_Contact[4] == 2)
				{
					T0_index = tri[2];
				}

				contact_force_all[obj1_index][P_index] += contact_force_PP[i][0];
				contact_force_all[obj2_index][T0_index] += contact_force_PP[i][1];

			}


			for (int i = 0; i < contact_pairs.Point_Triangle_PE_index.size(); i++)
			{
				int ct_index = contact_pairs.Point_Triangle_PE_index[i];
				Vector5i PE_Contact = contact_pairs.Point_Triangle[ct_index];
				
				int obj1_index = PE_Contact[0], obj2_index = PE_Contact[2];
				int P_index = PE_Contact[1], tri_index = PE_Contact[3];
				int T0_index = -99, T1_index = -99;
				Eigen::Vector3i tri = triSimMesh.allObjects[obj2_index].objectSurfaceMesh.faces[tri_index];
				if (PE_Contact[4] == 3)
				{
					T0_index = tri[0];
					T1_index = tri[1];
				}
				else if (PE_Contact[4] == 4)
				{
					T0_index = tri[1];
					T1_index = tri[2];
				}
				else if (PE_Contact[4] == 5)
				{
					T0_index = tri[2];
					T1_index = tri[0];
				}

				contact_force_all[obj1_index][P_index] += contact_force_PE[i][0];
				contact_force_all[obj2_index][T0_index] += contact_force_PE[i][1];
				contact_force_all[obj2_index][T1_index] += contact_force_PE[i][2];

			}


			for (int i = 0; i < contact_pairs.Point_Triangle_PT_index.size(); i++)
			{
				int ct_index = contact_pairs.Point_Triangle_PT_index[i];
				Vector5i PT_Contact = contact_pairs.Point_Triangle[ct_index];
				
				int obj1_index = PT_Contact[0], obj2_index = PT_Contact[2];
				int P_index = PT_Contact[1], tri_index = PT_Contact[3];				
				Eigen::Vector3i tri = triSimMesh.allObjects[obj2_index].objectSurfaceMesh.faces[tri_index];
				int T0_index = tri[0], T1_index = tri[1], T2_index = tri[2];

				contact_force_all[obj1_index][P_index] += contact_force_PT[i][0];
				contact_force_all[obj2_index][T0_index] += contact_force_PT[i][1];
				contact_force_all[obj2_index][T1_index] += contact_force_PT[i][2];
				contact_force_all[obj2_index][T2_index] += contact_force_PT[i][3];

			}


			for (int i = 0; i < contact_pairs.Edge_Edge.size(); i++)
			{
				Vector5i EE_Contact = contact_pairs.Edge_Edge[i];

				int obj1_index = EE_Contact[0], obj2_index = EE_Contact[2];
				int E1_index = EE_Contact[1], E2_index = EE_Contact[3];
				Eigen::Vector2i E1 = triSimMesh.allObjects[obj1_index].objectSurfaceMesh.edges[E1_index];
				Eigen::Vector2i E2 = triSimMesh.allObjects[obj2_index].objectSurfaceMesh.edges[E2_index];


				contact_force_all[obj1_index][E1[0]] += contact_force_EE[i][0];
				contact_force_all[obj1_index][E1[1]] += contact_force_EE[i][1];
				contact_force_all[obj2_index][E2[0]] += contact_force_EE[i][2];
				contact_force_all[obj2_index][E2[1]] += contact_force_EE[i][3];

			}

		}


		// check if any object will be broken
		if (iteration == 0 && timestep == 16)
		{
			if_start_fracture_sim(triSimMesh, broken_objects, contact_force_all);
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

	
	//std::cout << "			Newton solver iterations: " << solver.iterations() << "; error: " << solver.error() << std::endl;


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

	// update the actual vertex position on the surface
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int vI = 0; vI < triSimMesh.allObjects.size(); vI++)
	{
		// update each object's affine 
		triSimMesh.allObjects[vI].affine = triSimMesh.allObjects[vI].affine_last + step * ABD_direction[vI];

		// update each object's surface mesh
		for (int j = 0; j < triSimMesh.allObjects[vI].pos_node_surface.size(); j++)
		{
			Eigen::Matrix<double, 3, 12> Jx = build_Jx_matrix_for_ABD(triSimMesh.allObjects[vI].objectSurfaceMesh.vertices[j].transpose());
			triSimMesh.allObjects[vI].pos_node_surface[j] = triSimMesh.allObjects[vI].pos_node_surface_prev[j] +step * Jx * ABD_direction[vI];
		}

		// update each object's interior points
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

void if_start_fracture_sim(triMesh& triSimMesh, std::map<int, std::vector<Vector6d>>& broken_objects, std::map<int, std::map<int, Eigen::Vector3d>>& contact_force_all)
{
	broken_objects.clear();
	for (std::map<int, std::map<int, Eigen::Vector3d>>::iterator it = contact_force_all.begin(); it != contact_force_all.end(); it++)
	{
		int obj_index = it->first;
		std::map<int, Eigen::Vector3d> vert_force = it->second;

		if (triSimMesh.allObjects[obj_index].breakable == true)
		{
			double max_force = -999999999.0;
			std::vector<Vector6d> contactForce;

			for (std::map<int, Eigen::Vector3d>::iterator itf = vert_force.begin(); itf != vert_force.end(); itf++)
			{
				Vector6d contact_position_force;
				Eigen::Vector3d position_world = triSimMesh.allObjects[obj_index].pos_node_surface[itf->first];
				Eigen::Vector3d position_material = triSimMesh.allObjects[obj_index].objectSurfaceMesh.vertices[itf->first];
				Eigen::Vector3d force_world = itf->second;
				// the force is in the world space, convert it into material space
				Eigen::Vector3d t = triSimMesh.allObjects[obj_index].affine.segment<3>(0).transpose();
				Eigen::Matrix3d F = Eigen::Matrix3d::Zero();
				F.col(0) = triSimMesh.allObjects[obj_index].affine.segment<3>(3).transpose();
				F.col(1) = triSimMesh.allObjects[obj_index].affine.segment<3>(6).transpose();
				F.col(2) = triSimMesh.allObjects[obj_index].affine.segment<3>(9).transpose();
				Eigen::Vector3d force_material = F.inverse() * ((position_world + force_world) - t) - position_world;


				if (force_material.norm() > 0.0001) // to avoid numerical error
				{
					contact_position_force.block(0, 0, 3, 1) = position_material;
					contact_position_force.block(3, 0, 3, 1) = force_material * 100;
					contactForce.push_back(contact_position_force);
				}
				max_force = std::max(force_material.norm(), max_force);
			}
			std::cout << "			max_force = " << max_force << std::endl;
			if (max_force > triSimMesh.allObjects[obj_index].objectMaterial.fracture_start_force )
			{
				broken_objects[obj_index] = contactForce;
			}		
		}
	}

}

void fracture_sim(FEMParamters& parameters, triMesh& triSimMesh, std::map<int, std::vector<Vector6d>>& broken_objects, std::map<int, objMeshFormat>& crackSurface_object)
{
	crackSurface_object.clear(); // int: object's index

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

	for (auto it = crackSurface_object.begin(); it != crackSurface_object.end(); it++)
	{
		int objectIndex = it->first;
		objMeshFormat full_cut = it->second;
		float voxel_size = static_cast<float>(parameters.IPC_dis) * 3;
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

			ABD_Object child_object = parent_object;
			child_object.objectNote = parent_object.objectNote + "_F_" + std::to_string(cc);
			child_object.objectSurfaceMesh = childMesh;
			child_object.need_update = true;

			child_object.pos_node_surface = childMesh.vertices;
			child_object.pos_node_surface_prev = childMesh.vertices;
			child_object.pos_node_surface_direction.resize(childMesh.vertices.size());

			

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

