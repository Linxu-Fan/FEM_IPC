//
////
////Vector12d constructQVec(Eigen::Vector3d translation, Eigen::Matrix3d deformation) // construct q vector according translation and deformation
////{
////	Vector12d qb = Vector12d::Zero();
////	qb.block(0, 0, 3, 1) = translation;
////	qb.block(3, 0, 3, 1) = deformation.col(0);
////	qb.block(6, 0, 3, 1) = deformation.col(1);
////	qb.block(9, 0, 3, 1) = deformation.col(2);
////
////	return qb;
////}
////
////// compute the incremental potential energy
////double compute_IP_energy_ABD(Mesh_ABD& tetSimMesh, FEMParamters& parameters, int timestep)
////{
////	double energyVal = 0;	
////
////	// initeria energy contribution per affine body
////	std::vector<Vector12d> affine_force(tetSimMesh.num_meshes, Vector12d::Zero());
////	for (int i = 0; i < tetSimMesh.pos_node.size(); i++)
////	{
////		if (tetSimMesh.boundaryCondition_node[i].type == 2)
////		{
////			Eigen::Matrix<double, 3, 12> Jx;
////			Jx.block(0, 0, 3, 3) = Eigen::Matrix3d::Identity();
////			Jx.block(0, 3, 1, 3) = tetSimMesh.pos_node_Rest[i].transpose();
////			Jx.block(1, 6, 1, 3) = tetSimMesh.pos_node_Rest[i].transpose();
////			Jx.block(2, 9, 1, 3) = tetSimMesh.pos_node_Rest[i].transpose();
////
////			int AB_index = tetSimMesh.index_node[i][0];
////			Eigen::Vector3d ext_force = compute_external_force(tetSimMesh, i, timestep);
////			Vector12d affine_force_node = Jx.transpose() * ext_force;
////
////			affine_force[AB_index] += affine_force_node;		
////		}
////	}
////	std::vector<double> AB_ine_energy_vec(tetSimMesh.num_meshes,0);
////#pragma omp parallel for num_threads(parameters.numOfThreads)
////	for (int AI = 0; AI < AB_ine_energy_vec.size(); AI++)
////	{
////		// reconstruct q_prev vector
////		Vector12d qb_prev = constructQVec(tetSimMesh.translation_prev_ABD[AI], tetSimMesh.deformation_prev_ABD[AI]);
////		// reconstruct qb_vel vector
////		Vector12d qb_vel = constructQVec(tetSimMesh.translation_vel_ABD[AI], tetSimMesh.deformation_vel_ABD[AI]);	
////		// reconstruct q vector
////		Vector12d qb = constructQVec(tetSimMesh.translation_ABD[AI], tetSimMesh.deformation_ABD[AI]);
////
////		Vector12d qb_Minus_qbHat = qb - (qb_prev + parameters.dt * qb_vel + parameters.dt * parameters.dt * tetSimMesh.massMatrix_ABD[AI] * affine_force[AI]);
////		AB_ine_energy_vec[AI] = 0.5 * qb_Minus_qbHat.transpose() * tetSimMesh.massMatrix_ABD[AI] * qb_Minus_qbHat;
////	}	
////	energyVal += std::accumulate(AB_ine_energy_vec.begin(), AB_ine_energy_vec.end(), 0.0);
////
////
////
////	// affine energy contribution per affine body
////	std::vector<double> AB_aff_energy_vec(tetSimMesh.num_meshes,0);
////#pragma omp parallel for num_threads(parameters.numOfThreads)
////	for (int AI = 0; AI < AB_aff_energy_vec.size(); AI++)
////	{
////		// original implementation
////		//AB_aff_energy_vec[AI] = parameters.ABD_Coeff * tetSimMesh.volume_ABD[AI] * (tetSimMesh.deformation_ABD[AI].transpose() * tetSimMesh.deformation_ABD[AI] - Eigen::Matrix3d::Identity()).squaredNorm();
////
////		// neo-hookean implementation
////		Material mat;
////		mat.mu = parameters.ABD_Coeff;
////		AB_aff_energy_vec[AI] = ElasticEnergy::Val(mat, "ARAP_linear", tetSimMesh.deformation_ABD[AI], parameters.dt, tetSimMesh.volume_ABD[AI]);
////
////	}
////	energyVal += std::accumulate(AB_aff_energy_vec.begin(), AB_aff_energy_vec.end(), 0.0);
////
////
////
////	// energy contribution from barrier
////	energyVal += compute_Barrier_energy(tetSimMesh, parameters, timestep);
////
////
////	return energyVal;
////}
////
////// compute energy gradient and hessian of the linear system and solve the syetem
////std::vector<Eigen::Vector3d> solve_linear_system_ABD(Mesh_ABD& tetSimMesh, FEMParamters& parameters, int timestep)
////{
////	// we only need to update boundary vertices' position to calculate the barrier energy update
////	std::vector<Eigen::Vector3d> movingDir(tetSimMesh.boundaryVertices.size());
////
////
////	// initialize grad and hessian
////	Eigen::SparseMatrix<double> leftHandSide(12 * tetSimMesh.num_meshes, 12 * tetSimMesh.num_meshes);
////	leftHandSide.setZero();
////	Eigen::VectorXd rightHandSide = Eigen::VectorXd(12 * tetSimMesh.num_meshes);
////	rightHandSide.setZero();
////
////
////	// contribution from inertia term 
////	std::vector<Vector12d> affine_force(tetSimMesh.num_meshes, Vector12d::Zero());
////	for (int i = 0; i < tetSimMesh.pos_node.size(); i++)
////	{
////		if (tetSimMesh.boundaryCondition_node[i].type == 2)
////		{
////			Eigen::Matrix<double, 3, 12> Jx;
////			Jx.block(0, 0, 3, 3) = Eigen::Matrix3d::Identity();
////			Jx.block(0, 3, 1, 3) = tetSimMesh.pos_node_Rest[i].transpose();
////			Jx.block(1, 6, 1, 3) = tetSimMesh.pos_node_Rest[i].transpose();
////			Jx.block(2, 9, 1, 3) = tetSimMesh.pos_node_Rest[i].transpose();
////
////			int AB_index = tetSimMesh.index_node[i][0];
////			Eigen::Vector3d ext_force = compute_external_force(tetSimMesh, i, timestep);
////			Vector12d affine_force_node = Jx.transpose() * ext_force;
////
////			affine_force[AB_index] += affine_force_node;
////		}
////	}
////	std::vector<double> AB_ine_energy_vec(tetSimMesh.num_meshes, 0);
////#pragma omp parallel for num_threads(parameters.numOfThreads)
////	for (int AI = 0; AI < AB_ine_energy_vec.size(); AI++)
////	{
////		// reconstruct q_prev vector
////		Vector12d qb_prev = constructQVec(tetSimMesh.translation_prev_ABD[AI], tetSimMesh.deformation_prev_ABD[AI]);
////		// reconstruct qb_vel vector
////		Vector12d qb_vel = constructQVec(tetSimMesh.translation_vel_ABD[AI], tetSimMesh.deformation_vel_ABD[AI]);
////		// reconstruct q vector
////		Vector12d qb = constructQVec(tetSimMesh.translation_ABD[AI], tetSimMesh.deformation_ABD[AI]);
////
////		Vector12d qb_Minus_qbHat = qb - (qb_prev + parameters.dt * qb_vel + parameters.dt * parameters.dt * tetSimMesh.massMatrix_ABD[AI] * affine_force[AI]);
////		rightHandSide.block(12 * AI, 0, 12, 1) = tetSimMesh.massMatrix_ABD[AI] * qb_Minus_qbHat;
////		leftHandSide.block(12 * AI, 12 * AI, 12, 12) = tetSimMesh.massMatrix_ABD[AI];
////	}
////
////
////	// contribution from affine energy
////	for (int AI = 0; AI < tetSimMesh.num_meshes; AI++)
////	{
////		
////		// calculate the grad and hessian for each column
////		for (int i = 0; i < 3; i++)
////		{
////			Eigen::Vector3d ai = tetSimMesh.deformation_ABD[AI].col(i);
////			Eigen::Matrix3d ajaj = Eigen::Matrix3d::Zero();
////			for (int j = 0; j < 3; j++)
////			{
////				if (j != i)
////				{
////					Eigen::Vector3d aj = tetSimMesh.deformation_ABD[AI].col(j);
////					ajaj += aj.cross(aj);
////				}			
////			}
////			
////			Eigen::Vector3d ai_grad = 2 * parameters.ABD_Coeff * tetSimMesh.volume_ABD[AI] * (2 * (ai.dot(ai) - 1) * ai + 2 * ajaj * ai);
////			Eigen::Matrix3d aii_hess = 2 * parameters.ABD_Coeff * tetSimMesh.volume_ABD[AI] * (4 * ai.cross(ai) + 2 * (ai.squaredNorm() - 1) * Eigen::Matrix3d::Identity() + 2 * ajaj);
////
////			rightHandSide.block(12 * AI + 3 + 3 * i, 0, 3, 1) = ai_grad;
////			leftHandSide.block(12 * AI + 3 + 3 * i, 12 * AI + 3 + 3 * i, 3, 3) = aii_hess;
////
////			for (int j = 0; j < 3; j++)
////			{
////				if (j != i)
////				{
////					Eigen::Vector3d aj = tetSimMesh.deformation_ABD[AI].col(j);
////
////					Eigen::Matrix3d aij_hess = 2 * parameters.ABD_Coeff * tetSimMesh.volume_ABD[AI] * ( 2 * (aj * ai.transpose() + ai * aj.transpose()));
////					leftHandSide.block(12 * AI + 3 + 3 * i, 12 * AI + 3 + 3 * j, 3, 3) = aij_hess;
////				}
////			}
////
////		}
////
////
////	}
////
////
////	// contribution from contact
////	//		calculate contacts
////	std::vector<Vector5i> PG_PG, PT_PP, PT_PE, PT_PT, EE_EE;
////	calContactInfo(tetSimMesh, parameters, timestep, PG_PG, PT_PP, PT_PE, PT_PT, EE_EE);
////
////
////
////
////
////	int gradSize = tetSimMesh.pos_node.size() * 6 + tetSimMesh.tetrahedrals.size() * 12 + PG_PG.size() * 3
////		+ PT_PP.size() * 6 + PT_PE.size() * 9 + PT_PT.size() * 12 + EE_EE.size() * 12;
////	int hessSize = tetSimMesh.pos_node.size() * 3 + tetSimMesh.tetrahedrals.size() * 144 + PG_PG.size() * 9
////		+ PT_PP.size() * 36 + PT_PE.size() * 81 + PT_PT.size() * 144 + EE_EE.size() * 144;
////	std::vector<std::pair<int, double>> grad_triplet(gradSize);
////	std::vector<Eigen::Triplet<double>> hessian_triplet(hessSize);
////
////
////	// energy contribution per vertex
////#pragma omp parallel for num_threads(parameters.numOfThreads)
////	for (int vI = 0; vI < tetSimMesh.pos_node.size(); vI++)
////	{
////		double nodeMass = tetSimMesh.mass_node[vI];
////		Eigen::Vector3d xt = tetSimMesh.pos_node_prev[vI];
////		Eigen::Vector3d v = tetSimMesh.vel_node[vI];
////		Eigen::Vector3d x = tetSimMesh.pos_node[vI];
////		Eigen::Vector3d extForce = compute_external_force(tetSimMesh, vI, timestep);
////
////		// !!!! The order is very important. "InertiaEnergy" must be the first and "ExternalEnergy" the second
////
////		// the inertia energy contribution
////		int actGradIndex = vI * 6;
////		InertiaEnergy::Grad(grad_triplet, actGradIndex, nodeMass, parameters.dt, xt, v, x, extForce, vI, parameters, tetSimMesh.boundaryCondition_node, timestep);
////		int actHessIndex = vI * 3;
////		InertiaEnergy::Hess(hessian_triplet, actHessIndex, nodeMass, vI, tetSimMesh.boundaryCondition_node, timestep, parameters);
////
////		// the external energy contribution
////		actGradIndex = vI * 6 + 3;
////		ExternalEnergy::Grad(grad_triplet, actGradIndex, nodeMass, parameters.dt, x, parameters, extForce, vI, tetSimMesh.boundaryCondition_node[vI].type);
////
////	}
////
////
////	int startIndex_grad = tetSimMesh.pos_node.size() * 6, startIndex_hess = tetSimMesh.pos_node.size() * 3;
////	// energy contribution per element
////#pragma omp parallel for num_threads(parameters.numOfThreads)
////	for (int eI = 0; eI < tetSimMesh.tetrahedrals.size(); eI++)
////	{
////
////		Eigen::Vector4i tetVertInd_BC;
////		for (int hg = 0; hg < 4; hg++)
////		{
////			int vertInd = tetSimMesh.tetrahedrals[eI][hg];
////			tetVertInd_BC[hg] = tetSimMesh.boundaryCondition_node[vertInd].type;
////		}
////
////		// the internal elastic energy contribution
////		int matInd = tetSimMesh.materialInd[eI];
////		Eigen::Matrix<double, 9, 12> dFdx = ElasticEnergy::dF_wrt_dx(tetSimMesh.tetra_DM_inv[eI]);
////		int actualStartIndex_hess = startIndex_hess + eI * 144, actualStartIndex_grad = startIndex_grad + eI * 12;
////		ElasticEnergy::Grad(grad_triplet, actualStartIndex_grad, tetSimMesh.materialMesh[matInd], parameters.model, tetSimMesh.tetra_F[eI], parameters.dt, tetSimMesh.tetra_vol[eI], dFdx, tetSimMesh.tetrahedrals[eI], tetVertInd_BC);
////		ElasticEnergy::Hess(hessian_triplet, actualStartIndex_hess, tetSimMesh.materialMesh[matInd], parameters.model, tetSimMesh.tetra_F[eI], parameters.dt, tetSimMesh.tetra_vol[eI], dFdx, tetSimMesh.tetrahedrals[eI], tetVertInd_BC);
////
////
////	}
////
////
////	//double endTime2 = omp_get_wtime();
////	//std::cout << "	Element grad_hess Time is : " << endTime2 - endTime1 << "s" << std::endl;
////
////
////	
////
////
////	// calculate the barrier gradient and hessian
////	{
////		startIndex_grad += tetSimMesh.tetrahedrals.size() * 12,
////			startIndex_hess += tetSimMesh.tetrahedrals.size() * 144;
////#pragma omp parallel for num_threads(parameters.numOfThreads)
////		for (int i = 0; i < PG_PG.size(); i++)
////		{
////			int ptInd = PG_PG[i][1];
////			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
////			int actualStartIndex_hess = startIndex_hess + i * 9,
////				actualStartIndex_grad = startIndex_grad + i * 3;
////			double z2 = P[2] * P[2];
////			Ground::gradAndHess(hessian_triplet, grad_triplet, actualStartIndex_hess,
////				actualStartIndex_grad, ptInd,
////				z2, tetSimMesh.boundaryVertices_area[ptInd], parameters);
////		}
////
////		startIndex_grad += PG_PG.size() * 3, startIndex_hess += PG_PG.size() * 9;
////#pragma omp parallel for num_threads(parameters.numOfThreads)
////		for (int i = 0; i < PT_PP.size(); i++)
////		{
////			Vector5i cont_PT = PT_PP[i];
////
////			int ptInd = cont_PT[1];
////			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
////			Eigen::Vector3i tri = tetSimMesh.boundaryTriangles[cont_PT[2]];
////			Eigen::Vector3d A = tetSimMesh.pos_node[tri[0]];
////			Eigen::Vector3d B = tetSimMesh.pos_node[tri[1]];
////			Eigen::Vector3d C = tetSimMesh.pos_node[tri[2]];
////			int type = cont_PT[3];
////
////			double dis2 = 0;
////			DIS::computePointTriD(P, A, B, C, dis2);
////			Eigen::Vector4i ptIndices = { ptInd , tri[0] , tri[1] , tri[2] };
////			int actualStartIndex_hess = startIndex_hess + i * 36,
////				actualStartIndex_grad = startIndex_grad + i * 6;
////			BarrierEnergy::gradAndHess_PT(hessian_triplet, grad_triplet, actualStartIndex_hess,
////				actualStartIndex_grad, ptIndices, type,
////				dis2, tetSimMesh, parameters);
////
////		}
////
////		startIndex_grad += PT_PP.size() * 6, startIndex_hess += PT_PP.size() * 36;
////#pragma omp parallel for num_threads(parameters.numOfThreads)
////		for (int i = 0; i < PT_PE.size(); i++)
////		{
////			Vector5i cont_PT = PT_PE[i];
////
////			int ptInd = cont_PT[1];
////			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
////			Eigen::Vector3i tri = tetSimMesh.boundaryTriangles[cont_PT[2]];
////			Eigen::Vector3d A = tetSimMesh.pos_node[tri[0]];
////			Eigen::Vector3d B = tetSimMesh.pos_node[tri[1]];
////			Eigen::Vector3d C = tetSimMesh.pos_node[tri[2]];
////			int type = cont_PT[3];
////
////			double dis2 = 0;
////			DIS::computePointTriD(P, A, B, C, dis2);
////			Eigen::Vector4i ptIndices = { ptInd , tri[0] , tri[1] , tri[2] };
////			int actualStartIndex_hess = startIndex_hess + i * 81,
////				actualStartIndex_grad = startIndex_grad + i * 9;
////			BarrierEnergy::gradAndHess_PT(hessian_triplet, grad_triplet, actualStartIndex_hess,
////				actualStartIndex_grad, ptIndices, type,
////				dis2, tetSimMesh, parameters);
////
////		}
////
////		startIndex_grad += PT_PE.size() * 9, startIndex_hess += PT_PE.size() * 81;
////#pragma omp parallel for num_threads(parameters.numOfThreads)
////		for (int i = 0; i < PT_PT.size(); i++)
////		{
////			Vector5i cont_PT = PT_PT[i];
////
////			int ptInd = cont_PT[1];
////			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
////			Eigen::Vector3i tri = tetSimMesh.boundaryTriangles[cont_PT[2]];
////			Eigen::Vector3d A = tetSimMesh.pos_node[tri[0]];
////			Eigen::Vector3d B = tetSimMesh.pos_node[tri[1]];
////			Eigen::Vector3d C = tetSimMesh.pos_node[tri[2]];
////			int type = cont_PT[3];
////
////			double dis2 = 0;
////			DIS::computePointTriD(P, A, B, C, dis2);
////			Eigen::Vector4i ptIndices = { ptInd , tri[0] , tri[1] , tri[2] };
////			int actualStartIndex_hess = startIndex_hess + i * 144,
////				actualStartIndex_grad = startIndex_grad + i * 12;
////			BarrierEnergy::gradAndHess_PT(hessian_triplet, grad_triplet, actualStartIndex_hess,
////				actualStartIndex_grad, ptIndices, type,
////				dis2, tetSimMesh, parameters);
////
////		}
////
////		startIndex_grad += PT_PT.size() * 12, startIndex_hess += PT_PT.size() * 144;
////#pragma omp parallel for num_threads(parameters.numOfThreads)
////		for (int i = 0; i < EE_EE.size(); i++)
////		{
////			Vector5i cont_EE = EE_EE[i];
////			int E1 = cont_EE[1], E2 = cont_EE[2];
////			int e1p1 = tetSimMesh.index_boundaryEdge[E1][0], e1p2 = tetSimMesh.index_boundaryEdge[E1][1],
////				e2p1 = tetSimMesh.index_boundaryEdge[E2][0], e2p2 = tetSimMesh.index_boundaryEdge[E2][1];
////
////			Eigen::Vector3d P1 = tetSimMesh.pos_node[e1p1];
////			Eigen::Vector3d P2 = tetSimMesh.pos_node[e1p2];
////			Eigen::Vector3d Q1 = tetSimMesh.pos_node[e2p1];
////			Eigen::Vector3d Q2 = tetSimMesh.pos_node[e2p2];
////			int type = cont_EE[3];
////
////			double dis2 = 0;
////			DIS::computeEdgeEdgeD(P1, P2, Q1, Q2, dis2);
////			Eigen::Vector4i ptIndices = { e1p1 , e1p2 , e2p1 , e2p2 };
////			int actualStartIndex_hess = startIndex_hess + i * 144,
////				actualStartIndex_grad = startIndex_grad + i * 12;
////			BarrierEnergy::gradAndHess_EE(hessian_triplet, grad_triplet, actualStartIndex_hess,
////				actualStartIndex_grad, ptIndices, type,
////				dis2, tetSimMesh, parameters);
////
////		}
////
////	}
////
////
////
////
////	//double endTime3 = omp_get_wtime();
////	//std::cout << "	Cal grad_hess Time is : " << endTime3 - endTime2 << "s" << std::endl;
////
////	// assemable the left-hand side 
////	Eigen::SparseMatrix<double> leftHandSide(3 * tetSimMesh.pos_node.size(), 3 * tetSimMesh.pos_node.size());
////	leftHandSide.setFromTriplets(hessian_triplet.begin(), hessian_triplet.end());
////
////	// assemable the right-hand side 
////	Eigen::VectorXd rightHandSide = Eigen::VectorXd(tetSimMesh.pos_node.size() * 3);
////	rightHandSide.setZero();
////	for (int i = 0; i < grad_triplet.size(); i++)
////	{
////		std::pair<int, double> ele = grad_triplet[i];
////		rightHandSide[ele.first] += ele.second;
////	}
////
////	//double endTime4 = omp_get_wtime();
////	//std::cout << "	Assemble grad_hess Time is : " << endTime4 - endTime3 << "s" << std::endl;
////
////	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
////	solver.compute(leftHandSide);
////	Eigen::VectorXd result = solver.solve(-rightHandSide);
////
////
////
////#pragma omp parallel for num_threads(parameters.numOfThreads)
////	for (int i = 0; i < tetSimMesh.pos_node.size(); i++)
////	{
////		movingDir[i] = result.block<3, 1>(3 * i, 0);
////	}
////
////
////
////	return movingDir;
////}
////
//
//
//
//
//
//
//
//
//#include "simulator.h"
//
//
//Vector12d constructQVec(Eigen::Vector3d translation, Eigen::Matrix3d deformation) // construct q vector according translation and deformation
//{
//	Vector12d qb = Vector12d::Zero();
//	qb.block(0, 0, 3, 1) = translation;
//	qb.block(3, 0, 3, 1) = deformation.col(0);
//	qb.block(6, 0, 3, 1) = deformation.col(1);
//	qb.block(9, 0, 3, 1) = deformation.col(2);
//
//	return qb;
//}
//
//// compute the incremental potential energy
//double compute_IP_energy_ABD(Mesh_ABD& tetSimMesh, FEMParamters& parameters, int timestep)
//{
//	double energyVal = 0;	
//
//	// initeria energy contribution per affine body
//	std::vector<Vector12d> affine_force(tetSimMesh.num_meshes, Vector12d::Zero());
//	for (int i = 0; i < tetSimMesh.pos_node.size(); i++)
//	{
//		if (tetSimMesh.boundaryCondition_node[i].type == 2)
//		{
//			Eigen::Matrix<double, 3, 12> Jx;
//			Jx.block(0, 0, 3, 3) = Eigen::Matrix3d::Identity();
//			Jx.block(0, 3, 1, 3) = tetSimMesh.pos_node_Rest[i].transpose();
//			Jx.block(1, 6, 1, 3) = tetSimMesh.pos_node_Rest[i].transpose();
//			Jx.block(2, 9, 1, 3) = tetSimMesh.pos_node_Rest[i].transpose();
//
//			int AB_index = tetSimMesh.index_node[i][0];
//			Eigen::Vector3d ext_force = compute_external_force(tetSimMesh, i, timestep);
//			Vector12d affine_force_node = Jx.transpose() * ext_force;
//
//			affine_force[AB_index] += affine_force_node;		
//		}
//	}
//	std::vector<double> AB_ine_energy_vec(tetSimMesh.num_meshes,0);
//#pragma omp parallel for num_threads(parameters.numOfThreads)
//	for (int AI = 0; AI < AB_ine_energy_vec.size(); AI++)
//	{
//		// reconstruct q_prev vector
//		Vector12d qb_prev = constructQVec(tetSimMesh.translation_prev_ABD[AI], tetSimMesh.deformation_prev_ABD[AI]);
//		// reconstruct qb_vel vector
//		Vector12d qb_vel = constructQVec(tetSimMesh.translation_vel_ABD[AI], tetSimMesh.deformation_vel_ABD[AI]);	
//		// reconstruct q vector
//		Vector12d qb = constructQVec(tetSimMesh.translation_ABD[AI], tetSimMesh.deformation_ABD[AI]);
//
//		Vector12d qb_Minus_qbHat = qb - (qb_prev + parameters.dt * qb_vel + parameters.dt * parameters.dt * tetSimMesh.massMatrix_ABD[AI].inverse() * affine_force[AI]);
//		AB_ine_energy_vec[AI] = 0.5 * qb_Minus_qbHat.transpose() * tetSimMesh.massMatrix_ABD[AI] * qb_Minus_qbHat;
//	}	
//	energyVal += std::accumulate(AB_ine_energy_vec.begin(), AB_ine_energy_vec.end(), 0.0);
//
//
//
//	// affine energy contribution per affine body
//	std::vector<double> AB_aff_energy_vec(tetSimMesh.num_meshes,0);
//#pragma omp parallel for num_threads(parameters.numOfThreads)
//	for (int AI = 0; AI < AB_aff_energy_vec.size(); AI++)
//	{
//		// original implementation
//		//AB_aff_energy_vec[AI] = parameters.ABD_Coeff * tetSimMesh.volume_ABD[AI] * (tetSimMesh.deformation_ABD[AI].transpose() * tetSimMesh.deformation_ABD[AI] - Eigen::Matrix3d::Identity()).squaredNorm();
//
//		// neo-hookean implementation
//		double J = tetSimMesh.deformation_ABD[AI].determinant();
//		AB_aff_energy_vec[AI] = parameters.dt * parameters.dt * tetSimMesh.volume_ABD[AI] * parameters.ABD_Coeff / 2.0 * (0.5 * (tetSimMesh.deformation_ABD[AI]
//								+ tetSimMesh.deformation_ABD[AI].transpose()) - Eigen::Matrix3d::Identity()).squaredNorm();
//
//	}
//	energyVal += std::accumulate(AB_aff_energy_vec.begin(), AB_aff_energy_vec.end(), 0.0);
//
//
//
//	// energy contribution from barrier
//	energyVal += compute_Barrier_energy(tetSimMesh, parameters, timestep);
//
//
//	return energyVal;
//}
//
//// compute energy gradient and hessian of the linear system and solve the syetem
//std::vector<Eigen::Vector3d> solve_linear_system_ABD(Mesh_ABD& tetSimMesh, FEMParamters& parameters, int timestep)
//{
//	// we only need to update boundary vertices' position to calculate the barrier energy update
//	std::vector<Eigen::Vector3d> movingDir(tetSimMesh.boundaryVertices.size());
//
//
//
//	// contribution from contact
//	//		calculate contacts
//	std::vector<Vector5i> PG_PG, PT_PP, PT_PE, PT_PT, EE_EE;
//	calContactInfo(tetSimMesh, parameters, timestep, PG_PG, PT_PP, PT_PE, PT_PT, EE_EE);
//
//
//	int gradSize = tetSimMesh.num_meshes * 12 + tetSimMesh.num_meshes * 9 + PG_PG.size() * 12
//		+ PT_PP.size() * 24 + PT_PE.size() * 36 + PT_PT.size() * 48 + EE_EE.size() * 48;
//	int hessSize = tetSimMesh.num_meshes * 144 + tetSimMesh.num_meshes * 81 + PG_PG.size() * 144
//		+ PT_PP.size() * 144 * 4 + PT_PE.size() * 144 * 9 + PT_PT.size() * 144 * 16 + EE_EE.size() * 144 * 16;
//	std::vector<std::pair<int, double>> grad_triplet(gradSize, std::make_pair(0,0.0));
//	std::vector<Eigen::Triplet<double>> hessian_triplet(hessSize, Eigen::Triplet<double>(0, 0, 0.0));
//
//
//
//
//	// contribution from inertia term 
//	std::vector<Vector12d> affine_force(tetSimMesh.num_meshes, Vector12d::Zero());
//	for (int i = 0; i < tetSimMesh.pos_node.size(); i++)
//	{
//		if (tetSimMesh.boundaryCondition_node[i].type == 2)
//		{
//			Eigen::Matrix<double, 3, 12> Jx;
//			Jx.block(0, 0, 3, 3) = Eigen::Matrix3d::Identity();
//			Jx.block(0, 3, 1, 3) = tetSimMesh.pos_node_Rest[i].transpose();
//			Jx.block(1, 6, 1, 3) = tetSimMesh.pos_node_Rest[i].transpose();
//			Jx.block(2, 9, 1, 3) = tetSimMesh.pos_node_Rest[i].transpose();
//
//			int AB_index = tetSimMesh.index_node[i][0];
//			Eigen::Vector3d ext_force = compute_external_force(tetSimMesh, i, timestep);
//			Vector12d affine_force_node = Jx.transpose() * ext_force;
//
//			affine_force[AB_index] += affine_force_node;
//		}
//	}
//	
//	
//	int startIndex_grad = 0, startIndex_hess = 0;
//	std::vector<double> AB_ine_energy_vec(tetSimMesh.num_meshes, 0);
//#pragma omp parallel for num_threads(parameters.numOfThreads)
//	for (int AI = 0; AI < AB_ine_energy_vec.size(); AI++)
//	{
//		// reconstruct q_prev vector
//		Vector12d qb_prev = constructQVec(tetSimMesh.translation_prev_ABD[AI], tetSimMesh.deformation_prev_ABD[AI]);
//		// reconstruct qb_vel vector
//		Vector12d qb_vel = constructQVec(tetSimMesh.translation_vel_ABD[AI], tetSimMesh.deformation_vel_ABD[AI]);
//		// reconstruct q vector
//		Vector12d qb = constructQVec(tetSimMesh.translation_ABD[AI], tetSimMesh.deformation_ABD[AI]);
//
//		Vector12d qb_Minus_qbHat = qb - (qb_prev + parameters.dt * qb_vel + parameters.dt * parameters.dt * tetSimMesh.massMatrix_ABD[AI].inverse() * affine_force[AI]);
//
//		Vector12d grad = tetSimMesh.massMatrix_ABD[AI] * qb_Minus_qbHat;
//		Matrix12d hess = tetSimMesh.massMatrix_ABD[AI];
//
//		for (int i = 0; i < 12; i++)
//		{
//			grad_triplet[startIndex_grad + AI * 12 + i] = { AI * 12 + i, grad[i] };
//			for (int j = 0; j < 12; j++)
//			{
//				hessian_triplet[startIndex_grad + AI * 144 + i * 12 + j] = { 12 * AI + i, 12 * AI + j, hess(i,j)};
//			}
//		}
//	}
//
//
//	// contribution from affine energy
//	startIndex_grad = tetSimMesh.num_meshes * 12 , startIndex_hess = tetSimMesh.num_meshes * 144;
//#pragma omp parallel for num_threads(parameters.numOfThreads)
//	for (int AI = 0; AI < tetSimMesh.num_meshes; AI++)
//	{
//		
//		//// original implementation
//		//for (int i = 0; i < 3; i++)
//		//{
//		//	Eigen::Vector3d ai = tetSimMesh.deformation_ABD[AI].col(i);
//		//	Eigen::Matrix3d ajaj = Eigen::Matrix3d::Zero();
//		//	for (int j = 0; j < 3; j++)
//		//	{
//		//		if (j != i)
//		//		{
//		//			Eigen::Vector3d aj = tetSimMesh.deformation_ABD[AI].col(j);
//		//			ajaj += aj.cross(aj);
//		//		}			
//		//	}
//		//	
//		//	Eigen::Vector3d ai_grad = 2 * parameters.ABD_Coeff * tetSimMesh.volume_ABD[AI] * (2 * (ai.dot(ai) - 1) * ai + 2 * ajaj * ai);
//		//	Eigen::Matrix3d aii_hess = 2 * parameters.ABD_Coeff * tetSimMesh.volume_ABD[AI] * (4 * ai.cross(ai) + 2 * (ai.squaredNorm() - 1) * Eigen::Matrix3d::Identity() + 2 * ajaj);
//
//		//	rightHandSide.block(12 * AI + 3 + 3 * i, 0, 3, 1) = ai_grad;
//		//	leftHandSide.block(12 * AI + 3 + 3 * i, 12 * AI + 3 + 3 * i, 3, 3) = aii_hess;
//
//		//	for (int j = 0; j < 3; j++)
//		//	{
//		//		if (j != i)
//		//		{
//		//			Eigen::Vector3d aj = tetSimMesh.deformation_ABD[AI].col(j);
//
//		//			Eigen::Matrix3d aij_hess = 2 * parameters.ABD_Coeff * tetSimMesh.volume_ABD[AI] * ( 2 * (aj * ai.transpose() + ai * aj.transpose()));
//		//			leftHandSide.block(12 * AI + 3 + 3 * i, 12 * AI + 3 + 3 * j, 3, 3) = aij_hess;
//		//		}
//		//	}
//
//		//}
//
//
//		// neo-hookean implementation
//		Eigen::Matrix3d PK1 = parameters.dt * parameters.dt * tetSimMesh.volume_ABD[AI] * parameters.ABD_Coeff / 2 * (tetSimMesh.deformation_ABD[AI] 
//							  + tetSimMesh.deformation_ABD[AI].transpose() - 2 * Eigen::Matrix3d::Identity());;
//		Vector9d grad_E_wrt_F = flatenMatrix3d(PK1);
//		Eigen::Matrix<double, 9, 9> hess_E_wrt_F = parameters.dt * parameters.dt * tetSimMesh.volume_ABD[AI] * parameters.ABD_Coeff * Eigen::Matrix<double, 9, 9>::Identity();
//
//		for (int i = 0; i < 9; i++)
//		{
//			grad_triplet[startIndex_grad + AI * 9 + i] = { AI * 12 + 3 + i, grad_E_wrt_F[i] };
//			for (int j = 0; j < 9; j++)
//			{
//				hessian_triplet[startIndex_grad + AI * 81 + i * 9 + j] = { 12 * AI + 3 + i, 12 * AI + 3 + j, hess_E_wrt_F(i,j) };
//			}
//		}
//
//	}
//
//
//
//
//
//
//	// calculate the barrier gradient and hessian
//	{
//		startIndex_grad += tetSimMesh.num_meshes * 9,
//			startIndex_hess += tetSimMesh.num_meshes * 81;
//#pragma omp parallel for num_threads(parameters.numOfThreads)
//		for (int i = 0; i < PG_PG.size(); i++)
//		{
//			int ptInd = PG_PG[i][1];
//			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
//			int actualStartIndex_hess = startIndex_hess + i * 9,
//				actualStartIndex_grad = startIndex_grad + i * 3;
//			double z2 = P[2] * P[2];
//			Ground::gradAndHess(hessian_triplet, grad_triplet, actualStartIndex_hess,
//				actualStartIndex_grad, ptInd,
//				z2, tetSimMesh.boundaryVertices_area[ptInd], parameters);
//		}
//
//		startIndex_grad += PG_PG.size() * 3, startIndex_hess += PG_PG.size() * 9;
//#pragma omp parallel for num_threads(parameters.numOfThreads)
//		for (int i = 0; i < PT_PP.size(); i++)
//		{
//			Vector5i cont_PT = PT_PP[i];
//
//			int ptInd = cont_PT[1];
//			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
//			Eigen::Vector3i tri = tetSimMesh.boundaryTriangles[cont_PT[2]];
//			Eigen::Vector3d A = tetSimMesh.pos_node[tri[0]];
//			Eigen::Vector3d B = tetSimMesh.pos_node[tri[1]];
//			Eigen::Vector3d C = tetSimMesh.pos_node[tri[2]];
//			int type = cont_PT[3];
//
//			double dis2 = 0;
//			DIS::computePointTriD(P, A, B, C, dis2);
//			Eigen::Vector4i ptIndices = { ptInd , tri[0] , tri[1] , tri[2] };
//			int actualStartIndex_hess = startIndex_hess + i * 36,
//				actualStartIndex_grad = startIndex_grad + i * 6;
//			BarrierEnergy::gradAndHess_PT(hessian_triplet, grad_triplet, actualStartIndex_hess,
//				actualStartIndex_grad, ptIndices, type,
//				dis2, tetSimMesh, parameters);
//
//		}
//
//		startIndex_grad += PT_PP.size() * 6, startIndex_hess += PT_PP.size() * 36;
//#pragma omp parallel for num_threads(parameters.numOfThreads)
//		for (int i = 0; i < PT_PE.size(); i++)
//		{
//			Vector5i cont_PT = PT_PE[i];
//
//			int ptInd = cont_PT[1];
//			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
//			Eigen::Vector3i tri = tetSimMesh.boundaryTriangles[cont_PT[2]];
//			Eigen::Vector3d A = tetSimMesh.pos_node[tri[0]];
//			Eigen::Vector3d B = tetSimMesh.pos_node[tri[1]];
//			Eigen::Vector3d C = tetSimMesh.pos_node[tri[2]];
//			int type = cont_PT[3];
//
//			double dis2 = 0;
//			DIS::computePointTriD(P, A, B, C, dis2);
//			Eigen::Vector4i ptIndices = { ptInd , tri[0] , tri[1] , tri[2] };
//			int actualStartIndex_hess = startIndex_hess + i * 81,
//				actualStartIndex_grad = startIndex_grad + i * 9;
//			BarrierEnergy::gradAndHess_PT(hessian_triplet, grad_triplet, actualStartIndex_hess,
//				actualStartIndex_grad, ptIndices, type,
//				dis2, tetSimMesh, parameters);
//
//		}
//
//		startIndex_grad += PT_PE.size() * 9, startIndex_hess += PT_PE.size() * 81;
//#pragma omp parallel for num_threads(parameters.numOfThreads)
//		for (int i = 0; i < PT_PT.size(); i++)
//		{
//			Vector5i cont_PT = PT_PT[i];
//
//			int ptInd = cont_PT[1];
//			Eigen::Vector3d P = tetSimMesh.pos_node[ptInd];
//			Eigen::Vector3i tri = tetSimMesh.boundaryTriangles[cont_PT[2]];
//			Eigen::Vector3d A = tetSimMesh.pos_node[tri[0]];
//			Eigen::Vector3d B = tetSimMesh.pos_node[tri[1]];
//			Eigen::Vector3d C = tetSimMesh.pos_node[tri[2]];
//			int type = cont_PT[3];
//
//			double dis2 = 0;
//			DIS::computePointTriD(P, A, B, C, dis2);
//			Eigen::Vector4i ptIndices = { ptInd , tri[0] , tri[1] , tri[2] };
//			int actualStartIndex_hess = startIndex_hess + i * 144,
//				actualStartIndex_grad = startIndex_grad + i * 12;
//			BarrierEnergy::gradAndHess_PT(hessian_triplet, grad_triplet, actualStartIndex_hess,
//				actualStartIndex_grad, ptIndices, type,
//				dis2, tetSimMesh, parameters);
//
//		}
//
//		startIndex_grad += PT_PT.size() * 12, startIndex_hess += PT_PT.size() * 144;
//#pragma omp parallel for num_threads(parameters.numOfThreads)
//		for (int i = 0; i < EE_EE.size(); i++)
//		{
//			Vector5i cont_EE = EE_EE[i];
//			int E1 = cont_EE[1], E2 = cont_EE[2];
//			int e1p1 = tetSimMesh.index_boundaryEdge[E1][0], e1p2 = tetSimMesh.index_boundaryEdge[E1][1],
//				e2p1 = tetSimMesh.index_boundaryEdge[E2][0], e2p2 = tetSimMesh.index_boundaryEdge[E2][1];
//
//			Eigen::Vector3d P1 = tetSimMesh.pos_node[e1p1];
//			Eigen::Vector3d P2 = tetSimMesh.pos_node[e1p2];
//			Eigen::Vector3d Q1 = tetSimMesh.pos_node[e2p1];
//			Eigen::Vector3d Q2 = tetSimMesh.pos_node[e2p2];
//			int type = cont_EE[3];
//
//			double dis2 = 0;
//			DIS::computeEdgeEdgeD(P1, P2, Q1, Q2, dis2);
//			Eigen::Vector4i ptIndices = { e1p1 , e1p2 , e2p1 , e2p2 };
//			int actualStartIndex_hess = startIndex_hess + i * 144,
//				actualStartIndex_grad = startIndex_grad + i * 12;
//			BarrierEnergy::gradAndHess_EE(hessian_triplet, grad_triplet, actualStartIndex_hess,
//				actualStartIndex_grad, ptIndices, type,
//				dis2, tetSimMesh, parameters);
//
//		}
//
//	}
//
//
//
//
//	//double endTime3 = omp_get_wtime();
//	//std::cout << "	Cal grad_hess Time is : " << endTime3 - endTime2 << "s" << std::endl;
//
//	// assemable the left-hand side 
//	Eigen::SparseMatrix<double> leftHandSide(3 * tetSimMesh.pos_node.size(), 3 * tetSimMesh.pos_node.size());
//	leftHandSide.setFromTriplets(hessian_triplet.begin(), hessian_triplet.end());
//
//	// assemable the right-hand side 
//	Eigen::VectorXd rightHandSide = Eigen::VectorXd(tetSimMesh.pos_node.size() * 3);
//	rightHandSide.setZero();
//	for (int i = 0; i < grad_triplet.size(); i++)
//	{
//		std::pair<int, double> ele = grad_triplet[i];
//		rightHandSide[ele.first] += ele.second;
//	}
//
//	//double endTime4 = omp_get_wtime();
//	//std::cout << "	Assemble grad_hess Time is : " << endTime4 - endTime3 << "s" << std::endl;
//
//	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
//	solver.compute(leftHandSide);
//	Eigen::VectorXd result = solver.solve(-rightHandSide);
//
//
//
//#pragma omp parallel for num_threads(parameters.numOfThreads)
//	for (int i = 0; i < tetSimMesh.pos_node.size(); i++)
//	{
//		movingDir[i] = result.block<3, 1>(3 * i, 0);
//	}
//
//
//
//	return movingDir;
//}
//
//
