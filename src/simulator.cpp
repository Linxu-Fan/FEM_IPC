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
		std::vector<BarrierEnergyRes> pTeEBarrVec;



		double lastEnergyVal = compute_IP_energy(tetMesh, parameters, timestep, pTeEBarrVec);
		std::vector<Eigen::Vector3d> direction = solve_linear_system(tetMesh, parameters, timestep, pTeEBarrVec);
		double dist_to_converge = infiniteNorm(direction);

		int iteration = 0;
		while (dist_to_converge / parameters.dt > parameters.searchResidual)
		{
			iteration += 1;

			double step = 1.0;
			step_forward(tetMesh, currentPosition, direction, step);
			double newEnergyVal = compute_IP_energy(tetMesh, parameters, timestep, pTeEBarrVec);
			while (newEnergyVal > lastEnergyVal)
			{
				step /= 2.0;
				step_forward(tetMesh, currentPosition, direction, step);
				newEnergyVal = compute_IP_energy(tetMesh, parameters, timestep, pTeEBarrVec);
			}
			currentPosition = tetMesh.pos_node;
			lastEnergyVal = newEnergyVal;

			direction = solve_linear_system(tetMesh, parameters, timestep, pTeEBarrVec);
			dist_to_converge = infiniteNorm(direction);
		
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

	//if (tetMesh.pos_node[vertInd][0] <= 0.15)
	//{
	//	extForce = { -100,0, 0 };
	//}
	//if (tetMesh.pos_node[vertInd][0] >= 0.85)
	//{
	//	extForce = {100,0, 0 };
	//}

	return extForce;
}

// compute the incremental potential energy
double compute_IP_energy(Mesh& tetMesh, FEMParamters& parameters, int timestep, std::vector<BarrierEnergyRes>& pTeEBarrVec)
{
	double energyVal = 0;
	tetMesh.update_F();
	pTeEBarrVec.clear();
	
	// energy contribution per vertex
	for (int vI = 0; vI < tetMesh.pos_node.size(); vI++)
	{
		double nodeMass = tetMesh.mass_node[vI];
		Eigen::Vector3d xt = tetMesh.pos_node_prev[vI];
		Eigen::Vector3d v = tetMesh.vel_node[vI];
		Eigen::Vector3d x = tetMesh.pos_node[vI];
		Eigen::Vector3d extForce = compute_external_force(tetMesh, vI, timestep);

		// the external energy contribution
		energyVal += ExternalEnergy::Val(nodeMass, parameters.dt, x, parameters, extForce);
		// the inertia energy contribution
		energyVal += InertiaEnergy::Val(nodeMass, parameters.dt, xt, v, x, extForce, parameters);
	}

	// energy contribution per element
	for (int eI = 0; eI < tetMesh.tetra_F.size(); eI++)
	{
		// the internal elastic energy contribution
		int matInd = tetMesh.materialInd[eI];
		energyVal += ElasticEnergy::Val(tetMesh.materialMesh[matInd], parameters.model, tetMesh.tetra_F[eI], parameters.dt, tetMesh.tetra_vol[eI]);
	}	

	// energy contribution from barrier
	// point-triangle barrier
	for (std::map<int, std::set<int>>::iterator it = tetMesh.boundaryVertices.begin(); it != tetMesh.boundaryVertices.end(); it++)
	{
		int ptInd = it->first;
		Eigen::Vector3d P = tetMesh.pos_node[ptInd];
		for (int tI = 0; tI < tetMesh.boundaryTriangles.size(); tI++)
		{
			if (it->second.find(tI) == it->second.end()) // this triangle is not incident with the point
			{
				Eigen::Vector3i tri = tetMesh.boundaryTriangles[tI];
				Eigen::Vector3d A = tetMesh.pos_node[tri[0]];
				Eigen::Vector3d B = tetMesh.pos_node[tri[1]];
				Eigen::Vector3d C = tetMesh.pos_node[tri[2]];

				int type = pointTriangleDisType(P, A, B, C);
				double dis2 = pointTriangleDis2(type, P, A, B, C);
				BarrierEnergyRes res;
				if (dis2 <= squaredDouble(parameters.IPC_dis)) // only calculate the energy when the distance is smaller than the threshold
				{
					Eigen::Vector4i vtInd = { ptInd , tri[0] , tri[1] , tri[2] };
					Eigen::Vector4i vtInd_BC = { tetMesh.boundaryCondition_node[ptInd].type, tetMesh.boundaryCondition_node[tri[0]].type , tetMesh.boundaryCondition_node[tri[1]].type , tetMesh.boundaryCondition_node[tri[2]].type };
					res.Val = BarrierEnergy::Val(true, dis2, tetMesh, vtInd, parameters.IPC_dis, 1.0, parameters.dt);
					res.Grad = BarrierEnergy::Grad(true, dis2, type, tetMesh, vtInd, vtInd_BC, parameters.IPC_dis, 1.0, parameters.dt);
					res.Hess = BarrierEnergy::Hess(true, dis2, type, tetMesh, vtInd, vtInd_BC, parameters.IPC_dis, 1.0, parameters.dt);
				}
				pTeEBarrVec.push_back(res);
				energyVal += res.Val;
			}
		}
	}
	// edge-edge barrier
	for (std::map<int, std::map<int, Eigen::Vector2i>>::iterator ite11 = tetMesh.boundaryEdges.begin(); ite11 != tetMesh.boundaryEdges.end(); ite11++)
	{
		int e1p1 = ite11->first;
		for (std::map<int, Eigen::Vector2i>::iterator ite12 = ite11->second.begin(); ite12 != ite11->second.end(); ite12++)
		{
			int e1p2 = ite12->first;

			for (std::map<int, std::map<int, Eigen::Vector2i>>::iterator ite21 = tetMesh.boundaryEdges.begin(); ite21 != tetMesh.boundaryEdges.end(); ite21++)
			{
				int e2p1 = ite21->first;
				for (std::map<int, Eigen::Vector2i>::iterator ite22 = ite11->second.begin(); ite22 != ite11->second.end(); ite22++)
				{
					int e2p2 = ite22->first;

					if (e1p1 != e2p1 && e1p1 != e2p2 && e1p2 != e2p1 && e1p2 != e2p2) // not duplicated and incident edges
					{
						Eigen::Vector3i P1 = tetMesh.pos_node[e1p1];
						Eigen::Vector3d P2 = tetMesh.pos_node[e1p2];
						Eigen::Vector3d Q1 = tetMesh.pos_node[e2p1];
						Eigen::Vector3d Q2 = tetMesh.pos_node[e2p2];

						int type = edgeEdgeDisType(P1, P2, Q1, Q2);
						double dis2 = edgeEdgeDis2(type, P1, P2, Q1, Q2);
						BarrierEnergyRes res;
						if (dis2 <= squaredDouble(parameters.IPC_dis)) // only calculate the energy when the distance is smaller than the threshold
						{
							Eigen::Vector4i vtInd = { e1p1 , e1p2 , e2p1 , e2p2 };
							Eigen::Vector4i vtInd_BC = { tetMesh.boundaryCondition_node[e1p1].type, tetMesh.boundaryCondition_node[e1p2].type , tetMesh.boundaryCondition_node[e2p1].type , tetMesh.boundaryCondition_node[e2p2].type };
							res.Val = BarrierEnergy::Val(false, dis2, tetMesh, vtInd, parameters.IPC_dis, 1.0, parameters.dt);
							res.Grad = BarrierEnergy::Grad(false, dis2, type, tetMesh, vtInd, vtInd_BC, parameters.IPC_dis, 1.0, parameters.dt);
							res.Hess = BarrierEnergy::Hess(false, dis2, type, tetMesh, vtInd, vtInd_BC, parameters.IPC_dis, 1.0, parameters.dt);
						}
						pTeEBarrVec.push_back(res);
						energyVal += res.Val;
					}


				}
			}
		}
	}


	return energyVal;
}

// compute energy gradient and hessian of the linear system and solve the syetem
std::vector<Eigen::Vector3d> solve_linear_system(Mesh& tetMesh, FEMParamters& parameters, int timestep, std::vector<BarrierEnergyRes>& pTeEBarrVec)
{
	std::vector<Eigen::Vector3d> movingDir(tetMesh.pos_node.size());

	// update defromation gradient
	tetMesh.update_F();

	// hessian is the left-hand side, and grad is the right-hand side
	std::vector<Eigen::Triplet<double>> hessian_triplet;
	std::vector<std::pair<int, double>> grad_triplet;

	// energy contribution per vertex
	for (int vI = 0; vI < tetMesh.pos_node.size(); vI++)
	{
		double nodeMass = tetMesh.mass_node[vI];
		Eigen::Vector3d xt = tetMesh.pos_node_prev[vI];
		Eigen::Vector3d v = tetMesh.vel_node[vI];
		Eigen::Vector3d x = tetMesh.pos_node[vI];
		Eigen::Vector3d extForce = compute_external_force(tetMesh, vI, timestep);

		// the external energy contribution
		std::vector<std::pair<int, double>> extEngGrad = ExternalEnergy::Grad(nodeMass, parameters.dt, x, parameters, extForce, vI, tetMesh.boundaryCondition_node[vI].type);
		grad_triplet.insert(grad_triplet.end(), extEngGrad.begin(), extEngGrad.end());

		// the inertia energy contribution
		std::vector<std::pair<int, double>> inerEngGrad = InertiaEnergy::Grad(nodeMass, parameters.dt, xt, v, x, extForce, vI, parameters, tetMesh.boundaryCondition_node[vI].type);
		std::vector<Eigen::Triplet<double>> inerEngHess = InertiaEnergy::Hess(nodeMass, vI, tetMesh.boundaryCondition_node[vI].type);
		grad_triplet.insert(grad_triplet.end(), inerEngGrad.begin(), inerEngGrad.end());
		hessian_triplet.insert(hessian_triplet.end(), inerEngHess.begin(), inerEngHess.end());		

	}

	// energy contribution per element
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
		grad_triplet.insert(grad_triplet.end(), elasEngGrad.begin(), elasEngGrad.end());
		hessian_triplet.insert(hessian_triplet.end(), elasEngHess.begin(), elasEngHess.end());
	}

	// energy contribution from barrier
	for (int ptI = 0; ptI < pTeEBarrVec.size(); ptI++)
	{
		grad_triplet.insert(grad_triplet.end(), pTeEBarrVec[ptI].Grad.begin(), pTeEBarrVec[ptI].Grad.end());
		hessian_triplet.insert(hessian_triplet.end(), pTeEBarrVec[ptI].Hess.begin(), pTeEBarrVec[ptI].Hess.end());
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

	std::cout << "st" << std::endl;
	// apply the fixed boundary condition
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
	std::cout << "en" << std::endl;

	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
	solver.compute(leftHandSide);
	Eigen::VectorXd result = solver.solve(-rightHandSide);
	
	for (int i = 0; i < tetMesh.pos_node.size(); i++)
	{
		movingDir[i] = result.block<3, 1>(3 * i, 0);
	}

	return movingDir;
}


// move points' position according to the direction; Note this is a trial movement
void step_forward(Mesh& tetMesh, std::vector<Eigen::Vector3d>& currentPosition, std::vector<Eigen::Vector3d>& direction, double step)
{
	for (int vI = 0; vI < tetMesh.pos_node.size(); vI++)
	{
		tetMesh.pos_node[vI] = currentPosition[vI] + step * direction[vI];
	}
}


