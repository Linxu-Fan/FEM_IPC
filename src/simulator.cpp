#include "simulator.h" 
#include "ExternalEnergy.h" 
#include "InertiaEnergy.h" 
#include "ElasticEnergy.h" 

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
		std::vector<Eigen::Vector3d> direction = solve_linear_system(tetMesh, parameters, timestep);
		double dist_to_converge = infiniteNorm(direction);

		int iteration = 0;
		while (dist_to_converge / parameters.dt > parameters.searchResidual)
		{
			iteration += 1;

			double step = 1.0;
			step_forward(tetMesh, currentPosition, direction, step);
			double newEnergyVal = compute_IP_energy(tetMesh, parameters, timestep);			
			while (newEnergyVal > lastEnergyVal)
			{
				step /= 2.0;
				step_forward(tetMesh, currentPosition, direction, step);
				newEnergyVal = compute_IP_energy(tetMesh, parameters, timestep);
			}
			currentPosition = tetMesh.pos_node;
			lastEnergyVal = newEnergyVal;

			direction = solve_linear_system(tetMesh, parameters, timestep);
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
double compute_IP_energy(Mesh& tetMesh, FEMParamters& parameters, int timestep)
{
	double energyVal = 0;
	tetMesh.update_F();
	
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
		energyVal += ElasticEnergy::Val(tetMesh.materialMesh, parameters.model, tetMesh.tetra_F[eI], parameters.dt, tetMesh.tetra_vol[eI]);
	}	

	return energyVal;
}

// compute energy gradient and hessian of the linear system and solve the syetem
std::vector<Eigen::Vector3d> solve_linear_system(Mesh& tetMesh, FEMParamters& parameters, int timestep)
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
		Eigen::Matrix<double, 9, 12> dFdx = ElasticEnergy::dF_wrt_dx(tetMesh.tetra_DM_inv[eI]);
		std::vector<std::pair<int, double>> elasEngGrad = ElasticEnergy::Grad(tetMesh.materialMesh, parameters.model, tetMesh.tetra_F[eI], parameters.dt, tetMesh.tetra_vol[eI], dFdx, tetMesh.tetrahedrals[eI], tetVertInd_BC);
		std::vector<Eigen::Triplet<double>> elasEngHess = ElasticEnergy::Hess(tetMesh.materialMesh, parameters.model, tetMesh.tetra_F[eI], parameters.dt, tetMesh.tetra_vol[eI], dFdx, tetMesh.tetrahedrals[eI], tetVertInd_BC);
		grad_triplet.insert(grad_triplet.end(), elasEngGrad.begin(), elasEngGrad.end());
		hessian_triplet.insert(hessian_triplet.end(), elasEngHess.begin(), elasEngHess.end());
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


