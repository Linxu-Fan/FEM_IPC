#include "simulator.h" 
#include "ExternalEnergy.h" 
#include "InertiaEnergy.h" 
#include "ElasticEnergy.h" 

// implicit integration
void implicitFEM(Mesh& tetMesh, FEMParamters& parameters)
{
	for (int timestep = 0; timestep < parameters.num_timesteps; timestep++)
	{
		tetMesh.pos_node_prev = tetMesh.pos_node;
		std::vector<Eigen::Vector3d> currentPosition = tetMesh.pos_node;

		double dist_to_converge = 1.0e10;
		double lastEnergyVal = compute_IP_energy(tetMesh, parameters, timestep);
		do
		{
			// line search
			double step = 1.0;
			std::vector<Eigen::Vector3d> direction = solve_linear_system(tetMesh, parameters, timestep);
			step_forward(tetMesh, currentPosition, direction, step);
			double newEnergyVal = compute_IP_energy(tetMesh, parameters, timestep);
			while (newEnergyVal > lastEnergyVal)
			{
				step /= 2.0;
				step_forward(tetMesh, currentPosition, direction, step);		
				newEnergyVal = compute_IP_energy(tetMesh, parameters, timestep);
			}
			lastEnergyVal = newEnergyVal;

			dist_to_converge = infiniteNorm(direction);
		} while (dist_to_converge / parameters.dt > parameters.searchResidual);

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

	if (tetMesh.pos_node[vertInd][2] <= 0.15)
	{
		extForce = { -100,0, 0 };
	}
	if (tetMesh.pos_node[vertInd][2] >= 0.85)
	{
		extForce = {100,0, 0 };
	}

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
		energyVal += InertiaEnergy::Val(nodeMass, parameters.dt, xt, v, x, extForce);
	}


	// energy contribution per element
	for (int eI = 0; eI < tetMesh.tetra_F.size(); eI++)
	{
		// the internal elastic energy contribution
		Eigen::Matrix3d tmp = tetMesh.tetra_F[eI];
		energyVal += ElasticEnergy::Val(tetMesh.materialMesh, parameters.model, tetMesh.tetra_F[eI], parameters.dt);
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
		std::vector<std::pair<int, double>> extEngGrad = ExternalEnergy::Grad(nodeMass, parameters.dt, x, parameters, extForce, vI);
		grad_triplet.insert(grad_triplet.end(), extEngGrad.begin(), extEngGrad.end());

		// the inertia energy contribution
		std::vector<std::pair<int, double>> inerEngGrad = InertiaEnergy::Grad(nodeMass, parameters.dt, xt, v, x, extForce, vI);
		std::vector<Eigen::Triplet<double>> inerEngHess = InertiaEnergy::Hess(nodeMass, vI);
		grad_triplet.insert(grad_triplet.end(), inerEngGrad.begin(), inerEngGrad.end());
		hessian_triplet.insert(hessian_triplet.end(), inerEngHess.begin(), inerEngHess.end());			
	}

	// energy contribution per element
	for (int eI = 0; eI < tetMesh.tetrahedrals.size(); eI++)
	{
		// the internal elastic energy contribution
		Eigen::Matrix<double, 9, 12> dFdx = ElasticEnergy::dF_wrt_dx(tetMesh.tetra_DM_inv[eI]);
		std::vector<std::pair<int, double>> elasEngGrad = ElasticEnergy::Grad(tetMesh.materialMesh, parameters.model, tetMesh.tetra_F[eI], parameters.dt, tetMesh.tetra_vol[eI], dFdx, tetMesh.tetrahedrals[eI]);
		std::vector<Eigen::Triplet<double>> elasEngHess = ElasticEnergy::Hess(tetMesh.materialMesh, parameters.model, tetMesh.tetra_F[eI], parameters.dt, tetMesh.tetra_vol[eI], dFdx, tetMesh.tetrahedrals[eI]);
		grad_triplet.insert(grad_triplet.end(), elasEngGrad.begin(), elasEngGrad.end());
		hessian_triplet.insert(hessian_triplet.end(), elasEngHess.begin(), elasEngHess.end());
	}

	// assemable the left-hand side 
	Eigen::SparseMatrix<double> leftHandSide(3 * tetMesh.pos_node.size(), 3 * tetMesh.pos_node.size());
	leftHandSide.setFromTriplets(hessian_triplet.begin(), hessian_triplet.end());

	// assemable the right-hand side 
	Eigen::VectorXd rightHandSide = Eigen::VectorXd(tetMesh.pos_node.size() * 3);
	for (int i = 0; i < grad_triplet.size(); i++)
	{
		std::pair<int, double> ele = grad_triplet[i];
		rightHandSide[ele.first] = ele.second;
	}

	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
	solver.compute(-leftHandSide);
	Eigen::VectorXd result = solver.solve(rightHandSide);
	
	for (int i = 0; i < tetMesh.pos_node.size(); i++)
	{
		movingDir[i] = result.block<3, 1>(3 * i, 0);
	}

	return movingDir;
}

// do line search to find the optimal step
void lineSearch(Mesh& tetMesh, std::vector<Eigen::Vector3d>& direction, FEMParamters& parameters, double& lastEnergyVal, int timestep)
{
	std::vector<Eigen::Vector3d> currentPosition = tetMesh.pos_node;
	double step = 1.0;
	step_forward(tetMesh, currentPosition, direction, step);
	double newEnergyVal = compute_IP_energy(tetMesh, parameters, timestep);
	while (newEnergyVal > lastEnergyVal)
	{
		step_forward(tetMesh, currentPosition, direction, step);
		step /= 2.0;
		newEnergyVal = compute_IP_energy(tetMesh, parameters, timestep);
	}
	lastEnergyVal = newEnergyVal;
}

// move points' position according to the direction
void step_forward(Mesh& tetMesh, std::vector<Eigen::Vector3d>& currentPosition, std::vector<Eigen::Vector3d>& direction, double step)
{
	for (int vI = 0; vI < tetMesh.pos_node.size(); vI++)
	{
		tetMesh.pos_node[vI] = currentPosition[vI] + step * direction[vI];
	}
}
