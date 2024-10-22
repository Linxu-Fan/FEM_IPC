#include "simulator.h" 
#include "tools.h" 

// TODO
// 1. Accelerate and parallel openMP
// 2. If the code doesn't work, check line45 of InertiaEnergy.cpp


int main()
{

	if (1)
	{
		



		/*FEMParamters parameters;
		parameters.IPC_dis = 0.01;
		parameters.IPC_kStiffness = 1.0e16;
		parameters.numOfThreads = 24;
		Eigen::Vector3d bbx_min = {23.5, -99, -99}, bbx_max = {26.5, 99 , 99};


		std::vector<Eigen::Vector3d> forceFrame(500, Eigen::Vector3d::Zero());
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int frame = 0 ; frame < 500; frame++)
		{		
			objMeshFormat testMesh;
			testMesh.readObjFile("E:/hydroStatic_object/Libuipc/libuipc/output/tests/sim_case/25_linear_arap_beam.cpp/scene_surface" + std::to_string(frame) + ".obj");
			testMesh.sepConnectedComponents();

			Eigen::Vector3d force = compute_contact_force(testMesh.componentsSep[0], testMesh.componentsSep[3], bbx_min, bbx_max, parameters);
			forceFrame[frame] = force;

			std::cout<<"Frame = "<<frame<<"; Force = ("<< " " << force[0] << " " << force[1] << " " << force[2]<<")" << std::endl;
			std::cout << std::endl;
		}






		std::ofstream outfile9("./output/force.txt", std::ios::trunc);
		for (int frame = 0; frame < 500; frame++)
		{
			outfile9 << std::scientific << std::setprecision(8) << frame << " " << forceFrame[frame][0] << " " << forceFrame[frame][1] << " " << forceFrame[frame][2] << std::endl;
		}
		outfile9.close();*/




		//double pi = 3.141592653;

		//Eigen::Matrix3d F = Eigen::Matrix3d::Identity();


		//// Rotate along x,y,z axis
		//Eigen::Matrix3d RX = Eigen::Matrix3d::Zero();
		//Eigen::Matrix3d RY = Eigen::Matrix3d::Zero();
		//Eigen::Matrix3d RZ = Eigen::Matrix3d::Zero();
		//double tx = 200, ty = 142, tz = 1450;
		//RX(0, 0) = 1.0;
		//RX(1, 1) = cos(tx);
		//RX(1, 2) = -sin(tx);
		//RX(2, 1) = sin(tx);
		//RX(2, 2) = cos(tx);

		//RY(1, 1) = 1.0;
		//RY(0, 0) = cos(ty);
		//RY(0, 2) = -sin(ty);
		//RY(2, 0) = sin(ty);
		//RY(2, 2) = cos(ty);

		//RZ(2, 2) = 1.0;
		//RZ(1, 1) = cos(tz);
		//RZ(0, 1) = -sin(tz);
		//RZ(1, 0) = sin(tz);
		//RZ(0, 0) = cos(tz);

		//F = F * RX * RY * RZ;


		//// Rotate along arbitrary axis
		//double theta = 458;
		//double cosTheta = std::cos(theta);
		//double sinTheta = std::sin(theta);
		//Eigen::Vector3d axis = {1.5,1.9,74.4};
		//axis.normalize();
		//Eigen::Matrix3d outerProduct = axis * axis.transpose();
		//Eigen::Matrix3d K;
		//K << 0, -axis.z(), axis.y(),
		//	 axis.z(), 0, -axis.x(),
		//	 -axis.y(), axis.x(), 0;
		//Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
		//// Rodrigues' rotation formula
		//Eigen::Matrix3d R = cosTheta * I + (1 - cosTheta) * outerProduct + sinTheta * K;

		//F = F * R;



		//Eigen::HouseholderQR<Eigen::MatrixXd> qr(F);
		//Eigen::MatrixXd QW = qr.householderQ();
		//Eigen::MatrixXd RW = qr.matrixQR().triangularView<Eigen::Upper>();
		//Eigen::MatrixXd A_reconstructed = QW * RW;
		//std::cout << "Matrix F:\n" << F << "\n" << std::endl;
		////std::cout << "Matrix Q:\n" << QW << "\n" << std::endl;
		////std::cout << "Matrix R:\n" << RW << "\n" << std::endl;
		////std::cout << "Q * R:\n" << A_reconstructed << "\n" << std::endl;



		//// Add stretch matrix
		//double s = 0.5;
		//Eigen::Matrix3d stretch = Eigen::Matrix3d::Zero();
		//stretch(0, 0) = s;
		//stretch(1, 1) = s;
		//stretch(2, 2) = s;
		//F = F * stretch;


		//F = Eigen::Matrix3d::Zero();
		//F(1, 1) = sqrt(2);
		//F(0, 0) = sqrt(2);

		//double xv2 = F(0,0) * F(0, 0) * PI / 2.0 * PI + F(0, 1) * F(0, 1) * PI / 2.0 * PI + F(0, 2) * F(0, 2) * PI * PI + 2 * F(0, 0) * F(0, 2) * 2.0 * PI ;
		//double yv2 = F(1,0) * F(1, 0) * PI / 2.0 * PI + F(1, 1) * F(1, 1) * PI / 2.0 * PI + F(1, 2) * F(1, 2) * PI * PI + 2 * F(1, 0) * F(1, 2) * 2.0 * PI ;
		//double zv2 = F(2,0) * F(2, 0) * PI / 2.0 * PI + F(2, 1) * F(2, 1) * PI / 2.0 * PI + F(2, 2) * F(2, 2) * PI * PI + 2 * F(2, 0) * F(2, 2) * 2.0 * PI ;

		//std::cout << "f(phi) = " << (xv2 + yv2 + zv2) / 2 / PI / PI   << std::endl;

		//// Test axis change_y
		//double xv2_2 = F(0, 0) * F(0, 0) * PI / 2.0 * PI + F(0, 2) * F(0, 2) * PI / 2.0 * PI + F(0, 1) * F(0, 1) * PI * PI + 2 * F(0, 0) * F(0, 1) * 2.0 * PI;
		//double yv2_2 = F(1, 0) * F(1, 0) * PI / 2.0 * PI + F(1, 2) * F(1, 2) * PI / 2.0 * PI + F(1, 1) * F(1, 1) * PI * PI + 2 * F(1, 0) * F(1, 1) * 2.0 * PI;
		//double zv2_2 = F(2, 0) * F(2, 0) * PI / 2.0 * PI + F(2, 2) * F(2, 2) * PI / 2.0 * PI + F(2, 1) * F(2, 1) * PI * PI + 2 * F(2, 0) * F(2, 1) * 2.0 * PI;

		//std::cout << "f_2(phi) = " << (xv2_2 + yv2_2 + zv2_2) / 2 / PI / PI  << std::endl;

		//// Test axis change_z
		//double xv2_3 = F(0, 1) * F(0, 1) * PI / 2.0 * PI + F(0, 2) * F(0, 2) * PI / 2.0 * PI + F(0, 0) * F(0, 0) * PI * PI / 2.0;
		//double yv2_3 = F(1, 1) * F(1, 1) * PI / 2.0 * PI + F(1, 2) * F(1, 2) * PI / 2.0 * PI + F(1, 0) * F(1, 0) * PI * PI / 2.0;
		//double zv2_3 = F(2, 1) * F(2, 1) * PI / 2.0 * PI + F(2, 2) * F(2, 2) * PI / 2.0 * PI + F(2, 0) * F(2, 0) * PI * PI / 2.0;

		//std::cout << "f_3(phi) = " << (xv2_3 + yv2_3 + zv2_3) / 2 / PI / PI  << std::endl;

		//std::cout << "F_00 * F_02 + F_10 * F_12 + F_20 * F_22 = " << F(0, 0) * F(0, 2) + F(1, 0) * F(1, 2) + F(2, 0) * F(2, 2) << std::endl;
	}
	else
	{
		// Case:
		// 0: Three-point-bending beam test
		// 1. Two tetrahedrals collision test to verify the correctness of the solver
		// 2. Single tetrahedral test
		// 3. ABD test
		int caseNum = 0;
		if (caseNum == 0)
		{
			Material mat1;
			mat1.density = 2780;
			mat1.E = 7.26e7;
			mat1.updateDenpendecies();


			
			std::vector<meshConfiguration> config;			
			meshConfiguration m1, m2, m3, m4;
			m1.filePath = "../input/beam.msh";
			m1.mesh_material = mat1;
			m1.scale = {1000, 1000, 1000};
			m1.note = "beam";
			config.push_back(m1);

			m2 = m1;
			m2.filePath = "../input/left_support.msh";
			m2.note = "left_support";
			m2.translation = { 0, 0, 1.0e-4};
			//config.push_back(m2);

			m3 = m1;
			m3.filePath = "../input/middle_support.msh";
			m3.note = "middle_support";
			m3.translation = { 0, 0, -1.0e-4 };
			//config.push_back(m3);

			m4 = m1;
			m4.filePath = "../input/impactor.msh";
			m4.note = "impactor";
			m4.translation = { 0, 0, 1.5e-2 };
			config.push_back(m4);




			Mesh tetMesh;
			std::vector<std::string> objectNames;
			for (int i = 0; i < config.size(); i++)
			{
				objectNames.push_back(config[i].note);
			}
			tetMesh.readMeshes(config);
			tetMesh.initializeMesh();
			tetMesh.surfaceMesh.outputFile("surfMesh");


			std::cout << "tetMesh.pos_node.size() = " << tetMesh.pos_node.size() << std::endl;
			std::cout << "tetMesh.tetrahedrals.size() = " << tetMesh.tetrahedrals.size() << std::endl;


			FEMParamters parameters;
			parameters.gravity = { 0, 0, 0 };
			parameters.num_timesteps = 3000;
			parameters.numOfThreads = 12;
			parameters.dt = 5.0e-4;
			parameters.outputFrequency = 20;
			parameters.enableGround = false;
			parameters.searchResidual = 5.0;
			parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
			parameters.rigidMode = true;
			parameters.objectNames = objectNames;
			parameters.IPC_dis = 0.01;
			parameters.IPC_eta = 0.05;
			parameters.IPC_kStiffness = 1.0e14;
			parameters.IPC_hashSize = tetMesh.calLargestEdgeLength() * 1.1;
			parameters.IPC_B3Stiffness = 500;



			for (int p = 0; p < tetMesh.pos_node.size(); p++)
			{
				if (tetMesh.note_node[p] == "impactor")
				{
					tetMesh.boundaryCondition_node[p].type = 1;
					for (int fra = 0; fra < parameters.num_timesteps; fra++)
					{
						double incre = 0.1 / (double)parameters.num_timesteps * (double)fra;
						Eigen::Vector3d inc = {0,0,-incre };
						Eigen::Vector3d desiredPos = inc + tetMesh.pos_node[p];

						tetMesh.boundaryCondition_node[p].location.push_back(desiredPos);
					}
				}

				if (tetMesh.note_node[p] == "beam")
				{
					if (tetMesh.pos_node[p][0] <= -22)
					{
						tetMesh.boundaryCondition_node[p].type = 1;
						for (int fra = 0; fra < parameters.num_timesteps; fra++)
						{
							Eigen::Vector3d desiredPos = tetMesh.pos_node[p];
							tetMesh.boundaryCondition_node[p].location.push_back(desiredPos);
						}
					}
				}

			}



			implicitFEM(tetMesh, parameters);


		}
		else if (caseNum == 1)
		{
			Material mat1;
			mat1.density = 1000;
			mat1.E = 7.26e5;
			mat1.updateDenpendecies();



			std::vector<meshConfiguration> config;
			meshConfiguration m1,m2;
			m1.filePath = "../input/tet_origin.msh";
			m1.mesh_material = mat1;
			m1.translation = {0,0,1.1};
			m1.velocity = {0,0,1};
			m1.note = "tet_origin";
			config.push_back(m1);

			m2.filePath = "../input/tet_origin.msh";
			m2.mesh_material = mat1;
			m2.translation = { 0,0,4.3 };
			m2.velocity = { 0,0,-1 };
			m2.note = "tet_origin_top";
			config.push_back(m2);



			Mesh tetMesh;
			std::vector<std::string> objectNames;
			for (int i = 0; i < config.size(); i++)
			{
				objectNames.push_back(config[i].note);
			}
			tetMesh.readMeshes(config);
			tetMesh.initializeMesh();
			tetMesh.surfaceMesh.outputFile("surfMesh");


			std::cout << "tetMesh.pos_node.size() = " << tetMesh.pos_node.size() << std::endl;
			std::cout << "tetMesh.tetrahedrals.size() = " << tetMesh.tetrahedrals.size() << std::endl;


			FEMParamters parameters;
			parameters.gravity = { 0, 0, 0 };
			parameters.num_timesteps = 3;
			parameters.numOfThreads = 1;
			parameters.dt = 1.0e-2;
			parameters.outputFrequency = 10;
			parameters.enableGround = true;
			parameters.searchResidual = 0.001;
			parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
			parameters.rigidMode = true;
			parameters.objectNames = objectNames;
			parameters.IPC_dis = 0.001;
			parameters.IPC_eta = 0.01;
			parameters.IPC_kStiffness = 1.0e12;
			parameters.IPC_hashSize = tetMesh.calLargestEdgeLength() * 1.3;
			parameters.IPC_B3Stiffness = 500;




			implicitFEM(tetMesh, parameters);


		}
		else if (caseNum == 2)
		{
			Material mat1;
			mat1.density = 5000;
			mat1.E = 7.26e10;
			mat1.updateDenpendecies();



			std::vector<meshConfiguration> config;
			meshConfiguration m1, m2;
			m1.filePath = "../input/tet_origin.msh";
			m1.mesh_material = mat1;
			m1.translation = { 0,0,1.1 };
			m1.velocity = { 0,0,0 };
			m1.note = "tet_origin";
			config.push_back(m1);




			Mesh tetMesh;
			std::vector<std::string> objectNames;
			for (int i = 0; i < config.size(); i++)
			{
				objectNames.push_back(config[i].note);
			}
			tetMesh.readMeshes(config);
			tetMesh.initializeMesh();
			tetMesh.surfaceMesh.outputFile("surfMesh");


			std::cout << "tetMesh.pos_node.size() = " << tetMesh.pos_node.size() << std::endl;
			std::cout << "tetMesh.tetrahedrals.size() = " << tetMesh.tetrahedrals.size() << std::endl;


			FEMParamters parameters;
			parameters.gravity = { 0, 0, 0 };
			parameters.num_timesteps = 50000;
			parameters.numOfThreads = 1;
			parameters.dt = 1.0e-4;
			parameters.outputFrequency = 1000;
			parameters.enableGround = true;
			parameters.searchResidual = 0.001;
			parameters.model = "ARAP"; // neoHookean ARAP ARAP_linear ACAP
			parameters.rigidMode = true;
			parameters.objectNames = objectNames;
			parameters.IPC_dis = 0.001;
			parameters.IPC_eta = 0.05;
			parameters.IPC_kStiffness = 1.0e12;
			parameters.IPC_hashSize = tetMesh.calLargestEdgeLength() * 1.3;
			parameters.IPC_B3Stiffness = 500;


			//tetMesh.pos_node[3] = {0,0.05,1.15};
			//tetMesh.update_F(1);


			//// Perform SVD on the matrix A
			//Eigen::JacobiSVD<Eigen::MatrixXd> svd(tetMesh.tetra_F[0], Eigen::ComputeFullU | Eigen::ComputeFullV);
			//Eigen::MatrixXd U = svd.matrixU();
			//Eigen::VectorXd S = svd.singularValues();
			//Eigen::MatrixXd V = svd.matrixV();

			//// Compute R and Q
			//Eigen::MatrixXd R = U * V.transpose();
			//Eigen::MatrixXd Q = V * S.asDiagonal() * V.transpose();

			//std::cout << "Matrix F:\n" << tetMesh.tetra_F[0] << "\n\n";
			//std::cout << "Rotation Matrix R:\n" << R << "\n\n";
			//std::cout << "Stretch Matrix Q:\n" << Q << "\n\n";
			//


			//Eigen::Matrix<double, 9, 12> DFDX = ElasticEnergy::dF_wrt_dx(tetMesh.tetra_DM_inv[0]);
			//Vector9d grad_Energy_wrt_F = flatenMatrix3d(R);
			//Eigen::Matrix<double, 12, 1> engGrad = DFDX.transpose() * grad_Energy_wrt_F;

			//std::cout << "Force_1 = (" << engGrad(0, 0) << "," << engGrad(1, 0) << "," << engGrad(2, 0) << ")" << std::endl;
			//std::cout << "Force_2 = (" << engGrad(3, 0) << "," << engGrad(4, 0) << "," << engGrad(5, 0) << ")" << std::endl;
			//std::cout << "Force_3 = (" << engGrad(6, 0) << "," << engGrad(7, 0) << "," << engGrad(8, 0) << ")" << std::endl;
			//std::cout << "Force_4 = (" << engGrad(9, 0) << "," << engGrad(10, 0) << "," << engGrad(11, 0) << ")" << std::endl;

			//

			// Explicit FEM
			{
				std::vector<Eigen::Vector3d> nodeForce;

				// Simulation loop
				for (int step = 0; step < parameters.num_timesteps; ++step) 
				{
					if (step % 100 == 0)
					{
						std::cout << "Timestep = " << step << std::endl;
						tetMesh.exportSurfaceMesh("surfMesh", step);
					}


					tetMesh.elastForce_node.resize(tetMesh.pos_node.size());
					// Reset forces
					for (int nf = 0; nf < tetMesh.elastForce_node.size(); nf++) 
					{
						tetMesh.elastForce_node[nf] = Eigen::Vector3d::Zero();
					}
					tetMesh.elastForce_node[0] = {0,135,100};
					tetMesh.elastForce_node[1] = {145,47,28};
					tetMesh.elastForce_node[2] = {13,167,95};
					tetMesh.elastForce_node[3] = {0,0,100};

					// Apply gravity
					for (int nf = 0; nf < tetMesh.elastForce_node.size(); nf++)
					{
						tetMesh.elastForce_node[nf] += tetMesh.mass_node[nf] * parameters.gravity;
					}

					// Compute forces from elements
					for (int nf = 0; nf < tetMesh.tetrahedrals.size(); nf++) 
					{
						// Compute deformation gradient F
						tetMesh.update_F(1);


						// Compute First Piola-Kirchhoff stress tensor P
						Eigen::JacobiSVD<Eigen::MatrixXd> svd(tetMesh.tetra_F[nf], Eigen::ComputeFullU | Eigen::ComputeFullV);
						Eigen::MatrixXd U = svd.matrixU();
						Eigen::VectorXd S = svd.singularValues();
						Eigen::MatrixXd V = svd.matrixV();
						Eigen::MatrixXd R = U * V.transpose();
						Eigen::MatrixXd Q = V * S.asDiagonal() * V.transpose();
						
						Eigen::Matrix3d P = R;


						if (step % 1000 == 0)
						{
							std::cout << "Rotation Matrix R:\n" << R << "\n\n";
							std::cout << "Stretch Matrix S:\n" << Q << "\n\n";
						}


						// Compute forces on nodes
						Eigen::Matrix3d H = -tetMesh.tetra_vol[nf] * P * tetMesh.tetra_DM_inv[nf].transpose();

						// Distribute forces to nodes
						Eigen::Vector3d f0 = H.col(0);
						Eigen::Vector3d f1 = H.col(1);
						Eigen::Vector3d f2 = H.col(2);
						Eigen::Vector3d f3 = -f0 - f1 - f2;


						Eigen::Vector4i tet = tetMesh.tetrahedrals[nf];
						tetMesh.elastForce_node[tet[0]] += f3;
						tetMesh.elastForce_node[tet[1]] += f0;
						tetMesh.elastForce_node[tet[2]] += f1;
						tetMesh.elastForce_node[tet[3]] += f2;
					}

					// Time integration (Explicit Euler)
					for (int nf = 0; nf < tetMesh.elastForce_node.size(); nf++) 
					{
						// Update velocity and position
						Eigen::Vector3d acceleration = tetMesh.elastForce_node[nf] / tetMesh.mass_node[nf];
						tetMesh.vel_node[nf] += acceleration * parameters.dt;
						tetMesh.pos_node[nf] += tetMesh.vel_node[nf] * parameters.dt;
					}
					
				}


			}



		}
		else if (caseNum == 3)
		{
			Material mat1;
			mat1.density = 2780;
			mat1.E = 7.26e7;
			mat1.updateDenpendecies();



			std::vector<meshConfiguration> config;
			meshConfiguration m1, m2, m3, m4;
			m1.filePath = "../input/beam.msh";
			m1.mesh_material = mat1;
			m1.scale = { 1000, 1000, 1000 };
			m1.note = "beam";
			config.push_back(m1);

			m4 = m1;
			m4.filePath = "../input/impactor.msh";
			m4.note = "impactor";
			m4.translation = { 0, 0, 1.5e-2 };
			config.push_back(m4);




			Mesh tetMesh;
			std::vector<std::string> objectNames;
			for (int i = 0; i < config.size(); i++)
			{
				objectNames.push_back(config[i].note);
			}
			tetMesh.readMeshes(config);
			tetMesh.initializeMesh();
			tetMesh.surfaceMesh.outputFile("surfMesh");


			std::cout << "tetMesh.pos_node.size() = " << tetMesh.pos_node.size() << std::endl;
			std::cout << "tetMesh.tetrahedrals.size() = " << tetMesh.tetrahedrals.size() << std::endl;


			FEMParamters parameters;
			parameters.gravity = { 0, 0, 0 };
			parameters.num_timesteps = 3000;
			parameters.numOfThreads = 12;
			parameters.dt = 5.0e-4;
			parameters.outputFrequency = 20;
			parameters.simulation_Mode = "ABD"; // Normal, ABD, Coupling
			parameters.enableGround = false;
			parameters.searchResidual = 5.0;
			parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
			parameters.rigidMode = true;
			parameters.objectNames = objectNames;
			parameters.IPC_dis = 0.01;
			parameters.IPC_eta = 0.05;
			parameters.IPC_kStiffness = 1.0e14;
			parameters.IPC_hashSize = tetMesh.calLargestEdgeLength() * 1.1;
			parameters.IPC_B3Stiffness = 500;



			for (int p = 0; p < tetMesh.pos_node.size(); p++)
			{
				if (tetMesh.note_node[p] == "impactor")
				{
					tetMesh.boundaryCondition_node[p].type = 1;
					for (int fra = 0; fra < parameters.num_timesteps; fra++)
					{
						double incre = 0.1 / (double)parameters.num_timesteps * (double)fra;
						Eigen::Vector3d inc = { 0,0,-incre };
						Eigen::Vector3d desiredPos = inc + tetMesh.pos_node[p];

						tetMesh.boundaryCondition_node[p].location.push_back(desiredPos);
					}
				}

				if (tetMesh.note_node[p] == "beam")
				{
					if (tetMesh.pos_node[p][0] <= -22)
					{
						tetMesh.boundaryCondition_node[p].type = 1;
						for (int fra = 0; fra < parameters.num_timesteps; fra++)
						{
							Eigen::Vector3d desiredPos = tetMesh.pos_node[p];
							tetMesh.boundaryCondition_node[p].location.push_back(desiredPos);
						}
					}
				}

			}



			implicitFEM(tetMesh, parameters);


		}

	}


	return 0;

}
