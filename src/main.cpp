#include "simulator.h" 
#include "tools.h" 

// TODO
// 1. Accelerate and parallel openMP
// 2. If the code doesn't work, check line45 of InertiaEnergy.cpp


int main()
{

	if (0)
	{
		FEMParamters parameters;
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
		outfile9.close();




		/*double pi = 3.141592653;

		//Eigen::Matrix3d F = Eigen::Matrix3d::Identity();

		Eigen::Matrix3d RX = Eigen::Matrix3d::Zero();
		Eigen::Matrix3d RY = Eigen::Matrix3d::Zero();
		Eigen::Matrix3d RZ = Eigen::Matrix3d::Zero();

		double tx = 0, ty = 0, tz = 0;
		RX(0, 0) = 1.0;
		RX(1, 1) = cos(tx);
		RX(1, 2) = -sin(tx);
		RX(2, 1) = sin(tx);
		RX(2, 2) = cos(tx);

		F = F * RX;

		double xv2 = F(0,0) * PI / 2.0 * PI + F(0, 1) * PI / 2.0 * PI + F(0, 2) * PI * PI + 2 * F(0, 0) * F(0, 2) * 2.0 * PI ;
		double yv2 = F(1,0) * PI / 2.0 * PI + F(1, 1) * PI / 2.0 * PI + F(1, 2) * PI * PI + 2 * F(1, 0) * F(1, 2) * 2.0 * PI ;
		double zv2 = F(2,0) * PI / 2.0 * PI + F(2, 1) * PI / 2.0 * PI + F(2, 2) * PI * PI + 2 * F(2, 0) * F(2, 2) * 2.0 * PI ;

		std::cout << "(xv2 + yv2 + zv2) / 2 / PI /PI = " << (xv2 + yv2 + zv2) / 2 / PI / PI << std::endl;*/




	}
	else
	{
		int caseNum = 0;
		if (caseNum == 0)
		{
			Material mat1;
			mat1.density = 2780;
			mat1.E = 7.26e10;
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
			config.push_back(m2);

			m3 = m1;
			m3.filePath = "../input/middle_support.msh";
			m3.note = "middle_support";
			m3.translation = { 0, 0, -1.0e-4 };
			config.push_back(m3);

			m4 = m1;
			m4.filePath = "../input/impactor.msh";
			m4.note = "impactor";
			m4.translation = { 0, 0, 1.0e-4 };
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
			parameters.numOfThreads = 32;
			parameters.dt = 5.0e-5;
			parameters.outputFrequency = 20;
			parameters.enableGround = false;
			parameters.searchResidual = 5.0;
			parameters.model = "ARAP_linear"; // neoHookean ARAP ARAP_linear ACAP
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
						double incre = 2.0 / (double)parameters.num_timesteps * (double)fra;
						Eigen::Vector3d inc = {0,0,incre };
						Eigen::Vector3d desiredPos = inc + tetMesh.pos_node[p];

						tetMesh.boundaryCondition_node[p].location.push_back(desiredPos);
					}
				}
			}



			implicitFEM(tetMesh, parameters);


		}
		

	}



	return 0;
}
