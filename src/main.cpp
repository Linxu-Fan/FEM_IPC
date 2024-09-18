#include "simulator.h" 

// TODO
// 1. Accelerate and parallel openMP
// 2. If the code doesn't work, check line45 of InertiaEnergy.cpp


int main()
{

	if (1)
	{
		double pi = 3.141592653;

		Eigen::Matrix3d F = Eigen::Matrix3d::Identity();

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

		std::cout << "(xv2 + yv2 + zv2) / 2 / PI /PI = " << (xv2 + yv2 + zv2) / 2 / PI / PI << std::endl;




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


			//// cube collision test
			//std::vector<meshConfiguration> config;
			//meshConfiguration m1, m2;
			//m1.filePath = "../input/cube.msh";
			//m1.mesh_material = mat1;
			//m1.note = "cube1";
			//m1.velocity = { 0, 0 , 0 };
			////m1.translation = { 0, 0, 0.1 };
			////config.push_back(m1);
			//m2 = m1;
			//m2.mesh_material = mat1;
			//m2.note = "cube2";
			//m2.velocity = { 0, 0 , 0 };
			//m2.rotation_point = { 0.5, 0.5, 0.5 };
			//m2.rotation_angle = { PI / 4.0, PI / 4.0, PI / 4.0 };
			//m2.translation = { 1.80, 0.5, -0.3 };
			////config.push_back(m2);



			//// tet drop test
			//std::vector<meshConfiguration> config;
			//meshConfiguration m1, m2;
			//m1.filePath = "../input/tet.msh";
			//m1.mesh_material = mat1;
			//m1.note = "tet1";
			////m1.velocity = { 0, 0 , -2.0 };
			////m1.translation = { 0, 0, 0.1 };
			//config.push_back(m1);
			//m2 = m1;
			//m2.mesh_material = mat1;
			//m2.velocity = { 0, 0 , -2.0 };
			//m2.rotation_point = { 0, 0, 0 };
			////m2.translation = { 1.03, 0, 1.0 };
			//m2.rotation_angle = { PI , 0, 0 };
			//m2.translation = { 0, 0, 1.034 };
			////config.push_back(m2);



			std::vector<meshConfiguration> config;			
			meshConfiguration m1, m2, m3, m4, m5;
			m1.filePath = "./input/beam.msh";
			m1.mesh_material = mat1;
			m1.scale = {1000, 1000, 1000};
			m1.velocity = { 0, 0 , 0 };
			m1.note = "beam";
			config.push_back(m1);

			m2 = m1;
			m2.filePath = "../input/Left_bottom_fix.msh";
			m2.note = "Left_bottom_fix";
			//config.push_back(m2);

			m3 = m1;
			m3.filePath = "../input/Left_top_fix.msh";
			m3.note = "Left_top_fix";
			m3.translation = { 0, 0, -0.1 };
			m3.velocity = { 0, 0 , -2 };
			//config.push_back(m3);

			m4 = m1;
			m4.filePath = "../input/Middle_support.msh";
			m4.note = "Middle_support";
			m4.translation = { 0, 0, 0.197 };
			//config.push_back(m4);

			m5 = m1;
			m5.filePath = "../input/Right_top_move.msh";
			m5.translation = { 0, 0, -0.1 };
			m5.velocity = { 0, 0 , -2 };
			m5.note = "Right_top_move";
			//config.push_back(m5);


			Mesh tetMesh;
			std::vector<std::string> objectNames;
			for (int i = 0; i < config.size(); i++)
			{
				objectNames.push_back(config[i].note);
			}
			tetMesh.readMeshes(config);
			tetMesh.initializeMesh();
			tetMesh.surfaceMesh.outputFile("surfMesh");

			//for (int p = 0; p < tetMesh.pos_node.size(); p++)
			//{
			//	if (tetMesh.pos_node[p][0] < 0.3)
			//	{
			//		tetMesh.vel_node[p] = { 5.0 , 0 , 0 };

			//		tetMesh.boundaryCondition_node[p].type = 3;
			//		tetMesh.boundaryCondition_node[p].velocity = {5.0 , 0 , 0};
			//	}
			//}



			for (int p = 0; p < tetMesh.pos_node.size(); p++)
			{
				if (tetMesh.pos_node[p][0] <= 2)
				{
					tetMesh.boundaryCondition_node[p].type = 1;
				}

				if (tetMesh.note_node[p] == "Middle_support")
				{
					tetMesh.boundaryCondition_node[p].type = 1;
				}
			}


			//for (int p = 0; p < tetMesh.pos_node.size(); p++)
			//{
			//	if (tetMesh.note_node[p] == "tet1" )
			//	{
			//		tetMesh.vel_node[p] = { 0, 0 , -2 };
			//		tetMesh.boundaryCondition_node[p].type = 3;
			//		tetMesh.boundaryCondition_node[p].velocity = { 0, 0 , -2 };
			//	}
			//}




			std::cout << "tetMesh.pos_node.size() = " << tetMesh.pos_node.size() << std::endl;
			std::cout << "tetMesh.tetrahedrals.size() = " << tetMesh.tetrahedrals.size() << std::endl;


			FEMParamters parameters;
			parameters.gravity = { 0, 0, 0 };
			parameters.num_timesteps = 2000000;
			parameters.numOfThreads = 32;
			parameters.dt = 1.0e-4;
			parameters.outputFrequency = 20;
			parameters.enableGround = false;
			parameters.searchResidual = 0.001;
			parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
			parameters.rigidMode = true;
			parameters.objectNames = objectNames;
			parameters.IPC_dis = 0.001;
			parameters.IPC_eta = 0.05;
			parameters.IPC_kStiffness = 1.0e16;
			parameters.IPC_hashSize = tetMesh.calLargestEdgeLength() * 1.1;

			std::cout << "parameters.IPC_hashSize = " << parameters.IPC_hashSize << std::endl;

			implicitFEM(tetMesh, parameters);


		}
		

	}



	return 0;
}
