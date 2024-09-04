#include "simulator.h" 

// TODO
// 1. Accelerate and parallel openMP
// 2. If the code doesn't work, check line45 of InertiaEnergy.cpp


int main()
{

	if (0)
	{
		Eigen::Vector2i res = findIntersectionOfTwoNums(0, 3, 2, 5);
		std::cout << res;

	}
	else
	{
		int caseNum = 0;
		if (caseNum == 0)
		{
			Material mat1;
			mat1.density = 1000;
			mat1.E = 1.0e6;
			mat1.updateDenpendecies();


			// cube collision test
			std::vector<meshConfiguration> config;
			meshConfiguration m1, m2;
			m1.filePath = "../input/cube.msh";
			m1.mesh_material = mat1;
			m1.note = "cube1";
			m1.velocity = { 0, 0 , 0 };
			//m1.translation = { 0, 0, 0.1 };
			config.push_back(m1);
			m2 = m1;
			m2.mesh_material = mat1;
			m2.note = "cube2";
			m2.velocity = { 0, 0 , 0 };
			m2.rotation_point = { 0.5, 0.5, 0.5 };
			m2.rotation_angle = { PI / 4.0, PI / 4.0, PI / 4.0 };
			m2.translation = { 1.80, 0.5, -0.3 };
			config.push_back(m2);



			//// tet drop test
			//std::vector<meshConfiguration> config;
			//meshConfiguration m1, m2;
			//m1.filePath = "../input/tet.msh";
			//m1.mesh_material = mat1;
			//m1.velocity = { 0, 0 , -2.0 };
			//m1.translation = { 0, 0, 0.1 };
			////config.push_back(m1);
			//m2 = m1;
			//m2.mesh_material = mat1;
			//m2.velocity = { 0, 0 , -2.0 };
			//m2.rotation_point = { 0, 0, 0 };
			////m2.translation = { 1.03, 0, 1.0 };
			//m2.rotation_angle = { PI , 0, 0 };
			//m2.translation = { 0, 0, 1.034 };
			//config.push_back(m2);



			//std::vector<meshConfiguration> config;			
			//meshConfiguration m1, m2, m3, m4, m5;
			//m1.filePath = "./input/beam.msh";
			//m1.mesh_material = mat1;
			//m1.scale = {1000, 1000, 1000};
			//m1.velocity = { 0, 0 , 0 };
			//m1.note = "beam";
			//config.push_back(m1);

			//m2 = m1;
			//m2.filePath = "../input/Left_bottom_fix.msh";
			//m2.note = "Left_bottom_fix";
			//config.push_back(m2);

			//m3 = m1;
			//m3.filePath = "../input/Left_top_fix.msh";
			//m3.note = "Left_top_fix";
			//config.push_back(m3);

			//m4 = m1;
			//m4.filePath = "../input/Middle_support.msh";
			//m4.note = "Middle_support";
			//config.push_back(m4);

			//m5 = m1;
			//m5.filePath = "../input/Right_top_move.msh";
			//m5.translation = { 0, 0, -0.1 };
			//m5.velocity = { 0, 0 , -2 };
			//m5.note = "Right_top_move";
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

			for (int p = 0; p < tetMesh.pos_node.size(); p++)
			{
				if (tetMesh.pos_node[p][0] < 0.3)
				{
					tetMesh.vel_node[p] = { 5.0 , 0 , 0 };

					tetMesh.boundaryCondition_node[p].type = 3;
					tetMesh.boundaryCondition_node[p].velocity = {5.0 , 0 , 0};
				}
			}


			std::cout << "tetMesh.pos_node.size() = " << tetMesh.pos_node.size() << std::endl;
			std::cout << "tetMesh.tetrahedrals.size() = " << tetMesh.tetrahedrals.size() << std::endl;


			FEMParamters parameters;
			parameters.gravity = { 0, 0, 0 };
			parameters.num_timesteps = 200000;
			parameters.dt = 1.0e-3;
			parameters.outputFrequency = 1;
			parameters.enableGround = false;
			parameters.searchResidual = 0.001;
			parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
			parameters.rigidMode = true;
			parameters.objectNames = objectNames;
			parameters.IPC_dis = 0.001;
			parameters.IPC_eta = 0.05;
			parameters.IPC_kStiffness = 1.0e15;
			parameters.IPC_hashSize = tetMesh.calLargestEdgeLength() * 1.1;

			std::cout << "parameters.IPC_hashSize = " << parameters.IPC_hashSize << std::endl;

			implicitFEM(tetMesh, parameters);


		}
		

	}



	return 0;
}
