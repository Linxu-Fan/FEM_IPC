#include "simulator.h" 

// TODO
// 1. Accelerate and parallel openMP



int main()
{

	if (0)
	{
		Material mat1;
		mat1.updateDenpendecies();

		Material mat2;
		mat2.updateDenpendecies();



		std::vector<meshConfiguration> config;
		meshConfiguration m1, m2;
		m1.filePath = "./input/tet_neg.msh";
		m1.mesh_material = mat1;
		m1.velocity = { 1, 0 , 0 };
		m1.translation = { 0, 0, 0.1 };
		config.push_back(m1);

		m2 = m1;
		m2.mesh_material = mat2;
		m2.velocity = { -1, 0 , 0 };
		m2.translation = { 1.1, 0, 0.1 };
		config.push_back(m2);
		Mesh tetMesh;
		tetMesh.readMeshes(config);
		tetMesh.initializeMesh();
		tetMesh.surfaceMesh.outputFile("surfMesh");



		FEMParamters parameters;
		parameters.dt = 1.0E-3;
		parameters.gravity = { 0, 0, 0 };
		parameters.num_timesteps = 1511;
		parameters.outputFrequency = 1;
		parameters.enableGround = true;
		parameters.searchResidual = 0.003;
		parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
		parameters.IPC_dis = 0.001;
		parameters.IPC_eta = 0.1;
		parameters.IPC_kStiffness = 1.0e14;
		parameters.IPC_hashSize = tetMesh.calLargestEdgeLength() * 40.0;




		implicitFEM(tetMesh, parameters);
	}
	else
	{

		Material mat1;
		mat1.density = 1000;
		mat1.E = 1.0e6;
		mat1.updateDenpendecies();


		//// cube collision test
		//std::vector<meshConfiguration> config;
		//meshConfiguration m1, m2;
		//m1.filePath = "./input/cube.msh";
		//m1.mesh_material = mat1;
		//m1.velocity = { 5.1245678, 0 , 0 };
		//m1.translation = { 0, 0, 1.0 };
		//config.push_back(m1);
		//m2 = m1;
		//m2.mesh_material = mat1;
		//m2.velocity = { -5.1245678, 0 , 0 };
		//m2.rotation_point = { 0.5, 0.5, 0.5 };
		////m2.translation = { 1.03, 0, 1.0 };
		//m2.rotation_angle = { PI / 4.0, PI / 4.0, PI / 4.0 };
		//m2.translation = { 1.66, 0.5, 0.6 };
		//config.push_back(m2);



		std::vector<meshConfiguration> config;
		meshConfiguration m1, m2, m3, m4, m5;
		m1.filePath = "./input/beam.msh";
		m1.mesh_material = mat1;
		m1.scale = {1000, 1000, 1000};
		m1.note = "beam";
		config.push_back(m1);

		m2 = m1;
		m2.filePath = "./input/Left_bottom_fix.msh";
		m2.note = "Left_bottom_fix";
		config.push_back(m2);

		m3 = m1;
		m3.filePath = "./input/Left_top_fix.msh";
		m3.note = "Left_top_fix";
		config.push_back(m3);

		m4 = m1;
		m4.filePath = "./input/Middle_support.msh";
		m4.note = "Middle_support";
		config.push_back(m4);

		m5 = m1;
		m5.filePath = "./input/Right_top_move.msh";
		m5.note = "Right_top_move";
		config.push_back(m5);


		Mesh tetMesh;
		tetMesh.readMeshes(config);
		tetMesh.initializeMesh();
		tetMesh.surfaceMesh.outputFile("surfMesh");


		std::cout << "tetMesh.pos_node.size() = " << tetMesh.pos_node.size() << std::endl;
		std::cout << "tetMesh.tetrahedrals.size() = " << tetMesh.tetrahedrals.size() << std::endl;


		FEMParamters parameters;
		parameters.dt = 2.0E-4;
		parameters.gravity = { 0, 0, 0 };
		parameters.num_timesteps = 10001;
		parameters.outputFrequency = 10;
		parameters.enableGround = true;
		parameters.searchResidual = 0.001;
		parameters.model = "ARAP_linear"; // neoHookean ARAP ARAP_linear ACAP
		parameters.IPC_dis = 0.000001;
		parameters.IPC_eta = 0.05;
		parameters.IPC_kStiffness = 1.0e15;
		parameters.IPC_hashSize = tetMesh.calLargestEdgeLength() * 1.1;

		

		std::vector<Eigen::Vector3d> direction(tetMesh.pos_node.size());
		std::vector<spatialHashCellData> spatialHash_vec;


		double startTime, endTime;
		startTime = omp_get_wtime();
		initSpatialHash(parameters, tetMesh, direction, parameters.IPC_hashSize, spatialHash_vec);
		std::cout << "spatialHash = " << spatialHash_vec.size() << std::endl;
		endTime = omp_get_wtime();
		std::cout << "Hash time = " << endTime - startTime << std::endl;

		double startTime2, endTime2;
		startTime2 = omp_get_wtime();
		double step = calMaxStep_spatialHash(parameters, tetMesh, direction, parameters.IPC_hashSize, parameters.IPC_dis, parameters.IPC_eta);
		endTime2 = omp_get_wtime();
		std::cout << "Stepsize time = " << endTime2 - startTime2 << std::endl;

		//std::cout << "step = " << step << std::endl;

		implicitFEM(tetMesh, parameters);



	}



	return 0;
}
