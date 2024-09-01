#include "simulator.h" 

// TODO
// 1. Accelerate and parallel openMP
// 2. If the code doesn't work, check line45 of InertiaEnergy.cpp


int main()
{

	if (0)
	{
		std::vector<double> test;
		test.push_back(2.0);
		test.push_back(4.0);
		test.push_back(6.0);
		test.push_back(1.0);
		double value = std::accumulate(test.begin(), test.end(), 0.0);
		std::cout << "value = " << value << std::endl;

		//std::vector<Eigen::Triplet<double>> hessian_triplet;
		//hessian_triplet.emplace_back(1, 1, 2);
		//hessian_triplet.emplace_back(2, 3, 5);
		//hessian_triplet.emplace_back(3, 1, 4);
		//hessian_triplet.emplace_back(2, 2, 4);

		//std::vector<Eigen::Triplet<double>> hessian_triplet2;
		//hessian_triplet2.emplace_back(1, 1, 3);
		//hessian_triplet2.emplace_back(2, 3, 4);
		//hessian_triplet2.emplace_back(3, 1, 1);
		//hessian_triplet2.emplace_back(2, 2, 2);

		//Eigen::SparseMatrix<double> leftHandSide(4, 4), leftHandSide2(4, 4);
		//leftHandSide.setFromTriplets(hessian_triplet.begin(), hessian_triplet.end());
		//leftHandSide2.setFromTriplets(hessian_triplet2.begin(), hessian_triplet2.end());
		//leftHandSide += leftHandSide2;

		//std::cout << "leftHandSide = " << leftHandSide << std::endl;

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


			//// cube collision test
			//std::vector<meshConfiguration> config;
			//meshConfiguration m1, m2;
			//m1.filePath = "../input/cube.msh";
			//m1.mesh_material = mat1;
			//m1.velocity = { 5.0, 0 , 0 };
			//m1.translation = { 0, 0, 0.1 };
			//config.push_back(m1);
			//m2 = m1;
			//m2.mesh_material = mat1;
			//m2.velocity = { -5.0, 0 , 0 };
			//m2.rotation_point = { 0.5, 0.5, 0.5 };
			////m2.translation = { 1.03, 0, 1.0 };
			//m2.rotation_angle = { PI / 4.0, PI / 4.0, PI / 4.0 };
			//m2.translation = { 1.76, 0.5, -0.1 };
			//config.push_back(m2);



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




			std::vector<meshConfiguration> config;
			meshConfiguration m1, m2, m3, m4, m5;
			m1.filePath = "./input/beam.msh";
			m1.mesh_material = mat1;
			//m1.scale = {1000, 1000, 1000};
			m1.note = "beam";
			config.push_back(m1);

			m2 = m1;
			m2.filePath = "../input/Left_bottom_fix.msh";
			m2.note = "Left_bottom_fix";
			config.push_back(m2);

			m3 = m1;
			m3.filePath = "../input/Left_top_fix.msh";
			m3.note = "Left_top_fix";
			config.push_back(m3);

			m4 = m1;
			m4.filePath = "../input/Middle_support.msh";
			m4.note = "Middle_support";
			config.push_back(m4);

			m5 = m1;
			m5.filePath = "../input/Right_top_move.msh";
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
			parameters.enableGround = false;
			parameters.searchResidual = 0.001;
			parameters.model = "ARAP_linear"; // neoHookean ARAP ARAP_linear ACAP
			parameters.IPC_dis = 0.001;
			parameters.IPC_eta = 0.05;
			parameters.IPC_kStiffness = 1.0e15;
			parameters.IPC_hashSize = tetMesh.calLargestEdgeLength() * 1.3;



			implicitFEM(tetMesh, parameters);


		}
		

	}



	return 0;
}
