#include "simulator.h" 

// TODO
// 1. Compute the barrier distance, energy and their derivative
// 2. Accelerate and parallel openMP



int main()
{

	if (0)
	{
		std::cout << "PP" << std::endl;
		DIS::derivTest_PP();
		std::cout << std::endl << "PE" << std::endl;
		DIS::derivTest_PE();
		std::cout << std::endl << "PT" << std::endl;
		DIS::derivTest_PT();
		std::cout << std::endl << "EE" << std::endl;
		DIS::derivTest_EE();
		std::cout << std::endl << "DType" << std::endl;
		DIS::checkDType();
	}
	else
	{

		Material mat1;
		mat1.updateDenpendecies();

		Material mat2;
		mat2.density = 100;
		mat2.updateDenpendecies();



		//meshConfiguration m1;
		//m1.filePath = "./input/cube.msh";
		//m1.mesh_material = mat1;
		//m1.velocity = { 0, 0 , 0 };
		//Mesh tetMesh;
		//tetMesh.readMesh(m1);
		//tetMesh.initializeMesh();
		//tetMesh.surfaceMesh.outputFile("surfMesh");


		std::vector<meshConfiguration> config;
		meshConfiguration m1, m2;
		m1.filePath = "./input/cube.msh";
		m1.mesh_material = mat1;
		m1.velocity = { 10, 0 , 0 };
		config.push_back(m1);
		m2 = m1;
		m2.mesh_material = mat2;
		m2.velocity = { -10, 0 , 0 };
		m2.shift = { 2.157, 0.3458, -0.235 };
		config.push_back(m2);
		Mesh tetMesh;
		tetMesh.readMeshes(config);
		tetMesh.initializeMesh();
		tetMesh.surfaceMesh.outputFile("surfMesh");



		FEMParamters parameters;
		parameters.dt = 1.0E-3;
		parameters.gravity = { 0, 0, 0 };
		parameters.num_timesteps = 1001;
		parameters.outputFrequency = 1;
		parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
		parameters.IPC_dis = 0.001;
		parameters.IPC_eta = 0.1;
		parameters.IPC_hashSize = tetMesh.calLargestEdgeLength() * 2.0;



		//for (int vI = 0; vI < tetMesh.pos_node.size(); vI++)
		//{
		//	if (tetMesh.pos_node[vI][0] <= 0.15)
		//	{
		//		tetMesh.boundaryCondition_node[vI].type = 2;
		//		tetMesh.boundaryCondition_node[vI].force = {-100, 0, 0};
		//	}
		//	if (tetMesh.pos_node[vI][0] >= 0.85)
		//	{
		//		tetMesh.boundaryCondition_node[vI].type = 2;
		//		tetMesh.boundaryCondition_node[vI].force = { 100, 0, 0 };
		//	}
		//}





		implicitFEM(tetMesh, parameters);



	}



	return 0;
}
