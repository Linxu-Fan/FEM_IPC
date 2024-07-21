#include "simulator.h" 

// TODO
// 1. Compute the barrier distance, energy and their derivative
// 2. Accelerate and parallel openMP



int main()
{
	Material mat1;
	mat1.updateDenpendecies();

	Material mat2;
	mat2.density = 100;
	mat2.updateDenpendecies();



	//Mesh tetMesh;
	//tetMesh.readMesh("./input/cube.msh", mat1);
	//tetMesh.initializeMesh();
	//tetMesh.extractSurfaceMesh();
	//tetMesh.surfaceMesh.outputFile("surfMesh");


	std::vector<std::pair<std::string, Material>> filePath_and_mat;
	filePath_and_mat.push_back(std::make_pair("./input/cube.msh", mat1));
	filePath_and_mat.push_back(std::make_pair("./input/cube.msh", mat2));
	Mesh tetMesh;
	tetMesh.readMeshes(filePath_and_mat);
	tetMesh.initializeMesh();
	tetMesh.surfaceMesh.outputFile("surfMesh");



	FEMParamters parameters;
	parameters.dt = 1.0E-2;
	parameters.gravity = {0, 0, -9.8};
	parameters.num_timesteps = 1001;
	parameters.outputFrequency = 100;
	parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
	parameters.IPC_dis = 0.001;
	parameters.IPC_eta = 0.1;
	parameters.IPC_hashSize = tetMesh.calLargestEdgeLength();

	

	for (int vI = 0; vI < tetMesh.pos_node.size(); vI++)
	{
		if (tetMesh.pos_node[vI][0] <= 0.15)
		{
			tetMesh.boundaryCondition_node[vI].type = 1;
		}
		if (tetMesh.pos_node[vI][0] >= 0.85)
		{
			tetMesh.boundaryCondition_node[vI].type = 1;
		}
	}





	implicitFEM(tetMesh, parameters);




	return 0;
}
