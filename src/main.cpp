#include "simulator.h" 


int main()
{
	Material mat1;
	mat1.updateDenpendecies();

	Mesh tetMesh;
	tetMesh.readMesh("./input/cube.msh");
	tetMesh.materialMesh = mat1;
	tetMesh.initializeMesh();
	tetMesh.cal_DS_or_DM(false);
	tetMesh.extractSurfaceMesh();
	tetMesh.surfaceMesh.outputFile("surfMesh");


	FEMParamters parameters;
	parameters.dt = 1.0E-2;
	parameters.gravity = {0, 0, 0};
	parameters.num_timesteps = 2001;
	parameters.outputFrequency = 100;
	parameters.model = "neoHookean";
	// neoHookean ARAP ARAP_linear ACAP



	std::cout << tetMesh.pos_node.size() << std::endl;
	std::cout << tetMesh.tetrahedrals.size() << std::endl;

	//
	//implicitFEM(tetMesh, parameters);




	return 0;
}
