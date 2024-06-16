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
	parameters.dt = 1.0E-3;
	parameters.gravity = {0, 0, 0};
	parameters.num_timesteps = 10001;
	parameters.outputFrequency = 1000;
	parameters.model = "ARAP_linear";
	// neoHookean ARAP ARAP_linear ACAP




	implicitFEM(tetMesh, parameters);




	return 0;
}
