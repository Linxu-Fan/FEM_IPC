#include "simulator.h" 


int main()
{
	Material mat1;
	//mat1.density = 100;
	//mat1.E = 100;
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
	parameters.gravity = {0, 0, -9.8};
	parameters.num_timesteps = 1001;
	parameters.outputFrequency = 100;
	parameters.model = "neoHookean";
	// neoHookean ARAP ARAP_linear ACAP

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
