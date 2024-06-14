#include "FEM.h"






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





	std::ofstream outfile9("./output/left.obj", std::ios::trunc);
	std::ofstream outfile10("./output/right.obj", std::ios::trunc);
	tetMesh.boundaryCondition_node.resize(tetMesh.pos_node.size());
	//for (int k = 0; k < tetMesh.pos_node.size(); k++)
	//{
	//	if (tetMesh.pos_node[k][2] <= 10)
	//	{
	//		tetMesh.boundaryCondition_node[k].type = 1;
	//		outfile9 << std::scientific << std::setprecision(8) << "v " << tetMesh.pos_node[k][0] << " " << tetMesh.pos_node[k][1] << " " << tetMesh.pos_node[k][2] << std::endl;
	//	}


	//	if (tetMesh.pos_node[k][2] >= 125)
	//	{
	//		tetMesh.boundaryCondition_node[k].type = 2;
	//		tetMesh.boundaryCondition_node[k].forceMagnitude = {0, -1000, 0 };

	//		outfile10 << std::scientific << std::setprecision(8) << "v " << tetMesh.pos_node[k][0] << " " << tetMesh.pos_node[k][1] << " " << tetMesh.pos_node[k][2] << std::endl;
	//	}

	//}

	for (int k = 0; k < tetMesh.pos_node.size(); k++)
	{
		if (tetMesh.pos_node[k][0] <= 0.2)
		{
			tetMesh.boundaryCondition_node[k].type = 2;
			tetMesh.boundaryCondition_node[k].forceMagnitude = { -100, 0, 0 };
		}


		if (tetMesh.pos_node[k][0] >= 0.8)
		{
			tetMesh.boundaryCondition_node[k].type = 2;
			tetMesh.boundaryCondition_node[k].forceMagnitude = { 100, 0, 0 };
			
		}

	}
	outfile9.close();
	outfile10.close();

	FEMParamters parameters;
	parameters.dt = 1.0E-2;
	parameters.gravity = {0, -9.8, 0};
	parameters.num_timesteps = 201;
	parameters.outputFrequency = 10;
	parameters.model = "ACAP";
	// neoHookean ARAP ARAP_linear ACAP

	
	implicitFEM(tetMesh, mat1, parameters);



	std::cout << tetMesh.pos_node.size() << std::endl;
	std::cout << tetMesh.tetrahedrals.size() << std::endl;
	return 0;
}
