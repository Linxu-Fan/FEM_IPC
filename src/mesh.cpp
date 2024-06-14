#include "mesh.h"


void objMesh::clear()
{
	vertices.clear();
	faces.clear();
}


void objMesh::outputFile(std::string fileName, int timestep)
{
	std::ofstream outfile9("./output/" + fileName + "_" + std::to_string(timestep) + ".obj", std::ios::trunc);
	for (int k = 0; k < vertices.size(); k++)
	{
		Eigen::Vector3d scale = vertices[k];
		outfile9 << std::scientific << std::setprecision(8) << "v " << scale[0] << " " << scale[1] << " " << scale[2] << std::endl;
	}
	for (int k = 0; k < faces.size(); k++)
	{
		outfile9 << "f ";
		for (int m = 0; m < faces[k].size(); m++)
		{
			outfile9 << faces[k][m] + 1 << " ";
		}
		outfile9 << std::endl;
	}
	outfile9.close();
}


////////////////////////////////////////////////////////////////////////
// Is it possible that node and element are not placed in order? If possible, then the reading code may crash.
////////////////////////////////////////////////////////////////////////
// read .msh mesh file
void Mesh::readMesh(std::string filePath)
{
	std::ifstream in;
	in.open(filePath);
	std::string line;
	int nodeStart = 100000000000; // the line index where a node starts
	int numNodes = 100000000000; // number of nodes
	int elementStart = 100000000000; // the element index where an element starts
	int numElements = 100000000000; // number of elements
	int lineIndex = -1; // current line index

	while (getline(in, line))
	{
		if (line.size() > 0)
		{
			lineIndex += 1;
			
			std::vector<std::string> vecCoor = split(line, " ");
			if (vecCoor[0] == "$Nodes")
			{
				nodeStart = lineIndex + 2;
			}
			if (lineIndex == nodeStart - 1)
			{
				numNodes = std::stoi(vecCoor[0]);
			}

			if (vecCoor[0] == "$Elements")
			{
				elementStart = lineIndex + 2;
			}
			if (lineIndex == elementStart - 1)
			{
				numElements = std::stoi(vecCoor[0]);				
			}


			if (lineIndex >= nodeStart && lineIndex <= nodeStart + numNodes - 1)
			{
				Eigen::Vector3d nd_pos = { std::stod(vecCoor[1]) , std::stod(vecCoor[2]) , std::stod(vecCoor[3])};
				pos_node.push_back(nd_pos);
			}

			if (lineIndex >= elementStart && lineIndex <= elementStart + numElements - 1)
			{
				int numItemsLine = vecCoor.size(); // the number of items in a line
				if (vecCoor[1] == "4")
				{
					Eigen::Vector4i ele = { std::stoi(vecCoor[numItemsLine - 4]) - 1 ,std::stoi(vecCoor[numItemsLine - 3]) - 1 ,std::stoi(vecCoor[numItemsLine - 2]) - 1,std::stoi(vecCoor[numItemsLine - 1]) - 1 };
					tetrahedrals.push_back(ele);
				}
			}
		}
	}
	in.close();



}

// initialize the mesh after reading
void Mesh::initializeMesh() // initialize the mesh 
{
	for (int i = 0; i < pos_node.size(); i++)
	{
		Eigen::Vector3d vel = { 0,0,0 };
		vel_node.push_back(vel);
	}

	std::vector<std::vector<int>> nodePerElement_tmp(pos_node.size());
	nodeSharedByElement = nodePerElement_tmp;


	// find all elements that share a node
	for (int eleInd = 0; eleInd < tetrahedrals.size(); eleInd++)
	{
		for (int j = 0; j < 4; j++)
		{
			int nodeInd = tetrahedrals[eleInd][j];
			nodeSharedByElement[nodeInd].push_back(eleInd);
		}
	}

	cal_DS_or_DM(false);


	update_F();


	calculateNodeMass();


	pos_node_prev = pos_node;
}

// calculate the DM_inv or DS matrix
void Mesh::cal_DS_or_DM(bool DS)
{
	if (DS)
	{
		tetra_DS.clear();
	}

	for (int eleInd = 0; eleInd < tetrahedrals.size(); eleInd++)
	{
		int si = tetrahedrals[eleInd][0], ai = tetrahedrals[eleInd][1], bi = tetrahedrals[eleInd][2], ci = tetrahedrals[eleInd][3];
		Eigen::Matrix3d DMS = Eigen::Matrix3d::Zero();
		DMS << pos_node[ai] - pos_node[si], pos_node[bi] - pos_node[si], pos_node[ci] - pos_node[si];


		if (DS)
		{
			tetra_DS.push_back(DMS);
		}
		else
		{
			tetra_DM_inv.push_back(DMS.inverse());
			tetra_vol.push_back(-DMS.determinant() / 6.0);
			//std::cout << "DMS.determinant() / 6.0 = " << DMS.determinant() / 6.0 << std::endl;
		}
	}
}

// output the mesh
void Mesh::output(int timestep)
{
	std::ofstream outfile2("./output/" + std::to_string(timestep) + ".obj", std::ios::trunc);
	for (int vert = 0; vert < pos_node.size(); ++vert)
	{
		outfile2 << std::scientific << std::setprecision(8) << "v " << pos_node[vert][0] << " " << pos_node[vert][1] << " " << pos_node[vert][2] << " " << std::endl;
	}
	outfile2.close();
}

// extract the surface mesh
void Mesh::extractSurfaceMesh()
{
	surfaceMesh.clear();

	std::map<TriangleFace, int> faceCount;
	for (const auto& tetr : tetrahedrals)
	{
		std::vector<TriangleFace> faces = {
			TriangleFace(tetr[0], tetr[1], tetr[2]),
			TriangleFace(tetr[0], tetr[1], tetr[3]),
			TriangleFace(tetr[0], tetr[2], tetr[3]),
			TriangleFace(tetr[1], tetr[2], tetr[3])
		};
		for (const auto& face : faces) {
			faceCount[face]++;
		}
	}

	for (const auto& pair : faceCount)
	{
		if (pair.second == 1)
		{
			std::vector<int> fc;
			fc.push_back(pair.first.n1);
			fc.push_back(pair.first.n2);
			fc.push_back(pair.first.n3);
			surfaceMesh.faces.push_back(fc);
		}
	}
	surfaceMesh.vertices = pos_node;

}

// export surface mesh
void Mesh::exportSurfaceMesh(std::string fileName, int timestep)
{
	surfaceMesh.vertices = pos_node;
	surfaceMesh.outputFile(fileName, timestep);
}

// update each tetrahedral's deformation gradient
void Mesh::update_F()
{
	tetra_F.resize(tetrahedrals.size());
	for (int i = 0; i < tetrahedrals.size(); i++) 
	{
		Eigen::Vector3d x0 = pos_node[tetrahedrals[i][0]];
		Eigen::Vector3d x1 = pos_node[tetrahedrals[i][1]];
		Eigen::Vector3d x2 = pos_node[tetrahedrals[i][2]];
		Eigen::Vector3d x3 = pos_node[tetrahedrals[i][3]];
		Eigen::Matrix<double, 3, 3> Ds;
		Ds << x1 - x0, x2 - x0, x3 - x0;
		tetra_F[i] = Ds * tetra_DM_inv[i];
	}
}

// calculate the mass of each node
void Mesh::calculateNodeMass()
{
	mass_node.clear(); // mass of a node
	for (int nd = 0; nd < nodeSharedByElement.size(); nd++)
	{
		double mass_ = 0;
		for (int te = 0; te < nodeSharedByElement[nd].size(); te++)
		{
			int treInd = nodeSharedByElement[nd][te];
			mass_ += tetra_vol[treInd] * materialMesh.density / 4.0;
		}
		mass_node.push_back(mass_);
	}
}

// calculate the bounding box of the mesh
std::pair<Eigen::Vector3d, Eigen::Vector3d> Mesh::calculateBoundingBox()
{
	double maxx = -1.0E9, maxy = -1.0E9, maxz = -1.0E9;
	double minx = 1.0E9, miny = 1.0E9, minz = 1.0E9;
	for (auto ver : pos_node) 
	{
		maxx = std::max(maxx, ver.x());
		maxy = std::max(maxy, ver.y());
		maxz = std::max(maxz, ver.z());

		minx = std::min(minx, ver.x());
		miny = std::min(miny, ver.y());
		minz = std::min(minz, ver.z());
	}
	Eigen::Vector3d min = { minx, miny, minz};
	Eigen::Vector3d max = { maxx, maxy, maxz};
	return std::make_pair(min , max);
}





