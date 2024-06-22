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
void Mesh::readMesh(std::string filePath, Material mat)
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
					materialInd.push_back(0);
				}
			}
		}
	}
	in.close();

	materialMesh.push_back(mat);

}


void Mesh::readMeshes(std::vector<std::pair<std::string, Material>> filePath_and_mat)
{
	int prevNodesNum = 0;
	int currNodesNum = 0;
	for (int ii = 0; ii < filePath_and_mat.size(); ii++)
	{
		std::string filePath = filePath_and_mat[ii].first;
		Material mat = filePath_and_mat[ii].second;

		materialMesh.push_back(mat);
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
						Eigen::Vector3d nd_pos = { std::stod(vecCoor[1]) , std::stod(vecCoor[2]) , std::stod(vecCoor[3]) };
						pos_node.push_back(nd_pos);
						currNodesNum += 1;
					}

					if (lineIndex >= elementStart && lineIndex <= elementStart + numElements - 1)
					{
						int numItemsLine = vecCoor.size(); // the number of items in a line
						if (vecCoor[1] == "4")
						{
							Eigen::Vector4i ele = { prevNodesNum + std::stoi(vecCoor[numItemsLine - 4]) - 1 ,  prevNodesNum + std::stoi(vecCoor[numItemsLine - 3]) - 1 ,  prevNodesNum + std::stoi(vecCoor[numItemsLine - 2]) - 1 ,  prevNodesNum + std::stoi(vecCoor[numItemsLine - 1]) - 1 };
							tetrahedrals.push_back(ele);
							materialInd.push_back(ii);
						}
					}
				}
			}
			in.close();

			prevNodesNum = currNodesNum;
		}

	}
}


// initialize the mesh after reading
void Mesh::initializeMesh() // initialize the mesh 
{
	for (int i = 0; i < pos_node.size(); i++)
	{
		Eigen::Vector3d vel = { 0,0,0 };
		vel_node.push_back(vel);

		boundaryCondition BC;
		boundaryCondition_node.push_back(BC);
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


	findBoundaryElements();


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
			tetra_vol.push_back(std::abs(DMS.determinant() / 6.0));
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

// find boundary elements including vertices, edges and triangles
void Mesh::findBoundaryElements()
{
	surfaceMesh.clear();
	surfaceMesh.vertices = pos_node;

	std::map<TriangleFace, std::pair<int, std::vector<int>>> faceCount;
	for (const auto& tetr : tetrahedrals)
	{
		std::vector<TriangleFace> faces = {
			TriangleFace(tetr[0], tetr[1], tetr[2]),
			TriangleFace(tetr[0], tetr[3], tetr[1]),
			TriangleFace(tetr[3], tetr[2], tetr[1]),
			TriangleFace(tetr[0], tetr[2], tetr[3])
		};
		for (int i = 0; i < 4; i++)
		{
			faceCount[faces[i]].first++;
			if (i == 0)
			{
				faceCount[faces[i]].second.push_back(tetr[0]);
				faceCount[faces[i]].second.push_back(tetr[1]);
				faceCount[faces[i]].second.push_back(tetr[2]);
			}
			else if (i == 1)
			{
				faceCount[faces[i]].second.push_back(tetr[0]);
				faceCount[faces[i]].second.push_back(tetr[3]);
				faceCount[faces[i]].second.push_back(tetr[1]);
			}
			else if (i == 2)
			{
				faceCount[faces[i]].second.push_back(tetr[3]);
				faceCount[faces[i]].second.push_back(tetr[2]);
				faceCount[faces[i]].second.push_back(tetr[1]);
			}
			else if (i == 3)
			{
				faceCount[faces[i]].second.push_back(tetr[0]);
				faceCount[faces[i]].second.push_back(tetr[2]);
				faceCount[faces[i]].second.push_back(tetr[3]);
			}

		}
	}
	for (const auto& pair : faceCount)
	{
		if (pair.second.first == 1)
		{
			surfaceMesh.faces.push_back(pair.second.second);

			Eigen::Vector3i tri = { pair.second.second[0], pair.second.second[1], pair.second.second[2]};
			boundaryTriangles.push_back(tri);
		}
	}

	std::set<std::string> allEdges;
	// find boundary vertices
	for (int h = 0; h < boundaryTriangles.size(); h++)
	{
		Eigen::Vector3i tri = boundaryTriangles[h];
		boundaryVertices[tri[0]].insert(h);
		boundaryVertices[tri[1]].insert(h);
		boundaryVertices[tri[2]].insert(h);


		// find all edges
		std::string edge1 = std::to_string(std::min(tri[0], tri[1])) + "#" + std::to_string(std::max(tri[0], tri[1]));
		std::string edge2 = std::to_string(std::min(tri[1], tri[2])) + "#" + std::to_string(std::max(tri[1], tri[2]));
		std::string edge3 = std::to_string(std::min(tri[2], tri[0])) + "#" + std::to_string(std::max(tri[2], tri[0]));
		allEdges.insert(edge1);
		allEdges.insert(edge2);
		allEdges.insert(edge3);
	}


	// find boundary edges
	for (std::set<std::string>::iterator it = allEdges.begin(); it != allEdges.end(); it++)
	{
		std::string edge = *it;
		std::vector<std::string> seg = split(edge, "#");
		int v1 = std::stoi(seg[0]);
		int v2 = std::stoi(seg[1]);
		Eigen::Vector2i edg = {v1, v2};

		std::set<int> v1Tris = boundaryVertices[v1], v2Tris = boundaryVertices[v2];
		std::vector<int> edgeTris;
		std::set_intersection(v1Tris.begin(), v1Tris.end(), v2Tris.begin(), v2Tris.end(),std::back_inserter(edgeTris));
		Eigen::Vector2i tris = { edgeTris[0], edgeTris[1]};

		boundaryEdges.push_back(std::make_pair(edg, tris));
	}
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
			int matInd = materialInd[treInd];
			mass_ += tetra_vol[treInd] * materialMesh[matInd].density / 4.0;
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





