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
		if (1)
		{
			for (int m = 0; m < faces[k].size(); m++)
			{
				outfile9 << faces[k][m] + 1 << " ";
			}
		}
		else
		{
			for (int m = faces[k].size() - 1; m >= 0; m--)
			{
				outfile9 << faces[k][m] + 1 << " ";
			}
		}
		outfile9 << std::endl;
	}
	outfile9.close();
}


////////////////////////////////////////////////////////////////////////
// Is it possible that node and element are not placed in order? If possible, then the reading code may crash.
////////////////////////////////////////////////////////////////////////
// read .msh mesh file
void Mesh::readMesh(meshConfiguration& config)
{
	// !!!!!!!!!!!!!!!!! This function if not fully implemented. DONOT USE IT.

	//std::ifstream in;
	//in.open(config.filePath);
	//std::string line;
	//int nodeStart = 100000000000; // the line index where a node starts
	//int numNodes = 100000000000; // number of nodes
	//int elementStart = 100000000000; // the element index where an element starts
	//int numElements = 100000000000; // number of elements
	//int lineIndex = -1; // current line index

	//while (getline(in, line))
	//{
	//	if (line.size() > 0)
	//	{
	//		lineIndex += 1;
	//		
	//		std::vector<std::string> vecCoor = split(line, " ");
	//		if (vecCoor[0] == "$Nodes")
	//		{
	//			nodeStart = lineIndex + 2;
	//		}
	//		if (lineIndex == nodeStart - 1)
	//		{
	//			numNodes = std::stoi(vecCoor[0]);
	//		}

	//		if (vecCoor[0] == "$Elements")
	//		{
	//			elementStart = lineIndex + 2;
	//		}
	//		if (lineIndex == elementStart - 1)
	//		{
	//			numElements = std::stoi(vecCoor[0]);				
	//		}


	//		if (lineIndex >= nodeStart && lineIndex <= nodeStart + numNodes - 1)
	//		{
	//			Eigen::Vector3d nd_pos = { std::stod(vecCoor[1]) , std::stod(vecCoor[2]) , std::stod(vecCoor[3])};
	//			pos_node.push_back(nd_pos);
	//			vel_node.push_back(config.velocity);
	//		}

	//		if (lineIndex >= elementStart && lineIndex <= elementStart + numElements - 1)
	//		{
	//			int numItemsLine = vecCoor.size(); // the number of items in a line
	//			if (vecCoor[1] == "4")
	//			{
	//				Eigen::Vector4i ele = { std::stoi(vecCoor[numItemsLine - 4]) - 1 ,std::stoi(vecCoor[numItemsLine - 3]) - 1 ,std::stoi(vecCoor[numItemsLine - 2]) - 1,std::stoi(vecCoor[numItemsLine - 1]) - 1 };
	//				tetrahedrals.push_back(ele);
	//				materialInd.push_back(0);
	//			}
	//		}
	//	}
	//}
	//in.close();

	//materialMesh.push_back(config.mesh_material);

}


void Mesh::readMeshes(std::vector<meshConfiguration>& config)
{
	int prevNodesNum = 0;
	int currNodesNum = 0;
	for (int ii = 0; ii < config.size(); ii++)
	{
		std::string filePath = config[ii].filePath;
		Material mat = config[ii].mesh_material;
		Eigen::Vector3d scale = config[ii].scale;
		Eigen::Vector3d translation = config[ii].translation;
		// Create a transformation matrix
		Eigen::Affine3d rotation = Eigen::Affine3d::Identity();
		rotation.translate(-config[ii].rotation_point)
			.rotate(Eigen::AngleAxisd(config[ii].rotation_angle[0], Eigen::Vector3d::UnitX()))
			.rotate(Eigen::AngleAxisd(config[ii].rotation_angle[1], Eigen::Vector3d::UnitY()))
			.rotate(Eigen::AngleAxisd(config[ii].rotation_angle[2], Eigen::Vector3d::UnitZ()))
			.translate(config[ii].rotation_point);



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
						Eigen::Vector3d nd_pos = { std::stod(vecCoor[1]) * scale[0] , std::stod(vecCoor[2]) * scale[1] , std::stod(vecCoor[3]) * scale[2] };
						note_node.push_back(config[ii].note);
						pos_node.push_back(rotation * nd_pos + translation);
						vel_node.push_back(config[ii].velocity);
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
	update_F(1);
	calculateNodeMass();
	findBoundaryElements();
	updateBoundaryElementsInfo();

	pos_node_prev = pos_node;
	pos_node_Rest = pos_node;
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
	// give each edge an unique index
	int edgeIndex = 0;
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


		int emin = std::min(v1, v2), emax = std::max(v1, v2);
		boundaryEdges[emin][emax] = tris;
		boundaryVertices_egde[v1].insert(edgeIndex);
		boundaryVertices_egde[v2].insert(edgeIndex);
		boundaryEdge_index[emin][emax] = edgeIndex;
		edgeIndex += 1;
	}
	// reverse index
	for (std::map<int, std::map<int, int>>::iterator it1 = boundaryEdge_index.begin(); it1 != boundaryEdge_index.end(); it1++)
	{
		for (std::map<int, int>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++)
		{
			Eigen::Vector2i ed = {it1->first, it2->first};
			int index = it2->second;
			index_boundaryEdge[index] = ed;
			index_boundaryEdge_vec.push_back(index);
		}
	}


	for (std::map<int, std::set<int>>::iterator it = boundaryVertices.begin(); it != boundaryVertices.end(); it++)
	{
		boundaryVertices_vec.push_back(it->first);
	}

}

// update boundary elements' information: area
void Mesh::updateBoundaryElementsInfo()
{
	boundaryVertices_area.clear();
	boundaryEdges_area.clear();
	boundaryTriangles_area.clear();

	// calculate the area of each triangle
	for (int i = 0; i < boundaryTriangles.size(); i++)
	{
		Eigen::Vector3i triangle = boundaryTriangles[i];
		Eigen::Vector3d t1Coor = pos_node[triangle[0]], t2Coor = pos_node[triangle[1]], t3Coor = pos_node[triangle[2]];
		Eigen::Vector3d t1t2 = t2Coor - t1Coor, t1t3 = t3Coor - t1Coor;
		double triArea = 1.0 / 2.0 * (t1t2.cross(t1t3)).norm();
		boundaryTriangles_area.push_back(triArea);
	}

	// calculate the distributed area of each vertex
	for (std::map<int, std::set<int>>::iterator it = boundaryVertices.begin(); it != boundaryVertices.end(); it++)
	{
		std::set<int> incidentTriangles = it->second;
		double vertArea = 0;
		for (std::set<int>::iterator it = incidentTriangles.begin(); it != incidentTriangles.end(); it++)
		{
			vertArea += boundaryTriangles_area[*it];
		}
		boundaryVertices_area[it->first] = vertArea;
	}

	// calculate the distributed area of each edge
	for (std::map<int, std::map<int, Eigen::Vector2i>>::iterator it1 = boundaryEdges.begin(); it1 != boundaryEdges.end(); it1++)
	{
		int v1 = it1->first;
		std::map<int, Eigen::Vector2i> trisInd = it1->second;	
		for (std::map<int, Eigen::Vector2i>::iterator it2 = trisInd.begin(); it2 != trisInd.end(); it2++)
		{
			int v2 = it2->first;	
			double edgeArea = 0;
			edgeArea += boundaryTriangles_area[it2->second[0]];
			edgeArea += boundaryTriangles_area[it2->second[1]];
			boundaryEdges_area[v1][v2] = edgeArea;
		}
	}

}

// export surface mesh
void Mesh::exportSurfaceMesh(std::string fileName, int timestep)
{
	surfaceMesh.vertices = pos_node;
	surfaceMesh.outputFile(fileName, timestep);
}

// update each tetrahedral's deformation gradient
void Mesh::update_F(int numOfThreads)
{
	tetra_F.resize(tetrahedrals.size());
#pragma omp parallel for num_threads(numOfThreads)
	for (int i = 0; i < tetrahedrals.size(); i++) 
	{
		Eigen::Vector3d x0 = pos_node[tetrahedrals[i][0]];
		Eigen::Vector3d x1 = pos_node[tetrahedrals[i][1]];
		Eigen::Vector3d x2 = pos_node[tetrahedrals[i][2]];
		Eigen::Vector3d x3 = pos_node[tetrahedrals[i][3]];
		Eigen::Matrix<double, 3, 3> Ds;
		Ds << x1 - x0, x2 - x0, x3 - x0;
		tetra_F[i] = Ds * tetra_DM_inv[i];

		assert(tetra_F[i].determinant() > 0);
		if (tetra_F[i].determinant() < 0)
		{
			std::cout << "Tetrahedral reverse!" << std::endl;
		}

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

double Mesh::calLargestEdgeLength()
{
	double largestLength = -1.0E9;
	for (std::map<int, Eigen::Vector2i>::iterator it = index_boundaryEdge.begin(); it != index_boundaryEdge.end(); it++)
	{
		int v1_index = it->second[0], v2_index = it->second[1];
		Eigen::Vector3d v1 = pos_node[v1_index], v2 = pos_node[v2_index];
		double length = (v1 - v2).norm();
		if (largestLength < length)
		{
			largestLength = length;
		}
	}
	return largestLength;
}


double Mesh::calBBXDiagSize()
{
	std::pair<Eigen::Vector3d, Eigen::Vector3d> BBX = calculateBoundingBox();
	return (BBX.first - BBX.second).norm();
}


