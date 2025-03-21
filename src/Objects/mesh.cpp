#include "mesh.h"


void surface_Info::updateBEInfo(std::vector<Eigen::Vector3d>& vertices, std::vector<Eigen::Vector3i>& faces)
{
	boundaryTriangles = faces;


	std::set<std::string> allEdges;
	// find boundary vertices
	for (int h = 0; h < boundaryTriangles.size(); h++)
	{
		Eigen::Vector3i tri = boundaryTriangles[h];
		boundaryVertices[tri[0]].insert(h);
		boundaryVertices[tri[1]].insert(h);
		boundaryVertices[tri[2]].insert(h);


		// find all edges
		std::string edge1 = std::to_string(std::min(tri[0], tri[1])) + "#"
			+ std::to_string(std::max(tri[0], tri[1]));
		std::string edge2 = std::to_string(std::min(tri[1], tri[2])) + "#"
			+ std::to_string(std::max(tri[1], tri[2]));
		std::string edge3 = std::to_string(std::min(tri[2], tri[0])) + "#"
			+ std::to_string(std::max(tri[2], tri[0]));
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
		Eigen::Vector2i edg = { v1, v2 };

		std::set<int> v1Tris = boundaryVertices[v1], v2Tris = boundaryVertices[v2];
		std::vector<int> edgeTris;
		std::set_intersection(v1Tris.begin(), v1Tris.end(), v2Tris.begin(),
			v2Tris.end(), std::back_inserter(edgeTris));
		Eigen::Vector2i tris = { edgeTris[0], edgeTris[1] };


		int emin = std::min(v1, v2), emax = std::max(v1, v2);
		boundaryEdges[emin][emax] = tris;
		boundaryVertices_egde[v1].insert(edgeIndex);
		boundaryVertices_egde[v2].insert(edgeIndex);
		boundaryEdge_index[emin][emax] = edgeIndex;
		edgeIndex += 1;
	}
	// reverse index
	for (std::map<int, std::map<int, int>>::iterator it1 = boundaryEdge_index.begin();
		it1 != boundaryEdge_index.end(); it1++)
	{
		for (std::map<int, int>::iterator it2 = it1->second.begin();
			it2 != it1->second.end(); it2++)
		{
			Eigen::Vector2i ed = { it1->first, it2->first };
			int index = it2->second;
			index_boundaryEdge[index] = ed;
			index_boundaryEdge_vec.push_back(index);
		}
	}


	for (std::map<int, std::set<int>>::iterator it = boundaryVertices.begin();
		it != boundaryVertices.end(); it++)
	{
		boundaryVertices_vec.push_back(it->first);
	}


	// calculate the area of each triangle
	for (int i = 0; i < boundaryTriangles.size(); i++)
	{
		Eigen::Vector3i triangle = boundaryTriangles[i];
		Eigen::Vector3d t1Coor = vertices[triangle[0]], t2Coor = vertices[triangle[1]],
			t3Coor = vertices[triangle[2]];
		Eigen::Vector3d t1t2 = t2Coor - t1Coor, t1t3 = t3Coor - t1Coor;
		double triArea = 1.0 / 2.0 * (t1t2.cross(t1t3)).norm();
		boundaryTriangles_area.push_back(triArea);
	}

	// calculate the distributed area of each vertex
	for (std::map<int, std::set<int>>::iterator it = boundaryVertices.begin();
		it != boundaryVertices.end(); it++)
	{
		std::set<int> incidentTriangles = it->second;
		double vertArea = 0;
		for (std::set<int>::iterator it = incidentTriangles.begin();
			it != incidentTriangles.end(); it++)
		{
			vertArea += boundaryTriangles_area[*it];
		}
		boundaryVertices_area[it->first] = vertArea;
	}

	// calculate the distributed area of each edge
	for (std::map<int, std::map<int, Eigen::Vector2i>>::iterator it1 = boundaryEdges.begin();
		it1 != boundaryEdges.end(); it1++)
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


////////////////////////////////////////////////////////////////////////
// Is it possible that node and element are not placed in order? If possible, then the reading code may crash.
////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////
// Tetrahedral mesh for simulation
//////////////////////////////////////////////

// read .msh mesh file
void tetMesh::readMesh(meshConfiguration& config)
{
	materialTetMesh = config.mesh_material;
	tetMeshNote = config.note;
	Eigen::Vector3d scale = config.scale;
	Eigen::Vector3d translation = config.translation;
	// Create a transformation matrix
	Eigen::Affine3d rotation = Eigen::Affine3d::Identity();
	rotation.translate(-config.rotation_point)
		.rotate(Eigen::AngleAxisd(config.rotation_angle[0], Eigen::Vector3d::UnitX()))
		.rotate(Eigen::AngleAxisd(config.rotation_angle[1], Eigen::Vector3d::UnitY()))
		.rotate(Eigen::AngleAxisd(config.rotation_angle[2], Eigen::Vector3d::UnitZ()))
		.translate(config.rotation_point);

	
	{
		std::ifstream in;
		in.open(config.filePath);
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
					Eigen::Vector3d nd_pos = { std::stod(vecCoor[1]) * scale[0] ,
						std::stod(vecCoor[2]) * scale[1] , std::stod(vecCoor[3]) * scale[2] };
					pos_node.push_back(rotation * nd_pos + translation);
					pos_node_Rest.push_back(rotation * nd_pos + translation);
					vel_node.push_back(config.velocity);
				}

				if (lineIndex >= elementStart && lineIndex <= elementStart + numElements - 1)
				{
					int numItemsLine = vecCoor.size(); // the number of items in a line
					if (vecCoor[1] == "4")
					{
						Eigen::Vector4i ele = {std::stoi(vecCoor[numItemsLine - 4]) - 1 , 
							std::stoi(vecCoor[numItemsLine - 3]) - 1 ,  std::stoi(vecCoor[numItemsLine - 2]) - 1 ,  
							std::stoi(vecCoor[numItemsLine - 1]) - 1 };
						tetrahedrals.push_back(ele);
					}
				}
			}
		}
		in.close();

	}

	
}

// initialize the mesh after reading
void tetMesh::initializeTetMesh() // initialize the mesh 
{
	boundaryCondition_node.resize(pos_node.size());

	cal_DM_inv();
	update_F(1);
	calculateNodeMass();
	findSurfaceMesh();

}

// calculate the DM_inv or DS matrix
void tetMesh::cal_DM_inv()
{
	for (int eleInd = 0; eleInd < tetrahedrals.size(); eleInd++)
	{
		int si = tetrahedrals[eleInd][0], ai = tetrahedrals[eleInd][1],
			bi = tetrahedrals[eleInd][2], ci = tetrahedrals[eleInd][3];
		Eigen::Matrix3d DMS = Eigen::Matrix3d::Zero();
		DMS << pos_node[ai] - pos_node[si], pos_node[bi] - pos_node[si], pos_node[ci] - pos_node[si];

		tetra_DM_inv.push_back(DMS.inverse());
		tetra_vol.push_back(std::abs(DMS.determinant() / 6.0));		
	}

}

// update each tetrahedral's deformation gradient
void tetMesh::update_F(int numOfThreads)
{
	tetra_F.resize(tetrahedrals.size());
#pragma omp parallel for num_threads(numOfThreads)
	for (int i = 0; i < tetrahedrals.size(); i++)
	{
		tetra_F[i] = Eigen::Matrix3d::Identity();
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
			exit(0);
		}

	}

}

// calculate the mass of each node
void tetMesh::calculateNodeMass()
{
	std::vector<std::vector<int>> nodeSharedByElement(pos_node.size());
	// find all elements that share a node
	for (int eleInd = 0; eleInd < tetrahedrals.size(); eleInd++)
	{
		for (int j = 0; j < 4; j++)
		{
			int nodeInd = tetrahedrals[eleInd][j];
			nodeSharedByElement[nodeInd].push_back(eleInd);
		}
	}


	mass_node.clear(); // mass of a node
	for (int nd = 0; nd < nodeSharedByElement.size(); nd++)
	{
		double mass_ = 0;
		for (int te = 0; te < nodeSharedByElement[nd].size(); te++)
		{
			int treInd = nodeSharedByElement[nd][te];
			mass_ += tetra_vol[treInd] * materialTetMesh.density / 4.0;
		}
		mass_node.push_back(mass_);
	}
}

// find boundary elements including vertices, edges and triangles
void tetMesh::findSurfaceMesh()
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
			Eigen::Vector3i tri = { pair.second.second[0], pair.second.second[1], pair.second.second[2] };
			surfaceMesh.faces.push_back(tri);
		}
	}
}

// export surface mesh
void tetMesh::exportSurfaceMesh(std::string fileName, int timestep)
{
	surfaceMesh.vertices = pos_node;
	surfaceMesh.outputFile(fileName, timestep);
}

// output the mesh
void tetMesh::output(int timestep)
{
	std::ofstream outfile2("./output/" + std::to_string(timestep) + ".obj", std::ios::trunc);
	for (int vert = 0; vert < pos_node.size(); ++vert)
	{
		outfile2 << std::scientific << std::setprecision(8) << "v " 
			<< pos_node[vert][0] << " " << pos_node[vert][1] << " " 
			<< pos_node[vert][2] << " " << std::endl;
	}
	outfile2.close();
}

void tetMesh::exportEdges(std::string fileName)
{
	std::set<std::pair<int, int>> edges_set;
	for (int i = 0; i < tetrahedrals.size(); ++i) {
		int n[4] = { tetrahedrals[i][0], tetrahedrals[i][1], tetrahedrals[i][2], tetrahedrals[i][3] };

		// Generate all combinations of edges in a tetrahedron
		for (int a = 0; a < 4; ++a) {
			for (int b = a + 1; b < 4; ++b) {
				int idx1 = n[a];
				int idx2 = n[b];
				if (idx1 > idx2) std::swap(idx1, idx2);
				edges_set.insert(std::make_pair(idx1, idx2));
			}
		}
	}

	std::ofstream outfile2("./output/" + fileName + ".obj", std::ios::trunc);
	for (int vert = 0; vert < pos_node.size(); ++vert)
	{
		outfile2 << std::scientific << std::setprecision(8) << "v "
			<< pos_node[vert][0] << " " << pos_node[vert][1] << " "
			<< pos_node[vert][2] << " " << std::endl;
	}

	for (const auto& edge : edges_set) 
	{
		outfile2 << "l " << edge.first + 1 << " " << edge.second + 1 << std::endl;  // Convert back to one-based indexing
	}

	outfile2.close();

}




void Mesh::createGlobalSimulationMesh()
{
	// calculate the number of meshes in the scene
	num_meshes = objectsTetMesh.size();

	// assemble material vector
	std::map<std::string, int> materialNameIndex; 
	int matId = 0, meshId = 0;
	for (std::map<std::string, tetMesh>::iterator it = objectsTetMesh.begin(); 
		it != objectsTetMesh.end(); it++)
	{
		tetMesh obj_tet = it->second;
		if (materialNameIndex.find(obj_tet.materialTetMesh.name) == materialNameIndex.end())
		{
			materialNameIndex[obj_tet.materialTetMesh.name] = matId;
			materialMesh.push_back(obj_tet.materialTetMesh);
			matId += 1;
		}

		tetMeshIndex[obj_tet.tetMeshNote] = meshId;
		meshId += 1;
	}

	// assemble node vectors
	for (std::map<std::string, tetMesh>::iterator it = objectsTetMesh.begin(); 
		it != objectsTetMesh.end(); it++)
	{
		tetMesh obj_tet = it->second;
		int currentMeshID = tetMeshIndex[obj_tet.tetMeshNote];
		// per node 
		for (int n = 0; n < obj_tet.pos_node.size(); n++)
		{
			pos_node.push_back(obj_tet.pos_node[n]);
			pos_node_Rest.push_back(obj_tet.pos_node_Rest[n]);
			vel_node.push_back(obj_tet.vel_node[n]);
			mass_node.push_back(obj_tet.mass_node[n]);
			pos_node_prev.push_back(obj_tet.pos_node_Rest[n]);
			boundaryCondition_node.push_back(obj_tet.boundaryCondition_node[n]);
			note_node.push_back(obj_tet.tetMeshNote);

			Eigen::Vector2i index_vert = { currentMeshID , 0 };
			index_node.push_back(index_vert);
		}
	}
	elastForce_node.resize(pos_node.size(), Eigen::Vector3d::Zero());
	contactForce_node.resize(pos_node.size(), Eigen::Vector3d::Zero());

	// assemble tetrahedral vectors
	for (std::map<std::string, tetMesh>::iterator it = objectsTetMesh.begin(); 
		it != objectsTetMesh.end(); it++)
	{
		tetMesh obj_tet = it->second;
		// per node 
		for (int n = 0; n < obj_tet.tetrahedrals.size(); n++)
		{
			tetra_vol.push_back(obj_tet.tetra_vol[n]);
			tetra_DM_inv.push_back(obj_tet.tetra_DM_inv[n]);
			tetra_F.push_back(obj_tet.tetra_F[n]);
			materialInd.push_back(materialNameIndex[obj_tet.materialTetMesh.name]);
			Dp.push_back(0);
		}
	}

	// modify tetrahedral elements 
	int currNodesNum = 0;
	for (std::map<std::string, tetMesh>::iterator it = objectsTetMesh.begin(); 
		it != objectsTetMesh.end(); it++)
	{
		tetMesh obj_tet = it->second;
		// per node 
		for (int n = 0; n < obj_tet.tetrahedrals.size(); n++)
		{
			Eigen::Vector4i tetra = obj_tet.tetrahedrals[n];
			tetra += currNodesNum * Eigen::Vector4i::Ones();
			tetrahedrals.push_back(tetra);
		}
		currNodesNum += obj_tet.pos_node.size();
	}

	// find tetrahedrons that share a node
	tetrahedrals_node.resize(pos_node.size());
	for (int i = 0; i < tetrahedrals.size(); i++)
	{
		Eigen::Vector4i tet = tetrahedrals[i];
		for (int j = 0; j < 4; j++)
		{
			tetrahedrals_node[tet[j]].push_back(i);
		}

	}

	// find boundary
	findSurfaceMesh();
	surfaceInfo.updateBEInfo(surfaceMesh.vertices, surfaceMesh.faces);
	for (std::map<int, std::set<int>>::iterator it = surfaceInfo.boundaryVertices.begin();
		it != surfaceInfo.boundaryVertices.end(); it++)
	{
		index_node[it->first][1] = 1;
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
	for (std::map<int, Eigen::Vector2i>::iterator it = surfaceInfo.index_boundaryEdge.begin();
		it != surfaceInfo.index_boundaryEdge.end(); it++)
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

void Mesh::sample_MLS_points_inside_tetrahedral(objMeshFormat& crack, int& tetIndex, int& num_points, double& radius)
{
	Eigen::Vector4i tet = tetrahedrals[tetIndex];
	Eigen::Vector3d V1 = pos_node_Rest[tet[0]], V2 = pos_node_Rest[tet[1]], 
		V3 = pos_node_Rest[tet[2]], V4 = pos_node_Rest[tet[3]];



	for (int i = 0; i < num_points; i++)
	{
		Eigen::Vector3d pt = randomPointInTetrahedron(V1, V2, V3, V4);		
		std::vector<int> validNodes;

		int maxLayers = 2;	
		std::vector<int> neigNodes_tet = find_Neighbour_Nodes_tetrahedral(tetIndex, maxLayers);
		for (int h = 0; h < neigNodes_tet.size(); h++)
		{
			int ni = neigNodes_tet[h];
			Eigen::Vector3d vertPos = pos_node[ni];
			if ((vertPos - pt).norm() < radius )
			{
				if (!crack.checkIfMeshIntersectWithLine(vertPos, pt))
				{
					validNodes.push_back(ni);
				}
			}
		}
		
		if (validNodes.size() < 4)
		{
			std::cout << "Insufficient control points for MLS points! Num of ctr pts is "<< validNodes.size() << std::endl;
		}




		MLSPoints mpt;
		mpt.init_MLS(pt, tetra_vol[tetIndex] / (double)num_points, validNodes, "Gaussian", radius);

		MLSPoints_tet_map[tetIndex].push_back(mpt);

	}
	



}

std::vector<int> Mesh::find_Neighbour_Nodes_tetrahedral(int tetIndex, int maxLayers)
{
	std::set<int> visitedNodes;
	std::queue<std::pair<int, int>> nodeQueue; // pair of node index and current layer

	// Initialize the queue with the nodes of the given tetrahedral
	for (int i = 0; i < 4; ++i) 
	{
		int nodeIndex = tetrahedrals[tetIndex][i];
		nodeQueue.push({ nodeIndex, 0 });
		visitedNodes.insert(nodeIndex);
	}

	while (!nodeQueue.empty()) 
	{
		std::pair<int, int> node_and_layer = nodeQueue.front();
		int currentNode = node_and_layer.first, currentLayer = node_and_layer.second;
		nodeQueue.pop();

		if (currentLayer < maxLayers) 
		{
			// Find all tetrahedrals that share the current node
			for (int h = 0; h < tetrahedrals_node[currentNode].size(); h++) 
			{
				int tetIndex = tetrahedrals_node[currentNode][h];
				for (int i = 0; i < 4; ++i) 
				{
					int neighborNode = tetrahedrals[tetIndex][i];
					if (visitedNodes.find(neighborNode) == visitedNodes.end()) 
					{
						nodeQueue.push({ neighborNode, currentLayer + 1 });
						visitedNodes.insert(neighborNode);
					}
				}
			}
		}
	}

	return std::vector<int>(visitedNodes.begin(), visitedNodes.end());


}

void Mesh::sample_MLS_points(objMeshFormat& crack, int& num_points, double& radius, int& numOfThreads)
{
	std::vector<bool> tetsIntersected(tetrahedrals.size(), false);
#pragma omp parallel for num_threads(numOfThreads)
	for (int t = 0; t < tetrahedrals.size(); t++)
	{
		Eigen::Vector4i tet = tetrahedrals[t];
		Eigen::Vector3d tet_v0 = pos_node_Rest[tet[0]], tet_v1 = pos_node_Rest[tet[1]],
			tet_v2 = pos_node_Rest[tet[2]], tet_v3 = pos_node_Rest[tet[3]];

		if (crack.checkIfMeshIntersectWithTetrahedron(tet_v0, tet_v1, tet_v2, tet_v3))
		{
			tetsIntersected[t] = true;
		}
	}

	for (int h = 0; h < tetsIntersected.size(); h++)
	{
		if (tetsIntersected[h])
		{
			sample_MLS_points_inside_tetrahedral(crack, h, num_points, radius);
		}
	}


}



void Mesh_ABD::createGlobalSimulationMesh_ABD()
{
	createGlobalSimulationMesh();

	massMatrix_ABD.resize(num_meshes, Eigen::Matrix<double, 12, 12>::Zero());

	translation_prev_ABD.resize(num_meshes, Eigen::Vector3d::Zero());
	translation_ABD.resize(num_meshes, Eigen::Vector3d::Zero());
	translation_vel_ABD.resize(num_meshes, Eigen::Vector3d::Zero());

	deformation_prev_ABD.resize(num_meshes, Eigen::Matrix3d::Identity());
	deformation_ABD.resize(num_meshes, Eigen::Matrix3d::Identity());
	deformation_vel_ABD.resize(num_meshes, Eigen::Matrix3d::Zero());

	for (int i = 0; i < pos_node.size(); i++)
	{
		Eigen::Matrix<double, 3, 12> Jx = build_Jx_matrix_for_ABD(pos_node[i]);

		int AB_index = index_node[i][0];
		massMatrix_ABD[AB_index] += mass_node[i] * Jx.transpose() * Jx;
	}


	volume_ABD.resize(num_meshes, 0);
	for (int i = 0; i < tetrahedrals.size(); i++)
	{
		int vtInd = tetrahedrals[i][0];
		int ABInd = index_node[vtInd][0];
		volume_ABD[ABInd] += tetra_vol[i];
	}


}




//////////////////////////////////////////////
// Triangular mesh for simulation
//////////////////////////////////////////////

void triMesh::clear()
{
	triMeshIndex.clear(); 
	num_meshes = 0;


	note_node_surface.clear();
	contactForce_node_surface.clear();
	pos_node_surface.clear();
	pos_node_Rest_surface.clear();
	index_node_surface.clear();
	boundaryCondition_node_surface.clear();


	objMeshFormat surfaceMeshGlobal_;
	surface_Info surfaceInfo_;
	surfaceMeshGlobal = surfaceMeshGlobal_;
	surfaceInfo = surfaceInfo_;

	note_node_interior.clear();
	pos_node_interior.clear();
	pos_node_Rest_interior.clear();
	index_node_interior.clear();
	boundaryCondition_node_interior.clear();
	mass_node_interior.clear();
	vol_node_interior.clear();
}



void triMesh::createGlobalSimulationTriMesh_ABD(std::vector<meshConfiguration>& configs)
{
	readMeshes(configs);
	build_surface_mesh();
	sample_points_inside();
	update_ABD_info();

}


void triMesh::updateGlobalSimulationTriMesh_ABD()
{
	clear();


	num_meshes = allObjects.size();
	// read the 
	for (int i = 0; i < allObjects.size(); i++)
	{
		triMeshIndex[allObjects[i].objectNote] = i;
	}

	if (triMeshIndex.size() != allObjects.size()) // check if there exists identical object names
	{
		std::cout << "Identical object names exist! Please change the object's name!" << std::endl;
		exit(0);
	}



	build_surface_mesh();
	sample_points_inside();
	update_ABD_info();


	// !!!!!!!!!!!!!!!!!!!!!!!
// update points' position in the rest configuration
	std::vector<Eigen::Matrix3d> deformation_ABD_inverse(allObjects.size(), Eigen::Matrix3d::Identity());
#pragma omp parallel for
	for (int i = 0 ; i < allObjects.size(); i++)
	{
		deformation_ABD_inverse[i] = allObjects[i].deformation_ABD.inverse();
	}
#pragma omp parallel for
	for (int i = 0; i < pos_node_Rest_surface.size(); i++)
	{
		int obj_index = index_node_surface[i][0];
		if (allObjects[obj_index].need_update_rest_position == true)
		{
			Eigen::Vector3d curr_pos_minuse_trans = pos_node_Rest_surface[i] - allObjects[obj_index].translation_ABD;
			pos_node_Rest_surface[i] = deformation_ABD_inverse[obj_index] * curr_pos_minuse_trans;
		}
	}
	for (int i = 0; i < allObjects.size(); i++)
	{
		allObjects[i].need_update_rest_position = false;
	}

}


void triMesh::readMeshes(std::vector<meshConfiguration>& configs)
{
	num_meshes = configs.size();
	// read the 
	for (int i = 0; i < configs.size(); i++)
	{
		meshConfiguration config = configs[i];

		triMeshIndex[config.note] = i;
	
		Eigen::Vector3d scale = config.scale;
		Eigen::Vector3d translation = config.translation;
		// Create a transformation matrix
		Eigen::Affine3d rotation = Eigen::Affine3d::Identity();
		rotation.translate(-config.rotation_point)
			.rotate(Eigen::AngleAxisd(config.rotation_angle[0], Eigen::Vector3d::UnitX()))
			.rotate(Eigen::AngleAxisd(config.rotation_angle[1], Eigen::Vector3d::UnitY()))
			.rotate(Eigen::AngleAxisd(config.rotation_angle[2], Eigen::Vector3d::UnitZ()))
			.translate(config.rotation_point);

		objMeshFormat triMesh_;
		triMesh_.readObjFile(config.filePath, false, rotation, scale, translation);
		triMesh_.initialVelocity = config.velocity;


		ABD_Object obj_;
		obj_.objectMaterial = config.mesh_material;
		obj_.objectNote = config.note;
		obj_.breakable = config.breakable;
		obj_.objectSurfaceMesh = triMesh_;
		obj_.objectSurfaceMesh.updateVolume();
		obj_.translation_vel_ABD = config.velocity;
		obj_.per_point_volume = config.per_point_volume;
		allObjects.push_back(obj_);

	}

	if (triMeshIndex.size() != configs.size()) // check if there exists identical object names
	{
		std::cout << "Identical object names exist! Please change the object's name!" << std::endl;
		exit(0);
	}

}


void triMesh::build_surface_mesh()
{
	// surface points
	for (int i = 0; i < allObjects.size(); i++)
	{
		pos_node_surface.insert(pos_node_surface.end(), allObjects[i].objectSurfaceMesh.vertices.begin(), allObjects[i].objectSurfaceMesh.vertices.end());


		Eigen::Vector2i index_vert = { i , 1 };
		// store the note infor of each object
		for (int n = 0; n < allObjects[i].objectSurfaceMesh.vertices.size(); n++)
		{
			note_node_surface.push_back(allObjects[i].objectNote);
			index_node_surface.push_back(index_vert);
		}
	}
	pos_node_Rest_surface = pos_node_surface;

	// surface triangles
	int countNodeNum = 0;
	for (int i = 0; i < allObjects.size(); i++)
	{
		Eigen::Vector2i se_nodes = {0,0};
		se_nodes[0] = countNodeNum;
		for (int j = 0; j < allObjects[i].objectSurfaceMesh.faces.size(); j++)
		{
			Eigen::Vector3i face = allObjects[i].objectSurfaceMesh.faces[j];
			face += Eigen::Vector3i::Ones() * countNodeNum;
			surfaceMeshGlobal.faces.push_back(face);
		}
		countNodeNum += allObjects[i].objectSurfaceMesh.vertices.size();
		se_nodes[1] = countNodeNum;
		allObjects[i].objectSurfaceMeshes_node_start_end = se_nodes;
	}
	surfaceMeshGlobal.vertices = pos_node_surface;

	surfaceInfo.updateBEInfo(surfaceMeshGlobal.vertices, surfaceMeshGlobal.faces);

}


void triMesh::sample_points_inside()
{
	// sample points inside first
	for (int i = 0; i < allObjects.size(); i++)
	{

		int num_points = std::floor(allObjects[i].objectSurfaceMesh.volume / allObjects[i].per_point_volume);
		std::vector<Eigen::Vector3d> pts = allObjects[i].objectSurfaceMesh.sample_points_inside_mesh(num_points);
		pos_node_interior.insert(pos_node_interior.end(), pts.begin(), pts.end());


		Eigen::Vector2i index_vert = { i , 0 };
		// store the note infor of each object
		for (int n = 0; n < pts.size(); n++)
		{
			note_node_interior.push_back(allObjects[i].objectNote);
			index_node_interior.push_back(index_vert);
			vol_node_interior.push_back(allObjects[i].per_point_volume);
			mass_node_interior.push_back(allObjects[i].per_point_volume * allObjects[i].objectMaterial.density);
		}
	}
	pos_node_Rest_interior = pos_node_interior;



	boundaryCondition bc;
	// add boundary conditions
	for (int i = 0; i < pos_node_interior.size(); i++)
	{
		boundaryCondition_node_interior.push_back(bc);
	}
	for (int i = 0; i < pos_node_surface.size(); i++)
	{
		boundaryCondition_node_surface.push_back(bc);
	}
	contactForce_node_surface.resize(pos_node_surface.size(), Eigen::Vector3d::Zero());


}


void triMesh::update_ABD_info()
{

	massMatrix_ABD.resize(allObjects.size(), Eigen::Matrix<double, 12, 12>::Zero());
	for (int i = 0; i < pos_node_Rest_interior.size(); i++)
	{
		Eigen::Matrix<double, 3, 12> Jx = build_Jx_matrix_for_ABD(pos_node_Rest_interior[i]);

		int AB_index = index_node_interior[i][0];
		massMatrix_ABD[AB_index] += mass_node_interior[i] * Jx.transpose() * Jx;
	}
	for (int i = 0; i < allObjects.size(); i++)
	{
		translation_prev_ABD.push_back(allObjects[i].translation_prev_ABD);
		translation_vel_ABD.push_back(allObjects[i].translation_vel_ABD);
		translation_ABD.push_back(allObjects[i].translation_ABD);
		deformation_prev_ABD.push_back(allObjects[i].deformation_prev_ABD);
		deformation_vel_ABD.push_back(allObjects[i].deformation_vel_ABD);
		deformation_ABD.push_back(allObjects[i].deformation_ABD);

		volume_ABD.push_back(allObjects[i].objectSurfaceMesh.volume);
	}
	
}

void triMesh::exportSurfaceMesh(std::string fileName, int timestep)
{

	surfaceMeshGlobal.vertices = pos_node_surface;
	surfaceMeshGlobal.outputFile(fileName, timestep);
}

double triMesh::calLargestEdgeLength()
{
	double largestLength = -1.0E9;
	for (std::map<int, Eigen::Vector2i>::iterator it = surfaceInfo.index_boundaryEdge.begin();
		it != surfaceInfo.index_boundaryEdge.end(); it++)
	{
		int v1_index = it->second[0], v2_index = it->second[1];
		Eigen::Vector3d v1 = pos_node_surface[v1_index], v2 = pos_node_surface[v2_index];
		double length = (v1 - v2).norm();
		if (largestLength < length)
		{
			largestLength = length;
		}
	}
	return largestLength;
}

void triMesh::updateEachObjectSurfaceMesh()
{
	for (int i = 0; i < allObjects.size(); i++)
	{
		for (int j = allObjects[i].objectSurfaceMeshes_node_start_end[0]; j < allObjects[i].objectSurfaceMeshes_node_start_end[1]; j++)
		{
			int localIndex = j - allObjects[i].objectSurfaceMeshes_node_start_end[0];
			allObjects[i].objectSurfaceMesh.vertices[localIndex] = pos_node_surface[j];
		}
	}
}