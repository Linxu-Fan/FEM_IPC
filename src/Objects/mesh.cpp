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