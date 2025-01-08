#include "mesh.h"

void bounding_box::cal_min_max_ABD(const Vector12d& affine, const double& dilation)
{
	Eigen::Vector3d left_front_bottom_ = build_Jx_matrix_for_ABD(left_front_bottom) * affine;
	Eigen::Vector3d right_front_bottom_ = build_Jx_matrix_for_ABD(right_front_bottom) * affine;
	Eigen::Vector3d right_back_bottom_ = build_Jx_matrix_for_ABD(right_back_bottom) * affine;
	Eigen::Vector3d left_back_bottom_ = build_Jx_matrix_for_ABD(left_back_bottom) * affine;
	Eigen::Vector3d left_front_top_ = build_Jx_matrix_for_ABD(left_front_top) * affine;
	Eigen::Vector3d right_front_top_ = build_Jx_matrix_for_ABD(right_front_top) * affine;
	Eigen::Vector3d right_back_top_ = build_Jx_matrix_for_ABD(right_back_top) * affine;
	Eigen::Vector3d left_back_top_ = build_Jx_matrix_for_ABD(left_back_top) * affine;

	Eigen::Vector3d min1 = left_front_bottom_.cwiseMin(right_front_bottom_).cwiseMin(right_back_bottom_).cwiseMin(left_back_bottom_);
	Eigen::Vector3d max1 = left_front_bottom_.cwiseMax(right_front_bottom_).cwiseMax(right_back_bottom_).cwiseMax(left_back_bottom_);

	Eigen::Vector3d min2 = left_front_top_.cwiseMin(right_front_top_).cwiseMin(right_back_top_).cwiseMin(left_back_top_);
	Eigen::Vector3d max2 = left_front_top_.cwiseMax(right_front_top_).cwiseMax(right_back_top_).cwiseMax(left_back_top_);

	min = min1.cwiseMin(max1).cwiseMin(min2).cwiseMin(max2);
	max = max1.cwiseMax(max1).cwiseMax(min2).cwiseMax(max2);

	min -= Eigen::Vector3d(dilation, dilation, dilation);
	max += Eigen::Vector3d(dilation, dilation, dilation);
}

void bounding_box::merges(const bounding_box& other)
{
	min = min.cwiseMin(other.min).cwiseMin(other.max);
	max = max.cwiseMax(other.min).cwiseMax(other.max);
}

bool bounding_box::intersects(const bounding_box& other)
{
	return (min.x() <= other.max.x() && max.x() >= other.min.x() &&
		min.y() <= other.max.y() && max.y() >= other.min.y() &&
		min.z() <= other.max.z() && max.z() >= other.min.z());
}

void bounding_box::export_BBX_rest(std::string fileName)
{
	std::ofstream outfile9("./output/" + fileName + ".obj", std::ios::trunc);

	outfile9 << std::scientific << std::setprecision(8) << "v " << left_front_bottom[0] << " " << left_front_bottom[1] << " " << left_front_bottom[2] << std::endl;
	outfile9 << std::scientific << std::setprecision(8) << "v " << right_front_bottom[0] << " " << right_front_bottom[1] << " " << right_front_bottom[2] << std::endl;
	outfile9 << std::scientific << std::setprecision(8) << "v " << right_back_bottom[0] << " " << right_back_bottom[1] << " " << right_back_bottom[2] << std::endl;
	outfile9 << std::scientific << std::setprecision(8) << "v " << left_back_bottom[0] << " " << left_back_bottom[1] << " " << left_back_bottom[2] << std::endl;
	outfile9 << std::scientific << std::setprecision(8) << "v " << left_front_top[0] << " " << left_front_top[1] << " " << left_front_top[2] << std::endl;
	outfile9 << std::scientific << std::setprecision(8) << "v " << right_front_top[0] << " " << right_front_top[1] << " " << right_front_top[2] << std::endl;
	outfile9 << std::scientific << std::setprecision(8) << "v " << right_back_top[0] << " " << right_back_top[1] << " " << right_back_top[2] << std::endl;
	outfile9 << std::scientific << std::setprecision(8) << "v " << left_back_top[0] << " " << left_back_top[1] << " " << left_back_top[2] << std::endl;

	outfile9 << "l 1 2" << std::endl;
	outfile9 << "l 2 3" << std::endl;
	outfile9 << "l 3 4" << std::endl;
	outfile9 << "l 4 1" << std::endl;
	outfile9 << "l 5 6" << std::endl;
	outfile9 << "l 6 7" << std::endl;
	outfile9 << "l 7 8" << std::endl;
	outfile9 << "l 8 5" << std::endl;
	outfile9 << "l 1 5" << std::endl;
	outfile9 << "l 2 6" << std::endl;
	outfile9 << "l 3 7" << std::endl;
	outfile9 << "l 4 8" << std::endl;

	outfile9.close();
}

void bounding_box::export_BBX_world(std::string fileName)
{
	std::ofstream outfile9("./output/" + fileName + ".obj", std::ios::trunc);


	Eigen::Vector3d dxyz = max - min;
	double dx = dxyz[0], dy = dxyz[1], dz = dxyz[2];

	Eigen::Vector3d v1 = {min[0], min[1], min[2]};
	Eigen::Vector3d v2 = {min[0] + dx, min[1], min[2]};
	Eigen::Vector3d v3 = {min[0] + dx, min[1] + dy, min[2]};
	Eigen::Vector3d v4 = {min[0], min[1] + dy, min[2]};
	Eigen::Vector3d v5 = {min[0], min[1], min[2] + dz};
	Eigen::Vector3d v6 = {min[0] + dx, min[1], min[2] + dz };
	Eigen::Vector3d v7 = {min[0] + dx, min[1] + dy, min[2] + dz };
	Eigen::Vector3d v8 = {min[0], min[1] + dy, min[2] + dz };



	outfile9 << std::scientific << std::setprecision(8) << "v " << v1[0] << " " << v1[1] << " " << v1[2] << std::endl;
	outfile9 << std::scientific << std::setprecision(8) << "v " << v2[0] << " " << v2[1] << " " << v2[2] << std::endl;
	outfile9 << std::scientific << std::setprecision(8) << "v " << v3[0] << " " << v3[1] << " " << v3[2] << std::endl;
	outfile9 << std::scientific << std::setprecision(8) << "v " << v4[0] << " " << v4[1] << " " << v4[2] << std::endl;
	outfile9 << std::scientific << std::setprecision(8) << "v " << v5[0] << " " << v5[1] << " " << v5[2] << std::endl;
	outfile9 << std::scientific << std::setprecision(8) << "v " << v6[0] << " " << v6[1] << " " << v6[2] << std::endl;
	outfile9 << std::scientific << std::setprecision(8) << "v " << v7[0] << " " << v7[1] << " " << v7[2] << std::endl;
	outfile9 << std::scientific << std::setprecision(8) << "v " << v8[0] << " " << v8[1] << " " << v8[2] << std::endl;

	outfile9 << "l 1 2" << std::endl;
	outfile9 << "l 2 3" << std::endl;
	outfile9 << "l 3 4" << std::endl;
	outfile9 << "l 4 1" << std::endl;
	outfile9 << "l 5 6" << std::endl;
	outfile9 << "l 6 7" << std::endl;
	outfile9 << "l 7 8" << std::endl;
	outfile9 << "l 8 5" << std::endl;
	outfile9 << "l 1 5" << std::endl;
	outfile9 << "l 2 6" << std::endl;
	outfile9 << "l 3 7" << std::endl;
	outfile9 << "l 4 8" << std::endl;

	outfile9.close();
}

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

void triMesh::build_BVH_object(double dilation, int object_index)
{

	// initialize face BVH
	{
		std::vector<AABB> obj_AABBs_face(allObjects[object_index].objectSurfaceMesh.faces.size());
		for (int jj = 0; jj < allObjects[object_index].objectSurfaceMesh.faces.size(); ++jj)
		{
			AABB aabb_face;
			aabb_face.init(allObjects[object_index].pos_node_surface, allObjects[object_index].objectSurfaceMesh.faces[jj], jj, dilation);
			obj_AABBs_face[jj] = aabb_face;
			//std::cout << "aabb_face = " << aabb_face.vert_edge_face << std::endl;
		}
		allObjects[object_index].object_BVH_faces = buildBVH(obj_AABBs_face, 0, obj_AABBs_face.size(), 0);
	}

	// initialize edge BVH
	{
		std::vector<AABB> obj_AABBs_edge(allObjects[object_index].objectSurfaceMesh.edges.size());
		for (int jj = 0; jj < allObjects[object_index].objectSurfaceMesh.edges.size(); ++jj)
		{
			AABB aabb_edge;
			aabb_edge.init(allObjects[object_index].pos_node_surface, allObjects[object_index].objectSurfaceMesh.edges[jj], jj, dilation);
			obj_AABBs_edge[jj] = aabb_edge;
		}
		allObjects[object_index].object_BVH_edges = buildBVH(obj_AABBs_edge, 0, obj_AABBs_edge.size(), 0);
	}

	// initialize node BVH
	{
		std::vector<AABB> obj_AABBs_node(allObjects[object_index].objectSurfaceMesh.vertices.size());
		for (int jj = 0; jj < allObjects[object_index].objectSurfaceMesh.vertices.size(); ++jj)
		{
			AABB aabb_node;
			aabb_node.init(allObjects[object_index].pos_node_surface, jj, dilation);
			obj_AABBs_node[jj] = aabb_node;
		}
		allObjects[object_index].object_BVH_nodes = buildBVH(obj_AABBs_node, 0, obj_AABBs_node.size(), 0);
	}


}

void triMesh::build_BVH_object_advect(double dilation, int object_index)
{

	// initialize face BVH
	{
		std::vector<AABB> obj_AABBs_face(allObjects[object_index].objectSurfaceMesh.faces.size());
		for (int jj = 0; jj < allObjects[object_index].objectSurfaceMesh.faces.size(); ++jj)
		{
			AABB aabb_face;
			aabb_face.init_advect(allObjects[object_index].pos_node_surface, allObjects[object_index].objectSurfaceMesh.faces[jj], jj, dilation, allObjects[object_index].pos_node_surface_direction);
			obj_AABBs_face[jj] = aabb_face;
			//std::cout << "aabb_face = " << aabb_face.vert_edge_face << std::endl;
		}
		allObjects[object_index].object_BVH_faces = buildBVH(obj_AABBs_face, 0, obj_AABBs_face.size(), 0);
	}

	// initialize edge BVH
	{
		std::vector<AABB> obj_AABBs_edge(allObjects[object_index].objectSurfaceMesh.edges.size());
		for (int jj = 0; jj < allObjects[object_index].objectSurfaceMesh.edges.size(); ++jj)
		{
			AABB aabb_edge;
			aabb_edge.init_advect(allObjects[object_index].pos_node_surface, allObjects[object_index].objectSurfaceMesh.edges[jj], jj, dilation, allObjects[object_index].pos_node_surface_direction);
			obj_AABBs_edge[jj] = aabb_edge;
		}
		allObjects[object_index].object_BVH_edges = buildBVH(obj_AABBs_edge, 0, obj_AABBs_edge.size(), 0);
	}

	// initialize node BVH
	{
		std::vector<AABB> obj_AABBs_node(allObjects[object_index].objectSurfaceMesh.vertices.size());
		for (int jj = 0; jj < allObjects[object_index].objectSurfaceMesh.vertices.size(); ++jj)
		{
			AABB aabb_node;
			aabb_node.init_advect(allObjects[object_index].pos_node_surface, jj, dilation, allObjects[object_index].pos_node_surface_direction);
			obj_AABBs_node[jj] = aabb_node;
		}
		allObjects[object_index].object_BVH_nodes = buildBVH(obj_AABBs_node, 0, obj_AABBs_node.size(), 0);
	}


}

void triMesh::update_box_corner()
{
#pragma omp parallel for
	for (int ii = 0; ii < allObjects.size(); ii++)
	{
		build_BVH_object(0, ii);

		{
			Eigen::Vector3d dxyz = allObjects[ii].object_BVH_faces->box.max - allObjects[ii].object_BVH_faces->box.min;
			double dx = dxyz[0], dy = dxyz[1], dz = dxyz[2];

			//std::cout << "min = " << allObjects[ii].object_BVH_faces->box.min.transpose() << "; max = " << allObjects[ii].object_BVH_faces->box.max.transpose() << std::endl;

			allObjects[ii].BBX.left_front_bottom = allObjects[ii].object_BVH_faces->box.min;
			allObjects[ii].BBX.right_front_bottom = allObjects[ii].object_BVH_faces->box.min;
			allObjects[ii].BBX.right_back_bottom = allObjects[ii].object_BVH_faces->box.min;
			allObjects[ii].BBX.left_back_bottom = allObjects[ii].object_BVH_faces->box.min;
			allObjects[ii].BBX.left_front_top = allObjects[ii].object_BVH_faces->box.min;
			allObjects[ii].BBX.right_front_top = allObjects[ii].object_BVH_faces->box.min;
			allObjects[ii].BBX.right_back_top = allObjects[ii].object_BVH_faces->box.min;
			allObjects[ii].BBX.left_back_top = allObjects[ii].object_BVH_faces->box.min;

			allObjects[ii].BBX.right_front_bottom[0] += dx;

			allObjects[ii].BBX.right_back_bottom[0] += dx;
			allObjects[ii].BBX.right_back_bottom[1] += dy;

			allObjects[ii].BBX.left_back_bottom[1] += dy;


			allObjects[ii].BBX.left_front_top[2] += dz;

			allObjects[ii].BBX.right_front_top[0] += dx;
			allObjects[ii].BBX.right_front_top[2] += dz;

			allObjects[ii].BBX.right_back_top[0] += dx;
			allObjects[ii].BBX.right_back_top[1] += dy;
			allObjects[ii].BBX.right_back_top[2] += dz;

			allObjects[ii].BBX.left_back_top[1] += dy;
			allObjects[ii].BBX.left_back_top[2] += dz;
		}
	}
}

void triMesh::clear()
{
	triMeshIndex.clear(); 
	num_meshes = 0;
}

void triMesh::createGlobalSimulationTriMesh_ABD(std::vector<meshConfiguration>& configs)
{
	readMeshes(configs);
	build_surface_mesh();
	sample_points_inside();
	update_ABD_info();
	update_box_corner();
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
	
	//////////////////////////////
	// Need to update each object's surfac mesh
	//////////////////////////////

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
		triMesh_.updateMesh();
		triMesh_.initialVelocity = config.velocity;


		ABD_Object obj_;
		obj_.objectMaterial = config.mesh_material;
		obj_.objectNote = config.note;
		obj_.breakable = config.breakable;
		obj_.objectSurfaceMesh = triMesh_;
		obj_.objectSurfaceMesh.updateVolume();
		obj_.translation_vel_ABD = config.velocity;

		obj_.volume = obj_.objectSurfaceMesh.volume;

		obj_.surfaceInfo.updateBEInfo(obj_.objectSurfaceMesh.vertices, obj_.objectSurfaceMesh.faces);
		obj_.pos_node_surface = obj_.objectSurfaceMesh.vertices;
		obj_.pos_node_surface_prev = obj_.objectSurfaceMesh.vertices;
		obj_.pos_node_surface_direction.resize(obj_.pos_node_surface.size(), Eigen::Vector3d::Zero());
		obj_.contactForce_node_surface.resize(obj_.pos_node_surface.size(), Eigen::Vector3d::Zero());

		obj_.affine[3] = 1.0;
		obj_.affine[7] = 1.0;
		obj_.affine[11] = 1.0;

		obj_.affine_prev = obj_.affine;

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
		surfaceMeshGlobal.vertices.insert(surfaceMeshGlobal.vertices.end(), allObjects[i].objectSurfaceMesh.vertices.begin(), allObjects[i].objectSurfaceMesh.vertices.end());
	}
	

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
	}
}

void triMesh::sample_points_inside()
{
	// sample points inside first
	for (int i = 0; i < allObjects.size(); i++)
	{

		int num_points = std::floor(allObjects[i].objectSurfaceMesh.volume / allObjects[i].per_point_volume);
		std::vector<Eigen::Vector3d> pts = allObjects[i].objectSurfaceMesh.sample_points_inside_mesh(num_points);
		allObjects[i].pos_node_interior = pts;
		allObjects[i].pos_node_interior_prev = pts;
		allObjects[i].pos_node_Rest_interior = pts;


		// store the note infor of each object
		for (int n = 0; n < pts.size(); n++)
		{
			allObjects[i].vol_node_interior.push_back(allObjects[i].per_point_volume);
			allObjects[i].mass_node_interior.push_back(allObjects[i].per_point_volume * allObjects[i].objectMaterial.density);

		}
	}

}

void triMesh::update_ABD_info()
{
	for (int i = 0; i < allObjects.size(); i++)
	{
		for (int j = 0; j < allObjects[i].pos_node_interior.size(); j++)
		{
			Eigen::Matrix<double, 3, 12> Jx = build_Jx_matrix_for_ABD(allObjects[i].pos_node_interior[j]);

			allObjects[i].massMatrix_ABD += allObjects[i].mass_node_interior[j] * Jx.transpose() * Jx;
		}
	}

}

void triMesh::exportSurfaceMesh(std::string fileName, int timestep)
{
	surfaceMeshGlobal.vertices.clear();
	for (int i = 0; i < allObjects.size(); i++)
	{
		surfaceMeshGlobal.vertices.insert(surfaceMeshGlobal.vertices.end(), allObjects[i].pos_node_surface.begin(), allObjects[i].pos_node_surface.end());
	}
	surfaceMeshGlobal.outputFile(fileName, timestep);
}

double triMesh::calLargestEdgeLength()
{
	double largestLength = -1.0E9;
	//for (std::map<int, Eigen::Vector2i>::iterator it = surfaceInfo.index_boundaryEdge.begin();
	//	it != surfaceInfo.index_boundaryEdge.end(); it++)
	//{
	//	int v1_index = it->second[0], v2_index = it->second[1];
	//	Eigen::Vector3d v1 = pos_node_surface[v1_index], v2 = pos_node_surface[v2_index];
	//	double length = (v1 - v2).norm();
	//	if (largestLength < length)
	//	{
	//		largestLength = length;
	//	}
	//}
	return largestLength;
}

void triMesh::updateEachObjectSurfaceMesh()
{
	/*for (int i = 0; i < allObjects.size(); i++)
	{
		for (int j = allObjects[i].objectSurfaceMeshes_node_start_end[0]; j < allObjects[i].objectSurfaceMeshes_node_start_end[1]; j++)
		{
			int localIndex = j - allObjects[i].objectSurfaceMeshes_node_start_end[0];
			allObjects[i].objectSurfaceMesh.vertices[localIndex] = pos_node_surface[j];
		}
	}*/
}