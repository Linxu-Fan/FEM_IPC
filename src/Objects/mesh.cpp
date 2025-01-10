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
		if (allObjects[ii].need_update)
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
	reset_object_state();
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
	update_box_corner();
	//!!!!!!!!!!!!!!!!!
	//! Note update_position_world_after_cutting must be after update_box_corner
	update_position_world_after_cutting();
	reset_object_state();


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


		obj_.objectSurfaceMesh.updateBEInfo();
		obj_.pos_node_surface = obj_.objectSurfaceMesh.vertices;
		obj_.pos_node_surface_prev = obj_.objectSurfaceMesh.vertices;
		obj_.pos_node_surface_direction.resize(obj_.pos_node_surface.size(), Eigen::Vector3d::Zero());

		obj_.affine[3] = 1.0;
		obj_.affine[7] = 1.0;
		obj_.affine[11] = 1.0;

		obj_.affine_prev = obj_.affine;
		obj_.affine_last = obj_.affine;
		obj_.affine_vel[0] = config.velocity[0];
		obj_.affine_vel[1] = config.velocity[1];
		obj_.affine_vel[2] = config.velocity[2];

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
		if(allObjects[i].need_update)
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

}

void triMesh::update_ABD_info()
{
	for (int i = 0; i < allObjects.size(); i++)
	{
		if (allObjects[i].need_update)
		{
			for (int j = 0; j < allObjects[i].pos_node_interior.size(); j++)
			{
				Eigen::Matrix<double, 3, 12> Jx = build_Jx_matrix_for_ABD(allObjects[i].pos_node_interior[j]);

				allObjects[i].massMatrix_ABD += allObjects[i].mass_node_interior[j] * Jx.transpose() * Jx;
			}
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

void triMesh::reset_object_state()
{
	for (int i = 0; i < allObjects.size(); i++)
	{
		allObjects[i].need_update = false;
	}

}

void triMesh::update_position_world_after_cutting()
{
	for (int i = 0; i < allObjects.size(); i++)
	{
		if (allObjects[i].need_update)
		{
			// update interior
#pragma omp parallel for 
			for (int iv = 0; iv < allObjects[i].pos_node_interior.size(); iv++)
			{
				Eigen::Matrix<double, 3, 12> Jx = build_Jx_matrix_for_ABD(allObjects[i].pos_node_interior[iv]);
				allObjects[i].pos_node_interior[iv] = Jx * allObjects[i].affine;
			}

			// update surface
#pragma omp parallel for 
			for (int sv = 0; sv < allObjects[i].pos_node_surface.size(); sv++)
			{
				Eigen::Matrix<double, 3, 12> Jx = build_Jx_matrix_for_ABD(allObjects[i].pos_node_surface[sv]);
				allObjects[i].pos_node_surface[sv] = Jx * allObjects[i].affine;
			}
		}
	}

}
