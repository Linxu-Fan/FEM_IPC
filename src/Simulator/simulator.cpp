#include "simulator.h"


void contact_Info::clear()
{
	PT.clear();
	EE.clear();
	PG.clear();


	PG_PG.clear();
	PT_PP.clear();
	PT_PE.clear();
	PT_PT.clear();
	EE_EE.clear();
}


// compute external force excluding gravity
Eigen::Vector3d compute_external_force(
	std::vector<boundaryCondition>& boundaryCondition_node, 
	int vertInd, 
	int timestep)
{
	Eigen::Vector3d extForce = Eigen::Vector3d::Zero();

	if (boundaryCondition_node[vertInd].type == 2)
	{
		if (timestep >= boundaryCondition_node[vertInd].appliedTime[0]
			&& timestep <= boundaryCondition_node[vertInd].appliedTime[1])
		{
			extForce = boundaryCondition_node[vertInd].force[timestep];
		}
	}

	return extForce;
}

double compute_Barrier_energy(FEMParamters& parameters,
	surface_Info& surfaceInfo,
	const std::vector<Eigen::Vector3d>& pos_node,
	const std::vector<Eigen::Vector3d>& pos_node_Rest,
	contact_Info& contact_pairs,
	const int timestep)
{
	double threshold2 = parameters.IPC_dis * parameters.IPC_dis;

	double barrierEnergy = 0;
	std::vector<double> energy_PT_PP(contact_pairs.PT_PP.size(), 0), energy_PT_PE(contact_pairs.PT_PE.size(), 0),
		energy_PT_PT(contact_pairs.PT_PT.size(), 0), energy_EE_EE(contact_pairs.EE_EE.size(), 0), energy_PG_PG(contact_pairs.PG_PG.size(), 0);
	if (contact_pairs.PG_PG.size() + contact_pairs.PT_PP.size() + contact_pairs.PT_PE.size() + contact_pairs.PT_PT.size() + contact_pairs.EE_EE.size() != 0)
	{
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < contact_pairs.PG_PG.size(); i++)
		{
			int ptInd = contact_pairs.PG_PG[i][1];
			Eigen::Vector3d P = pos_node[ptInd];
			double z2 = P[2] * P[2];
			if (z2 < threshold2)
			{
				energy_PG_PG[i] = Ground::val(z2, surfaceInfo.boundaryVertices_area[ptInd], parameters);

			}
		}
		barrierEnergy += std::accumulate(energy_PG_PG.begin(), energy_PG_PG.end(), 0.0);



#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < contact_pairs.PT_PP.size(); i++)
		{
			Vector5i cont_PT = contact_pairs.PT_PP[i];

			int ptInd = cont_PT[1];
			Eigen::Vector3d P = pos_node[ptInd];
			Eigen::Vector3i tri = surfaceInfo.boundaryTriangles[cont_PT[2]];
			Eigen::Vector3d A = pos_node[tri[0]];
			Eigen::Vector3d B = pos_node[tri[1]];
			Eigen::Vector3d C = pos_node[tri[2]];

			double dis2 = 0;
			DIS::computePointTriD(P, A, B, C, dis2, cont_PT[3]);
			if (dis2 < threshold2)
			{
				energy_PT_PP[i] = BarrierEnergy::val_PT(surfaceInfo.boundaryVertices_area[ptInd], dis2, parameters);
			}
		}
		barrierEnergy += std::accumulate(energy_PT_PP.begin(), energy_PT_PP.end(), 0.0);


#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < contact_pairs.PT_PE.size(); i++)
		{
			Vector5i cont_PT = contact_pairs.PT_PE[i];

			int ptInd = cont_PT[1];
			Eigen::Vector3d P = pos_node[ptInd];
			Eigen::Vector3i tri = surfaceInfo.boundaryTriangles[cont_PT[2]];
			Eigen::Vector3d A = pos_node[tri[0]];
			Eigen::Vector3d B = pos_node[tri[1]];
			Eigen::Vector3d C = pos_node[tri[2]];

			double dis2 = 0;
			DIS::computePointTriD(P, A, B, C, dis2, cont_PT[3]);
			if (dis2 < threshold2)
			{
				energy_PT_PE[i] = BarrierEnergy::val_PT(surfaceInfo.boundaryVertices_area[ptInd], dis2, parameters);
			}
		}
		barrierEnergy += std::accumulate(energy_PT_PE.begin(), energy_PT_PE.end(), 0.0);


#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < contact_pairs.PT_PT.size(); i++)
		{
			Vector5i cont_PT = contact_pairs.PT_PT[i];

			int ptInd = cont_PT[1];
			Eigen::Vector3d P = pos_node[ptInd];
			Eigen::Vector3i tri = surfaceInfo.boundaryTriangles[cont_PT[2]];
			Eigen::Vector3d A = pos_node[tri[0]];
			Eigen::Vector3d B = pos_node[tri[1]];
			Eigen::Vector3d C = pos_node[tri[2]];

			double dis2 = 0;
			DIS::computePointTriD(P, A, B, C, dis2, cont_PT[3]);
			if (dis2 < threshold2)
			{
				energy_PT_PT[i] = BarrierEnergy::val_PT(surfaceInfo.boundaryVertices_area[ptInd], dis2, parameters);
			}

		}
		barrierEnergy += std::accumulate(energy_PT_PT.begin(), energy_PT_PT.end(), 0.0);


#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < contact_pairs.EE_EE.size(); i++)
		{
			Vector5i cont_EE = contact_pairs.EE_EE[i];
			int E1 = cont_EE[1], E2 = cont_EE[2];
			int e1p1 = surfaceInfo.index_boundaryEdge[E1][0], e1p2 = surfaceInfo.index_boundaryEdge[E1][1],
				e2p1 = surfaceInfo.index_boundaryEdge[E2][0], e2p2 = surfaceInfo.index_boundaryEdge[E2][1];

			Eigen::Vector3d P1 = pos_node[e1p1];
			Eigen::Vector3d P2 = pos_node[e1p2];
			Eigen::Vector3d Q1 = pos_node[e2p1];
			Eigen::Vector3d Q2 = pos_node[e2p2];

			double dis2 = 0;
			DIS::computeEdgeEdgeD(P1, P2, Q1, Q2, dis2, cont_EE[3]);
			if (dis2 < threshold2)
			{
				Eigen::Vector4i ptIndices = { e1p1 , e1p2 , e2p1 , e2p2 };
				energy_EE_EE[i] = BarrierEnergy::val_EE(surfaceInfo.boundaryEdges_area[e1p1][e1p2],
					dis2, pos_node, pos_node_Rest, ptIndices, parameters);
			}
		}
		barrierEnergy += std::accumulate(energy_EE_EE.begin(), energy_EE_EE.end(), 0.0);



	}


	return barrierEnergy;
}



std::vector<std::pair<int, int>> find_contact_pair_BBX_level(
	const double& dilation,
	triMesh& triSimMesh,
	const std::vector<Vector12d>& moving_direction_ABD)
{

#pragma omp parallel for 
	for (int ii = 0; ii < triSimMesh.allObjects.size(); ii++)
	{
		Vector12d AFFINE;
		AFFINE.block(0, 0, 3, 1) = triSimMesh.translation_ABD[ii];
		AFFINE.block(3, 0, 3, 1) = triSimMesh.deformation_ABD[ii].col(0);
		AFFINE.block(6, 0, 3, 1) = triSimMesh.deformation_ABD[ii].col(1);
		AFFINE.block(9, 0, 3, 1) = triSimMesh.deformation_ABD[ii].col(2);


		triSimMesh.allObjects[ii].BBX.cal_min_max_ABD(AFFINE, dilation);
	}


	if (moving_direction_ABD.size() != 0)
	{
#pragma omp parallel for 
		for (int ii = 0; ii < triSimMesh.allObjects.size(); ii++)
		{

			Vector12d AFFINE;
			AFFINE.block(0, 0, 3, 1) = triSimMesh.translation_ABD[ii];
			AFFINE.block(3, 0, 3, 1) = triSimMesh.deformation_ABD[ii].col(0);
			AFFINE.block(6, 0, 3, 1) = triSimMesh.deformation_ABD[ii].col(1);
			AFFINE.block(9, 0, 3, 1) = triSimMesh.deformation_ABD[ii].col(2);
			AFFINE += moving_direction_ABD[ii];


			bounding_box BBX_increment = triSimMesh.allObjects[ii].BBX;
			BBX_increment.cal_min_max_ABD(AFFINE, dilation);

			triSimMesh.allObjects[ii].BBX.merges(BBX_increment);

		}
	}



	std::vector<std::pair<int, int>> contact_objects_BBX_level;
	for (int i = 0; i < triSimMesh.allObjects.size() - 1; i++)
	{
		for (int j = i + 1; j < triSimMesh.allObjects.size(); j++)
		{
			if (triSimMesh.allObjects[i].BBX.intersects(triSimMesh.allObjects[j].BBX))
			{
				contact_objects_BBX_level.push_back(std::make_pair(i, j));
			}
		}
	}

	return contact_objects_BBX_level;

}



std::vector<std::pair<int, int>> find_contact_pair_BVH_level(
	const double& dilation,
	const std::vector<std::pair<int, int>>& BBX_pair,
	triMesh& triSimMesh,
	std::vector<int>& BVH_objects,
	const std::vector<Eigen::Vector3d>& advection_direction)
{
	std::set<int> need_BVH_objects; // objects that need bounding box for a more detailed contact detection
	for (int i = 0; i < BBX_pair.size(); i++)
	{
		need_BVH_objects.insert(BBX_pair[i].first);
		need_BVH_objects.insert(BBX_pair[i].second);
	}

	BVH_objects.clear();
	BVH_objects.insert(BVH_objects.end(), need_BVH_objects.begin(), need_BVH_objects.end());
	if (advection_direction.size() == 0)
	{
#pragma omp parallel for 
		for (int k = 0; k < BVH_objects.size(); k++)
		{
			triSimMesh.build_BVH_object(dilation, BVH_objects[k]);
		}
	}
	else
	{
#pragma omp parallel for 
		for (int k = 0; k < BVH_objects.size(); k++)
		{
			//std::cout << "build advect BVH = "<< need_BVH_objects_vec[k] << "; advection_direction = " << advection_direction.size() << std::endl;
			triSimMesh.build_BVH_object_advect(dilation, advection_direction, BVH_objects[k]);
		}
	}

	// find contact pairs between objects
	std::vector<std::pair<int, int>> result;
	for (int i = 0; i < BBX_pair.size(); i++)
	{
		int obj1 = BBX_pair[i].first, obj2 = BBX_pair[i].second;
		bool intersect = triSimMesh.allObjects[obj1].object_BVH_faces->box.intersects(triSimMesh.allObjects[obj2].object_BVH_faces->box);
		if (intersect)
		{
			
			result.push_back(std::make_pair(obj1, obj2));
		}

	}


	return result;

}



void find_contact_pair_element_level(
	const std::vector<std::pair<int, int>>& BVH_pair,
	triMesh& triSimMesh,
	FEMParamters& parameters,
	contact_Info& contact_pairs)
{

	// find potential contact face pair
	std::vector<std::vector<std::pair<int, int>>> contact_pointTriangle_pair(BVH_pair.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int i = 0; i < BVH_pair.size(); i++)
	{
		int obj_1 = BVH_pair[i].first, obj_2 = BVH_pair[i].second;
		int obj_1_node_start = triSimMesh.allObjects[obj_1].objectSurfaceMeshes_node_start_end[0];
		int obj_2_face_start = triSimMesh.allObjects[obj_2].objectSurfaceMeshes_face_start_end[0];


		std::vector<std::pair<int, int>> cont_pair;

		// Contact between obj_1's point and obj_2's faces
		{
			std::vector<std::pair<int, int>> result_;
			queryBVH(triSimMesh.allObjects[obj_1].object_BVH_nodes, triSimMesh.allObjects[obj_2].object_BVH_faces, result_);

			// compute the actual contact
			for (int k = 0; k < result_.size(); k++)
			{
				int vert_global = result_[k].first + obj_1_node_start;
				int face_global = result_[k].second + obj_2_face_start;
				cont_pair.push_back(std::make_pair(vert_global, face_global));
			}
		}


		//// Contact between obj_1's faces and obj_2's points
		//{
		//	std::vector<std::pair<int, int>> result_;
		//	queryBVH(triSimMesh.allObjects[obj_1].object_BVH_faces, triSimMesh.allObjects[obj_2].object_BVH_nodes, result_);

		//	// compute the actual contact
		//	for (int k = 0; k < result_.size(); k++)
		//	{
		//		int vert_global = result_[k].first + obj_1_node_start;
		//		int face_global = result_[k].second + obj_2_face_start;
		//		cont_pair.push_back(std::make_pair(vert_global, face_global));
		//	}
		//}


		contact_pointTriangle_pair[i] = cont_pair;
	}


	// find potential contact edge pair
	std::vector<std::vector<std::pair<int, int>>> contact_edgeEdge_pair(BVH_pair.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int i = 0; i < BVH_pair.size(); i++)
	{
		int obj_1 = BVH_pair[i].first, obj_2 = BVH_pair[i].second;
		std::vector<std::pair<int, int>> result_;
		queryBVH(triSimMesh.allObjects[obj_1].object_BVH_edges, triSimMesh.allObjects[obj_2].object_BVH_edges, result_);
		// convert the local mesh info to global mesh info
		int obj_1_node_start = triSimMesh.allObjects[obj_1].objectSurfaceMeshes_node_start_end[0];
		int obj_2_node_start = triSimMesh.allObjects[obj_2].objectSurfaceMeshes_node_start_end[0];

		// compute the actual contact
		std::vector<std::pair<int, int>> cont_pair;
		for (int k = 0; k < result_.size(); k++)
		{
			Eigen::Vector2i edge1 = triSimMesh.allObjects[obj_1].objectSurfaceMesh.edges[result_[k].first] + Eigen::Vector2i::Ones() * obj_1_node_start;
			Eigen::Vector2i edge2 = triSimMesh.allObjects[obj_2].objectSurfaceMesh.edges[result_[k].second] + Eigen::Vector2i::Ones() * obj_2_node_start;

			int edge1_index = triSimMesh.surfaceInfo.boundaryEdge_index[edge1[0]][edge1[1]];
			int edge2_index = triSimMesh.surfaceInfo.boundaryEdge_index[edge2[0]][edge2[1]];

			cont_pair.push_back(std::make_pair(edge1_index, edge2_index));
		}

		contact_edgeEdge_pair[i] = cont_pair;
	}



	// find final actual contact pairs
	for (int i = 0; i < contact_pointTriangle_pair.size(); i++)
	{
		if (contact_pointTriangle_pair[i].size() != 0)
		{
			contact_pairs.PT.insert(contact_pairs.PT.end(), contact_pointTriangle_pair[i].begin(), contact_pointTriangle_pair[i].end());
		}
	}


	for (int i = 0; i < contact_edgeEdge_pair.size(); i++)
	{
		if (contact_edgeEdge_pair[i].size() != 0)
		{
			contact_pairs.EE.insert(contact_pairs.EE.end(), contact_edgeEdge_pair[i].begin(), contact_edgeEdge_pair[i].end());
		}
	}



	// Step 3: find ground contact pairs if any
	if (parameters.enableGround == true)
	{
		for (int ft = 0; ft < triSimMesh.surfaceInfo.boundaryVertices_vec.size(); ft++)
		{
			int ptInd = triSimMesh.surfaceInfo.boundaryVertices_vec[ft];
			if (triSimMesh.pos_node_surface[ptInd][2] <= parameters.IPC_dis)
			{
				contact_pairs.PG.push_back(ptInd);
			}
		}

	}


	std::cout << "			PT.size() = " << contact_pairs.PT.size();
	std::cout << "; EE.size() = " << contact_pairs.EE.size();
	std::cout << "; PG.size() = " << contact_pairs.PG.size() << std::endl;


}



void find_contact_pair_IPC_level(
	triMesh& triSimMesh,
	FEMParamters& parameters,
	contact_Info& contact_pairs)
{
	// find PT pair
	std::vector<bool> cal_PT(contact_pairs.PT.size());
	std::vector<Vector5i> ct_PT_pair(contact_pairs.PT.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int i = 0; i < contact_pairs.PT.size(); i++)
	{
		int vert_global = contact_pairs.PT[i].first;
		int face_global = contact_pairs.PT[i].second;

		Eigen::Vector3i triVerts = triSimMesh.surfaceInfo.boundaryTriangles[face_global];
		Eigen::Vector3d P = triSimMesh.pos_node_surface[vert_global];
		Eigen::Vector3d A = triSimMesh.pos_node_surface[triVerts[0]];
		Eigen::Vector3d B = triSimMesh.pos_node_surface[triVerts[1]];
		Eigen::Vector3d C = triSimMesh.pos_node_surface[triVerts[2]];


		int type = DIS::dType_PT(P, A, B, C);
		double dis2 = 0;
		DIS::computePointTriD(P, A, B, C, dis2);

		Vector5i ct = Vector5i::Zero();
		bool inside = false;
		if (dis2 <= squaredDouble(parameters.IPC_dis)) // only calculate the energy when the distance is smaller than the threshold
		{
			if (type <= 2)
			{
				ct = { 0 , vert_global , face_global , type, 2 };
			}
			else if (type > 2 && type <= 5)
			{
				ct = { 0 , vert_global , face_global , type, 3 };
			}
			else if (type == 6)
			{
				ct = { 0 , vert_global , face_global , type, 4 };
			}
			inside = true;
		}

		ct_PT_pair[i] = ct;
		cal_PT[i] = inside;
	}


	// find potential contact edge pair
	std::vector<bool> cal_EE(contact_pairs.EE.size());
	std::vector<Vector5i> ct_EE_pair(contact_pairs.EE.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int i = 0; i < contact_pairs.EE.size(); i++)
	{

		int edge1_index = contact_pairs.EE[i].first;
		int edge2_index = contact_pairs.EE[i].second;

		Eigen::Vector2i edge1 = triSimMesh.surfaceInfo.index_boundaryEdge[edge1_index];
		Eigen::Vector2i edge2 = triSimMesh.surfaceInfo.index_boundaryEdge[edge2_index];


		int P1I = edge1[0], P2I = edge1[1], Q1I = edge2[0], Q2I = edge2[1];
		Eigen::Vector3d P1 = triSimMesh.pos_node_surface[P1I];
		Eigen::Vector3d P2 = triSimMesh.pos_node_surface[P2I];
		Eigen::Vector3d Q1 = triSimMesh.pos_node_surface[Q1I];
		Eigen::Vector3d Q2 = triSimMesh.pos_node_surface[Q2I];

		int type = DIS::dType_EE(P1, P2, Q1, Q2);
		double dis2 = 0;
		DIS::computeEdgeEdgeD(P1, P2, Q1, Q2, dis2);

		Vector5i ct = Vector5i::Zero();
		bool inside = false;
		if (dis2 <= squaredDouble(parameters.IPC_dis)) // only calculate the energy when the distance is smaller than the threshold
		{
			ct = { 1 , edge1_index , edge2_index , type , 4 };
			inside = true;
		}
		ct_EE_pair[i] = ct;
		cal_EE[i] = inside;
	}



	// find final actual contact pairs
	for (int i = 0; i < cal_PT.size(); i++)
	{
		if (cal_PT[i])
		{
			if (ct_PT_pair[i][4] == 2)
			{
				contact_pairs.PT_PP.push_back(ct_PT_pair[i]);
			}
			else if (ct_PT_pair[i][4] == 3)
			{
				contact_pairs.PT_PE.push_back(ct_PT_pair[i]);
			}
			else
			{
				contact_pairs.PT_PT.push_back(ct_PT_pair[i]);
			}

		}
	}


	for (int i = 0; i < cal_EE.size(); i++)
	{
		if (cal_EE[i])
		{
			contact_pairs.EE_EE.push_back(ct_EE_pair[i]);
		}
	}



	// Step 3: find ground contact pairs if any
	if (parameters.enableGround == true)
	{
		for (int ft = 0; ft < contact_pairs.PG.size(); ft++)
		{
			Vector5i ct = { 2 , contact_pairs.PG[ft] , 0 , 0 , 0};
			contact_pairs.PG_PG.push_back(ct);
		}

	}


}



void find_contact(
	triMesh& triSimMesh,
	FEMParamters& parameters,
	contact_Info& contact_pairs,
	const std::vector<Vector12d>& moving_direction_ABD,
	const std::vector<Eigen::Vector3d>& moving_direction_pos)
{
	contact_pairs.clear();
	std::vector<int> BVH_objects;

	if (moving_direction_ABD.size() != 0)
	{

		if (moving_direction_pos.size() != 0) // detect the contact before calculating the maximum step size
		{
			std::vector<std::pair<int, int>> BBX_contact = find_contact_pair_BBX_level(parameters.IPC_dis / 2.0, triSimMesh, moving_direction_ABD);
			if (BBX_contact.size() != 0 || parameters.enableGround == true)
			{
				std::vector<std::pair<int, int>> BVH_contact = find_contact_pair_BVH_level(parameters.IPC_dis / 2.0, BBX_contact, triSimMesh, BVH_objects, moving_direction_pos);
				if (BVH_contact.size() != 0 || parameters.enableGround == true)
				{
					find_contact_pair_element_level(BVH_contact, triSimMesh, parameters, contact_pairs);


				}
			}
						
		}
	}
	else  // detect the contact at the begining of a timestep
	{
		std::vector<std::pair<int, int>> BBX_contact = find_contact_pair_BBX_level(parameters.IPC_dis / 2.0, triSimMesh);
		if (BBX_contact.size() != 0 || parameters.enableGround == true)
		{
			std::vector<std::pair<int, int>> BVH_contact = find_contact_pair_BVH_level(parameters.IPC_dis / 2.0, BBX_contact, triSimMesh, BVH_objects);
			if (BVH_contact.size() != 0 || parameters.enableGround == true)
			{
				find_contact_pair_element_level(BVH_contact, triSimMesh, parameters,contact_pairs);
				find_contact_pair_IPC_level(triSimMesh, parameters, contact_pairs);
			}
		}

		
	}



#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int k = 0; k < BVH_objects.size(); k++)
	{
		int index = BVH_objects[k];

		deleteBVH(triSimMesh.allObjects[index].object_BVH_nodes);
		triSimMesh.allObjects[index].object_BVH_nodes = nullptr;


		deleteBVH(triSimMesh.allObjects[index].object_BVH_edges);
		triSimMesh.allObjects[index].object_BVH_edges = nullptr;

		deleteBVH(triSimMesh.allObjects[index].object_BVH_faces);
		triSimMesh.allObjects[index].object_BVH_faces = nullptr;
	

	}


}




// calculate the maximum feasible step size
double calMaxStepSize(
	FEMParamters& parameters,
	surface_Info& surfaceInfo,
	const std::vector<Eigen::Vector3d>& direction,
	const std::vector<Eigen::Vector3d>& pos_node,
	const std::vector<std::string>& note_node,
	std::map<std::string, int>& tetMeshIndex,
	const int timestep)
{
	////std::cout << "Calculating maximum step!" << std::endl;
	//std::set<int> culledSet; // vertices who are in a contact

	//// Step 1: calculate the culled constraint
	//for (int i  = 0; i < pTeEBarrVec.PT_Indices.size(); i++)
	//{
	//	culledSet.insert(pTeEBarrVec.PT_Indices[i].begin(), pTeEBarrVec.PT_Indices[i].end());
	//}
	//for (int i = 0; i < pTeEBarrVec.EE_Indices.size(); i++)
	//{
	//	culledSet.insert(pTeEBarrVec.EE_Indices[i].begin(), pTeEBarrVec.EE_Indices[i].end());
	//}
	////std::cout << "culledSet.size() = " << culledSet.size() << std::endl;
	//
	//// Step 2: calculate alpha_F
	//double alpha_F = 1.0; // Eq.3 in IPC paper's supplementary document
	//for (int i = 0; i < direction.size(); i++)
	//{
	//	if (culledSet.find(i) == culledSet.end())
	//	{
	//		alpha_F = std::min(alpha_F, parameters.IPC_dis / 2.0 / direction[i].norm());
	//	}
	//}
	////std::cout << "alpha_F = " << alpha_F <<"; pTeEBarrVec.PT_Indices.size() = "<< pTeEBarrVec.PT_Indices.size() << std::endl;

	//// Step 3: calculate alpha_C_hat
	//double alpha_C_hat = 1.0;
	//for (int i = 0; i < pTeEBarrVec.PT_Indices.size(); i++)
	//{
	//	Eigen::Vector4i PP_Index = pTeEBarrVec.PT_Indices[i];
	//	int p0 = PP_Index[0], p1 = PP_Index[1], p2 = PP_Index[2], p3 = PP_Index[3];
	//	//std::cout << "	PT " << i << " ";
	//	bool intersect = pointTriangleCCDBroadphase(tetSimMesh.pos_node[p0], direction[p0], tetSimMesh.pos_node[p1], 
	//		direction[p1], tetSimMesh.pos_node[p2], direction[p2], tetSimMesh.pos_node[p3], direction[p3], parameters.IPC_dis);
	//	//std::cout << "	intersect = " << intersect << std::endl;
	//	if (intersect)
	//	{
	//		//std::cout << "		tetSimMesh.pos_node[p0] = "<< tetSimMesh.pos_node[p0] << std::endl;
	//		//std::cout << "		tetSimMesh.pos_node[p1] = "<< tetSimMesh.pos_node[p1] << std::endl;
	//		//std::cout << "		tetSimMesh.pos_node[p2] = "<< tetSimMesh.pos_node[p2] << std::endl;
	//		//std::cout << "		tetSimMesh.pos_node[p3] = "<< tetSimMesh.pos_node[p3] << std::endl;
	//		//std::cout << "		direction[p0] = "<< direction[p0] << std::endl;
	//		//std::cout << "		direction[p1] = "<< direction[p1] << std::endl;
	//		//std::cout << "		direction[p2] = "<< direction[p2] << std::endl;
	//		//std::cout << "		direction[p3] = "<< direction[p3] << std::endl;
	//		double alpha_tmp = pointTriangleCCDNarrowphase(tetSimMesh.pos_node[p0], direction[p0], tetSimMesh.pos_node[p1],
	//			direction[p1], tetSimMesh.pos_node[p2], direction[p2], tetSimMesh.pos_node[p3], direction[p3], parameters.IPC_eta);
	//		//std::cout << "		xxx2" << std::endl;
	//		alpha_C_hat = std::min(alpha_C_hat, alpha_tmp);
	//	}
	//	
	//}
	//for (int i = 0; i < pTeEBarrVec.EE_Indices.size(); i++)
	//{
	//	Eigen::Vector4i PP_Index = pTeEBarrVec.EE_Indices[i];
	//	int p0 = PP_Index[0], p1 = PP_Index[1], p2 = PP_Index[2], p3 = PP_Index[3];
	//	//std::cout << "	EE " << i << " ";
	//	bool intersect = edgeEdgeCCDBroadphase(tetSimMesh.pos_node[p0], direction[p0], tetSimMesh.pos_node[p1],
	//		direction[p1], tetSimMesh.pos_node[p2], direction[p2], tetSimMesh.pos_node[p3], direction[p3], parameters.IPC_dis);
	//	if (intersect)
	//	{
	//		double alpha_tmp = edgeEdgeCCDNarrowphase(tetSimMesh.pos_node[p0], direction[p0], tetSimMesh.pos_node[p1],
	//			direction[p1], tetSimMesh.pos_node[p2], direction[p2], tetSimMesh.pos_node[p3], direction[p3], parameters.IPC_eta);
	//		alpha_C_hat = std::min(alpha_C_hat, alpha_tmp);
	//	}
	//	
	//}
	////std::cout << "alpha_C_hat = " << alpha_C_hat << std::endl;

	//// Step 4: calculate full CCD if necessary
	////if (alpha_F >= 0.5 * alpha_C_hat)
	//if(0)
	//{
	//	//std::cout << "Partial calculation!" << std::endl;
	//	if (parameters.enableGround == true)
	//	{
	//		double stepGround = 1.0;
	//		for (std::map<int, std::set<int>>::iterator it = tetSimMesh.boundaryVertices.begin(); it != tetSimMesh.boundaryVertices.end(); it++)
	//		{
	//			int ptInd = it->first;
	//			if (direction[ptInd][2] < 0)
	//			{
	//				double coor_z = tetSimMesh.pos_node[ptInd][2];
	//				stepGround = std::min(stepGround, coor_z * (1.0 - parameters.IPC_eta) / std::abs(direction[ptInd][2]));
	//			}
	//			
	//		}
	//		alpha_C_hat = std::min(stepGround, alpha_C_hat);
	//	}
	//	return std::min(alpha_F, alpha_C_hat);
	//}
	//else // use spatial hash to calculate the actual full CCD
	{
		//double CCD_step = calMaxStep_spatialHash(
		//	parameters,
		//	surfaceInfo,
		//	direction,
		//	pos_node,
		//	note_node,
		//	tetMeshIndex,
		//	timestep);
		


		return 0;
	}

}


