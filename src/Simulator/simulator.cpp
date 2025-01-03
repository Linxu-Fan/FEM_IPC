#include "simulator.h"


void contact_Info::clear()
{
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
			energy_PG_PG[i] = Ground::val(z2, surfaceInfo.boundaryVertices_area[ptInd], parameters);		
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
			energy_PT_PP[i] = BarrierEnergy::val_PT(surfaceInfo.boundaryVertices_area[ptInd], dis2, parameters);
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
			energy_PT_PE[i] = BarrierEnergy::val_PT(surfaceInfo.boundaryVertices_area[ptInd], dis2, parameters);

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
			energy_PT_PT[i] = BarrierEnergy::val_PT(surfaceInfo.boundaryVertices_area[ptInd], dis2, parameters);

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
			Eigen::Vector4i ptIndices = { e1p1 , e1p2 , e2p1 , e2p2 };
			energy_EE_EE[i] = BarrierEnergy::val_EE(surfaceInfo.boundaryEdges_area[e1p1][e1p2],
				dis2, pos_node, pos_node_Rest, ptIndices, parameters);

		}
		barrierEnergy += std::accumulate(energy_EE_EE.begin(), energy_EE_EE.end(), 0.0);



	}


	return barrierEnergy;
}

void calContactInfo(
	triMesh& triSimMesh,
	FEMParamters& parameters,
	surface_Info& surfaceInfo,
	const std::vector<Eigen::Vector3d>& pos_node,
	const std::vector<std::string>& note_node,
	std::map<std::string, int>& tetMeshIndex,
	const int timestep,
	contact_Info& contact_pairs)
{

	double startTime1, endTime1;
	startTime1 = omp_get_wtime();

	contact_pairs.clear();

	// update the BVH of each object
	triSimMesh.build_BVH_object(parameters.IPC_dis / 2.0);


	// find contact pairs between objects
	std::vector<std::pair<int, int>> contact_objects;
	for (int i = 0; i < triSimMesh.allObjects.size() - 1; i++)
	{
		for (int j = i + 1; j < triSimMesh.allObjects.size(); j++)
		{
			bool intersect = triSimMesh.allObjects[i].object_BVH_faces->box.intersects(triSimMesh.allObjects[j].object_BVH_faces->box);
			if (intersect)
			{
				contact_objects.push_back(std::make_pair(i, j));
			}
		}
	}


	// find potential contact face pair
	std::vector<std::vector<Vector5i>> contact_pointTriangle_pair(contact_objects.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int i = 0; i < contact_objects.size(); i++)
	{
		int obj_1 = contact_objects[i].first, obj_2 = contact_objects[i].second;
		std::vector<std::pair<int, int>> result_;
		queryBVH(triSimMesh.allObjects[obj_1].object_BVH_nodes, triSimMesh.allObjects[obj_2].object_BVH_faces, result_);
		// convert the local mesh info to global mesh info
		int obj_1_node_start = triSimMesh.allObjects[obj_1].objectSurfaceMeshes_node_start_end[0];
		int obj_2_face_start = triSimMesh.allObjects[obj_2].objectSurfaceMeshes_face_start_end[0];

		// compute the actual contact
		std::vector<Vector5i> cont_pair;
		for (int k = 0; k < result_.size(); k++)
		{
			int vert_global = result_[k].first + obj_1_node_start;
			int face_global = result_[k].second + obj_2_face_start;

			Eigen::Vector3i triVerts = surfaceInfo.boundaryTriangles[face_global];
			Eigen::Vector3d P = pos_node[vert_global];
			Eigen::Vector3d A = pos_node[triVerts[0]];
			Eigen::Vector3d B = pos_node[triVerts[1]];
			Eigen::Vector3d C = pos_node[triVerts[2]];


			int type = DIS::dType_PT(P, A, B, C);
			double dis2 = 0;
			DIS::computePointTriD(P, A, B, C, dis2);

			if (dis2 <= squaredDouble(parameters.IPC_dis)) // only calculate the energy when the distance is smaller than the threshold
			{
				if (type <= 2)
				{
					Vector5i ct = { 0 , vert_global , face_global , type, 2 };
					cont_pair.push_back(ct);
				}
				else if (type > 2 && type <= 5)
				{
					Vector5i ct = { 0 , vert_global , face_global , type, 3 };
					cont_pair.push_back(ct);
				}
				else if (type == 6)
				{
					Vector5i ct = { 0 , vert_global , face_global , type, 4 };
					cont_pair.push_back(ct);
				}
			}
			

		}

		contact_pointTriangle_pair[i] = cont_pair;
	}


	// find potential contact edge pair
	std::vector<std::vector<Vector5i>> contact_edgeEdge_pair(contact_objects.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int i = 0; i < contact_objects.size(); i++)
	{
		int obj_1 = contact_objects[i].first, obj_2 = contact_objects[i].second;
		std::vector<std::pair<int, int>> result_;
		queryBVH(triSimMesh.allObjects[obj_1].object_BVH_edges, triSimMesh.allObjects[obj_2].object_BVH_edges, result_);
		// convert the local mesh info to global mesh info
		int obj_1_node_start = triSimMesh.allObjects[obj_1].objectSurfaceMeshes_node_start_end[0];
		int obj_2_node_start = triSimMesh.allObjects[obj_2].objectSurfaceMeshes_node_start_end[0];

		// compute the actual contact
		std::vector<Vector5i> cont_pair;
		for (int k = 0; k < result_.size(); k++)
		{
			Eigen::Vector2i edge1 = triSimMesh.allObjects[obj_1].objectSurfaceMesh.edges[result_[k].first] + Eigen::Vector2i::Ones() * obj_1_node_start;
			Eigen::Vector2i edge2 = triSimMesh.allObjects[obj_2].objectSurfaceMesh.edges[result_[k].second] + Eigen::Vector2i::Ones() * obj_2_node_start;

			int edge1_index = triSimMesh.surfaceInfo.boundaryEdge_index[edge1[0]][edge1[1]];
			int edge2_index = triSimMesh.surfaceInfo.boundaryEdge_index[edge2[0]][edge2[1]];

			int P1I = edge1[0], P2I = edge1[1], Q1I = edge2[0], Q2I = edge2[1];
			Eigen::Vector3d P1 = pos_node[P1I];
			Eigen::Vector3d P2 = pos_node[P2I];
			Eigen::Vector3d Q1 = pos_node[Q1I];
			Eigen::Vector3d Q2 = pos_node[Q2I];

			int type = DIS::dType_EE(P1, P2, Q1, Q2);
			double dis2 = 0;
			DIS::computeEdgeEdgeD(P1, P2, Q1, Q2, dis2);

			if (dis2 <= squaredDouble(parameters.IPC_dis)) // only calculate the energy when the distance is smaller than the threshold
			{
				Vector5i ct = { 1 , edge1_index , edge2_index , type , 4 };
				cont_pair.push_back(ct);
			}
		
		}

		contact_edgeEdge_pair[i] = cont_pair;
	}



	// find final actual contact pairs
	for (int i = 0; i < contact_pointTriangle_pair.size(); i++)
	{
		if (contact_pointTriangle_pair[i].size() != 0)
		{
			for (int j = 0; j < contact_pointTriangle_pair[i].size(); j++)
			{
				if (contact_pointTriangle_pair[i][j][4] == 2)
				{
					contact_pairs.PT_PP.push_back(contact_pointTriangle_pair[i][j]);
				}
				else if (contact_pointTriangle_pair[i][j][4] == 3)
				{
					contact_pairs.PT_PE.push_back(contact_pointTriangle_pair[i][j]);
				}
				else
				{
					contact_pairs.PT_PT.push_back(contact_pointTriangle_pair[i][j]);
				}
			}
		}		
	}


	for (int i = 0; i < contact_edgeEdge_pair.size(); i++)
	{
		if (contact_edgeEdge_pair[i].size() != 0)
		{
			for (int j = 0; j < contact_edgeEdge_pair[i].size(); j++)
			{
				contact_pairs.EE_EE.push_back(contact_edgeEdge_pair[i][j]);
			}
		}
	}


	// Step 3: find ground contact pairs if any
	if (parameters.enableGround == true)
	{
		for (int ft = 0; ft < surfaceInfo.boundaryVertices_vec.size(); ft++)
		{
			int ptInd = surfaceInfo.boundaryVertices_vec[ft];
			if (pos_node[ptInd][2] <= parameters.IPC_dis)
			{
				Vector5i ct = { 2 , ptInd , 0 , 0 , 0 };
				contact_pairs.PG_PG.push_back(ct);
			}
		}

	}


	std::cout << "			PT_PP.size() = " << contact_pairs.PT_PP.size();
	std::cout << "; PT_PE.size() = " << contact_pairs.PT_PE.size();
	std::cout << "; PT_PT.size() = " << contact_pairs.PT_PT.size();
	std::cout << "; EE_EE.size() = " << contact_pairs.EE_EE.size() << std::endl;


	//if (timestep == 87)
	//{
	//	for (int h = 0; h < triSimMesh.allObjects.size(); h++)
	//	{
	//		triSimMesh.allObjects[h].object_BVH_faces->export_bounding_box_mesh(true, triSimMesh.allObjects[h].objectNote);
	//	}
	//}

	endTime1 = omp_get_wtime();
	std::cout << "			Contact Pairs Time is : " << endTime1 - startTime1 << "s" << std::endl;

}




void calContactInfo_advect(
	triMesh& triSimMesh,
	FEMParamters& parameters,
	surface_Info& surfaceInfo,
	const std::vector<Eigen::Vector3d>& pos_node,
	const std::vector<std::string>& note_node,
	std::map<std::string, int>& tetMeshIndex,
	const int timestep,
	contact_Info& contact_pairs,
	const std::vector<Eigen::Vector3d>& moving_direction)
{

	double startTime1, endTime1;
	startTime1 = omp_get_wtime();

	contact_pairs.clear();

	// update the BVH of each object
	triSimMesh.build_BVH_object_advect(parameters.IPC_dis / 2.0, moving_direction);



	// find contact pairs between objects
	std::vector<std::pair<int, int>> contact_objects;
	for (int i = 0; i < triSimMesh.allObjects.size() - 1; i++)
	{
		for (int j = i + 1; j < triSimMesh.allObjects.size(); j++)
		{
			bool intersect = triSimMesh.allObjects[i].object_BVH_faces->box.intersects(triSimMesh.allObjects[j].object_BVH_faces->box);
			if (intersect)
			{
				contact_objects.push_back(std::make_pair(i, j));
			}
		}
	}


	// find potential contact face pair
	std::vector<std::vector<Vector5i>> contact_pointTriangle_pair(contact_objects.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int i = 0; i < contact_objects.size(); i++)
	{
		int obj_1 = contact_objects[i].first, obj_2 = contact_objects[i].second;
		std::vector<std::pair<int, int>> result_;
		queryBVH(triSimMesh.allObjects[obj_1].object_BVH_nodes, triSimMesh.allObjects[obj_2].object_BVH_faces, result_);
		// convert the local mesh info to global mesh info
		int obj_1_node_start = triSimMesh.allObjects[obj_1].objectSurfaceMeshes_node_start_end[0];
		int obj_2_face_start = triSimMesh.allObjects[obj_2].objectSurfaceMeshes_face_start_end[0];

		// compute the actual contact
		std::vector<Vector5i> cont_pair;
		for (int k = 0; k < result_.size(); k++)
		{
			int vert_global = result_[k].first + obj_1_node_start;
			int face_global = result_[k].second + obj_2_face_start;

			Vector5i ct = { 0 , vert_global , face_global , -99, -99 };
			cont_pair.push_back(ct);
		}

		contact_pointTriangle_pair[i] = cont_pair;
	}


	// find potential contact edge pair
	std::vector<std::vector<Vector5i>> contact_edgeEdge_pair(contact_objects.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int i = 0; i < contact_objects.size(); i++)
	{
		int obj_1 = contact_objects[i].first, obj_2 = contact_objects[i].second;
		std::vector<std::pair<int, int>> result_;
		queryBVH(triSimMesh.allObjects[obj_1].object_BVH_edges, triSimMesh.allObjects[obj_2].object_BVH_edges, result_);
		// convert the local mesh info to global mesh info
		int obj_1_node_start = triSimMesh.allObjects[obj_1].objectSurfaceMeshes_node_start_end[0];
		int obj_2_node_start = triSimMesh.allObjects[obj_2].objectSurfaceMeshes_node_start_end[0];

		// compute the actual contact
		std::vector<Vector5i> cont_pair;
		for (int k = 0; k < result_.size(); k++)
		{
			Eigen::Vector2i edge1 = triSimMesh.allObjects[obj_1].objectSurfaceMesh.edges[result_[k].first] + Eigen::Vector2i::Ones() * obj_1_node_start;
			Eigen::Vector2i edge2 = triSimMesh.allObjects[obj_2].objectSurfaceMesh.edges[result_[k].second] + Eigen::Vector2i::Ones() * obj_2_node_start;

			int edge1_index = triSimMesh.surfaceInfo.boundaryEdge_index[edge1[0]][edge1[1]];
			int edge2_index = triSimMesh.surfaceInfo.boundaryEdge_index[edge2[0]][edge2[1]];

			Vector5i ct = { 1 , edge1_index , edge2_index , -99 , -99 };
			cont_pair.push_back(ct);
		}

		contact_edgeEdge_pair[i] = cont_pair;
	}



	// find final actual contact pairs
	for (int i = 0; i < contact_pointTriangle_pair.size(); i++)
	{
		if (contact_pointTriangle_pair[i].size() != 0)
		{
			for (int j = 0; j < contact_pointTriangle_pair[i].size(); j++)
			{
				if (contact_pointTriangle_pair[i][j][4] == 2)
				{
					contact_pairs.PT_PP.push_back(contact_pointTriangle_pair[i][j]);
				}
				else if (contact_pointTriangle_pair[i][j][4] == 3)
				{
					contact_pairs.PT_PE.push_back(contact_pointTriangle_pair[i][j]);
				}
				else
				{
					contact_pairs.PT_PT.push_back(contact_pointTriangle_pair[i][j]);
				}
			}
		}
	}


	for (int i = 0; i < contact_edgeEdge_pair.size(); i++)
	{
		if (contact_edgeEdge_pair[i].size() != 0)
		{
			for (int j = 0; j < contact_edgeEdge_pair[i].size(); j++)
			{
				contact_pairs.EE_EE.push_back(contact_edgeEdge_pair[i][j]);
			}
		}
	}


	// Step 3: find ground contact pairs if any
	if (parameters.enableGround == true)
	{
		for (int ft = 0; ft < surfaceInfo.boundaryVertices_vec.size(); ft++)
		{
			int ptInd = surfaceInfo.boundaryVertices_vec[ft];
			if (pos_node[ptInd][2] <= parameters.IPC_dis)
			{
				Vector5i ct = { 2 , ptInd , 0 , 0 , 0 };
				contact_pairs.PG_PG.push_back(ct);
			}
		}

	}


	std::cout << "			PT_PP.size() = " << contact_pairs.PT_PP.size();
	std::cout << "; PT_PE.size() = " << contact_pairs.PT_PE.size();
	std::cout << "; PT_PT.size() = " << contact_pairs.PT_PT.size();
	std::cout << "; EE_EE.size() = " << contact_pairs.EE_EE.size() << std::endl;



	endTime1 = omp_get_wtime();
	std::cout << "			Contact Pairs Time is : " << endTime1 - startTime1 << "s" << std::endl;

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


