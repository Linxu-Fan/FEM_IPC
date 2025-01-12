#include "simulator.h"


void contact_Info::clear()
{
	//PT.clear();
	//EE.clear();
	Point_Ground.clear();


	//PG_PG.clear();
	//PT_PP.clear();
	//PT_PE.clear();
	//PT_PT.clear();
	//EE_EE.clear();


	Point_Triangle.clear();
	Point_Triangle_PP_Dis.clear();
	Point_Triangle_PE_Dis.clear();
	Point_Triangle_PT_Dis.clear();
	Point_Triangle_PP_index.clear();
	Point_Triangle_PE_index.clear();
	Point_Triangle_PT_index.clear();

	Edge_Edge.clear();
	Edge_Edge_Dis.clear();

}



double compute_Barrier_energy(
	triMesh& triSimMesh,
	FEMParamters& parameters,
	contact_Info& contact_pairs,
	const int timestep)
{
	double threshold2 = parameters.IPC_dis * parameters.IPC_dis;

	double barrierEnergy = 0;
	if (contact_pairs.Point_Triangle_PP_index.size() + contact_pairs.Point_Triangle_PE_index.size() + 
		contact_pairs.Point_Triangle_PT_index.size() + contact_pairs.Edge_Edge.size() + contact_pairs.Point_Ground.size() != 0)
	{
		std::vector<double> energy_PT_PP(contact_pairs.Point_Triangle_PP_index.size(), 0), energy_PT_PE(contact_pairs.Point_Triangle_PE_index.size(), 0), energy_PT_PT(contact_pairs.Point_Triangle_PT_index.size(), 0);
		std::vector<double> energy_EE(contact_pairs.Edge_Edge.size(), 0), energy_PG(contact_pairs.Point_Ground.size(), 0);


#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < contact_pairs.Point_Ground.size(); i++)
		{
			Vector2i ct = contact_pairs.Point_Ground[i];
			Eigen::Vector3d P = triSimMesh.allObjects[ct[0]].pos_node_surface[ct[1]];
			double z2 = P[2] * P[2];
			if (z2 < threshold2)
			{
				energy_PG[i] = Ground::val(z2, triSimMesh.allObjects[ct[0]].objectSurfaceMesh.boundaryVertices_area[ct[1]], parameters);

			}
		}
		barrierEnergy += std::accumulate(energy_PG.begin(), energy_PG.end(), 0.0);



#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < contact_pairs.Point_Triangle_PP_index.size(); i++)
		{
			int ppindex = contact_pairs.Point_Triangle_PP_index[i];
			Vector5i cont_PT = contact_pairs.Point_Triangle[ppindex];

			double dis2 = contact_pairs.Point_Triangle_PP_Dis[i];

			int obj_1 = cont_PT[0], vert = cont_PT[1];

			double surface_area = triSimMesh.allObjects[obj_1].objectSurfaceMesh.boundaryVertices_area[vert];
			energy_PT_PP[i] = BarrierEnergy::val_PT(surface_area, dis2, parameters);

		}
		barrierEnergy += std::accumulate(energy_PT_PP.begin(), energy_PT_PP.end(), 0.0);

#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < contact_pairs.Point_Triangle_PE_index.size(); i++)
		{
			int ppindex = contact_pairs.Point_Triangle_PE_index[i];
			Vector5i cont_PT = contact_pairs.Point_Triangle[ppindex];

			double dis2 = contact_pairs.Point_Triangle_PE_Dis[i];

			int obj_1 = cont_PT[0], vert = cont_PT[1];

			double surface_area = triSimMesh.allObjects[obj_1].objectSurfaceMesh.boundaryVertices_area[vert];
			energy_PT_PE[i] = BarrierEnergy::val_PT(surface_area, dis2, parameters);

		}
		barrierEnergy += std::accumulate(energy_PT_PE.begin(), energy_PT_PE.end(), 0.0);

#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < contact_pairs.Point_Triangle_PT_index.size(); i++)
		{
			int ppindex = contact_pairs.Point_Triangle_PT_index[i];
			Vector5i cont_PT = contact_pairs.Point_Triangle[ppindex];

			double dis2 = contact_pairs.Point_Triangle_PT_Dis[i];

			int obj_1 = cont_PT[0], vert = cont_PT[1];

			double surface_area = triSimMesh.allObjects[obj_1].objectSurfaceMesh.boundaryVertices_area[vert];
			energy_PT_PT[i] = BarrierEnergy::val_PT(surface_area, dis2, parameters);

		}
		barrierEnergy += std::accumulate(energy_PT_PT.begin(), energy_PT_PT.end(), 0.0);

		//std::cout << "	Energy after PT contact = " << barrierEnergy << std::endl;


		//int count = 0;
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < contact_pairs.Edge_Edge.size(); i++)
		{
			Vector5i cont_PT = contact_pairs.Edge_Edge[i];
			{
				int obj_1 = cont_PT[0], obj_2 = cont_PT[2];
				int edge1 = cont_PT[1], edge2 = cont_PT[3];

		
				Eigen::Vector2i edge1_vert = triSimMesh.allObjects[obj_1].objectSurfaceMesh.index_boundaryEdge[edge1];
				Eigen::Vector2i edge2_vert = triSimMesh.allObjects[obj_2].objectSurfaceMesh.index_boundaryEdge[edge2];


				int P1I = edge1_vert[0], P2I = edge1_vert[1], Q1I = edge2_vert[0], Q2I = edge2_vert[1];
				double surface_area = triSimMesh.allObjects[obj_1].objectSurfaceMesh.boundaryEdges_area[P1I][P2I];

				Eigen::Vector3d P1 = triSimMesh.allObjects[obj_1].pos_node_surface[P1I];
				Eigen::Vector3d P2 = triSimMesh.allObjects[obj_1].pos_node_surface[P2I];
				Eigen::Vector3d Q1 = triSimMesh.allObjects[obj_2].pos_node_surface[Q1I];
				Eigen::Vector3d Q2 = triSimMesh.allObjects[obj_2].pos_node_surface[Q2I];

				Eigen::Vector3d P1_Rest = triSimMesh.allObjects[obj_1].objectSurfaceMesh.vertices[P1I];
				Eigen::Vector3d P2_Rest = triSimMesh.allObjects[obj_1].objectSurfaceMesh.vertices[P2I];
				Eigen::Vector3d Q1_Rest = triSimMesh.allObjects[obj_2].objectSurfaceMesh.vertices[Q1I];
				Eigen::Vector3d Q2_Rest = triSimMesh.allObjects[obj_2].objectSurfaceMesh.vertices[Q2I];


				double dis2 = contact_pairs.Edge_Edge_Dis[i];

				double eps_x = DIS::cal_EEM_eps_x(P1_Rest, P2_Rest, Q1_Rest, Q2_Rest);
				double val_ek = 0;
				DIS::compute_e(P1, P2, Q1, Q2, eps_x, val_ek);


				energy_EE[i] = BarrierEnergy::val_EE(surface_area,
					dis2, val_ek, parameters);

				//count += 1;
	
			}

		}
		barrierEnergy += std::accumulate(energy_EE.begin(), energy_EE.end(), 0.0);

		//std::cout << "	Energy after EE contact = " << barrierEnergy << "; Actual EE = " << count << std::endl;

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
		triSimMesh.allObjects[ii].BBX.cal_min_max_ABD(triSimMesh.allObjects[ii].affine, dilation);


		//triSimMesh.allObjects[ii].BBX.export_BBX_world(triSimMesh.allObjects[ii].objectNote+"_" + std::to_string(ii));
	}


	if (moving_direction_ABD.size() != 0)
	{
#pragma omp parallel for 
		for (int ii = 0; ii < triSimMesh.allObjects.size(); ii++)
		{
			bounding_box BBX_increment = triSimMesh.allObjects[ii].BBX;
			BBX_increment.cal_min_max_ABD(triSimMesh.allObjects[ii].affine + moving_direction_ABD[ii], dilation);
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
	bool advection)
{
	std::set<int> need_BVH_objects; // objects that need bounding box for a more detailed contact detection
	for (int i = 0; i < BBX_pair.size(); i++)
	{
		need_BVH_objects.insert(BBX_pair[i].first);
		need_BVH_objects.insert(BBX_pair[i].second);
	}

	BVH_objects.clear();
	BVH_objects.insert(BVH_objects.end(), need_BVH_objects.begin(), need_BVH_objects.end());
	if (!advection)
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
			triSimMesh.build_BVH_object_advect(dilation, BVH_objects[k]);
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


void find_contact_pair_element_IPC(
	const std::vector<std::pair<int, int>>& BVH_pair,
	triMesh& triSimMesh,
	FEMParamters& parameters,
	contact_Info& contact_pairs)
{
	// find potential contact face pair
	std::vector<std::vector<Vector5i>> contact_pointTriangle_pair(BVH_pair.size());
	std::vector<std::vector<double>> contact_pointTriangle_dis2(BVH_pair.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int i = 0; i < BVH_pair.size(); i++)
	{
		int obj_1 = BVH_pair[i].first, obj_2 = BVH_pair[i].second;
		std::vector<Vector5i> cont_pair;
		std::vector<double> ct_dis2;


		std::vector<std::pair<int, int>> result_;
		queryBVH(triSimMesh.allObjects[obj_1].object_BVH_nodes, triSimMesh.allObjects[obj_2].object_BVH_faces, result_);
		for (int k = 0; k < result_.size(); k++)
		{
			int vert = result_[k].first, face = result_[k].second;


			Eigen::Vector3d P = triSimMesh.allObjects[obj_1].pos_node_surface[vert];
			Eigen::Vector3i triVerts = triSimMesh.allObjects[obj_2].objectSurfaceMesh.faces[face];
			Eigen::Vector3d A = triSimMesh.allObjects[obj_2].pos_node_surface[triVerts[0]];
			Eigen::Vector3d B = triSimMesh.allObjects[obj_2].pos_node_surface[triVerts[1]];
			Eigen::Vector3d C = triSimMesh.allObjects[obj_2].pos_node_surface[triVerts[2]];


			int type = DIS::dType_PT(P, A, B, C);
			double dis2 = 0;
			DIS::computePointTriD(P, A, B, C, dis2);

			if (dis2 <= squaredDouble(parameters.IPC_dis)) // only calculate the energy when the distance is smaller than the threshold
			{
				Vector5i ct = { obj_1 ,result_[k].first ,obj_2,result_[k].second, type };
				cont_pair.push_back(ct);
				ct_dis2.push_back(dis2);
			}
		}
		contact_pointTriangle_pair[i] = cont_pair;
		contact_pointTriangle_dis2[i] = ct_dis2;
	}


	// find potential contact edge pair
	std::vector<std::vector<Vector5i>> contact_edgeEdge_pair(BVH_pair.size());
	std::vector<std::vector<double>> contact_edgeEdge_dis2(BVH_pair.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int i = 0; i < BVH_pair.size(); i++)
	{
		int obj_1 = BVH_pair[i].first, obj_2 = BVH_pair[i].second;
		std::vector<Vector5i> cont_pair;
		std::vector<double> ct_dis2;


		std::vector<std::pair<int, int>> result_;
		queryBVH(triSimMesh.allObjects[obj_1].object_BVH_edges, triSimMesh.allObjects[obj_2].object_BVH_edges, result_);


		// compute the actual contact	
		for (int k = 0; k < result_.size(); k++)
		{
			int edge1 = result_[k].first, edge2 = result_[k].second;
			Eigen::Vector2i edge1_vert = triSimMesh.allObjects[obj_1].objectSurfaceMesh.index_boundaryEdge[edge1];
			Eigen::Vector2i edge2_vert = triSimMesh.allObjects[obj_2].objectSurfaceMesh.index_boundaryEdge[edge2];

			int P1I = edge1_vert[0], P2I = edge1_vert[1], Q1I = edge2_vert[0], Q2I = edge2_vert[1];

			Eigen::Vector3d P1 = triSimMesh.allObjects[obj_1].pos_node_surface[P1I];
			Eigen::Vector3d P2 = triSimMesh.allObjects[obj_1].pos_node_surface[P2I];
			Eigen::Vector3d Q1 = triSimMesh.allObjects[obj_2].pos_node_surface[Q1I];
			Eigen::Vector3d Q2 = triSimMesh.allObjects[obj_2].pos_node_surface[Q2I];

			int type = DIS::dType_EE(P1, P2, Q1, Q2);
			double dis2 = 0;
			DIS::computeEdgeEdgeD(P1, P2, Q1, Q2, dis2);
			if (dis2 <= squaredDouble(parameters.IPC_dis)) // only calculate the energy when the distance is smaller than the threshold
			{

				Vector5i ct = { obj_1 ,result_[k].first ,obj_2,result_[k].second, type };
				cont_pair.push_back(ct);
				ct_dis2.push_back(dis2);
			}
		}
		contact_edgeEdge_pair[i] = cont_pair;
		contact_edgeEdge_dis2[i] = ct_dis2;
	}


	std::vector<double> PT_DIS;
	// find final actual contact pairs
	for (int i = 0; i < contact_pointTriangle_pair.size(); i++)
	{
		if (contact_pointTriangle_pair[i].size() != 0)
		{
			contact_pairs.Point_Triangle.insert(contact_pairs.Point_Triangle.end(), contact_pointTriangle_pair[i].begin(), contact_pointTriangle_pair[i].end());
			PT_DIS.insert(PT_DIS.end(), contact_pointTriangle_dis2[i].begin(), contact_pointTriangle_dis2[i].end());
		}
	}


	for (int i = 0; i < contact_edgeEdge_pair.size(); i++)
	{
		if (contact_edgeEdge_pair[i].size() != 0)
		{
			contact_pairs.Edge_Edge.insert(contact_pairs.Edge_Edge.end(), contact_edgeEdge_pair[i].begin(), contact_edgeEdge_pair[i].end());
			contact_pairs.Edge_Edge_Dis.insert(contact_pairs.Edge_Edge_Dis.end(), contact_edgeEdge_dis2[i].begin(), contact_edgeEdge_dis2[i].end());
		}
	}




	// find final actual contact pairs
	for (int i = 0; i < contact_pairs.Point_Triangle.size(); i++)
	{
		if (contact_pairs.Point_Triangle[i][4] >= 0)
		{
			if (contact_pairs.Point_Triangle[i][4] <= 2)
			{
				contact_pairs.Point_Triangle_PP_index.push_back(i);
				contact_pairs.Point_Triangle_PP_Dis.push_back(PT_DIS[i]);
			}
			else if (contact_pairs.Point_Triangle[i][4] > 2 && contact_pairs.Point_Triangle[i][4] <= 5)
			{
				contact_pairs.Point_Triangle_PE_index.push_back(i);
				contact_pairs.Point_Triangle_PE_Dis.push_back(PT_DIS[i]);
			}
			else if (contact_pairs.Point_Triangle[i][4] == 6)
			{
				contact_pairs.Point_Triangle_PT_index.push_back(i);
				contact_pairs.Point_Triangle_PT_Dis.push_back(PT_DIS[i]);
			}

		}
	}



	// Step 3: find ground contact pairs if any
	if (parameters.enableGround == true)
	{
		for (int ft = 0; ft < triSimMesh.allObjects.size(); ft++)
		{
			for (int vt = 0; vt < triSimMesh.allObjects[ft].pos_node_surface.size(); vt++)
			{
				if (triSimMesh.allObjects[ft].pos_node_surface[vt][2] <= parameters.IPC_dis)
				{
					Vector2i gt = { ft,vt };
					contact_pairs.Point_Ground.push_back(gt);
				}
			}

		}

	}


}


void find_contact_pair_element_advect(
	const std::vector<std::pair<int, int>>& BVH_pair,
	triMesh& triSimMesh,
	FEMParamters& parameters,
	contact_Info& contact_pairs)
{

	// find potential contact face pair
	std::vector<std::vector<Vector5i>> contact_pointTriangle_pair(BVH_pair.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int i = 0; i < BVH_pair.size(); i++)
	{
		int obj_1 = BVH_pair[i].first, obj_2 = BVH_pair[i].second;
		std::vector<Vector5i> cont_pair;


		std::vector<std::pair<int, int>> result_;
		queryBVH(triSimMesh.allObjects[obj_1].object_BVH_nodes, triSimMesh.allObjects[obj_2].object_BVH_faces, result_);
		// compute the actual contact
		for (int k = 0; k < result_.size(); k++)
		{
			Vector5i ct = { obj_1 ,result_[k].first ,obj_2,result_[k].second, -1};
			cont_pair.push_back(ct);
		}
		contact_pointTriangle_pair[i] = cont_pair;
	}


	// find potential contact edge pair
	std::vector<std::vector<Vector5i>> contact_edgeEdge_pair(BVH_pair.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
	for (int i = 0; i < BVH_pair.size(); i++)
	{
		int obj_1 = BVH_pair[i].first, obj_2 = BVH_pair[i].second;
		std::vector<Vector5i> cont_pair;


		std::vector<std::pair<int, int>> result_;
		queryBVH(triSimMesh.allObjects[obj_1].object_BVH_edges, triSimMesh.allObjects[obj_2].object_BVH_edges, result_);


		// compute the actual contact	
		for (int k = 0; k < result_.size(); k++)
		{
			Vector5i ct = { obj_1 ,result_[k].first ,obj_2, result_[k].second ,-1};
			cont_pair.push_back(ct);

			//std::cout << "ct = (" << result_[k].first << "," << result_[k].second << ")" << std::endl;
		}
		contact_edgeEdge_pair[i] = cont_pair;
	}



	// find final actual contact pairs
	for (int i = 0; i < contact_pointTriangle_pair.size(); i++)
	{
		if (contact_pointTriangle_pair[i].size() != 0)
		{
			contact_pairs.Point_Triangle.insert(contact_pairs.Point_Triangle.end(), contact_pointTriangle_pair[i].begin(), contact_pointTriangle_pair[i].end());
		}
	}


	for (int i = 0; i < contact_edgeEdge_pair.size(); i++)
	{
		if (contact_edgeEdge_pair[i].size() != 0)
		{
			contact_pairs.Edge_Edge.insert(contact_pairs.Edge_Edge.end(), contact_edgeEdge_pair[i].begin(), contact_edgeEdge_pair[i].end());
		}
	}



	// Step 3: find ground contact pairs if any
	if (parameters.enableGround == true)
	{
		for (int ft = 0; ft < triSimMesh.allObjects.size(); ft++)
		{
			for (int vt = 0; vt < triSimMesh.allObjects[ft].pos_node_surface.size(); vt++)
			{
				if (triSimMesh.allObjects[ft].pos_node_surface[vt][2] <= parameters.IPC_dis)
				{
					Vector2i gt = {ft,vt};
					contact_pairs.Point_Ground.push_back(gt);
				}
			}

		}

	}


	std::cout << "			PT = " << contact_pairs.Point_Triangle.size() << "; EE = " << contact_pairs.Edge_Edge.size() << "; PG = " << contact_pairs.Point_Ground.size() << std::endl;



	//std::cout << "			PT.size() = " << contact_pairs.PT.size();
	//std::cout << "; EE.size() = " << contact_pairs.EE.size();
	//std::cout << "; PG.size() = " << contact_pairs.PG.size() << std::endl;


}




void find_contact(
	triMesh& triSimMesh,
	FEMParamters& parameters,
	contact_Info& contact_pairs,
	const std::vector<Vector12d>& moving_direction_ABD)
{
	contact_pairs.clear();
	std::vector<int> BVH_objects;

	if (moving_direction_ABD.size() != 0)
	{
		std::vector<std::pair<int, int>> BBX_contact = find_contact_pair_BBX_level(parameters.IPC_dis / 2.0, triSimMesh, moving_direction_ABD);
		if (BBX_contact.size() != 0 || parameters.enableGround == true)
		{
			std::vector<std::pair<int, int>> BVH_contact = find_contact_pair_BVH_level(parameters.IPC_dis / 2.0, BBX_contact, triSimMesh, BVH_objects, true);
			if (BVH_contact.size() != 0 || parameters.enableGround == true)
			{
				find_contact_pair_element_advect(BVH_contact, triSimMesh, parameters, contact_pairs);
			}
		}
	}
	else  // detect the contact at the begining of a timestep
	{
		std::vector<std::pair<int, int>> BBX_contact = find_contact_pair_BBX_level(parameters.IPC_dis / 2.0, triSimMesh);
		if (BBX_contact.size() != 0 || parameters.enableGround == true)
		{
			std::vector<std::pair<int, int>> BVH_contact = find_contact_pair_BVH_level(parameters.IPC_dis / 2.0, BBX_contact, triSimMesh, BVH_objects, false);
			if (BVH_contact.size() != 0 || parameters.enableGround == true)
			{
				find_contact_pair_element_IPC(BVH_contact, triSimMesh, parameters,contact_pairs);
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





double calMaxStep(
	FEMParamters& parameters,
	triMesh& triSimMesh,
	contact_Info& contact_pairs,
	const int timestep)
{
	double dist_threshold = parameters.IPC_dis;
	double eta = parameters.IPC_eta;
	double step = 1.0;

	std::cout << "			Step size calculation: PT = " << contact_pairs.Point_Triangle.size() << "; EE = " << contact_pairs.Edge_Edge.size() << "; PG = " << contact_pairs.Point_Ground.size() << std::endl;



	// PT_PT maximum step
	if (contact_pairs.Point_Triangle.size() != 0)
	{
		std::vector<double> PT_Step(contact_pairs.Point_Triangle.size(), 1.0);
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < contact_pairs.Point_Triangle.size(); i++)
		{
			int obj_1 = contact_pairs.Point_Triangle[i][0], obj_2 = contact_pairs.Point_Triangle[i][2];
			int vert = contact_pairs.Point_Triangle[i][1], face = contact_pairs.Point_Triangle[i][3];
	
			Eigen::Vector3i triVerts = triSimMesh.allObjects[obj_2].objectSurfaceMesh.faces[face];

			Eigen::Vector3d P = triSimMesh.allObjects[obj_1].pos_node_surface[vert];
			Eigen::Vector3d A = triSimMesh.allObjects[obj_2].pos_node_surface[triVerts[0]];
			Eigen::Vector3d B = triSimMesh.allObjects[obj_2].pos_node_surface[triVerts[1]];
			Eigen::Vector3d C = triSimMesh.allObjects[obj_2].pos_node_surface[triVerts[2]];

			Eigen::Vector3d dP = triSimMesh.allObjects[obj_1].pos_node_surface_direction[vert];
			Eigen::Vector3d dA = triSimMesh.allObjects[obj_2].pos_node_surface_direction[triVerts[0]];
			Eigen::Vector3d dB = triSimMesh.allObjects[obj_2].pos_node_surface_direction[triVerts[1]];
			Eigen::Vector3d dC = triSimMesh.allObjects[obj_2].pos_node_surface_direction[triVerts[2]];


			PT_Step[i] = pointTriangleCCDNarrowphase(P, dP, A,dA, B, dB, C, dC, eta);

		}
		double min_PT_Step = *std::min_element(PT_Step.begin(), PT_Step.end());
		step = std::min(step, min_PT_Step);
	}


	// EE_EE maximum step
	if (contact_pairs.Edge_Edge.size() != 0)
	{
		std::vector<double> EE_Step(contact_pairs.Edge_Edge.size(), 1.0);
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < contact_pairs.Edge_Edge.size(); i++)
		{
			int obj_1 = contact_pairs.Edge_Edge[i][0], obj_2 = contact_pairs.Edge_Edge[i][2];
			int edge1 = contact_pairs.Edge_Edge[i][1], edge2 = contact_pairs.Edge_Edge[i][3];

			Eigen::Vector2i edge1_vert = triSimMesh.allObjects[obj_1].objectSurfaceMesh.index_boundaryEdge[edge1];
			Eigen::Vector2i edge2_vert = triSimMesh.allObjects[obj_2].objectSurfaceMesh.index_boundaryEdge[edge2];


			int P1I = edge1_vert[0], P2I = edge1_vert[1], Q1I = edge2_vert[0], Q2I = edge2_vert[1];

			Eigen::Vector3d P1 = triSimMesh.allObjects[obj_1].pos_node_surface[P1I];
			Eigen::Vector3d P2 = triSimMesh.allObjects[obj_1].pos_node_surface[P2I];
			Eigen::Vector3d Q1 = triSimMesh.allObjects[obj_2].pos_node_surface[Q1I];
			Eigen::Vector3d Q2 = triSimMesh.allObjects[obj_2].pos_node_surface[Q2I];


			Eigen::Vector3d dP1 = triSimMesh.allObjects[obj_1].pos_node_surface_direction[P1I];
			Eigen::Vector3d dP2 = triSimMesh.allObjects[obj_1].pos_node_surface_direction[P2I];
			Eigen::Vector3d dQ1 = triSimMesh.allObjects[obj_2].pos_node_surface_direction[Q1I];
			Eigen::Vector3d dQ2 = triSimMesh.allObjects[obj_2].pos_node_surface_direction[Q2I];


			EE_Step[i] = edgeEdgeCCDNarrowphase(P1, P2, dP1,
				dP2, Q1, Q2, dQ1, dQ2, eta);

		}

		double min_EE_Step = *std::min_element(EE_Step.begin(), EE_Step.end());
		step = std::min(step, min_EE_Step);
	}


	// PG_PG maximum step
	if (contact_pairs.Point_Ground.size() != 0)
	{
		std::vector<double> PG_Step(contact_pairs.Point_Ground.size(), 1.0);
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int i = 0; i < contact_pairs.Point_Ground.size(); i++)
		{
			Vector2i ct = contact_pairs.Point_Ground[i];
			int objInd = ct[0];
			int ptInd = ct[1];

			Eigen::Vector3d P = triSimMesh.allObjects[objInd].pos_node_surface[ptInd];
			double coor_z = P[2];		
			PG_Step[i] = coor_z * (1.0 - parameters.IPC_eta) / std::abs(triSimMesh.allObjects[objInd].pos_node_surface_direction[ptInd][2]);
		
		}

		double min_PG_Step = *std::min_element(PG_Step.begin(), PG_Step.end());
		step = std::min({ step, min_PG_Step });
	}




	return step;

}




// calculate the maximum feasible step size
double calMaxStepSize(
	FEMParamters& parameters,
	objMeshFormat& surfaceInfo,
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


