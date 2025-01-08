#include "AABB.h"


void AABB::init(const std::vector<Eigen::Vector3d>& pos_node_surface, Eigen::Vector3i& face_vertices_, int index_, double dx)
{
	vert_edge_face = 2;
	index = index_;
	face_vertices = face_vertices_;
	compute_min_max(pos_node_surface, dx);
	
}

void AABB::init(const std::vector<Eigen::Vector3d>& pos_node_surface, Eigen::Vector2i& edge_vertices_, int index_, double dx)
{
	vert_edge_face = 1;
	index = index_;
	edge_vertices = edge_vertices_;
	compute_min_max(pos_node_surface, dx);

}

void AABB::init(const std::vector<Eigen::Vector3d>& pos_node_surface, int index_, double dx)
{
	vert_edge_face = 0;
	index = index_;
	compute_min_max(pos_node_surface, dx);
}

void AABB::init_advect(const std::vector<Eigen::Vector3d>& pos_node_surface, Eigen::Vector3i& face_vertices_, int index_, double dx, const std::vector<Eigen::Vector3d>& direction)
{
	vert_edge_face = 2;
	index = index_;
	face_vertices = face_vertices_;
	compute_min_max_advect(pos_node_surface, dx, direction);

}

void AABB::init_advect(const std::vector<Eigen::Vector3d>& pos_node_surface, Eigen::Vector2i& edge_vertices_, int index_, double dx, const std::vector<Eigen::Vector3d>& direction)
{
	vert_edge_face = 1;
	index = index_;
	edge_vertices = edge_vertices_;
	compute_min_max_advect(pos_node_surface, dx, direction);

}

void AABB::init_advect(const std::vector<Eigen::Vector3d>& pos_node_surface, int index_, double dx, const std::vector<Eigen::Vector3d>& direction)
{
	vert_edge_face = 0;
	index = index_;
	compute_min_max_advect(pos_node_surface, dx, direction);
}


void AABB::compute_min_max(const std::vector<Eigen::Vector3d>& pos_node_surface, double dx)
{
	if (vert_edge_face == 2)
	{
		Eigen::Vector3d v0 = pos_node_surface[face_vertices[0]];
		Eigen::Vector3d v1 = pos_node_surface[face_vertices[1]];
		Eigen::Vector3d v2 = pos_node_surface[face_vertices[2]];

		min = v0.cwiseMin(v1).cwiseMin(v2);
		max = v0.cwiseMax(v1).cwiseMax(v2);

		dilate(dx);
	}
	else if (vert_edge_face == 1)
	{
		Eigen::Vector3d v0 = pos_node_surface[edge_vertices[0]];
		Eigen::Vector3d v1 = pos_node_surface[edge_vertices[1]];

		min = v0.cwiseMin(v1);
		max = v0.cwiseMax(v1);
		dilate(dx);
	}
	else
	{
		Eigen::Vector3d v0 = pos_node_surface[index];

		min = v0;
		max = v0;
		dilate(dx);
	}

}

void AABB::compute_min_max_advect(const std::vector<Eigen::Vector3d>& pos_node_surface, double dx, const std::vector<Eigen::Vector3d>& direction)
{
	if (vert_edge_face == 2)
	{
		Eigen::Vector3d v0 = pos_node_surface[face_vertices[0]];
		Eigen::Vector3d v1 = pos_node_surface[face_vertices[1]];
		Eigen::Vector3d v2 = pos_node_surface[face_vertices[2]];

		Eigen::Vector3d v3 = v0 + direction[face_vertices[0]];
		Eigen::Vector3d v4 = v1 + direction[face_vertices[1]];
		Eigen::Vector3d v5 = v2 + direction[face_vertices[2]];

		min = v0.cwiseMin(v1).cwiseMin(v2).cwiseMin(v3).cwiseMin(v4).cwiseMin(v5);
		max = v0.cwiseMax(v1).cwiseMax(v2).cwiseMax(v3).cwiseMax(v4).cwiseMax(v5);

		dilate(dx);
	}
	else if (vert_edge_face == 1)
	{
		Eigen::Vector3d v0 = pos_node_surface[edge_vertices[0]];
		Eigen::Vector3d v1 = pos_node_surface[edge_vertices[1]];

		Eigen::Vector3d v2 = v0 + direction[edge_vertices[0]];
		Eigen::Vector3d v3 = v1 + direction[edge_vertices[1]];

		min = v0.cwiseMin(v1).cwiseMin(v2).cwiseMin(v3);
		max = v0.cwiseMax(v1).cwiseMax(v2).cwiseMax(v3);
		dilate(dx);
	}
	else
	{
		Eigen::Vector3d v0 = pos_node_surface[index];

		Eigen::Vector3d v1 = v0 + direction[index];

		min = v0.cwiseMin(v1);
		max = v0.cwiseMax(v1);
		dilate(dx);
	}

}


bool AABB::intersects(const AABB& other)
{
	return (min.x() <= other.max.x() && max.x() >= other.min.x() &&
		min.y() <= other.max.y() && max.y() >= other.min.y() &&
		min.z() <= other.max.z() && max.z() >= other.min.z());
}

void AABB::dilate(double dx) {
	min -= Eigen::Vector3d(dx, dx, dx);
	max += Eigen::Vector3d(dx, dx, dx);
}

void BVHNode::find_bounding_box_vertices_leaf(std::vector<Eigen::Vector3d>& vertices)
{
	if (isLeaf())
	{
		Eigen::Vector3d v0 = { box.min[0], box.min[1], box.min[2] };
		Eigen::Vector3d v1 = { box.max[0], box.min[1], box.min[2] };
		Eigen::Vector3d v2 = { box.max[0], box.max[1], box.min[2] };
		Eigen::Vector3d v3 = { box.min[0], box.max[1], box.min[2] };
		Eigen::Vector3d v4 = { box.min[0], box.min[1], box.max[2] };
		Eigen::Vector3d v5 = { box.max[0], box.min[1], box.max[2] };
		Eigen::Vector3d v6 = { box.max[0], box.max[1], box.max[2] };
		Eigen::Vector3d v7 = { box.min[0], box.max[1], box.max[2] };
		vertices.push_back(v0);
		vertices.push_back(v1);
		vertices.push_back(v2);
		vertices.push_back(v3);
		vertices.push_back(v4);
		vertices.push_back(v5);
		vertices.push_back(v6);
		vertices.push_back(v7);

		return;
	}
	else
	{
		left->find_bounding_box_vertices_leaf(vertices);
		right->find_bounding_box_vertices_leaf(vertices);
	}

}

void BVHNode::find_bounding_box_vertices(bool one_layer, std::vector<Eigen::Vector3d>& vertices)
{
	Eigen::Vector3d v0 = { box.min[0], box.min[1], box.min[2] };
	Eigen::Vector3d v1 = { box.max[0], box.min[1], box.min[2] };
	Eigen::Vector3d v2 = { box.max[0], box.max[1], box.min[2] };
	Eigen::Vector3d v3 = { box.min[0], box.max[1], box.min[2] };
	Eigen::Vector3d v4 = { box.min[0], box.min[1], box.max[2] };
	Eigen::Vector3d v5 = { box.max[0], box.min[1], box.max[2] };
	Eigen::Vector3d v6 = { box.max[0], box.max[1], box.max[2] };
	Eigen::Vector3d v7 = { box.min[0], box.max[1], box.max[2] };
	vertices.push_back(v0);
	vertices.push_back(v1);
	vertices.push_back(v2);
	vertices.push_back(v3);
	vertices.push_back(v4);
	vertices.push_back(v5);
	vertices.push_back(v6);
	vertices.push_back(v7);

	if (one_layer || isLeaf())
	{
		return;
	}
	else
	{
		left->find_bounding_box_vertices(one_layer, vertices);
		right->find_bounding_box_vertices(one_layer, vertices);
	}

}

void BVHNode::export_bounding_box_mesh(bool one_layer, std::string fileName, bool leaf_only)
{

	std::vector<Eigen::Vector3d> vertices;
	if (leaf_only)
	{
		find_bounding_box_vertices_leaf(vertices);		
	}
	else
	{
		find_bounding_box_vertices(one_layer, vertices);
	}
	

	std::ofstream outfile9("./output/" + fileName + ".obj", std::ios::trunc);
	for (int k = 0; k < vertices.size(); k++)
	{
		Eigen::Vector3d scale = vertices[k];
		outfile9 << std::scientific << std::setprecision(8) << "v " << scale[0] << " " << scale[1] << " " << scale[2] << std::endl;
	}
	for (int k = 0; k < vertices.size() / 8; k++)
	{
		outfile9 << "l " << 8 * k + 0 + 1 << " " << 8 * k + 1 + 1 << std::endl;
		outfile9 << "l " << 8 * k + 1 + 1 << " " << 8 * k + 2 + 1 << std::endl;
		outfile9 << "l " << 8 * k + 2 + 1 << " " << 8 * k + 3 + 1 << std::endl;
		outfile9 << "l " << 8 * k + 3 + 1 << " " << 8 * k + 0 + 1 << std::endl;
		outfile9 << "l " << 8 * k + 0 + 1 << " " << 8 * k + 4 + 1 << std::endl;
		outfile9 << "l " << 8 * k + 1 + 1 << " " << 8 * k + 5 + 1 << std::endl;
		outfile9 << "l " << 8 * k + 2 + 1 << " " << 8 * k + 6 + 1 << std::endl;
		outfile9 << "l " << 8 * k + 3 + 1 << " " << 8 * k + 7 + 1 << std::endl;
		outfile9 << "l " << 8 * k + 4 + 1 << " " << 8 * k + 5 + 1 << std::endl;
		outfile9 << "l " << 8 * k + 5 + 1 << " " << 8 * k + 6 + 1 << std::endl;
		outfile9 << "l " << 8 * k + 6 + 1 << " " << 8 * k + 7 + 1 << std::endl;
		outfile9 << "l " << 8 * k + 7 + 1 << " " << 8 * k + 4 + 1 << std::endl;
	}
	outfile9.close();

}

BVHNode* buildBVH(const std::vector<AABB>& triangles, int start, int end, size_t depth) {
	BVHNode* node = new BVHNode();
	node->depth = depth;


	if (end - start == 1) {
		node->box = triangles[start];
		return node;
	}

	//std::cout << "111.0" << std::endl;

	//std::cout << "triangles.size() = " << triangles.size() << "; start = " << start << "; end = " << end << std::endl;


	AABB box;
	for (int i = start; i < end; ++i) {
		//std::cout << "    vef = "<< triangles[i].vert_edge_face<<"; index = " << triangles[i].index << std::endl;
		box.min = box.min.cwiseMin(triangles[i].min);
		box.max = box.max.cwiseMax(triangles[i].max);
	}
	node->box = box;

	//std::cout << "111.1" << std::endl;


	int axis = 0;
	Eigen::Vector3d size = box.max - box.min;
	if (size.y() > size.x() && size.y() > size.z()) axis = 1; 
	if (size.z() > size.x() && size.z() > size.y()) axis = 2; 

	//std::cout << "axis = " << axis << std::endl;

	std::vector<AABB> sortedTriangles(triangles.begin() + start, triangles.begin() + end);
	std::sort(sortedTriangles.begin(), sortedTriangles.end(), [axis](const AABB& a, const AABB& b) {
		return a.min[axis] < b.min[axis];
		});





	int mid = sortedTriangles.size() / 2;

	/*std::cout << "111.2" << std::endl;
	std::cout << "sortedTriangles.size() = " << sortedTriangles.size() <<"; mid = "<<mid<<"; end = "<< sortedTriangles.size() <<std::endl << std::endl;*/
	node->left = buildBVH(sortedTriangles, 0, mid, depth + 1);
	node->right = buildBVH(sortedTriangles, mid, sortedTriangles.size(), depth + 1);




	return node;
}

void deleteBVH(BVHNode* node) 
{
	if (node == nullptr) {
		return; 
	}

	deleteBVH(node->left);
	deleteBVH(node->right);

	delete node;
}

void updateBVHLeafNodes(BVHNode* node, const std::vector<Eigen::Vector3d>& pos_node_surface, double dilation) {
	if (!node) return;


	if (node->isLeaf())
	{
		node->box.compute_min_max(pos_node_surface, dilation);
		return;
	}


	updateBVHLeafNodes(node->left, pos_node_surface, dilation);
	updateBVHLeafNodes(node->right, pos_node_surface, dilation);


	node->box.min = node->left->box.min.cwiseMin(node->right->box.min);
	node->box.max = node->left->box.max.cwiseMax(node->right->box.max);
}

void parallelUpdateBVH(BVHNode* root, const std::vector<Eigen::Vector3d>& pos_node_surface, double dilation)
{
	if (!root) return;

	std::vector<std::thread> threads;
	std::mutex mutex;


	auto task = [&](BVHNode* node) {
		updateBVHLeafNodes(node, pos_node_surface, dilation);
		};


	threads.emplace_back(task, root->left);
	threads.emplace_back(task, root->right);


	for (auto& thread : threads) {
		thread.join();
	}
}

void queryBVH(BVHNode* nodeA, BVHNode* nodeB, std::vector<std::pair<int, int>>& results) {
	if (!nodeA || !nodeB) return;


	if (!nodeA->box.intersects(nodeB->box))
	{
		//std::cout << "Level = " << nodeA->depth << "; " << nodeB->depth << std::endl;
		return;
	}
	else
	{
		if (nodeA->isLeaf())
		{
			if (nodeB->isLeaf())
			{
				//std::cout << "intersect = " << "(" << nodeA->triangleIndex << "," << nodeB->triangleIndex << ")";
				//std::cout << ";add = " << "(" << nodeA << "," << nodeB << ")" << std::endl;
				results.emplace_back(nodeA->box.index, nodeB->box.index);
				return;
			}
			else
			{
				queryBVH(nodeA, nodeB->right, results);
				queryBVH(nodeA, nodeB->left, results);
			}
		}
		else
		{
			if (nodeB->isLeaf())
			{
				queryBVH(nodeA->right, nodeB, results);
				queryBVH(nodeA->left, nodeB, results);
			}
			else
			{
				queryBVH(nodeA->left, nodeB->left, results);
				queryBVH(nodeA->left, nodeB->right, results);
				queryBVH(nodeA->right, nodeB->left, results);
				queryBVH(nodeA->right, nodeB->right, results);
			}
		}

	}

}

void count_leaf_node(BVHNode* nodeA, int& num_leaves)
{
	if (nodeA->isLeaf())
	{
		num_leaves += 1;

		if (nodeA->box.index == 0)
		{
			std::cout << "min = " << nodeA->box.min.transpose() << "; max = " << nodeA->box.max.transpose() << std::endl;
		}
	}
	else
	{
		count_leaf_node(nodeA->left, num_leaves);
		count_leaf_node(nodeA->right, num_leaves);
	}
}

