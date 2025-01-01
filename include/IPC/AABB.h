//#ifndef AABB_H
//#define AABB_H
//
//#include "utils.h"
//
//
//
//struct AABB
//{
//	Eigen::Vector3d min = Eigen::Vector3d::Ones() * 1.0e9;
//	Eigen::Vector3d max = -Eigen::Vector3d::Ones() * 1.0e9;
//	int face_index = -99;
//	Eigen::Vector3i face_vertices = Eigen::Vector3i::Ones() * (-99);
//
//
//	void init(const std::vector<Eigen::Vector3d>& pos_node_surface, Eigen::Vector3i& face_vertices_, int index_, double dx);
//
//	void compute_min_max(const std::vector<Eigen::Vector3d>& pos_node_surface, double dx);
//
//	bool intersects(const AABB& other);
//
//	void dilate(double dx);
//};
//
//
//struct BVHNode {
//	AABB box;
//	BVHNode* left = nullptr;
//	BVHNode* right = nullptr;
//	size_t depth = 0;
//
//
//	bool isLeaf() {
//		return left == nullptr && right == nullptr;
//	}
//};
//
//
//
//BVHNode* buildBVH(const std::vector<AABB>& triangles, int start, int end, size_t depth) {
//	BVHNode* node = new BVHNode();
//	node->depth = depth;
//
//
//	if (end - start == 1) {
//		node->box = triangles[start];
//		return node;
//	}
//
//	//std::cout << "111.0" << std::endl;
//
//	//std::cout << "triangles.size() = " << triangles.size() << "; start = " << start << "; end = " << end << std::endl;
//
//
//	AABB box;
//	for (int i = start; i < end; ++i) {
//		//std::cout << "    i = " << i << std::endl;
//		box.min = box.min.cwiseMin(triangles[i].min);
//		box.max = box.max.cwiseMax(triangles[i].max);
//	}
//	node->box = box;
//
//	//std::cout << "111.1" << std::endl;
//
//
//	int axis = 0;
//	Eigen::Vector3d size = box.max - box.min;
//	if (size.y() > size.x() && size.y() > size.z()) axis = 1; // Y ��
//	if (size.z() > size.x() && size.z() > size.y()) axis = 2; // Z ��
//
//	//std::cout << "axis = " << axis << std::endl;
//
//	std::vector<AABB> sortedTriangles(triangles.begin() + start, triangles.begin() + end);
//	std::sort(sortedTriangles.begin(), sortedTriangles.end(), [axis](const AABB& a, const AABB& b) {
//		return a.min[axis] < b.min[axis];
//		});
//
//
//
//
//
//	int mid = sortedTriangles.size() / 2;
//
//	/*std::cout << "111.2" << std::endl;
//	std::cout << "sortedTriangles.size() = " << sortedTriangles.size() <<"; mid = "<<mid<<"; end = "<< sortedTriangles.size() <<std::endl << std::endl;*/
//	node->left = buildBVH(sortedTriangles, 0, mid, depth + 1);
//	node->right = buildBVH(sortedTriangles, mid, sortedTriangles.size(), depth + 1);
//
//
//
//
//	return node;
//}
//
//
//void updateBVHLeafNodes(BVHNode* node, const std::vector<Eigen::Vector3d>& pos_node_surface, double dilation) {
//	if (!node) return;
//
//
//	if (node->isLeaf())
//	{
//		node->box.compute_min_max(pos_node_surface, dilation);
//		return;
//	}
//
//
//	updateBVHLeafNodes(node->left, pos_node_surface, dilation);
//	updateBVHLeafNodes(node->right, pos_node_surface, dilation);
//
//
//	node->box.min = node->left->box.min.cwiseMin(node->right->box.min);
//	node->box.max = node->left->box.max.cwiseMax(node->right->box.max);
//}
//
//
//void parallelUpdateBVH(BVHNode* root, const std::vector<Eigen::Vector3d>& pos_node_surface, double dilation) {
//	if (!root) return;
//
//	std::vector<std::thread> threads;
//	std::mutex mutex;
//
//
//	auto task = [&](BVHNode* node) {
//		updateBVHLeafNodes(node, pos_node_surface, dilation);
//		};
//
//
//	threads.emplace_back(task, root->left);
//	threads.emplace_back(task, root->right);
//
//
//	for (auto& thread : threads) {
//		thread.join();
//	}
//}
//
//
//void queryBVH(BVHNode* nodeA, BVHNode* nodeB, std::vector<std::pair<int, int>>& results) {
//	if (!nodeA || !nodeB) return;
//
//
//	if (!nodeA->box.intersects(nodeB->box))
//	{
//		//std::cout << "Level = " << nodeA->depth << "; " << nodeB->depth << std::endl;
//		return;
//	}
//	else
//	{
//		if (nodeA->isLeaf())
//		{
//			if (nodeB->isLeaf())
//			{
//				//std::cout << "intersect = " << "(" << nodeA->triangleIndex << "," << nodeB->triangleIndex << ")";
//				//std::cout << ";add = " << "(" << nodeA << "," << nodeB << ")" << std::endl;
//				results.emplace_back(nodeA->box.face_index, nodeB->box.face_index);
//				return;
//			}
//			else
//			{
//				queryBVH(nodeA, nodeB->right, results);
//				queryBVH(nodeA, nodeB->left, results);
//			}
//		}
//		else
//		{
//			if (nodeB->isLeaf())
//			{
//				queryBVH(nodeA->right, nodeB, results);
//				queryBVH(nodeA->left, nodeB, results);
//			}
//			else
//			{
//				queryBVH(nodeA->left, nodeB->left, results);
//				queryBVH(nodeA->left, nodeB->right, results);
//				queryBVH(nodeA->right, nodeB->left, results);
//				queryBVH(nodeA->right, nodeB->right, results);
//			}
//		}
//
//	}
//
//}
//
//void count_leaf_node(BVHNode* nodeA, int& num_leaves)
//{
//	if (nodeA->isLeaf())
//	{
//		num_leaves += 1;
//
//		if (nodeA->box.face_index == 0)
//		{
//			std::cout << "min = " << nodeA->box.min.transpose() << "; max = " << nodeA->box.max.transpose() << std::endl;
//		}
//	}
//	else
//	{
//		count_leaf_node(nodeA->left, num_leaves);
//		count_leaf_node(nodeA->right, num_leaves);
//	}
//}
//
//
//
//
//
//
//#endif