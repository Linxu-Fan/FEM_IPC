#ifndef AABB_H
#define AABB_H

#include "utils.h"


struct AABB
{
	Eigen::Vector3d min = Eigen::Vector3d::Ones() * 1.0e9;
	Eigen::Vector3d max = -Eigen::Vector3d::Ones() * 1.0e9;

	int vert_edge_face = 0; // indicate the AABB type. 0: node; 1: edge; 2: face
	int index = -99; // index of vertex/edge/face in the local mesh
	Eigen::Vector2i edge_vertices = Eigen::Vector2i::Ones() * (-99);
	Eigen::Vector3i face_vertices = Eigen::Vector3i::Ones() * (-99);

	// initialize face
	void init(const std::vector<Eigen::Vector3d>& pos_node_surface, Eigen::Vector3i& face_vertices_, int index_, double dx);

	// initialize edge
	void init(const std::vector<Eigen::Vector3d>& pos_node_surface, Eigen::Vector2i& edge_vertices_, int index_, double dx);

	// initialize node
	void init(const std::vector<Eigen::Vector3d>& pos_node_surface, int index_, double dx);

	// initialize face
	void init_advect(const std::vector<Eigen::Vector3d>& pos_node_surface, Eigen::Vector3i& face_vertices_, int index_, double dx, const std::vector<Eigen::Vector3d>& direction);

	// initialize edge
	void init_advect(const std::vector<Eigen::Vector3d>& pos_node_surface, Eigen::Vector2i& edge_vertices_, int index_, double dx, const std::vector<Eigen::Vector3d>& direction);

	// initialize node
	void init_advect(const std::vector<Eigen::Vector3d>& pos_node_surface, int index_, double dx, const std::vector<Eigen::Vector3d>& direction);



	void compute_min_max(const std::vector<Eigen::Vector3d>& pos_node_surface, double dx);

	void compute_min_max_advect(const std::vector<Eigen::Vector3d>& pos_node_surface, double dx, const std::vector<Eigen::Vector3d>& direction);

	bool intersects(const AABB& other);

	void dilate(double dx);
};

struct BVHNode {
	AABB box;
	BVHNode* left = nullptr;
	BVHNode* right = nullptr;
	size_t depth = 0;


	bool isLeaf() {
		return left == nullptr && right == nullptr;
	}

	/**
	 * @brief find the bounding box vertices of this BVH node
	 *
	 * @param one_layer true: only export this layer; false export children layers as well
	 */
	void find_bounding_box_vertices(bool one_layer, std::vector<Eigen::Vector3d>& vertices);

	/**
	 * @brief find the bounding box vertices of leaf node
	 *
	 */
	void find_bounding_box_vertices_leaf(std::vector<Eigen::Vector3d>& vertices);

	/**
	 * @brief export the bounding box of this BVH node
	 *
	 * @param one_layer true: only export this layer; false export children layers as well
	 */
	void export_bounding_box_mesh(bool one_layer, std::string fileName, bool leaf_only = false);

};

void deleteBVH(BVHNode* node);

BVHNode* buildBVH(const std::vector<AABB>& triangles, int start, int end, size_t depth);

void updateBVHLeafNodes(BVHNode* node, const std::vector<Eigen::Vector3d>& pos_node_surface, double dilation);

void parallelUpdateBVH(BVHNode* root, const std::vector<Eigen::Vector3d>& pos_node_surface, double dilation);

void queryBVH(BVHNode* nodeA, BVHNode* nodeB, std::vector<std::pair<int, int>>& results);

void count_leaf_node(BVHNode* nodeA, int& num_leaves);



#endif