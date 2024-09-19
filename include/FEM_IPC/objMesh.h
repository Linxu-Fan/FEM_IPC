﻿#ifndef OBJMESH_H
#define OBJMESH_H

#include "utils.h"


// 用于表示三角形面，包含三个节点的ID
struct TriangleFace
{
	int n1, n2, n3;

	TriangleFace(int v1, int v2, int v3) {
		// 确保节点ID始终以相同的顺序存储
		std::vector<int> nodes = { v1, v2, v3 };
		std::sort(nodes.begin(), nodes.end());
		n1 = nodes[0];
		n2 = nodes[1];
		n3 = nodes[2];
	}

	// 为 map 容器提供 < 操作符
	bool operator<(const TriangleFace& other) const {
		return std::tie(n1, n2, n3) < std::tie(other.n1, other.n2, other.n3);
	}
};

void depthFirstSearch(int v, const std::vector<std::set<int>>& adjList, std::vector<bool>& visited, std::vector<int>& component);

struct objMeshFormat
{
    std::vector<Eigen::Vector3d> vertices;
	std::vector<Eigen::Vector2i> edges; // 1st int: smaller vertex index; 2nd int: larger vertex index
    std::vector<Eigen::Vector3i> faces;
	std::vector<std::vector<int>> vertFaces; // faces that share this vertex
	std::vector<objMeshFormat> componentsSep; // separated connected-components
	


	void clear();

	void readObjFile(std::string fileName);

    void sepConnectedComponents(); // separate connected comonents

    void outputFile(std::string fileName, int timestep = -99);

	void findVertFaces_Edges();

};


#endif


