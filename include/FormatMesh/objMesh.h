#ifndef OBJMESH_H
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

void depthFirstSearch(int v, const std::vector<std::set<int>>& adjList, 
	std::vector<bool>& visited, std::vector<int>& component);



// Function to compute face normal
Eigen::Vector3d computeNormal(const Eigen::Vector3d& v0,
	const Eigen::Vector3d& v1,
	const Eigen::Vector3d& v2);

// Function to compute projection of points onto an axis
void projectOntoAxis(const std::vector<Eigen::Vector3d>& points,
	const Eigen::Vector3d& axis,
	double& minProj,
	double& maxProj);

// Triangle-Tetrahedron intersection test using SAT
bool triangleTetrahedronIntersect(const Eigen::Vector3d& tri_v0,
	const Eigen::Vector3d& tri_v1,
	const Eigen::Vector3d& tri_v2,
	const Eigen::Vector3d& tet_v0,
	const Eigen::Vector3d& tet_v1,
	const Eigen::Vector3d& tet_v2,
	const Eigen::Vector3d& tet_v3);



// Möller–Trumbore intersection algorithm for line segment and triangle
bool lineSegmentIntersectsTriangle(const Eigen::Vector3d& orig,
	const Eigen::Vector3d& dest,
	const Eigen::Vector3d& vert0,
	const Eigen::Vector3d& vert1,
	const Eigen::Vector3d& vert2,
	Eigen::Vector3d& intersectPoint);



struct objMeshFormat 
{
    std::vector<Eigen::Vector3d> vertices;
	std::vector<Eigen::Vector2i> edges; // 1st int: smaller vertex index; 2nd int: larger vertex index
    std::vector<Eigen::Vector3i> faces;
	std::vector<std::vector<int>> facesPolygonal; // polygonal face
	std::vector<std::vector<int>> vertFaces; // faces that share this vertex
	std::vector<objMeshFormat> componentsSep; // separated connected-components
	

	// initial velocity of this mesh. It is useful when the mesh is used for simulation
	Eigen::Vector3d initialVelocity = Eigen::Vector3d::Zero();
	double volume = 0; // watertight volume of this object
	


	//  data structure of boundary elements
	std::map<int, std::set<int>> boundaryVertices; // int: vertex's index in the original mesh; set<int>: neighbour triangles of this vertex	
	std::map<int, std::set<int>> boundaryVertices_egde; // int: vertex's index in the original mesh; set<int>: neighbour edges of this vertex	
	std::vector<int> boundaryVertices_vec; // int: vertex's index in the original mesh
	std::map<int, double> boundaryVertices_area; // boundary vertex's area (distributed area of this vertex)	
	std::map<int, std::map<int, Eigen::Vector2i>> boundaryEdges; // 1st (smaller one) & 2nd (larger one) int: edge index containing two vertices in the ORIGINAL mesh; Eigen::Vector2i: triangle indices
	std::map<int, std::map<int, double>> boundaryEdges_area; // boundary edge's area (distributed area of this edge)
	std::map<int, std::map<int, int>> boundaryEdge_index; // 1st (smaller one) & 2nd (larger one) int: edge index containing two vertices in the ORIGINAL mesh; 3rd int: index of this edge
	std::map<int, Eigen::Vector2i> index_boundaryEdge; // 1st int: index of this edge; Eigen::Vector2i two vertices in the ORIGINAL mesh
	std::vector<int> index_boundaryEdge_vec; // 1st int: index of this edge
	std::vector<Eigen::Vector3i> boundaryTriangles;
	std::vector<double> boundaryTriangles_area; // boundary triangle's area

	// update the boundary information given a surface mesh
	void updateBEInfo();



	/**
	 * @brief update the mesh as far as possible
	 *
	 * @note Only work for triangular mesh
	 */
	void updateMesh();

	void clear();

	void readObjFile(std::string fileName, bool polygonal = false, Eigen::Affine3d rotation = Eigen::Affine3d::Identity(),
		Eigen::Vector3d scale = Eigen::Vector3d::Ones(), Eigen::Vector3d translation = Eigen::Vector3d::Zero()); // if polygonal is true, read polygonal mesh

    void sepConnectedComponents(); // separate connected comonents. ONLY FOR TRIANGULAR MESH

    void outputFile(std::string fileName, int timestep = -99, bool polygonal = false);

	bool checkIfMeshIntersectWithTetrahedron(const Eigen::Vector3d tet_v0, const Eigen::Vector3d tet_v1, 
		const Eigen::Vector3d tet_v2, const Eigen::Vector3d tet_v3); // check if this mesh intersect with a tetrahedron (v0, v1, v2, v3)

	bool checkIfMeshIntersectWithLine(const Eigen::Vector3d line_pt1, const Eigen::Vector3d line_pt2); // check if this mesh intersect with a line segment (line_pt1, line_pt2)

	// convert to a pure triangle mesh
	void triangulate();

	// convert to openvdb format
	void to_openVDB_format(std::vector<openvdb::Vec3s>& verticesVdb, std::vector<openvdb::Vec3I>& trianglesVdb);


	/**
	 * @brief update the volume of the object
	 *
	 */
	void updateVolume();

	/**
	 * @brief convert the vertices and faces into libigl format for further processing
	 *
	 */
	std::pair<Eigen::MatrixXd, Eigen::MatrixXi> to_libigl_mesh();

	/**
	 * @brief convert the vertices and faces into CGAL format for further processing
	 *
	 */
	CGAL_Surface_mesh to_CGAL_mesh();

	/**
	 * @brief Generate random points inside a closed triangular mesh.
	 *
	 * @param num_samples Number of points to sample
	 */
	std::vector<Eigen::Vector3d> sample_points_inside_mesh(int num_samples);

	/**
	 * @brief Reconstruct the object with openvdb
	 *
	 * @param voxel_size openVDB voxel size
	 */
	objMeshFormat reconstruct_with_vdb(float& voxel_size);

	/**
	 * @brief Boolean operation of difference with another mesh B
	 *
	 */
	objMeshFormat boolean_difference_with_mesh(objMeshFormat& B);

};






#endif


