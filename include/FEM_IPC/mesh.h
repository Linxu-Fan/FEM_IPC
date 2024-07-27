#ifndef MESH_H
#define MESH_H

#include "utils.h"


struct objMesh 
{
	std::vector<Eigen::Vector3d> vertices;
	std::vector<std::vector<int>> faces;

	void clear();

	void outputFile(std::string fileName, int timestep = -99);

};


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


// boundary condition of each vertex in the mesh
struct boundaryCondition
{
	// vertex's type: 0) default: without constraint; 1) fixed points: velocity = 0; 2) external force, f_ext = xxx
	int type = 0; 
	// 1st and 2nd element are the starting and ending timestep when a boundary condition is applied
	Eigen::Vector2i appliedTime = {0, 100000};
	// for type2 particles, the applied force magnitude
	Eigen::Vector3d force = {0,0,0};
};
 

// store key information of simulating mesh like velocity, material etc.
struct meshConfiguration
{
	std::string filePath = "";
	Material mesh_material;
	Eigen::Vector3d velocity = {0, 0, 0};
	Eigen::Vector3d shift = {0, 0, 0};
};


struct Mesh 
{
    std::vector<Eigen::Vector3d> pos_node; // position of each node
	std::vector<Eigen::Vector3d> vel_node; // velocity of each node
	std::vector<Eigen::Vector3d> elastForce_node; // elastic force of each node
	std::vector<boundaryCondition> boundaryCondition_node; // the respective boundary condition of each node
	std::vector<double> mass_node; // mass of each node
    std::vector<Eigen::Vector4i> tetrahedrals;
    std::vector<double> tetra_vol; // volume of the tetrahedral
	std::vector<Eigen::Matrix3d> tetra_DM_inv; // inverse of matrix DM
	std::vector<Eigen::Matrix3d> tetra_DS; // DS matrix
	std::vector<Eigen::Matrix3d> tetra_F; // deformation gradient of each tetrahedral
    std::vector<std::vector<int>> nodeSharedByElement; // indices of element that shares the node

	std::vector<Eigen::Vector3d> pos_node_prev; // tetrahedral node's previous position

	std::vector<int> materialInd; // index of the materials(materialMesh) used in this tetrahedral
	std::vector<Material> materialMesh; // all materials used in the simulation


	objMesh surfaceMesh;


	//  data structure of boundary elements
	std::map<int, std::set<int>> boundaryVertices; // int: vertex's index in the original mesh; set<int>: neighbour triangles of this vertex	
	std::map<int, double> boundaryVertices_area; // boundary vertex's area (distributed area of this vertex)	
	std::map<int, std::map<int, Eigen::Vector2i>> boundaryEdges; // 1st (smaller one) & 2nd (larger one) int: edge index containing two vertices in the ORIGINAL mesh; Eigen::Vector2i: triangle indices
	std::map<int, std::map<int, double>> boundaryEdges_area; // boundary edge's area (distributed area of this edge)
	std::map<int, std::map<int, int>> boundaryEdge_index; // 1st (smaller one) & 2nd (larger one) int: edge index containing two vertices in the ORIGINAL mesh; 3rd int: index of this edge
	std::map<int, Eigen::Vector2i> index_boundaryEdge; // 1st int: index of this edge; Eigen::Vector2i two vertices in the ORIGINAL mesh
	std::vector<Eigen::Vector3i> boundaryTriangles;
	std::vector<double> boundaryTriangles_area; // boundary triangle's area
	
	


	////////////////////////////////////////////////////////////////////////
	// Is it possible that node and element are not placed in order? If possible, then the reading code may crash.
	////////////////////////////////////////////////////////////////////////
	// read .msh mesh file
	void readMesh(meshConfiguration& config);
	// read multiple meshes at the same time
	void readMeshes(std::vector<meshConfiguration>& config);
	// initialize the mesh after reading
	void initializeMesh(); // initialize the mesh 
	// calculate the DM_inv or DS matrix
	void cal_DS_or_DM(bool DS);
	// output the mesh
	void output(int timestep);
	// export surface mesh
	void exportSurfaceMesh(std::string fileName, int timestep = -99);
	// update each tetrahedral's deformation gradient
	void update_F();
	// calculate the mass of each node
	void calculateNodeMass();
	// calculate the bounding box of the mesh
	std::pair<Eigen::Vector3d, Eigen::Vector3d> calculateBoundingBox();
	// find boundary elements including vertices, edges and triangles
	void findBoundaryElements();
	// update boundary elements' information: area
	void updateBoundaryElementsInfo();
	// check the largest edge length
	double calLargestEdgeLength();
	

};





#endif


