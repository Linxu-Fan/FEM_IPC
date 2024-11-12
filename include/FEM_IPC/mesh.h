#ifndef MESH_H
#define MESH_H

#include "utils.h"
#include "objMesh.h"

// store key information of simulating mesh like velocity, material etc.
struct meshConfiguration
{
	std::string filePath = "";
	Material mesh_material;
	Eigen::Vector3d velocity = {0, 0, 0};
	Eigen::Vector3d scale = {1.0, 1.0, 1.0 };
	Eigen::Vector3d translation = {0, 0, 0};
	Eigen::Vector3d rotation_angle = {0, 0, 0};
	Eigen::Vector3d rotation_point = {0, 0, 0};

	std::string note = "";
};


class MLSPoints
{
public:
	double volume = 0; // volume of the MLS point
	Eigen::Vector3d pos = { 0,0,0 }; // position of the MLS point
	Eigen::Matrix3d F = Eigen::Matrix3d::Identity(); // deformation gradient of the MLS point
	std::vector<int> index_node; // index of the nodes in the MLS point's support domain; Index is the index of the node in the original simulation mesh
	std::vector<double> weight; // weight of each node in the MLS point's support domain
	std::vector<Eigen::Matrix<double, 9, 3>> dFdx; // the partial derivative of deformation gradient wrt the node position of the tetrahedral

	// start index of the global gradient vector and hessian matrix
	int startIndex_grad = 0;
	int startIndex_hess = 0;

	

};


// tetrahedral mesh class
class tetMesh
{
public:
	
	std::string tetMeshNote = "";
	std::vector<Eigen::Vector3d> pos_node; // position of each node
	std::vector<Eigen::Vector3d> pos_node_Rest; // tetrahedral node's position at the rest configuration. Note it is different from pos_node_prev which is the position at timestep = n - 1
	std::vector<Eigen::Vector3d> vel_node; // velocity of each node
	std::vector<Eigen::Vector4i> tetrahedrals;
	std::vector<boundaryCondition> boundaryCondition_node; // the respective boundary condition of each node
	std::vector<double> mass_node; // mass of each node	
	std::vector<double> tetra_vol; // volume of the tetrahedral
	std::vector<Eigen::Matrix3d> tetra_DM_inv; // inverse of matrix DM
	std::vector<Eigen::Matrix3d> tetra_F; // deformation gradient of each tetrahedral

	objMeshFormat surfaceMesh;
	Material materialTetMesh;


	////////////////////////////////////////////////////////////////////////
	// Is it possible that node and element are not placed in order? If possible, then the reading code may crash.
	////////////////////////////////////////////////////////////////////////
	// read .msh mesh file
	void readMesh(meshConfiguration& config);
	// initialize the mesh after reading
	void initializeTetMesh(); // initialize the mesh 
	// calculate the DM_inv or DS matrix
	void cal_DM_inv();
	// export surface mesh
	void exportSurfaceMesh(std::string fileName, int timestep = -99);
	// update each tetrahedral's deformation gradient
	void update_F(int numOfThreads);
	// calculate the mass of each node
	void calculateNodeMass();
	// find boundary elements including vertices, edges and triangles
	void findSurfaceMesh();
	// output the mesh
	void output(int timestep);



};


class Mesh : public tetMesh
{
public:
	std::map<std::string, tetMesh> objectsTetMesh; // store all objects' tetrahedral meshes in the scene

	int num_meshes = 1; // number of independant tetrahedral meshes inputed into the simulator
	std::map<std::string, int> tetMeshIndex; // std::string: mesh name; int: mesh index
	std::vector<Eigen::Vector2i> index_node; // index of each node; 1st int: index of the tetrhedral mesh, 2nd int: 0(interior vertex), 1(surface vertex) 

	std::vector<double> Dp; // damage value of the tetrahedral!!!!!!!
	std::vector<std::string> note_node; // note of each node
	std::vector<Eigen::Vector3d> elastForce_node; // elastic force of each node
	std::vector<Eigen::Vector3d> pos_node_prev; // tetrahedral node's previous position
	std::vector<int> materialInd; // index of the materials(materialMesh) used in this tetrahedral
	std::vector<Material> materialMesh; // all materials used in the simulation

	std::map<int, std::vector<MLSPoints>> MLSPoints_tet_map; // int: index of the tetrahedral; std::vector<MLSPoints>: MLS points in this tetrahedral

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


	// create a global mesh that is suitable for simulation
	void createGlobalSimulationMesh(); 
	// calculate the bounding box of the mesh
	std::pair<Eigen::Vector3d, Eigen::Vector3d> calculateBoundingBox();
	// find boundary elements including vertices, edges and triangles
	void findBoundaryElements();
	// update boundary elements' information: area
	void updateBoundaryElementsInfo();
	// check the largest edge length
	double calLargestEdgeLength();
	// calculate the boundingbox's diagonal size
	double calBBXDiagSize();


};


class Mesh_ABD : public Mesh
{
public: 
	// data structure specialized for ABD
	std::vector<Matrix12d> massMatrix_ABD; // the mass matrix of each mesh if in ABD mode
	std::vector<Eigen::Vector3d> translation_prev_ABD; // translation of each mesh in previous timestep if in ABD mode
	std::vector<Eigen::Vector3d> translation_ABD; // translation of each mesh if in ABD mode
	std::vector<Eigen::Vector3d> translation_vel_ABD; // translation velocity
	std::vector<Eigen::Matrix3d> deformation_prev_ABD; // deformation of each mesh in previous timestep if in ABD mode
	std::vector<Eigen::Matrix3d> deformation_ABD; // deformation of each mesh if in ABD mode
	std::vector<Eigen::Matrix3d> deformation_vel_ABD; // deformation velocity
	std::vector<double> volume_ABD; // volume of each mesh if in ABD mode


	void createGlobalSimulationMesh_ABD();
	void updateTranslationDeformation_ABD(std::vector<Eigen::Vector3d>& translation_ABD_Incre, std::vector<Eigen::Matrix3d>& deformation_ABD_Incre);
};





#endif


