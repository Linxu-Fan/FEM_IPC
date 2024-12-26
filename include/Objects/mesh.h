#ifndef MESH_H
#define MESH_H

#include "utils.h"
#include "objMesh.h"
#include "MLS.h"

// store key information of simulating mesh like velocity, material etc.
struct meshConfiguration
{
	std::string filePath = "";
	bool breakable = false;
	Material mesh_material;
	Eigen::Vector3d velocity = {0, 0, 0};
	Eigen::Vector3d scale = {1.0, 1.0, 1.0 };
	Eigen::Vector3d translation = {0, 0, 0};
	Eigen::Vector3d rotation_angle = {0, 0, 0};
	Eigen::Vector3d rotation_point = {0, 0, 0};

	std::string note = "";
};

class ABD_Info
{
public:
	std::vector<bool> breakable; // if this object is breakable or not

	// data structure specialized for ABD
	std::vector<Matrix12d> massMatrix_ABD; // the mass matrix of each mesh if in ABD mode
	std::vector<Eigen::Vector3d> translation_prev_ABD; // translation of each mesh in previous timestep if in ABD mode
	std::vector<Eigen::Vector3d> translation_vel_ABD; // translation velocity
	std::vector<Eigen::Vector3d> translation_ABD; // translation of each mesh if in ABD mode
	std::vector<Eigen::Matrix3d> deformation_prev_ABD; // deformation of each mesh in previous timestep if in ABD mode
	std::vector<Eigen::Matrix3d> deformation_vel_ABD; // deformation velocity
	std::vector<Eigen::Matrix3d> deformation_ABD; // deformation of each mesh if in ABD mode
	std::vector<double> volume_ABD; // volume of each mesh if in ABD mode

};

class surface_Info
{
public:

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
	void updateBEInfo(std::vector<Eigen::Vector3d>& vertices, std::vector<Eigen::Vector3i>& faces);

};

//////////////////////////////////////////////
// Tetrahedral mesh for simulation
//////////////////////////////////////////////

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
	// export the mesh's edges
	void exportEdges(std::string fileName);



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
	std::vector<Eigen::Vector3d> contactForce_node; // contact force applied to the surface node
	std::vector<Eigen::Vector3d> pos_node_prev; // tetrahedral node's previous position
	std::vector<std::vector<int>> tetrahedrals_node; // tetrahedrals that contains this node
	std::vector<int> materialInd; // index of the materials(materialMesh) used in this tetrahedral
	std::vector<Material> materialMesh; // all materials used in the simulation

	std::map<int, std::vector<MLSPoints>> MLSPoints_tet_map; // int: index of the tetrahedral; std::vector<MLSPoints>: MLS points in this tetrahedral


	surface_Info surfaceInfo; // store the surface information of the mesh


	// create a global mesh that is suitable for simulation
	void createGlobalSimulationMesh(); 
	// calculate the bounding box of the mesh
	std::pair<Eigen::Vector3d, Eigen::Vector3d> calculateBoundingBox();
	// check the largest edge length
	double calLargestEdgeLength();
	// calculate the boundingbox's diagonal size
	double calBBXDiagSize();


	// Given a crack surface, check intersected tetrahedrons and sample MLS points
	void sample_MLS_points(objMeshFormat& crack, int& num_points, double& radius, int& numOfThreads);
	// sample MLS points inside a tetrahedral
	void sample_MLS_points_inside_tetrahedral(objMeshFormat& crack, int& tetIndex, int& num_points, double& radius);
	// using BSF to find neighbouring nodes of a tetrahedron
	std::vector<int> find_Neighbour_Nodes_tetrahedral(int tetIndex, int maxLayers); // maxLayers: the maximum number of layers to search



};

class Mesh_ABD : public Mesh , public ABD_Info
{
public:
	void createGlobalSimulationMesh_ABD();
};






//////////////////////////////////////////////
// Triangular mesh for simulation
//////////////////////////////////////////////

class triMesh : public ABD_Info
{
public:
	std::vector<std::string> triMeshNote;                                                                                      // **********
	std::map<std::string, int> triMeshIndex; // std::string: mesh name; int: mesh index
	std::vector<Material> materialMesh; // materials(materialMesh) used in each ABD body                                       // **********
	std::vector<objMeshFormat> objectSurfaceMeshes; // store all objects' triangular meshes in the scene                       // **********
	std::vector<Eigen::Vector2i> objectSurfaceMeshes_node_start_end; // start and end node of this surface mesh in pos_node_surface (end node index is one smaller)   // **********
	int num_meshes = 0; // number of independant ABD objects                                                                   // **********


	std::vector<std::string> note_node_surface;
	std::vector<Eigen::Vector3d> contactForce_node_surface; // contact force applied to the surface node                       // **********
	std::vector<Eigen::Vector3d> pos_node_surface; // position of each point on the surface                                    // **********
	std::vector<Eigen::Vector3d> pos_node_Rest_surface; // rest position of each point on the surface                          // **********
	std::vector<Eigen::Vector3d> pos_node_prev_surface; // previous position of each point on the surface                      // **********
	std::vector<Eigen::Vector2i> index_node_surface; // To reuse the code of previous implementation, we keep the same data structure of Class Mesh                    // **********
	std::vector<boundaryCondition> boundaryCondition_node_surface;  // the respective boundary condition of each node on the surface    // **********
	objMeshFormat surfaceMeshGlobal; // use global index of vertex                                                             // **********
	
	surface_Info surfaceInfo; // store the surface information of the mesh


	std::vector<std::string> note_node_interior;
	std::vector<Eigen::Vector3d> pos_node_interior; // position of each point in the interior                                  // **********
	std::vector<Eigen::Vector3d> pos_node_Rest_interior; // rest position of each point in the interiorc                       // **********
	std::vector<Eigen::Vector2i> index_node_interior; // To reuse the code of previous implementation, we keep the same data structure of Class Mesh                   // **********
	std::vector<boundaryCondition> boundaryCondition_node_interior;  // the respective boundary condition of each node on the interior    // **********
	std::vector<double> mass_node_interior; // mass of each node in the interior	                                           // **********
	std::vector<double> vol_node_interior; // volume of each node in the interior	                                           // **********


	/**
	 * @brief create the simulation mesh
	 *
	 */
	void createGlobalSimulationTriMesh_ABD(std::vector<meshConfiguration>& configs, const double& per_point_volume);


	// read meshes from file
	void readMeshes(std::vector<meshConfiguration>& configs);

	// build up the surface information 
	void build_surface_mesh();

	/**
	 * @brief sample points inside of each ABD body
	 *
	 * @param per_point_volume  the volume of each sampled points. 
	 * @return void
	 * @note The number of points sampled in each object equals to vol_obj / per_point_volume
	 */
	void sample_points_inside(double per_point_volume); 

	/**
	 * @brief update the ABD system information
	 *
	 */
	void update_ABD_info();

	/**
	 * @brief export the surface mesh
	 *
	 */
	void exportSurfaceMesh(std::string fileName, int timestep);

	/**
	 * @brief largest edge length of the triangular surface mesh
	 *
	 */
	double calLargestEdgeLength();

	/**
	 * @brief update each ABD object's surface mesh
	 *
	 */
	void updateEachObjectSurfaceMesh();

};




#endif


