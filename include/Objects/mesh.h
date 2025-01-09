#ifndef MESH_H
#define MESH_H

#include "utils.h"
#include "objMesh.h"
#include "MLS.h"
#include "AABB.h"

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
	double per_point_volume = 0.01; // the volume of each sampled point
	std::string note = "";
};

class ABD_Info
{
public:
	//std::vector<bool> breakable; // if this object is breakable or not

	// data structure specialized for ABD
	std::vector<Matrix12d> massMatrix_ABD; // the mass matrix of each mesh if in ABD mode
	std::vector<Vector12d> affine;
	std::vector<Eigen::Vector3d> translation_prev_ABD; // translation of each mesh in previous timestep if in ABD mode
	std::vector<Eigen::Vector3d> translation_vel_ABD; // translation velocity
	std::vector<Eigen::Vector3d> translation_ABD; // translation of each mesh if in ABD mode
	std::vector<Eigen::Matrix3d> deformation_prev_ABD; // deformation of each mesh in previous timestep if in ABD mode
	std::vector<Eigen::Matrix3d> deformation_vel_ABD; // deformation velocity
	std::vector<Eigen::Matrix3d> deformation_ABD; // deformation of each mesh if in ABD mode
	std::vector<double> volume_ABD; // volume of each mesh if in ABD mode

};


struct bounding_box
{
	Eigen::Vector3d left_front_bottom = Eigen::Vector3d::Zero();
	Eigen::Vector3d right_front_bottom = Eigen::Vector3d::Zero();
	Eigen::Vector3d right_back_bottom = Eigen::Vector3d::Zero();
	Eigen::Vector3d left_back_bottom = Eigen::Vector3d::Zero();
	Eigen::Vector3d left_front_top = Eigen::Vector3d::Zero();
	Eigen::Vector3d right_front_top = Eigen::Vector3d::Zero();
	Eigen::Vector3d right_back_top = Eigen::Vector3d::Zero();
	Eigen::Vector3d left_back_top = Eigen::Vector3d::Zero();

	Eigen::Vector3d min = Eigen::Vector3d::Ones() * 1.0e9;
	Eigen::Vector3d max = -Eigen::Vector3d::Ones() * 1.0e9;


	void cal_min_max_ABD(const Vector12d& affine, const double& dilation);

	void merges(const bounding_box& other);

	bool intersects(const bounding_box& other);

	void export_BBX_rest(std::string fileName);

	void export_BBX_world(std::string fileName);

};


//////////////////////////////////////////////
// Triangular mesh for simulation
//////////////////////////////////////////////

struct ABD_Object
{
	std::string objectNote = "";              //
	Material objectMaterial;                  //
	objMeshFormat objectSurfaceMesh;          //
	bool breakable = false;                   //
	double volume = 0;                        //
	double per_point_volume = 0.01;


	Vector12d affine_last = Vector12d::Zero();
	Vector12d affine_prev = Vector12d::Zero();
	Vector12d affine = Vector12d::Zero();	
	Matrix12d massMatrix_ABD = Matrix12d::Zero();

	Eigen::Vector3d translation_last_ABD = Eigen::Vector3d::Zero();
	Eigen::Vector3d translation_prev_ABD = Eigen::Vector3d::Zero();
	Eigen::Vector3d translation_vel_ABD = Eigen::Vector3d::Zero();
	Eigen::Vector3d translation_ABD = Eigen::Vector3d::Zero();

	Eigen::Matrix3d deformation_last_ABD = Eigen::Matrix3d::Identity();
	Eigen::Matrix3d deformation_prev_ABD = Eigen::Matrix3d::Identity();
	Eigen::Matrix3d deformation_vel_ABD = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d deformation_ABD = Eigen::Matrix3d::Identity();
	bool need_update_rest_position = false; // update the position in the rest configuration


	// sampled points inside of the object
	std::vector<Eigen::Vector3d> pos_node_interior; // position of each point in the interior                                  
	std::vector<Eigen::Vector3d> pos_node_interior_prev; // position of each point in the interior                                  
	std::vector<Eigen::Vector3d> pos_node_Rest_interior; // rest position of each point in the interior        
	std::vector<double> mass_node_interior; // mass of each node in the interior	    
	std::vector<double> vol_node_interior; // volume of each node in the interior


	std::vector<Eigen::Vector3d> pos_node_surface;                    
	std::vector<Eigen::Vector3d> pos_node_surface_prev;                    
	std::vector<Eigen::Vector3d> pos_node_surface_direction; 
	std::vector<Eigen::Vector3d> contactForce_node_surface;


	// Bounding box in the rest configuration
	bounding_box BBX;

	BVHNode* object_BVH_nodes = nullptr;
	BVHNode* object_BVH_edges = nullptr;
	BVHNode* object_BVH_faces = nullptr;

};

class triMesh 
{
public:
	std::vector<ABD_Object> allObjects; // all ABD objects in the simulation                                                   // **********
	std::map<std::string, int> triMeshIndex; // std::string: mesh name; int: mesh index
	int num_meshes = 0; // number of independant ABD objects                                                                   // **********


	objMeshFormat surfaceMeshGlobal;
	//surface_Info surfaceInfo; // store the surface information of the mesh







	/**
	 * @brief build the BVH data for this ABD object in advection mode, i.e., calculate the step size
	 *
	 * @param direction moving direction of the surface mesh
	 */
	void build_BVH_object_advect(double dilation, const int object_index);


	/**
	 * @brief build the BVH data for this ABD object
	 *
	 */
	void build_BVH_object(double dilation, int object_index);


	/**
	 * @brief update mesh's box corner in the rest configuration
	 *
	 */
	void update_box_corner();


	/**
	 * @brief clear all information except allObjects
	 *
	 */
	void clear();


	/**
	 * @brief update the simulation mesh after generating fragments
	 *
	 */
	void updateGlobalSimulationTriMesh_ABD();

	/**
	 * @brief create the simulation mesh the first time, i.e. from configuration file
	 *
	 */
	void createGlobalSimulationTriMesh_ABD(std::vector<meshConfiguration>& configs);


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
	void sample_points_inside(); 

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


