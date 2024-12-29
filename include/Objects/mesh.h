﻿#ifndef MESH_H
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
	double per_point_volume = 0.01; // the volume of each sampled point
	std::string note = "";
};

class ABD_Info
{
public:
	//std::vector<bool> breakable; // if this object is breakable or not

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
	Eigen::Vector2i objectSurfaceMeshes_node_start_end;  //
	Eigen::Vector3d translation_prev_ABD = Eigen::Vector3d::Zero();
	Eigen::Vector3d translation_vel_ABD = Eigen::Vector3d::Zero();
	Eigen::Vector3d translation_ABD = Eigen::Vector3d::Zero();
	Eigen::Matrix3d deformation_prev_ABD = Eigen::Matrix3d::Identity();
	Eigen::Matrix3d deformation_vel_ABD = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d deformation_ABD = Eigen::Matrix3d::Identity();
	bool need_update_rest_position = false; // update the position in the rest configuration

};

class triMesh : public ABD_Info
{
public:
	std::vector<ABD_Object> allObjects; // all ABD objects in the simulation                                                   // **********
	std::map<std::string, int> triMeshIndex; // std::string: mesh name; int: mesh index
	int num_meshes = 0; // number of independant ABD objects                                                                   // **********


	std::vector<std::string> note_node_surface;
	std::vector<Eigen::Vector3d> contactForce_node_surface; // contact force applied to the surface node                       // **********
	std::vector<Eigen::Vector3d> pos_node_surface; // position of each point on the surface                                    // **********
	std::vector<Eigen::Vector3d> pos_node_Rest_surface; // rest position of each point on the surface                          // **********
	std::vector<Eigen::Vector2i> index_node_surface; // To reuse the code of previous implementation, we keep the same data structure of Class Mesh                    // **********
	std::vector<boundaryCondition> boundaryCondition_node_surface;  // the respective boundary condition of each node on the surface    // **********
	

	objMeshFormat surfaceMeshGlobal;
	surface_Info surfaceInfo; // store the surface information of the mesh


	std::vector<std::string> note_node_interior;
	std::vector<Eigen::Vector3d> pos_node_interior; // position of each point in the interior                                  // **********
	std::vector<Eigen::Vector3d> pos_node_Rest_interior; // rest position of each point in the interiorc                       // **********
	std::vector<Eigen::Vector2i> index_node_interior; // To reuse the code of previous implementation, we keep the same data structure of Class Mesh                   // **********
	std::vector<boundaryCondition> boundaryCondition_node_interior;  // the respective boundary condition of each node on the interior    // **********
	std::vector<double> mass_node_interior; // mass of each node in the interior	                                           // **********
	std::vector<double> vol_node_interior; // volume of each node in the interior	                                           // **********


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


