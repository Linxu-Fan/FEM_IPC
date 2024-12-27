#ifndef MPMSIMULATOR_H

#define MPMSIMULATOR_H

#include "extractCrack.h"



namespace mpmSimulator
{

	// Struct of grid
	struct Grid {

		double mass = 0; // each node's mass
		Eigen::Vector3d momentum = { 0 , 0 , 0 }; // each node's momentum
		Eigen::Vector3d velocity = { 0 , 0 , 0 }; // each node's velocity
		Eigen::Vector3d force = { 0 , 0 , 0 }; // each node's force


		Eigen::Vector3i posIndex = { 0 , 0 , 0 }; // node's position in the grid
		Eigen::Vector3d position = { 0 , 0 , 0 }; // node's real coordinate
		Eigen::Vector3d positionUpdated = { 0 , 0 , 0 }; // node's updated coordinate


		// particle index in the support radius of this node. The order of the vector is important
		std::vector<int> supportParticles; // store the position of the particle in vector "particles"; int is the position of the individual particle vector 
		std::vector<double> supportParticlesWeight; // store the weight of particle to the grid node 
		std::vector<Eigen::Vector3d> supportParticlesDeltaWeight; // store the weight of node to the particle


		// if this node contacins two kinds of particles
		Eigen::Vector2i twoPoints = { -99 , -99 };

		double Di = 0; // value of damage field
		double sw = 0; // sum of particle-grid weight
		Eigen::Vector3d deltaDi = { 0, 0, 0 };

	};


	// Struct of particle type-1
	struct MPMParticle {

		double mass = 0; // each particle's mass
		double volume = 0; // each particle's volume

		Eigen::Vector3d velocity = { 0 , 0 , 0 }; // each particle's velocity
		Eigen::Vector3i posIndex = { 0 , 0 , 0 }; // particle base index
		Eigen::Vector3d position = { 0 , 0 , 0 }; // each particle's position 

		Eigen::Matrix3d weight = Eigen::Matrix3d::Zero();
		Eigen::Matrix3d deltaWeight = Eigen::Matrix3d::Zero();
		Eigen::Matrix3d F = Eigen::Matrix3d::Identity(); // each particle's deformation gradient
		Eigen::Matrix3d affine = Eigen::Matrix3d::Zero(); // each particle's affine term
		Eigen::Matrix3d cauchyStress = Eigen::Matrix3d::Zero(); // each particle's internal cauchy stress


		// node index in the support radius of this particle. The order of the vector is important
		std::vector<int> supportNodes; // store the position of the node in vector "nodeVec" 
		std::vector<double> supportNodeWeight; // store the weight of node to the particle
		std::vector<Eigen::Vector3d> supportNodeDeltaWeight; // store the weight of node to the particle


		// material 
		int color = 0;
		double dp = 0;

		Eigen::Vector3d deltaD = { 0 , 0 , 0 }; // particle's damage gradient
		double Dg = 0;

	};


	struct MPMParamters {

		// number of parallel threads
		int numOfThreads = 6;

		Eigen::Vector3d gravity = { 0 , 0 , 0 };
		double dt = 1.0E-6; // timestep

		// Bcakground Eulerian grid
		double dx = 2.0E-3; // grid space

		// damping coefficient
		double nu = 0;
		double damageThreshold = 0.97;

		Material mat_mpm;

	};


	std::string calculateID_string(int x, int y, int z); // string id of the grid node


	// find the surrounding support nodes of each particle and calculate the weights of the particle
	void calWeightsAndNodes(std::vector<MPMParticle>& particles, MPMParamters& param, std::vector<Grid>& nodesVec, std::map<std::string, int>& gridMap);


	// particle to grid transfer
	void particle2Grid(std::vector<MPMParticle>& particles, MPMParamters& param, std::vector<Grid>& nodesVec);


	// update each particle's cauchy stress
	void updateParInternalForce(std::vector<MPMParticle>& particles, MPMParamters& param);


	// calculate the grid node's internal force induced by particles
	void calculateNodeForce(std::vector<MPMParticle>& particles, MPMParamters& param, std::vector<Grid>& nodesVec);


	// particle to grid transfer
	void gridUpdate(std::vector<Grid>& nodesVec, MPMParamters& param);


	// grid to particle transfer
	void grid2Particle(std::vector<MPMParticle>& particles, MPMParamters& param, std::vector<Grid>& nodesVec);


	
	/**
	 * @brief apply contact force from IPC
	 *
	 * @param contactForce  the first .
	 */
	void applyPointForce(MPMParamters& param, std::vector<Grid>& nodesVec, std::map<std::string, int>& gridMap, std::vector<Vector6d>& contactForce);


	// the calculation of each timestep
	void advanceStep(std::vector<MPMParticle>& particles, MPMParamters& param, std::vector<Vector6d>& contactForce, int timestep); // simulation parameters, particle vector, the total number of type-1,2,3 particles


	// extract crack surface
	std::tuple<bool, objMeshFormat, objMeshFormat, std::vector<objMeshFormat>> tryToExtractCracks(std::vector<MPMParticle>& particles, MPMParamters& param, int timestep);


	std::tuple<bool, objMeshFormat, objMeshFormat, std::vector<objMeshFormat>> crackSimulation(
		const std::vector<Eigen::Vector3d>& points, 
		const double& volume, 
		const Material& mat_mpm,
		MPMParamters& param, 
		std::vector<Vector6d>& contactForce, 
		int num_timestep);


	/**
	 * @brief initialize mpm particles from std::vector<Eigen::Vector3d>
	 *
	 * @volume volume of each point
	 * @density density of each point
	 */
	void initialize_mpm_particles(std::vector<MPMParticle>& particles, const std::vector<Eigen::Vector3d>& points, const double& volume, const Material& mat_mpm);
    
}

#endif
