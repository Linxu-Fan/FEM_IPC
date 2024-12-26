#ifndef EXTRACTCRACK_H

#define EXTRACTCRACK_H

#include "utils.h"
#include "objMesh.h"
#include "voro++.hh"
#include <algorithm>
#include "omp.h"


namespace extractCrackSurface
{

    // MPM weights 
    // Struct of weights
    struct weightAndDreri {
        Eigen::Vector3i ppIndex = { 0 , 0 , 0 }; // each particle's weight
        Eigen::Vector3d space = { 0 , 0 , 0 }; // each particle's weight derivative
        Eigen::MatrixXd weight;
        Eigen::MatrixXd deltaWeight;

        weightAndDreri(Eigen::Vector3i ippIndex, Eigen::Vector3d ispace, Eigen::MatrixXd iweight, Eigen::MatrixXd ideltaWeight) :
            ppIndex(ippIndex),
            space(ispace),
            weight(iweight),
            deltaWeight(ideltaWeight) {}

    };

    // Struct of particles
    struct Particle 
    {
        Eigen::Vector3d pos = { 0, 0, 0 }, vel = { 0, 0, 0 }; // each particle's x and y position , velocity, and momentum
        double m = 0; // each particle's mass
        Eigen::Vector3i ppIndex = { 0, 0, 0 }; // particle base index
        Eigen::MatrixXd weight;
        Eigen::MatrixXd deltaWeight;
        double Dp = 0; // particle's scalar damage value
        Eigen::Vector3d deltaD = { 0, 0, 0 }; // particle's damage gradient
        int color = 0; // boundary that traction is applied

        Particle(Eigen::Vector3d ipos, Eigen::Vector3d ivel, double im, int ic, double iDp)
            : pos(ipos)
            , vel(ivel)
            , m(im)
            , color(ic)
            , Dp(iDp)
        {
        }
    };

    // crack extraction parameters
    struct CRACKParamters {

        // computational domain
        Eigen::Vector3d length = { 1, 1, 1 }; // computation cube lengths of three dimensions (x , y , z). The origin point is (0 , 0 , 0)
        Eigen::Vector3d minCoordinate = { 0, 0, 0 }; // the minimum coordinate of the computation domain

        // Bcakground Eulerian grid
        double dx = 2E-2;

        // openvdb voxel size
        double vdbVoxelSize = 0.0005;

        // applied force
        std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> appliedForce;

        double damageThreshold = 0.97; // after this threshold,
    };

    // Struct of grid
    struct Grid {
        // property of velocity-field 0
        double m = 0; // each node's mass
        Eigen::Vector3d mom = { 0, 0, 0 }; // each node's momentum
        Eigen::Vector3d velocity = { 0, 0, 0 }; // each node's velocity
        Eigen::Vector3d force = { 0, 0, 0 }; // each node's force

        // general grid node property
        Eigen::Vector3i posIndex = { 0, 0, 0 };
        Eigen::Vector3d deltaDi = { 0, 0, 0 }; // gradient of damage field
        double Di = 0; // value of damage field
        double sw = 0; // sum of particle-grid weight

        // particle index in the support radius of this node. The order of the vector is important
        std::vector<int> supportParticles; // store the position of the particle in vector "particles";
        std::vector<double> supportParticlesWeight; // store the weight of particle to the grid node

        // set of crack surface points withing the grid cell
        std::vector<int> crackPoints;
        int nearestPoint = -1000; // (nearestPoint < 0) means it is far away from the crack surface
        Eigen::Vector3d crackSurfaceNormal = { 0, 0, 0 }; // the vector pointing from the nearest point on the crack surface to the grid node

        // parameters of contact algorithm
        double mass_0 = 0;
        Eigen::Vector3d mom_0 = { 0, 0, 0 };
        Eigen::Vector3d velocity_0 = { 0, 0, 0 };
        Eigen::Vector3d force_0 = { 0, 0, 0 };

        double mass_1 = 0;
        Eigen::Vector3d mom_1 = { 0, 0, 0 };
        Eigen::Vector3d velocity_1 = { 0, 0, 0 };
        Eigen::Vector3d force_1 = { 0, 0, 0 };

        Grid(double im)
            : m(im)
        {
        }
    };

    // Struct of particles
    struct Point {
        int index = 0; // index of each point
        Eigen::Vector3d pos = { 0, 0, 0 }; // each point's position
        int numVertices = 0; // number of vertices
        std::vector<Eigen::Vector3d> verticsCoor; // vertices' coordinate
        int numFaces = 0; // number of faces
        std::vector<std::vector<int>> verticesFace; // vertices of each face
        std::vector<Eigen::Vector3d> surfaceNormal; // vertices' coordinate
        std::vector<int> neighbour; // neighbour points that share common faces with this point
        std::vector<int> neighbourCalculated; // neighbour points that have already find shared face

        std::vector<int> neighbourSameSide; // neighbour points that are on the same side with this point
        std::vector<int> neighbourOtherSide; // neighbour points that are on the other side with this point


        std::vector<int> globalIndexed; // if a vertex finds a global index or not. If yes, the index is a value lager than 0; If not, it is -99;
        std::vector<bool> faceIndexed; // if a face is indexed or not;


        Point(int iindex, Eigen::Vector3d ipos, int inumVertices, std::vector<Eigen::Vector3d> iverticsCoor, int inumFaces, std::vector<std::vector<int>> iverticesFace, std::vector<Eigen::Vector3d> isurfaceNormal, std::vector<int> ineighbour)
            : index(iindex)
            , pos(ipos)
            , numVertices(inumVertices)
            , verticsCoor(iverticsCoor)
            , numFaces(inumFaces)
            , verticesFace(iverticesFace)
            , surfaceNormal(isurfaceNormal)
            , neighbour(ineighbour)
        {
        }
    };

    // calculate MPM weights
    weightAndDreri calWeight(double, Eigen::Vector3d);

    // calculate the MPM grid node hash index
    int calculateID(int x, int y, int z, Eigen::Vector3d len, double dx);

    // calculate the damage gradient of all particles and grid nodes.
    void calDamageGradient(std::vector<Particle>* particles, CRACKParamters param, double dx, std::map<int, int>* gridMap, std::vector<Grid>* nodesVec);
    
    // calculate the damage gradient of any give point
    Eigen::Vector3d calDamageGradientPoint(Eigen::Vector3d pos, CRACKParamters param, double dx, std::map<int, int>* gridMap, std::vector<Grid>* nodesVec);

    // calculate the damage value of any give point
    double calDamageValuePoint(Eigen::Vector3d pos, CRACKParamters param, double dx, std::map<int, int>* gridMap, std::vector<Grid>* nodesVec);

    // Set the damage phase of a grid node vector into a specific value
    void setNodeValue(std::vector<Grid>* nodesVec, int va);

    // Find the bounding box boundary nodes and set its damage phase into a specific value
    void findBoundaryNodes(std::vector<Particle>* particles, std::vector<Grid>* nodesVec, std::map<int, int>* gridMap, struct CRACKParamters parti, int va);

    // Read particles' positions and damage phases
    void readParticles(std::vector<Particle>* particlesRaw, std::vector<Particle>* particles, bool ifFully, struct CRACKParamters param);

    // Calculate the damage value of any point and return the value
    double ifFullyDamaged(Eigen::Vector3d pos, CRACKParamters param, std::map<int, int>* gridMap, std::vector<Grid>* nodesVec);

    // Store the index of each vertex and its position
    int findIndexVertex(Eigen::Vector3d pos, std::vector<Eigen::Vector3d>* vertexIndex);

    // Read all structured nodes and calculate the damage gradient
    void readParticlesAndCalGradient(std::vector<Grid>* fullyDamagedParticlesNodesVec, std::vector<Particle>* particles, CRACKParamters param, std::map<int, int>* gridMap, std::vector<Grid>* nodesVec);

    // Find paths between two nodes
    bool findPath(Eigen::Vector3d startNode, Eigen::Vector3d stopNode, CRACKParamters param, std::vector<Grid>* fullyDamagedParticlesNodesVec, std::map<int, int>* fullyDamagedParticlesGridMap, std::vector<int>* surfaceNodesID);

    // Find if a pair of nodes belong to critical nodes. The function return true if one node is a critical node
    bool ifCriticalNode(Eigen::Vector3d node1, CRACKParamters param, std::vector<int>* criticalNodeIndex);

    // Find the nearest boundary node of a critical node
    Eigen::Vector3i findNearestBoundaryNode(int nodeIDPoint, std::vector<Point>* points, std::vector<int>* boundaryNodesID, std::vector<Eigen::Vector3i>* boundaryNodesPosIndex, std::map<int, int>* pointIndexFind, CRACKParamters param, std::vector<int>* criticalNodeIndex);

    // Judge if a pair of points are on different sides of a crack
    bool ifTwoSides(int startNode, int stopNode, std::vector<Point>* points, std::vector<int>* boundaryNodesID, std::vector<Eigen::Vector3i>* boundaryNodesPosIndex, std::map<int, int>* pointIndexFind, CRACKParamters param, std::vector<int>* criticalNodeIndex);

    // Extract the crack surface
    std::tuple<bool, objMeshFormat, objMeshFormat, std::vector<objMeshFormat>> extractCrackSurf(std::vector<Particle>* particlesRaw, CRACKParamters param);

}

#endif
