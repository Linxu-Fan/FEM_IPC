#ifndef UTILS_H
#define UTILS_H

#include <cassert> 
#include "omp.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <cmath>
#include <set>
#include <queue>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <map>
#include <unordered_map>
#include <string.h>
#include <cassert>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry> // For affine transformations
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/AutoDiff>
#include <unsupported/Eigen/KroneckerProduct>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <sstream>
#include <cstdio>
#include <cassert>
#include <assert.h>
#include <cstring>
#include <cfloat>
#include <complex>
#include <cstdlib>
#include <optional> 



#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_stop_predicate.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_distance_placement.h>
#include <CGAL/property_map.h>



typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 CGAL_Point_3;
typedef CGAL::Surface_mesh<CGAL_Point_3> CGAL_Surface_mesh;



namespace PMP = CGAL::Polygon_mesh_processing;
namespace SMS = CGAL::Surface_mesh_simplification;



#include <openvdb/openvdb.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/LevelSetUtil.h>
#include <openvdb/tools/MeshToVolume.h>


#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/decimate.h>
#include <igl/winding_number.h>
#include <igl/centroid.h>


#define PI_Value 3.14159265358979323846 

const double PI = 3.141592653589793238463L;
const double M_2PI = 2 * PI;
const double eps = 1e-12;

typedef std::complex<double> DComplex;


typedef Eigen::Matrix<int, 7, 1> Vector7i;
typedef Eigen::Matrix<int, 5, 1> Vector5i;

typedef Eigen::Matrix<double, 3, 1> Vector3d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 9, 1> Vector9d;
typedef Eigen::Matrix<double, 12, 1> Vector12d;
typedef Eigen::Matrix<double, 3, 3> Matrix3d;
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<double, 9, 9> Matrix9d;
typedef Eigen::Matrix<double, 12, 12> Matrix12d;





struct FEMParamters
{
    double dt = 1.0E-2; // simulation timestep size
    int num_timesteps = 2000; // the number of simulation timesteps
    Eigen::Vector3d gravity = { 0,0,0 }; // the gravity magnitude and direction
    int outputFrequency = 200; // export the simulation result per interval
    std::string model = "neoHookean";
	int numOfThreads = 5; // number of openMP threads

    bool enableGround = true;
	double searchResidual = 1.0e-4;

    bool rigidMode = false; // if rigidMode, IPC will only check contact pairs that belongs to two different objects. It is suitable for small deformation
    std::string simulation_Mode = "Normal"; // Normal: normal simulation; ABD: ABD simulation; Coupling: ABD and deformable
    //std::vector<std::string> objectNames; // name of objects in the scene
	double IPC_dis = 0.01; // the distance of IPC gap, i.e., barrier energy emerges if the distance is smaller than this threshold
	double IPC_eta = 0.1; // the distance ratio that bring two edges or point-triangle to
	double IPC_hashSize = 0.1; // the spatial hash's grid cell size
    double IPC_kStiffness = 1.0e14;
    double IPC_B3Stiffness = 200; // in order to move a point to the desired position, we define an energy

    double MLS_radius = 1.5; // the support radius of MLS points
    int MLS_num_MLS_Pts = 6; // number of MLS points inside a tetrahedron

    double ABD_Coeff = 1.0e10; // penalty coefficient for ABD

};


// Struct of material
struct Material
{
    std::string name = "mat"; // material's name
	double density = 1000; // particle density
	double E = 1.0E4; // Young's modulus
	double nu = 0.3; //Poisson ratio

	double mu = E / (2.0 * (1.0 + nu)); // lame parameter mu / shear modulus  ::::only for isotropic material
	double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)); // lame parameter lambda  ::::only for isotropic material
	double K = 2.0 / 3.0 * mu + lambda; // bulk modulus  ::::only for isotropic material


	// local damage field parameters
	double thetaf = 480; // the reference stress for conventional phase field method
	double Gf = 3;
	double lch = sqrt(2) * 1.0;
	double HsBar = thetaf * thetaf / 2.0 / E / Gf;
	double Hs = HsBar * lch / (1.0 - HsBar * lch);

	double thetaF = 480; // the reference stress for no-weakening phase field


	// return mapping stress threshold(only for bending stress)
	double bendingStressThreshold = 1.0E6;

	// The minimum force magnitude required to break this material
	double fracture_start_force = 1.0E5;


	void updateDenpendecies()
	{
		mu = E / (2.0 * (1.0 + nu)); // lame parameter mu / shear modulus  ::::only for isotropic material
		lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)); // lame parameter lambda  ::::only for isotropic material
		K = 2.0 / 3.0 * mu + lambda; // bulk modulus  ::::only for isotropic material

		HsBar = thetaf * thetaf / 2.0 / E / Gf;
		Hs = HsBar * lch / (1.0 - HsBar * lch);
	}

    double calHsBar()
    {
        return thetaf * thetaf / 2.0 / E / Gf;
    }

};


// split a line from a text file
std::vector<std::string> split(const std::string& s, const std::string& seperator);

// flaten a 3x3 matrix into a vector
Vector9d flatenMatrix3d(const Eigen::Matrix3d & matrix);

// convert a vector to a cross-product matrix. Eq.4.23 in Kim_Eberle_2022
Eigen::Matrix3d vector2CrossProductMatrix(const Eigen::Vector3d& vec);

// calculate the infinite norm of a vector contanining Eigen::Vector3d
double infiniteNorm(std::vector<Eigen::Vector3d>& vec3d);


Eigen::Vector2i findIntersectionOfTwoNums(int min1, int max1, int min2, int max2);
Eigen::Vector2d findIntersectionOfTwoNums(double min1, double max1, double min2, double max2);

 
bool insideBoundingBox(Eigen::Vector3d pt, Eigen::Vector3d bbx_min, Eigen::Vector3d bbx_max);

bool findIntersectionOfTwoVector3i(Eigen::Vector3i& minCoor1, Eigen::Vector3i& maxCoor1, Eigen::Vector3i& minCoor2, Eigen::Vector3i& maxCoor2, Eigen::Vector3i& intersectMin, Eigen::Vector3i& intersectMax);
bool findIntersectionOfTwoVector3d(Eigen::Vector3d& minCoor1, Eigen::Vector3d& maxCoor1, Eigen::Vector3d& minCoor2, Eigen::Vector3d& maxCoor2, Eigen::Vector3d& intersectMin, Eigen::Vector3d& intersectMax);


// find if an element exists in a multimap
bool existInMultiMap_2(int x, int y, std::map<int, std::map<int, int>>& gridMap);


inline double squaredDouble(double a)
{
	return a * a;
}


// find the bounding box of a vector of Eigen::Vector3d
std::pair<Eigen::Vector3d, Eigen::Vector3d> findBoundingBox_vec(std::vector<Eigen::Vector3d>& pos_node);

static std::string calculateID(Eigen::Vector3i index) 
{
    std::string ID = std::to_string(index[0]) + "#" + std::to_string(index[1]) + "#" + std::to_string(index[2]);
	return ID;
};


// boundary condition of each vertex in the mesh
struct boundaryCondition
{
    // vertex's type: 0) default: without constraint; 1) specify location in each timestep(include fixed point); 2) external force, f_ext = xxx; 
    int type = 0;
    // for tyepe1 particles, the location of each timestep 
    std::vector<Eigen::Vector3d> location;
    // for type2 particles, the applied force magnitude
    std::vector<Eigen::Vector3d> force;
    // for type1 and type2 particles, the minimum and maximum applied timesteps
    Eigen::Vector2i appliedTime = {-99, 9999999};


};


struct contact_Info
{
	// 1st int: 0(PT), 1(EE), 2(PG); 
	// 2nd int: index of P(E1)(P: for ground contact case); 
	// 3rd int: index of T(E2); 
	// 4th int: type; 
	// 5th int: actual involved elements, i.e. PP(2), PE(3), PT(4) and EE(4)  
	std::vector<Vector5i> PG_PG;
	std::vector<Vector5i> PT_PP;
	std::vector<Vector5i> PT_PE;
	std::vector<Vector5i> PT_PT;
	std::vector<Vector5i> EE_EE;

	void clear();
};



template <typename Scalar, int size>
static void makePD(Eigen::Matrix<Scalar, size, size>& symMtr)
{
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, size, size>> eigenSolver(symMtr);
	if (eigenSolver.eigenvalues()[0] >= 0.0) {
		return;
	}
	Eigen::DiagonalMatrix<Scalar, size> D(eigenSolver.eigenvalues());
	int rows = ((size == Eigen::Dynamic) ? symMtr.rows() : size);
	for (int i = 0; i < rows; i++) {
		if (D.diagonal()[i] < 0.0) {
			D.diagonal()[i] = 0.0;
		}
		else {
			break;
		}
	}
	symMtr = eigenSolver.eigenvectors() * D * eigenSolver.eigenvectors().transpose();
}



Eigen::Vector3d randomPointInTetrahedron(const Eigen::Vector3d& V1, const Eigen::Vector3d& V2, const Eigen::Vector3d& V3, const Eigen::Vector3d& V4);


// Build the constant Jx matrix for ABD
// Jx is defined under Eq.2 in "Affine body dynamics fast, stable and intersection-free simulation of stiff materials"
Eigen::Matrix<double, 3, 12> build_Jx_matrix_for_ABD(const Eigen::Vector3d& pos);



// kroneckerProduct of two 3x3 matrices
Matrix9d kroneckerProduct_matrices(const Eigen::Matrix3d& A, const Eigen::Matrix3d& B);














void ComputeEigenvector0(double a00, double a01, double a02, double a11, double a12, double a22, double eval0, Eigen::Vector3d& evec0);

void ComputeOrthogonalComplement(Eigen::Vector3d const& W, Eigen::Vector3d& U, Eigen::Vector3d& V);

void ComputeEigenvector1(double a00, double a01, double a02, double a11, double a12, double a22, Eigen::Vector3d const& evec0, double eval1, Eigen::Vector3d& evec1);




//---------------------------------------------------------------------------
// useful for testing
inline DComplex polinom_2(DComplex x, double a, double b)
{
	//Horner's scheme for x*x + a*x + b
	return x * (x + a) + b;
}

//---------------------------------------------------------------------------
// useful for testing
inline DComplex polinom_3(DComplex x, double a, double b, double c)
{
	//Horner's scheme for x*x*x + a*x*x + b*x + c;
	return x * (x * (x + a) + b) + c;
}

//---------------------------------------------------------------------------
// useful for testing
inline DComplex polinom_4(DComplex x, double a, double b, double c, double d)
{
	//Horner's scheme for x*x*x*x + a*x*x*x + b*x*x + c*x + d;
	return x * (x * (x * (x + a) + b) + c) + d;
}

//---------------------------------------------------------------------------
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] � i*x[2], return 1
unsigned int solveP3(double* x, double a, double b, double c);

//---------------------------------------------------------------------------
// Solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d
// (attention - this function returns dynamically allocated array. It has to be released afterwards)
DComplex* solve_quartic(double a, double b, double c, double d);



#endif