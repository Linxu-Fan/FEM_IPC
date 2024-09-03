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
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <map>
#include <unordered_map>
#include <string.h>
#include <cassert>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <sstream>
#include <cstdio>
#include <cassert>
#include <assert.h>
#include <cstring>
#include <cfloat>
#include <complex>

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

    bool rigidMode = false; // if rigidMode, IPC will only check contact pairs that belongs to two different objects
    std::vector<std::string> objectNames; // name of objects in the scene
	double IPC_dis = 0.01; // the distance of IPC gap, i.e., barrier energy emerges if the distance is smaller than this threshold
	double IPC_eta = 0.1; // the distance ratio that bring two edges or point-triangle to
	double IPC_hashSize = 0.1; // the spatial hash's grid cell size
    double IPC_kStiffness = 1.0e14;

};



// Struct of material
struct Material
{
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


	void updateDenpendecies()
	{
		mu = E / (2.0 * (1.0 + nu)); // lame parameter mu / shear modulus  ::::only for isotropic material
		lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)); // lame parameter lambda  ::::only for isotropic material
		K = 2.0 / 3.0 * mu + lambda; // bulk modulus  ::::only for isotropic material

		HsBar = thetaf * thetaf / 2.0 / E / Gf;
		Hs = HsBar * lch / (1.0 - HsBar * lch);
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
    // vertex's type: 0) default: without constraint; 1) fixed points: velocity = 0; 2) external force, f_ext = xxx
    int type = 0;
    // 1st and 2nd element are the starting and ending timestep when a boundary condition is applied
    Eigen::Vector2i appliedTime = { 0, 100000 };
    // for type2 particles, the applied force magnitude
    Eigen::Vector3d force = { 0,0,0 };
};



struct BE_Grad_Hess
{
    std::vector<Eigen::Triplet<double>> hessian_triplet;
    std::vector<std::pair<int, double>> grad_triplet;
};

void BE_to_triplet(std::vector<Eigen::Triplet<double>>& hessian_triplet, std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_hess, int& startIndex_grad, std::vector<boundaryCondition>& boundaryCondition_node, int& D1Index, Vector3d& V3, Matrix3d& H3x3);
void BE_to_triplet(std::vector<Eigen::Triplet<double>>& hessian_triplet, std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_hess, int& startIndex_grad, std::vector<boundaryCondition>& boundaryCondition_node, Eigen::Vector2i& D2Index, Vector6d& V6, Matrix6d& H6x6);
void BE_to_triplet(std::vector<Eigen::Triplet<double>>& hessian_triplet, std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_hess, int& startIndex_grad, std::vector<boundaryCondition>& boundaryCondition_node, Eigen::Vector3i& D3Index, Vector9d& V9, Matrix9d& H9x9);
void BE_to_triplet(std::vector<Eigen::Triplet<double>>& hessian_triplet, std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_hess, int& startIndex_grad, std::vector<boundaryCondition>& boundaryCondition_node, Eigen::Vector4i& D4Index, Vector12d& V12, Matrix12d& H12x12);

struct BarrierEnergyElement
{
    int PT_EE_or_Ground = 0; // PT: 0; EE: 1; Ground: 2
    int size = 1;

    Eigen::Vector4i Indices;

    int D1Index; 
    Eigen::Vector2i D2Index;
    Eigen::Vector3i D3Index;
    Eigen::Vector4i D4Index;

    Vector3d V3;
    Vector6d V6;
    Vector9d V9;
    Vector12d V12;

    Matrix3d H3x3; 
    Matrix6d H6x6; 
    Matrix9d H9x9; 
    Matrix12d H12x12; 

};


class BarrierEnergyRes
{
public:
    std::vector<std::vector<Eigen::Triplet<double>>> hessian_triplet_vec;
    std::vector<std::vector<std::pair<int, double>>> grad_triplet_vec;

    void clear();
};




// project a matrix into semi-positive definite
Matrix3d projToSPD(Matrix3d& top);

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


// The following code is generated by Huancheng
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