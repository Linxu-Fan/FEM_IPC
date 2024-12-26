#include "utils.h"


// calculate the infinite norm of a vector contanining Eigen::Vector3d
double infiniteNorm(std::vector<Eigen::Vector3d>& vec3d)
{
	double norm = 0;
	for (int ii = 0; ii < vec3d.size(); ii++)
	{
		for (int jj = 0; jj < 3; jj++)
		{
			if (norm < abs(vec3d[ii][jj]))
			{
				norm = abs(vec3d[ii][jj]);
			}
		}
	}
	return norm;
}

// split a line from a text file
std::vector<std::string> split(const std::string& s, const std::string& seperator) 
{
	std::vector<std::string> result;
	typedef std::string::size_type string_size;
	string_size i = 0;

	while (i != s.size()) {
		int flag = 0;
		while (i != s.size() && flag == 0) {
			flag = 1;
			for (string_size x = 0; x < seperator.size(); ++x)
				if (s[i] == seperator[x]) {
					++i;
					flag = 0;
					break;
				}
		}

		flag = 0;
		string_size j = i;
		while (j != s.size() && flag == 0) {
			for (string_size x = 0; x < seperator.size(); ++x)
				if (s[j] == seperator[x]) {
					flag = 1;
					break;
				}
			if (flag == 0)
				++j;
		}
		if (i != j) {
			result.push_back(s.substr(i, j - i));
			i = j;
		}
	}
	return result;
}

// flaten a 3x3 matrix into a vector
Vector9d flatenMatrix3d(const Eigen::Matrix3d& matrix)
{
	Vector9d result;
	result[0] = matrix(0, 0);
	result[1] = matrix(1, 0);
	result[2] = matrix(2, 0);
	result[3] = matrix(0, 1);
	result[4] = matrix(1, 1);
	result[5] = matrix(2, 1);
	result[6] = matrix(0, 2);
	result[7] = matrix(1, 2);
	result[8] = matrix(2, 2);

	return result;
}

// convert a vector to a cross-product matrix. Eq.4.23 in Kim_Eberle_2022
Eigen::Matrix3d vector2CrossProductMatrix(const Eigen::Vector3d& vec)
{
	Eigen::Matrix3d result = Eigen::Matrix3d::Zero();
	result(2, 1) = vec[0];
	result(1, 2) = -vec[0];
	result(0, 2) = vec[1];
	result(2, 0) = -vec[1];
	result(1, 0) = vec[2];
	result(0, 1) = -vec[2];
	return result;
}

Eigen::Vector2d findIntersectionOfTwoNums(double min1, double max1, double min2, double max2)
{
	Eigen::Vector2d res = { -999999999999, -999999999999 };
	if (max2 <= min1)
	{
		res = { -999999999999, -999999999999 };
	}
	else if (max2 > min1 && max2 <= max1)
	{
		if (min2 <= min1)
		{
			res = { min1, max2 };
		}
		else
		{
			res = { min2, max2 };
		}
	}
	else
	{
		if (min2 <= min1)
		{
			res = { min1, max1 };
		}
		else if (min2 > min1 && min2 <= max1)
		{
			res = { min2, max1 };
		}
		else
		{
			res = { -999999999999, -999999999999 };
		}
	}

	if (res[0] == res[1])
	{
		res = { -999999999999, -999999999999 };
	}

	return res;
}

Eigen::Vector2i findIntersectionOfTwoNums(int min1, int max1, int min2, int max2)
{
	Eigen::Vector2i res = {-999999999999, -999999999999 };
	if (max2 <= min1)
	{
		res = { -999999999999, -999999999999 };
	}
	else if (max2 > min1 && max2 <= max1)
	{
		if (min2 <= min1)
		{
			res = { min1, max2 };
		}
		else
		{
			res = { min2, max2 };
		}
	}
	else
	{
		if (min2 <= min1)
		{
			res = { min1, max1 };
		}
		else if (min2 > min1 && min2 <= max1)
		{
			res = { min2, max1 };
		}
		else
		{
			res = { -999999999999, -999999999999 };
		}
	}

	if (res[0] == res[1])
	{
		res = { -999999999999, -999999999999 };
	}

	return res;
}

bool findIntersectionOfTwoVector3i(Eigen::Vector3i& minCoor1, Eigen::Vector3i& maxCoor1, 
	Eigen::Vector3i& minCoor2, Eigen::Vector3i& maxCoor2, Eigen::Vector3i& intersectMin, 
	Eigen::Vector3i& intersectMax)
{
	bool success = true;
	for (int i = 0; i < 3; i++)
	{
		Eigen::Vector2i xyz = findIntersectionOfTwoNums(minCoor1[i], maxCoor1[i], minCoor2[i], maxCoor2[i]);
		if (xyz[0] == -999999999999)
		{
			success = false;
			break;
		}
		intersectMin[i] = xyz[0];
		intersectMax[i] = xyz[1];
	}
	
	return success;
}


bool findIntersectionOfTwoVector3d(Eigen::Vector3d& minCoor1, Eigen::Vector3d& maxCoor1, 
	Eigen::Vector3d& minCoor2, Eigen::Vector3d& maxCoor2, Eigen::Vector3d& intersectMin,
	Eigen::Vector3d& intersectMax)
{
	bool success = true;
	for (int i = 0; i < 3; i++)
	{
		Eigen::Vector2d xyz = findIntersectionOfTwoNums(minCoor1[i], maxCoor1[i], minCoor2[i], maxCoor2[i]);
		if (xyz[0] == -999999999999)
		{
			success = false;
			break;
		}
		intersectMin[i] = xyz[0];
		intersectMax[i] = xyz[1];
	}

	return success;
}

// find if an element exists in a multimap
bool existInMultiMap_2(int x, int y, std::map<int, std::map<int, int>>& gridMap)
{
	if (gridMap.find(x) != gridMap.end())
	{
		if (gridMap[x].find(y) != gridMap[x].end())
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}


std::pair<Eigen::Vector3d, Eigen::Vector3d> findBoundingBox_vec(std::vector<Eigen::Vector3d>& pos_node)
{
	double minx = 1.0E10, miny = 1.0E10, minz = 1.0E10, maxx = -1.0E10, maxy = -1.0E10, maxz = -1.0E10;
	for (int i = 0; i < pos_node.size(); i++)
	{
		minx = std::min(minx, pos_node[i][0]);
		miny = std::min(miny, pos_node[i][1]);
		minz = std::min(minz, pos_node[i][2]);

		maxx = std::max(maxx, pos_node[i][0]);
		maxy = std::max(maxy, pos_node[i][1]);
		maxz = std::max(maxz, pos_node[i][2]);
	}
	Eigen::Vector3d min = { minx , miny , minz };
	Eigen::Vector3d max = { maxx , maxy , maxz };

	return std::make_pair(min, max);
}




bool insideBoundingBox(Eigen::Vector3d pt, Eigen::Vector3d bbx_min, Eigen::Vector3d bbx_max)
{
	bool res = false;
	if (pt[0] >= bbx_min[0] && pt[0] <= bbx_max[0])
	{
		if (pt[1] >= bbx_min[1] && pt[1] <= bbx_max[1])
		{
			if (pt[2] >= bbx_min[2] && pt[2] <= bbx_max[2])
			{
				res = true;
			}
		}
	}
	return res;
}




Eigen::Vector3d randomPointInTetrahedron(const Eigen::Vector3d& V1, const Eigen::Vector3d& V2, 
	const Eigen::Vector3d& V3, const Eigen::Vector3d& V4)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0.0, 1.0);

	double r1 = dis(gen);
	double r2 = dis(gen);
	double r3 = dis(gen);
	double r4 = dis(gen);

	double sum = r1 + r2 + r3 + r4;
	double b1 = r1 / sum;
	double b2 = r2 / sum;
	double b3 = r3 / sum;
	double b4 = r4 / sum;

	return b1 * V1 + b2 * V2 + b3 * V3 + b4 * V4;
}





Eigen::Matrix<double, 3, 12> build_Jx_matrix_for_ABD(const Eigen::Vector3d& pos)
{
	Eigen::Matrix<double, 3, 12> Jx = Eigen::Matrix<double, 3, 12>::Zero();
	Jx.block(0, 0, 3, 3) = Eigen::Matrix3d::Identity();
	Jx.block(0, 3, 1, 3) = pos.transpose();
	Jx.block(1, 6, 1, 3) = pos.transpose();
	Jx.block(2, 9, 1, 3) = pos.transpose();
	return Jx;
}



Matrix9d kroneckerProduct_matrices(const Eigen::Matrix3d& A, const Eigen::Matrix3d& B)
{
	Matrix9d result = Matrix9d::Zero();
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			result.block(i * 3, j * 3, 3, 3) = A(i, j) * B;
		}
	}
	return result;
}








// The following code is generated by Huancheng
void ComputeEigenvector0(double a00, double a01, double a02, double a11, 
	double a12, double a22, double eval0, Eigen::Vector3d& evec0)
{
	Eigen::Vector3d row0 = { a00 - eval0, a01, a02 };
	Eigen::Vector3d row1 = { a01, a11 - eval0, a12 };
	Eigen::Vector3d row2 = { a02, a12, a22 - eval0 };
	Eigen::Vector3d r0xr1 = row0.cross(row1);
	Eigen::Vector3d r0xr2 = row0.cross(row2);
	Eigen::Vector3d r1xr2 = row1.cross(row2);
	double d0 = r0xr1.dot(r0xr1);
	double d1 = r0xr2.dot(r0xr2);
	double d2 = r1xr2.dot(r1xr2);

	double dmax = d0;
	int32_t imax = 0;
	if (d1 > dmax) {
		dmax = d1;
		imax = 1;
	}
	if (d2 > dmax) {
		imax = 2;
	}

	if (imax == 0) {
		evec0 = r0xr1 / std::sqrt(d0);
	}
	else if (imax == 1) {
		evec0 = r0xr2 / std::sqrt(d1);
	}
	else {
		evec0 = r1xr2 / std::sqrt(d2);
	}
}


void ComputeOrthogonalComplement(Eigen::Vector3d const& W, Eigen::Vector3d& U, Eigen::Vector3d& V)
{
	double invLength;
	if (std::fabs(W[0]) > std::fabs(W[1])) {
		invLength = (double)1 / std::sqrt(W[0] * W[0] + W[2] * W[2]);
		U = { -W[2] * invLength, (double)0, +W[0] * invLength };
	}
	else {
		invLength = (double)1 / std::sqrt(W[1] * W[1] + W[2] * W[2]);
		U = { (double)0, +W[2] * invLength, -W[1] * invLength };
	}
	V = W.cross(U);
}


void ComputeEigenvector1(double a00, double a01, double a02, double a11, 
	double a12, double a22, Eigen::Vector3d const& evec0, double eval1, Eigen::Vector3d& evec1)
{
	Eigen::Vector3d U, V;
	ComputeOrthogonalComplement(evec0, U, V);

	Eigen::Vector3d AU =
	{
			a00 * U[0] + a01 * U[1] + a02 * U[2],
			a01 * U[0] + a11 * U[1] + a12 * U[2],
			a02 * U[0] + a12 * U[1] + a22 * U[2]
	};

	Eigen::Vector3d AV =
	{
			a00 * V[0] + a01 * V[1] + a02 * V[2],
			a01 * V[0] + a11 * V[1] + a12 * V[2],
			a02 * V[0] + a12 * V[1] + a22 * V[2]
	};

	double m00 = U[0] * AU[0] + U[1] * AU[1] + U[2] * AU[2] - eval1;
	double m01 = U[0] * AV[0] + U[1] * AV[1] + U[2] * AV[2];
	double m11 = V[0] * AV[0] + V[1] * AV[1] + V[2] * AV[2] - eval1;

	double absM00 = std::fabs(m00);
	double absM01 = std::fabs(m01);
	double absM11 = std::fabs(m11);
	double maxAbsComp;
	if (absM00 >= absM11) {
		maxAbsComp = std::max(absM00, absM01);
		if (maxAbsComp > (double)0) {
			if (absM00 >= absM01) {
				m01 /= m00;
				m00 = (double)1 / std::sqrt((double)1 + m01 * m01);
				m01 *= m00;
			}
			else {
				m00 /= m01;
				m01 = (double)1 / std::sqrt((double)1 + m00 * m00);
				m00 *= m01;
			}
			evec1 = m01 * U - m00 * V;
		}
		else {
			evec1 = U;
		}
	}
	else {
		maxAbsComp = std::max(absM11, absM01);
		if (maxAbsComp > (double)0) {
			if (absM11 >= absM01) {
				m01 /= m11;
				m11 = (double)1 / std::sqrt((double)1 + m01 * m01);
				m01 *= m11;
			}
			else {
				m11 /= m01;
				m01 = (double)1 / std::sqrt((double)1 + m11 * m11);
				m11 *= m01;
			}
			evec1 = m11 * U - m01 * V;
		}
		else {
			evec1 = U;
		}
	}
}








//---------------------------------------------------------------------------
// solve cubic equation x^3 + a*x^2 + b*x + c
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] ± i*x[2], return 1
unsigned int solveP3(double* x, double a, double b, double c) {
	double a2 = a * a;
	double q = (a2 - 3 * b) / 9;
	double r = (a * (2 * a2 - 9 * b) + 27 * c) / 54;
	double r2 = r * r;
	double q3 = q * q * q;
	double A, B;
	if (r2 < q3)
	{
		double t = r / sqrt(q3);
		if (t < -1) t = -1;
		if (t > 1) t = 1;
		t = acos(t);
		a /= 3; q = -2 * sqrt(q);
		x[0] = q * cos(t / 3) - a;
		x[1] = q * cos((t + M_2PI) / 3) - a;
		x[2] = q * cos((t - M_2PI) / 3) - a;
		return 3;
	}
	else
	{
		A = -pow(fabs(r) + sqrt(r2 - q3), 1. / 3);
		if (r < 0) A = -A;
		B = (0 == A ? 0 : q / A);

		a /= 3;
		x[0] = (A + B) - a;
		x[1] = -0.5 * (A + B) - a;
		x[2] = 0.5 * sqrt(3.) * (A - B);
		if (fabs(x[2]) < eps) { x[2] = x[1]; return 2; }

		return 1;
	}
}

//---------------------------------------------------------------------------
// Solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d
// (attention - this function returns dynamically allocated array. It has to be released afterwards)
DComplex* solve_quartic(double a, double b, double c, double d)
{
	double a3 = -b;
	double b3 = a * c - 4. * d;
	double c3 = -a * a * d - c * c + 4. * b * d;

	// cubic resolvent
	// y^3 − b*y^2 + (ac−4d)*y − a^2*d−c^2+4*b*d = 0

	double x3[3];
	unsigned int iZeroes = solveP3(x3, a3, b3, c3);

	double q1, q2, p1, p2, D, sqD, y;

	y = x3[0];
	// THE ESSENCE - choosing Y with maximal absolute value !
	if (iZeroes != 1)
	{
		if (fabs(x3[1]) > fabs(y)) y = x3[1];
		if (fabs(x3[2]) > fabs(y)) y = x3[2];
	}

	// h1+h2 = y && h1*h2 = d  <=>  h^2 -y*h + d = 0    (h === q)

	D = y * y - 4 * d;
	if (fabs(D) < eps) //in other words - D==0
	{
		q1 = q2 = y * 0.5;
		// g1+g2 = a && g1+g2 = b-y   <=>   g^2 - a*g + b-y = 0    (p === g)
		D = a * a - 4 * (b - y);
		if (fabs(D) < eps) //in other words - D==0
			p1 = p2 = a * 0.5;

		else
		{
			sqD = sqrt(D);
			p1 = (a + sqD) * 0.5;
			p2 = (a - sqD) * 0.5;
		}
	}
	else
	{
		sqD = sqrt(D);
		q1 = (y + sqD) * 0.5;
		q2 = (y - sqD) * 0.5;
		// g1+g2 = a && g1*h2 + g2*h1 = c       ( && g === p )  Krammer
		p1 = (a * q1 - c) / (q1 - q2);
		p2 = (c - a * q2) / (q1 - q2);
	}

	DComplex* retval = new DComplex[4];

	// solving quadratic eq. - x^2 + p1*x + q1 = 0
	D = p1 * p1 - 4 * q1;
	if (D < 0.0)
	{
		retval[0].real(-p1 * 0.5);
		retval[0].imag(sqrt(-D) * 0.5);
		retval[1] = std::conj(retval[0]);
	}
	else
	{
		sqD = sqrt(D);
		retval[0].real((-p1 + sqD) * 0.5);
		retval[1].real((-p1 - sqD) * 0.5);
	}

	// solving quadratic eq. - x^2 + p2*x + q2 = 0
	D = p2 * p2 - 4 * q2;
	if (D < 0.0)
	{
		retval[2].real(-p2 * 0.5);
		retval[2].imag(sqrt(-D) * 0.5);
		retval[3] = std::conj(retval[2]);
	}
	else
	{
		sqD = sqrt(D);
		retval[2].real((-p2 + sqD) * 0.5);
		retval[3].real((-p2 - sqD) * 0.5);
	}

	return retval;
}




