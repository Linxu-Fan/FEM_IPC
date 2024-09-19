#include "ElasticEnergy.h"


// compute df/dx: the partial derivative of deformation gradient wrt the node position per tetrahedral
Eigen::Matrix<double, 9, 12> ElasticEnergy::dF_wrt_dx(Eigen::Matrix3d& DmInv)
{
	Eigen::Matrix<double, 9, 12> PFPx = Eigen::Matrix<double, 9, 12>::Zero(); // the partial derivative of deformation gradient wrt the node position per tetrahedral. It is vectorized.

	const double m = DmInv(0, 0);
	const double n = DmInv(0, 1);
	const double o = DmInv(0, 2);
	const double p = DmInv(1, 0);
	const double q = DmInv(1, 1);
	const double r = DmInv(1, 2);
	const double s = DmInv(2, 0);
	const double t = DmInv(2, 1);
	const double u = DmInv(2, 2);
	const double t1 = -m - p - s;
	const double t2 = -n - q - t;
	const double t3 = -o - r - u;

	PFPx(0, 0) = t1;
	PFPx(0, 3) = m;
	PFPx(0, 6) = p;
	PFPx(0, 9) = s;
	PFPx(1, 1) = t1;
	PFPx(1, 4) = m;
	PFPx(1, 7) = p;
	PFPx(1, 10) = s;
	PFPx(2, 2) = t1;
	PFPx(2, 5) = m;
	PFPx(2, 8) = p;
	PFPx(2, 11) = s;
	PFPx(3, 0) = t2;
	PFPx(3, 3) = n;
	PFPx(3, 6) = q;
	PFPx(3, 9) = t;
	PFPx(4, 1) = t2;
	PFPx(4, 4) = n;
	PFPx(4, 7) = q;
	PFPx(4, 10) = t;
	PFPx(5, 2) = t2;
	PFPx(5, 5) = n;
	PFPx(5, 8) = q;
	PFPx(5, 11) = t;
	PFPx(6, 0) = t3;
	PFPx(6, 3) = o;
	PFPx(6, 6) = r;
	PFPx(6, 9) = u;
	PFPx(7, 1) = t3;
	PFPx(7, 4) = o;
	PFPx(7, 7) = r;
	PFPx(7, 10) = u;
	PFPx(8, 2) = t3;
	PFPx(8, 5) = o;
	PFPx(8, 8) = r;
	PFPx(8, 11) = u;
	return PFPx;
}

// compute the elastic energy
double ElasticEnergy::Val(Material& mat, std::string model, Eigen::Matrix3d& F, double dt, double vol)
{
	double energy = 0;
	if (model == "neoHookean")
	{
		double J = F.determinant();
		energy = mat.mu / 2.0 * (F.squaredNorm() - 3.0) - mat.mu * log(J) + mat.lambda / 2.0 * log(J) * log(J);
	}
	else if (model == "ARAP")
	{
		double i1 = F.squaredNorm();
		double i2 = (F.transpose() * F).squaredNorm();
		double i3 = (F.transpose() * F).determinant();
		double J = F.determinant();
		double a = 0;
		double b = -2 * i1;
		double c = -8 * J;
		double d = i1 * i1 - 2 * (i1 * i1 - i2);
		std::complex<double>* solutions = solve_quartic(a, b, c, d);
		double f = solutions[2].real();
		delete[] solutions;

		energy = mat.mu * 0.5 * (i1 - 2 * f + 3);
	}
	else if (model == "ARAP_linear")
	{
		double J = F.determinant();
		energy = mat.mu / 2.0 * (13.0 / 8.0) * (13.0 / 8.0) * (F.transpose() - Eigen::Matrix3d::Identity()).squaredNorm();
	}
	else if (model == "ACAP")
	{
		double J = F.determinant();
		energy = mat.mu / 2.0 * (F - Eigen::Matrix3d::Identity()).squaredNorm();
	}

	return dt * dt * energy * vol;

}

// compute the energy gradient dPHI/dF (PK1 stress). Return a vectorized gradient which is a 9x1 matrix
void ElasticEnergy::Grad(std::vector<std::pair<int, double>>& grad_triplet, int& startIndex_grad, Material& mat, std::string model, Eigen::Matrix3d& F, double dt, double vol, Eigen::Matrix<double, 9, 12>& dFdx, Eigen::Vector4i& tetVertInd, Eigen::Vector4i& tetVertInd_BC)
{
	Vector9d grad_tmp;
	if (model == "neoHookean")
	{
		double J = F.determinant();
		Eigen::Matrix3d pJpF = J * F.inverse().transpose();

		pJpF.setZero();
		pJpF.col(0) = F.col(1).cross(F.col(2));
		pJpF.col(1) = F.col(2).cross(F.col(0));
		pJpF.col(2) = F.col(0).cross(F.col(1));

		Eigen::Matrix3d PK1 = mat.mu * (F - 1.0 / J * pJpF) + mat.lambda * log(J) / J * pJpF;
		grad_tmp = flatenMatrix3d(PK1);
		
	}
	else if (model == "ARAP")
	{
		Eigen::Matrix<double, 9, 1> g1;
		g1.block<3, 1>(0, 0) = 2 * F.col(0);
		g1.block<3, 1>(3, 0) = 2 * F.col(1);
		g1.block<3, 1>(6, 0) = 2 * F.col(2);
		Eigen::Matrix3d mat_g2 = 4 * F * F.transpose() * F;
		Eigen::Matrix<double, 9, 1> g2;
		g2.block<3, 1>(0, 0) = mat_g2.col(0);
		g2.block<3, 1>(3, 0) = mat_g2.col(1);
		g2.block<3, 1>(6, 0) = mat_g2.col(2);
		double J = F.determinant();
		Eigen::Matrix<double, 9, 1> gJ;
		gJ.block<3, 1>(0, 0) = F.col(1).cross(F.col(2));
		gJ.block<3, 1>(3, 0) = F.col(2).cross(F.col(0));
		gJ.block<3, 1>(6, 0) = F.col(0).cross(F.col(1));
		Eigen::Matrix<double, 9, 1> g3 = 2 * J * gJ;
		double i1 = F.squaredNorm();
		double i2 = (F.transpose() * F).squaredNorm();
		double i3 = (F.transpose() * F).determinant();
		double a = 0;
		double b = -2 * i1;
		double c = -8 * J;
		double d = i1 * i1 - 2 * (i1 * i1 - i2);

		std::complex<double>* solutions = solve_quartic(a, b, c, d);
		double f = solutions[2].real();
		double f1 = (2 * f * f + 2 * i1) / (4 * pow(f, 3) - 4 * i1 * f - 8 * J);
		double f2 = -2 / (4 * pow(f, 3) - 4 * i1 * f - 8 * J);
		double fJ = (8 * f) / (4 * pow(f, 3) - 4 * i1 * f - 8 * J);
		delete[] solutions;
		grad_tmp = mat.mu * 0.5 * (g1 - 2 * (f1 * g1 + f2 * g2 + fJ * gJ));
	}
	else if (model == "ARAP_linear")
	{
		Eigen::Matrix3d PK1 = mat.mu * (13.0 / 8.0) * (13.0 / 8.0) * (F - Eigen::Matrix3d::Identity());
		grad_tmp = flatenMatrix3d(PK1);
	}
	else if (model == "ACAP")
	{
		Eigen::Matrix3d PK1 = mat.mu * (F - Eigen::Matrix3d::Identity());
		grad_tmp = flatenMatrix3d(PK1);
	}

	Eigen::Matrix<double, 12, 1> engGrad = dt * dt * vol * dFdx.transpose() * grad_tmp;
	for (int m = 0; m < 4; m++)
	{
		int x1_Ind = tetVertInd[m]; // the first vertex index
		for (int xd = 0; xd < 3; xd++)
		{
			grad_triplet[startIndex_grad + m * 3 + xd] = { x1_Ind * 3 + xd, engGrad(m * 3 + xd, 1) };
		}
	}

}

// compute the energy hessian dPHI2/d2F. Return a vectorized gradient which is a 9x9 matrix
void ElasticEnergy::Hess(std::vector<Eigen::Triplet<double>>& hessian_triplet, int& startIndex_hess, Material& mat, std::string model, Eigen::Matrix3d& F, double dt, double vol, Eigen::Matrix<double, 9, 12>& dFdx, Eigen::Vector4i& tetVertInd, Eigen::Vector4i& tetVertInd_BC)
{

	Eigen::Matrix<double, 9, 9> hessian = Eigen::Matrix<double, 9, 9>::Zero();
	if (model == "neoHookean")
	{
		double J = F.determinant();
		Eigen::Matrix3d pJpF = J * F.inverse().transpose();
		pJpF.setZero();
		pJpF.col(0) = F.col(1).cross(F.col(2));
		pJpF.col(1) = F.col(2).cross(F.col(0));
		pJpF.col(2) = F.col(0).cross(F.col(1));
		Vector9d pJpF_vec = flatenMatrix3d(pJpF);



		// 1st term in the last Equation of page_40
		Eigen::Matrix<double, 9, 9> I9 = Eigen::Matrix<double, 9, 9>::Identity();
		hessian += mat.mu * I9;

		// 2nd term in the last Equation of page_40
		double coef_2 = (mat.mu + mat.lambda * (1.0 - log(J))) / J / J;
		hessian += coef_2 * pJpF_vec * pJpF_vec.transpose();


		//std::cout << "coef_2 * pJpF_vec * pJpF_vec.transpose() = "<< coef_2 * pJpF_vec * pJpF_vec.transpose() << std::endl;

		// 3rd term in the last Equation of page_40
		double coef_3 = (mat.lambda * log(J) - mat.mu) / J;
		Eigen::Matrix<double, 9, 9> HJ = Eigen::Matrix<double, 9, 9>::Zero();
		Eigen::Vector3d f0 = F.col(0);
		Eigen::Vector3d f1 = F.col(1);
		Eigen::Vector3d f2 = F.col(2);
		Eigen::Matrix3d f0_mat = vector2CrossProductMatrix(f0);
		Eigen::Matrix3d f1_mat = vector2CrossProductMatrix(f1);
		Eigen::Matrix3d f2_mat = vector2CrossProductMatrix(f2);
		HJ.block<3, 3>(6, 3) = f0_mat;
		HJ.block<3, 3>(3, 6) = -f0_mat;
		HJ.block<3, 3>(0, 6) = f1_mat;
		HJ.block<3, 3>(6, 0) = -f1_mat;
		HJ.block<3, 3>(3, 0) = f2_mat;
		HJ.block<3, 3>(0, 3) = -f2_mat;
		hessian += coef_3 * HJ;

	}
	else if (model == "ARAP")
	{
		double i1 = F.squaredNorm();
		double i2 = (F.transpose() * F).squaredNorm();
		double i3 = (F.transpose() * F).determinant();
		double J = F.determinant();
		double a = 0;
		double b = -2 * i1;
		double c = -8 * J;
		double d = i1 * i1 - 2 * (i1 * i1 - i2);

		std::complex<double>* solutions = solve_quartic(a, b, c, d);
		double x1 = solutions[0].real();
		double x2 = solutions[1].real();
		double x3 = solutions[3].real();
		double x4 = solutions[2].real();
		double sig1 = (x1 + x4) / 2;
		double sig2 = (x2 + x4) / 2;
		double sig3 = (x3 + x4) / 2;
		double f = solutions[2].real();

		delete[] solutions;
		if (sig1 > sig2) {
			std::swap(sig1, sig2);
		}
		if (sig1 > sig3) {
			std::swap(sig1, sig3);
		}
		if (sig2 > sig3) {
			std::swap(sig2, sig3);
		}
		double f1 = (2 * f * f + 2 * i1) / (4 * f * f * f - 4 * i1 * f - 8 * J);
		double f2 = -2 / (4 * f * f * f - 4 * i1 * f - 8 * J);
		double fJ = (8 * f) / (4 * f * f * f - 4 * i1 * f - 8 * J);
		Eigen::Matrix3d g1 = 2 * F;
		Eigen::Matrix3d g2 = 4 * F * F.transpose() * F;
		Eigen::Matrix3d gJ;
		gJ.col(0) = F.col(1).cross(F.col(2));
		gJ.col(1) = F.col(2).cross(F.col(0));
		gJ.col(2) = F.col(0).cross(F.col(1));
		Eigen::Matrix3d R = f1 * g1 + f2 * g2 + fJ * gJ;
		Eigen::Matrix3d S = R.transpose() * F;
		Eigen::Vector3d V0, V1, V2;

		Eigen::Matrix3d V;
		double norm = S(0, 1) * S(0, 1) + S(0, 2) * S(0, 2) + S(1, 2) * S(1, 2);
		if (norm > 0) {
			double q = (S(0, 0) + S(1, 1) + S(2, 2)) / (double)3;
			double b00 = S(0, 0) - q;
			double b11 = S(1, 1) - q;
			double b22 = S(2, 2) - q;
			double p = std::sqrt((b00 * b00 + b11 * b11 + b22 * b22 + norm * (double)2) / (double)6);
			double c00 = b11 * b22 - S(1, 2) * S(1, 2);
			double c01 = S(0, 1) * b22 - S(1, 2) * S(0, 2);
			double c02 = S(0, 1) * S(1, 2) - b11 * S(0, 2);
			double det = (b00 * c00 - S(0, 1) * c01 + S(0, 2) * c02) / (p * p * p);

			double halfDet = det * (double)0.5;
			halfDet = std::min(std::max(halfDet, (double)-1), (double)1);

			if (halfDet >= (double)0) {
				ComputeEigenvector0(S(0, 0), S(0, 1), S(0, 2), S(1, 1), S(1, 2), S(2, 2), sig3, V2);
				ComputeEigenvector1(S(0, 0), S(0, 1), S(0, 2), S(1, 1), S(1, 2), S(2, 2), V2, sig2, V1);
				V0 = V1.cross(V2);
			}
			else {
				ComputeEigenvector0(S(0, 0), S(0, 1), S(0, 2), S(1, 1), S(1, 2), S(2, 2), sig1, V0);
				ComputeEigenvector1(S(0, 0), S(0, 1), S(0, 2), S(1, 1), S(1, 2), S(2, 2), V0, sig2, V1);
				V2 = V0.cross(V1);
			}
			V.col(2) = V0;
			V.col(1) = V1;
			V.col(0) = V2;
		}
		else
		{
			V = Eigen::Matrix3d::Identity();
		}
		Eigen::Matrix3d Sigma = Eigen::Matrix3d::Identity();
		Sigma(0, 0) = sig3;
		Sigma(1, 1) = sig2;
		Sigma(2, 2) = sig1;

		Eigen::Matrix3d U = R * V;

		Eigen::Matrix3d T0{
			{0, -1, 0},
			{ 1, 0,  0 },
			{ 0, 0,  0 }
		};
		T0 = (1 / sqrt(2)) * U * T0 * V.transpose();
		Eigen::Matrix3d T1{
			{0, 0, 0},
			{ 0, 0,  1 },
			{ 0, -1, 0 }
		};
		T1 = (1 / sqrt(2)) * U * T1 * V.transpose();
		Eigen::Matrix3d T2{
			{0, 0, 1},
			{ 0,  0, 0 },
			{ -1, 0, 0 }
		};
		T2 = (1 / sqrt(2)) * U * T2 * V.transpose();
		Eigen::Matrix<double, 9, 1> t0, t1, t2;
		t0.block<3, 1>(0, 0) = T0.col(0);
		t0.block<3, 1>(3, 0) = T0.col(1);
		t0.block<3, 1>(6, 0) = T0.col(2);
		t1.block<3, 1>(0, 0) = T1.col(0);
		t1.block<3, 1>(3, 0) = T1.col(1);
		t1.block<3, 1>(6, 0) = T1.col(2);
		t2.block<3, 1>(0, 0) = T2.col(0);
		t2.block<3, 1>(3, 0) = T2.col(1);
		t2.block<3, 1>(6, 0) = T2.col(2);
		double sx = Sigma(0, 0);
		double sy = Sigma(1, 1);
		double sz = Sigma(2, 2);
		double lambda0 = 2 / (sx + sy);
		double lambda1 = 2 / (sz + sy);
		double lambda2 = 2 / (sx + sz);

		if (sx + sy < 2)lambda0 = 1;
		if (sz + sy < 2)lambda1 = 1;
		if (sx + sz < 2)lambda2 = 1;
		Eigen::Matrix<double, 9, 9> SH = Eigen::Matrix<double, 9, 9>::Identity();
		SH -= lambda0 * (t0 * t0.transpose());
		SH -= lambda1 * (t1 * t1.transpose());
		SH -= lambda2 * (t2 * t2.transpose());
		hessian = mat.mu * SH;
	}
	else if (model == "ARAP_linear")
	{
		Eigen::Matrix<double, 9, 9> I9 = Eigen::Matrix<double, 9, 9>::Identity();
		hessian = mat.mu * (13.0 / 8.0) * (13.0 / 8.0) * I9;
	}
	else if (model == "ACAP")
	{
		Eigen::Matrix<double, 9, 9> I9 = Eigen::Matrix<double, 9, 9>::Identity();
		hessian = mat.mu * I9;
	}

	Eigen::Matrix<double, 12, 12> engHess = dt * dt * vol * dFdx.transpose() * hessian * dFdx;
	makePD<double, 12>(engHess);
	for (int m = 0; m < 4; m++)
	{
		engHess.col(m * 3).setZero();
		engHess.col(m * 3 + 1).setZero();
		engHess.col(m * 3 + 2).setZero();

		engHess.row(m * 3).setZero();
		engHess.row(m * 3 + 1).setZero();
		engHess.row(m * 3 + 2).setZero();
	}

	int ac = 0;
	for (int m = 0; m < 4; m++)
	{
		int x1_Ind = tetVertInd[m]; // the first vertex index
		for (int n = 0; n < 4; n++)
		{
			int x2_Ind = tetVertInd[n]; // the second vertex index	
			for (int xd = 0; xd < 3; xd++)
			{
				for (int yd = 0; yd < 3; yd++)
				{
					double value = engHess(m * 3 + xd, n * 3 + yd);
					Eigen::Triplet<double> pa = { x1_Ind * 3 + xd, x2_Ind * 3 + yd, value };
					hessian_triplet[startIndex_hess + ac] = pa;
					ac += 1;
				}
			}			
		}
	}


}
