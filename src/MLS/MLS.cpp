#include "mesh.h"


double MLSPoints::computeWeight(const Eigen::Vector3d& x_s, 
    const Eigen::Vector3d& x_i, double h)
{
    double r = (x_s - x_i).norm();
    return exp(-(r / h) * (r / h));
}

Eigen::Vector3d MLSPoints::computeWeightDerivative(const Eigen::Vector3d& x_s, 
    const Eigen::Vector3d& x_i, double h)
{
    Eigen::Vector3d diff = x_s - x_i;
    double r = diff.norm();
    if (r == 0)
    {
        return Eigen::Vector3d::Zero();
    }
    double w = computeWeight(x_s, x_i, h);
    double w_prime = (-2.0 * w) / (h * h);
    return w_prime * diff;
}

void MLSPoints::MLS_approximation(const std::vector<Vector3d>& pos_node_Rest, 
    const std::vector<Vector3d>& pos_node, double& radius)
{

    // old approximation
    {
        //const size_t N = index_node.size();
        //assert(N == x_i.size());

        //// Initialize weights
        //std::vector<double> w_i(N);
        //double w_sum = 0.0;

        //// Compute weights
        //for (size_t i = 0; i < N; ++i)
        //{
        //    w_i[i] = computeWeight(pos_node_Rest[index_node[i]], pos_Rest, radius);
        //    w_sum += w_i[i];
        //}

        //// Compute weighted centroids
        //Eigen::Vector3d p_bar = Eigen::Vector3d::Zero();
        //Eigen::Vector3d x_bar = Eigen::Vector3d::Zero();

        //for (size_t i = 0; i < N; ++i)
        //{
        //    p_bar += w_i[i] * pos_node_Rest[index_node[i]];
        //    x_bar += w_i[i] * pos_node[index_node[i]];
        //}
        //p_bar /= w_sum;
        //x_bar /= w_sum;

        //// Compute centered coordinates
        //std::vector<Eigen::Vector3d> p_tilde(N);
        //std::vector<Eigen::Vector3d> x_tilde(N);
        //for (size_t i = 0; i < N; ++i)
        //{
        //    p_tilde[i] = pos_node_Rest[index_node[i]] - p_bar;
        //    x_tilde[i] = pos_node[index_node[i]] - x_bar;
        //}
        //Eigen::Vector3d p_s_tilde = pos_Rest - p_bar;

        //// Assemble moment matrices
        //Eigen::Matrix3d M_pp = Eigen::Matrix3d::Zero();
        //Eigen::Matrix3d M_xp = Eigen::Matrix3d::Zero();

        //for (size_t i = 0; i < N; ++i)
        //{
        //    M_pp += w_i[i] * p_tilde[i] * p_tilde[i].transpose();
        //    M_xp += w_i[i] * x_tilde[i] * p_tilde[i].transpose();
        //}

        //// Compute deformation gradient F
        //Eigen::Matrix3d M_pp_inv = M_pp.inverse();
        //F = M_xp * M_pp_inv;

        //// Approximate deformed position x_s
        //pos = F * p_s_tilde + x_bar;

        //// Compute derivative of F with respect to x_i
        //std::vector<std::vector<Eigen::Matrix3d>> dF_dxi(N, std::vector<Eigen::Matrix3d>(3)); // 3 components for each x_i

        //// Precompute M_pp_inv * p_tilde[i] for each i
        //std::vector<Eigen::Vector3d> Mpp_inv_p_tilde(N);
        //for (size_t i = 0; i < N; ++i)
        //{
        //    Mpp_inv_p_tilde[i] = M_pp_inv * p_tilde[i];
        //}

        //dFdx.clear();
        //dFdx.resize(index_node.size());

        //for (size_t i = 0; i < N; ++i)
        //{
        //    // For each component alpha (0, 1, 2 corresponding to x, y, z)
        //    for (int alpha = 0; alpha < 3; ++alpha)
        //    {
        //        Eigen::Matrix3d dF_dxi_alpha = Eigen::Matrix3d::Zero();

        //        // Compute w_i * ¦Ä_{a¦Á} * [M_pp_inv * p_tilde_i]_b
        //        // ¦Ä_{a¦Á} is 1 if a == alpha, 0 otherwise
        //        for (int a = 0; a < 3; ++a)
        //        {
        //            double delta_a_alpha = (a == alpha) ? 1.0 : 0.0;
        //            double factor = w_i[i] * delta_a_alpha;

        //            // Compute the derivative column vector
        //            Eigen::Vector3d derivative_column = factor * Mpp_inv_p_tilde[i];

        //            // Assign to the matrix
        //            for (int b = 0; b < 3; ++b)
        //            {
        //                dF_dxi_alpha(a, b) = derivative_column(b);
        //            }
        //        }




        //        // Store the derivative tensor component
        //        Vector9d fv = flatenMatrix3d(dF_dxi_alpha);

        //        //std::cout << "fv: " << fv.transpose() << std::endl;


        //        dFdx[i].col(alpha) = fv;

        //        //std::cout << "dFdx[i]: " << dFdx[i] << std::endl;
        //        //std::cout << "dFdx[i]: " << dFdx[i] << std::endl;
        //    }
        //}

    }

    // new approximation
    {


        // Number of points
        int N = index_node.size();
        // Dimension
        int d = 3;




        // Reference positions p_i
        std::vector<Eigen::Vector3d> p_list;
        // Deformed positions x_i
        std::vector<Eigen::Vector3d> x_list;

        for (int h = 0; h < index_node.size(); h++)
        {
            p_list.push_back(pos_node_Rest[index_node[h]]);
            x_list.push_back(pos_node[index_node[h]]);
        }


        // Fill p_list and x_list with your data

        // Position of point s in the reference configuration
        Eigen::VectorXd p_s = pos_Rest;

        // Smoothing length for the weight function
        double h = radius;

        // Step 1: Compute weights
        Eigen::VectorXd w(N);
        for (int i = 0; i < N; ++i) {
            double r = (p_list[i] - p_s).norm();
            w(i) = computeWeight(p_list[i], p_s, radius);
        }

        // Step 2: Compute shifted basis functions
        int m = d + 1; // Shifted linear basis
        Eigen::MatrixXd P(N, m);
        for (int i = 0; i < N; ++i) {
            P(i, 0) = 1.0;
            for (int j = 0; j < d; ++j) {
                P(i, j + 1) = p_list[i](j) - p_s(j);
            }
        }
        // Basis at p_s (shifted basis at the evaluation point)
        Eigen::VectorXd P_s(m);
        P_s(0) = 1.0;
        for (int j = 0; j < d; ++j) {
            P_s(j + 1) = 0.0;
        }

        // Step 3: Compute matrices A and B
        Eigen::MatrixXd W = w.asDiagonal();
        Eigen::MatrixXd A = P.transpose() * W * P;
        // Construct X matrix
        Eigen::MatrixXd X(N, d);
        for (int i = 0; i < N; ++i) {
            X.row(i) = x_list[i].transpose();
        }
        Eigen::MatrixXd B = P.transpose() * W * X;

        // Step 4: Solve for a
        Eigen::MatrixXd a = A.ldlt().solve(B); // Solving A * a = B

        // Step 5: Approximate x_s
        // Since P_s = [1; 0; ..., 0], x_s = a.row(0)
        Eigen::VectorXd x_s = a.row(0); // First row of a


        ////////////////////
        pos = x_s;
        ////////////////////



        // Step 6: Compute deformation gradient F_s
        Eigen::MatrixXd F_s = a.block(1, 0, d, d).transpose(); // F_s = a_{1:d}^T



        ////////////////////
        F = F_s;
        ////////////////////


        // Step 7: Compute derivative of F_s with respect to each x_i
        std::vector<Eigen::MatrixXd> dFdx_list(N); // List to store derivatives (size d^2 x d)
        Eigen::MatrixXd A_inv = A.inverse();
        // For performance, consider solving systems instead of computing the inverse

        for (int i = 0; i < N; ++i) {
            Eigen::VectorXd P_i = P.row(i).transpose(); // m x 1 vector
            Eigen::VectorXd v = w(i) * A_inv * P_i;     // m x 1 vector
            // Exclude the first component (corresponding to a_0)
            Eigen::VectorXd v_p = v.segment(1, d);     // v_p in the derivation (size d)

            // Initialize derivative matrix D_i of size d^2 x d
            Eigen::MatrixXd D_i(d * d, d);
            D_i.setZero();
            // Fill D_i according to the formula
            for (int n = 0; n < d; ++n) {
                for (int p = 0; p < d; ++p) {
                    int row = p * d + n;
                    D_i(row, n) = v_p(p); // Note the absence of negative sign since derivatives are directly v_p
                }
            }
            dFdx_list[i] = D_i;

            ////////////////////
            dFdx[i] = D_i;
            ////////////////////

        }


    }

    // quadratic approximation
    {
        //// Number of points
        //int N = index_node.size();
        //// Dimension
        //int d = 3;

        //// Reference positions p_i
        //std::vector<Eigen::Vector3d> p_list(N);
        //// Deformed positions x_i
        //std::vector<Eigen::Vector3d> x_list(N);

        //for (int h = 0; h < index_node.size(); h++)
        //{
        //    p_list.push_back(pos_node_Rest[index_node[h]]);
        //    x_list.push_back(pos_node[index_node[h]]);
        //}



        //// Fill p_list and x_list with your data

        //// Position of point s in the reference configuration
        //Eigen::VectorXd p_s = pos_Rest;

        //// Smoothing length for the weight function
        //double h = radius;

        //// Step 1: Compute weights
        //Eigen::VectorXd w(N);
        //for (int i = 0; i < N; ++i) {
        //    double r = (p_list[i] - p_s).norm();
        //    w(i) = computeWeight(p_list[i], p_s, radius);
        //}

        //// Step 2: Compute quadratic shifted basis functions
        //int m = 1 + d + d * (d + 1) / 2; // Number of basis functions for quadratic basis
        //Eigen::MatrixXd P(N, m);
        //for (int i = 0; i < N; ++i) {
        //    Eigen::VectorXd P_i = computeQuadraticShiftedBasis(p_list[i], p_s);
        //    P.row(i) = P_i.transpose();
        //}
        //// Basis at p_s (shifted basis at the evaluation point)
        //Eigen::VectorXd P_s = Eigen::VectorXd::Zero(m);
        //P_s(0) = 1.0; // Since (p_i - p_s) = 0 at p_s, only the constant term is non-zero

        //// Step 3: Compute matrices A and B
        //Eigen::MatrixXd W = w.asDiagonal();
        //Eigen::MatrixXd A = P.transpose() * W * P;
        //// Construct X matrix
        //Eigen::MatrixXd X(N, d);
        //for (int i = 0; i < N; ++i) {
        //    X.row(i) = x_list[i].transpose();
        //}
        //Eigen::MatrixXd B = P.transpose() * W * X;

        //// Step 4: Solve for a
        //Eigen::MatrixXd a = A.ldlt().solve(B); // Solving A * a = B

        //// Step 5: Approximate x_s
        //// Since P_s = [1; 0; ..., 0], x_s = a.row(0)
        //Eigen::VectorXd x_s = a.row(0); // First row of a


        ////////////////////////
        //pos = x_s;
        ////////////////////////



        //// Step 6: Compute deformation gradient F_s
        //// F_s = (G_s^T) * a
        //Eigen::MatrixXd G_s_T = computeBasisGradientAtPs(d, m).transpose(); // d x m
        //Eigen::MatrixXd F_s = G_s_T * a; // d x d


        ////////////////////////
        //F = F_s;
        ////////////////////////


        //// Step 7: Compute derivative of F_s with respect to each x_i
        //std::vector<Eigen::MatrixXd> dFdx_list(N); // List to store derivatives (each is d^2 x d)
        //Eigen::MatrixXd A_inv = A.inverse();
        //// For better performance, consider solving A * y = P_i instead of computing A_inv

        //for (int i = 0; i < N; ++i) {
        //    Eigen::VectorXd P_i = P.row(i).transpose(); // m x 1 vector
        //    // Compute v = w_i * A_inv * P_i
        //    Eigen::VectorXd y = A.ldlt().solve(P_i);
        //    Eigen::VectorXd v = w(i) * y; // m x 1 vector

        //    // Compute the product G_s^T * v
        //    Eigen::VectorXd GsT_v = G_s_T * v; // d x 1 vector

        //    // Initialize derivative matrix D_i of size d^2 x d
        //    Eigen::MatrixXd D_i(d * d, d);
        //    D_i.setZero();
        //    // Fill D_i according to the derivative formula
        //    for (int p = 0; p < d; ++p) {
        //        for (int q = 0; q < d; ++q) {
        //            for (int n = 0; n < d; ++n) {
        //                if (q == n) {
        //                    int row = p * d + q;
        //                    D_i(row, n) = GsT_v(p);
        //                }
        //            }
        //        }
        //    }
        //    dFdx_list[i] = D_i;

        //    ////////////////////
        //    dFdx[i] = D_i;
        //    ////////////////////

        //}

    }
    

}

void MLSPoints::init_MLS(Eigen::Vector3d& pos_Rest_, double volume_, 
    std::vector<int>& index_node_, std::string kernel, double& h)
{
    pos_Rest = pos_Rest_;
    pos = pos_Rest_;
	volume = volume_;
	index_node = index_node_;
    dFdx.resize(index_node_.size());
}














// new MLS approximation
// Compute the derivative of the basis functions (shifted linear basis)
Eigen::MatrixXd MLSPoints::computeBasisGradient(int d) 
{
    // For shifted linear basis, gradient is [0; I_d]
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(d + 1, d);
    G.block(1, 0, d, d) = Eigen::MatrixXd::Identity(d, d);
    return G;
}




// Function to construct the quadratic shifted basis for a point
Eigen::VectorXd MLSPoints::computeQuadraticShiftedBasis(const Eigen::VectorXd& p_i, const Eigen::VectorXd& p_s) {
    int d = p_i.size(); // Dimension (2D or 3D)
    int m = 1 + d + d * (d + 1) / 2;
    Eigen::VectorXd P(m);
    P(0) = 1.0; // Constant term

    // Linear terms
    Eigen::VectorXd p_shifted = p_i - p_s; // Shifted position
    for (int i = 0; i < d; ++i) {
        P(1 + i) = p_shifted(i);
    }

    // Quadratic terms
    int index = 1 + d;
    for (int i = 0; i < d; ++i) {
        for (int j = i; j < d; ++j) {
            P(index++) = p_shifted(i) * p_shifted(j);
        }
    }

    return P;
}

// Function to compute the gradient of the quadratic shifted basis at p_s
Eigen::MatrixXd MLSPoints::computeBasisGradientAtPs(int d, int m) {
    // Gradient of basis functions w.r.t p at p = p_s
    // Only linear terms have non-zero gradients
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(m, d);
    // Linear terms
    G.block(1, 0, d, d) = Eigen::MatrixXd::Identity(d, d);
    // Quadratic terms
    // At p = p_s, the derivative of quadratic terms is zero
    return G;
}
