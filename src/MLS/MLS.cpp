#include "mesh.h"


double MLSPoints::computeWeight(const Eigen::Vector3d& x_s, const Eigen::Vector3d& x_i, double h)
{
    double r = (x_s - x_i).norm();
    return exp(-(r / h) * (r / h));
}

Eigen::Vector3d MLSPoints::computeWeightDerivative(const Eigen::Vector3d& x_s, const Eigen::Vector3d& x_i, double h)
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

void MLSPoints::MLS_approximation(const std::vector<Vector3d>& pos_node_Rest, const std::vector<Vector3d>& pos_node, double h)
{

    const size_t N = index_node.size();
    assert(N == x_i.size());

    // Initialize weights
    std::vector<double> w_i(N);
    double w_sum = 0.0;

    // Compute weights
    for (size_t i = 0; i < N; ++i) 
    {
        w_i[i] = computeWeight(pos_node_Rest[index_node[i]], pos_Rest, h);
        w_sum += w_i[i];
    }

    // Compute weighted centroids
    Eigen::Vector3d p_bar = Eigen::Vector3d::Zero();
    Eigen::Vector3d x_bar = Eigen::Vector3d::Zero();

    for (size_t i = 0; i < N; ++i) 
    {
        p_bar += w_i[i] * pos_node_Rest[index_node[i]];
        x_bar += w_i[i] * pos_node[index_node[i]];
    }
    p_bar /= w_sum;
    x_bar /= w_sum;

    // Compute centered coordinates
    std::vector<Eigen::Vector3d> p_tilde(N);
    std::vector<Eigen::Vector3d> x_tilde(N);
    for (size_t i = 0; i < N; ++i) 
    {
        p_tilde[i] = pos_node_Rest[index_node[i]] - p_bar;
        x_tilde[i] = pos_node[index_node[i]] - x_bar;
    }
    Eigen::Vector3d p_s_tilde = pos_Rest - p_bar;

    // Assemble moment matrices
    Eigen::Matrix3d M_pp = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d M_xp = Eigen::Matrix3d::Zero();

    for (size_t i = 0; i < N; ++i) 
    {
        M_pp += w_i[i] * p_tilde[i] * p_tilde[i].transpose();
        M_xp += w_i[i] * x_tilde[i] * p_tilde[i].transpose();
    }

    // Compute deformation gradient F
    Eigen::Matrix3d M_pp_inv = M_pp.inverse();
    F = M_xp * M_pp_inv;

    // Approximate deformed position x_s
    pos = F * p_s_tilde + x_bar;

    // Compute derivative of F with respect to x_i
    std::vector<std::vector<Eigen::Matrix3d>> dF_dxi(N, std::vector<Eigen::Matrix3d>(3)); // 3 components for each x_i

    // Precompute M_pp_inv * p_tilde[i] for each i
    std::vector<Eigen::Vector3d> Mpp_inv_p_tilde(N);
    for (size_t i = 0; i < N; ++i) 
    {
        Mpp_inv_p_tilde[i] = M_pp_inv * p_tilde[i];
    }

    dFdx.clear();
    dFdx.resize(index_node.size());

    for (size_t i = 0; i < N; ++i) 
    {
        // For each component alpha (0, 1, 2 corresponding to x, y, z)
        for (int alpha = 0; alpha < 3; ++alpha) 
        {
            Eigen::Matrix3d dF_dxi_alpha = Eigen::Matrix3d::Zero();

            // Compute w_i * ¦Ä_{a¦Á} * [M_pp_inv * p_tilde_i]_b
            // ¦Ä_{a¦Á} is 1 if a == alpha, 0 otherwise
            for (int a = 0; a < 3; ++a) 
            {
                double delta_a_alpha = (a == alpha) ? 1.0 : 0.0;
                double factor = w_i[i] * delta_a_alpha;

                // Compute the derivative column vector
                Eigen::Vector3d derivative_column = factor * Mpp_inv_p_tilde[i];

                // Assign to the matrix
                for (int b = 0; b < 3; ++b) 
                {
                    dF_dxi_alpha(a, b) = derivative_column(b);
                }
            }

            // Store the derivative tensor component
            Vector9d fv = flatenMatrix3d(dF_dxi_alpha);
            dFdx[i].col(alpha) = fv;

            //dF_dxi[i][alpha] = dF_dxi_alpha;
        }
    }


}

void MLSPoints::init_MLS(Eigen::Vector3d& pos_Rest_, double volume_, std::vector<int>& index_node_, std::string kernel, double h)
{
    pos_Rest = pos_Rest_;
    pos = pos_Rest_;
	volume = volume_;
	index_node = index_node_;
    dFdx.resize(index_node_.size());
}



