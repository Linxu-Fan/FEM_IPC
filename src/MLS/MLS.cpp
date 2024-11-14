#include "mesh.h"

Eigen::Vector3d MLSPoints::cal_MLS_point_pos(const std::vector<Eigen::Vector3d>& pos_node)
{
    Eigen::Vector3d pos = { 0, 0, 0 };
    for (int i = 0; i < index_node.size(); i++)
    {
        pos += weight[i] * pos_node[index_node[i]];
    }

    return pos;
}

double MLSPoints::ComputeGaussianWeight(double distance_squared, double kernelSize)
{
	return exp(-distance_squared / (kernelSize * kernelSize));
}

Eigen::Vector3d MLSPoints::ComputeWeightedCentroid(const std::vector<Eigen::Vector3d>& pos_node)
{
    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
    for (int i = 0; i < index_node.size(); ++i)
    {
        centroid += weight[i] * pos_node[index_node[i]];
    }
    return centroid;
}

Eigen::Matrix3d MLSPoints::ComputeDeformationGradient(const std::vector<Eigen::Vector3d>& pos_node)
{
    int N = index_node.size();

    // Step 2: Compute weighted centroids
    Eigen::Vector3d centroid_def = ComputeWeightedCentroid(pos_node);
    Eigen::Vector3d pos_Defo = cal_MLS_point_pos(pos_node);

    // Step 3: Compute weighted deviations
    std::vector<Eigen::Vector3d> q_ref(N, Eigen::Vector3d::Zero());
    std::vector<Eigen::Vector3d> q_def(N, Eigen::Vector3d::Zero());
    for (int i = 0; i < N; ++i) 
    {
        q_ref[i] = pos_Rest - centroid_ref;
        q_def[i] = pos_Defo - centroid_def;
    }

    // Step 4: Construct weighted matrices for least squares
    Eigen::MatrixXd WQ_ref(N, 3);
    Eigen::MatrixXd WQ_def(N, 3);
    for (int i = 0; i < N; ++i) 
    {
        double sqrt_w = std::sqrt(weight[i]);
        WQ_ref.row(i) = sqrt_w * q_ref[i].transpose();
        WQ_def.row(i) = sqrt_w * q_def[i].transpose();
    }

    // Step 5: Compute deformation gradient F using least squares
    // F = (WQ_def^T * WQ_ref) * (WQ_ref^T * WQ_ref)^-1
    // Alternatively, use SVD for better numerical stability

    // Compute matrices for normal equations
    Eigen::Matrix3d A = WQ_ref.adjoint() * WQ_ref; // 3x3
    Eigen::Matrix3d B = WQ_def.adjoint() * WQ_ref; // 3x3

    // Solve for F
    return B * A.inverse();

}

void MLSPoints::init_MLS(Eigen::Vector3d& pos_Rest_, double volume_, std::vector<int>& index_node_, std::vector<Eigen::Vector3d>& pos_node_Rest, std::string kernel, double kernelSize)
{
    pos_Rest = pos_Rest_;
    pos = pos_Rest_;
	volume = volume_;
	index_node = index_node_;
    centroid_ref = ComputeWeightedCentroid(pos_node_Rest);
	cal_MLS_weights(pos_node_Rest, kernel, kernelSize);

}

void MLSPoints::cal_MLS_weights(std::vector<Eigen::Vector3d>& pos_node_Rest, std::string kernel, double kernelSize)
{
	if (kernel == "Gaussian")
	{
        int N = index_node.size();
        std::vector<int> weights_unnormalized(N, 0.0);
        weight.resize(N, 0.0);
        double weight_sum = 0.0;

        // Compute raw weights
        for (int i = 0; i < N; ++i)
        {
            double dist_sq = (pos_node_Rest[index_node[i]] - pos_Rest).squaredNorm();
            weights_unnormalized[i] = ComputeGaussianWeight(dist_sq, kernelSize);
            weight_sum += weights_unnormalized[i];
        }

        // Normalize weights
        if (weight_sum > 1e-8) 
        { // Prevent division by zero
            for (size_t i = 0; i < N; ++i) 
            {
                weight[i] = weights_unnormalized[i] / weight_sum;
            }
        }
	}
}

void MLSPoints::update_MLS(std::vector<Eigen::Vector3d>& pos_node)
{
    pos = cal_MLS_point_pos(pos_node);
    F = ComputeDeformationGradient(pos_node);
}



