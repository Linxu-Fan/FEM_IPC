#include "CCD.h"

// Check if two edges' boundingbox intersect or not
bool edgeEdgeCCDBroadphase(const Eigen::Vector3d& P1, const Eigen::Vector3d& P2, const Eigen::Vector3d& P1_Next, const Eigen::Vector3d& P2_Next, const Eigen::Vector3d& Q1, const Eigen::Vector3d& Q2, const Eigen::Vector3d& Q1_Next, const Eigen::Vector3d& Q2_Next, double dist_threshold)
{
    auto max_1 = P1.array().max(P2.array()).max((P1_Next).array()).max((P2_Next).array());
    auto min_1 = P1.array().min(P2.array()).min((P1_Next).array()).min((P2_Next).array());
    auto max_2 = Q1.array().max(Q2.array()).max((Q1_Next).array()).max((Q2_Next).array());
    auto min_2 = Q1.array().min(Q2.array()).min((Q1_Next).array()).min((Q2_Next).array());
    if ((min_1 - max_2 > dist_threshold).any() || (min_2 - max_1 > dist_threshold).any())
    {
        return false; // two bounding boxes don't intersect
    } 
    else
    {
        return true; // two bounding boxes intersect
    }
        
}


// check if a point and a triangle's boundingboxes intersect or not
bool pointTriangleCCDBroadphase(const Eigen::Vector3d& P, const Eigen::Vector3d& P_Next, const Eigen::Vector3d& A, const Eigen::Vector3d& A_Next, const Eigen::Vector3d& B, const Eigen::Vector3d& B_Next, const Eigen::Vector3d& C, const Eigen::Vector3d& C_Next, double dist_threshold)
{
    auto max_P = P.array().max((P_Next).array());
    auto min_P = P.array().min((P_Next).array());
    auto max_T = A.array().max(B.array()).max(C.array()).max((A_Next).array()).max((B_Next).array()).max((C_Next).array());
    auto min_T = A.array().min(B.array()).min(C.array()).min((A_Next).array()).min((B_Next).array()).min((C_Next).array());
    if ((min_P - max_T > dist_threshold).any() || (min_T - max_P > dist_threshold).any())
    {
        return false; // two bounding boxes don't intersect
    }
    else
    {
        return true; // two bounding boxes intersect
    }
}


bool edgeEdgeCCDNarrowphase(const Eigen::Vector3d& P1, const Eigen::Vector3d& P2, const Eigen::Vector3d& dP1, const Eigen::Vector3d& dP2, const Eigen::Vector3d& Q1, const Eigen::Vector3d& Q2, const Eigen::Vector3d& dQ1, const Eigen::Vector3d& dQ2, double eta)
{
    Eigen::Vector3d P1_ = P1, P2_ = P2, Q1_ = Q1, Q2_ = Q2, dP1_ = dP1, dP2_ = dP2, dQ1_ = dQ1, dQ2_ = dQ2;
    Eigen::Vector3d mov = (dP1_ + dP2_ + dQ1_ + dQ2_) / 4.0; // Use relative displacement for better convergence
    dP1_ -= mov, dP2_ -= mov, dQ1_ -= mov, dQ2_ -= mov;
    // Suppose these two edges move towards each other 
    double max_disp_mag = std::sqrt(std::max(dP1_.squaredNorm(), dP2_.squaredNorm())) + std::sqrt(std::max(dQ1_.squaredNorm(), dQ2_.squaredNorm()));
    if (max_disp_mag == 0) // No motion
    {
        return 1.0;
    }
    double dist2_cur = edgeEdgeDis2(P1_, P2_, Q1_, Q2_);
    double dFunc = dist2_cur;
    if (dFunc <= 0) {
        std::vector<double> dists{ (P1_ - Q1_).squaredNorm(), (P1_ - Q2_).squaredNorm(), (P2_ - Q1_).squaredNorm(), (P2_ - Q2_).squaredNorm() };
        dist2_cur = *std::min_element(dists.begin(), dists.end());
        dFunc = dist2_cur;
    }
    double dist_cur = std::sqrt(dist2_cur);
    double gap = eta * dFunc / (dist_cur);
    double toc = 0.0;
    while (true) 
    {
        double toc_lower_bound = (1 - eta) * dFunc / ((dist_cur) * max_disp_mag);
        P1_ += toc_lower_bound * dP1_;
        P2_ += toc_lower_bound * dP2_;
        Q1_ += toc_lower_bound * dQ1_;
        Q2_ += toc_lower_bound * dQ2_;
        dist2_cur = edgeEdgeDis2(P1_, P2_, Q1_, Q2_);
        dFunc = dist2_cur;
        if (dFunc <= 0) {
            std::vector<double> dists{ (P1_ - Q1_).squaredNorm(), (P1_ - Q2_).squaredNorm(), (P2_ - Q1_).squaredNorm(), (P2_ - Q2_).squaredNorm() };
            dist2_cur = *std::min_element(dists.begin(), dists.end());
            dFunc = dist2_cur;
        }
        dist_cur = std::sqrt(dist2_cur);
        if (toc && (dFunc / (dist_cur) < gap)) {
            break;
        }
        toc += toc_lower_bound;
        if (toc > 1.0)
        {
            return 1.0;
        }          
    }
    return toc;
}

bool pointTriangleCCDNarrowphase(const Eigen::Vector3d& P, const Eigen::Vector3d& dP, const Eigen::Vector3d& A, const Eigen::Vector3d& dA, const Eigen::Vector3d& B, const Eigen::Vector3d& dB, const Eigen::Vector3d& C, const Eigen::Vector3d& dC, double eta)
{
    Eigen::Vector3d P_ = P, A_ = A, B_ = B, C_ = C, dP_ = dP, dA_ = dA, dB_ = dB, dC_ = dC;
    Eigen::Vector3d mov = (dA_ + dB_ + dC_ + dP_) / 4;
    dA_ -= mov, dB_ -= mov, dC_ -= mov, dP_ -= mov;
    std::vector<double> disp_mag2_vec{ dA_.squaredNorm(), dB_.squaredNorm(), dC_.squaredNorm() };
    double max_disp_mag = dP_.norm() + std::sqrt(*std::max_element(disp_mag2_vec.begin(), disp_mag2_vec.end()));
    if (max_disp_mag == 0)
    {
        return 1.0;
    }
        
    double dist_cur = std::sqrt(pointTriangleDis2(P_, A_, B_, C_));
    double gap = eta * dist_cur;
    double toc = 0.0;
    while (true) 
    {
        double toc_lower_bound = (1 - eta) * dist_cur / max_disp_mag;
        P_ += toc_lower_bound * dP_;
        A_ += toc_lower_bound * dA_;
        B_ += toc_lower_bound * dB_;
        C_ += toc_lower_bound * dC_;
        dist_cur = std::sqrt(pointTriangleDis2(P_, A_, B_, C_));
        if (toc && (dist_cur < gap))
        {
            break;
        }

        toc += toc_lower_bound;
        if (toc > 1.0) 
        {
            return 1.0;
        }
    }
    return toc;
    
}

