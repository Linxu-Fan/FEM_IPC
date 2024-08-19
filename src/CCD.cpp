#include "CCD.h"

// Check if two edges' boundingbox intersect or not
bool edgeEdgeCCDBroadphase(const Eigen::Vector3d& P1, const Eigen::Vector3d& P2, const Eigen::Vector3d& dP1, 
    const Eigen::Vector3d& dP2, const Eigen::Vector3d& Q1, const Eigen::Vector3d& Q2, const Eigen::Vector3d& dQ1, const Eigen::Vector3d& dQ2, double dist_threshold)
{
    auto max_1 = P1.array().max(P2.array()).max((P1 + dP1).array()).max((P2 + dP2).array());
    auto min_1 = P1.array().min(P2.array()).min((P1 + dP1).array()).min((P2 + dP2).array());
    auto max_2 = Q1.array().max(Q2.array()).max((Q1 + dQ1).array()).max((Q2 + dQ2).array());
    auto min_2 = Q1.array().min(Q2.array()).min((Q1 + dQ1).array()).min((Q2 + dQ2).array());
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
bool pointTriangleCCDBroadphase(const Eigen::Vector3d& P, const Eigen::Vector3d& dP, const Eigen::Vector3d& A, 
    const Eigen::Vector3d& dA, const Eigen::Vector3d& B, const Eigen::Vector3d& dB, const Eigen::Vector3d& C, const Eigen::Vector3d& dC, double dist_threshold)
{
    auto max_P = P.array().max((P + dP).array());
    auto min_P = P.array().min((P + dP).array());
    auto max_T = A.array().max(B.array()).max(C.array()).max((A + dA).array()).max((B + dB).array()).max((C + dC).array());
    auto min_T = A.array().min(B.array()).min(C.array()).min((A + dA).array()).min((B + dB).array()).min((C + dC).array());
    if ((min_P - max_T > dist_threshold).any() || (min_T - max_P > dist_threshold).any())
    {
        return false; // two bounding boxes don't intersect
    }
    else
    {
        return true; // two bounding boxes intersect
    }
}


double edgeEdgeCCDNarrowphase(const Eigen::Vector3d& P1, const Eigen::Vector3d& P2, const Eigen::Vector3d& dP1, 
    const Eigen::Vector3d& dP2, const Eigen::Vector3d& Q1, const Eigen::Vector3d& Q2, const Eigen::Vector3d& dQ1, const Eigen::Vector3d& dQ2, double eta)
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
    double dist2_cur = 0;
    DIS::computeEdgeEdgeD(P1_, P2_, Q1_, Q2_, dist2_cur);
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
        
        DIS::computeEdgeEdgeD(P1_, P2_, Q1_, Q2_, dist2_cur);

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


double pointTriangleCCDNarrowphase(const Eigen::Vector3d& P, const Eigen::Vector3d& dP, const Eigen::Vector3d& A, 
    const Eigen::Vector3d& dA, const Eigen::Vector3d& B, const Eigen::Vector3d& dB, const Eigen::Vector3d& C, const Eigen::Vector3d& dC, double eta)
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
    double dist_cur2 = 0;
    DIS::computePointTriD(P_, A_, B_, C_, dist_cur2);
    double dist_cur = std::sqrt(dist_cur2);
    double gap = eta * dist_cur;
    double toc = 0.0;
    while (true) 
    {
        double toc_lower_bound = (1 - eta) * dist_cur / max_disp_mag;
        P_ += toc_lower_bound * dP_;
        A_ += toc_lower_bound * dA_;
        B_ += toc_lower_bound * dB_;
        C_ += toc_lower_bound * dC_;
        DIS::computePointTriD(P_, A_, B_, C_, dist_cur2);
        dist_cur = std::sqrt(dist_cur2);
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



void initSpatialHash(Mesh& tetMesh, std::vector<Eigen::Vector3d>& direction, double cellSize, std::unordered_map<int, spatialHashCellData>& spatialHash)
{
    spatialHash.clear();

    std::vector<Eigen::Vector3d> pos_node_current_and_next;
    for (std::map<int, std::set<int>>::iterator it = tetMesh.boundaryVertices.begin(); it != tetMesh.boundaryVertices.end(); it++)
    {
        pos_node_current_and_next.push_back(tetMesh.pos_node[it->first]);
        pos_node_current_and_next.push_back(tetMesh.pos_node[it->first] + direction[it->first]);
    }
    std::pair<Eigen::Vector3d, Eigen::Vector3d> bbx = findBoundingBox_vec(pos_node_current_and_next);
    Eigen::Vector3i minFloor = { static_cast<int>(std::floor(bbx.first[0] / cellSize)) ,
        static_cast<int>(std::floor(bbx.first[1] / cellSize)) , static_cast<int>(std::floor(bbx.first[2] / cellSize)) };
    Eigen::Vector3i maxFloor = { static_cast<int>(std::floor(bbx.second[0] / cellSize)) ,
        static_cast<int>(std::floor(bbx.second[1] / cellSize)) , static_cast<int>(std::floor(bbx.second[2] / cellSize)) };
    Eigen::Vector3i numCell = maxFloor - minFloor + Eigen::Vector3i::Constant(1);

    // find spatial hash of boundary vertices
    for (std::map<int, std::set<int>>::iterator it = tetMesh.boundaryVertices.begin(); it != tetMesh.boundaryVertices.end(); it++)
    {
        std::vector<Eigen::Vector3d> currNext;
        currNext.push_back(tetMesh.pos_node[it->first]);
        currNext.push_back(tetMesh.pos_node[it->first] + direction[it->first]);

        std::pair<Eigen::Vector3d, Eigen::Vector3d> bbx_cn = findBoundingBox_vec(currNext);
        Eigen::Vector3i minFloor_cn = { static_cast<int>(std::floor((bbx_cn.first[0] - minFloor[0]) / cellSize)) , 
            static_cast<int>(std::floor((bbx_cn.first[1] - minFloor[1]) / cellSize)) , static_cast<int>(std::floor((bbx_cn.first[2] - minFloor[2]) / cellSize)) };
        Eigen::Vector3i maxFloor_cn = { static_cast<int>(std::floor((bbx_cn.second[0] - minFloor[0]) / cellSize)) ,
            static_cast<int>(std::floor((bbx_cn.second[1] - minFloor[1]) / cellSize)) , static_cast<int>(std::floor((bbx_cn.second[2] - minFloor[2]) / cellSize)) };
        for (int mnx = minFloor_cn[0]; mnx < maxFloor_cn[0] + 1; mnx++)
        {
            for (int mny = minFloor_cn[1]; mny < maxFloor_cn[1] + 1; mny++)
            {
                for (int mnz = minFloor_cn[2]; mnz < maxFloor_cn[2] + 1; mnz++)
                {
                    Eigen::Vector3i index = { mnx , mny , mnz };
                    int ID = calculateID(numCell, index);
                    spatialHash[ID].vertIndices.insert(it->first);
                }
            }
        }
    }

    // find spatial hash of edges
    for (std::map<int, Eigen::Vector2i>::iterator it1 = tetMesh.index_boundaryEdge.begin(); it1 != tetMesh.index_boundaryEdge.end(); it1++)
    {
        int p1 = it1->second[0], p2 = it1->second[1];
        std::vector<Eigen::Vector3d> currNext;
        currNext.push_back(tetMesh.pos_node[p1]);
        currNext.push_back(tetMesh.pos_node[p1] + direction[p1]);
        currNext.push_back(tetMesh.pos_node[p2]);
        currNext.push_back(tetMesh.pos_node[p2] + direction[p2]);

        std::pair<Eigen::Vector3d, Eigen::Vector3d> bbx_cn = findBoundingBox_vec(currNext);
        Eigen::Vector3i minFloor_cn = { static_cast<int>(std::floor((bbx_cn.first[0] - minFloor[0]) / cellSize)) ,
            static_cast<int>(std::floor((bbx_cn.first[1] - minFloor[1]) / cellSize)) , static_cast<int>(std::floor((bbx_cn.first[2] - minFloor[2]) / cellSize)) };
        Eigen::Vector3i maxFloor_cn = { static_cast<int>(std::floor((bbx_cn.second[0] - minFloor[0]) / cellSize)) ,
            static_cast<int>(std::floor((bbx_cn.second[1] - minFloor[1]) / cellSize)) , static_cast<int>(std::floor((bbx_cn.second[2] - minFloor[2]) / cellSize)) };
        for (int mnx = minFloor_cn[0]; mnx < maxFloor_cn[0] + 1; mnx++)
        {
            for (int mny = minFloor_cn[1]; mny < maxFloor_cn[1] + 1; mny++)
            {
                for (int mnz = minFloor_cn[2]; mnz < maxFloor_cn[2] + 1; mnz++)
                {
                    Eigen::Vector3i index = { mnx , mny , mnz };
                    int ID = calculateID(numCell, index);
                    spatialHash[ID].edgeIndices.insert(it1->first);
                }
            }
        }
    }

    // find spatial hash of triangles
    for (int i = 0; i < tetMesh.boundaryTriangles.size(); i++)
    {
        int p1 = tetMesh.boundaryTriangles[i][0], p2 = tetMesh.boundaryTriangles[i][1], p3 = tetMesh.boundaryTriangles[i][2];
        std::vector<Eigen::Vector3d> currNext;
        currNext.push_back(tetMesh.pos_node[p1]);
        currNext.push_back(tetMesh.pos_node[p1] + direction[p1]);
        currNext.push_back(tetMesh.pos_node[p2]);
        currNext.push_back(tetMesh.pos_node[p2] + direction[p2]);
        currNext.push_back(tetMesh.pos_node[p3]);
        currNext.push_back(tetMesh.pos_node[p3] + direction[p3]);

        std::pair<Eigen::Vector3d, Eigen::Vector3d> bbx_cn = findBoundingBox_vec(currNext);
        Eigen::Vector3i minFloor_cn = { static_cast<int>(std::floor((bbx_cn.first[0] - minFloor[0]) / cellSize)) ,
            static_cast<int>(std::floor((bbx_cn.first[1] - minFloor[1]) / cellSize)) , static_cast<int>(std::floor((bbx_cn.first[2] - minFloor[2]) / cellSize)) };
        Eigen::Vector3i maxFloor_cn = { static_cast<int>(std::floor((bbx_cn.second[0] - minFloor[0]) / cellSize)) ,
            static_cast<int>(std::floor((bbx_cn.second[1] - minFloor[1]) / cellSize)) , static_cast<int>(std::floor((bbx_cn.second[2] - minFloor[2]) / cellSize)) };
        for (int mnx = minFloor_cn[0]; mnx < maxFloor_cn[0] + 1; mnx++)
        {
            for (int mny = minFloor_cn[1]; mny < maxFloor_cn[1] + 1; mny++)
            {
                for (int mnz = minFloor_cn[2]; mnz < maxFloor_cn[2] + 1; mnz++)
                {
                    Eigen::Vector3i index = { mnx , mny , mnz };
                    int ID = calculateID(numCell, index);
                    spatialHash[ID].triaIndices.insert(i);
                }
            }
        }
    }

}


double calMaxStep_spatialHash(Mesh& tetMesh, std::vector<Eigen::Vector3d>& direction, double cellSize, double dist_threshold, double eta)
{
    // build and initialize a spatial hash
    std::unordered_map<int, spatialHashCellData> spatialHash;
    initSpatialHash(tetMesh, direction, cellSize, spatialHash);


    double step = 1.0;

    std::map<int, std::map<int, bool>> PTCal; // if vert (1st int) and triangle (2nd int) pair has been calculated or not
    std::map<int, std::map<int, bool>> EECal; // if edge1 (1st int) and edge2 (2nd int) pair has been calculated or not
    for (std::unordered_map<int, spatialHashCellData>::iterator itCell = spatialHash.begin(); itCell != spatialHash.end(); itCell++)
    {
        int cellIndex = itCell->first;


        // calcualte the PT pair
        for (std::set<int>::iterator itP = itCell->second.vertIndices.begin(); itP != itCell->second.vertIndices.end(); itP++)
        {
            for (std::set<int>::iterator itT = itCell->second.triaIndices.begin(); itT != itCell->second.triaIndices.end(); itT++)
            {
                int vert = *itP, tri = *itT;
                if (PTCal[vert][tri] == false)
                {
                    if (tetMesh.boundaryVertices[vert].find(tri) == tetMesh.boundaryVertices[vert].end()) // this triangle is not incident with the point
                    {
                        Eigen::Vector3i triVerts = tetMesh.boundaryTriangles[tri];

                        Eigen::Vector3d P = tetMesh.pos_node[vert];
                        Eigen::Vector3d dP = direction[vert];
                        Eigen::Vector3d A = tetMesh.pos_node[triVerts[0]];
                        Eigen::Vector3d B = tetMesh.pos_node[triVerts[1]];
                        Eigen::Vector3d C = tetMesh.pos_node[triVerts[2]];
                        Eigen::Vector3d dA = direction[triVerts[0]];
                        Eigen::Vector3d dB = direction[triVerts[1]];
                        Eigen::Vector3d dC = direction[triVerts[2]];

                        if (pointTriangleCCDBroadphase(P, dP, A, dA, B, dB, C, dC, dist_threshold))
                        {
                            double step_tmp = pointTriangleCCDNarrowphase(P, dP, A, dA, B, dB, C, dC, dist_threshold);
                            step = std::min(step , step_tmp);
                        }
                        PTCal[vert][tri] = true;
                    }                
                }
            }
        }
    
        // calcualte the EE pair
        for (std::set<int>::iterator itE1 = itCell->second.edgeIndices.begin(); itE1 != itCell->second.edgeIndices.end(); itE1++)
        {
            for (std::set<int>::iterator itE2 = itCell->second.edgeIndices.begin(); itE2 != itCell->second.edgeIndices.end(); itE2++)
            {
                if (*itE1 != *itE2)
                {
                    if (EECal[*itE1][*itE2] == false)
                    {
                        Eigen::Vector2i E1 = tetMesh.index_boundaryEdge[*itE1], E2 = tetMesh.index_boundaryEdge[*itE2];
                        int P1I = E1[0], P2I = E1[1], Q1I = E2[0], Q2I = E2[1];
                        if (P1I != Q1I && P1I != Q2I && P2I != Q1I && P2I != Q2I) // not duplicated and incident edges
                        {
                            Eigen::Vector3d P1 = tetMesh.pos_node[P1I];
                            Eigen::Vector3d P2 = tetMesh.pos_node[P2I];
                            Eigen::Vector3d Q1 = tetMesh.pos_node[Q1I];
                            Eigen::Vector3d Q2 = tetMesh.pos_node[Q2I];
                            Eigen::Vector3d dP1 = direction[P1I];
                            Eigen::Vector3d dP2 = direction[P2I];
                            Eigen::Vector3d dQ1 = direction[Q1I];
                            Eigen::Vector3d dQ2 = direction[Q2I];

                            if (edgeEdgeCCDBroadphase(P1, P2, dP1, dP2, Q1, Q2, dQ1, dQ2, eta))
                            {
                                double step_tmp = edgeEdgeCCDNarrowphase(P1, P2, dP1, dP2, Q1, Q2, dQ1, dQ2, eta);
                                step = std::min(step, step_tmp);
                            }
                        }

                        EECal[*itE1][*itE2] = true;
                    }
                }
            }
        }

    }

    return step;

}

