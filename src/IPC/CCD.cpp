#include "CCD.h"

// Check if two edges' boundingbox intersect or not
bool edgeEdgeCCDBroadphase(const Eigen::Vector3d& P1, const Eigen::Vector3d& P2, 
    const Eigen::Vector3d& dP1, const Eigen::Vector3d& dP2, const Eigen::Vector3d& Q1, 
    const Eigen::Vector3d& Q2, const Eigen::Vector3d& dQ1, const Eigen::Vector3d& dQ2, 
    double dist_threshold)
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

// Check if two edges' boundingbox intersect or not (no advection)
bool edgeEdgeCCDBroadphase(const Eigen::Vector3d& P1, const Eigen::Vector3d& P2, 
    const Eigen::Vector3d& Q1, const Eigen::Vector3d& Q2, double dist_threshold)
{
    auto max_1 = P1.array().max(P2.array()).max((P1).array()).max((P2).array());
    auto min_1 = P1.array().min(P2.array()).min((P1).array()).min((P2).array());
    auto max_2 = Q1.array().max(Q2.array()).max((Q1).array()).max((Q2).array());
    auto min_2 = Q1.array().min(Q2.array()).min((Q1).array()).min((Q2).array());
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
bool pointTriangleCCDBroadphase(const Eigen::Vector3d& P, const Eigen::Vector3d& dP, 
    const Eigen::Vector3d& A, const Eigen::Vector3d& dA, const Eigen::Vector3d& B, 
    const Eigen::Vector3d& dB, const Eigen::Vector3d& C, const Eigen::Vector3d& dC, 
    double dist_threshold)
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


// check if a point and a triangle's boundingboxes intersect or not (no advection)
bool pointTriangleCCDBroadphase(const Eigen::Vector3d& P, const Eigen::Vector3d& A, 
    const Eigen::Vector3d& B, const Eigen::Vector3d& C, double dist_threshold)
{
    auto max_P = P.array().max((P).array());
    auto min_P = P.array().min((P).array());
    auto max_T = A.array().max(B.array()).max(C.array()).max((A).array()).max((B).array()).max((C).array());
    auto min_T = A.array().min(B.array()).min(C.array()).min((A).array()).min((B).array()).min((C).array());
    if ((min_P - max_T > dist_threshold).any() || (min_T - max_P > dist_threshold).any())
    {
        return false; // two bounding boxes don't intersect
    }
    else
    {
        return true; // two bounding boxes intersect
    }
}



double edgeEdgeCCDNarrowphase(const Eigen::Vector3d& P1, const Eigen::Vector3d& P2, 
    const Eigen::Vector3d& dP1, const Eigen::Vector3d& dP2, const Eigen::Vector3d& Q1, 
    const Eigen::Vector3d& Q2, const Eigen::Vector3d& dQ1, const Eigen::Vector3d& dQ2, 
    double eta)
{
    Eigen::Vector3d P1_ = P1, P2_ = P2, Q1_ = Q1, Q2_ = Q2, dP1_ = dP1, dP2_ = dP2, dQ1_ = dQ1, dQ2_ = dQ2;
    Eigen::Vector3d mov = (dP1_ + dP2_ + dQ1_ + dQ2_) / 4.0; // Use relative displacement for better convergence
    dP1_ -= mov, dP2_ -= mov, dQ1_ -= mov, dQ2_ -= mov;
    // Suppose these two edges move towards each other 
    double max_disp_mag = std::sqrt(std::max(dP1_.squaredNorm(), dP2_.squaredNorm())) 
        + std::sqrt(std::max(dQ1_.squaredNorm(), dQ2_.squaredNorm()));
    if (max_disp_mag == 0) // No motion
    {
        return 1.0;
    }
    double dist2_cur = 0;
    DIS::computeEdgeEdgeD(P1_, P2_, Q1_, Q2_, dist2_cur);
    double dFunc = dist2_cur;
    if (dFunc <= 0) {
        std::vector<double> dists{ (P1_ - Q1_).squaredNorm(), (P1_ - Q2_).squaredNorm(), 
            (P2_ - Q1_).squaredNorm(), (P2_ - Q2_).squaredNorm() };
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
            std::vector<double> dists{ (P1_ - Q1_).squaredNorm(), 
                (P1_ - Q2_).squaredNorm(), (P2_ - Q1_).squaredNorm(), (P2_ - Q2_).squaredNorm() };
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


double pointTriangleCCDNarrowphase(const Eigen::Vector3d& P, const Eigen::Vector3d& dP, 
    const Eigen::Vector3d& A, const Eigen::Vector3d& dA, const Eigen::Vector3d& B, 
    const Eigen::Vector3d& dB, const Eigen::Vector3d& C, const Eigen::Vector3d& dC, 
    double eta)
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



// get the minimum and maximum corner of a vector of points
// bool advected: if consider the advection of points. If not, direction is ZERO
std::pair<Eigen::Vector3d, Eigen::Vector3d> getMinMaxCorner_boundaryVertices(bool advected, 
    FEMParamters& parameters, Mesh& tetSimMesh, std::vector<Eigen::Vector3d>& direction, 
    double cellSize)
{
    int numOfPts = tetSimMesh.boundaryVertices.size();
    if (advected)
    {
        numOfPts = tetSimMesh.boundaryVertices.size() * 2;
    }

    std::vector<Eigen::Vector3d> pos_node_current_and_next(numOfPts);
    if (advected)
    {
#pragma omp parallel for num_threads(parameters.numOfThreads)
        for (int ft = 0; ft < tetSimMesh.boundaryVertices_vec.size(); ft++)
        {
            int ptInd = tetSimMesh.boundaryVertices_vec[ft];
            pos_node_current_and_next[ft * 2] = tetSimMesh.pos_node[ptInd];
            pos_node_current_and_next[ft * 2 + 1] = tetSimMesh.pos_node[ptInd] + direction[ptInd];
        }
    }
    else
    {
#pragma omp parallel for num_threads(parameters.numOfThreads)
        for (int ft = 0; ft < tetSimMesh.boundaryVertices_vec.size(); ft++)
        {
            int ptInd = tetSimMesh.boundaryVertices_vec[ft];
            pos_node_current_and_next[ft] = tetSimMesh.pos_node[ptInd];
        }
    }

    


    return findBoundingBox_vec(pos_node_current_and_next);
}



// get the minimum and maximum corner of a vector of points of each object
// bool advected: if consider the advection of points. If not, direction is ZERO
// we add another layer of cells outside of the original cells
std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> getMinMaxCorner_boundaryVertices_eachObject(bool advected, 
    FEMParamters& parameters, Mesh& tetSimMesh, std::vector<Eigen::Vector3d>& direction, double cellSize)
{
    
    int numOfPts = tetSimMesh.boundaryVertices.size();
    if (advected)
    {
        numOfPts = tetSimMesh.boundaryVertices.size() * 2;
    }

    std::vector<std::vector<Eigen::Vector3d>> pos_node_current_and_next(tetSimMesh.num_meshes);
    std::map<std::string, int> objectNameIndex = tetSimMesh.tetMeshIndex;


    if (advected)
    {
        for (int ft = 0; ft < tetSimMesh.boundaryVertices_vec.size(); ft++)
        {
            int ptInd = tetSimMesh.boundaryVertices_vec[ft];
            int indexPt = objectNameIndex[tetSimMesh.note_node[ptInd]];
            pos_node_current_and_next[indexPt].push_back(tetSimMesh.pos_node[ptInd]);
            pos_node_current_and_next[indexPt].push_back(tetSimMesh.pos_node[ptInd] + direction[ptInd]);
        }
    }
    else
    {
        for (int ft = 0; ft < tetSimMesh.boundaryVertices_vec.size(); ft++)
        {
            int ptInd = tetSimMesh.boundaryVertices_vec[ft];
            int indexPt = objectNameIndex[tetSimMesh.note_node[ptInd]];
            pos_node_current_and_next[indexPt].push_back(tetSimMesh.pos_node[ptInd]);
        }
    }


    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> res;
    for (int i = 0; i < pos_node_current_and_next.size(); i++)
    {
        std::pair<Eigen::Vector3d, Eigen::Vector3d> bbx = findBoundingBox_vec(pos_node_current_and_next[i]);
        res.push_back(bbx);
    }
    return res;
}



// find the intersection of objects' hash
std::set<std::string> getTheIntersectedHash(bool advected, FEMParamters& parameters, 
    Mesh& tetSimMesh, std::vector<Eigen::Vector3d>& direction, double cellSize, 
    Eigen::Vector3d& minFloorGlobal)
{
    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> intersectVec = 
        getMinMaxCorner_boundaryVertices_eachObject(advected,parameters, tetSimMesh, direction, cellSize);
    std::set<std::string> hashCellNames;
    for (int i = 0; i < intersectVec.size(); i++)
    {
        std::pair<Eigen::Vector3d, Eigen::Vector3d> I1 = intersectVec[i];
        I1.first -= cellSize * Eigen::Vector3d::Ones();
        I1.second += cellSize * Eigen::Vector3d::Ones();

        for (int j = i + 1; j < intersectVec.size(); j++)
        {
            std::pair<Eigen::Vector3d, Eigen::Vector3d> I2 = intersectVec[j];
            I2.first -= cellSize * Eigen::Vector3d::Ones();
            I2.second += cellSize * Eigen::Vector3d::Ones();

            Eigen::Vector3d intersectMin, intersectMax;
            bool success = findIntersectionOfTwoVector3d(I1.first , 
                I1.second, I2.first, I2.second, intersectMin, intersectMax);

            Eigen::Vector3i minFloor = { static_cast<int>(std::floor((intersectMin[0] - minFloorGlobal[0]) / cellSize)) , 
                static_cast<int>(std::floor((intersectMin[1] - minFloorGlobal[1]) / cellSize)) , 
                static_cast<int>(std::floor((intersectMin[2] - minFloorGlobal[2]) / cellSize)) };
            Eigen::Vector3i maxFloor = { static_cast<int>(std::floor((intersectMax[0] - minFloorGlobal[0]) / cellSize)) , 
                static_cast<int>(std::floor((intersectMax[1] - minFloorGlobal[1]) / cellSize)) , 
                static_cast<int>(std::floor((intersectMax[2] - minFloorGlobal[2]) / cellSize)) };

            if (success)
            {
                for (int mnx = minFloor[0]; mnx < maxFloor[0] + 1; mnx++)
                {
                    for (int mny = minFloor[1]; mny < maxFloor[1] + 1; mny++)
                    {
                        for (int mnz = minFloor[2]; mnz < maxFloor[2] + 1; mnz++)
                        {
                            Eigen::Vector3i index = { mnx , mny , mnz };
                            std::string ID = calculateID(index);
                            hashCellNames.insert(ID);

                        }
                    }
                }

            }
        }
    }
    return hashCellNames;
}


void initSpatialHash(bool advected, FEMParamters& parameters, Mesh& tetSimMesh, 
    std::vector<Eigen::Vector3d>& direction, double cellSize, 
    std::vector<spatialHashCellData>& spatialHash_vec, 
    std::map<std::string, int>& hashNameIndex, int timestep = 0)
{
    spatialHash_vec.clear();
    hashNameIndex.clear();

    std::unordered_map<std::string, spatialHashCellData> spatialHash;


    std::pair<Eigen::Vector3d, Eigen::Vector3d> minmax = getMinMaxCorner_boundaryVertices(advected, 
        parameters, tetSimMesh, direction, cellSize);
    Eigen::Vector3d minFloor = minmax.first;
    Eigen::Vector3d maxFloor = minmax.second;

    std::set<std::string> hashCellNames = getTheIntersectedHash(advected, parameters, tetSimMesh, 
        direction, cellSize, minFloor);



    // find spatial hash of boundary vertices
    for (std::map<int, std::set<int>>::iterator it = tetSimMesh.boundaryVertices.begin(); 
        it != tetSimMesh.boundaryVertices.end(); it++)
    {
        std::vector<Eigen::Vector3d> currNext;
        currNext.push_back(tetSimMesh.pos_node[it->first]);
        currNext.push_back(tetSimMesh.pos_node[it->first] + direction[it->first]);


        std::pair<Eigen::Vector3d, Eigen::Vector3d> bbx_cn = findBoundingBox_vec(currNext);
        Eigen::Vector3i minFloor_cn = { static_cast<int>(std::floor((bbx_cn.first[0] - minFloor[0]) / cellSize)) , 
            static_cast<int>(std::floor((bbx_cn.first[1] - minFloor[1]) / cellSize)) , 
            static_cast<int>(std::floor((bbx_cn.first[2] - minFloor[2]) / cellSize)) };
        Eigen::Vector3i maxFloor_cn = { static_cast<int>(std::floor((bbx_cn.second[0] - minFloor[0]) / cellSize)) ,
            static_cast<int>(std::floor((bbx_cn.second[1] - minFloor[1]) / cellSize)) , 
            static_cast<int>(std::floor((bbx_cn.second[2] - minFloor[2]) / cellSize)) };       
        for (int mnx = minFloor_cn[0]; mnx < maxFloor_cn[0] + 1; mnx++)
        {
            for (int mny = minFloor_cn[1]; mny < maxFloor_cn[1] + 1; mny++)
            {
                for (int mnz = minFloor_cn[2]; mnz < maxFloor_cn[2] + 1; mnz++)
                {
                    Eigen::Vector3i index = { mnx , mny , mnz };
                    
                    std::string ID = calculateID(index);

                    spatialHash[ID].bottomLeftCorner = index;

                    spatialHash[ID].vertIndices.insert(it->first);

                    for (std::set<int>::iterator ite = tetSimMesh.boundaryVertices_egde[it->first].begin(); 
                        ite != tetSimMesh.boundaryVertices_egde[it->first].end(); ite++)
                    {
                        spatialHash[ID].edgeIndices.insert(*ite);
                    }

                    for (std::set<int>::iterator itt = tetSimMesh.boundaryVertices[it->first].begin(); 
                        itt != tetSimMesh.boundaryVertices[it->first].end(); itt++)
                    {
                        spatialHash[ID].triaIndices.insert(*itt);
                    }

                }
            }
        }
    }

   
    // give spatial hash string name
    if (parameters.rigidMode == true)
    {
        int indexStart = 0;
        for (std::unordered_map<std::string, spatialHashCellData>::iterator itCell = spatialHash.begin(); 
            itCell != spatialHash.end(); itCell++)
        {
            if (hashCellNames.find(itCell->first) != hashCellNames.end())
            {
                //std::vector<int> edgeIndices_vec(itCell->second.edgeIndices.begin(), itCell->second.edgeIndices.end());
                //itCell->second.edgeIndices_vec = edgeIndices_vec;

                spatialHash_vec.push_back(itCell->second);
                hashNameIndex[itCell->first] = indexStart;               
                indexStart += 1;                
            }
        }
    }
    else
    {
        int indexStart = 0;
        for (std::unordered_map<std::string, spatialHashCellData>::iterator itCell = spatialHash.begin(); 
            itCell != spatialHash.end(); itCell++)
        {
            //std::vector<int> edgeIndices_vec(itCell->second.edgeIndices.begin(), itCell->second.edgeIndices.end());
            //itCell->second.edgeIndices_vec = edgeIndices_vec;

            spatialHash_vec.push_back(itCell->second);
            hashNameIndex[itCell->first] = indexStart;
            indexStart += 1;
        }
    }



    {

        //std::ofstream outfile2("./output/hash.obj", std::ios::trunc);
        //for (int k = 0; k < spatialHash_vec.size(); k++)
        //{
        //    Eigen::Vector3d P1 = (spatialHash_vec[k].bottomLeftCorner.cast<double>()) * cellSize + minFloor;
        //    //if (k == 8)
        //    //{
        //    //    std::cout << "P1 = " << P1 << std::endl;
        //    //    std::cout << "minFloor = " << minFloor << std::endl;
        //    //}

        //    //if (k == 18)
        //    //{
        //    //    std::cout << "P2 = " << P1 << std::endl;
        //    //    //std::cout << "P2 = " << P1 << std::endl;
        //    //}
        //    Eigen::Vector3d P2 = P1;
        //    P2[0] += cellSize;
        //    Eigen::Vector3d P3 = P1;
        //    P3[1] += cellSize;
        //    Eigen::Vector3d P4 = P1;
        //    P4[0] += cellSize;
        //    P4[1] += cellSize;

        //    Eigen::Vector3d P1Q = P1;
        //    P1Q[2] += cellSize;
        //    Eigen::Vector3d P2Q = P2;
        //    P2Q[2] += cellSize;
        //    Eigen::Vector3d P3Q = P3;
        //    P3Q[2] += cellSize;
        //    Eigen::Vector3d P4Q = P4;
        //    P4Q[2] += cellSize;


        //    outfile2 << std::scientific << std::setprecision(8) << "v " << P1[0] << " " << P1[1] << " " << P1[2] << std::endl;
        //    outfile2 << std::scientific << std::setprecision(8) << "v " << P2[0] << " " << P2[1] << " " << P2[2] << std::endl;
        //    outfile2 << std::scientific << std::setprecision(8) << "v " << P3[0] << " " << P3[1] << " " << P3[2] << std::endl;
        //    outfile2 << std::scientific << std::setprecision(8) << "v " << P4[0] << " " << P4[1] << " " << P4[2] << std::endl;
        //    outfile2 << std::scientific << std::setprecision(8) << "v " << P1Q[0] << " " << P1Q[1] << " " << P1Q[2] << std::endl;
        //    outfile2 << std::scientific << std::setprecision(8) << "v " << P2Q[0] << " " << P2Q[1] << " " << P2Q[2] << std::endl;
        //    outfile2 << std::scientific << std::setprecision(8) << "v " << P3Q[0] << " " << P3Q[1] << " " << P3Q[2] << std::endl;
        //    outfile2 << std::scientific << std::setprecision(8) << "v " << P4Q[0] << " " << P4Q[1] << " " << P4Q[2] << std::endl;
        //}
        //outfile2.close();
    }




}


double calMaxStep_spatialHash(FEMParamters& parameters, Mesh& tetSimMesh, 
    std::vector<Eigen::Vector3d>& direction, double cellSize, 
    double dist_threshold, double eta)
{


    //double startTime1, endTime1;
    //startTime1 = omp_get_wtime();


    // build and initialize a spatial hash
    std::vector<spatialHashCellData> spatialHash_vec;
    std::map<std::string, int> hashNameIndex;
    initSpatialHash(true, parameters, tetSimMesh, direction, 
        cellSize, spatialHash_vec, hashNameIndex);


    //endTime1 = omp_get_wtime();
    //std::cout << "Initialization Time is : " << endTime1 - startTime1 << "s" << std::endl;



    //std::cout << "              spatialHash_vec.size() = " << spatialHash_vec.size() << std::endl;

    std::vector<double> stepValue_perHash(spatialHash_vec.size());
#pragma omp parallel for num_threads(parameters.numOfThreads)
    for (int y = 0; y < spatialHash_vec.size(); y++)
    {
        double step = 1.0;
        for (std::set<int>::iterator itP = spatialHash_vec[y].vertIndices.begin(); 
            itP != spatialHash_vec[y].vertIndices.end(); itP++)
        {
            Eigen::Vector3i bottomLeftCorner = spatialHash_vec[y].bottomLeftCorner;
            for (int xx = bottomLeftCorner[0] - 1; xx <= bottomLeftCorner[0] + 1; xx++)
            {
                for (int yy = bottomLeftCorner[1] - 1; yy <= bottomLeftCorner[1] + 1; yy++)
                {
                    for (int zz = bottomLeftCorner[2] - 1; zz <= bottomLeftCorner[2] + 1; zz++)
                    {
                        Eigen::Vector3i index = { xx , yy , zz };
                        std::string ID = calculateID(index);
                        if (hashNameIndex.find(ID) != hashNameIndex.end())
                        {
                            int neigHashIndex = hashNameIndex[ID];
                            for (std::set<int>::iterator itT = spatialHash_vec[neigHashIndex].triaIndices.begin(); 
                                itT != spatialHash_vec[neigHashIndex].triaIndices.end(); itT++)
                            {
                                int vert = *itP, tri = *itT;
                                if (tetSimMesh.boundaryVertices[vert].find(tri) == tetSimMesh.boundaryVertices[vert].end()) // this triangle is not incident with the point
                                {                                   
                                    Eigen::Vector3i triVerts = tetSimMesh.boundaryTriangles[tri];
                                    Eigen::Vector3d P = tetSimMesh.pos_node[vert];
                                    Eigen::Vector3d dP = direction[vert];
                                    Eigen::Vector3d A = tetSimMesh.pos_node[triVerts[0]];
                                    Eigen::Vector3d B = tetSimMesh.pos_node[triVerts[1]];
                                    Eigen::Vector3d C = tetSimMesh.pos_node[triVerts[2]];
                                    Eigen::Vector3d dA = direction[triVerts[0]];
                                    Eigen::Vector3d dB = direction[triVerts[1]];
                                    Eigen::Vector3d dC = direction[triVerts[2]];

                                    if (pointTriangleCCDBroadphase(P, dP, A, dA, B, dB, C, dC, dist_threshold))
                                    {

                                        double step_tmp = pointTriangleCCDNarrowphase(P, dP, A, 
                                            dA, B, dB, C, dC, eta);
                                        step = std::min(step, step_tmp);



                                        //std::cout << "Step VT: V = " << vert << "; T = " << tri << "; Z = " << A[2] <<"; Step = "<< step_tmp << std::endl;


                                    }

                                }
                            }
                        }
                    }
                }
            }
        }

        for (std::set<int>::iterator itE1 = spatialHash_vec[y].edgeIndices.begin(); 
            itE1 != spatialHash_vec[y].edgeIndices.end(); itE1++)
        {
            Eigen::Vector3i bottomLeftCorner = spatialHash_vec[y].bottomLeftCorner;
            for (int xx = bottomLeftCorner[0]; xx <= bottomLeftCorner[0] + 1; xx++)
            {
                for (int yy = bottomLeftCorner[1]; yy <= bottomLeftCorner[1] + 1; yy++)
                {
                    for (int zz = bottomLeftCorner[2]; zz <= bottomLeftCorner[2] + 1; zz++)
                    {
                        Eigen::Vector3i index = { xx , yy , zz };
                        std::string ID = calculateID(index);
                        if (hashNameIndex.find(ID) != hashNameIndex.end())
                        {
                            int neigHashIndex = hashNameIndex[ID];
                            for (std::set<int>::iterator itE2 = spatialHash_vec[neigHashIndex].edgeIndices.begin(); 
                                itE2 != spatialHash_vec[neigHashIndex].edgeIndices.end(); itE2++)
                            {
                                if (*itE1 != *itE2)
                                {
                                    Eigen::Vector2i E1 = tetSimMesh.index_boundaryEdge[*itE1], 
                                        E2 = tetSimMesh.index_boundaryEdge[*itE2];
                                    int P1I = E1[0], P2I = E1[1], Q1I = E2[0], Q2I = E2[1];
                                    if (P1I != Q1I && P1I != Q2I && P2I != Q1I && P2I != Q2I) // not duplicated and incident edges
                                    {
                                        

                                        Eigen::Vector3d P1 = tetSimMesh.pos_node[P1I];
                                        Eigen::Vector3d P2 = tetSimMesh.pos_node[P2I];
                                        Eigen::Vector3d Q1 = tetSimMesh.pos_node[Q1I];
                                        Eigen::Vector3d Q2 = tetSimMesh.pos_node[Q2I];
                                        Eigen::Vector3d dP1 = direction[P1I];
                                        Eigen::Vector3d dP2 = direction[P2I];
                                        Eigen::Vector3d dQ1 = direction[Q1I];
                                        Eigen::Vector3d dQ2 = direction[Q2I];

                                        if (edgeEdgeCCDBroadphase(P1, P2, dP1, dP2, 
                                            Q1, Q2, dQ1, dQ2, dist_threshold))
                                        {

                                            double step_tmp = edgeEdgeCCDNarrowphase(P1, P2, dP1, 
                                                dP2, Q1, Q2, dQ1, dQ2, eta);
                                            step = std::min(step, step_tmp);


                                            //std::cout<<"Hash index = "<< y << "; Step EE: E1 = " << *itE1 << "; E2 = " << *itE2 << "; Step = " << step_tmp << std::endl;


                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        stepValue_perHash[y] = step;
    }

    double minStep = 1.0;
    for (int h = 0; h < stepValue_perHash.size(); h++)
    {
        minStep = std::min(minStep, stepValue_perHash[h]);
    }

    //double endTime2 = omp_get_wtime();
    //std::cout << "Cal step Time is : " << endTime2 - endTime1 << "s" << std::endl;




    return minStep;

}

