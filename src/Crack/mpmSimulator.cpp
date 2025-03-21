﻿//#include "mpmSimulator.h"
//
//
//namespace mpmSimulator
//{
//	
//	std::string calculateID_string(int x, int y, int z) // string id of the particle
//	{
//		return std::to_string(x) + "#" + std::to_string(y) + "#" + std::to_string(z);
//	};
//
//	// calculate particles' weights and find neighbouring nodes
//	void calWeightsAndNodes(std::vector<MPMParticle>& particles, MPMParamters& param, std::vector<Grid>& nodesVec, std::map<std::string, int>& gridMap)
//	{
//		// store the key and value of gridMap: active grid
//
//		Eigen::Vector3i minCellIndex = { 100000000, 100000000, 100000000 };
//		Eigen::Vector3i maxCellIndex = { -100000000, -100000000, -100000000 };
//#pragma omp parallel for num_threads(param.numOfThreads)
//		for (int f = 0; f < particles.size(); f++)
//		{
//			extractCrackSurface::weightAndDreri  WD = extractCrackSurface::calWeight(param.dx, particles[f].position);
//
//			particles[f].posIndex = WD.ppIndex;
//			particles[f].weight = WD.weight;
//			particles[f].deltaWeight = WD.deltaWeight;
//
//			particles[f].supportNodes.clear();
//			particles[f].supportNodeWeight.clear();
//			particles[f].supportNodeDeltaWeight.clear();
//
//
//			minCellIndex[0] = std::min(WD.ppIndex[0], minCellIndex[0]);
//			maxCellIndex[0] = std::max(WD.ppIndex[0], maxCellIndex[0]);
//
//			minCellIndex[1] = std::min(WD.ppIndex[1], minCellIndex[1]);
//			maxCellIndex[1] = std::max(WD.ppIndex[1], maxCellIndex[1]);
//
//			minCellIndex[2] = std::min(WD.ppIndex[2], minCellIndex[2]);
//			maxCellIndex[2] = std::max(WD.ppIndex[2], maxCellIndex[2]);
//		};
//
//		// find parameters that used to devide the background grid into disconnected patches
//		int blockSize = 6; // the size of each block. The minimum size for quadratic kernel is 5!!!!!!!!!!!!!!!!
//		std::map<std::string, std::vector<int>> bottomLeftBloack;
//		std::map<std::string, std::vector<int>> bottomRightBloack;
//		std::map<std::string, std::vector<int>> topLeftBloack;
//		std::map<std::string, std::vector<int>> topRightBloack;
//		Eigen::Vector3i span = { maxCellIndex[0] - minCellIndex[0], maxCellIndex[1] - minCellIndex[1], maxCellIndex[2] - minCellIndex[2] };
//		int minSpan = std::min(span[0], std::min(span[1], span[2]));
//		if (minSpan == span[0])
//		{
//			for (int f = 0; f < particles.size(); f++)
//			{
//				int remainZ = (particles[f].posIndex[2] - minCellIndex[2]) % (2 * blockSize);
//				int remainY = (particles[f].posIndex[1] - minCellIndex[1]) % (2 * blockSize);
//
//				int numZ = (particles[f].posIndex[2] - minCellIndex[2]) / (2 * blockSize);
//				int numY = (particles[f].posIndex[1] - minCellIndex[1]) / (2 * blockSize);
//				std::string blockID = std::to_string(numZ) + "#" + std::to_string(numY);
//				if (remainZ >= 0 && remainZ < blockSize) // left 
//				{
//					if (remainY >= 0 && remainY < blockSize) // bottom left
//					{
//						bottomLeftBloack[blockID].push_back(f);
//					}
//					else // top left
//					{
//						topLeftBloack[blockID].push_back(f);
//					}
//				}
//				else // right
//				{
//					if (remainY >= 0 && remainY < blockSize) // top right
//					{
//						bottomRightBloack[blockID].push_back(f);
//					}
//					else // bottm right
//					{
//						topRightBloack[blockID].push_back(f);
//					}
//				}
//			}
//		}
//		else if (minSpan == span[1])
//		{
//			for (int f = 0; f < particles.size(); f++)
//			{
//				int remainX = (particles[f].posIndex[0] - minCellIndex[0]) % (2 * blockSize);
//				int remainZ = (particles[f].posIndex[2] - minCellIndex[2]) % (2 * blockSize);
//
//				int numX = (particles[f].posIndex[0] - minCellIndex[0]) / (2 * blockSize);
//				int numZ = (particles[f].posIndex[2] - minCellIndex[2]) / (2 * blockSize);
//				std::string blockID = std::to_string(numX) + "#" + std::to_string(numZ);
//				if (remainX >= 0 && remainX < blockSize) // left 
//				{
//					if (remainZ >= 0 && remainZ < blockSize) // bottom left
//					{
//						bottomLeftBloack[blockID].push_back(f);
//					}
//					else // top left
//					{
//						topLeftBloack[blockID].push_back(f);
//					}
//				}
//				else // right
//				{
//					if (remainZ >= 0 && remainZ < blockSize) // top right
//					{
//						bottomRightBloack[blockID].push_back(f);
//					}
//					else // bottm right
//					{
//						topRightBloack[blockID].push_back(f);
//					}
//				}
//			}
//		}
//		else
//		{
//			for (int f = 0; f < particles.size(); f++)
//			{
//				int remainX = (particles[f].posIndex[0] - minCellIndex[0]) % (2 * blockSize);
//				int remainY = (particles[f].posIndex[1] - minCellIndex[1]) % (2 * blockSize);
//
//				int numX = (particles[f].posIndex[0] - minCellIndex[0]) / (2 * blockSize);
//				int numY = (particles[f].posIndex[1] - minCellIndex[1]) / (2 * blockSize);
//				std::string blockID = std::to_string(numX) + "#" + std::to_string(numY);
//				if (remainX >= 0 && remainX < blockSize) // left 
//				{
//					if (remainY >= 0 && remainY < blockSize) // bottom left
//					{
//						bottomLeftBloack[blockID].push_back(f);
//					}
//					else // top left
//					{
//						topLeftBloack[blockID].push_back(f);
//					}
//				}
//				else // right
//				{
//					if (remainY >= 0 && remainY < blockSize) // bottom right
//					{
//						bottomRightBloack[blockID].push_back(f);
//					}
//					else // top right
//					{
//						topRightBloack[blockID].push_back(f);
//					}
//				}
//			}
//
//
//
//		}
//
//
//		int count = nodesVec.size() - 1; // count the number of active grid node
//		for (int f = 0; f < particles.size(); f++)
//		{
//			for (int i = 0; i < 3; i++)
//			{
//				for (int j = 0; j < 3; j++)
//				{
//					for (int k = 0; k < 3; k++)
//					{
//						std::string ID = calculateID_string(particles[f].posIndex[0] + i, particles[f].posIndex[1] + j, particles[f].posIndex[2] + k);
//
//						if (gridMap.find(ID) == gridMap.end())
//						{
//							count += 1;
//							gridMap[ID] = count;
//
//							Grid node;
//							nodesVec.push_back(node);
//							nodesVec[count].posIndex = { particles[f].posIndex[0] + i , particles[f].posIndex[1] + j , particles[f].posIndex[2] + k };
//							nodesVec[count].position = nodesVec[count].posIndex.cast<double>() * param.dx;
//
//						}
//
//					};
//				};
//			};
//
//		}
//
//
//		std::vector<std::vector<int>> bottomLeftBloackVec;
//		std::map<std::string, std::vector<int>>::iterator it;
//		for (it = bottomLeftBloack.begin(); it != bottomLeftBloack.end(); it++)
//		{
//			bottomLeftBloackVec.push_back(it->second);
//		}
//		std::vector<std::vector<int>> topLeftBloackVec;
//		for (it = topLeftBloack.begin(); it != topLeftBloack.end(); it++)
//		{
//			topLeftBloackVec.push_back(it->second);
//		}
//		std::vector<std::vector<int>> topRightBloackVec;
//		for (it = topRightBloack.begin(); it != topRightBloack.end(); it++)
//		{
//			topRightBloackVec.push_back(it->second);
//		}
//		std::vector<std::vector<int>> bottomRightBloackVec;
//		for (it = bottomRightBloack.begin(); it != bottomRightBloack.end(); it++)
//		{
//			bottomRightBloackVec.push_back(it->second);
//		}
//
//
//#pragma omp parallel for num_threads(param.numOfThreads)
//		for (int n = 0; n < bottomLeftBloackVec.size(); n++)
//		{
//			for (int m = 0; m < bottomLeftBloackVec[n].size(); m++)
//			{
//				int f = bottomLeftBloackVec[n][m];
//				for (int i = 0; i < 3; i++)
//				{
//					for (int j = 0; j < 3; j++)
//					{
//						for (int k = 0; k < 3; k++)
//						{
//							std::string ID = calculateID_string(particles[f].posIndex[0] + i, particles[f].posIndex[1] + j, particles[f].posIndex[2] + k);
//							double weight = particles[f].weight(0, i) * particles[f].weight(1, j) * particles[f].weight(2, k);
//							Eigen::Vector3d deltaWeight = { particles[f].deltaWeight(0, i) * particles[f].weight(1, j) * particles[f].weight(2, k),   particles[f].weight(0, i) * particles[f].deltaWeight(1, j) * particles[f].weight(2, k),  particles[f].weight(0, i) * particles[f].weight(1, j) * particles[f].deltaWeight(2, k) };
//
//							int nodeIndex = gridMap[ID];
//
//							nodesVec[nodeIndex].supportParticles.push_back(f);
//							nodesVec[nodeIndex].supportParticlesWeight.push_back(weight);
//							nodesVec[nodeIndex].supportParticlesDeltaWeight.push_back(deltaWeight);
//
//							particles[f].supportNodes.push_back(nodeIndex);
//							particles[f].supportNodeWeight.push_back(weight);
//							particles[f].supportNodeDeltaWeight.push_back(deltaWeight);
//
//						};
//					};
//				};
//			}
//		}
//
//
//#pragma omp parallel for num_threads(param.numOfThreads)
//		for (int n = 0; n < topLeftBloackVec.size(); n++)
//		{
//			for (int m = 0; m < topLeftBloackVec[n].size(); m++)
//			{
//				int f = topLeftBloackVec[n][m];
//				for (int i = 0; i < 3; i++)
//				{
//					for (int j = 0; j < 3; j++)
//					{
//						for (int k = 0; k < 3; k++)
//						{
//							std::string ID = calculateID_string(particles[f].posIndex[0] + i, particles[f].posIndex[1] + j, particles[f].posIndex[2] + k);
//							double weight = particles[f].weight(0, i) * particles[f].weight(1, j) * particles[f].weight(2, k);
//							Eigen::Vector3d deltaWeight = { particles[f].deltaWeight(0, i) * particles[f].weight(1, j) * particles[f].weight(2, k),   particles[f].weight(0, i) * particles[f].deltaWeight(1, j) * particles[f].weight(2, k),  particles[f].weight(0, i) * particles[f].weight(1, j) * particles[f].deltaWeight(2, k) };
//
//							int nodeIndex = gridMap[ID];
//
//							nodesVec[nodeIndex].supportParticles.push_back(f);
//							nodesVec[nodeIndex].supportParticlesWeight.push_back(weight);
//							nodesVec[nodeIndex].supportParticlesDeltaWeight.push_back(deltaWeight);
//
//							particles[f].supportNodes.push_back(nodeIndex);
//							particles[f].supportNodeWeight.push_back(weight);
//							particles[f].supportNodeDeltaWeight.push_back(deltaWeight);
//
//						};
//					};
//				};
//			}
//		}
//
//
//#pragma omp parallel for num_threads(param.numOfThreads)
//		for (int n = 0; n < topRightBloackVec.size(); n++)
//		{
//			for (int m = 0; m < topRightBloackVec[n].size(); m++)
//			{
//				int f = topRightBloackVec[n][m];
//				for (int i = 0; i < 3; i++)
//				{
//					for (int j = 0; j < 3; j++)
//					{
//						for (int k = 0; k < 3; k++)
//						{
//							std::string ID = calculateID_string(particles[f].posIndex[0] + i, particles[f].posIndex[1] + j, particles[f].posIndex[2] + k);
//							double weight = particles[f].weight(0, i) * particles[f].weight(1, j) * particles[f].weight(2, k);
//							Eigen::Vector3d deltaWeight = { particles[f].deltaWeight(0, i) * particles[f].weight(1, j) * particles[f].weight(2, k),   particles[f].weight(0, i) * particles[f].deltaWeight(1, j) * particles[f].weight(2, k),  particles[f].weight(0, i) * particles[f].weight(1, j) * particles[f].deltaWeight(2, k) };
//
//							int nodeIndex = gridMap[ID];
//
//							Eigen::Vector3d pos = nodesVec[nodeIndex].position;
//							nodesVec[nodeIndex].supportParticles.push_back(f);
//							nodesVec[nodeIndex].supportParticlesWeight.push_back(weight);
//							nodesVec[nodeIndex].supportParticlesDeltaWeight.push_back(deltaWeight);
//
//							particles[f].supportNodes.push_back(nodeIndex);
//							particles[f].supportNodeWeight.push_back(weight);
//							particles[f].supportNodeDeltaWeight.push_back(deltaWeight);
//
//						};
//					};
//				};
//			}
//		}
//
//
//#pragma omp parallel for num_threads(param.numOfThreads)
//		for (int n = 0; n < bottomRightBloackVec.size(); n++)
//		{
//			for (int m = 0; m < bottomRightBloackVec[n].size(); m++)
//			{
//				int f = bottomRightBloackVec[n][m];
//				for (int i = 0; i < 3; i++)
//				{
//					for (int j = 0; j < 3; j++)
//					{
//						for (int k = 0; k < 3; k++)
//						{
//							std::string ID = calculateID_string(particles[f].posIndex[0] + i, particles[f].posIndex[1] + j, particles[f].posIndex[2] + k);
//							double weight = particles[f].weight(0, i) * particles[f].weight(1, j) * particles[f].weight(2, k);
//							Eigen::Vector3d deltaWeight = { particles[f].deltaWeight(0, i) * particles[f].weight(1, j) * particles[f].weight(2, k),   particles[f].weight(0, i) * particles[f].deltaWeight(1, j) * particles[f].weight(2, k),  particles[f].weight(0, i) * particles[f].weight(1, j) * particles[f].deltaWeight(2, k) };
//
//
//							int nodeIndex = gridMap[ID];
//
//							nodesVec[nodeIndex].supportParticles.push_back(f);
//							nodesVec[nodeIndex].supportParticlesWeight.push_back(weight);
//							nodesVec[nodeIndex].supportParticlesDeltaWeight.push_back(deltaWeight);
//
//							particles[f].supportNodes.push_back(nodeIndex);
//							particles[f].supportNodeWeight.push_back(weight);
//							particles[f].supportNodeDeltaWeight.push_back(deltaWeight);
//
//						};
//					};
//				};
//			}
//		}
//
//
//	}
//
//
//	// particle to grid transfer
//	void particle2Grid(std::vector<MPMParticle>& particles, MPMParamters& param, std::vector<Grid>& nodesVec)
//	{
//		// transfer particle's mass and momentum to grid nodes
//#pragma omp parallel for num_threads(param.numOfThreads)
//		for (int g = 0; g < nodesVec.size(); g++)
//		{
//			for (int p = 0; p < nodesVec[g].supportParticles.size(); p++)
//			{
//				int parPosInParticleVec = nodesVec[g].supportParticles[p];
//				double weight = nodesVec[g].supportParticlesWeight[p];
//
//				double mass = particles[parPosInParticleVec].mass;
//				Eigen::Vector3d CMultiPos = nodesVec[g].position - particles[parPosInParticleVec].position;
//				Eigen::Vector3d affineContribution = particles[parPosInParticleVec].affine * CMultiPos;
//
//				// transfer mass and momentum
//				nodesVec[g].mass += mass * weight;
//				//nodesVec[g].momentum += mass * weight * particles[parPosInParticleVec].velocity;
//
//				// MLS-MPM
//				nodesVec[g].momentum += mass * weight * (particles[parPosInParticleVec].velocity + affineContribution);
//			}
//		}
//
//	}
//
//
//	// update each particle's cauchy stress
//	void updateParInternalForce(std::vector<MPMParticle>& particles, MPMParamters& param)
//	{
//		// calculate each particle's internal cauchy stress
//#pragma omp parallel for num_threads(param.numOfThreads)
//		for (int f = 0; f < particles.size(); f++)
//		{
//			Eigen::Matrix3d F = particles[f].F;
//			double J = F.determinant();
//			Eigen::Matrix3d cauchyStressE = (param.mat_mpm.lambda * log(J) / J - param.mat_mpm.mu / J) * Eigen::Matrix3d::Identity() + param.mat_mpm.mu / J * F * F.transpose();
//
//
//			// compute eigenvalue and eigenvector
//			Eigen::EigenSolver<Eigen::MatrixXd> es(cauchyStressE);
//			Eigen::Vector3d eigenValues = { es.eigenvalues()[0].real() ,  es.eigenvalues()[1].real() ,  es.eigenvalues()[2].real() };
//			Eigen::Matrix3d eigenVectors;
//			eigenVectors << es.eigenvectors().col(0)[0].real(), es.eigenvectors().col(1)[0].real(), es.eigenvectors().col(2)[0].real(),
//				es.eigenvectors().col(0)[1].real(), es.eigenvectors().col(1)[1].real(), es.eigenvectors().col(2)[1].real(),
//				es.eigenvectors().col(0)[2].real(), es.eigenvectors().col(1)[2].real(), es.eigenvectors().col(2)[2].real();
//			double maxEigenValue = std::max(std::max(eigenValues[0], eigenValues[1]), eigenValues[2]);
//
//
//			if (maxEigenValue > param.mat_mpm.thetaf)
//			{
//				double tempDp = (1 + param.mat_mpm.Hs) * (1 - param.mat_mpm.thetaf / maxEigenValue);
//				if (maxEigenValue > (1 + 1 / param.mat_mpm.Hs) * param.mat_mpm.thetaf)
//				{
//					particles[f].dp = 1.0;
//				}
//				else
//				{
//					if (tempDp > particles[f].dp)
//					{
//						particles[f].dp = tempDp;
//					};
//				};
//			};
//
//			Eigen::Vector3d sigmaPlus = { 0 , 0 , 0 };
//			for (int i = 0; i < 3; i++)
//			{
//				if (eigenValues[i] > 0)
//				{
//					if (particles[f].dp >= param.damageThreshold)
//					{
//						sigmaPlus[i] = 0;
//					}
//					else
//					{
//						sigmaPlus[i] = (1 - particles[f].dp) * eigenValues[i];
//					};
//				}
//				else
//				{
//					sigmaPlus[i] = eigenValues[i];
//				};
//			};
//
//			Eigen::Matrix3d sigma = Eigen::Matrix3d::Zero();
//			for (int i = 0; i < 3; i++)
//			{
//				sigma = sigma + sigmaPlus[i] * eigenVectors.col(i) * (eigenVectors.col(i).transpose());
//			};
//
//			cauchyStressE = sigma;
//
//			particles[f].cauchyStress = cauchyStressE;
//		}
//
//		
//
//
//
//	}
//
//
//	// calculate the grid node's internal force induced by particles
//	void calculateNodeForce(std::vector<MPMParticle>& particles, MPMParamters& param, std::vector<Grid>& nodesVec)
//	{
//		// transfer particle's interal force to grid nodes
//#pragma omp parallel for num_threads(param.numOfThreads)
//		for (int g = 0; g < nodesVec.size(); g++)
//		{
//			for (int p = 0; p < nodesVec[g].supportParticles.size(); p++)
//			{
//				int parPosInParticleVec = nodesVec[g].supportParticles[p];
//				double weight = nodesVec[g].supportParticlesWeight[p];
//				Eigen::Vector3d deltaWeight = nodesVec[g].supportParticlesDeltaWeight[p];
//				{
//					// // APIC-MPM implementation
//					// nodesVec[g].force += -weight * param.dt * particles[parPosInParticleVec].volume * (particles[parPosInParticleVec].F).determinant() * particles[parPosInParticleVec].cauchyStress
//					// 	* nodesVec[g].supportParticlesDeltaWeight[p];
//
//
//					//nodesVec[g].force += - param.dt * particles[parPosInParticleVec].volume * particles[parPosInParticleVec].cauchyStress * deltaWeight;
//
//
//					//MLS-MPM implementation
//					nodesVec[g].force += - param.dt * weight / (param.dx * param.dx / 4.0) *  particles[parPosInParticleVec].volume * (particles[parPosInParticleVec].F).determinant() * (particles[parPosInParticleVec].cauchyStress * (nodesVec[g].position - particles[parPosInParticleVec].position)).transpose();
//				}
//			}
//		}
//
//	}
//
//
//	// grid momentum update
//	void gridUpdate(std::vector<Grid>& nodesVec, MPMParamters& param)
//	{
//		// calculate nodes' force, solve the momentum equation and update node's velocity
//		// add gravity and boundary condition
//#pragma omp parallel for num_threads(param.numOfThreads)
//		for (int g = 0; g < nodesVec.size(); g++)
//		{
//			if (nodesVec[g].mass > 0)
//			{
//				Eigen::Vector3d velocity = nodesVec[g].momentum / nodesVec[g].mass; // node velcity of timestep n
//				Eigen::Vector3d acceleration = nodesVec[g].force / nodesVec[g].mass + param.gravity;
//				nodesVec[g].velocity = velocity + param.dt * acceleration; // node velocity of Lagrangian phase
//				//nodesVec[g].positionUpdated = nodesVec[g].position + nodesVec[g].velocity * param.dt;
//			};
//		};
//
//
//
//
//
//	}
//
//
//	// grid to particle transfer
//	void grid2Particle(std::vector<MPMParticle>& particles, MPMParamters& param, std::vector<Grid>& nodesVec)
//	{
//
//#pragma omp parallel for num_threads(param.numOfThreads)
//		for (int f = 0; f < particles.size(); f++)
//		{
//
//			particles[f].velocity = Eigen::Vector3d::Zero(); // MPMParticle-in-cell method
//			Eigen::Matrix3d affine = Eigen::Matrix3d::Zero();
//			for (int g = 0; g < particles[f].supportNodes.size(); g++)
//			{
//				int nodePosInNodesVec = particles[f].supportNodes[g];
//				Eigen::Vector3d CMultiPos = nodesVec[nodePosInNodesVec].position - particles[f].position;
//
//				double weight = particles[f].supportNodeWeight[g];
//
//				particles[f].velocity += weight * nodesVec[nodePosInNodesVec].velocity;
//				affine += weight / (param.dx * param.dx / 4.0) * nodesVec[nodePosInNodesVec].velocity * CMultiPos.transpose();
//			}
//			particles[f].F = (Eigen::Matrix3d::Identity() + param.dt * affine) * particles[f].F;
//			particles[f].affine = affine;
//			particles[f].position += param.dt * particles[f].velocity;
//
//
//			//particles[f].F = (Eigen::Matrix3d::Identity() + param.dt * deltaVel) * particles[f].F;
//
//
//		};
//
//	}
//
//
//	// apply point force
//	void applyPointForce(MPMParamters& param, std::vector<Grid>& nodesVec, std::map<std::string, int>& gridMap, std::vector<Vector6d>& contactForce)
//	{
//		for (int i = 0; i < contactForce.size(); i++)
//		{		
//			Eigen::Vector3d forcePosition = contactForce[i].segment(0,3);
//			Eigen::Vector3d forceMagnitude = contactForce[i].segment(3, 3);
//
//			extractCrackSurface::weightAndDreri WD = extractCrackSurface::calWeight(param.dx, forcePosition);
//			Eigen::Vector3i ppIndex = WD.ppIndex;
//			Eigen::MatrixXd weightForce = WD.weight;
//
//			for (int i = 0; i < 3; i++) {
//				for (int j = 0; j < 3; j++) {
//					for (int k = 0; k < 3; k++) {
//						std::string ID = calculateID_string(ppIndex[0] + i, ppIndex[1] + j, ppIndex[2] + k);
//						double weight = weightForce(0, i) * weightForce(1, j) * weightForce(2, k);
//						if (gridMap.find(ID) != gridMap.end())
//						{
//							int nodeIndex = gridMap[ID];
//							if (weight != 0)
//							{
//								int eid = gridMap[ID];
//								nodesVec[nodeIndex].force += weight * forceMagnitude;
//							};
//						}
//						else
//						{
//							std::cout << "Warning: the conatct force is not applied to mpm particles!" << std::endl;
//						}
//
//					};
//				};
//			};
//		}
//		
//	}
//
//
//	void advanceStep(std::vector<MPMParticle>& particles, MPMParamters& param, std::vector<Vector6d>& contactForce, int timestep) // prticle vector, timestep
//	{
//		// // initialize background grid nodes
//		std::vector<Grid> nodesVec;
//		std::map<std::string, int> gridMap;
//
//
//		// calculate the reationreship between particles and background grid 
//		calWeightsAndNodes(particles, param, nodesVec, gridMap);
//		// transfer information from particle to grdi nodes
//		particle2Grid(particles, param, nodesVec);
//		// update each material particle's cauchy stress
//		updateParInternalForce(particles, param);
//		// apply point force
//		//applyPointForce(param, nodesVec, gridMap, contactForce);
//		// calculate the grid node's internal force induced by particles
//		calculateNodeForce(particles, param, nodesVec);
//		// grid nodes momentum update
//		gridUpdate(nodesVec, param);
//		// transfer information back form grid to particles
//		grid2Particle(particles, param, nodesVec);
//
//		{
//			std::ofstream outfile9("./output/velocity_" + std::to_string(timestep) + ".obj", std::ios::trunc);
//			for (int k = 0; k < nodesVec.size(); k++)
//			{
//				Eigen::Vector3d scale = nodesVec[k].position;
//				outfile9 << std::scientific << std::setprecision(8) << scale[0] << " " << scale[1] << " " << scale[2] << " " << nodesVec[k].velocity.norm() << std::endl;
//			}
//			outfile9.close();
//		}
//
//		{
//			std::ofstream outfile9("./output/force_" + std::to_string(timestep) + ".obj", std::ios::trunc);
//			for (int k = 0; k < nodesVec.size(); k++)
//			{
//				Eigen::Vector3d scale = nodesVec[k].position;
//				outfile9 << std::scientific << std::setprecision(8) << scale[0] << " " << scale[1] << " " << scale[2] << " " << nodesVec[k].force.norm() << std::endl;
//			}
//			outfile9.close();
//		}
//
//
//		{
//			std::ofstream outfile9("./output/velocity_par_" + std::to_string(timestep) + ".obj", std::ios::trunc);
//			for (int k = 0; k < particles.size(); k++)
//			{
//				Eigen::Vector3d scale = particles[k].position;
//				outfile9 << std::scientific << std::setprecision(8) << scale[0] << " " << scale[1] << " " << scale[2] << " " << particles[k].velocity.norm() << std::endl;
//			}
//			outfile9.close();
//		}
//
//
//
//
//	};
//	
//
//	// extract crack surface
//	std::tuple<bool, objMeshFormat, objMeshFormat, std::vector<objMeshFormat>> tryToExtractCracks(std::vector<MPMParticle>& particles, MPMParamters& param, int timestep)
//	{
//		std::vector<extractCrackSurface::Particle> particlesRaw;
//		extractCrackSurface::CRACKParamters paramCrack;
//		for (int i = 0; i < particles.size(); i++)
//		{
//			particlesRaw.push_back(extractCrackSurface::Particle(particles[i].position, particles[i].velocity, particles[i].mass, particles[i].color, particles[i].dp));
//		}
//		paramCrack.dx = param.dx;
//
//
//		std::tuple<bool, objMeshFormat, objMeshFormat, std::vector<objMeshFormat>> cracks = extractCrackSurface::extractCrackSurf(&particlesRaw, paramCrack);
//
//		bool findCrackSurface = std::get<0>(cracks);
//
//		objMeshFormat crackSurfacePartialCut;
//		crackSurfacePartialCut.vertices = std::get<1>(cracks).vertices;
//		crackSurfacePartialCut.faces = std::get<1>(cracks).faces;
//
//		objMeshFormat crackSurfaceFullCut;
//		crackSurfaceFullCut.vertices = std::get<2>(cracks).vertices;
//		crackSurfaceFullCut.faces = std::get<2>(cracks).faces;
//
//		std::vector<objMeshFormat> allFragmentsObj;
//		for (int k = 0; k < std::get<3>(cracks).size(); k++)
//		{
//			objMeshFormat frag;
//			frag.vertices = std::get<3>(cracks)[k].vertices;
//			frag.faces = std::get<3>(cracks)[k].faces;
//			allFragmentsObj.push_back(frag);
//		}
//
//		std::tuple<bool, objMeshFormat, objMeshFormat, std::vector<objMeshFormat>> resultReturn(findCrackSurface, crackSurfacePartialCut, crackSurfaceFullCut, allFragmentsObj);
//		return resultReturn;
//	}
//
//
//	std::tuple<bool, objMeshFormat, objMeshFormat, std::vector<objMeshFormat>> crackSimulation(
//		const std::vector<Eigen::Vector3d>& points, 
//		const double& volume, 
//		const Material& mat_mpm,
//		MPMParamters& param,
//		std::vector<Vector6d>& contactForce, 
//		int num_timestep)
//	{
//		
//		std::vector<MPMParticle> particles;
//		//initialize_mpm_particles(particles, points, volume, mat_mpm);
//
//		for (int i = 0; i < points.size(); i++)
//		{
//			MPMParticle pt;
//			pt.volume = volume;
//			pt.mass = volume * mat_mpm.density;
//			Eigen::Vector3d dev = {0,-4.2,0};
//			pt.position = points[i] + dev;
//			pt.velocity = {0, 10, 0};
//			particles.push_back(pt);
//		}
//
//
//
//		bool generateFragments = false;
//		bool simulationFailed = false;
//		std::tuple<bool, objMeshFormat, objMeshFormat, std::vector<objMeshFormat>> crackSurfs;
//		std::cout << "MPM timestep = ";
//		for (int timestep = 0; timestep <= 2000 && generateFragments == false && !simulationFailed; timestep++)
//		{
//			std::cout << timestep<<" ";
//			advanceStep(particles, param, contactForce, timestep);
//
//
//			{
//				std::ofstream outfile9("./output/MPM_particles_" + std::to_string(timestep) + ".obj", std::ios::trunc);
//				for (int k = 0; k < particles.size(); k++)
//				{
//					Eigen::Vector3d scale = particles[k].position;
//					outfile9 << std::scientific << std::setprecision(8) << "v " << scale[0] << " " << scale[1] << " " << scale[2] << std::endl;
//				}
//				outfile9.close();
//			}
//
//			// The simulation fails due to unexpected reasons. Test
//			for (int k = 0; k < particles.size(); k++)
//			{
//				Eigen::Vector3d scale = particles[k].position;
//				if (scale.hasNaN() || scale.array().isInf().any())
//				{
//					std::cout << "Conataining inf or NaN, exit!" << std::endl << std::endl << std::endl;
//					simulationFailed = true;
//					break;
//				}
//			}
//
//			if (timestep % 20 == 0 && !simulationFailed)
//			{
//				crackSurfs = tryToExtractCracks(particles, param, timestep);
//				if (std::get<0>(crackSurfs) == true && std::get<3>(crackSurfs).size() > 1)
//				{
//					generateFragments = true;
//				}
//			}
//
//		}
//		std::cout << std::endl;
//
//		return crackSurfs;
//
//	}
//
//
//	void initialize_mpm_particles(std::vector<MPMParticle>& particles, const std::vector<Eigen::Vector3d>& points, const double& volume, const Material& mat_mpm)
//	{
//		particles.clear();
//		for (int i = 0; i < points.size(); i++)
//		{
//			MPMParticle pt;
//			pt.volume = volume;
//			pt.mass = volume * mat_mpm.density;
//			pt.position = points[i];
//			particles.push_back(pt);
//		}
//		
//	}
//
//}
//
//
//
//
//




#include "mpmSimulator.h"


namespace mpmSimulator
{
	
	std::string calculateID_string(int x, int y, int z) // string id of the particle
	{
		return std::to_string(x) + "#" + std::to_string(y) + "#" + std::to_string(z);
	};


	// calculate particles' weights and find neighbouring nodes
	void calWeightsAndNodes(std::vector<MPMParticle>& particles, MPMParamters& param, std::vector<Grid>& nodesVec, std::map<std::string, int>& gridMap)
	{
		// store the key and value of gridMap: active grid

		Eigen::Vector3i minCellIndex = { 100000000, 100000000, 100000000 };
		Eigen::Vector3i maxCellIndex = { -100000000, -100000000, -100000000 };
#pragma omp parallel for num_threads(param.numOfThreads)
		for (int f = 0; f < particles.size(); f++)
		{
			extractCrackSurface::weightAndDreri  WD = extractCrackSurface::calWeight(param.dx, particles[f].position);

			particles[f].posIndex = WD.ppIndex;
			particles[f].weight = WD.weight;
			particles[f].deltaWeight = WD.deltaWeight;

			particles[f].supportNodes.clear();
			particles[f].supportNodeWeight.clear();
			particles[f].supportNodeDeltaWeight.clear();


			minCellIndex[0] = std::min(WD.ppIndex[0], minCellIndex[0]);
			maxCellIndex[0] = std::max(WD.ppIndex[0], maxCellIndex[0]);

			minCellIndex[1] = std::min(WD.ppIndex[1], minCellIndex[1]);
			maxCellIndex[1] = std::max(WD.ppIndex[1], maxCellIndex[1]);

			minCellIndex[2] = std::min(WD.ppIndex[2], minCellIndex[2]);
			maxCellIndex[2] = std::max(WD.ppIndex[2], maxCellIndex[2]);
		};

		// find parameters that used to devide the background grid into disconnected patches
		int blockSize = 6; // the size of each block. The minimum size for quadratic kernel is 5!!!!!!!!!!!!!!!!
		std::map<std::string, std::vector<int>> bottomLeftBloack;
		std::map<std::string, std::vector<int>> bottomRightBloack;
		std::map<std::string, std::vector<int>> topLeftBloack;
		std::map<std::string, std::vector<int>> topRightBloack;
		Eigen::Vector3i span = { maxCellIndex[0] - minCellIndex[0], maxCellIndex[1] - minCellIndex[1], maxCellIndex[2] - minCellIndex[2] };
		int minSpan = std::min(span[0], std::min(span[1], span[2]));
		if (minSpan == span[0])
		{
			for (int f = 0; f < particles.size(); f++)
			{
				int remainZ = (particles[f].posIndex[2] - minCellIndex[2]) % (2 * blockSize);
				int remainY = (particles[f].posIndex[1] - minCellIndex[1]) % (2 * blockSize);

				int numZ = (particles[f].posIndex[2] - minCellIndex[2]) / (2 * blockSize);
				int numY = (particles[f].posIndex[1] - minCellIndex[1]) / (2 * blockSize);
				std::string blockID = std::to_string(numZ) + "#" + std::to_string(numY);
				if (remainZ >= 0 && remainZ < blockSize) // left 
				{
					if (remainY >= 0 && remainY < blockSize) // bottom left
					{
						bottomLeftBloack[blockID].push_back(f);
					}
					else // top left
					{
						topLeftBloack[blockID].push_back(f);
					}
				}
				else // right
				{
					if (remainY >= 0 && remainY < blockSize) // top right
					{
						bottomRightBloack[blockID].push_back(f);
					}
					else // bottm right
					{
						topRightBloack[blockID].push_back(f);
					}
				}
			}
		}
		else if (minSpan == span[1])
		{
			for (int f = 0; f < particles.size(); f++)
			{
				int remainX = (particles[f].posIndex[0] - minCellIndex[0]) % (2 * blockSize);
				int remainZ = (particles[f].posIndex[2] - minCellIndex[2]) % (2 * blockSize);

				int numX = (particles[f].posIndex[0] - minCellIndex[0]) / (2 * blockSize);
				int numZ = (particles[f].posIndex[2] - minCellIndex[2]) / (2 * blockSize);
				std::string blockID = std::to_string(numX) + "#" + std::to_string(numZ);
				if (remainX >= 0 && remainX < blockSize) // left 
				{
					if (remainZ >= 0 && remainZ < blockSize) // bottom left
					{
						bottomLeftBloack[blockID].push_back(f);
					}
					else // top left
					{
						topLeftBloack[blockID].push_back(f);
					}
				}
				else // right
				{
					if (remainZ >= 0 && remainZ < blockSize) // top right
					{
						bottomRightBloack[blockID].push_back(f);
					}
					else // bottm right
					{
						topRightBloack[blockID].push_back(f);
					}
				}
			}
		}
		else
		{
			for (int f = 0; f < particles.size(); f++)
			{
				int remainX = (particles[f].posIndex[0] - minCellIndex[0]) % (2 * blockSize);
				int remainY = (particles[f].posIndex[1] - minCellIndex[1]) % (2 * blockSize);

				int numX = (particles[f].posIndex[0] - minCellIndex[0]) / (2 * blockSize);
				int numY = (particles[f].posIndex[1] - minCellIndex[1]) / (2 * blockSize);
				std::string blockID = std::to_string(numX) + "#" + std::to_string(numY);
				if (remainX >= 0 && remainX < blockSize) // left 
				{
					if (remainY >= 0 && remainY < blockSize) // bottom left
					{
						bottomLeftBloack[blockID].push_back(f);
					}
					else // top left
					{
						topLeftBloack[blockID].push_back(f);
					}
				}
				else // right
				{
					if (remainY >= 0 && remainY < blockSize) // bottom right
					{
						bottomRightBloack[blockID].push_back(f);
					}
					else // top right
					{
						topRightBloack[blockID].push_back(f);
					}
				}
			}



		}


		int count = nodesVec.size() - 1; // count the number of active grid node
		for (int f = 0; f < particles.size(); f++)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						std::string ID = calculateID_string(particles[f].posIndex[0] + i, particles[f].posIndex[1] + j, particles[f].posIndex[2] + k);

						if (gridMap.find(ID) == gridMap.end())
						{
							count += 1;
							gridMap[ID] = count;

							Grid node;
							nodesVec.push_back(node);
							nodesVec[count].posIndex = { particles[f].posIndex[0] + i , particles[f].posIndex[1] + j , particles[f].posIndex[2] + k };
							nodesVec[count].position = nodesVec[count].posIndex.cast<double>() * param.dx;

						}

					};
				};
			};

		}


		std::vector<std::vector<int>> bottomLeftBloackVec;
		std::map<std::string, std::vector<int>>::iterator it;
		for (it = bottomLeftBloack.begin(); it != bottomLeftBloack.end(); it++)
		{
			bottomLeftBloackVec.push_back(it->second);
		}
		std::vector<std::vector<int>> topLeftBloackVec;
		for (it = topLeftBloack.begin(); it != topLeftBloack.end(); it++)
		{
			topLeftBloackVec.push_back(it->second);
		}
		std::vector<std::vector<int>> topRightBloackVec;
		for (it = topRightBloack.begin(); it != topRightBloack.end(); it++)
		{
			topRightBloackVec.push_back(it->second);
		}
		std::vector<std::vector<int>> bottomRightBloackVec;
		for (it = bottomRightBloack.begin(); it != bottomRightBloack.end(); it++)
		{
			bottomRightBloackVec.push_back(it->second);
		}


#pragma omp parallel for num_threads(param.numOfThreads)
		for (int n = 0; n < bottomLeftBloackVec.size(); n++)
		{
			for (int m = 0; m < bottomLeftBloackVec[n].size(); m++)
			{
				int f = bottomLeftBloackVec[n][m];
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						for (int k = 0; k < 3; k++)
						{
							std::string ID = calculateID_string(particles[f].posIndex[0] + i, particles[f].posIndex[1] + j, particles[f].posIndex[2] + k);
							double weight = particles[f].weight(0, i) * particles[f].weight(1, j) * particles[f].weight(2, k);
							Eigen::Vector3d deltaWeight = { particles[f].deltaWeight(0, i) * particles[f].weight(1, j) * particles[f].weight(2, k),   particles[f].weight(0, i) * particles[f].deltaWeight(1, j) * particles[f].weight(2, k),  particles[f].weight(0, i) * particles[f].weight(1, j) * particles[f].deltaWeight(2, k) };

							int nodeIndex = gridMap[ID];

							nodesVec[nodeIndex].supportParticles.push_back(f);
							nodesVec[nodeIndex].supportParticlesWeight.push_back(weight);
							nodesVec[nodeIndex].supportParticlesDeltaWeight.push_back(deltaWeight);

							particles[f].supportNodes.push_back(nodeIndex);
							particles[f].supportNodeWeight.push_back(weight);
							particles[f].supportNodeDeltaWeight.push_back(deltaWeight);

						};
					};
				};
			}
		}


#pragma omp parallel for num_threads(param.numOfThreads)
		for (int n = 0; n < topLeftBloackVec.size(); n++)
		{
			for (int m = 0; m < topLeftBloackVec[n].size(); m++)
			{
				int f = topLeftBloackVec[n][m];
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						for (int k = 0; k < 3; k++)
						{
							std::string ID = calculateID_string(particles[f].posIndex[0] + i, particles[f].posIndex[1] + j, particles[f].posIndex[2] + k);
							double weight = particles[f].weight(0, i) * particles[f].weight(1, j) * particles[f].weight(2, k);
							Eigen::Vector3d deltaWeight = { particles[f].deltaWeight(0, i) * particles[f].weight(1, j) * particles[f].weight(2, k),   particles[f].weight(0, i) * particles[f].deltaWeight(1, j) * particles[f].weight(2, k),  particles[f].weight(0, i) * particles[f].weight(1, j) * particles[f].deltaWeight(2, k) };

							int nodeIndex = gridMap[ID];

							nodesVec[nodeIndex].supportParticles.push_back(f);
							nodesVec[nodeIndex].supportParticlesWeight.push_back(weight);
							nodesVec[nodeIndex].supportParticlesDeltaWeight.push_back(deltaWeight);

							particles[f].supportNodes.push_back(nodeIndex);
							particles[f].supportNodeWeight.push_back(weight);
							particles[f].supportNodeDeltaWeight.push_back(deltaWeight);

						};
					};
				};
			}
		}


#pragma omp parallel for num_threads(param.numOfThreads)
		for (int n = 0; n < topRightBloackVec.size(); n++)
		{
			for (int m = 0; m < topRightBloackVec[n].size(); m++)
			{
				int f = topRightBloackVec[n][m];
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						for (int k = 0; k < 3; k++)
						{
							std::string ID = calculateID_string(particles[f].posIndex[0] + i, particles[f].posIndex[1] + j, particles[f].posIndex[2] + k);
							double weight = particles[f].weight(0, i) * particles[f].weight(1, j) * particles[f].weight(2, k);
							Eigen::Vector3d deltaWeight = { particles[f].deltaWeight(0, i) * particles[f].weight(1, j) * particles[f].weight(2, k),   particles[f].weight(0, i) * particles[f].deltaWeight(1, j) * particles[f].weight(2, k),  particles[f].weight(0, i) * particles[f].weight(1, j) * particles[f].deltaWeight(2, k) };

							int nodeIndex = gridMap[ID];

							Eigen::Vector3d pos = nodesVec[nodeIndex].position;
							nodesVec[nodeIndex].supportParticles.push_back(f);
							nodesVec[nodeIndex].supportParticlesWeight.push_back(weight);
							nodesVec[nodeIndex].supportParticlesDeltaWeight.push_back(deltaWeight);

							particles[f].supportNodes.push_back(nodeIndex);
							particles[f].supportNodeWeight.push_back(weight);
							particles[f].supportNodeDeltaWeight.push_back(deltaWeight);

						};
					};
				};
			}
		}


#pragma omp parallel for num_threads(param.numOfThreads)
		for (int n = 0; n < bottomRightBloackVec.size(); n++)
		{
			for (int m = 0; m < bottomRightBloackVec[n].size(); m++)
			{
				int f = bottomRightBloackVec[n][m];
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						for (int k = 0; k < 3; k++)
						{
							std::string ID = calculateID_string(particles[f].posIndex[0] + i, particles[f].posIndex[1] + j, particles[f].posIndex[2] + k);
							double weight = particles[f].weight(0, i) * particles[f].weight(1, j) * particles[f].weight(2, k);
							Eigen::Vector3d deltaWeight = { particles[f].deltaWeight(0, i) * particles[f].weight(1, j) * particles[f].weight(2, k),   particles[f].weight(0, i) * particles[f].deltaWeight(1, j) * particles[f].weight(2, k),  particles[f].weight(0, i) * particles[f].weight(1, j) * particles[f].deltaWeight(2, k) };


							int nodeIndex = gridMap[ID];

							nodesVec[nodeIndex].supportParticles.push_back(f);
							nodesVec[nodeIndex].supportParticlesWeight.push_back(weight);
							nodesVec[nodeIndex].supportParticlesDeltaWeight.push_back(deltaWeight);

							particles[f].supportNodes.push_back(nodeIndex);
							particles[f].supportNodeWeight.push_back(weight);
							particles[f].supportNodeDeltaWeight.push_back(deltaWeight);

						};
					};
				};
			}
		}


	}


	// particle to grid transfer
	void particle2Grid(std::vector<MPMParticle>& particles, MPMParamters& param, std::vector<Grid>& nodesVec)
	{
		// transfer particle's mass and momentum to grid nodes
#pragma omp parallel for num_threads(param.numOfThreads)
		for (int g = 0; g < nodesVec.size(); g++)
		{
			for (int p = 0; p < nodesVec[g].supportParticles.size(); p++)
			{
				int parPosInParticleVec = nodesVec[g].supportParticles[p];
				double weight = nodesVec[g].supportParticlesWeight[p];

				double mass = particles[parPosInParticleVec].mass;
				Eigen::Vector3d CMultiPos = nodesVec[g].position - particles[parPosInParticleVec].position;
				Eigen::Vector3d affineContribution = particles[parPosInParticleVec].affine * CMultiPos;

				// transfer mass and momentum
				nodesVec[g].mass += mass * weight;
				nodesVec[g].momentum += mass * weight * (particles[parPosInParticleVec].velocity + affineContribution);
			}
		}

	}


	// update each particle's cauchy stress
	void updateParInternalForce(std::vector<MPMParticle>& particles, MPMParamters& param)
	{
		// calculate each particle's internal cauchy stress
#pragma omp parallel for num_threads(param.numOfThreads)
		for (int f = 0; f < particles.size(); f++)
		{
			Eigen::Matrix3d F = particles[f].F;
			double J = F.determinant();
			Eigen::Matrix3d cauchyStressE = (param.mat_mpm.lambda * log(J) / J - param.mat_mpm.mu / J) * Eigen::Matrix3d::Identity() + param.mat_mpm.mu / J * F * F.transpose();


			{
				// compute eigenvalue and eigenvector
				Eigen::EigenSolver<Eigen::MatrixXd> es(cauchyStressE);
				Eigen::Vector3d eigenValues = { es.eigenvalues()[0].real() ,  es.eigenvalues()[1].real() ,  es.eigenvalues()[2].real() };
				Eigen::Matrix3d eigenVectors;
				eigenVectors << es.eigenvectors().col(0)[0].real(), es.eigenvectors().col(1)[0].real(), es.eigenvectors().col(2)[0].real(),
					es.eigenvectors().col(0)[1].real(), es.eigenvectors().col(1)[1].real(), es.eigenvectors().col(2)[1].real(),
					es.eigenvectors().col(0)[2].real(), es.eigenvectors().col(1)[2].real(), es.eigenvectors().col(2)[2].real();
				double maxEigenValue = std::max(std::max(eigenValues[0], eigenValues[1]), eigenValues[2]);


				if (maxEigenValue > param.mat_mpm.thetaf)
				{
					double tempDp = (1 + param.mat_mpm.Hs) * (1 - param.mat_mpm.thetaf / maxEigenValue);
					if (maxEigenValue > (1 + 1 / param.mat_mpm.Hs) * param.mat_mpm.thetaf)
					{
						particles[f].dp = 1.0;
					}
					else
					{
						if (tempDp > particles[f].dp)
						{
							particles[f].dp = tempDp;
						};
					};
				};

				Eigen::Vector3d sigmaPlus = { 0 , 0 , 0 };
				for (int i = 0; i < 3; i++)
				{
					if (eigenValues[i] > 0)
					{
						if (particles[f].dp >= param.damageThreshold)
						{
							sigmaPlus[i] = 0;
						}
						else
						{
							sigmaPlus[i] = (1 - particles[f].dp) * eigenValues[i];
						};
					}
					else
					{
						sigmaPlus[i] = eigenValues[i];
					};
				};

				Eigen::Matrix3d sigma = Eigen::Matrix3d::Zero();
				for (int i = 0; i < 3; i++)
				{
					sigma = sigma + sigmaPlus[i] * eigenVectors.col(i) * (eigenVectors.col(i).transpose());
				};

				cauchyStressE = sigma;
			}

			particles[f].cauchyStress = cauchyStressE;



		}


	}


	// calculate the grid node's internal force induced by particles
	void calculateNodeForce(std::vector<MPMParticle>& particles, MPMParamters& param, std::vector<Grid>& nodesVec)
	{

		// transfer particle's interal force to grid nodes
#pragma omp parallel for num_threads(param.numOfThreads)
		for (int g = 0; g < nodesVec.size(); g++)
		{
			for (int p = 0; p < nodesVec[g].supportParticles.size(); p++)
			{
				int parPosInParticleVec = nodesVec[g].supportParticles[p];
				double weight = nodesVec[g].supportParticlesWeight[p];
				{
					// // APIC-MPM implementation
					// nodesVec[g].force += -weight * param.dt * particles[parPosInParticleVec].volume * (particles[parPosInParticleVec].F).determinant() * particles[parPosInParticleVec].cauchyStress
					// 	* nodesVec[g].supportParticlesDeltaWeight[p];




					Eigen::Vector3d CMultiPos = nodesVec[g].position - particles[parPosInParticleVec].position;
					nodesVec[g].force += -weight / (param.dx * param.dx / 4.0) * param.dt * particles[parPosInParticleVec].volume * (particles[parPosInParticleVec].F).determinant() * (particles[parPosInParticleVec].cauchyStress * CMultiPos).transpose();



					//MLS-MPM implementation
					// nodesVec[g].force += -weight / param.DP * param.dt * particles[parPosInParticleVec].volume * (particles[parPosInParticleVec].F).determinant() * (particles[parPosInParticleVec].cauchyStress * (nodesVec[g].position - particles[parPosInParticleVec].position)).transpose();
				}
			}
		}

	}


	// grid momentum update
	void gridUpdate(std::vector<Grid>& nodesVec, MPMParamters& param)
	{
		// calculate nodes' force, solve the momentum equation and update node's velocity
		// add gravity and boundary condition
#pragma omp parallel for num_threads(param.numOfThreads)
		for (int g = 0; g < nodesVec.size(); g++)
		{
			if (nodesVec[g].mass > 0)
			{

				Eigen::Vector3d velocity = nodesVec[g].momentum / nodesVec[g].mass; // node velcity of timestep n
				Eigen::Vector3d acceleration = nodesVec[g].force / nodesVec[g].mass + param.gravity;
				nodesVec[g].velocity = velocity + param.dt * acceleration; // node velocity of Lagrangian phase
				//nodesVec[g].positionUpdated = nodesVec[g].position + nodesVec[g].velocity * param.dt;
			};
		};

	}


	// grid to particle transfer
	void grid2Particle(std::vector<MPMParticle>& particles, MPMParamters& param, std::vector<Grid>& nodesVec)
	{

#pragma omp parallel for num_threads(param.numOfThreads)
		for (int f = 0; f < particles.size(); f++)
		{
			particles[f].velocity = Eigen::Vector3d::Zero(); // MPMParticle-in-cell method
			Eigen::Matrix3d affine = Eigen::Matrix3d::Zero();
			for (int g = 0; g < particles[f].supportNodes.size(); g++)
			{
				//std::cout << "g = " << g << "; weight = " << particles[f].supportNodeWeight[g] << std::endl;
				int nodePosInNodesVec = particles[f].supportNodes[g];
				Eigen::Vector3d CMultiPos = nodesVec[nodePosInNodesVec].position - particles[f].position;

				double weight = particles[f].supportNodeWeight[g];
				particles[f].velocity += weight * nodesVec[nodePosInNodesVec].velocity;
				affine += weight / (param.dx * param.dx / 4.0) * nodesVec[nodePosInNodesVec].velocity * CMultiPos.transpose();
			}
			particles[f].F = (Eigen::Matrix3d::Identity() + param.dt * affine) * particles[f].F;
			particles[f].affine = affine;



			particles[f].position += param.dt * particles[f].velocity;


		};



	}


	// apply point force
	void applyPointForce(MPMParamters& param, std::vector<Grid>& nodesVec, std::map<std::string, int>& gridMap, std::vector<Vector6d>& contactForce)
	{
		for (int i = 0; i < contactForce.size(); i++)
		{		
			Eigen::Vector3d forcePosition = contactForce[i].segment(0,3);
			Eigen::Vector3d forceMagnitude = contactForce[i].segment(3, 3);

			extractCrackSurface::weightAndDreri WD = extractCrackSurface::calWeight(param.dx, forcePosition);
			Eigen::Vector3i ppIndex = WD.ppIndex;
			Eigen::MatrixXd weightForce = WD.weight;

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					for (int k = 0; k < 3; k++) {
						std::string ID = calculateID_string(ppIndex[0] + i, ppIndex[1] + j, ppIndex[2] + k);
						double weight = weightForce(0, i) * weightForce(1, j) * weightForce(2, k);
						if (gridMap.find(ID) != gridMap.end())
						{
							int nodeIndex = gridMap[ID];
							if (weight != 0)
							{
								int eid = gridMap[ID];
								nodesVec[nodeIndex].force += weight * forceMagnitude;
							};
						}
						else
						{
							std::cout << "Warning: the conatct force is not applied to mpm particles!" << std::endl;
						}

					};
				};
			};
		}
		
	}


	void advanceStep(std::vector<MPMParticle>& particles, MPMParamters& param, std::vector<Vector6d>& contactForce, int timestep) // prticle vector, timestep
	{
		// // initialize background grid nodes
		std::vector<Grid> nodesVec;
		std::map<std::string, int> gridMap;


		// calculate the reationreship between particles and background grid 
		calWeightsAndNodes(particles, param, nodesVec, gridMap);
		// transfer information from particle to grdi nodes
		particle2Grid(particles, param, nodesVec);
		// update each material particle's cauchy stress
		updateParInternalForce(particles, param);
		// apply point force
		applyPointForce(param, nodesVec, gridMap, contactForce);
		// calculate the grid node's internal force induced by particles
		calculateNodeForce(particles, param, nodesVec);
		// grid nodes momentum update
		gridUpdate(nodesVec, param);
		// transfer information back form grid to particles
		grid2Particle(particles, param, nodesVec);

	};
	

	// extract crack surface
	std::tuple<bool, objMeshFormat, objMeshFormat, std::vector<objMeshFormat>> tryToExtractCracks(std::vector<MPMParticle>& particles, MPMParamters& param, int timestep)
	{
		std::vector<extractCrackSurface::Particle> particlesRaw;
		extractCrackSurface::CRACKParamters paramCrack;
		for (int i = 0; i < particles.size(); i++)
		{
			particlesRaw.push_back(extractCrackSurface::Particle(particles[i].position, particles[i].velocity, particles[i].mass, particles[i].color, particles[i].dp));
		}
		paramCrack.dx = param.dx;


		std::tuple<bool, objMeshFormat, objMeshFormat, std::vector<objMeshFormat>> cracks = extractCrackSurface::extractCrackSurf(&particlesRaw, paramCrack);

		bool findCrackSurface = std::get<0>(cracks);

		objMeshFormat crackSurfacePartialCut;
		crackSurfacePartialCut.vertices = std::get<1>(cracks).vertices;
		crackSurfacePartialCut.facesPolygonal = std::get<1>(cracks).facesPolygonal;

		objMeshFormat crackSurfaceFullCut;
		crackSurfaceFullCut.vertices = std::get<2>(cracks).vertices;
		crackSurfaceFullCut.facesPolygonal = std::get<2>(cracks).facesPolygonal;

		std::vector<objMeshFormat> allFragmentsObj;
		for (int k = 0; k < std::get<3>(cracks).size(); k++)
		{
			objMeshFormat frag;
			frag.vertices = std::get<3>(cracks)[k].vertices;
			frag.facesPolygonal = std::get<3>(cracks)[k].facesPolygonal;
			allFragmentsObj.push_back(frag);
		}

		std::tuple<bool, objMeshFormat, objMeshFormat, std::vector<objMeshFormat>> resultReturn(findCrackSurface, crackSurfacePartialCut, crackSurfaceFullCut, allFragmentsObj);
		return resultReturn;
	}


	std::tuple<bool, objMeshFormat, objMeshFormat, std::vector<objMeshFormat>> crackSimulation(
		const std::vector<Eigen::Vector3d>& points, 
		const double& volume, 
		const Material& mat_mpm,
		MPMParamters& param,
		std::vector<Vector6d>& contactForce, 
		int num_timestep)
	{
		
		std::vector<MPMParticle> particles;
		initialize_mpm_particles(particles, points, volume, mat_mpm);


		bool generateFragments = false;
		bool simulationFailed = false;
		std::tuple<bool, objMeshFormat, objMeshFormat, std::vector<objMeshFormat>> crackSurfs;
		
		for (int timestep = 0; timestep <= 2000 && generateFragments == false && !simulationFailed; timestep++)
		{
			std::cout << "MPM timestep = "<< timestep<<" "<<std::endl;
			advanceStep(particles, param, contactForce, timestep);




			// The simulation fails due to unexpected reasons. Test
			for (int k = 0; k < particles.size(); k++)
			{
				Eigen::Vector3d scale = particles[k].position;
				if (scale.hasNaN() || scale.array().isInf().any())
				{
					std::cout << "Conataining inf or NaN, exit!" << std::endl << std::endl << std::endl;
					simulationFailed = true;
					break;
				}
			}

			if (timestep % 20 == 0 && !simulationFailed)
			{
				crackSurfs = tryToExtractCracks(particles, param, timestep);
				if (std::get<0>(crackSurfs) == true && std::get<3>(crackSurfs).size() > 1)
				{
					generateFragments = true;

					{
						std::ofstream outfile9("./output/MPM_particles_" + std::to_string(timestep) + ".txt", std::ios::trunc);
						for (int k = 0; k < particles.size(); k++)
						{
							Eigen::Vector3d scale = particles[k].position;
							outfile9 << std::scientific << std::setprecision(8) << scale[0] << " " << scale[1] << " " << scale[2]<<" "<< particles[k].dp << std::endl;
						}
						outfile9.close();
					}
				}
			}

		}
		std::cout << std::endl;

		return crackSurfs;

	}


	void initialize_mpm_particles(std::vector<MPMParticle>& particles, const std::vector<Eigen::Vector3d>& points, const double& volume, const Material& mat_mpm)
	{
		particles.clear();
		for (int i = 0; i < points.size(); i++)
		{
			MPMParticle pt;
			pt.volume = volume;
			pt.mass = volume * mat_mpm.density;
			pt.position = points[i];
			particles.push_back(pt);
		}
		
	}

}

















