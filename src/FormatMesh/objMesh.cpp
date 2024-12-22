#include "objMesh.h"



// Function to compute face normal
Eigen::Vector3d computeNormal(const Eigen::Vector3d& v0,
	const Eigen::Vector3d& v1,
	const Eigen::Vector3d& v2)
{
	return (v1 - v0).cross(v2 - v0).normalized();
}

// Function to compute projection of points onto an axis
void projectOntoAxis(const std::vector<Eigen::Vector3d>& points,
	const Eigen::Vector3d& axis,
	double& minProj,
	double& maxProj)
{
	minProj = std::numeric_limits<double>::infinity();
	maxProj = -std::numeric_limits<double>::infinity();
	for (const auto& p : points)
	{
		double proj = axis.dot(p);
		if (proj < minProj) minProj = proj;
		if (proj > maxProj) maxProj = proj;
	}
}

// Triangle-Tetrahedron intersection test using SAT
bool triangleTetrahedronIntersect(const Eigen::Vector3d& tri_v0,
	const Eigen::Vector3d& tri_v1,
	const Eigen::Vector3d& tri_v2,
	const Eigen::Vector3d& tet_v0,
	const Eigen::Vector3d& tet_v1,
	const Eigen::Vector3d& tet_v2,
	const Eigen::Vector3d& tet_v3)
{
	// Collect vertices
	std::vector<Eigen::Vector3d> triVerts = { tri_v0, tri_v1, tri_v2 };
	std::vector<Eigen::Vector3d> tetVerts = { tet_v0, tet_v1, tet_v2, tet_v3 };

	// Edges of the triangle
	std::vector<Eigen::Vector3d> triEdges = {
		tri_v1 - tri_v0,
		tri_v2 - tri_v1,
		tri_v0 - tri_v2
	};

	// Edges of the tetrahedron
	std::vector<Eigen::Vector3d> tetEdges = {
		tet_v1 - tet_v0,
		tet_v2 - tet_v0,
		tet_v3 - tet_v0,
		tet_v2 - tet_v1,
		tet_v3 - tet_v1,
		tet_v3 - tet_v2
	};

	// Normals of the tetrahedron faces
	std::vector<Eigen::Vector3d> tetNormals = {
		(tet_v1 - tet_v0).cross(tet_v2 - tet_v0).normalized(),
		(tet_v2 - tet_v0).cross(tet_v3 - tet_v0).normalized(),
		(tet_v3 - tet_v0).cross(tet_v1 - tet_v0).normalized(),
		(tet_v1 - tet_v2).cross(tet_v3 - tet_v2).normalized()
	};

	// Triangle normal
	Eigen::Vector3d triNormal = computeNormal(tri_v0, tri_v1, tri_v2);

	// Axes to test
	std::vector<Eigen::Vector3d> axes;

	// Add triangle normal
	axes.push_back(triNormal);

	// Add tetrahedron face normals
	axes.insert(axes.end(), tetNormals.begin(), tetNormals.end());

	// Add cross products of edges
	for (const auto& e1 : triEdges)
	{
		for (const auto& e2 : tetEdges)
		{
			Eigen::Vector3d axis = e1.cross(e2);
			if (axis.norm() > 1e-8)  // Avoid zero-length axes
				axes.push_back(axis.normalized());
		}
	}

	// Test axes
	for (const auto& axis : axes)
	{
		double triMin, triMax;
		double tetMin, tetMax;
		projectOntoAxis(triVerts, axis, triMin, triMax);
		projectOntoAxis(tetVerts, axis, tetMin, tetMax);

		// Check for separation
		if (triMax < tetMin || tetMax < triMin)
		{
			// No overlap on this axis, shapes are separated
			return false;
		}
	}

	// No separating axis found, shapes intersect
	return true;
}


// Möller–Trumbore intersection algorithm for line segment and triangle
bool lineSegmentIntersectsTriangle(const Eigen::Vector3d& orig,
	const Eigen::Vector3d& dest,
	const Eigen::Vector3d& vert0,
	const Eigen::Vector3d& vert1,
	const Eigen::Vector3d& vert2,
	Eigen::Vector3d& intersectPoint)
{
	const double EPSILON = 1e-8;

	Eigen::Vector3d edge1 = vert1 - vert0;
	Eigen::Vector3d edge2 = vert2 - vert0;
	Eigen::Vector3d direction = dest - orig;
	Eigen::Vector3d pvec = direction.cross(edge2);
	double det = edge1.dot(pvec);

	// If the determinant is near zero, the line lies in the plane of the triangle or is parallel to it
	if (std::abs(det) < EPSILON)
		return false;

	double invDet = 1.0 / det;
	Eigen::Vector3d tvec = orig - vert0;
	double u = tvec.dot(pvec) * invDet;
	if (u < 0.0 || u > 1.0)
		return false;

	Eigen::Vector3d qvec = tvec.cross(edge1);
	double v = direction.dot(qvec) * invDet;
	if (v < 0.0 || u + v > 1.0)
		return false;

	double t = edge2.dot(qvec) * invDet;
	if (t < 0.0 || t > 1.0)
		return false;

	// Intersection point
	intersectPoint = orig + t * direction;
	return true;
}


void objMeshFormat::readObjFile(std::string fileName, bool polygonal)
{
	clear();

	//int ct = 0;
	// read vertices and faces
	std::ifstream in;
	in.open(fileName);
	std::string line;
	if (!polygonal)
	{
		while (getline(in, line))
		{
			if (line.size() > 0)
			{
				std::vector<std::string> vecCoor = split(line, " ");
				if (vecCoor[0] == "v")
				{
					Eigen::Vector3d vt = { std::stod(vecCoor[1]) , std::stod(vecCoor[2]) , std::stod(vecCoor[3]) };
					vertices.push_back(vt);
				}

				if (vecCoor[0] == "f")
				{
					Eigen::Vector3i ele = { std::stoi(vecCoor[1]) - 1 ,  std::stoi(vecCoor[2]) - 1 ,  std::stoi(vecCoor[3]) - 1 };
					faces.push_back(ele);
				}

			}
		}
	}
	else
	{
		while (getline(in, line))
		{
			if (line.size() > 0)
			{
				std::vector<std::string> vecCoor = split(line, " ");
				if (vecCoor[0] == "v")
				{
					Eigen::Vector3d vt = { std::stod(vecCoor[1]) , std::stod(vecCoor[2]) , std::stod(vecCoor[3]) };
					vertices.push_back(vt);
				}

				if (vecCoor[0] == "f")
				{
					//std::cout << "vecCoor = " << line << std::endl;
					//std::cout << "vecCoor.size() = " << vecCoor.size() << std::endl;
					std::vector<int> facePolygonal;
					for (int h = 1; h < vecCoor.size(); h++)
					{
						//std::cout << "	ele = " << std::stoi(vecCoor[h]) - 1 << std::endl;
						facePolygonal.push_back(std::stoi(vecCoor[h]) - 1 );				
					}
					facesPolygonal.push_back(facePolygonal);

					//ct += 1;
				}

			}
		}
	}

	in.close();

	//std::cout << "ct = " << ct << std::endl;

	findVertFaces_Edges();

}

void objMeshFormat::findVertFaces_Edges()
{
	std::set<std::string> edgesStored;
	for (int i = 0; i < faces.size(); i++)
	{
		int v1 = faces[i][0], v2 = faces[i][1], v3 = faces[i][2];
		std::string v1v2 = std::to_string(std::min(v1, v2)) + std::to_string(std::max(v1, v2));
		std::string v1v3 = std::to_string(std::min(v1, v3)) + std::to_string(std::max(v1, v3));
		std::string v2v3 = std::to_string(std::min(v2, v3)) + std::to_string(std::max(v2, v3));
		if (edgesStored.find(v1v2) == edgesStored.end())
		{
			edgesStored.insert(v1v2);
			Eigen::Vector2i eg = { std::min(v1, v2) ,  std::max(v1, v2) };
			edges.push_back(eg);
		}

		if (edgesStored.find(v1v3) == edgesStored.end())
		{
			edgesStored.insert(v1v3);
			Eigen::Vector2i eg = { std::min(v1, v3) ,  std::max(v1, v3) };
			edges.push_back(eg);
		}

		if (edgesStored.find(v2v3) == edgesStored.end())
		{
			edgesStored.insert(v2v3);
			Eigen::Vector2i eg = { std::min(v2, v3) ,  std::max(v2, v3) };
			edges.push_back(eg);
		}
	}


	// find vertFaces
	vertFaces.resize(vertices.size());
	for (int fi = 0; fi < faces.size(); fi++)
	{
		Eigen::Vector3i fc = faces[fi];
		vertFaces[fc[0]].push_back(fi);
		vertFaces[fc[1]].push_back(fi);
		vertFaces[fc[2]].push_back(fi);
	}


}

void depthFirstSearch(int v, const std::vector<std::set<int>>& adjList, 
	std::vector<bool>& visited, std::vector<int>& component) 
{
	visited[v] = true;
	component.push_back(v);
	for (int neighbor : adjList[v]) 
	{
		if (!visited[neighbor]) 
		{
			depthFirstSearch(neighbor, adjList, visited, component);
		}
	}
}

void objMeshFormat::sepConnectedComponents()
{

	// find connected components
	std::vector<std::vector<int>> components;
	std::vector<std::set<int>> adjList(vertices.size());
	for (const auto& f : faces) 
	{
		adjList[f[0]].insert(f[1]);
		adjList[f[0]].insert(f[2]);
		adjList[f[1]].insert(f[0]);
		adjList[f[1]].insert(f[2]);
		adjList[f[2]].insert(f[0]);
		adjList[f[2]].insert(f[1]);
	}
	std::vector<bool> visited(vertices.size(), false);
	for (size_t i = 0; i < vertices.size(); ++i) 
	{
		if (!visited[i]) 
		{
			std::vector<int> component;
			depthFirstSearch(i, adjList, visited, component);
			components.push_back(component);
		}
	}

	// remove unnecessary vertices
	for (size_t i = 0; i < components.size(); ++i) 
	{
		objMeshFormat cop;
		std::unordered_map<int, int> vertexMap;

		for (int idx : components[i]) {
			vertexMap[idx] = cop.vertices.size();
			cop.vertices.push_back(vertices[idx]);
		}

		for (const auto& f : faces) 
		{
			if (vertexMap.count(f[0]) && vertexMap.count(f[1]) && vertexMap.count(f[2]))
			{
				cop.faces.push_back({ vertexMap[f[0]], vertexMap[f[1]], vertexMap[f[2]] });
			}
		}
		
		cop.findVertFaces_Edges();
		componentsSep.push_back(cop);
	}

}

void objMeshFormat::clear()
{
	vertices.clear();
	faces.clear();
	vertFaces.clear();
	componentsSep.clear();
	facesPolygonal.clear();
}

void objMeshFormat::outputFile(std::string fileName, int timestep, bool polygonal)
{
	std::ofstream outfile9("./output/" + fileName + "_" + std::to_string(timestep) + ".obj", std::ios::trunc);
	for (int k = 0; k < vertices.size(); k++)
	{
		Eigen::Vector3d scale = vertices[k];
		outfile9 << std::scientific << std::setprecision(8) << "v " << scale[0] << " " << scale[1] << " " << scale[2] << std::endl;
	}
	if (!polygonal)
	{
		for (int k = 0; k < faces.size(); k++)
		{
			outfile9 << "f ";
			if (0)
			{
				for (int m = 0; m < faces[k].size(); m++)
				{
					outfile9 << faces[k][m] + 1 << " ";
				}
			}
			else
			{
				for (int m = faces[k].size() - 1; m >= 0; m--)
				{
					outfile9 << faces[k][m] + 1 << " ";
				}
			}
			outfile9 << std::endl;
		}
	}
	else
	{
		for (int k = 0; k < facesPolygonal.size(); k++)
		{
			outfile9 << "f ";
			if (0)
			{
				for (int m = 0; m < facesPolygonal[k].size(); m++)
				{
					outfile9 << facesPolygonal[k][m] + 1 << " ";
				}
			}
			else
			{
				for (int m = facesPolygonal[k].size() - 1; m >= 0; m--)
				{
					outfile9 << facesPolygonal[k][m] + 1 << " ";
				}
			}
			outfile9 << std::endl;
		}
	}	
	outfile9.close();
}

bool objMeshFormat::checkIfMeshIntersectWithTetrahedron(const Eigen::Vector3d tet_v0, const Eigen::Vector3d tet_v1,
	const Eigen::Vector3d tet_v2, const Eigen::Vector3d tet_v3)
{
	// Iterate over each face
	for (const auto& face : facesPolygonal)
	{
		// Triangulate face if it has more than 3 vertices
		for (size_t i = 1; i + 1 < face.size(); ++i)
		{
			Eigen::Vector3d tri_v0 = vertices[face[0]];
			Eigen::Vector3d tri_v1 = vertices[face[i]];
			Eigen::Vector3d tri_v2 = vertices[face[i + 1]];

			// Check for intersection
			if (triangleTetrahedronIntersect(tri_v0, tri_v1, tri_v2,
				tet_v0, tet_v1, tet_v2, tet_v3))
			{
				// Intersection found
				return true;
			}
		}
	}

	// No intersection found
	return false;
}

bool objMeshFormat::checkIfMeshIntersectWithLine(const Eigen::Vector3d line_pt1, const Eigen::Vector3d line_pt2)
{
	// Iterate over each face
	for (const auto& face : facesPolygonal)
	{
		// Triangulate face if it has more than 3 vertices
		for (size_t i = 1; i + 1 < face.size(); ++i)
		{
			Eigen::Vector3d tri_v0 = vertices[face[0]];
			Eigen::Vector3d tri_v1 = vertices[face[i]];
			Eigen::Vector3d tri_v2 = vertices[face[i + 1]];

			// Check for intersection
			Eigen::Vector3d intersectPoint;
			if (lineSegmentIntersectsTriangle(line_pt1, line_pt2,
				tri_v0, tri_v1, tri_v2,
				intersectPoint))
			{
				// Intersection found
				return true;
			}
		}
	}

	// No intersection found
	return false;
}

void objMeshFormat::triangulate()
{
	for (int i = 0; i < facesPolygonal.size(); i++)
	{
		std::vector<int> facePolygonal = facesPolygonal[i];
		if (facePolygonal.size() > 3)
		{
			for (int j = 1; j < facePolygonal.size() - 1; j++)
			{
				Eigen::Vector3i faceTri = { facePolygonal[0], facePolygonal[j], facePolygonal[j + 1] };
				faces.push_back(faceTri);
			}
		}
		else
		{
			Eigen::Vector3i faceTri = { facePolygonal[0], facePolygonal[1], facePolygonal[2] };
			faces.push_back(faceTri);
		}
	}
}

void objMeshFormat::to_openVDB_format(std::vector<openvdb::Vec3s>& verticesVdb, std::vector<openvdb::Vec3I>& trianglesVdb)
{
	verticesVdb.clear();
	trianglesVdb.clear();
	for (int i = 0 ; i < vertices.size(); i++)
	{
		verticesVdb.push_back({ static_cast<float>(vertices[i][0]),static_cast<float>(vertices[i][1]),static_cast<float>(vertices[i][2]) });
	}
	for (int i = 0; i < faces.size(); i++)
	{
		trianglesVdb.push_back({ static_cast<uint32_t>(faces[i][0]), static_cast<uint32_t>(faces[i][1]), static_cast<uint32_t>(faces[i][2]) });
	}
}