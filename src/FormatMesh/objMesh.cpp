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


void objMeshFormat::updateMesh()
{
	// check if the mesh is outward or not. If not, flip the surface
	std::pair<Eigen::MatrixXd, Eigen::MatrixXi> libigl_mesh = to_libigl_mesh();
	Eigen::Vector3d center = Eigen::Vector3d::Zero();
	double volume_ = 0;
	igl::centroid(libigl_mesh.first, libigl_mesh.second, center, volume_);
	if (volume_ < 0.0)
	{
		std::cout << "Warning: the normal is not outward. flip it!" << std::endl;
		for (int f = 0; f < faces.size(); f++)
		{
			Eigen::Vector3i tmp = { faces[f][2],faces[f][1],faces[f][0] };
			faces[f] = tmp;
		}	
	}
	volume = volume_;
	//std::cout << "ct = " << ct << std::endl;
	vertFaces.clear();
	edges.clear();
	findVertFaces_Edges();
}


void objMeshFormat::readObjFile(std::string fileName, bool polygonal, Eigen::Affine3d rotation,
	Eigen::Vector3d scale, Eigen::Vector3d translation)
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
					Eigen::Vector3d vt = { std::stod(vecCoor[1]) * scale[0] , std::stod(vecCoor[2]) * scale[1] , std::stod(vecCoor[3]) * scale[2] };
					vt = rotation * vt + translation;
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
					Eigen::Vector3d vt = { std::stod(vecCoor[1]) * scale[0] , std::stod(vecCoor[2]) * scale[1] , std::stod(vecCoor[3]) * scale[2] };
					vt = rotation * vt + translation;
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
						facePolygonal.push_back(std::stoi(vecCoor[h]) - 1);
					}
					facesPolygonal.push_back(facePolygonal);

					//ct += 1;
				}

			}
		}
	}

	in.close();


	// check if the mesh is outward or not. If not, flip the surface
	std::pair<Eigen::MatrixXd, Eigen::MatrixXi> libigl_mesh = to_libigl_mesh();
	Eigen::Vector3d center = Eigen::Vector3d::Zero();
	double volume_ = 0;
	igl::centroid(libigl_mesh.first, libigl_mesh.second, center, volume_);
	if (volume_ < 0.0)
	{
		std::cout << "Warning: the normal is not outward. flip it!" << std::endl;
		if (!polygonal)
		{
			for (int f = 0; f < faces.size(); f++)
			{
				Eigen::Vector3i tmp = { faces[f][2],faces[f][1],faces[f][0] };
				faces[f] = tmp;
			}
		}
		else
		{
			for (int f = 0; f < facesPolygonal.size(); f++)
			{
				std::reverse(facesPolygonal[f].begin(), facesPolygonal[f].end());
			}
		}

	}

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


	// check if the mesh is outward or not. If not, flip the surface
	std::pair<Eigen::MatrixXd, Eigen::MatrixXi> libigl_mesh = to_libigl_mesh();
	Eigen::Vector3d center = Eigen::Vector3d::Zero();
	double volume_ = 0;
	igl::centroid(libigl_mesh.first, libigl_mesh.second, center, volume_);
	bool flip = false;
	if (volume_ < 0.0)
	{
		flip = true;
	}

	if (!polygonal)
	{
		for (int k = 0; k < faces.size(); k++)
		{
			outfile9 << "f ";
			if (!flip)
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
			if (!flip)
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
	if (faces.size() == 0 )
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

}

void objMeshFormat::to_openVDB_format(std::vector<openvdb::Vec3s>& verticesVdb, std::vector<openvdb::Vec3I>& trianglesVdb)
{
	if (faces.size() == 0)
	{
		triangulate();
	}

	verticesVdb.clear();
	trianglesVdb.clear();
	for (int i = 0; i < vertices.size(); i++)
	{
		verticesVdb.push_back({ static_cast<float>(vertices[i][0]),static_cast<float>(vertices[i][1]),static_cast<float>(vertices[i][2]) });
	}
	for (int i = 0; i < faces.size(); i++)
	{
		trianglesVdb.push_back({ static_cast<uint32_t>(faces[i][0]), static_cast<uint32_t>(faces[i][1]), static_cast<uint32_t>(faces[i][2]) });
	}
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXi> objMeshFormat::to_libigl_mesh()
{
	if (faces.size() == 0)
	{
		triangulate();
	}

	Eigen::MatrixXd V = Eigen::MatrixXd::Zero(vertices.size(), 3);
	Eigen::MatrixXi F = Eigen::MatrixXi::Zero(faces.size(), 3);

	for (int k = 0; k < vertices.size(); k++)
	{
		V.row(k) = vertices[k];
	}
	for (int k = 0; k < faces.size(); k++)
	{
		F.row(k) = faces[k];
	}

	return std::make_pair(V, F);
}

std::vector<Eigen::Vector3d> objMeshFormat::sample_points_inside_mesh(int num_samples)
{

	std::pair<Eigen::MatrixXd, Eigen::MatrixXi> libigl_mesh = to_libigl_mesh();
	Eigen::MatrixXd V = libigl_mesh.first;
	Eigen::MatrixXi F = libigl_mesh.second;


	// Initialize a random number generator
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis_x(V.col(0).minCoeff(), V.col(0).maxCoeff());
	std::uniform_real_distribution<> dis_y(V.col(1).minCoeff(), V.col(1).maxCoeff());
	std::uniform_real_distribution<> dis_z(V.col(2).minCoeff(), V.col(2).maxCoeff());

	// Store sampled points
	std::vector<Eigen::Vector3d> inside_points;


	bool sufficient = false;
	do
	{
		Eigen::VectorXd inside_check = Eigen::VectorXd::Zero(num_samples);
		std::vector<Eigen::Vector3d> inside_points_tmp(num_samples);

#pragma omp parallel for 
		for (int i = 0; i < num_samples; ++i)
		{
			Eigen::Vector3d point;
			point[0] = dis_x(gen);
			point[1] = dis_y(gen);
			point[2] = dis_z(gen);

			// Compute the winding number to check if the point is inside the mesh
			double winding_number;
			winding_number = igl::winding_number(V, F, point.transpose());

			// If winding number is close to 1, the point is inside
			if (std::abs(winding_number - 1.0) < 1e-6)
			{
				inside_points_tmp[i] = point;
				inside_check[i] = 1;
			}
		}

		for (int i = 0; i < num_samples; i++)
		{
			if (inside_check[i] == 1)
			{
				inside_points.push_back(inside_points_tmp[i]);
			}
		}

		if (inside_points.size() >= num_samples)
		{
			sufficient = true;
		}


	} while (!sufficient);


	inside_points.resize(num_samples);

	return inside_points;
}

void objMeshFormat::updateVolume()
{
	std::pair<Eigen::MatrixXd, Eigen::MatrixXi> libigl_mesh = to_libigl_mesh();
	Eigen::Vector3d center = Eigen::Vector3d::Zero();
	double volume_ = 0;
	igl::centroid(libigl_mesh.first, libigl_mesh.second, center, volume_);
	volume = std::abs(volume_);
}

objMeshFormat objMeshFormat::reconstruct_with_vdb(float& voxel_size)
{
	objMeshFormat result;

	if (faces.size() == 0)
	{
		triangulate();
	}

	std::vector<openvdb::Vec3s> vertices;
	std::vector<openvdb::Vec3I> triangles;
	to_openVDB_format(vertices, triangles);


	// define openvdb linear transformation
	openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(voxel_size);
	openvdb::FloatGrid::Ptr crackLevelSetGrid = openvdb::tools::meshToUnsignedDistanceField<openvdb::FloatGrid>(
		*transform,
		vertices,
		triangles,
		std::vector<openvdb::Vec4I>(),
		3);

	for (openvdb::FloatGrid::ValueOnIter iter = crackLevelSetGrid->beginValueOn(); iter; ++iter) {
		float dist = iter.getValue();
		float value = dist - std::sqrt(3 * std::pow(voxel_size, 2));
		iter.setValue(value);
	}
	crackLevelSetGrid->setGridClass(openvdb::GRID_LEVEL_SET);



	{
		openvdb::tools::VolumeToMesh volumeToMeshHandle;
		volumeToMeshHandle(*crackLevelSetGrid);
		openvdb::tools::PointList* verts = &volumeToMeshHandle.pointList();
		openvdb::tools::PolygonPoolList* polys = &volumeToMeshHandle.polygonPoolList();

		for (size_t i = 0; i < volumeToMeshHandle.pointListSize(); i++)
		{
			openvdb::Vec3s v = (*verts)[i];
			result.vertices.push_back(Eigen::Vector3d(v.x(), v.y(), v.z()));
		}
		for (size_t i = 0; i < volumeToMeshHandle.polygonPoolListSize(); i++) 
		{
			for (size_t ndx = 0; ndx < (*polys)[i].numQuads(); ndx++) 
			{
				openvdb::Vec4I* p = &((*polys)[i].quad(ndx));
				openvdb::Vec3I f0 = { p->z() ,p->y() ,p->x() };
				openvdb::Vec3I f1 = { p->w() ,p->z() ,p->x() };

				result.faces.push_back(Eigen::Vector3i(p->z(), p->y(), p->x()));
				result.faces.push_back(Eigen::Vector3i(p->w(), p->z(), p->x()));
			}
		}

	}

	return result;
}


CGAL_Surface_mesh objMeshFormat::to_CGAL_mesh()
{
	CGAL_Surface_mesh mesh;

	if (faces.size() == 0)
	{
		triangulate();
	}


	std::vector<CGAL_Surface_mesh::Vertex_index> vertex_indices;
	for (int i = 0; i < vertices.size(); i++)
	{
		CGAL_Surface_mesh::Vertex_index vi = mesh.add_vertex(CGAL_Point_3(vertices[i][0], vertices[i][1], vertices[i][2]));
		vertex_indices.push_back(vi);
	}

	
	for (int j = 0; j < faces.size(); j++)
	{
		std::vector<CGAL_Surface_mesh::Vertex_index> face_vertices;
		face_vertices.push_back(vertex_indices[faces[j][0]]);
		face_vertices.push_back(vertex_indices[faces[j][1]]);
		face_vertices.push_back(vertex_indices[faces[j][2]]);

		mesh.add_face(face_vertices);
	}

	return mesh;
}

objMeshFormat objMeshFormat::boolean_difference_with_mesh(objMeshFormat& B_)
{
	objMeshFormat mesh_diff;

	CGAL_Surface_mesh A = to_CGAL_mesh();
	CGAL_Surface_mesh B = B_.to_CGAL_mesh();
	CGAL_Surface_mesh result;

	// 执行布尔操作 A - B
	if (!PMP::corefine_and_compute_difference(A, B, result)) 
	{
		std::cerr << "Boolean operation failed!" << std::endl;
		std::exit(0);
	}


	// 映射顶点描述符到索引
	std::map<CGAL_Surface_mesh::Vertex_index, int> vertex_index_map;

	// 遍历顶点并存储到 Eigen::Vector3d
	int index = 0;
	for (auto vd : CGAL::vertices(result)) {
		CGAL_Point_3 p = result.point(vd);
		mesh_diff.vertices.emplace_back(p.x(), p.y(), p.z());
		vertex_index_map[vd] = index++; // 将顶点描述符映射到索引
	}

	// 遍历面并存储顶点索引
	for (auto fd : CGAL::faces(result)) {
		std::vector<int> face_indices;
		for (auto vd : vertices_around_face(result.halfedge(fd), result)) {
			face_indices.push_back(vertex_index_map[vd]); // 获取顶点索引
		}
		mesh_diff.facesPolygonal.push_back(face_indices);
	}

	mesh_diff.triangulate();

	return mesh_diff;

}
