#include "simulator.h" 
#include "tools.h"
#include "AABB.h"



// TODO
// 1. Accelerate and parallel openMP
// 2. Decimate the generated mesh for performance
// 3. triSimMesh.boundaryCondition_node_surface (& boundaryCondition_node_interior) is not properly handled
// 4. There are penetration occasionally
// 5. No friction implementation
// 6. contact_Info can be simplified
// 7. Remove redundant functions in simulator.h
// 8. surface_Info may not be necessary
// 9. change the step update rule to per object base
// 10. change the timestep size according to the solution automatically





int main()
{

	if (0)
	{
		
		objMeshFormat bunny;
		bunny.readObjFile("D:/Research/Hydrostatic_object/code/FEM_IPC/build/output/bunny.obj");
		bunny.updateBEInfo();

		objMeshFormat crack;
		crack.readObjFile("D:/Research/Hydrostatic_object/code/FEM_IPC/build/output/crack.obj");
		crack.updateBEInfo();


		float voxel_size = 0.01f / 2.0 / std::sqrt(3);

		objMeshFormat crack_surf = crack.reconstruct_with_vdb(voxel_size);
		crack_surf.outputFile("crack_volumetric");
		crack_surf.decimate_mesh(50000);



		//std::pair<Eigen::MatrixXd, Eigen::MatrixXi> VF = crack_surf.to_libigl_mesh();
		//Eigen::MatrixXd V = VF.first;
		//Eigen::MatrixXi F = VF.second; 
		//Eigen::MatrixXd U; // 简化后的顶点
		//Eigen::MatrixXi G; // 简化后的面
		//Eigen::VectorXi J; // 面的映射
		//Eigen::VectorXi I; // 顶点的映射
		//bool block_intersections = false;
		//if (!igl::decimate(V, F, 50000, block_intersections, U, G, J, I)) {
		//	std::cerr << "Decimation failed!" << std::endl;
		//	return 1;
		//}
		//std::cout << "Simplified mesh: " << G.rows() << " faces, " << U.rows() << " vertices." << std::endl;
		//// 保存简化后的网格
		//if (!igl::writeOBJ("./output/reconstructed_surface_decimated.obj", U, G)) {
		//	std::cerr << "Failed to save simplified mesh to " << std::endl;
		//	return 1;
		//}
		//objMeshFormat crack_surf_decimate;
		//crack_surf_decimate.build_with_libigl(U, G);




		objMeshFormat children = bunny.boolean_difference_with_mesh(crack_surf);
		children.triangulate();
		children.sepConnectedComponents();
		for (int i = 0; i < children.componentsSep.size(); i++)
		{
			children.componentsSep[i].outputFile("child_"+std::to_string(i));
		}



		//// 定义物体 A 的顶点和面片
		//std::vector<Eigen::Vector3d> verticesA = {
		//	{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}
		//};
		//std::vector<Eigen::Vector3i> facesA = {
		//	{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}
		//};

		//// 定义物体 B 的顶点和面片
		//std::vector<Eigen::Vector3d> verticesB = {
		//	{0.5, 0.5, 0.5}, {1.5, 0.5, 0.5}, {0.5, 1.5, 0.5}, {0.5, 0.5, 1.5}
		//};
		//std::vector<Eigen::Vector3i> facesB = {
		//	{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}
		//};

		//// 构建 BVH
		//RigidBody bodyA, bodyB;
		//bodyA.bvh = buildBVH(verticesA, facesA, 0, facesA.size());
		//bodyB.bvh = buildBVH(verticesB, facesB, 0, facesB.size());

		//// 设置仿射变换
		//bodyA.transform = Eigen::Affine3d(Eigen::Translation3d(0, 0, 0));
		//bodyA.inverseTransform = bodyA.transform.inverse();

		//bodyB.transform = Eigen::Affine3d(Eigen::Translation3d(0.25, 0.25, 0.25));
		//bodyB.inverseTransform = bodyB.transform.inverse();

		//// 碰撞检测
		//if (checkCollision(bodyA, bodyB)) {
		//	std::cout << "Collision detected!" << std::endl;
		//}
		//else {
		//	std::cout << "No collision." << std::endl;
		//}








		//std::string bunny = "../input/bunny_highRes.obj";
		//std::string cube = "../input/bunny_highRes.obj";


		//objMeshFormat bunnyMesh;
		//bunnyMesh.readObjFile(bunny, true);


		//objMeshFormat cubeMesh;
		//Eigen::Vector3d trans = {0,1.0,0};
		//cubeMesh.readObjFile(cube, true, Eigen::Affine3d::Identity(), Eigen::Vector3d::Ones(), trans);

		//objMeshFormat diff = bunnyMesh.boolean_difference_with_mesh(cubeMesh);

		//diff.outputFile("../output/bunny_cube_diff.obj", -99, false);
		//










		//// 为网格 A 和 B 的顶点添加来源属性
		//std::map<SurfaceMesh::Vertex_index, int> vertex_source_map_A;
		//std::map<SurfaceMesh::Vertex_index, int> vertex_source_map_B;

		//for (auto v : bunny.vertices()) {
		//	vertex_source_map_A[v] = 1; // 标记 A 的顶点来源为 1
		//}

		//for (auto v : cube.vertices()) {
		//	vertex_source_map_B[v] = 2; // 标记 B 的顶点来源为 2
		//}

		//// 执行布尔操作 A - B
		//if (!PMP::corefine_and_compute_difference(bunny, cube, result,
		//	CGAL::parameters::vertex_point_map(get(CGAL::vertex_point, bunny)),
		//	CGAL::parameters::vertex_point_map(get(CGAL::vertex_point, cube)))) {
		//	std::cerr << "Boolean operation failed!" << std::endl;
		//	return 1;
		//}

		//// 遍历结果网格的顶点，检查顶点的来源
		//for (auto v : result.vertices()) {
		//	Point p = result.point(v);

		//	// 检查顶点是否在原始网格 A 的顶点映射中
		//	auto it = vertex_source_map_A.find(v);
		//	if (it != vertex_source_map_A.end() && it->second == 1) {
		//		std::cout << "Vertex " << p << " comes from mesh A." << std::endl;
		//	}
		//	else {
		//		std::cout << "Vertex " << p << " does not come from mesh A." << std::endl;
		//	}
		//}




		//// Input .obj file
		//std::string input_filename = "D:/Research/Hydrostatic_object/code/FEM_IPC/input/bunnyCrack.obj";


		//objMeshFormat crackMesh;
		//crackMesh.readObjFile(input_filename, true);
		//crackMesh.triangulate();

		//std::vector<openvdb::Vec3s> vertices;
		//std::vector<openvdb::Vec3I> triangles;
		//crackMesh.to_openVDB_format(vertices, triangles);

		//{
		//	std::ofstream outfile9("./output/original_surface.obj", std::ios::trunc);
		//	for (int k = 0; k < vertices.size(); k++)
		//	{
		//		openvdb::Vec3s scale = vertices[k];
		//		outfile9 << std::scientific << std::setprecision(8) << "v " << scale.x() << " " << scale.y() << " " << scale.z() << std::endl;
		//	}
		//	for (int k = 0; k < triangles.size(); k++)
		//	{
		//		openvdb::Vec3I scale = triangles[k];
		//		outfile9 << std::scientific << std::setprecision(8) << "f " << scale[0] + 1 << " " << scale[1] + 1 << " " << scale[2] + 1 << std::endl;
		//	}
		//	outfile9.close();
		//}




		//float voxel_size = 0.04f;
		//// define openvdb linear transformation
		//openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(voxel_size);
		//openvdb::FloatGrid::Ptr crackLevelSetGrid = openvdb::tools::meshToUnsignedDistanceField<openvdb::FloatGrid>(
		//	*transform,
		//	vertices,
		//	triangles,
		//	std::vector<openvdb::Vec4I>(),
		//	3);

		//for (openvdb::FloatGrid::ValueOnIter iter = crackLevelSetGrid->beginValueOn(); iter; ++iter) {
		//	float dist = iter.getValue();
		//	float value = dist - std::sqrt(3 * std::pow(voxel_size, 2));
		//	iter.setValue(value);
		//}
		//crackLevelSetGrid->setGridClass(openvdb::GRID_LEVEL_SET);


		//std::vector<openvdb::Vec3s> surfaceVertices; // List of surface vertices
		//std::vector<openvdb::Vec3I> surfaceTriangles;    // List of surface quads (faces)




		//{
		//	openvdb::tools::VolumeToMesh volumeToMeshHandle;
		//	volumeToMeshHandle(*crackLevelSetGrid);


		//	openvdb::tools::PointList* verts = &volumeToMeshHandle.pointList();
		//	openvdb::tools::PolygonPoolList* polys = &volumeToMeshHandle.polygonPoolList();





		//	for (size_t i = 0; i < volumeToMeshHandle.pointListSize(); i++)
		//	{
		//		openvdb::Vec3s v = (*verts)[i];
		//		surfaceVertices.push_back(v);
		//	}

		//	for (size_t i = 0; i < volumeToMeshHandle.polygonPoolListSize(); i++) {

		//		for (size_t ndx = 0; ndx < (*polys)[i].numTriangles(); ndx++) {
		//			openvdb::Vec3I* p = &((*polys)[i].triangle(ndx));
		//		}

		//		for (size_t ndx = 0; ndx < (*polys)[i].numQuads(); ndx++) {
		//			openvdb::Vec4I* p = &((*polys)[i].quad(ndx));

		//			openvdb::Vec3I f0 = { p->z() ,p->y() ,p->x() };
		//			openvdb::Vec3I f1 = { p->w() ,p->z() ,p->x() };
		//			surfaceTriangles.push_back(f0);
		//			surfaceTriangles.push_back(f1);
		//		}
		//	}


		//}



		//objMeshFormat rms;
		//{
		//	std::ofstream outfile9("./output/reconstructed_surface.obj", std::ios::trunc);
		//	for (int k = 0; k < surfaceVertices.size(); k++)
		//	{
		//		openvdb::Vec3s scale = surfaceVertices[k];
		//		outfile9 << std::scientific << std::setprecision(8) << "v " << scale.x() << " " << scale.y() << " " << scale.z() << std::endl;

		//		Eigen::Vector3d vt = { scale.x() , scale.y() , scale.z() };
		//		rms.vertices.push_back(vt);

		//	}
		//	for (int k = 0; k < surfaceTriangles.size(); k++)
		//	{
		//		openvdb::Vec3I scale = surfaceTriangles[k];
		//		outfile9 << std::scientific << std::setprecision(8) << "f " << scale[0] + 1 << " " << scale[1] + 1 << " " << scale[2] + 1 << std::endl;


		//		Eigen::Vector3i ft = { static_cast<int>(scale.x()) , static_cast<int>(scale.y()) ,static_cast<int>(scale.z()) };
		//		rms.faces.push_back(ft);

		//	}
		//	outfile9.close();
		//}



		//Eigen::MatrixXd V(surfaceVertices.size(), 3);
		//Eigen::MatrixXi F(surfaceTriangles.size(), 3);
		//{
		//	for (int k = 0; k < surfaceVertices.size(); k++)
		//	{
		//		openvdb::Vec3s scale = surfaceVertices[k];
		//		V(k, 0) = scale.x();
		//		V(k, 1) = scale.y();
		//		V(k, 2) = scale.z();
		//	}
		//	for (int k = 0; k < surfaceTriangles.size(); k++)
		//	{
		//		openvdb::Vec3I scale = surfaceTriangles[k];
		//		F(k, 0) = scale.x();
		//		F(k, 1) = scale.y();
		//		F(k, 2) = scale.z();
		//	}
		//}
		//// 简化网格


		//std::vector<Eigen::Vector3d> pts = rms.sample_points_inside_mesh(100000);
		//{
		//	std::ofstream outfile9("./output/sampledPoints.obj", std::ios::trunc);
		//	for (int k = 0; k < pts.size(); k++)
		//	{
		//		outfile9 << std::scientific << std::setprecision(8) << "v " << pts[k][0] << " " << pts[k][1] << " " << pts[k][2] << std::endl;
		//	}
		//	outfile9.close();
		//}


		//Eigen::MatrixXd U; // 简化后的顶点
		//Eigen::MatrixXi G; // 简化后的面
		//Eigen::VectorXi J; // 面的映射
		//Eigen::VectorXi I; // 顶点的映射
		//bool block_intersections = false;

		//if (!igl::decimate(V, F, 50000, block_intersections, U, G, J, I)) {
		//	std::cerr << "Decimation failed!" << std::endl;
		//	return 1;
		//}

		//std::cout << "Simplified mesh: " << G.rows() << " faces, " << U.rows() << " vertices." << std::endl;

		//// 保存简化后的网格
		//if (!igl::writeOBJ("./output/reconstructed_surface_decimated.obj", U, G)) {
		//	std::cerr << "Failed to save simplified mesh to " << std::endl;
		//	return 1;
		//}























		//float voxel_size = 0.01;
		//openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(voxel_size);



		//openvdb::initialize();
		//std::vector<openvdb::Vec4I> quads;
		//// Step 3: Convert the triangular mesh into a signed distance field (SDF)
		//openvdb::FloatGrid::Ptr sdfGrid = openvdb::tools::meshToSignedDistanceField<openvdb::FloatGrid>
		//	(*transform,vertices,triangles,quads,1.0,1.0);



		//// Step 4: Extract the surface mesh from the voxel grid
		//std::vector<openvdb::Vec3s> surfaceVertices; // List of surface vertices
		//std::vector<openvdb::Vec3I> surfaceTriangles;    // List of surface quads (faces)
		//std::vector<openvdb::Vec4I> quads_out;

		//openvdb::tools::volumeToMesh(
		//	*sdfGrid,         // Input SDF grid
		//	surfaceVertices,  // Output vertices of the surface mesh
		//	surfaceTriangles,     // Output quads of the surface mesh
		//	quads_out,
		//	0.0,              // Iso-value (surface level set, 0.0 for the zero level set)
		//	1.0               // Adaptivity (0 = full detail, higher = more simplified mesh)
		//);



		//// Step 5: Output the surface mesh (e.g., to a file or console)
		//std::cout << "Extracted vertices: "<< surfaceVertices.size() << std::endl;
		//std::cout << "Extracted triangles: "<< surfaceTriangles.size() << std::endl;
		//std::cout << "Extracted quads_out: "<< quads_out.size() << std::endl;


		//{
		//	std::ofstream outfile9("./output/reconstructed_surface.obj", std::ios::trunc);
		//	for (int k = 0; k < surfaceVertices.size(); k++)
		//	{
		//		openvdb::Vec3s scale = surfaceVertices[k];
		//		outfile9 << std::scientific << std::setprecision(8) << "v " << scale.x() << " " << scale.y() << " " << scale.z() << std::endl;
		//	}
		//	//for (int k = 0; k < surfaceTriangles.size(); k++)
		//	//{
		//	//	openvdb::Vec3I scale = surfaceTriangles[k];
		//	//	outfile9 << std::scientific << std::setprecision(8) << "f " << scale[0] + 1 << " " << scale[1] + 1 << " " << scale[2] + 1 << std::endl;
		//	//}
		//	for (int k = 0; k < quads_out.size(); k++)
		//	{
		//		openvdb::Vec4I scale = quads_out[k];
		//		outfile9 << std::scientific << std::setprecision(8) << "f " << scale[0] + 1 << " " << scale[1] + 1 << " " << scale[2] + 1 << " " << scale[3] + 1 << std::endl;
		//	}
		//	outfile9.close();
		//}



































		//// Load meshes
		//Eigen::MatrixXd VA, VB; // Vertices of meshes A and B
		//Eigen::MatrixXi FA, FB; // Faces of meshes A and B


		//if (!igl::readOBJ("D:/Research/Hydrostatic_object/code/FEM_IPC/input/cube.obj", VB, FB)) {
		//	std::cerr << "Failed to load meshB.obj" << std::endl;
		//	return -1;
		//}

		//if (!igl::readOBJ("D:/Research/Hydrostatic_object/code/FEM_IPC/input/bunny_highRes.obj", VA, FA)) {
		//	std::cerr << "Failed to load meshA.obj" << std::endl;
		//	return -1;
		//}


		//// Output result mesh
		//Eigen::MatrixXd VC; // Resulting vertices
		//Eigen::MatrixXi FC; // Resulting faces

		//// Perform Boolean operation (e.g., UNION)
		//igl::copyleft::cgal::mesh_boolean(VA, FA, VB, FB, igl::MESH_BOOLEAN_TYPE_MINUS, VC, FC);

		//std::cout << "Boolean operation completed!" << std::endl;

		//// Save the result into an OBJ file
		//if (!igl::writeOBJ("./output/result.obj", VC, FC)) {
		//	std::cerr << "Failed to save result.obj" << std::endl;
		//	return -1;
		//}

		//std::cout << "Result saved to result.obj" << std::endl;















		//objMeshFormat start;
		//start.readObjFile("./input/surfMesh_0.obj", false);

		//objMeshFormat end;
		//end.readObjFile("./input/surfMesh_2900.obj", false);
		//std::cout << "start.vertices.size() = " << start.vertices.size() << std::endl;

		//std::vector<std::pair<double, double>> displacement;
		//for (int i = 0; i < start.vertices.size(); i++)
		//{
		//	Eigen::Vector3d pos = start.vertices[i];
		//	if (std::abs(pos[1] + 10) <= 0.0001 && std::abs(pos[2] - 1) <= 0.0001 && pos[0] >= -22)
		//	{
		//		std::pair<double, double> pa = { pos[0] + 22, -(end.vertices[i][2] - start.vertices[i][2])};
		//		displacement.push_back(pa);
		//	}
		//}


		//// 使用 std::sort 并提供一个自定义的比较函数
		//std::sort(displacement.begin(), displacement.end(),
		//	[](const std::pair<double, double>& a, const std::pair<double, double>& b) {
		//		return a.first < b.first; // 按照第一个元素排序
		//	});

		//// 输出排序后的结果
		//std::cout << "Sorted displacement:" << std::endl;
		//for (const auto& p : displacement) {
		//	std::cout << "(" << p.first << ", " << p.second << ")" << std::endl;
		//}



		//std::ofstream outfile9("./output/comp.txt", std::ios::trunc);
		//for (int k = 0; k < displacement.size(); k++)
		//{
		//	outfile9 << std::scientific << std::setprecision(8) << displacement[k].first << " " << displacement[k].second << std::endl;
		//}
		//outfile9.close();















		//Eigen::Vector3d tet_v0 = { 0, 2, 2.8284271247 };
		//Eigen::Vector3d tet_v1 = { 1.7320508076, - 1, 2.8284271247 };
		//Eigen::Vector3d tet_v2 = { -1.7320508076, - 1, 2.8284271247 };
		//Eigen::Vector3d tet_v3 = { 0, 0, 0 };



		//objMeshFormat testMesh;
		////Eigen::Vector3d p1 = {-4,-4,1};
		////Eigen::Vector3d p2 = {4,-4,1};
		////Eigen::Vector3d p3 = {4,4,1};
		////Eigen::Vector3d p4 = {-4,4,1};

		////Eigen::Vector3d p1 = { -0.1,-0.1,1 };
		////Eigen::Vector3d p2 = { 0.1,-0.1,1 };
		////Eigen::Vector3d p3 = { 0.1,0.1,1 };
		////Eigen::Vector3d p4 = { -0.1,0.1,1 };

		////Eigen::Vector3d p1 = { -0.1,-0.1,1 };
		////Eigen::Vector3d p2 = { 0.1,-0.1,1 };
		////Eigen::Vector3d p3 = { 0.1,1,1 };
		////Eigen::Vector3d p4 = { -0.1,1,1 };

		////Eigen::Vector3d p1 = { -0.1,5,1 };
		////Eigen::Vector3d p2 = { 0.1,5,1 };
		////Eigen::Vector3d p3 = { 0.1,1,1 };
		////Eigen::Vector3d p4 = { -0.1,1,1 };

		//Eigen::Vector3d p1 = { -0.1,5,0 };
		//Eigen::Vector3d p2 = { 0.1,5,0 };
		//Eigen::Vector3d p3 = { 0.1,1,0 };
		//Eigen::Vector3d p4 = { -0.1,1,0 };

		//testMesh.vertices.push_back(p1);
		//testMesh.vertices.push_back(p2);
		//testMesh.vertices.push_back(p3);
		//testMesh.vertices.push_back(p4);



		//std::vector<int> fc;
		//fc.push_back(0);
		//fc.push_back(1);
		//fc.push_back(2);
		//fc.push_back(3);
		//testMesh.facesPolygonal.push_back(fc);


		//testMesh.outputFile("faces",-99);

		//bool intersect = testMesh.checkIfMeshIntersectWithTetrahedron(tet_v0, tet_v1, tet_v2, tet_v3);
		//std::cout << "intersect = " << intersect << std::endl;








		/*FEMParamters parameters;
		parameters.IPC_dis = 0.01;
		parameters.IPC_kStiffness = 1.0e16;
		parameters.numOfThreads = 24;
		Eigen::Vector3d bbx_min = {23.5, -99, -99}, bbx_max = {26.5, 99 , 99};


		std::vector<Eigen::Vector3d> forceFrame(500, Eigen::Vector3d::Zero());
#pragma omp parallel for num_threads(parameters.numOfThreads)
		for (int frame = 0 ; frame < 500; frame++)
		{		
			objMeshFormat testMesh;
			testMesh.readObjFile("E:/hydroStatic_object/Libuipc/libuipc/output/tests/sim_case/25_linear_arap_beam.cpp/scene_surface" + std::to_string(frame) + ".obj");
			testMesh.sepConnectedComponents();

			Eigen::Vector3d force = compute_contact_force(testMesh.componentsSep[0], testMesh.componentsSep[3], bbx_min, bbx_max, parameters);
			forceFrame[frame] = force;

			std::cout<<"Frame = "<<frame<<"; Force = ("<< " " << force[0] << " " << force[1] << " " << force[2]<<")" << std::endl;
			std::cout << std::endl;
		}






		std::ofstream outfile9("./output/force.txt", std::ios::trunc);
		for (int frame = 0; frame < 500; frame++)
		{
			outfile9 << std::scientific << std::setprecision(8) << frame << " " << forceFrame[frame][0] << " " << forceFrame[frame][1] << " " << forceFrame[frame][2] << std::endl;
		}
		outfile9.close();*/




		//double pi = 3.141592653;

		//Eigen::Matrix3d F = Eigen::Matrix3d::Identity();


		//// Rotate along x,y,z axis
		//Eigen::Matrix3d RX = Eigen::Matrix3d::Zero();
		//Eigen::Matrix3d RY = Eigen::Matrix3d::Zero();
		//Eigen::Matrix3d RZ = Eigen::Matrix3d::Zero();
		//double tx = 200, ty = 142, tz = 1450;
		//RX(0, 0) = 1.0;
		//RX(1, 1) = cos(tx);
		//RX(1, 2) = -sin(tx);
		//RX(2, 1) = sin(tx);
		//RX(2, 2) = cos(tx);

		//RY(1, 1) = 1.0;
		//RY(0, 0) = cos(ty);
		//RY(0, 2) = -sin(ty);
		//RY(2, 0) = sin(ty);
		//RY(2, 2) = cos(ty);

		//RZ(2, 2) = 1.0;
		//RZ(1, 1) = cos(tz);
		//RZ(0, 1) = -sin(tz);
		//RZ(1, 0) = sin(tz);
		//RZ(0, 0) = cos(tz);

		//F = F * RX * RY * RZ;


		//// Rotate along arbitrary axis
		//double theta = 458;
		//double cosTheta = std::cos(theta);
		//double sinTheta = std::sin(theta);
		//Eigen::Vector3d axis = {1.5,1.9,74.4};
		//axis.normalize();
		//Eigen::Matrix3d outerProduct = axis * axis.transpose();
		//Eigen::Matrix3d K;
		//K << 0, -axis.z(), axis.y(),
		//	 axis.z(), 0, -axis.x(),
		//	 -axis.y(), axis.x(), 0;
		//Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
		//// Rodrigues' rotation formula
		//Eigen::Matrix3d R = cosTheta * I + (1 - cosTheta) * outerProduct + sinTheta * K;

		//F = F * R;



		//Eigen::HouseholderQR<Eigen::MatrixXd> qr(F);
		//Eigen::MatrixXd QW = qr.householderQ();
		//Eigen::MatrixXd RW = qr.matrixQR().triangularView<Eigen::Upper>();
		//Eigen::MatrixXd A_reconstructed = QW * RW;
		//std::cout << "Matrix F:\n" << F << "\n" << std::endl;
		////std::cout << "Matrix Q:\n" << QW << "\n" << std::endl;
		////std::cout << "Matrix R:\n" << RW << "\n" << std::endl;
		////std::cout << "Q * R:\n" << A_reconstructed << "\n" << std::endl;



		//// Add stretch matrix
		//double s = 0.5;
		//Eigen::Matrix3d stretch = Eigen::Matrix3d::Zero();
		//stretch(0, 0) = s;
		//stretch(1, 1) = s;
		//stretch(2, 2) = s;
		//F = F * stretch;


		//F = Eigen::Matrix3d::Zero();
		//F(1, 1) = sqrt(2);
		//F(0, 0) = sqrt(2);

		//double xv2 = F(0,0) * F(0, 0) * PI / 2.0 * PI + F(0, 1) * F(0, 1) * PI / 2.0 * PI + F(0, 2) * F(0, 2) * PI * PI + 2 * F(0, 0) * F(0, 2) * 2.0 * PI ;
		//double yv2 = F(1,0) * F(1, 0) * PI / 2.0 * PI + F(1, 1) * F(1, 1) * PI / 2.0 * PI + F(1, 2) * F(1, 2) * PI * PI + 2 * F(1, 0) * F(1, 2) * 2.0 * PI ;
		//double zv2 = F(2,0) * F(2, 0) * PI / 2.0 * PI + F(2, 1) * F(2, 1) * PI / 2.0 * PI + F(2, 2) * F(2, 2) * PI * PI + 2 * F(2, 0) * F(2, 2) * 2.0 * PI ;

		//std::cout << "f(phi) = " << (xv2 + yv2 + zv2) / 2 / PI / PI   << std::endl;

		//// Test axis change_y
		//double xv2_2 = F(0, 0) * F(0, 0) * PI / 2.0 * PI + F(0, 2) * F(0, 2) * PI / 2.0 * PI + F(0, 1) * F(0, 1) * PI * PI + 2 * F(0, 0) * F(0, 1) * 2.0 * PI;
		//double yv2_2 = F(1, 0) * F(1, 0) * PI / 2.0 * PI + F(1, 2) * F(1, 2) * PI / 2.0 * PI + F(1, 1) * F(1, 1) * PI * PI + 2 * F(1, 0) * F(1, 1) * 2.0 * PI;
		//double zv2_2 = F(2, 0) * F(2, 0) * PI / 2.0 * PI + F(2, 2) * F(2, 2) * PI / 2.0 * PI + F(2, 1) * F(2, 1) * PI * PI + 2 * F(2, 0) * F(2, 1) * 2.0 * PI;

		//std::cout << "f_2(phi) = " << (xv2_2 + yv2_2 + zv2_2) / 2 / PI / PI  << std::endl;

		//// Test axis change_z
		//double xv2_3 = F(0, 1) * F(0, 1) * PI / 2.0 * PI + F(0, 2) * F(0, 2) * PI / 2.0 * PI + F(0, 0) * F(0, 0) * PI * PI / 2.0;
		//double yv2_3 = F(1, 1) * F(1, 1) * PI / 2.0 * PI + F(1, 2) * F(1, 2) * PI / 2.0 * PI + F(1, 0) * F(1, 0) * PI * PI / 2.0;
		//double zv2_3 = F(2, 1) * F(2, 1) * PI / 2.0 * PI + F(2, 2) * F(2, 2) * PI / 2.0 * PI + F(2, 0) * F(2, 0) * PI * PI / 2.0;

		//std::cout << "f_3(phi) = " << (xv2_3 + yv2_3 + zv2_3) / 2 / PI / PI  << std::endl;

		//std::cout << "F_00 * F_02 + F_10 * F_12 + F_20 * F_22 = " << F(0, 0) * F(0, 2) + F(1, 0) * F(1, 2) + F(2, 0) * F(2, 2) << std::endl;
	}
	else
	{
		// Case:
		// 0. Cube tower stress test for ABD of triMesh ABD implementation
		// 1. Bunny test to verify the correctness of mpm simulation
		int caseNum = 1;
		if (caseNum == 0)
		{

			double IPC_dis = 0.01;

			Material mat1;
			mat1.density = 4000;
			mat1.E = 7.26e12;
			mat1.updateDenpendecies();

			Material mat2;
			mat2.density = 80;
			mat2.E = 7.26e12;
			mat2.updateDenpendecies();



			std::vector<meshConfiguration> config;
			meshConfiguration m1;
			m1.filePath = "../input/cube_eq.obj";
			m1.mesh_material = mat1;
			m1.note = "cube_0";
			Eigen::Vector3d trans = { -3, 2.4, 2.7 };
			m1.velocity = { 5,0,0 };
			m1.translation = trans;
			m1.per_point_volume = 0.01;
			config.push_back(m1);




			int count = 1;
			for (int z = 0; z < 4; z++)
			{
				for (int x = 0; x < 4; x++)
				{
					for (int y = 0; y < 4; y++)
					{
						count += 1;
						m1.mesh_material = mat2;
						m1.note = "cube_" + std::to_string(count);
						Eigen::Vector3d trans = { (double)x * 1.03 , (double)y * 1.03 , (double)z * 1.03 + 2 };
						m1.translation = trans;
						m1.velocity = { 0,0,0 };
						config.push_back(m1);
					}
				}
			}



			triMesh triSimMesh;
			triSimMesh.createGlobalSimulationTriMesh_ABD(config);




			FEMParamters parameters;
			parameters.gravity = { 0, 0, -9.8};
			parameters.num_timesteps = 10000;
			parameters.numOfThreads = 20;
			parameters.dt = 1.0e-3;
			parameters.outputFrequency = 20;
			parameters.simulation_Mode = "ABD"; // Normal, ABD, Coupling
			parameters.enableGround = true;
			parameters.searchResidual = 0.01;
			parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
			parameters.rigidMode = true;
			parameters.IPC_dis = IPC_dis;
			parameters.IPC_eta = 0.05;
			parameters.IPC_kStiffness = 1.0e12;
			parameters.IPC_B3Stiffness = 500;
			parameters.ABD_Coeff = 1.0e12;



			implicitFEM_ABD_triMesh(triSimMesh, parameters);

			

		}
		else if (caseNum == 1)
		{
			double IPC_dis = 0.01;

			Material mat1;
			mat1.density = 3000;
			mat1.E = 3.26e12;
			mat1.thetaF = 8.0e9;
			mat1.fracture_start_force = 1.0E1;
			mat1.updateDenpendecies();



			std::vector<meshConfiguration> config;
			meshConfiguration m1;
			m1.filePath = "../input/bunny_highRes.obj";
			m1.mesh_material = mat1;
			m1.note = "bunny1";
			m1.breakable = true;
			m1.translation = { 0, 0, 6 };
			m1.rotation_angle = {0,20,45};
			m1.velocity = { 0,10,0 };
			//m1.velocity = { 0,0,-10 };
			m1.per_point_volume = 0.01;
			config.push_back(m1);

			meshConfiguration m2;
			m2.filePath = "../input/bunny_highRes.obj";
			m2.mesh_material = mat1;
			m2.note = "bunny2";
			m2.translation = { 1, 4, 4 };
			m2.velocity = {0,0,0};
			//m2.velocity = {0,0,-10};
			m2.breakable = false;
			m2.per_point_volume = 0.01;
			config.push_back(m2);






			triMesh triSimMesh;
			triSimMesh.createGlobalSimulationTriMesh_ABD(config);




			FEMParamters parameters;
			parameters.gravity = { 0, 0, -9.8 };
			parameters.num_timesteps = 10000;
			parameters.numOfThreads = 20;
			parameters.dt = 1.0e-2;
			parameters.outputFrequency = 1;
			parameters.simulation_Mode = "ABD"; // Normal, ABD, Coupling
			parameters.enableGround = true;
			parameters.searchResidual = 0.1;
			parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
			parameters.rigidMode = true;
			parameters.IPC_dis = IPC_dis;
			parameters.IPC_eta = 0.05;
			parameters.IPC_kStiffness = 1.0e12;
			parameters.IPC_B3Stiffness = 500;
			parameters.ABD_Coeff = 1.0e12;



			implicitFEM_ABD_triMesh(triSimMesh, parameters);

		}

	}


	return 0;

}
