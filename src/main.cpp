#include "simulator.h" 
#include "tools.h"



// TODO
// 1. Accelerate and parallel openMP
// 2. Decimate the generated mesh for performance
// 3. triSimMesh.boundaryCondition_node_surface (& boundaryCondition_node_interior) is not properly handled






int main()
{

	if (0)
	{
	
		std::string bunny = "../input/bunny_highRes.obj";
		std::string cube = "../input/bunny_highRes.obj";


		objMeshFormat bunnyMesh;
		bunnyMesh.readObjFile(bunny, true);


		objMeshFormat cubeMesh;
		Eigen::Vector3d trans = {0,1.0,0};
		cubeMesh.readObjFile(cube, true, Eigen::Affine3d::Identity(), Eigen::Vector3d::Ones(), trans);

		objMeshFormat diff = bunnyMesh.boolean_difference_with_mesh(cubeMesh);

		diff.outputFile("../output/bunny_cube_diff.obj", -99, false);
		










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
		// 0: Three-point-bending beam test
		// 1. Two tetrahedrals collision test to verify the correctness of the solver
		// 2. Single tetrahedral test
		// 3. ABD tes: cube collision
		// 4. Explicit test
		// 5. Object hanging test to verify element sampling algorithm
		// 6. Two tetrahedrals collision test to verify the correctness of the solver in ABD 
		// 7. Cube tower stress test for ABD 
		// 8. Cube tower stress test for ABD of triMesh ABD implementation
		// 9. Bunny test to verify the correctness of mpm simulation
		int caseNum = 9;
		if (caseNum == 0)
		{
			//Material mat1;
			//mat1.density = 7880;
			//mat1.E = 7.26e10;
			//mat1.updateDenpendecies();


			//
			//std::vector<meshConfiguration> config;			
			//meshConfiguration m1, m2, m3, m4;
			//m1.filePath = "../input/beam.msh";
			//m1.mesh_material = mat1;
			//m1.scale = {1000, 1000, 1000};
			//m1.note = "beam";
			//config.push_back(m1);

			//m2 = m1;
			//m2.filePath = "../input/left_support.msh";
			//m2.note = "left_support";
			//m2.translation = { 0, 0, 1.0e-4};
			////config.push_back(m2);

			//m3 = m1;
			//m3.filePath = "../input/middle_support.msh";
			//m3.note = "middle_support";
			//m3.translation = { 0, 0, -1.0e-4 };
			////config.push_back(m3);

			//m4 = m1;
			//m4.filePath = "../input/impactor.msh";
			//m4.note = "impactor";
			//m4.translation = { 0, 0, 1.5e-2 };
			//config.push_back(m4);




			//Mesh tetSimMesh;
			//for (int i = 0; i < config.size(); i++)
			//{
			//	tetMesh msh_tmp;
			//	msh_tmp.readMesh(config[i]);
			//	msh_tmp.initializeTetMesh();
			//	tetSimMesh.objectsTetMesh[msh_tmp.tetMeshNote] = msh_tmp;
			//}
			//tetSimMesh.createGlobalSimulationMesh();



			//



			//std::cout << "tetSimMesh.pos_node.size() = " << tetSimMesh.pos_node.size() << std::endl;
			//std::cout << "tetSimMesh.tetrahedrals.size() = " << tetSimMesh.tetrahedrals.size() << std::endl;


			//FEMParamters parameters;
			//parameters.gravity = { 0, 0, 0 };
			//parameters.num_timesteps = 500;
			//parameters.numOfThreads = 12;
			//parameters.dt = 5.0e-4;
			//parameters.outputFrequency = 20;
			//parameters.enableGround = false;
			//parameters.searchResidual = 5.0;
			//parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
			//parameters.rigidMode = true;
			////parameters.objectNames = objectNames;
			//parameters.IPC_dis = 0.01;
			//parameters.IPC_eta = 0.05;
			//parameters.IPC_kStiffness = 1.0e14;
			//parameters.IPC_hashSize = tetSimMesh.calLargestEdgeLength() * 1.1;
			//parameters.IPC_B3Stiffness = 500;



			//for (int p = 0; p < tetSimMesh.pos_node.size(); p++)
			//{
			//	if (tetSimMesh.note_node[p] == "impactor")
			//	{
			//		tetSimMesh.boundaryCondition_node[p].type = 1;
			//		for (int fra = 0; fra < parameters.num_timesteps; fra++)
			//		{
			//			double incre = 0.1 / (double)parameters.num_timesteps * (double)fra;
			//			Eigen::Vector3d inc = {0,0,-incre };
			//			Eigen::Vector3d desiredPos = inc + tetSimMesh.pos_node[p];

			//			tetSimMesh.boundaryCondition_node[p].location.push_back(desiredPos);
			//		}
			//	}

			//	if (tetSimMesh.note_node[p] == "beam")
			//	{
			//		if (tetSimMesh.pos_node[p][0] <= -22)
			//		{
			//			tetSimMesh.boundaryCondition_node[p].type = 1;
			//			for (int fra = 0; fra < parameters.num_timesteps; fra++)
			//			{
			//				Eigen::Vector3d desiredPos = tetSimMesh.pos_node[p];
			//				tetSimMesh.boundaryCondition_node[p].location.push_back(desiredPos);
			//			}
			//		}
			//	}

			//}



			//implicitFEM(tetSimMesh, parameters);


		}
		else if (caseNum == 1)
		{
			//Material mat1;
			//mat1.density = 1000;
			//mat1.E = 7.26e5;
			//mat1.updateDenpendecies();



			//std::vector<meshConfiguration> config;
			//meshConfiguration m1,m2;
			//m1.filePath = "../input/tet_origin.msh";
			//m1.mesh_material = mat1;
			//m1.translation = {0,0,1.1};
			//m1.velocity = {0,0,1};
			//m1.note = "tet_origin";
			//config.push_back(m1);

			//m2.filePath = "../input/tet_origin.msh";
			//m2.mesh_material = mat1;
			//m2.translation = { 0,0,4.3 };
			//m2.velocity = { 0,0,-1 };
			//m2.note = "tet_origin_top";
			//config.push_back(m2);



			//Mesh tetSimMesh;
			//for (int i = 0; i < config.size(); i++)
			//{
			//	tetMesh msh_tmp;
			//	msh_tmp.readMesh(config[i]);
			//	msh_tmp.initializeTetMesh();
			//	tetSimMesh.objectsTetMesh[msh_tmp.tetMeshNote] = msh_tmp;
			//}
			//tetSimMesh.createGlobalSimulationMesh();


			//std::cout << "tetSimMesh.pos_node.size() = " << tetSimMesh.pos_node.size() << std::endl;
			//std::cout << "tetSimMesh.tetrahedrals.size() = " << tetSimMesh.tetrahedrals.size() << std::endl;



			//FEMParamters parameters;
			//parameters.gravity = { 0, 0, 0 };
			//parameters.num_timesteps = 3;
			//parameters.numOfThreads = 1;
			//parameters.dt = 1.0e-2;
			//parameters.outputFrequency = 10;
			//parameters.enableGround = true;
			//parameters.searchResidual = 0.001;
			//parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
			//parameters.rigidMode = true;
			////parameters.objectNames = objectNames;
			//parameters.IPC_dis = 0.001;
			//parameters.IPC_eta = 0.01;
			//parameters.IPC_kStiffness = 1.0e12;
			//parameters.IPC_hashSize = tetSimMesh.calLargestEdgeLength() * 1.3;
			//parameters.IPC_B3Stiffness = 500;




			//implicitFEM(tetSimMesh, parameters);


		}
		else if (caseNum == 2)
		{
			//Material mat1;
			//mat1.density = 5000;
			//mat1.E = 7.26e10;
			//mat1.updateDenpendecies();



			//std::vector<meshConfiguration> config;
			//meshConfiguration m1, m2;
			//m1.filePath = "../input/tet_origin.msh";
			//m1.mesh_material = mat1;
			//m1.translation = { 0,0,1.1 };
			//m1.velocity = { 0,0,0 };
			//m1.note = "tet_origin";
			//config.push_back(m1);




			//Mesh tetSimMesh;
			//for (int i = 0; i < config.size(); i++)
			//{
			//	tetMesh msh_tmp;
			//	msh_tmp.readMesh(config[i]);
			//	msh_tmp.initializeTetMesh();
			//	tetSimMesh.objectsTetMesh[msh_tmp.tetMeshNote] = msh_tmp;
			//}
			//tetSimMesh.createGlobalSimulationMesh();


			//std::cout << "tetSimMesh.pos_node.size() = " << tetSimMesh.pos_node.size() << std::endl;
			//std::cout << "tetSimMesh.tetrahedrals.size() = " << tetSimMesh.tetrahedrals.size() << std::endl;



			//FEMParamters parameters;
			//parameters.gravity = { 0, 0, 0 };
			//parameters.num_timesteps = 50000;
			//parameters.numOfThreads = 1;
			//parameters.dt = 1.0e-4;
			//parameters.outputFrequency = 1000;
			//parameters.enableGround = true;
			//parameters.searchResidual = 0.001;
			//parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
			//parameters.rigidMode = true;
			////parameters.objectNames = objectNames;
			//parameters.IPC_dis = 0.001;
			//parameters.IPC_eta = 0.05;
			//parameters.IPC_kStiffness = 1.0e12;
			//parameters.IPC_hashSize = tetSimMesh.calLargestEdgeLength() * 1.3;
			//parameters.IPC_B3Stiffness = 500;


			////tetMesh.pos_node[3] = {0,0.05,1.15};
			////tetMesh.update_F(1);


			////// Perform SVD on the matrix A
			////Eigen::JacobiSVD<Eigen::MatrixXd> svd(tetMesh.tetra_F[0], Eigen::ComputeFullU | Eigen::ComputeFullV);
			////Eigen::MatrixXd U = svd.matrixU();
			////Eigen::VectorXd S = svd.singularValues();
			////Eigen::MatrixXd V = svd.matrixV();

			////// Compute R and Q
			////Eigen::MatrixXd R = U * V.transpose();
			////Eigen::MatrixXd Q = V * S.asDiagonal() * V.transpose();

			////std::cout << "Matrix F:\n" << tetMesh.tetra_F[0] << "\n\n";
			////std::cout << "Rotation Matrix R:\n" << R << "\n\n";
			////std::cout << "Stretch Matrix Q:\n" << Q << "\n\n";
			////


			////Eigen::Matrix<double, 9, 12> DFDX = ElasticEnergy::dF_wrt_dx(tetMesh.tetra_DM_inv[0]);
			////Vector9d grad_Energy_wrt_F = flatenMatrix3d(R);
			////Eigen::Matrix<double, 12, 1> engGrad = DFDX.transpose() * grad_Energy_wrt_F;

			////std::cout << "Force_1 = (" << engGrad(0, 0) << "," << engGrad(1, 0) << "," << engGrad(2, 0) << ")" << std::endl;
			////std::cout << "Force_2 = (" << engGrad(3, 0) << "," << engGrad(4, 0) << "," << engGrad(5, 0) << ")" << std::endl;
			////std::cout << "Force_3 = (" << engGrad(6, 0) << "," << engGrad(7, 0) << "," << engGrad(8, 0) << ")" << std::endl;
			////std::cout << "Force_4 = (" << engGrad(9, 0) << "," << engGrad(10, 0) << "," << engGrad(11, 0) << ")" << std::endl;

			////

			//// Explicit FEM
			//{
			//	std::vector<Eigen::Vector3d> nodeForce;

			//	// Simulation loop
			//	for (int step = 0; step < parameters.num_timesteps; ++step) 
			//	{
			//		if (step % 100 == 0)
			//		{
			//			std::cout << "Timestep = " << step << std::endl;
			//			tetSimMesh.exportSurfaceMesh("surfMesh", step);
			//		}


			//		tetSimMesh.elastForce_node.resize(tetSimMesh.pos_node.size());
			//		// Reset forces
			//		for (int nf = 0; nf < tetSimMesh.elastForce_node.size(); nf++)
			//		{
			//			tetSimMesh.elastForce_node[nf] = Eigen::Vector3d::Zero();
			//		}
			//		tetSimMesh.elastForce_node[0] = {0,135,100};
			//		tetSimMesh.elastForce_node[1] = {145,47,28};
			//		tetSimMesh.elastForce_node[2] = {13,167,95};
			//		tetSimMesh.elastForce_node[3] = {0,0,100};

			//		// Apply gravity
			//		for (int nf = 0; nf < tetSimMesh.elastForce_node.size(); nf++)
			//		{
			//			tetSimMesh.elastForce_node[nf] += tetSimMesh.mass_node[nf] * parameters.gravity;
			//		}

			//		// Compute forces from elements
			//		for (int nf = 0; nf < tetSimMesh.tetrahedrals.size(); nf++)
			//		{
			//			// Compute deformation gradient F
			//			tetSimMesh.update_F(1);


			//			// Compute First Piola-Kirchhoff stress tensor P
			//			Eigen::JacobiSVD<Eigen::MatrixXd> svd(tetSimMesh.tetra_F[nf], Eigen::ComputeFullU | Eigen::ComputeFullV);
			//			Eigen::MatrixXd U = svd.matrixU();
			//			Eigen::VectorXd S = svd.singularValues();
			//			Eigen::MatrixXd V = svd.matrixV();
			//			Eigen::MatrixXd R = U * V.transpose();
			//			Eigen::MatrixXd Q = V * S.asDiagonal() * V.transpose();
			//			
			//			Eigen::Matrix3d P = R;


			//			if (step % 1000 == 0)
			//			{
			//				std::cout << "Rotation Matrix R:\n" << R << "\n\n";
			//				std::cout << "Stretch Matrix S:\n" << Q << "\n\n";
			//			}


			//			// Compute forces on nodes
			//			Eigen::Matrix3d H = -tetSimMesh.tetra_vol[nf] * P * tetSimMesh.tetra_DM_inv[nf].transpose();

			//			// Distribute forces to nodes
			//			Eigen::Vector3d f0 = H.col(0);
			//			Eigen::Vector3d f1 = H.col(1);
			//			Eigen::Vector3d f2 = H.col(2);
			//			Eigen::Vector3d f3 = -f0 - f1 - f2;


			//			Eigen::Vector4i tet = tetSimMesh.tetrahedrals[nf];
			//			tetSimMesh.elastForce_node[tet[0]] += f3;
			//			tetSimMesh.elastForce_node[tet[1]] += f0;
			//			tetSimMesh.elastForce_node[tet[2]] += f1;
			//			tetSimMesh.elastForce_node[tet[3]] += f2;
			//		}

			//		// Time integration (Explicit Euler)
			//		for (int nf = 0; nf < tetSimMesh.elastForce_node.size(); nf++)
			//		{
			//			// Update velocity and position
			//			Eigen::Vector3d acceleration = tetSimMesh.elastForce_node[nf] / tetSimMesh.mass_node[nf];
			//			tetSimMesh.vel_node[nf] += acceleration * parameters.dt;
			//			tetSimMesh.pos_node[nf] += tetSimMesh.vel_node[nf] * parameters.dt;
			//		}
			//		
			//	}


			//}



		}
		else if (caseNum == 3)
		{
			//Material mat1;
			//mat1.density = 2780;
			//mat1.E = 7.26e12;
			//mat1.updateDenpendecies();



			//std::vector<meshConfiguration> config;
			//meshConfiguration m1, m2;
			//m1.filePath = "../input/cube.msh";
			//m1.mesh_material = mat1;
			////m1.rotation_angle = {-90, -19, 1.7};
			//m1.translation = {0, 0, 2.0};
			//m1.note = "cube1";
			//config.push_back(m1);


			//m2.filePath = "../input/cube.msh";
			//m2.mesh_material = mat1;
			//m2.note = "cube2";
			//m2.translation = { 1.5, 0, 2.0 };
			//config.push_back(m2);




			//Mesh_ABD tetSimMesh;
			//for (int i = 0; i < config.size(); i++)
			//{
			//	tetMesh msh_tmp;
			//	msh_tmp.readMesh(config[i]);
			//	msh_tmp.initializeTetMesh();
			//	tetSimMesh.objectsTetMesh[msh_tmp.tetMeshNote] = msh_tmp;
			//}
			//tetSimMesh.createGlobalSimulationMesh_ABD();
			//tetSimMesh.translation_vel_ABD[0] = { 1,0,-3 };
			//tetSimMesh.translation_vel_ABD[1] = { 0,0,-3 };


			//std::cout << "tetSimMesh.pos_node.size() = " << tetSimMesh.pos_node.size() << std::endl;
			//std::cout << "tetSimMesh.tetrahedrals.size() = " << tetSimMesh.tetrahedrals.size() << std::endl;



			//FEMParamters parameters;
			//parameters.gravity = { 0, 0, 0 };
			//parameters.num_timesteps = 10000;
			//parameters.numOfThreads = 12;
			//parameters.dt = 1.0e-3;
			//parameters.outputFrequency = 10;
			//parameters.simulation_Mode = "ABD"; // Normal, ABD, Coupling
			//parameters.enableGround = true;
			//parameters.searchResidual = 6.0;
			//parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
			//parameters.rigidMode = true;
			//parameters.IPC_dis = 0.01;
			//parameters.IPC_eta = 0.05;
			//parameters.IPC_kStiffness = 1.0e12;
			//parameters.IPC_hashSize = tetSimMesh.calLargestEdgeLength() * 1.1;
			//parameters.IPC_B3Stiffness = 500;
			//parameters.ABD_Coeff = 1.0e9;



			//implicitFEM_ABD(tetSimMesh, parameters);

		}
		else if (caseNum == 4)
		{
			//Material mat1;
			//mat1.density = 2700;
			//mat1.E = 1.0e12;
			//mat1.nu = 0.3;

			//mat1.thetaf = 6.0E5;
			//mat1.Gf = 300;
			//double hsTarget = 0.45;
			////mat1.lch = hsTarget / (1 + hsTarget) / mat1.calHsBar() * 2;
			//mat1.lch = 1.0;
			//mat1.updateDenpendecies();
			//std::cout << "material1.lch = " << mat1.lch << std::endl;
			//std::cout << "material1.Hs = " << mat1.Hs << std::endl;




			//std::vector<meshConfiguration> config;
			//meshConfiguration m1, m2, m3, m4;
			//m1.filePath = "../input/impactor.msh";
			//m1.mesh_material = mat1;
			//m1.scale = { 1000, 1000, 1000 };
			//m1.note = "impactor";
			//config.push_back(m1);



			//Mesh tetSimMesh;
			//for (int i = 0; i < config.size(); i++)
			//{
			//	tetMesh msh_tmp;
			//	msh_tmp.readMesh(config[i]);
			//	msh_tmp.initializeTetMesh();
			//	tetSimMesh.objectsTetMesh[msh_tmp.tetMeshNote] = msh_tmp;
			//}
			//tetSimMesh.createGlobalSimulationMesh();


			//std::cout << "tetSimMesh.pos_node.size() = " << tetSimMesh.pos_node.size() << std::endl;
			//std::cout << "tetSimMesh.tetrahedrals.size() = " << tetSimMesh.tetrahedrals.size() << std::endl;


			//FEMParamters parameters;
			//parameters.gravity = { 0, 0, 0 };
			//parameters.num_timesteps = 100000;
			//parameters.numOfThreads = 12;
			//parameters.dt = 1.0e-6;
			//parameters.outputFrequency = 100;
			//parameters.enableGround = false;
			//parameters.searchResidual = 5.0;
			//parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
			//parameters.rigidMode = true;
			//parameters.IPC_dis = 0.01;
			//parameters.IPC_eta = 0.05;
			//parameters.IPC_kStiffness = 1.0e14;
			//parameters.IPC_hashSize = tetSimMesh.calLargestEdgeLength() * 1.1;
			//parameters.IPC_B3Stiffness = 500;



			//for (int p = 0; p < tetSimMesh.pos_node.size(); p++)
			//{
			//	if (tetSimMesh.note_node[p] == "impactor")
			//	{
			//		if (tetSimMesh.pos_node[p][1] >= 9.0)
			//		{
			//			tetSimMesh.boundaryCondition_node[p].type = 2;
			//			for (int fra = 0; fra < parameters.num_timesteps; fra++)
			//			{
			//				Eigen::Vector3d inc = { 0,1000,0 };
			//				tetSimMesh.boundaryCondition_node[p].force.push_back(inc);
			//			}
			//		}

			//		/*if (tetSimMesh.pos_node[p][1] < -9.95)
			//		{
			//			tetSimMesh.boundaryCondition_node[p].type = 2;
			//			for (int fra = 0; fra < parameters.num_timesteps; fra++)
			//			{
			//				Eigen::Vector3d inc = { 0,-500,0 };
			//				tetSimMesh.boundaryCondition_node[p].force.push_back(inc);
			//			}
			//		}*/

			//	}

			//	if (std::abs(tetSimMesh.elastForce_node[p].norm() - 0) > 0.0001)
			//	{
			//		//std::cout << "fsdjf" << std::endl;
			//	}

			//}



			//explicitFEM(tetSimMesh, parameters);


		}
		else if (caseNum == 5)
		{
			//Material mat1;
			//mat1.density = 2700;
			//mat1.E = 1.0e7;
			//mat1.nu = 0.3;

			//mat1.thetaf = 6.0E5;
			//mat1.Gf = 300;
			//double hsTarget = 0.45;
			//mat1.lch = 1.0;
			//mat1.updateDenpendecies();



			//std::vector<meshConfiguration> config;
			//meshConfiguration m1, m2, m3, m4;
			//m1.filePath = "../input/bunny_highRes_mesh.msh";
			//m1.mesh_material = mat1;
			//m1.note = "bunny";
			////m1.scale = {1000,1000,1000};
			////m1.rotation_angle = {PI_Value / 2.0,0,0};
			//m1.translation = { 0,0,0 };
			//config.push_back(m1);



			////std::vector<meshConfiguration> config;
			////meshConfiguration m1, m2, m3, m4;
			////m1.filePath = "../input/Right_top_move.msh";
			////m1.mesh_material = mat1;
			////m1.note = "cube";
			////m1.scale = { 1000,1000,1000 };
			////m1.rotation_angle = { PI_Value / 2.0,0,0 };
			////m1.translation = { -48,4,5.3 };
			////config.push_back(m1);



			//Mesh tetSimMesh;
			//for (int i = 0; i < config.size(); i++)
			//{
			//	tetMesh msh_tmp;
			//	msh_tmp.readMesh(config[i]);
			//	msh_tmp.initializeTetMesh();
			//	tetSimMesh.objectsTetMesh[msh_tmp.tetMeshNote] = msh_tmp;
			//}
			//tetSimMesh.createGlobalSimulationMesh();

			//tetSimMesh.exportEdges("bunnyEdges");




			//FEMParamters parameters;
			//parameters.gravity = { 0, 0, -9.8 };
			//parameters.num_timesteps = 1000000;
			//parameters.numOfThreads = 16;
			//parameters.dt = 1.0e-4;
			//parameters.outputFrequency = 100;
			//parameters.enableGround = false;
			//parameters.searchResidual = 7.0;
			//parameters.model = "ARAP"; // neoHookean ARAP ARAP_linear ACAP
			//parameters.rigidMode = true;
			//parameters.IPC_dis = 0.01;
			//parameters.IPC_eta = 0.05;
			//parameters.IPC_kStiffness = 1.0e12;
			//parameters.IPC_hashSize = tetSimMesh.calLargestEdgeLength() * 1.1;
			//parameters.IPC_B3Stiffness = 500;

			//parameters.MLS_radius = 0.5;

			//std::cout << "parameters.IPC_hashSize = " << parameters.IPC_hashSize << std::endl;




			///*{
			//	objMeshFormat crack;
			//	Eigen::Vector3d p1 = { -4,-4,15 };
			//	Eigen::Vector3d p2 = { 4,-4,15 };
			//	Eigen::Vector3d p3 = { 4,4,15 };
			//	Eigen::Vector3d p4 = { -4,4,15 };
			//	crack.vertices.push_back(p1);
			//	crack.vertices.push_back(p2);
			//	crack.vertices.push_back(p3);
			//	crack.vertices.push_back(p4);
			//	std::vector<int> fc;
			//	fc.push_back(0);
			//	fc.push_back(1);
			//	fc.push_back(2);
			//	fc.push_back(3);
			//	crack.vertFaces.push_back(fc);

			//	crack.outputFile("crack", -99);

			//	tetSimMesh.sample_MLS_points(crack, parameters.MLS_num_MLS_Pts, parameters.MLS_radius, parameters.numOfThreads);

			//}*/


			//{
			//	//objMeshFormat crack;	
			//	//crack.readObjFile("../input/bunnyCrack.obj", true);
			//	////crack.readObjFile("../input/bunnyCrack_Plane.obj", true);
			//	//crack.outputFile("bunnyCrack", -99, true);
			//	//tetSimMesh.sample_MLS_points(crack, parameters.MLS_num_MLS_Pts, parameters.MLS_radius, parameters.numOfThreads);
			//}
			//


			//std::cout << "tetSimMesh.pos_node.size() = " << tetSimMesh.pos_node.size() << std::endl;
			//std::cout << "tetSimMesh.tetrahedrals.size() = " << tetSimMesh.tetrahedrals.size() << std::endl;




			//for (int p = 0; p < tetSimMesh.pos_node.size(); p++)
			//{
			//	/*if (tetSimMesh.note_node[p] == "cube")
			//	{
			//		if (tetSimMesh.pos_node[p][2] >= 25)
			//		{
			//			tetSimMesh.boundaryCondition_node[p].type = 1;
			//			for (int fra = 0; fra < parameters.num_timesteps; fra++)
			//			{
			//				tetSimMesh.boundaryCondition_node[p].location.push_back(tetSimMesh.pos_node[p]);
			//			}
			//		}

			//		if (tetSimMesh.pos_node[p][2] <= 6)
			//		{
			//			tetSimMesh.boundaryCondition_node[p].type = 1;
			//			Eigen::Vector3d vel = { 0,0,-1.0 };
			//			for (int fra = 0; fra < parameters.num_timesteps; fra++)
			//			{				
			//				tetSimMesh.boundaryCondition_node[p].location.push_back(tetSimMesh.pos_node[p] + vel * parameters.dt * (double)fra);
			//			}
			//		}


			//	}*/


			//}


			//{
			//	//// export all tetrahedron that has MLS points
			//	//{
			//	//	std::ofstream outfile9("./output/tetMLS.obj", std::ios::trunc);
			//	//	for (std::map<int, std::vector<MLSPoints>>::iterator it = tetSimMesh.MLSPoints_tet_map.begin(); it != tetSimMesh.MLSPoints_tet_map.end(); it++)
			//	//	{
			//	//		int indTet = it->first;
			//	//		for (int i = 0; i < 4; i++)
			//	//		{
			//	//			int indVt = tetSimMesh.tetrahedrals[indTet][i];
			//	//			Eigen::Vector3d pos = tetSimMesh.pos_node[indVt];
			//	//			outfile9 << std::scientific << std::setprecision(8) << "v " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
			//	//		}
			//	//	}
			//	//	for (int j = 0; j < tetSimMesh.MLSPoints_tet_map.size(); j++)
			//	//	{
			//	//		outfile9 << std::scientific << std::setprecision(8) << "f " << j * 4 + 1 << " " << j * 4 + 2 << " " << j * 4 + 3 << std::endl;
			//	//		outfile9 << std::scientific << std::setprecision(8) << "f " << j * 4 + 2 << " " << j * 4 + 3 << " " << j * 4 + 4 << std::endl;
			//	//		outfile9 << std::scientific << std::setprecision(8) << "f " << j * 4 + 3 << " " << j * 4 + 4 << " " << j * 4 + 1 << std::endl;
			//	//		outfile9 << std::scientific << std::setprecision(8) << "f " << j * 4 + 4 << " " << j * 4 + 1 << " " << j * 4 + 2 << std::endl;
			//	//	}
			//	//	outfile9.close();
			//	//}

			//	//// export tetrahedron-74
			//	//{
			//	//	std::ofstream outfile9("./output/tetMLS-74.obj", std::ios::trunc);
			//	//	for (int i = 0; i < 4; i++)
			//	//	{
			//	//		int indVt = tetSimMesh.tetrahedrals[74][i];
			//	//		Eigen::Vector3d pos = tetSimMesh.pos_node[indVt];
			//	//		outfile9 << std::scientific << std::setprecision(8) << "v " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
			//	//	}
			//	//	outfile9 << std::scientific << std::setprecision(8) << "f " << 1 << " " << 2 << " " << 3 << std::endl;
			//	//	outfile9 << std::scientific << std::setprecision(8) << "f " << 2 << " " << 3 << " " << 4 << std::endl;
			//	//	outfile9 << std::scientific << std::setprecision(8) << "f " << 3 << " " << 4 << " " << 1 << std::endl;
			//	//	outfile9 << std::scientific << std::setprecision(8) << "f " << 4 << " " << 1 << " " << 2 << std::endl;
			//	//	outfile9.close();
			//	//}

			//	//// export tetrahedron-74's neighbouring tets
			//	//{
			//	//	std::ofstream outfile9("./output/tetMLS-74-neigs.obj", std::ios::trunc);
			//	//	std::set<int> allTets;
			//	//	for (int i = 0; i < 4; i++)
			//	//	{
			//	//		int indVt = tetSimMesh.tetrahedrals[74][i];
			//	//		std::vector<int> tetsThisNode = tetSimMesh.tetrahedrals_node[indVt];
			//	//		for (int j = 0; j < tetsThisNode.size(); j++)
			//	//		{
			//	//			int indTet = tetsThisNode[j];
			//	//			allTets.insert(indTet);
			//	//		}
			//	//	}


			//	//	for (std::set<int>::iterator it = allTets.begin(); it != allTets.end(); it++)
			//	//	{
			//	//		int indTet = *it;
			//	//		for (int i = 0; i < 4; i++)
			//	//		{
			//	//			int indVt = tetSimMesh.tetrahedrals[indTet][i];
			//	//			Eigen::Vector3d pos = tetSimMesh.pos_node[indVt];
			//	//			outfile9 << std::scientific << std::setprecision(8) << "v " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
			//	//		}
			//	//	}
			//	//	for (int j = 0; j < allTets.size(); j++)
			//	//	{
			//	//		outfile9 << std::scientific << std::setprecision(8) << "f " << j * 4 + 1 << " " << j * 4 + 2 << " " << j * 4 + 3 << std::endl;
			//	//		outfile9 << std::scientific << std::setprecision(8) << "f " << j * 4 + 2 << " " << j * 4 + 3 << " " << j * 4 + 4 << std::endl;
			//	//		outfile9 << std::scientific << std::setprecision(8) << "f " << j * 4 + 3 << " " << j * 4 + 4 << " " << j * 4 + 1 << std::endl;
			//	//		outfile9 << std::scientific << std::setprecision(8) << "f " << j * 4 + 4 << " " << j * 4 + 1 << " " << j * 4 + 2 << std::endl;
			//	//	}

			//	//	outfile9.close();
			//	//}

			//	//// export tetrahedron-74 MLS-0
			//	//{
			//	//	std::ofstream outfile9("./output/tetMLS-74-MLS-0.obj", std::ios::trunc);
			//	//	Eigen::Vector3d pos = tetSimMesh.MLSPoints_tet_map[74][0].pos_Rest;
			//	//	outfile9 << std::scientific << std::setprecision(8) << "v " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
			//	//	outfile9.close();
			//	//}

			//	//// export tetrahedron-74 MLS-0 control points
			//	//{
			//	//	std::ofstream outfile9("./output/tetMLS-74-MLS-0-ctrPts.obj", std::ios::trunc);
			//	//	for (auto value : tetSimMesh.MLSPoints_tet_map[74][0].index_node)
			//	//	{
			//	//		Eigen::Vector3d pos = tetSimMesh.pos_node[value];
			//	//		outfile9 << std::scientific << std::setprecision(8) << "v " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;						
			//	//	}
			//	//	outfile9.close();

			//	//}



			//	//// export all MLS points that have insufficient control points
			//	//{
			//	//	std::ofstream outfile9("./output/tetMLS_insuff.obj", std::ios::trunc);
			//	//	for (std::map<int, std::vector<MLSPoints>>::iterator it = tetSimMesh.MLSPoints_tet_map.begin(); it != tetSimMesh.MLSPoints_tet_map.end(); it++)
			//	//	{
			//	//		int indTet = it->first;
			//	//		for (int m = 0; m < tetSimMesh.MLSPoints_tet_map[indTet].size(); m++)
			//	//		{
			//	//			//int indVt = tetSimMesh.tetrahedrals[indTet][m];
			//	//			if (tetSimMesh.MLSPoints_tet_map[indTet][m].index_node.size() < 4)
			//	//			{
			//	//				Eigen::Vector3d pos = tetSimMesh.MLSPoints_tet_map[indTet][m].pos_Rest;
			//	//				outfile9 << std::scientific << std::setprecision(8) << "v " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
			//	//			}					
			//	//		}
			//	//	}
			//	//	outfile9.close();
			//	//}



			//}

			//
			//


			//implicitFEM(tetSimMesh, parameters);
		}
		else if (caseNum == 6)
		{
			//Material mat1;
			//mat1.density = 1000;
			//mat1.E = 7.26e5;
			//mat1.updateDenpendecies();



			//std::vector<meshConfiguration> config;
			//meshConfiguration m1,m2;
			//m1.filePath = "../input/tet_origin.msh";
			//m1.mesh_material = mat1;
			//m1.translation = {0,0,1.1};
			////m1.velocity = {0,0,1};
			//m1.note = "tet_origin";
			//config.push_back(m1);

			//m2.filePath = "../input/tet_origin.msh";
			//m2.mesh_material = mat1;
			//m2.translation = { 0,0,4.3 };
			////m2.velocity = { 0,0,-1 };
			//m2.note = "tet_origin_top";
			//config.push_back(m2);



			//Mesh_ABD tetSimMesh;
			//for (int i = 0; i < config.size(); i++)
			//{
			//	tetMesh msh_tmp;
			//	msh_tmp.readMesh(config[i]);
			//	msh_tmp.initializeTetMesh();
			//	tetSimMesh.objectsTetMesh[msh_tmp.tetMeshNote] = msh_tmp;
			//}
			//tetSimMesh.createGlobalSimulationMesh_ABD();
			//tetSimMesh.translation_vel_ABD[0] = { 0,0,1 };
			//tetSimMesh.translation_vel_ABD[1] = { 0,0,-1 };


			//std::cout << "tetSimMesh.pos_node.size() = " << tetSimMesh.pos_node.size() << std::endl;
			//std::cout << "tetSimMesh.tetrahedrals.size() = " << tetSimMesh.tetrahedrals.size() << std::endl;



			//FEMParamters parameters;
			//parameters.gravity = { 0, 0, 0 };
			//parameters.num_timesteps = 1000;
			//parameters.numOfThreads = 12;
			//parameters.dt = 1.0e-3;
			//parameters.outputFrequency = 10;
			//parameters.simulation_Mode = "ABD"; // Normal, ABD, Coupling
			//parameters.enableGround = false;
			//parameters.searchResidual = 6.0;
			//parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
			//parameters.rigidMode = true;
			//parameters.IPC_dis = 0.01;
			//parameters.IPC_eta = 0.05;
			//parameters.IPC_kStiffness = 1.0e12;
			//parameters.IPC_hashSize = tetSimMesh.calLargestEdgeLength() * 1.1;
			//parameters.IPC_B3Stiffness = 500;
			//parameters.ABD_Coeff = 1.0e9;




			//implicitFEM_ABD(tetSimMesh, parameters);


		}
		else if (caseNum == 7)
		{
			//Material mat1;
			//mat1.density = 8000;
			//mat1.E = 7.26e12;
			//mat1.updateDenpendecies();

			//Material mat2;
			//mat2.density = 800;
			//mat2.E = 7.26e12;
			//mat2.updateDenpendecies();



			//std::vector<meshConfiguration> config;
			//meshConfiguration m1;
			//m1.filePath = "../input/cube.msh";
			//m1.mesh_material = mat1;


			//m1.note = "cube_0";
			//Eigen::Vector3d trans = { -7.5, 2.0, 2.0 };
			//m1.translation = trans;
			//config.push_back(m1);




			//int count = 1;
			//for (int z = 0; z < 3; z++)
			//{
			//	for (int x = 0; x < 3; x++)
			//	{
			//		for (int y = 0; y < 3; y++)
			//		{
			//			count += 1;
			//			m1.mesh_material = mat2;
			//			m1.note = "cube_" + std::to_string(count);
			//			Eigen::Vector3d trans = { (double)x * 1.1, (double)y * 1.1, (double)z * 1.1 + 0.8 };
			//			m1.translation = trans;
			//			config.push_back(m1);
			//		}
			//	}
			//}



			//Mesh_ABD tetSimMesh;
			//for (int i = 0; i < config.size(); i++)
			//{
			//	tetMesh msh_tmp;
			//	msh_tmp.readMesh(config[i]);
			//	msh_tmp.initializeTetMesh();
			//	tetSimMesh.objectsTetMesh[msh_tmp.tetMeshNote] = msh_tmp;
			//}
			//tetSimMesh.createGlobalSimulationMesh_ABD();
			//for (int num = 1; num < tetSimMesh.translation_vel_ABD.size(); num++)
			//{
			//	tetSimMesh.translation_vel_ABD[num] = { 0,0,-3 };
			//}
			//tetSimMesh.translation_vel_ABD[0] = { 3,0,0 };



			//std::cout << "tetSimMesh.pos_node.size() = " << tetSimMesh.pos_node.size() << std::endl;
			//std::cout << "tetSimMesh.tetrahedrals.size() = " << tetSimMesh.tetrahedrals.size() << std::endl;



			//FEMParamters parameters;
			//parameters.gravity = { 0, 0, -9.8};
			//parameters.num_timesteps = 10000;
			//parameters.numOfThreads = 12;
			//parameters.dt = 1.0e-2;
			//parameters.outputFrequency = 20;
			//parameters.simulation_Mode = "ABD"; // Normal, ABD, Coupling
			//parameters.enableGround = true;
			//parameters.searchResidual = 0.5;
			//parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
			//parameters.rigidMode = true;
			//parameters.IPC_dis = 0.01;
			//parameters.IPC_eta = 0.05;
			//parameters.IPC_kStiffness = 1.0e9;
			//parameters.IPC_hashSize = tetSimMesh.calLargestEdgeLength() * 1.1;
			//parameters.IPC_B3Stiffness = 500;
			//parameters.ABD_Coeff = 1.0e10;



			//implicitFEM_ABD(tetSimMesh, parameters);

		}
		else if (caseNum == 8)
		{

			Material mat1;
			mat1.density = 8000;
			mat1.E = 7.26e12;
			mat1.updateDenpendecies();

			Material mat2;
			mat2.density = 800;
			mat2.E = 7.26e12;
			mat2.updateDenpendecies();



			std::vector<meshConfiguration> config;
			meshConfiguration m1;
			m1.filePath = "../input/cube_eq.obj";
			m1.mesh_material = mat1;
			m1.note = "cube_0";
			Eigen::Vector3d trans = { -7.5, 2.0, 2.0 };
			m1.translation = trans;
			m1.per_point_volume = 0.01;
			config.push_back(m1);




			int count = 1;
			for (int z = 0; z < 3; z++)
			{
				for (int x = 0; x < 3; x++)
				{
					for (int y = 0; y < 3; y++)
					{
						count += 1;
						m1.mesh_material = mat2;
						m1.note = "cube_" + std::to_string(count);
						Eigen::Vector3d trans = { (double)x * 1.1, (double)y * 1.1, (double)z * 1.1 + 0.8 };
						m1.translation = trans;
						config.push_back(m1);
					}
				}
			}



			triMesh triSimMesh;
			triSimMesh.createGlobalSimulationTriMesh_ABD(config);
			for (int num = 1; num < triSimMesh.translation_vel_ABD.size(); num++)
			{
				triSimMesh.translation_vel_ABD[num] = { 0,0,-3 };
			}
			triSimMesh.translation_vel_ABD[0] = { 3,0,0 };


			std::cout << "tetSimMesh.pos_node_surface.size() = " << triSimMesh.pos_node_surface.size() << std::endl;


			FEMParamters parameters;
			parameters.gravity = { 0, 0, -9.8};
			parameters.num_timesteps = 5000;
			parameters.numOfThreads = 12;
			parameters.dt = 1.0e-2;
			parameters.outputFrequency = 20;
			parameters.simulation_Mode = "ABD"; // Normal, ABD, Coupling
			parameters.enableGround = true;
			parameters.searchResidual = 0.5;
			parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
			parameters.rigidMode = true;
			parameters.IPC_dis = 0.01;
			parameters.IPC_eta = 0.05;
			parameters.IPC_kStiffness = 1.0e9;
			parameters.IPC_hashSize = triSimMesh.calLargestEdgeLength() * 1.1;
			parameters.IPC_B3Stiffness = 500;
			parameters.ABD_Coeff = 1.0e10;



			implicitFEM_ABD_triMesh(triSimMesh, parameters);

		}
		else if (caseNum == 9)
		{

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
			m1.note = "bunny";
			m1.breakable = true;
			m1.per_point_volume = 0.01;
			config.push_back(m1);

			meshConfiguration m2;
			m2.filePath = "../input/cube_eq.obj";
			m2.mesh_material = mat1;
			m2.note = "cube";
			Eigen::Vector3d trans = { 0, 1.5, 0 };
			m2.translation = trans;
			m2.velocity = {0,-20,0};
			m2.breakable = false;
			m2.per_point_volume = 0.01;
			config.push_back(m2);






			triMesh triSimMesh;
			triSimMesh.createGlobalSimulationTriMesh_ABD(config);


			std::cout << "tetSimMesh.pos_node_surface.size() = " << triSimMesh.pos_node_surface.size() << std::endl;


			FEMParamters parameters;
			parameters.gravity = { 0, 0, 0 };
			parameters.num_timesteps = 10000;
			parameters.numOfThreads = 12;
			parameters.dt = 1.0e-2;
			parameters.outputFrequency = 1;
			parameters.simulation_Mode = "ABD"; // Normal, ABD, Coupling
			parameters.enableGround = false;
			parameters.searchResidual = 0.5;
			parameters.model = "neoHookean"; // neoHookean ARAP ARAP_linear ACAP
			parameters.rigidMode = true;
			parameters.IPC_dis = 0.05;
			parameters.IPC_eta = 0.05;
			parameters.IPC_kStiffness = 1.0e12;
			parameters.IPC_hashSize = triSimMesh.calLargestEdgeLength() * 1.1;
			parameters.IPC_B3Stiffness = 500;
			parameters.ABD_Coeff = 1.0e10;



			implicitFEM_ABD_triMesh(triSimMesh, parameters);

		}

	}


	return 0;

}
