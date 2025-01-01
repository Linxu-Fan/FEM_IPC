#include "simulator.h" 
#include "tools.h"



// TODO
// 1. Accelerate and parallel openMP
// 2. Decimate the generated mesh for performance
// 3. triSimMesh.boundaryCondition_node_surface (& boundaryCondition_node_interior) is not properly handled




class ThreadPool {
public:
	ThreadPool(size_t numThreads) : stop(false) {
		for (size_t i = 0; i < numThreads; ++i) {
			workers.emplace_back([this]() {
				while (true) {
					std::function<void()> task;
					{
						std::unique_lock<std::mutex> lock(queueMutex);
						condition.wait(lock, [this]() { return stop || !tasks.empty(); });
						if (stop && tasks.empty()) return;
						task = std::move(tasks.front());
						tasks.pop();
					}
					task();
				}
				});
		}
	}

	template <typename F>
	void enqueue(F&& task) {
		{
			std::unique_lock<std::mutex> lock(queueMutex);
			if (stop) throw std::runtime_error("ThreadPool is stopped.");
			tasks.emplace(std::forward<F>(task));
		}
		condition.notify_one();
	}

	~ThreadPool() {
		{
			std::unique_lock<std::mutex> lock(queueMutex);
			stop = true;
		}
		condition.notify_all();
		for (std::thread& worker : workers) {
			if (worker.joinable()) worker.join();
		}
	}

private:
	std::vector<std::thread> workers;
	std::queue<std::function<void()>> tasks;
	std::mutex queueMutex;
	std::condition_variable condition;
	std::atomic<bool> stop;
};





struct AABB 
{
	Eigen::Vector3d min = Eigen::Vector3d::Ones() * 1.0e9;
	Eigen::Vector3d max = -Eigen::Vector3d::Ones() * 1.0e9;
	int face_index = -99;
	Eigen::Vector3d v0 = Eigen::Vector3d::Zero();
	Eigen::Vector3d v1 = Eigen::Vector3d::Zero();
	Eigen::Vector3d v2 = Eigen::Vector3d::Zero();


	void init(const std::vector<Eigen::Vector3d>& pos_node_surface, Eigen::Vector3i& face_vertices_, int index_, double dx)
	{
		face_index = index_;
		v0 = pos_node_surface[face_vertices_[0]];
		v1 = pos_node_surface[face_vertices_[0]];
		v2 = pos_node_surface[face_vertices_[0]];
		compute_min_max(pos_node_surface, dx);
	}

	void compute_min_max(const std::vector<Eigen::Vector3d>& pos_node_surface, double dx)
	{
		min = v0.cwiseMin(v1).cwiseMin(v2);
		max = v0.cwiseMax(v1).cwiseMax(v2);
		dilate(dx);
	}


	bool intersects(const AABB& other)
	{
		return (min.x() <= other.max.x() && max.x() >= other.min.x() &&
			min.y() <= other.max.y() && max.y() >= other.min.y() &&
			min.z() <= other.max.z() && max.z() >= other.min.z());
	}

	void dilate(double dx) {
		min -= Eigen::Vector3d(dx, dx, dx);
		max += Eigen::Vector3d(dx, dx, dx);
	}
};


struct BVHNode {
	AABB box;                
	BVHNode* left = nullptr; 
	BVHNode* right = nullptr; 
	size_t depth = 0;

	
	bool isLeaf() {
		return left == nullptr && right == nullptr;
	}
};



BVHNode* buildBVH(const std::vector<AABB>& triangles, int start, int end, size_t depth) {
	BVHNode* node = new BVHNode();
	node->depth = depth;


	if (end - start == 1) {
		node->box = triangles[start];
		return node;
	}

	//std::cout << "111.0" << std::endl;

	//std::cout << "triangles.size() = " << triangles.size() << "; start = " << start << "; end = " << end << std::endl;


	AABB box;
	for (int i = start; i < end; ++i) {
		//std::cout << "    i = " << i << std::endl;
		box.min = box.min.cwiseMin(triangles[i].min);
		box.max = box.max.cwiseMax(triangles[i].max);
	}
	node->box = box;

	//std::cout << "111.1" << std::endl;


	int axis = 0; 
	Eigen::Vector3d size = box.max - box.min;
	if (size.y() > size.x() && size.y() > size.z()) axis = 1; // Y ��
	if (size.z() > size.x() && size.z() > size.y()) axis = 2; // Z ��

	//std::cout << "axis = " << axis << std::endl;

	std::vector<AABB> sortedTriangles(triangles.begin() + start, triangles.begin() + end);
	std::sort(sortedTriangles.begin(), sortedTriangles.end(), [axis](const AABB& a, const AABB& b) {
		return a.min[axis] < b.min[axis];
		});





	int mid = sortedTriangles.size() / 2;

	/*std::cout << "111.2" << std::endl;
	std::cout << "sortedTriangles.size() = " << sortedTriangles.size() <<"; mid = "<<mid<<"; end = "<< sortedTriangles.size() <<std::endl << std::endl;*/
	node->left = buildBVH(sortedTriangles, 0, mid, depth + 1);
	node->right = buildBVH(sortedTriangles, mid, sortedTriangles.size(), depth + 1);




	return node;
}


void updateBVHLeafNodes(BVHNode* node, const std::vector<Eigen::Vector3d>& pos_node_surface, double dilation) {
	if (!node) return;


	if (node->isLeaf()) 
	{
		node->box.compute_min_max(pos_node_surface, dilation);
		return;
	}


	updateBVHLeafNodes(node->left, pos_node_surface, dilation);
	updateBVHLeafNodes(node->right, pos_node_surface, dilation);


	node->box.min = node->left->box.min.cwiseMin(node->right->box.min);
	node->box.max = node->left->box.max.cwiseMax(node->right->box.max);
}


void parallelUpdateBVH(BVHNode* root, const std::vector<Eigen::Vector3d>& pos_node_surface, double dilation) {
	if (!root) return;

	std::vector<std::thread> threads;
	std::mutex mutex; 


	auto task = [&](BVHNode* node) {
		updateBVHLeafNodes(node, pos_node_surface, dilation);
		};


	threads.emplace_back(task, root->left);
	threads.emplace_back(task, root->right);


	for (auto& thread : threads) {
		thread.join();
	}
}


void queryBVH(BVHNode* nodeA, BVHNode* nodeB, std::vector<std::pair<int, int>>& results) {
	if (!nodeA || !nodeB) return;


	if (!nodeA->box.intersects(nodeB->box))
	{
		//std::cout << "Level = " << nodeA->depth << "; " << nodeB->depth << std::endl;
		return;
	}
	else
	{
		if (nodeA->isLeaf())
		{
			if (nodeB->isLeaf())
			{
				//std::cout << "intersect = " << "(" << nodeA->triangleIndex << "," << nodeB->triangleIndex << ")";
				//std::cout << ";add = " << "(" << nodeA << "," << nodeB << ")" << std::endl;
				results.emplace_back(nodeA->box.face_index, nodeB->box.face_index);
				return;
			}
			else
			{
				queryBVH(nodeA, nodeB->right, results);
				queryBVH(nodeA, nodeB->left, results);
			}
		}
		else
		{
			if (nodeB->isLeaf())
			{
				queryBVH(nodeA->right, nodeB, results);
				queryBVH(nodeA->left, nodeB, results);
			}
			else
			{
				queryBVH(nodeA->left, nodeB->left, results);
				queryBVH(nodeA->left, nodeB->right, results);
				queryBVH(nodeA->right, nodeB->left, results);
				queryBVH(nodeA->right, nodeB->right, results);
			}
		}

	}

}

void count_leaf_node(BVHNode* nodeA, int& num_leaves)
{
	if (nodeA->isLeaf())
	{
		num_leaves += 1;

		if (nodeA->box.face_index == 0)
		{
			std::cout << "min = " << nodeA->box.min.transpose() << "; max = " << nodeA->box.max.transpose() << std::endl;
		}
	}
	else
	{
		count_leaf_node(nodeA->left, num_leaves);
		count_leaf_node(nodeA->right, num_leaves);
	}
}







void parallelQueryBVH(BVHNode* nodeA, BVHNode* nodeB, ThreadPool& threadPool,
	std::vector<std::vector<std::pair<int, int>>>& localResults, int threadIndex, int maxDepth) {
	if (!nodeA || !nodeB) return;

	if (!nodeA->box.intersects(nodeB->box)) return;

	if (nodeA->isLeaf() && nodeB->isLeaf()) {
		localResults[threadIndex].emplace_back(nodeA->box.face_index, nodeB->box.face_index);
		return;
	}

	if (maxDepth <= 0) {

		if (!nodeA->isLeaf()) {
			parallelQueryBVH(nodeA->left, nodeB, threadPool, localResults, threadIndex, maxDepth - 1);
			parallelQueryBVH(nodeA->right, nodeB, threadPool, localResults, threadIndex, maxDepth - 1);
		}
		if (!nodeB->isLeaf()) {
			parallelQueryBVH(nodeA, nodeB->left, threadPool, localResults, threadIndex, maxDepth - 1);
			parallelQueryBVH(nodeA, nodeB->right, threadPool, localResults, threadIndex, maxDepth - 1);
		}
		return;
	}


	if (!nodeA->isLeaf()) {
		threadPool.enqueue([=, &localResults, &threadPool]() {
			parallelQueryBVH(nodeA->left, nodeB, threadPool, localResults, threadIndex, maxDepth - 1);
			});
		threadPool.enqueue([=, &localResults, &threadPool]() {
			parallelQueryBVH(nodeA->right, nodeB, threadPool, localResults, threadIndex, maxDepth - 1);
			});
	}

	if (!nodeB->isLeaf()) {
		threadPool.enqueue([=, &localResults, &threadPool]() {
			parallelQueryBVH(nodeA, nodeB->left, threadPool, localResults, threadIndex, maxDepth - 1);
			});
		threadPool.enqueue([=, &localResults, &threadPool]() {
			parallelQueryBVH(nodeA, nodeB->right, threadPool, localResults, threadIndex, maxDepth - 1);
			});
	}
}














int main()
{

	if (1)
	{
		
		double dilation = 0;

		std::string bunny = "D:/Research/Hydrostatic_object/code/FEM_IPC/input/bunny_highRes.obj";
		objMeshFormat bunny1;
		bunny1.readObjFile(bunny, false);
		std::vector<AABB> bunny1_AABBs;
		for (size_t i = 0; i < bunny1.faces.size(); ++i) 
		{
			AABB aabb_;
			aabb_.init(bunny1.vertices, bunny1.faces[i],i, dilation);
			bunny1_AABBs.push_back(aabb_);
		}


		std::cout << "111" << std::endl;
		BVHNode* bunny1_root = buildBVH(bunny1_AABBs, 0, bunny1_AABBs.size(), 0);

		int num_leaves = 0;
		count_leaf_node(bunny1_root, num_leaves);
		std::cout << "num_leaves = " << num_leaves << std::endl;

		std::cout << "222" << std::endl;



		objMeshFormat bunny2 = bunny1;
		Eigen::Vector3d translation = { 3.0066,0,0 };
		for (int i = 0; i < bunny1.vertices.size(); i++)
		{
			bunny2.vertices[i] = bunny1.vertices[i] + translation;
		}
		std::vector<AABB> bunny2_AABBs;
		for (size_t i = 0; i < bunny2.faces.size(); ++i)
		{
			AABB aabb_;
			aabb_.init(bunny2.vertices, bunny2.faces[i], i, dilation);
			bunny2_AABBs.push_back(aabb_);
		}
		BVHNode* bunny2_root = buildBVH(bunny2_AABBs, 0, bunny2_AABBs.size(), 0);

		std::cout << "333" << std::endl;


		bunny1.outputFile("bunny1", -99);
		bunny2.outputFile("bunny2", -99);

		std::vector<std::pair<int, int>> results;


		double startTime1, endTime1;
		startTime1 = omp_get_wtime();

		queryBVH(bunny1_root, bunny2_root, results);

		endTime1 = omp_get_wtime();
		std::cout << "Query Time is : " << endTime1 - startTime1 << "s" << std::endl;


		{
			std::ofstream outfile9("./output/collision.obj", std::ios::trunc);
			for (int k = 0; k < bunny1.vertices.size(); k++)
			{
				Eigen::Vector3d scale = bunny1.vertices[k];
				outfile9 << std::scientific << std::setprecision(8) << "v " << scale[0] << " " << scale[1] << " " << scale[2] << std::endl;
			}

			for (int km = 0; km < results.size(); km++)
			{
				int k = results[km].first;
				outfile9 << "f ";
				for (int m = 0; m < bunny1.faces[k].size(); m++)
				{
					outfile9 << bunny1.faces[k][m] + 1 << " ";
				}
				outfile9 << std::endl;
			}


			for (int k = 0; k < bunny2.vertices.size(); k++)
			{
				Eigen::Vector3d scale = bunny2.vertices[k];
				outfile9 << std::scientific << std::setprecision(8) << "v " << scale[0] << " " << scale[1] << " " << scale[2] << std::endl;
			}

			for (int km = 0; km < results.size(); km++)
			{
				int k = results[km].second;
				outfile9 << "f ";
				for (int m = 0; m < bunny2.faces[k].size(); m++)
				{
					outfile9 << bunny2.faces[k][m] + 1 + bunny1.vertices.size() << " ";
				}
				outfile9 << std::endl;
			}

		
			outfile9.close();
		}


		std::cout << "results = " << results.size() << std::endl;
		{
			std::ofstream outfile9("./output/collision.txt", std::ios::trunc);
			for (int km = 0; km < results.size(); km++)
			{
				outfile9 << results[km].first << " - " << results[km].second << std::endl;
			}
			outfile9.close();
		}





		Eigen::Vector3d translation_2 = { -1.0,0,0 };
		for (int i = 0; i < bunny2.vertices.size(); i++)
		{
			bunny2.vertices[i] = bunny2.vertices[i] + translation_2;
		}



		parallelUpdateBVH(bunny2_root, bunny2.vertices, dilation);

		std::vector<std::pair<int, int>> results2;
		queryBVH(bunny1_root, bunny2_root, results2);
		std::cout << "results2 = " << results2.size() << std::endl;


		{
			std::ofstream outfile9("./output/collision2.obj", std::ios::trunc);
			for (int k = 0; k < bunny1.vertices.size(); k++)
			{
				Eigen::Vector3d scale = bunny1.vertices[k];
				outfile9 << std::scientific << std::setprecision(8) << "v " << scale[0] << " " << scale[1] << " " << scale[2] << std::endl;
			}

			for (int km = 0; km < results2.size(); km++)
			{
				int k = results2[km].first;
				outfile9 << "f ";
				for (int m = 0; m < bunny1.faces[k].size(); m++)
				{
					outfile9 << bunny1.faces[k][m] + 1 << " ";
				}
				outfile9 << std::endl;
			}


			for (int k = 0; k < bunny2.vertices.size(); k++)
			{
				Eigen::Vector3d scale = bunny2.vertices[k];
				outfile9 << std::scientific << std::setprecision(8) << "v " << scale[0] << " " << scale[1] << " " << scale[2] << std::endl;
			}

			for (int km = 0; km < results2.size(); km++)
			{
				int k = results2[km].second;
				outfile9 << "f ";
				for (int m = 0; m < bunny2.faces[k].size(); m++)
				{
					outfile9 << bunny2.faces[k][m] + 1 + bunny1.vertices.size() << " ";
				}
				outfile9 << std::endl;
			}


			outfile9.close();
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

			Material mat1;
			mat1.density = 800;
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
			Eigen::Vector3d trans = { 1.5, 1.5, 4.8 };
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
						Eigen::Vector3d trans = { (double)x * 1.03, (double)y * 1.03, (double)z * 1.03 + 0.5 };
						m1.translation = trans;
						config.push_back(m1);
					}
				}
			}



			triMesh triSimMesh;
			triSimMesh.createGlobalSimulationTriMesh_ABD(config);
			for (int num = 1; num < triSimMesh.translation_vel_ABD.size(); num++)
			{
				triSimMesh.translation_vel_ABD[num] = { 0,0,0 };
			}
			triSimMesh.translation_vel_ABD[0] = { 0,0,-1 };


			std::cout << "tetSimMesh.pos_node_surface.size() = " << triSimMesh.pos_node_surface.size() << std::endl;


			FEMParamters parameters;
			parameters.gravity = { 0, 0, -9.8};
			parameters.num_timesteps = 5000;
			parameters.numOfThreads = 12;
			parameters.dt = 1.0e-2;
			parameters.outputFrequency = 5;
			parameters.simulation_Mode = "ABD"; // Normal, ABD, Coupling
			parameters.enableGround = true;
			parameters.searchResidual = 0.05;
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
		else if (caseNum == 1)
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
			m1.velocity = { 0,1.0,0 };
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
