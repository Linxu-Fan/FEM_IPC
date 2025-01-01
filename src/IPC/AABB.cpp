//#include "AABB.h"
//
//void AABB::init(const std::vector<Eigen::Vector3d>& pos_node_surface, Eigen::Vector3i& face_vertices_, int index_, double dx)
//{
//	face_index = index_;
//	face_vertices = face_vertices_;
//	compute_min_max(pos_node_surface, dx);
//}
//
//void AABB::compute_min_max(const std::vector<Eigen::Vector3d>& pos_node_surface, double dx)
//{
//	min = pos_node_surface[face_vertices[0]].cwiseMin(pos_node_surface[face_vertices[1]]).cwiseMin(pos_node_surface[face_vertices[2]]);
//	max = pos_node_surface[face_vertices[0]].cwiseMax(pos_node_surface[face_vertices[1]]).cwiseMax(pos_node_surface[face_vertices[2]]);
//	dilate(dx);
//}
//
//
//bool AABB::intersects(const AABB& other)
//{
//	return (min.x() <= other.max.x() && max.x() >= other.min.x() &&
//		min.y() <= other.max.y() && max.y() >= other.min.y() &&
//		min.z() <= other.max.z() && max.z() >= other.min.z());
//}
//
//void AABB::dilate(double dx) {
//	min -= Eigen::Vector3d(dx, dx, dx);
//	max += Eigen::Vector3d(dx, dx, dx);
//}