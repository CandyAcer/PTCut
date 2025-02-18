// $Intro: 本文件是切割算法的第一部分, 负责计算交点, 构建细化和提取所需的信息.
// $License: GPL-3.0-or-later
// 
// 
// $Author: Liu Zhiqiang



#ifndef MCAL_ALGORITHM_INTERSECTION_IMPL_H
#define MCAL_ALGORITHM_INTERSECTION_IMPL_H


#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/iterator.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Kernel/global_functions_3.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/functional/hash.hpp>

#include <unordered_set>
#include <unordered_map>

#include "utils/sorted_pair.h"
#include "intersection_nodes.h"
#include "intersection_callback.h"
#include "intersection_type.h"

extern int iExactComputeCount;
extern int iIntervalComputeCount;

namespace MCAL    // Mesh Cut Algorithm Libaray, 存放多面体网格切割算法的所有代码
{

	// 阅读代码需要注意的一个地方:
	// filter_intersections(), compute_intersection_points()等函数都会有两个版本, 
	// 这是由于polyhedron grid和surface mesh的行为不同导致的, 为了区分, 我们在函数后面加了后缀
	// 
	// sepf: sm edge to pg face
	// pesf: pg edge to sm face
	//


	/*
	* IntersectionImpl是切割算法的第一部分, 大体流程如下:
	*   1) 检测表面网格的边与哪些多面体网格面相交, 对可能相交的多面体网格面三角化,
	*      方便几何计算; 检测多面体网格的边与哪些表面网格面可能相交.
	*   2) 对这些筛选出的信息, 计算边与三角形的交点, 同时记录相关信息, 主要有以下几类：
	*        2.1) 交点的坐标信息, 交由IntersectionNodes统一管理, 由NodeId进行标识.
	*        2.2) 点落在哪个面, 哪条边或哪个点上, 记录这类对应关系.
	*        2.3) 相邻信息, 两个三角形交线上的两个点相邻.
	*   3) 以cell为单位构建交线.
	*
	* 2记录的信息对细化有用, 故记录在细化类的data member中.
	* 3记录的信息则是用来建立最终结果, 故记录在OutputBuilder的data member中.
	*
	* @param IntersectionVisitor: 负责对网格进行细化的类, 作为本类的Visitor.
	*/
	template <typename SurfaceMesh,
		typename PolyhedronGrid,
		typename VPMSM,
		typename VPMPG,
		typename IntersectionVisitor>
	class IntersectionImpl
	{
		// typedefs
		typedef boost::graph_traits<SurfaceMesh>             GT;
		typedef typename GT::vertex_descriptor               SMVertex;
		typedef typename GT::edge_descriptor                 SMEdge;
		typedef typename GT::halfedge_descriptor             SMHalfedge;
		typedef typename GT::face_descriptor                 SMFace;

		typedef typename PolyhedronGrid::VertexIndex         PGVertex;
		typedef typename PolyhedronGrid::EdgeIndex           PGEdge;
		typedef typename PolyhedronGrid::FaceIndex           PGFace;
		typedef typename PolyhedronGrid::CellIndex           PGCell;
		typedef typename PolyhedronGrid::Halfedge            PGHalfedge;

		typedef typename boost::property_traits<VPMSM>::value_type  Point3;
		static_assert((std::is_same<Point3, typename PolyhedronGrid::Point>::value),
			"输入网格的点类型不一致");

		// Axis aligned bounding box (with info).
		typedef BoxInfo<SurfaceMesh, PolyhedronGrid>         Info;
		typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, Info>  Box;

		typedef std::unordered_set<PGFace> PGFaceSet;
		typedef std::unordered_map<SMEdge, PGFaceSet> SMEdgeToPGFaces;
		typedef std::unordered_set<SMFace> SMFaceSet;
		typedef std::unordered_map<PGEdge, SMFaceSet> PGEdgeToSMFaces;

		typedef IntersectionNodes<SurfaceMesh, PolyhedronGrid, VPMSM, VPMPG> NodeVector;

		typedef std::size_t NodeId;

		// adjacent information.
		typedef std::pair<SMFace, PGFace> FacePair;
		typedef std::unordered_map<FacePair, SortedPair<NodeId>, boost::hash<FacePair>> FacesToNodes;

		typedef std::set<FacePair> CoplanarFaceSet;

		// data members
		SMEdgeToPGFaces sm_edge_to_pg_faces;  // intersected surface mesh edge and polyhedron grid faces.
		PGEdgeToSMFaces pg_edge_to_sm_faces;  // intersected polyhedron grid edge and surface mesh faces.
		CoplanarFaceSet coplanar_faces;  // coplanar surface mesh face and polyhedron grid face. we don't handle their intersections.
		NodeVector nodes; // manager all intersection nodes, store exact points to offer exact predications.
		IntersectionVisitor visitor;   // intersection visitor to collect information used in triangulation.
		FacesToNodes f_to_node;   // Associate a pair of triangles to their intersection points.

		// member functions
	private:

		// 检测surface mesh的edge与哪些polyhedron grid的face相交(指包围盒相交, 非实际相交).
		// 过滤掉surface mesh中包围盒不干涉的edge.
		void filter_intersections_sepf(SurfaceMesh& sm,
			PolyhedronGrid& pg,
			VPMSM& vpm_sm,
			VPMPG& vpm_pg)
		{
			std::vector<Box> sm_edge_boxes;
			std::vector<Box> pg_face_boxes;

			// 为Surface mesh所有的edge套一个AABB包围盒.
			sm_edge_boxes.reserve(num_edges(sm));
			for (SMEdge e_sm : edges(sm))
			{
				SMHalfedge h_sm = halfedge(e_sm, sm);
				sm_edge_boxes.push_back(
					Box(get(vpm_sm, source(h_sm, sm)).bbox() +
						get(vpm_sm, target(h_sm, sm)).bbox(),
						Info(h_sm))
				);
			}

			// 为polyhedron grid的所有face套一个AABB包围盒.
			pg_face_boxes.reserve(pg.num_faces());
			for (PGFace face : pg.faces())
			{
				std::vector<PGVertex> vset;
				pg.incident_vertices(face, std::back_inserter(vset));
				CGAL::Bbox_3 face_bbox;
				for (PGVertex v_pg : vset)
					face_bbox += vpm_pg[v_pg].bbox();

				PGHalfedge h_pg = pg.halfedge(face);
				pg_face_boxes.push_back(Box(face_bbox, Info(h_pg)));
			}

			// callback在检测到包围盒相交时进行记录.
			typedef CollectSMEdgeToPGFaces<SurfaceMesh, PolyhedronGrid, SMEdgeToPGFaces>
				Callback;
			Callback callback(sm, pg, sm_edge_to_pg_faces);

			// 调用相交检测算法.
			CGAL::box_intersection_d(sm_edge_boxes.begin(), sm_edge_boxes.end(),
				pg_face_boxes.begin(), pg_face_boxes.end(),
				callback);

			// 由于多面体网格的面是空间多边形, 相较于三角形更为复杂, 需要在此将检测到的面进行三角化.
			std::unordered_map<PGFace, std::vector<PGFace>> polygon_to_triangles;
			for (typename SMEdgeToPGFaces::iterator it = sm_edge_to_pg_faces.begin();
				it != sm_edge_to_pg_faces.end();
				++it)
			{
				SMEdge e_sm = it->first;
				// 此处不能是引用, 由于后续会往set中添加元素, 引用会导致set中动态增加, 
				// 我们要遍历的就是原来的那些面.
				PGFaceSet polygon_set = it->second;

				for (PGFace polygon : polygon_set)
				{
					// 对干涉的多面体网格面进行三角化, 
					// 三角化后的面记录在传入的容器中, 其首元素一定是polygon.
					std::vector<PGFace> refined_triangles;
					if (!polygon_to_triangles.count(polygon))
					{
						pg.triangulate_a_polygon(polygon, std::back_inserter(refined_triangles));
						assert(polygon == refined_triangles.front());
						polygon_to_triangles.insert(std::make_pair(polygon, refined_triangles));
					}
					else
						refined_triangles = polygon_to_triangles[polygon];

					// 记录好三角化后的所有面.
					sm_edge_to_pg_faces[e_sm].insert(refined_triangles.begin(), refined_triangles.end());
				}
			}
		}

		// 检测polyhedron grid的edge与哪些Surface mesh的face相交(指包围盒相交, 非实际相交).
		// 过滤掉polyhedron grid中包围盒不干涉的edge.
		void filter_intersections_pesf(PolyhedronGrid& pg,
			SurfaceMesh& sm,
			VPMPG& vpm_pg,
			VPMSM& vpm_sm)
		{
			std::vector<Box> pg_cell_boxes;
			std::vector<Box> sm_face_boxes;

			// 不能为edge套包围盒, 因为它关联的面可能会存在非三角形面, 这样判断共面就会很麻烦.
			// 既然如此, 是否可以调用box_intersection_d筛出可能相交的edge, 然后三角化其关联的面?
			// 可以. 但这会导致多一次filter_pg_intersections(), 效率很低.
			// 
			// Note: 
			// 上面的思路三角化后新产生的edge包围盒可能更大(如果是斜边的话), 但不会有错误, 
			// 因为新产生的edge一定在多边形的内部, 它关联的面都是三角形面, 不会陷入死循环.

			// 为polyhedron grid的所有cell套上AABB包围盒.
			pg_cell_boxes.reserve(pg.num_cells());
			for (PGCell c : pg.cells())
			{
				std::vector<PGVertex> vset;
				pg.incident_vertices(c, std::back_inserter(vset));
				CGAL::Bbox_3 cell_bbox;
				for (PGVertex v_pg : vset)
					cell_bbox += vpm_pg[v_pg].bbox();

				PGHalfedge h_pg = pg.halfedge(c);
				pg_cell_boxes.push_back(Box(cell_bbox, Info(h_pg)));
			}

			// 为surface mesh所有的face套上AABB包围盒.
			sm_face_boxes.reserve(num_faces(sm));
			for (SMFace f_sm : faces(sm))
			{
				SMHalfedge h_sm = halfedge(f_sm, sm);
				sm_face_boxes.push_back(
					Box(get(vpm_sm, source(h_sm, sm)).bbox() +
						get(vpm_sm, target(h_sm, sm)).bbox() +
						get(vpm_sm, target(next(h_sm, sm), sm)).bbox(),
						Info(h_sm))
				);
			}

			// callback负责检测到包围盒相交时进行记录.
			typedef std::unordered_map<PGCell, std::unordered_set<SMFace>> PGCellToSMFaces;
			PGCellToSMFaces pg_cell_to_sm_faces;

			typedef CollectPGCellToSMFaces<PolyhedronGrid, SurfaceMesh, PGCellToSMFaces>
				Callback;
			Callback callback(pg, sm, pg_cell_to_sm_faces);

			// 调用相交检测算法.
			CGAL::box_intersection_d(pg_cell_boxes.begin(), pg_cell_boxes.end(),
				sm_face_boxes.begin(), sm_face_boxes.end(),
				callback);

			// 之所以使用多面体网格的单元而不是边进行相交检测, 是因为它的包围盒比包含的所有边都大,
			for (typename PGCellToSMFaces::iterator it = pg_cell_to_sm_faces.begin();
				it != pg_cell_to_sm_faces.end();
				++it)
			{
				PGCell c = it->first;
				SMFaceSet& fset_sm = it->second;  // 此处可以是引用, 因为不会insert, size固定不变.

				std::vector<PGFace> cell_faces;
				pg.incident_faces(c, std::back_inserter(cell_faces));

				for (PGFace polygon : cell_faces)
					pg.triangulate_a_polygon(polygon, CGAL::Emptyset_iterator());

				// 收集多边形polygon三角化后的所有edge.
				std::vector<PGEdge> all_edges;
				pg.incident_edges(c, std::back_inserter(all_edges));
				// 对于三角化后的所有边, 由于边的包围盒一般比cell的小,
				// 所以不一定与fset_sm中的所有表面网格面干涉, 需要剔除.
				for (PGEdge e_pg : all_edges)
				{
					PGHalfedge h_pg = pg.halfedge(e_pg);
					CGAL::Bbox_3 edge_bbox(
						vpm_pg[pg.source(h_pg)].bbox() +
						vpm_pg[pg.target(h_pg)].bbox());

					for (SMFace f_sm : fset_sm)
					{
						SMHalfedge h_sm = halfedge(f_sm, sm);
						CGAL::Bbox_3 face_bbox(
							get(vpm_sm, source(h_sm, sm)).bbox() +
							get(vpm_sm, target(h_sm, sm)).bbox() +
							get(vpm_sm, target(next(h_sm, sm), sm)).bbox());

						// 若pg edge的包围盒和sm face的包围盒确实干涉, 记录.
						if (CGAL::do_overlap(edge_bbox, face_bbox))
							pg_edge_to_sm_faces[e_pg].insert(f_sm);
					}
				}
			}
		}

		// 判断sm face和pg face是否共面.
		bool are_two_face_coplanar(SMFace& f_sm,
			PGFace& f_pg,
			SurfaceMesh& sm,
			PolyhedronGrid& pg,
			VPMSM& vpm_sm,
			VPMPG& vpm_pg)
		{
			assert(pg.is_triangle(f_pg));

			SMHalfedge h_sm = halfedge(f_sm, sm);
			Point3 a = get(vpm_sm, source(h_sm, sm));
			Point3 b = get(vpm_sm, target(h_sm, sm));
			Point3 c = get(vpm_sm, target(next(h_sm, sm), sm));

			std::vector<PGVertex> vset;
			pg.incident_vertices(f_pg, std::back_inserter(vset));
			Point3 p = vpm_pg[vset[0]];
			Point3 q = vpm_pg[vset[1]];
			Point3 r = vpm_pg[vset[2]];

			CGAL::Orientation abcp = CGAL::orientation(a, b, c, p);
			CGAL::Orientation abcq = CGAL::orientation(a, b, c, q);
			CGAL::Orientation abcr = CGAL::orientation(a, b, c, r);

			// 此处应该做剔除, 暂不处理.
			if (abcp != CGAL::COPLANAR || abcq != CGAL::COPLANAR || abcr != CGAL::COPLANAR)
				return false;
			else
				return true;
		}

		// filter coplanar sm face and pg face, we don't handle their intersections.
		void filter_coplanar_case(SurfaceMesh& sm,
			PolyhedronGrid& pg,
			VPMSM& vpm_sm,
			VPMPG& vpm_pg)
		{
			// 遍历sm_edge_to_pg_faces, 检测共面的情况.
			for (typename SMEdgeToPGFaces::iterator it = sm_edge_to_pg_faces.begin();
				it != sm_edge_to_pg_faces.end();
				++it)
			{
				SMEdge e_sm = it->first;
				PGFaceSet& fset_pg = it->second;

				for (PGFace f_pg : fset_pg)
				{
					SMHalfedge h_sm = halfedge(e_sm, sm);
					if (is_border(h_sm, sm))
						h_sm = opposite(h_sm, sm);

					SMFace f_sm = face(h_sm, sm);
					if (are_two_face_coplanar(f_sm, f_pg, sm, pg, vpm_sm, vpm_pg))
						coplanar_faces.insert(std::make_pair(f_sm, f_pg));

					SMHalfedge opp_h_sm = opposite(h_sm, sm);
					if (is_border(opp_h_sm, sm))
						continue;
					SMFace opp_f_sm = face(opp_h_sm, sm);
					if (are_two_face_coplanar(opp_f_sm, f_pg, sm, pg, vpm_sm, vpm_pg))
						coplanar_faces.insert(std::make_pair(opp_f_sm, f_pg));
				}
			}

			// 遍历pg_edge_to_sm_faces, 检测共面的情况.
			for (typename PGEdgeToSMFaces::iterator it = pg_edge_to_sm_faces.begin();
				it != pg_edge_to_sm_faces.end();
				++it)
			{
				PGEdge e_pg = it->first;
				SMFaceSet& fset_sm = it->second;

				std::vector<PGFace> fset_pg;
				pg.incident_faces(e_pg, std::back_inserter(fset_pg));

				for (SMFace f_sm : fset_sm)
				{
					for (PGFace f_pg : fset_pg)
					{
						if (are_two_face_coplanar(f_sm, f_pg, sm, pg, vpm_sm, vpm_pg))
							coplanar_faces.insert(std::make_pair(f_sm, f_pg));
					}
				}
			}
		}

		std::tuple<IntersectionType, PGHalfedge, bool, bool>
			find_intersection_sepf(Point3& p, Point3& q,  // segment
				Point3& a, Point3& b, Point3& c,   // triangle
				PGHalfedge h_pg, // halfedge of the triangle face, its target is a
				PolyhedronGrid& pg,
				bool is_src_coplanar = false,
				bool is_tgt_coplanar = false)
		{
			typedef std::tuple<IntersectionType, PGHalfedge, bool, bool>     result_type;

			CGAL::Orientation ab = CGAL::orientation(p, q, a, b);
			CGAL::Orientation bc = CGAL::orientation(p, q, b, c);
			CGAL::Orientation ca = CGAL::orientation(p, q, c, a);

			using CGAL::POSITIVE;
			using CGAL::NEGATIVE;
			using CGAL::COPLANAR;

			if (ab == POSITIVE || bc == POSITIVE || ca == POSITIVE)
				return result_type(EMPTY, pg.null_halfedge(), false, false);

			int nb_coplanar = (ab == COPLANAR ? 1 : 0) + (bc == COPLANAR ? 1 : 0) + (ca == COPLANAR ? 1 : 0);

			if (nb_coplanar == 0)
				return result_type(ON_FACE, h_pg, is_src_coplanar, is_tgt_coplanar);

			if (nb_coplanar == 1)
			{
				if (ab == COPLANAR)      // intersection is ab
					return result_type(ON_EDGE, pg.next(h_pg), is_src_coplanar, is_tgt_coplanar);
				if (bc == COPLANAR)      // intersection is bc
					return result_type(ON_EDGE, pg.prev(h_pg), is_src_coplanar, is_tgt_coplanar);
				assert(ca == COPLANAR);
				// intersection is ca
				return result_type(ON_EDGE, h_pg, is_src_coplanar, is_tgt_coplanar);
			}

			assert(nb_coplanar == 2);

			if (ab != COPLANAR)		// intersection is c
				return result_type(ON_VERTEX, pg.prev(h_pg), is_src_coplanar, is_tgt_coplanar);
			if (bc != COPLANAR)		// intersection is a
				return result_type(ON_VERTEX, h_pg, is_src_coplanar, is_tgt_coplanar);
			assert(ca != COPLANAR);
			// intersection is b
			return result_type(ON_VERTEX, pg.next(h_pg), is_src_coplanar, is_tgt_coplanar);
		}

		std::tuple<IntersectionType, PGHalfedge, bool, bool>
			compute_intersection_type_sepf(SMHalfedge h_sm,
				PGFace f_pg,
				SurfaceMesh& sm,
				PolyhedronGrid& pg,
				VPMSM& vpm_sm,
				VPMPG& vpm_pg)
		{
			typedef std::tuple<IntersectionType, PGHalfedge, bool, bool>     result_type;
			typedef typename CGAL::Kernel_traits<Point3>::Kernel             Kernel;

			assert(pg.is_triangle(f_pg));
			PGHalfedge h_pg = pg.halfedge(f_pg);

			Point3 a = vpm_pg[pg.target(h_pg)];
			Point3 b = vpm_pg[pg.target(pg.next(h_pg))];
			Point3 c = vpm_pg[pg.source(h_pg)];
			Point3 p = get(vpm_sm, source(h_sm, sm));
			Point3 q = get(vpm_sm, target(h_sm, sm));

			const CGAL::Orientation abcp = CGAL::orientation(a, b, c, p);
			const CGAL::Orientation abcq = CGAL::orientation(a, b, c, q);

			switch (abcp)
			{
				// 按照编程规范来说, using语句是禁用的, 我们仅在switch域内使用, 保证影响范围尽可能小.
				// 目的是为了使得代码尽可能简洁.
				using CGAL::POSITIVE;
				using CGAL::NEGATIVE;
				using CGAL::COPLANAR;
			case POSITIVE:
				switch (abcq)
				{
				case POSITIVE:
					// the segment lies in the positive open halfspaces defined by the
					// triangle's supporting plane
					return result_type(EMPTY, pg.null_halfedge(), false, false);
				case NEGATIVE:
					// p sees the triangle in counterclockwise order
					return find_intersection_sepf(p, q, a, b, c, h_pg, pg);
				case COPLANAR:
					// q belongs to the triangle's supporting plane
					// p sees the triangle in counterclockwise order
					return find_intersection_sepf(p, q, a, b, c, h_pg, pg, false, true);
				default: // should not happen.
					assert(false);
					return result_type(EMPTY, pg.null_halfedge(), false, false);
				}
			case NEGATIVE:
				switch (abcq)
				{
				case POSITIVE:
					// q sees the triangle in counterclockwise order
					return find_intersection_sepf(q, p, a, b, c, h_pg, pg);
				case NEGATIVE:
					// the segment lies in the negative open halfspaces defined by the
					// triangle's supporting plane
					return result_type(EMPTY, pg.null_halfedge(), false, false);
				case COPLANAR:
					// q belongs to the triangle's supporting plane
					// p sees the triangle in clockwise order
					return find_intersection_sepf(q, p, a, b, c, h_pg, pg, false, true);
				default:
					assert(false);
					return result_type(EMPTY, pg.null_halfedge(), false, true);
				}
			case COPLANAR: // p belongs to the triangle's supporting plane
				switch (abcq)
				{
				case POSITIVE:
					// q sees the triangle in counterclockwise order
					return find_intersection_sepf(q, p, a, b, c, h_pg, pg, true);
				case NEGATIVE:
					// q sees the triangle in clockwise order
					return find_intersection_sepf(p, q, a, b, c, h_pg, pg, true);
				case COPLANAR:
					// the segment is coplanar with the triangle's supporting plane
					// we test whether the segment intersects the triangle in the common
					// supporting plane
					if (::CGAL::Intersections::internal::do_intersect_coplanar(a, b, c, p, q, Kernel()))
						return result_type(COPLANAR_TRIANGLES, pg.null_halfedge(), true, true);
					return result_type(EMPTY, pg.null_halfedge(), true, true);

				default: // should not happen.
					assert(false);
					return result_type(EMPTY, pg.null_halfedge(), false, false);
				}
			default: // should not happen.
				assert(false);
				return result_type(EMPTY, pg.null_halfedge(), false, false);
			}
		}

		// 计算sm edge和pg face的交点坐标, 以精确形式记录在nodes中, 由node id唯一标识.
		void add_new_node_sepf(SMHalfedge h_sm,
			PGFace f_pg,
			SurfaceMesh& sm,
			PolyhedronGrid& pg,
			VPMSM& vpm_sm,
			VPMPG& vpm_pg,
			std::tuple<IntersectionType, PGHalfedge, bool, bool> res)
		{
			if (std::get<3>(res))  // 交点是h_sm的终点
				nodes.add_new_node(get(vpm_sm, target(h_sm, sm)));
			else if (std::get<2>(res))  // 交点是h_sm的起点
				nodes.add_new_node(get(vpm_sm, source(h_sm, sm)));
			else
				nodes.add_new_node(h_sm, f_pg, sm, pg, vpm_sm, vpm_pg);
		}

		void fill_face_pair_to_node_case_face_sepf(SMHalfedge h_sm,
			PGFace f_pg,
			SurfaceMesh& sm,
			PolyhedronGrid& pg,
			NodeId node_id)
		{
			if (!is_border(h_sm, sm))
			{
				SMFace f_sm = face(h_sm, sm);
				f_to_node[std::make_pair(f_sm, f_pg)].insert(node_id);
			}
			h_sm = opposite(h_sm, sm);
			if (!is_border(h_sm, sm))
			{
				SMFace f_sm = face(h_sm, sm);
				f_to_node[std::make_pair(f_sm, f_pg)].insert(node_id);
			}
		}

		void fill_face_pair_to_node_case_edge_sepf(SMHalfedge h_sm,
			PGHalfedge h_pg,
			PGFaceSet* fset_pg,
			SurfaceMesh& sm,
			PolyhedronGrid& pg,
			NodeId node_id)
		{
			std::vector<PGFace> faces_around_edge;
			pg.incident_faces(pg.edge(h_pg), std::back_inserter(faces_around_edge));

			for (PGFace f_pg : faces_around_edge)
			{
				fill_face_pair_to_node_case_face_sepf(h_sm, f_pg, sm, pg, node_id);
				if (fset_pg != nullptr)
					fset_pg->erase(f_pg);
			}

			typename PGEdgeToSMFaces::iterator it_pge_smfs = pg_edge_to_sm_faces.find(pg.edge(h_pg));

			if (it_pge_smfs == pg_edge_to_sm_faces.end())
				return;
			else
			{
				SMFaceSet& fset_sm = it_pge_smfs->second;
				if (!is_border(h_sm, sm))
					fset_sm.erase(face(h_sm, sm));

				h_sm = opposite(h_sm, sm);
				if (!is_border(h_sm, sm))
					fset_sm.erase(face(h_sm, sm));
			}
		}

		void fill_face_pair_to_node_case_vertex_sepf(SMHalfedge h_sm,
			PGHalfedge h_pg,
			PGFaceSet* fset_pg,
			SurfaceMesh& sm,
			PolyhedronGrid& pg,
			NodeId node_id)
		{
			std::vector<PGEdge> edges_around_vertex;
			pg.incident_edges(pg.target(h_pg), std::back_inserter(edges_around_vertex));

			for (PGEdge e_pg : edges_around_vertex)
				fill_face_pair_to_node_case_edge_sepf(h_sm, pg.halfedge(e_pg), fset_pg, sm, pg, node_id);
		}

		// 计算表面网格的边与多面体网格的面的交点并记录相关的信息.
		void compute_intersection_points_sepf(SMEdgeToPGFaces& sm_edge_to_pg_faces,
			SurfaceMesh& sm,
			PolyhedronGrid& pg,
			VPMSM& vpm_sm,
			VPMPG& vpm_pg,
			NodeId& current_node)
		{
			/*
			* 在此解释一下std::tuple<IntersectionType, Halfedge, bool, bool>的含义.
			* tuple[0]: ON_FACE/ON_EDGE/ON_VERTEX, 说明相交的类型.
			* tuple[1]: 三角形面上的一条halfedge. 注意它和线段边无关.
			*			若相交类型为ON_FACE, halfedge是face的拓扑中记录的半边;
			*			若相交类型为ON_EDGE, halfedge是相交边对应的半边;
			*			若相交类型为ON_VERTEX, halfedge指向那个vertex;
			* tuple[2]: 线段边的起始点是否交在三角面上.
			* tuple[3]: 线段边的终点是否交在三角面上.
			*
			* 之所以有这个类型, 是因为后续要收集每个face/edge上的交点, 必须在tuple中返回1~3的信息, 方便通过拓扑查询记录相关信息.
			*/
			typedef std::tuple<IntersectionType, PGHalfedge, bool, bool>  InterInfo;

			// 遍历相交检测出的所有元素.
			for (typename SMEdgeToPGFaces::iterator it = sm_edge_to_pg_faces.begin();
				it != sm_edge_to_pg_faces.end();
				++it)
			{
				SMEdge e_sm = it->first;
				SMHalfedge h_sm = halfedge(e_sm, sm);
				PGFaceSet& fset_pg = it->second;

				// 处理每个e_sm干涉的面.
				while (!fset_pg.empty())
				{
					PGFace f_pg = *fset_pg.begin();
					// Note: e_sm和f_pg通过halfedge()方法取到的半边具有不确定性, 
					// 即h_sm的方向有两个, f_pg也有两个相反的旋向, 这导致计算相交类型的不确定性, 
					// 不过我们是以物理面/边为key记录信息的, 只需保证按照res记录的信息来进行后续的操作即可.
					InterInfo res = compute_intersection_type_sepf(h_sm, f_pg, sm, pg, vpm_sm, vpm_pg);
					IntersectionType type = std::get<0>(res);

					// 注意, 交点是在遍历edge_to_faces的过程中, 由边与三角形计算得出的,
					// 我们给每个交点分配一个NodeId, 唯一地标识它, 这样方便后续查找. 
					// 如果交点是边的一个端点, 这个交点会关联很多条边, 它们与三角形的交点都是这个交点.
					// 这样就会重复计算, 导致明明一个交点, 却分配了多个NodeId, 这就事与愿违了.
					// 所以需要从关联的边记录的FaceSet中, 删除这个面, 避免多次计算.
					// 
					// all_halfedges就是用来收集这些关联的边的.
					//
					std::vector<SMHalfedge> all_halfedges;
					if (std::get<3>(res)) // is edge target in triangle plane
					{
						std::copy(halfedges_around_target(h_sm, sm).first,
							halfedges_around_target(h_sm, sm).second,
							std::back_inserter(all_halfedges));
					}
					else if (std::get<2>(res))  // is edge source in triangle plane
					{
						std::copy(halfedges_around_source(h_sm, sm).first,
							halfedges_around_source(h_sm, sm).second,
							std::back_inserter(all_halfedges));
					}
					else
						all_halfedges.push_back(h_sm);

					assert(all_halfedges[0] == h_sm || all_halfedges[0] == opposite(h_sm, sm));

					typename std::vector<SMHalfedge>::iterator it_halfedge = all_halfedges.begin();
					switch (type)
					{
					case COPLANAR_TRIANGLES:
						// 我们不处理这种情况, 之前的步骤已经过滤掉了.
						/*std::cout << "Coplanar triangles: this point should never be reached!";*/
						break;

					case EMPTY:  // 不相交, 删掉接着处理下个面即可.
						fset_pg.erase(f_pg);
						break;

						// Case when the edge pierces the face in its interior.
					case ON_FACE: // 交在了面的内部.
					{
						assert(f_pg == pg.face(std::get<1>(res)));

						NodeId node_id = ++current_node;  // 分配node id
						// 第一类信息: 计算出交点坐标, 存放在nodes中.
						add_new_node_sepf(h_sm, f_pg, sm, pg, vpm_sm, vpm_pg, res);
						PGHalfedge h_pg = std::get<1>(res);

						// 第二类信息: 细化所需的信息, 由visitor保管.
						visitor.new_node_added_sepf(node_id, ON_FACE, h_sm, h_pg, sm, pg, std::get<2>(res), std::get<3>(res));

						// 第三类信息: polyline相关信息, 用于构建output_bulider需要的数据.
						for (; it_halfedge != all_halfedges.end(); ++it_halfedge)
						{
							fill_face_pair_to_node_case_face_sepf(*it_halfedge, f_pg, sm, pg, node_id);

							if (it_halfedge == all_halfedges.begin())
								fset_pg.erase(f_pg);
							else
							{
								// 从关联的边记录的FaceSet中, 删除这个面, 避免多次计算.
								typename SMEdgeToPGFaces::iterator it_sme_pgfs =
									sm_edge_to_pg_faces.find(edge(*it_halfedge, sm));
								if (it_sme_pgfs != sm_edge_to_pg_faces.end())
									it_sme_pgfs->second.erase(f_pg);
							}
						}
					}	// end case ON_FACE
					break;

					// Case when the edge intersect one edge of the face.
					case ON_EDGE:
					{
						NodeId node_id = ++current_node;
						// 第一类信息: 计算出交点坐标, 存放在nodes中.
						add_new_node_sepf(h_sm, f_pg, sm, pg, vpm_sm, vpm_pg, res);
						PGHalfedge h_pg = std::get<1>(res);

						// 第二类信息: 细化所需的信息, 由visitor保管.
						visitor.new_node_added_sepf(node_id, ON_EDGE, h_sm, h_pg, sm, pg, std::get<2>(res), std::get<3>(res));

						// 第三类信息: polyline相关信息, 用于构建output_bulider需要的数据.
						for (; it_halfedge != all_halfedges.end(); ++it_halfedge)
						{
							if (it_halfedge == all_halfedges.begin())
								fill_face_pair_to_node_case_edge_sepf(*it_halfedge, h_pg, &fset_pg, sm, pg, node_id);
							else
							{
								typename SMEdgeToPGFaces::iterator it_sme_pgfs =
									sm_edge_to_pg_faces.find(edge(*it_halfedge, sm));

								PGFaceSet* fset;
								if (it_sme_pgfs != sm_edge_to_pg_faces.end())
									fset = &(it_sme_pgfs->second);
								else
									fset = nullptr;
								fill_face_pair_to_node_case_edge_sepf(*it_halfedge, h_pg, fset, sm, pg, node_id);
							}
						}
					}	// end case ON_EDGE
					break;

					case ON_VERTEX:
					{
						NodeId node_id = ++current_node;
						PGHalfedge h_pg = std::get<1>(res);

						// 第一类信息: 计算出交点坐标, 存放在nodes中.
						nodes.add_new_node(vpm_pg[pg.target(h_pg)]);

						// 第二类信息: 细化所需的信息, 由visitor保管.
						visitor.new_node_added_sepf(node_id, ON_VERTEX, h_sm, h_pg, sm, pg, std::get<2>(res), std::get<3>(res));

						// 第三类信息: polyline相关信息, 用于构建output_bulider需要的数据.
						for (; it_halfedge != all_halfedges.end(); ++it_halfedge)
						{
							if (it_halfedge == all_halfedges.begin())
								fill_face_pair_to_node_case_vertex_sepf(*it_halfedge, h_pg, &fset_pg, sm, pg, node_id);
							else
							{
								typename SMEdgeToPGFaces::iterator it_sme_pgfs =
									sm_edge_to_pg_faces.find(edge(*it_halfedge, sm));

								PGFaceSet* fset;
								if (it_sme_pgfs != sm_edge_to_pg_faces.end())
									fset = &(it_sme_pgfs->second);
								else
									fset = nullptr;
								fill_face_pair_to_node_case_vertex_sepf(*it_halfedge, h_pg, fset, sm, pg, node_id);
							}
						}
					} // end case ON_VERTEX
					break;

					} // end switch on the type of the intersection
				} // end loop on all faces that intersect the edge
			} // end loop on all entries in 'sm_edge_to_pg_faces'

			assert(nodes.size() == unsigned(current_node + 1));
		}

		//================ 以下代码与上面差异不大, 可参考上面的注释 =======================

		std::tuple<IntersectionType, SMHalfedge, bool, bool>
			find_intersection_pesf(Point3& p, Point3& q,
				Point3& a, Point3& b, Point3& c,
				SMHalfedge h_sm,
				SurfaceMesh& sm,
				bool is_src_coplanar = false,
				bool is_tgt_coplanar = false)
		{
			typedef std::tuple<IntersectionType, SMHalfedge, bool, bool> result_type;

			CGAL::Orientation ab = CGAL::orientation(p, q, a, b);
			CGAL::Orientation bc = CGAL::orientation(p, q, b, c);
			CGAL::Orientation ca = CGAL::orientation(p, q, c, a);

			using CGAL::POSITIVE;
			using CGAL::NEGATIVE;
			using CGAL::COPLANAR;

			if (ab == POSITIVE || bc == POSITIVE || ca == POSITIVE)
				return result_type(EMPTY, GT::null_halfedge(), false, false);

			int nb_coplanar = (ab == COPLANAR ? 1 : 0) + (bc == COPLANAR ? 1 : 0) + (ca == COPLANAR ? 1 : 0);

			if (nb_coplanar == 0)
				return result_type(ON_FACE, h_sm, is_src_coplanar, is_tgt_coplanar);

			if (nb_coplanar == 1)
			{
				if (ab == COPLANAR)
					return result_type(ON_EDGE, next(h_sm, sm), is_src_coplanar, is_tgt_coplanar);
				if (bc == COPLANAR)
					return result_type(ON_EDGE, prev(h_sm, sm), is_src_coplanar, is_tgt_coplanar);

				assert(ca == COPLANAR);
				return result_type(ON_EDGE, h_sm, is_src_coplanar, is_tgt_coplanar);
			}

			assert(nb_coplanar == 2);

			if (ab != COPLANAR)
				return result_type(ON_VERTEX, prev(h_sm, sm), is_src_coplanar, is_tgt_coplanar);
			if (bc != COPLANAR)
				return result_type(ON_VERTEX, h_sm, is_src_coplanar, is_tgt_coplanar);
			assert(ca != COPLANAR);
			return result_type(ON_VERTEX, next(h_sm, sm), is_src_coplanar, is_tgt_coplanar);
		}

		std::tuple<IntersectionType, SMHalfedge, bool, bool>
			compute_intersection_type_pesf(PGHalfedge h_pg,
				SMFace f_sm,
				PolyhedronGrid& pg,
				SurfaceMesh& sm,
				VPMPG& vpm_pg,
				VPMSM& vpm_sm)
		{
			typedef std::tuple<IntersectionType, SMHalfedge, bool, bool>    result_type;
			typedef typename CGAL::Kernel_traits<Point3>::Kernel            Kernel;

			SMHalfedge h_sm = halfedge(f_sm, sm);

			Point3 a = get(vpm_sm, target(h_sm, sm));
			Point3 b = get(vpm_sm, target(next(h_sm, sm), sm));
			Point3 c = get(vpm_sm, source(h_sm, sm));
			Point3 p = vpm_pg[pg.source(h_pg)];
			Point3 q = vpm_pg[pg.target(h_pg)];

			CGAL::Orientation abcp = CGAL::orientation(a, b, c, p);
			CGAL::Orientation abcq = CGAL::orientation(a, b, c, q);

			switch (abcp)
			{
				using CGAL::POSITIVE;
				using CGAL::NEGATIVE;
				using CGAL::COPLANAR;
			case POSITIVE:
				switch (abcq)
				{
				case POSITIVE:
					return result_type(EMPTY, GT::null_halfedge(), false, false);
				case NEGATIVE:
					return find_intersection_pesf(p, q, a, b, c, h_sm, sm);
				case COPLANAR:
					return find_intersection_pesf(p, q, a, b, c, h_sm, sm, false, true);
				default:
					assert(false);
					return result_type(EMPTY, GT::null_halfedge(), false, false);
				}
			case NEGATIVE:
				switch (abcq)
				{
				case POSITIVE:
					return find_intersection_pesf(q, p, a, b, c, h_sm, sm);
				case NEGATIVE:
					return result_type(EMPTY, GT::null_halfedge(), false, false);
				case COPLANAR:
					return find_intersection_pesf(q, p, a, b, c, h_sm, sm, false, true);
				default:
					assert(false);
					return result_type(EMPTY, GT::null_halfedge(), false, true);
				}
			case COPLANAR:
				switch (abcq)
				{
				case POSITIVE:
					return find_intersection_pesf(q, p, a, b, c, h_sm, sm, true);
				case NEGATIVE:
					return find_intersection_pesf(p, q, a, b, c, h_sm, sm, true);
				case COPLANAR:
					if (::CGAL::Intersections::internal::do_intersect_coplanar(a, b, c, p, q, Kernel()))
						return result_type(COPLANAR_TRIANGLES, GT::null_halfedge(), true, true);
					return result_type(EMPTY, GT::null_halfedge(), true, true);
				default:
					assert(false);
					return result_type(EMPTY, GT::null_halfedge(), false, false);
				}
			default:
				assert(false);
				return result_type(EMPTY, GT::null_halfedge(), false, false);
			}
		}

		void add_new_node_pesf(PGHalfedge h_pg,
			SMFace f_sm,
			PolyhedronGrid& pg,
			SurfaceMesh& sm,
			VPMPG& vpm_pg,
			VPMSM& vpm_sm,
			std::tuple<IntersectionType, SMHalfedge, bool, bool> res)
		{
			if (std::get<3>(res))
				nodes.add_new_node(vpm_pg[pg.target(h_pg)]);
			else if (std::get<2>(res))
				nodes.add_new_node(vpm_pg[pg.source(h_pg)]);
			else
				nodes.add_new_node(h_pg, f_sm, pg, sm, vpm_pg, vpm_sm);
		}

		void fill_face_pair_to_node_case_face_pesf(PGEdge e_pg,
			SMFace f_sm,
			PolyhedronGrid& pg,
			SurfaceMesh& sm,
			NodeId node_id)
		{
			std::vector<PGFace> faces_around_edge;
			pg.incident_faces(e_pg, std::back_inserter(faces_around_edge));

			for (PGFace f_pg : faces_around_edge)
			{
				FacePair face_pair(f_sm, f_pg);
				f_to_node[face_pair].insert(node_id);
			}
		}

		void fill_face_pair_to_node_case_edge_pesf(PGEdge e_pg,
			SMHalfedge h_sm,
			SMFaceSet* fset_sm,
			PolyhedronGrid& pg,
			SurfaceMesh& sm,
			NodeId node_id)
		{
			std::vector<PGFace> faces_around_edge;
			pg.incident_faces(e_pg, std::back_inserter(faces_around_edge));

			for (PGFace f_pg : faces_around_edge)
			{
				if (!is_border(h_sm, sm))
				{
					SMFace f_sm = face(h_sm, sm);
					FacePair face_pair(f_sm, f_pg);
					f_to_node[face_pair].insert(node_id);

					if (fset_sm != nullptr)
						fset_sm->erase(f_sm);
				}
				h_sm = opposite(h_sm, sm);
				if (!is_border(h_sm, sm))
				{
					SMFace f_sm = face(h_sm, sm);
					FacePair face_pair(f_sm, f_pg);
					f_to_node[face_pair].insert(node_id);

					if (fset_sm != nullptr)
						fset_sm->erase(f_sm);
				}
			}

			typename SMEdgeToPGFaces::iterator it_sme_pgfs =
				sm_edge_to_pg_faces.find(edge(h_sm, sm));

			if (it_sme_pgfs == sm_edge_to_pg_faces.end())
				return;
			else
			{
				PGFaceSet& fset_pg = it_sme_pgfs->second;

				for (PGFace f_pg : faces_around_edge)
					fset_pg.erase(f_pg);
			}
		}

		void fill_face_pair_to_node_case_vertex_pesf(PGEdge e_pg,
			SMHalfedge h_sm,
			SMFaceSet* fset_sm,
			PolyhedronGrid& pg,
			SurfaceMesh& sm,
			NodeId node_id)
		{
			for (SMHalfedge h : halfedges_around_target(h_sm, sm))
				fill_face_pair_to_node_case_edge_pesf(e_pg, h, fset_sm, pg, sm, node_id);
		}

		void compute_intersection_points_pesf(PGEdgeToSMFaces& pg_edge_to_sm_faces,
			PolyhedronGrid& pg,
			SurfaceMesh& sm,
			VPMPG& vpm_pg,
			VPMSM& vpm_sm,
			NodeId& current_node)
		{
			typedef std::tuple<IntersectionType, SMHalfedge, bool, bool> InterInfo;

			for (typename PGEdgeToSMFaces::iterator it = pg_edge_to_sm_faces.begin();
				it != pg_edge_to_sm_faces.end();
				++it)
			{
				PGEdge e_pg = it->first;
				PGHalfedge h_pg = pg.halfedge(e_pg);
				SMFaceSet& fset_sm = it->second;

				while (!fset_sm.empty())
				{
					SMFace f_sm = *fset_sm.begin();

					InterInfo res = compute_intersection_type_pesf(h_pg, f_sm, pg, sm, vpm_pg, vpm_sm);
					IntersectionType type = std::get<0>(res);

					// 由于我们没有实现像halfedges_around_target这样的迭代器, 所以使用已有接口.
					// 而且多面体的incident halfedges of v是没有规律的, 不像表面网格可以转回来, 所以用edge更合适.
					std::vector<PGEdge> all_edges;
					if (std::get<3>(res))
						pg.incident_edges(pg.target(h_pg), std::back_inserter(all_edges));
					else if (std::get<2>(res))
						pg.incident_edges(pg.source(h_pg), std::back_inserter(all_edges));
					else
						all_edges.push_back(e_pg);

					typename std::vector<PGEdge>::iterator it_edge = all_edges.begin();
					switch (type)
					{
					case COPLANAR_TRIANGLES:
						std::cout << "Coplanar triangles: this point should never be reached!";
						break;

					case EMPTY:
						fset_sm.erase(f_sm);
						break;

					case ON_FACE:
					{
						assert(f_sm == face(std::get<1>(res), sm));

						NodeId node_id = ++current_node;
						add_new_node_pesf(h_pg, f_sm, pg, sm, vpm_pg, vpm_sm, res);
						SMHalfedge h_sm = std::get<1>(res);

						visitor.new_node_added_pesf(node_id, ON_FACE, h_pg, h_sm, pg, sm, std::get<2>(res), std::get<3>(res));
						for (; it_edge != all_edges.end(); ++it_edge)
						{
							fill_face_pair_to_node_case_face_pesf(*it_edge, f_sm, pg, sm, node_id);

							if (*it_edge == e_pg)
								fset_sm.erase(f_sm);
							else
							{
								typename PGEdgeToSMFaces::iterator it_pge_smfs =
									pg_edge_to_sm_faces.find(*it_edge);
								if (it_pge_smfs != pg_edge_to_sm_faces.end())
									it_pge_smfs->second.erase(f_sm);
							}
						}
					} // end case ON_FACE
					break;

					case ON_EDGE:
					{
						NodeId node_id = ++current_node;
						add_new_node_pesf(h_pg, f_sm, pg, sm, vpm_pg, vpm_sm, res);
						SMHalfedge h_sm = std::get<1>(res);

						visitor.new_node_added_pesf(node_id, ON_EDGE, h_pg, h_sm, pg, sm, std::get<2>(res), std::get<3>(res));
						for (; it_edge != all_edges.end(); ++it_edge)
						{
							if (*it_edge == e_pg)
								fill_face_pair_to_node_case_edge_pesf(*it_edge, h_sm, &fset_sm, pg, sm, node_id);
							else
							{
								typename PGEdgeToSMFaces::iterator it_pge_smfs =
									pg_edge_to_sm_faces.find(*it_edge);

								SMFaceSet* fset;
								if (it_pge_smfs != pg_edge_to_sm_faces.end())
									fset = &(it_pge_smfs->second);
								else
									fset = nullptr;
								fill_face_pair_to_node_case_edge_pesf(*it_edge, h_sm, fset, pg, sm, node_id);
							}
						}
					} // end case ON_EDGE
					break;

					case ON_VERTEX:
					{
						NodeId node_id = ++current_node;
						SMHalfedge h_sm = std::get<1>(res);
						nodes.add_new_node(get(vpm_sm, target(h_sm, sm)));

						visitor.new_node_added_pesf(node_id, ON_VERTEX, h_pg, h_sm, pg, sm, std::get<2>(res), std::get<3>(res));
						for (; it_edge != all_edges.end(); ++it_edge)
						{
							if (*it_edge == e_pg)
								fill_face_pair_to_node_case_vertex_pesf(*it_edge, h_sm, &fset_sm, pg, sm, node_id);
							else
							{
								typename PGEdgeToSMFaces::iterator it_pge_smfs =
									pg_edge_to_sm_faces.find(*it_edge);

								SMFaceSet* fset;
								if (it_pge_smfs != pg_edge_to_sm_faces.end())
									fset = &(it_pge_smfs->second);
								else
									fset = nullptr;
								fill_face_pair_to_node_case_vertex_pesf(*it_edge, h_sm, fset, pg, sm, node_id);
							}
						}
					} // end case ON_VERTEX
					break;
					} // end switch on the type of the intersection
				} // end loop on all faces that intersect the edge
			} // end loop on all entries (edges) in 'pg_edge_to_sm_faces'

			assert(nodes.size() == unsigned(current_node + 1));
		}

		struct GraphNode
		{
			boost::container::flat_set<NodeId> neighbors;
			unsigned degree;

			GraphNode() :degree(0) {}

			void insert(NodeId node_id)
			{
				++degree;
				assert(!neighbors.count(node_id));
				neighbors.insert(node_id);
			}

			// degree不减一.
			void erase(NodeId node_id)
			{
				assert(neighbors.count(node_id) != 0);
				neighbors.erase(node_id);
			}

			void make_terminal() { if (degree == 2) degree = 45; }
			bool is_terminal() const { return degree != 2; }
			bool empty() const { return neighbors.empty(); }
			NodeId top() const { return *neighbors.begin(); }
			void pop()
			{
				assert(!neighbors.empty());
				neighbors.erase(neighbors.begin());
			}
		};

		// CGAL构建交线基于两个事实: 
		// 1) 不共面的两个三角形对应最多两个交点.
		// 2) surface mesh一条边最多关联两个面, 所以一个交点的邻居最多两个.
		// 
		// 算法可以随意找一个交点i, 找到i的一个邻居j, 从j记录的相邻交点中删掉i, 由于j最多记录两个邻居,
		// 剩下的那个自然是polyline的下一个点. 如此循环下去, 就可以找到polyline的所有交点.
		// 
		// 其中一个网格变为多面体网格后, 事实2不再成立.
		// 当一条pg edge 交在sm face上时, 由于pg edge关联了多个面, 邻居不止两个.
		// 
		// 对于细化来说, 我们需要知道一个交点的所有邻居, 以方便插入约束边/交边. 这个很简单, 全部记录就好.
		// 但对于提取来说, 事实2不再成立造成很多麻烦. 
		// 
		// 一个解决方法是, 以cell为单位construct_polylines, 此时cell可看作一个surface mesh, 1和2都成立.
		//
		void construct_polylines(PolyhedronGrid& pg)
		{
			// 细化的信息.
			std::vector<boost::container::flat_set<NodeId>> node_neighbors(nodes.size());
			// 提取的信息.
			// Note: 此处记录cell的polyline信息不应使用std::map<cell, std::vector<GraphNode>>.
			//       因为只使用nodes中的很少一部分点, 如果node_id偏大, 使用vector会占用大量内存.
			typedef std::map<PGCell, std::map<NodeId, GraphNode>> GraphPerCell;
			GraphPerCell cell_graph;

			for (typename FacesToNodes::iterator it_fnode = f_to_node.begin();
				it_fnode != f_to_node.end();
				++it_fnode)
			{
				FacePair face_pair = it_fnode->first;
				SortedPair<NodeId> segment = it_fnode->second;
				assert(segment.dimension() == 2);

				NodeId i = segment.first;
				NodeId j = segment.second;
				node_neighbors[i].insert(j);
				node_neighbors[j].insert(i);

				// 记录每个cell上交点的相邻信息.
				PGFace f_pg = face_pair.second;
				std::vector<PGCell> cset;
				pg.incident_cells(f_pg, std::back_inserter(cset));
				for (PGCell c : cset)
				{
					cell_graph[c][i].insert(j);
					cell_graph[c][j].insert(i);

					visitor.record_segment_cell(segment, c);
				}
			}

			visitor.annotate_graph(node_neighbors);

			// 以cell为单位, 构建每个cell的polyline信息
			for (typename GraphPerCell::iterator it = cell_graph.begin();
				it != cell_graph.end();
				++it)
			{
				PGCell c = it->first;
				std::map<NodeId, GraphNode>& graph_nodes = it->second;

				std::size_t nb_nodes = graph_nodes.size();

				// dynamic_bitset要求索引是连续的, 但graph_nodes中的node id并不连续.
				// map an index to a node id.
				std::vector<NodeId> idx_to_node;
				// map a node id to an index;
				boost::container::flat_map<NodeId, std::size_t> node_to_idx;

				for (std::pair<NodeId, GraphNode> node : graph_nodes)
				{
					std::size_t i = idx_to_node.size();
					idx_to_node.push_back(node.first);
					node_to_idx.insert(std::make_pair(node.first, i));
				}

				boost::dynamic_bitset<> terminal_nodes(nb_nodes), interior_nodes(nb_nodes);
				for (std::size_t i = 0; i < nb_nodes; ++i)
				{
					if (graph_nodes[idx_to_node[i]].is_terminal())
						terminal_nodes.set(i);
					else
						interior_nodes.set(i);
				}

				// 先处理不闭合的交线.
				while (terminal_nodes.any())
				{
					// 随机找到不闭合交线的一个端点, 从它开始摸出一条交线.
					NodeId i = idx_to_node[terminal_nodes.find_first()];
					GraphNode& node_i = graph_nodes[i];

					NodeId j = node_i.top();
					visitor.start_new_polyline(c, i, j);
					node_i.pop();
					if (node_i.empty())
						terminal_nodes.reset(node_to_idx[i]);

					while (true)
					{
						GraphNode& node_j = graph_nodes[j];
						assert(!node_j.empty());
						node_j.erase(i);
						i = j;
						if (node_j.is_terminal() && node_j.empty())
						{
							terminal_nodes.reset(node_to_idx[j]);
							break;
						}
						else
						{
							j = node_j.top();
							visitor.add_node_to_polyline(c, j);
							node_j.pop();
							assert(node_j.empty());
							interior_nodes.reset(node_to_idx[i]);
						}
					}
				}

				// 处理闭合的交线.
				while (interior_nodes.any())
				{
					NodeId i = idx_to_node[interior_nodes.find_first()];
					GraphNode& node_i = graph_nodes[i];

					NodeId j = node_i.top();
					visitor.start_new_polyline(c, i, j);
					interior_nodes.reset(node_to_idx[i]);

					NodeId first = i;
					do
					{
						GraphNode& node_j = graph_nodes[j];
						interior_nodes.reset(node_to_idx[j]);
						node_j.erase(i);
						i = j;
						j = node_j.top();
						visitor.add_node_to_polyline(c, j);
					} while (j != first);
				}
			}
		}

	public:
		IntersectionImpl(SurfaceMesh& sm,
			PolyhedronGrid& pg,
			VPMSM& vpm_sm,
			VPMPG& vpm_pg,
			const IntersectionVisitor& v = IntersectionVisitor())
			: nodes(sm, pg, vpm_sm, vpm_pg)
			, visitor(v)
		{
		}

		// Intersection functor. 由仿函数组织相交部分的算法.
		void operator()()
		{
			SurfaceMesh& sm = nodes.sm;
			PolyhedronGrid& pg = nodes.pg;
			VPMSM& vpm_sm = nodes.vpm_sm;
			VPMPG& vpm_pg = nodes.vpm_pg;

			// Step1 相交检测.
			filter_intersections_sepf(sm, pg, vpm_sm, vpm_pg);
			filter_intersections_pesf(pg, sm, vpm_pg, vpm_sm);
			//filter_coplanar_case(sm, pg, vpm_sm, vpm_pg);

			std::cout << "相交检测过程-----------------------------" << std::endl;
			std::cout << "精确计算的次数：" << iExactComputeCount << std::endl;
			std::cout << "区间代数的次数：" << iIntervalComputeCount << std::endl;
			std::cout << "----------------------------------------\n\n";


			std::cout << "筛选出的sm_edge数量：" << sm_edge_to_pg_faces.size() << std::endl;
			std::unordered_set<PGFace> pg_face_set;
			std::size_t cnt = 0;
			for (auto it = sm_edge_to_pg_faces.begin(); it != sm_edge_to_pg_faces.end(); ++it)
			{
				pg_face_set.insert(it->second.begin(), it->second.end());
				cnt += it->second.size();
			}
			std::cout << "筛选出的pg_face数量（没有去重）：" << cnt << std::endl;
			std::cout << "筛选出的pg_face数量（去重后）：" << pg_face_set.size() << "\n\n";


			std::cout << "筛选出的pg_edge数量：" << pg_edge_to_sm_faces.size() << std::endl;
			std::unordered_set<SMFace> sm_face_set;
			cnt = 0;
			for (auto it = pg_edge_to_sm_faces.begin(); it != pg_edge_to_sm_faces.end(); ++it)
			{
				sm_face_set.insert(it->second.begin(), it->second.end());
				cnt += it->second.size();
			}
			std::cout << "筛选出的sm_face数量（没有去重）：" << cnt << std::endl;
			std::cout << "筛选出的sm_face数量（去重）：" << sm_face_set.size() << "\n\n";

			// Step2 计算交点并记录相关信息.
			NodeId current_node(std::numeric_limits<NodeId>::max());
			assert(current_node + 1 == 0);

			int cur_exact_cnt = iExactComputeCount;
			int cur_interval_cnt = iIntervalComputeCount;

			compute_intersection_points_sepf(sm_edge_to_pg_faces, sm, pg, vpm_sm, vpm_pg, current_node);
			compute_intersection_points_pesf(pg_edge_to_sm_faces, pg, sm, vpm_pg, vpm_sm, current_node);

			std::cout << "计算交点过程-----------------------------" << std::endl;
			std::cout << "精确计算的次数：" << iExactComputeCount - cur_exact_cnt << std::endl;
			std::cout << "区间代数的次数：" << iIntervalComputeCount - cur_interval_cnt << std::endl;
			std::cout << "----------------------------------------\n\n";

			// Step3 构建提取所需的交线信息.
			construct_polylines(pg);

			// 转入算法的第二部分, 交给visitor, 由其负责对网格进行细化
			visitor.finalize(nodes, sm, pg, vpm_sm, vpm_pg);
		}

	};

}	// namespace MCAL


#endif