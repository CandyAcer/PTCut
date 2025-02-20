// $License: GPL-3.0
// $Author: Liu Zhiqiang(from Beihang University)
// $Email: lzhq0930@gmail.com(liuzhiqiang@buaa.edu.cn)
// 
// Part 1 of PTCut, compute intersection between surface and polyhedral mesh.
// 


#ifndef MCAL_ALGO_INTERSECTION_H
#define MCAL_ALGO_INTERSECTION_H

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

#include "util/sorted_pair.h"
#include "intersection_nodes.h"
#include "intersection_callback.h"
#include "intersection_type.h"


namespace MCAL    // Mesh Cut Algorithm Libaray, 存放多面体网格切割算法的所有代码
{

// 阅读代码需要注意的一个地方:
// filter_intersections(), compute_intersection_points()等函数都会有两个版本, 
// 这是由于polyhedral mesh和surface mesh的行为不同导致的, 
// 为了区分, 我们在函数后面加了后缀
// 
// sepf: sm edge to pm face
// pesf: pm edge to sm face
//


/*
* Intersection是切割算法的第一部分, 大体流程如下:
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
* @param Corefine: 负责对网格进行细化的类, 作为本类的Visitor.
*/
template <typename SurfaceMesh,
		  typename PolyhedralMesh,
		  typename VPMSM,
		  typename VPMPM,
		  typename Corefine>
class Intersection
{
	// typedefs
	typedef boost::graph_traits<SurfaceMesh>             GT;
	typedef typename GT::vertex_descriptor               SMVertex;
	typedef typename GT::edge_descriptor                 SMEdge;
	typedef typename GT::halfedge_descriptor             SMHalfedge;
	typedef typename GT::face_descriptor                 SMFace;

	typedef typename PolyhedralMesh::VertexIndex         PMVertex;
	typedef typename PolyhedralMesh::EdgeIndex           PMEdge;
	typedef typename PolyhedralMesh::FaceIndex           PMFace;
	typedef typename PolyhedralMesh::CellIndex           PMCell;
	typedef typename PolyhedralMesh::Halfedge            PMHalfedge;

	typedef typename boost::property_traits<VPMSM>::value_type  Point3;
	static_assert((std::is_same<Point3, typename PolyhedralMesh::Point>::value),
				   "输入网格的点类型不一致");

	// Axis aligned bounding box (with info).
	typedef BoxInfo<SurfaceMesh, PolyhedralMesh>         Info;
	typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, Info>  Box;

	typedef std::unordered_set<PMFace> PMFaceSet;
	typedef std::unordered_map<SMEdge, PMFaceSet> SMEdgeToPMFaces;
	typedef std::unordered_set<SMFace> SMFaceSet;
	typedef std::unordered_map<PMEdge, SMFaceSet> PMEdgeToSMFaces;

	typedef IntersectionNodes<SurfaceMesh, PolyhedralMesh, VPMSM, VPMPM> NodeVector;

	typedef std::size_t NodeId;

	// adjacent information.
	typedef std::pair<SMFace, PMFace> FacePair;
	typedef std::unordered_map<FacePair, SortedPair<NodeId>, boost::hash<FacePair>> FacesToNodes;

	typedef std::set<FacePair> CoplanarFaceSet;

	// data members
	SMEdgeToPMFaces sm_edge_to_pm_faces;  // intersected surface mesh edge and polyhedral mesh faces.
	PMEdgeToSMFaces pm_edge_to_sm_faces;  // intersected polyhedral mesh edge and surface mesh faces.
	CoplanarFaceSet coplanar_faces;  // coplanar surface mesh face and polyhedral mesh face. we don't handle their intersections.
	NodeVector nodes; // manager all intersection nodes, store exact points to offer exact predications.
	Corefine visitor;   // intersection visitor to collect information used in triangulation.
	FacesToNodes f_to_node;   // Associate a pair of triangles to their intersection points.

	// member functions
private:

	// 检测surface mesh的edge与哪些polyhedral mesh的face相交(指包围盒相交, 非实际相交).
	// 过滤掉surface mesh中包围盒不干涉的edge.
	void filter_intersections_sepf(SurfaceMesh& sm,
								   PolyhedralMesh& pm,
								   VPMSM& vpm_sm,
								   VPMPM& vpm_pm)
	{
		std::vector<Box> sm_edge_boxes;
		std::vector<Box> pm_face_boxes;

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

		// 为polyhedral mesh的所有face套一个AABB包围盒.
		pm_face_boxes.reserve(pm.num_faces());
		for (PMFace face : pm.faces())
		{
			std::vector<PMVertex> vset;
			pm.incident_vertices(face, std::back_inserter(vset));
			CGAL::Bbox_3 face_bbox;
			for (PMVertex v_pm : vset)
				face_bbox += vpm_pm[v_pm].bbox();

			PMHalfedge h_pm = pm.halfedge(face);
			pm_face_boxes.push_back(Box(face_bbox, Info(h_pm)));
		}

		// callback在检测到包围盒相交时进行记录.
		typedef CollectSMEdgeToPMFaces<SurfaceMesh, PolyhedralMesh, SMEdgeToPMFaces>
			Callback;
		Callback callback(sm, pm, sm_edge_to_pm_faces);

		// 调用相交检测算法.
		CGAL::box_intersection_d(sm_edge_boxes.begin(), sm_edge_boxes.end(),
								 pm_face_boxes.begin(), pm_face_boxes.end(),
								 callback);

		// 由于多面体网格的面是空间多边形, 相较于三角形更为复杂, 需要在此将检测到的面进行三角化.
		std::unordered_map<PMFace, std::vector<PMFace>> polygon_to_triangles;
		for (typename SMEdgeToPMFaces::iterator it = sm_edge_to_pm_faces.begin();
												it != sm_edge_to_pm_faces.end();
												++it)
		{
			SMEdge e_sm = it->first;
			// 此处不能是引用, 由于后续会往set中添加元素, 引用会导致set中动态增加, 
			// 我们要遍历的就是原来的那些面.
			PMFaceSet polygon_set = it->second;

			for (PMFace polygon : polygon_set)
			{
				// 对干涉的多面体网格面进行三角化, 
				// 三角化后的面记录在传入的容器中, 其首元素一定是polygon.
				std::vector<PMFace> refined_triangles;
				if (!polygon_to_triangles.count(polygon))
				{
					pm.triangulate_a_polygon(polygon, std::back_inserter(refined_triangles));
					assert(polygon == refined_triangles.front());
					polygon_to_triangles.insert(std::make_pair(polygon, refined_triangles));
				}
				else
					refined_triangles = polygon_to_triangles[polygon];

				// 记录好三角化后的所有面.
				sm_edge_to_pm_faces[e_sm].insert(refined_triangles.begin(), refined_triangles.end());
			}
		}
	}

	// 检测polyhedral mesh的edge与哪些Surface mesh的face相交(指包围盒相交, 非实际相交).
	// 过滤掉polyhedral mesh中包围盒不干涉的edge.
	void filter_intersections_pesf(PolyhedralMesh& pm,
								   SurfaceMesh& sm,
								   VPMPM& vpm_pm,
								   VPMSM& vpm_sm)
	{
		std::vector<Box> pm_cell_boxes;
		std::vector<Box> sm_face_boxes;

		// 不能为edge套包围盒, 因为它关联的面可能会存在非三角形面, 这样判断共面就会很麻烦.
		// 既然如此, 是否可以调用box_intersection_d筛出可能相交的edge, 然后三角化其关联的面?
		// 可以. 但这会导致多一次filter_pm_intersections(), 效率很低.
		// 
		// Note: 
		// 上面的思路三角化后新产生的edge包围盒可能更大(如果是斜边的话), 但不会有错误, 
		// 因为新产生的edge一定在多边形的内部, 它关联的面都是三角形面, 不会陷入死循环.

		// 为polyhedral mesh的所有cell套上AABB包围盒.
		pm_cell_boxes.reserve(pm.num_cells());
		for (PMCell c : pm.cells())
		{
			std::vector<PMVertex> vset;
			pm.incident_vertices(c, std::back_inserter(vset));
			CGAL::Bbox_3 cell_bbox;
			for (PMVertex v_pm : vset)
				cell_bbox += vpm_pm[v_pm].bbox();

			PMHalfedge h_pm = pm.halfedge(c);
			pm_cell_boxes.push_back(Box(cell_bbox, Info(h_pm)));
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
		typedef std::unordered_map<PMCell, std::unordered_set<SMFace>> PMCellToSMFaces;
		PMCellToSMFaces pm_cell_to_sm_faces;

		typedef CollectPMCellToSMFaces<PolyhedralMesh, SurfaceMesh, PMCellToSMFaces>
			Callback;
		Callback callback(pm, sm, pm_cell_to_sm_faces);

		// 调用相交检测算法.
		CGAL::box_intersection_d(pm_cell_boxes.begin(), pm_cell_boxes.end(),
								 sm_face_boxes.begin(), sm_face_boxes.end(),
								 callback);

		// 之所以使用多面体网格的单元而不是边进行相交检测, 是因为它的包围盒比包含的所有边都大,
		for (typename PMCellToSMFaces::iterator it = pm_cell_to_sm_faces.begin();
												it != pm_cell_to_sm_faces.end();
												++it)
		{
			PMCell c = it->first;
			SMFaceSet& fset_sm = it->second;  // 此处可以是引用, 因为不会insert, size固定不变.

			std::vector<PMFace> cell_faces;
			pm.incident_faces(c, std::back_inserter(cell_faces));

			for (PMFace polygon : cell_faces)
				pm.triangulate_a_polygon(polygon, CGAL::Emptyset_iterator());

			// 收集多边形polygon三角化后的所有edge.
			std::vector<PMEdge> all_edges;
			pm.incident_edges(c, std::back_inserter(all_edges));
			// 对于三角化后的所有边, 由于边的包围盒一般比cell的小,
			// 所以不一定与fset_sm中的所有表面网格面干涉, 需要剔除.
			for (PMEdge e_pm : all_edges)
			{
				PMHalfedge h_pm = pm.halfedge(e_pm);
				CGAL::Bbox_3 edge_bbox(
					vpm_pm[pm.source(h_pm)].bbox() +
					vpm_pm[pm.target(h_pm)].bbox());

				for (SMFace f_sm : fset_sm)
				{
					SMHalfedge h_sm = halfedge(f_sm, sm);
					CGAL::Bbox_3 face_bbox(
						get(vpm_sm, source(h_sm, sm)).bbox() +
						get(vpm_sm, target(h_sm, sm)).bbox() +
						get(vpm_sm, target(next(h_sm, sm), sm)).bbox());

					// 若pm edge的包围盒和sm face的包围盒确实干涉, 记录.
					if (CGAL::do_overlap(edge_bbox, face_bbox))
						pm_edge_to_sm_faces[e_pm].insert(f_sm);
				}
			}
		}
	}

	std::tuple<IntersectionType, PMHalfedge, bool, bool>
	find_intersection_sepf(Point3& p, Point3& q,  // segment
						   Point3& a, Point3& b, Point3& c,   // triangle
						   PMHalfedge h_pm, // halfedge of the triangle face, its target is a
						   PolyhedralMesh& pm,
						   bool is_src_coplanar = false,
						   bool is_tgt_coplanar = false)
	{
		typedef std::tuple<IntersectionType, PMHalfedge, bool, bool>     result_type;

		CGAL::Orientation ab = CGAL::orientation(p, q, a, b);
		CGAL::Orientation bc = CGAL::orientation(p, q, b, c);
		CGAL::Orientation ca = CGAL::orientation(p, q, c, a);

		using CGAL::POSITIVE;
		using CGAL::NEGATIVE;
		using CGAL::COPLANAR;

		if (ab == POSITIVE || bc == POSITIVE || ca == POSITIVE)
			return result_type(EMPTY, pm.null_halfedge(), false, false);

		int nb_coplanar = (ab == COPLANAR ? 1 : 0) + (bc == COPLANAR ? 1 : 0) + (ca == COPLANAR ? 1 : 0);

		if (nb_coplanar == 0)
			return result_type(ON_FACE, h_pm, is_src_coplanar, is_tgt_coplanar);

		if (nb_coplanar == 1)
		{
			if (ab == COPLANAR)      // intersection is ab
				return result_type(ON_EDGE, pm.next(h_pm), is_src_coplanar, is_tgt_coplanar);
			if (bc == COPLANAR)      // intersection is bc
				return result_type(ON_EDGE, pm.prev(h_pm), is_src_coplanar, is_tgt_coplanar);
			assert(ca == COPLANAR);
			// intersection is ca
			return result_type(ON_EDGE, h_pm, is_src_coplanar, is_tgt_coplanar);
		}

		assert(nb_coplanar == 2);

		if (ab != COPLANAR)		// intersection is c
			return result_type(ON_VERTEX, pm.prev(h_pm), is_src_coplanar, is_tgt_coplanar);
		if (bc != COPLANAR)		// intersection is a
			return result_type(ON_VERTEX, h_pm, is_src_coplanar, is_tgt_coplanar);
		assert(ca != COPLANAR);
		// intersection is b
		return result_type(ON_VERTEX, pm.next(h_pm), is_src_coplanar, is_tgt_coplanar);
	}

	std::tuple<IntersectionType, PMHalfedge, bool, bool>
	compute_intersection_type_sepf(SMHalfedge h_sm,
								   PMFace f_pm,
								   SurfaceMesh& sm,
								   PolyhedralMesh& pm,
								   VPMSM& vpm_sm,
								   VPMPM& vpm_pm)
	{
		typedef std::tuple<IntersectionType, PMHalfedge, bool, bool>     result_type;
		typedef typename CGAL::Kernel_traits<Point3>::Kernel             Kernel;

		assert(pm.is_triangle(f_pm));
		PMHalfedge h_pm = pm.halfedge(f_pm);

		Point3 a = vpm_pm[pm.target(h_pm)];
		Point3 b = vpm_pm[pm.target(pm.next(h_pm))];
		Point3 c = vpm_pm[pm.source(h_pm)];
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
				return result_type(EMPTY, pm.null_halfedge(), false, false);
			case NEGATIVE:
				// p sees the triangle in counterclockwise order
				return find_intersection_sepf(p, q, a, b, c, h_pm, pm);
			case COPLANAR:
				// q belongs to the triangle's supporting plane
				// p sees the triangle in counterclockwise order
				return find_intersection_sepf(p, q, a, b, c, h_pm, pm, false, true);
			default: // should not happen.
				assert(false);
				return result_type(EMPTY, pm.null_halfedge(), false, false);
			}
		case NEGATIVE:
			switch (abcq)
			{
			case POSITIVE:
				// q sees the triangle in counterclockwise order
				return find_intersection_sepf(q, p, a, b, c, h_pm, pm);
			case NEGATIVE:
				// the segment lies in the negative open halfspaces defined by the
				// triangle's supporting plane
				return result_type(EMPTY, pm.null_halfedge(), false, false);
			case COPLANAR:
				// q belongs to the triangle's supporting plane
				// p sees the triangle in clockwise order
				return find_intersection_sepf(q, p, a, b, c, h_pm, pm, false, true);
			default:
				assert(false);
				return result_type(EMPTY, pm.null_halfedge(), false, true);
			}
		case COPLANAR: // p belongs to the triangle's supporting plane
			switch (abcq)
			{
			case POSITIVE:
				// q sees the triangle in counterclockwise order
				return find_intersection_sepf(q, p, a, b, c, h_pm, pm, true);
			case NEGATIVE:
				// q sees the triangle in clockwise order
				return find_intersection_sepf(p, q, a, b, c, h_pm, pm, true);
			case COPLANAR:
				// the segment is coplanar with the triangle's supporting plane
				// we test whether the segment intersects the triangle in the common
				// supporting plane
				if (::CGAL::Intersections::internal::do_intersect_coplanar(a, b, c, p, q, Kernel()))
					return result_type(COPLANAR_TRIANGLES, pm.null_halfedge(), true, true);
				return result_type(EMPTY, pm.null_halfedge(), true, true);

			default: // should not happen.
				assert(false);
				return result_type(EMPTY, pm.null_halfedge(), false, false);
			}
		default: // should not happen.
			assert(false);
			return result_type(EMPTY, pm.null_halfedge(), false, false);
		}
	}

	// 计算sm edge和pm face的交点坐标, 以精确形式记录在nodes中, 由node id唯一标识.
	void add_new_node_sepf(SMHalfedge h_sm,
						   PMFace f_pm,
						   SurfaceMesh& sm,
						   PolyhedralMesh& pm,
						   VPMSM& vpm_sm,
						   VPMPM& vpm_pm,
						   std::tuple<IntersectionType, PMHalfedge, bool, bool> res)
	{
		if (std::get<3>(res))  // 交点是h_sm的终点
			nodes.add_new_node(get(vpm_sm, target(h_sm, sm)));
		else if (std::get<2>(res))  // 交点是h_sm的起点
			nodes.add_new_node(get(vpm_sm, source(h_sm, sm)));
		else
			nodes.add_new_node(h_sm, f_pm, sm, pm, vpm_sm, vpm_pm);
	}

	void fill_face_pair_to_node_case_face_sepf(SMHalfedge h_sm,
											   PMFace f_pm,
											   SurfaceMesh& sm,
											   PolyhedralMesh& pm,
											   NodeId node_id)
	{
		if (!is_border(h_sm, sm))
		{
			SMFace f_sm = face(h_sm, sm);
			f_to_node[std::make_pair(f_sm, f_pm)].insert(node_id);
		}
		h_sm = opposite(h_sm, sm);
		if (!is_border(h_sm, sm))
		{
			SMFace f_sm = face(h_sm, sm);
			f_to_node[std::make_pair(f_sm, f_pm)].insert(node_id);
		}
	}

	void fill_face_pair_to_node_case_edge_sepf(SMHalfedge h_sm,
											   PMHalfedge h_pm,
											   PMFaceSet* fset_pm,
											   SurfaceMesh& sm,
											   PolyhedralMesh& pm,
											   NodeId node_id)
	{
		std::vector<PMFace> faces_around_edge;
		pm.incident_faces(pm.edge(h_pm), std::back_inserter(faces_around_edge));

		for (PMFace f_pm : faces_around_edge)
		{
			fill_face_pair_to_node_case_face_sepf(h_sm, f_pm, sm, pm, node_id);
			if (fset_pm != nullptr)
				fset_pm->erase(f_pm);
		}

		typename PMEdgeToSMFaces::iterator it_pme_smfs = pm_edge_to_sm_faces.find(pm.edge(h_pm));

		if (it_pme_smfs == pm_edge_to_sm_faces.end())
			return;
		else
		{
			SMFaceSet& fset_sm = it_pme_smfs->second;
			if (!is_border(h_sm, sm))
				fset_sm.erase(face(h_sm, sm));

			h_sm = opposite(h_sm, sm);
			if (!is_border(h_sm, sm))
				fset_sm.erase(face(h_sm, sm));
		}
	}

	void fill_face_pair_to_node_case_vertex_sepf(SMHalfedge h_sm,
												 PMHalfedge h_pm,
												 PMFaceSet* fset_pm,
												 SurfaceMesh& sm,
												 PolyhedralMesh& pm,
												 NodeId node_id)
	{
		std::vector<PMEdge> edges_around_vertex;
		pm.incident_edges(pm.target(h_pm), std::back_inserter(edges_around_vertex));

		for (PMEdge e_pm : edges_around_vertex)
			fill_face_pair_to_node_case_edge_sepf(h_sm, pm.halfedge(e_pm), fset_pm, sm, pm, node_id);
	}

	// 计算表面网格的边与多面体网格的面的交点并记录相关的信息.
	void compute_intersection_points_sepf(SMEdgeToPMFaces& sm_edge_to_pm_faces,
										  SurfaceMesh& sm,
										  PolyhedralMesh& pm,
										  VPMSM& vpm_sm,
										  VPMPM& vpm_pm,
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
		typedef std::tuple<IntersectionType, PMHalfedge, bool, bool>  InterInfo;

		// 遍历相交检测出的所有元素.
		for (typename SMEdgeToPMFaces::iterator it = sm_edge_to_pm_faces.begin();
												it != sm_edge_to_pm_faces.end();
												++it)
		{
			SMEdge e_sm = it->first;
			SMHalfedge h_sm = halfedge(e_sm, sm);
			PMFaceSet& fset_pm = it->second;

			// 处理每个e_sm干涉的面.
			while (!fset_pm.empty())
			{
				PMFace f_pm = *fset_pm.begin();
				// Note: e_sm和f_pm通过halfedge()方法取到的半边具有不确定性, 
				// 即h_sm的方向有两个, f_pm也有两个相反的旋向, 这导致计算相交类型的不确定性, 
				// 不过我们是以物理面/边为key记录信息的, 只需保证按照res记录的信息来进行后续的操作即可.
				InterInfo res = compute_intersection_type_sepf(h_sm, f_pm, sm, pm, vpm_sm, vpm_pm);
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
					fset_pm.erase(f_pm);
					break;

					// Case when the edge pierces the face in its interior.
				case ON_FACE: // 交在了面的内部.
				{
					assert(f_pm == pm.face(std::get<1>(res)));

					NodeId node_id = ++current_node;  // 分配node id
					// 第一类信息: 计算出交点坐标, 存放在nodes中.
					add_new_node_sepf(h_sm, f_pm, sm, pm, vpm_sm, vpm_pm, res);
					PMHalfedge h_pm = std::get<1>(res);

					// 第二类信息: 细化所需的信息, 由visitor保管.
					visitor.new_node_added_sepf(node_id, ON_FACE, h_sm, h_pm, sm, pm, std::get<2>(res), std::get<3>(res));

					// 第三类信息: polyline相关信息, 用于构建output_bulider需要的数据.
					for (; it_halfedge != all_halfedges.end(); ++it_halfedge)
					{
						fill_face_pair_to_node_case_face_sepf(*it_halfedge, f_pm, sm, pm, node_id);

						if (it_halfedge == all_halfedges.begin())
							fset_pm.erase(f_pm);
						else
						{
							// 从关联的边记录的FaceSet中, 删除这个面, 避免多次计算.
							typename SMEdgeToPMFaces::iterator it_sme_pmfs =
								sm_edge_to_pm_faces.find(edge(*it_halfedge, sm));
							if (it_sme_pmfs != sm_edge_to_pm_faces.end())
								it_sme_pmfs->second.erase(f_pm);
						}
					}
				}	// end case ON_FACE
				break;

				// Case when the edge intersect one edge of the face.
				case ON_EDGE:
				{
					NodeId node_id = ++current_node;
					// 第一类信息: 计算出交点坐标, 存放在nodes中.
					add_new_node_sepf(h_sm, f_pm, sm, pm, vpm_sm, vpm_pm, res);
					PMHalfedge h_pm = std::get<1>(res);

					// 第二类信息: 细化所需的信息, 由visitor保管.
					visitor.new_node_added_sepf(node_id, ON_EDGE, h_sm, h_pm, sm, pm, std::get<2>(res), std::get<3>(res));

					// 第三类信息: polyline相关信息, 用于构建output_bulider需要的数据.
					for (; it_halfedge != all_halfedges.end(); ++it_halfedge)
					{
						if (it_halfedge == all_halfedges.begin())
							fill_face_pair_to_node_case_edge_sepf(*it_halfedge, h_pm, &fset_pm, sm, pm, node_id);
						else
						{
							typename SMEdgeToPMFaces::iterator it_sme_pmfs =
								sm_edge_to_pm_faces.find(edge(*it_halfedge, sm));

							PMFaceSet* fset;
							if (it_sme_pmfs != sm_edge_to_pm_faces.end())
								fset = &(it_sme_pmfs->second);
							else
								fset = nullptr;
							fill_face_pair_to_node_case_edge_sepf(*it_halfedge, h_pm, fset, sm, pm, node_id);
						}
					}
				}	// end case ON_EDGE
				break;

				case ON_VERTEX:
				{
					NodeId node_id = ++current_node;
					PMHalfedge h_pm = std::get<1>(res);

					// 第一类信息: 计算出交点坐标, 存放在nodes中.
					nodes.add_new_node(vpm_pm[pm.target(h_pm)]);

					// 第二类信息: 细化所需的信息, 由visitor保管.
					visitor.new_node_added_sepf(node_id, ON_VERTEX, h_sm, h_pm, sm, pm, std::get<2>(res), std::get<3>(res));

					// 第三类信息: polyline相关信息, 用于构建output_bulider需要的数据.
					for (; it_halfedge != all_halfedges.end(); ++it_halfedge)
					{
						if (it_halfedge == all_halfedges.begin())
							fill_face_pair_to_node_case_vertex_sepf(*it_halfedge, h_pm, &fset_pm, sm, pm, node_id);
						else
						{
							typename SMEdgeToPMFaces::iterator it_sme_pmfs =
								sm_edge_to_pm_faces.find(edge(*it_halfedge, sm));

							PMFaceSet* fset;
							if (it_sme_pmfs != sm_edge_to_pm_faces.end())
								fset = &(it_sme_pmfs->second);
							else
								fset = nullptr;
							fill_face_pair_to_node_case_vertex_sepf(*it_halfedge, h_pm, fset, sm, pm, node_id);
						}
					}
				} // end case ON_VERTEX
				break;

				} // end switch on the type of the intersection
			} // end loop on all faces that intersect the edge
		} // end loop on all entries in 'sm_edge_to_pm_faces'

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
	compute_intersection_type_pesf(PMHalfedge h_pm,
								   SMFace f_sm,
								   PolyhedralMesh& pm,
								   SurfaceMesh& sm,
								   VPMPM& vpm_pm,
								   VPMSM& vpm_sm)
	{
		typedef std::tuple<IntersectionType, SMHalfedge, bool, bool>    result_type;
		typedef typename CGAL::Kernel_traits<Point3>::Kernel            Kernel;

		SMHalfedge h_sm = halfedge(f_sm, sm);

		Point3 a = get(vpm_sm, target(h_sm, sm));
		Point3 b = get(vpm_sm, target(next(h_sm, sm), sm));
		Point3 c = get(vpm_sm, source(h_sm, sm));
		Point3 p = vpm_pm[pm.source(h_pm)];
		Point3 q = vpm_pm[pm.target(h_pm)];

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

	void add_new_node_pesf(PMHalfedge h_pm,
						   SMFace f_sm,
						   PolyhedralMesh& pm,
						   SurfaceMesh& sm,
						   VPMPM& vpm_pm,
						   VPMSM& vpm_sm,
						   std::tuple<IntersectionType, SMHalfedge, bool, bool> res)
	{
		if (std::get<3>(res))
			nodes.add_new_node(vpm_pm[pm.target(h_pm)]);
		else if (std::get<2>(res))
			nodes.add_new_node(vpm_pm[pm.source(h_pm)]);
		else
			nodes.add_new_node(h_pm, f_sm, pm, sm, vpm_pm, vpm_sm);
	}

	void fill_face_pair_to_node_case_face_pesf(PMEdge e_pm,
											   SMFace f_sm,
											   PolyhedralMesh& pm,
											   SurfaceMesh& sm,
											   NodeId node_id)
	{
		std::vector<PMFace> faces_around_edge;
		pm.incident_faces(e_pm, std::back_inserter(faces_around_edge));

		for (PMFace f_pm : faces_around_edge)
		{
			FacePair face_pair(f_sm, f_pm);
			f_to_node[face_pair].insert(node_id);
		}
	}

	void fill_face_pair_to_node_case_edge_pesf(PMEdge e_pm,
											   SMHalfedge h_sm,
											   SMFaceSet* fset_sm,
											   PolyhedralMesh& pm,
											   SurfaceMesh& sm,
											   NodeId node_id)
	{
		std::vector<PMFace> faces_around_edge;
		pm.incident_faces(e_pm, std::back_inserter(faces_around_edge));

		for (PMFace f_pm : faces_around_edge)
		{
			if (!is_border(h_sm, sm))
			{
				SMFace f_sm = face(h_sm, sm);
				FacePair face_pair(f_sm, f_pm);
				f_to_node[face_pair].insert(node_id);

				if (fset_sm != nullptr)
					fset_sm->erase(f_sm);
			}
			h_sm = opposite(h_sm, sm);
			if (!is_border(h_sm, sm))
			{
				SMFace f_sm = face(h_sm, sm);
				FacePair face_pair(f_sm, f_pm);
				f_to_node[face_pair].insert(node_id);

				if (fset_sm != nullptr)
					fset_sm->erase(f_sm);
			}
		}

		typename SMEdgeToPMFaces::iterator it_sme_pmfs =
			sm_edge_to_pm_faces.find(edge(h_sm, sm));

		if (it_sme_pmfs == sm_edge_to_pm_faces.end())
			return;
		else
		{
			PMFaceSet& fset_pm = it_sme_pmfs->second;

			for (PMFace f_pm : faces_around_edge)
				fset_pm.erase(f_pm);
		}
	}

	void fill_face_pair_to_node_case_vertex_pesf(PMEdge e_pm,
												 SMHalfedge h_sm,
												 SMFaceSet* fset_sm,
												 PolyhedralMesh& pm,
												 SurfaceMesh& sm,
												 NodeId node_id)
	{
		for (SMHalfedge h : halfedges_around_target(h_sm, sm))
			fill_face_pair_to_node_case_edge_pesf(e_pm, h, fset_sm, pm, sm, node_id);
	}

	void compute_intersection_points_pesf(PMEdgeToSMFaces& pm_edge_to_sm_faces,
										  PolyhedralMesh& pm,
										  SurfaceMesh& sm,
										  VPMPM& vpm_pm,
										  VPMSM& vpm_sm,
										  NodeId& current_node)
	{
		typedef std::tuple<IntersectionType, SMHalfedge, bool, bool> InterInfo;

		for (typename PMEdgeToSMFaces::iterator it = pm_edge_to_sm_faces.begin();
												it != pm_edge_to_sm_faces.end();
												++it)
		{
			PMEdge e_pm = it->first;
			PMHalfedge h_pm = pm.halfedge(e_pm);
			SMFaceSet& fset_sm = it->second;

			while (!fset_sm.empty())
			{
				SMFace f_sm = *fset_sm.begin();

				InterInfo res = compute_intersection_type_pesf(h_pm, f_sm, pm, sm, vpm_pm, vpm_sm);
				IntersectionType type = std::get<0>(res);

				// 由于我们没有实现像halfedges_around_target这样的迭代器, 所以使用已有接口.
				// 而且多面体的incident halfedges of v是没有规律的, 不像表面网格可以转回来, 所以用edge更合适.
				std::vector<PMEdge> all_edges;
				if (std::get<3>(res))
					pm.incident_edges(pm.target(h_pm), std::back_inserter(all_edges));
				else if (std::get<2>(res))
					pm.incident_edges(pm.source(h_pm), std::back_inserter(all_edges));
				else
					all_edges.push_back(e_pm);

				typename std::vector<PMEdge>::iterator it_edge = all_edges.begin();
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
					add_new_node_pesf(h_pm, f_sm, pm, sm, vpm_pm, vpm_sm, res);
					SMHalfedge h_sm = std::get<1>(res);

					visitor.new_node_added_pesf(node_id, ON_FACE, h_pm, h_sm, pm, sm, std::get<2>(res), std::get<3>(res));
					for (; it_edge != all_edges.end(); ++it_edge)
					{
						fill_face_pair_to_node_case_face_pesf(*it_edge, f_sm, pm, sm, node_id);

						if (*it_edge == e_pm)
							fset_sm.erase(f_sm);
						else
						{
							typename PMEdgeToSMFaces::iterator it_pme_smfs =
								pm_edge_to_sm_faces.find(*it_edge);
							if (it_pme_smfs != pm_edge_to_sm_faces.end())
								it_pme_smfs->second.erase(f_sm);
						}
					}
				} // end case ON_FACE
				break;

				case ON_EDGE:
				{
					NodeId node_id = ++current_node;
					add_new_node_pesf(h_pm, f_sm, pm, sm, vpm_pm, vpm_sm, res);
					SMHalfedge h_sm = std::get<1>(res);

					visitor.new_node_added_pesf(node_id, ON_EDGE, h_pm, h_sm, pm, sm, std::get<2>(res), std::get<3>(res));
					for (; it_edge != all_edges.end(); ++it_edge)
					{
						if (*it_edge == e_pm)
							fill_face_pair_to_node_case_edge_pesf(*it_edge, h_sm, &fset_sm, pm, sm, node_id);
						else
						{
							typename PMEdgeToSMFaces::iterator it_pme_smfs =
								pm_edge_to_sm_faces.find(*it_edge);

							SMFaceSet* fset;
							if (it_pme_smfs != pm_edge_to_sm_faces.end())
								fset = &(it_pme_smfs->second);
							else
								fset = nullptr;
							fill_face_pair_to_node_case_edge_pesf(*it_edge, h_sm, fset, pm, sm, node_id);
						}
					}
				} // end case ON_EDGE
				break;

				case ON_VERTEX:
				{
					NodeId node_id = ++current_node;
					SMHalfedge h_sm = std::get<1>(res);
					nodes.add_new_node(get(vpm_sm, target(h_sm, sm)));

					visitor.new_node_added_pesf(node_id, ON_VERTEX, h_pm, h_sm, pm, sm, std::get<2>(res), std::get<3>(res));
					for (; it_edge != all_edges.end(); ++it_edge)
					{
						if (*it_edge == e_pm)
							fill_face_pair_to_node_case_vertex_pesf(*it_edge, h_sm, &fset_sm, pm, sm, node_id);
						else
						{
							typename PMEdgeToSMFaces::iterator it_pme_smfs =
								pm_edge_to_sm_faces.find(*it_edge);

							SMFaceSet* fset;
							if (it_pme_smfs != pm_edge_to_sm_faces.end())
								fset = &(it_pme_smfs->second);
							else
								fset = nullptr;
							fill_face_pair_to_node_case_vertex_pesf(*it_edge, h_sm, fset, pm, sm, node_id);
						}
					}
				} // end case ON_VERTEX
				break;
				} // end switch on the type of the intersection
			} // end loop on all faces that intersect the edge
		} // end loop on all entries (edges) in 'pm_edge_to_sm_faces'

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
	
	// 构建交线基于两个事实: 
	// 1) 不共面的两个三角形对应最多两个交点.
	// 2) surface mesh一条边最多关联两个面, 所以一个交点的邻居最多两个.
	// 
	// 算法可以随意找一个交点i, 找到i的一个邻居j, 从j记录的相邻交点中删掉i, 由于j最多记录两个邻居,
	// 剩下的那个自然是polyline的下一个点. 如此循环下去, 就可以找到polyline的所有交点.
	// 
	// 其中一个网格变为多面体网格后, 事实2不再成立.
	// 当一条pm edge 交在sm face上时, 由于pm edge关联了多个面, 邻居不止两个.
	// 
	// 对于细化来说, 我们需要知道一个交点的所有邻居, 以方便插入约束边/交边. 这个很简单, 全部记录就好.
	// 但对于提取来说, 事实2不再成立造成很多麻烦. 
	// 
	// 一个解决方法是, 以cell为单位construct_polylines, 此时cell可看作一个surface mesh, 1和2都成立.
	//
	void construct_polylines(PolyhedralMesh& pm)
	{
		// 细化的信息.
		std::vector<boost::container::flat_set<NodeId>> node_neighbors(nodes.size());
		// 提取的信息.
		// Note: 此处记录cell的polyline信息不应使用std::map<cell, std::vector<GraphNode>>.
		//       因为只使用nodes中的很少一部分点, 如果node_id偏大, 使用vector会占用大量内存.
		typedef std::map<PMCell, std::map<NodeId, GraphNode>> GraphPerCell;
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
			PMFace f_pm = face_pair.second;
			std::vector<PMCell> cset;
			pm.incident_cells(f_pm, std::back_inserter(cset));
			for (PMCell c : cset)
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
			PMCell c = it->first;
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
	Intersection(SurfaceMesh& sm,
				 PolyhedralMesh& pm,
				 VPMSM& vpm_sm,
				 VPMPM& vpm_pm,
				 const Corefine& v = Corefine())
		: nodes(sm, pm, vpm_sm, vpm_pm)
		, visitor(v)
	{}

	// Intersection functor. 由仿函数组织相交部分的算法.
	void operator()()
	{
		SurfaceMesh& sm = nodes.sm;
		PolyhedralMesh& pm = nodes.pm;
		VPMSM& vpm_sm = nodes.vpm_sm;
		VPMPM& vpm_pm = nodes.vpm_pm;

		// Step1 相交检测.
		filter_intersections_sepf(sm, pm, vpm_sm, vpm_pm);
		filter_intersections_pesf(pm, sm, vpm_pm, vpm_sm);

		// Step2 计算交点并记录相关信息.
		NodeId current_node(std::numeric_limits<NodeId>::max());
		assert(current_node + 1 == 0);

		compute_intersection_points_sepf(sm_edge_to_pm_faces, sm, pm, vpm_sm, vpm_pm, current_node);
		compute_intersection_points_pesf(pm_edge_to_sm_faces, pm, sm, vpm_pm, vpm_sm, current_node);

		// Step3 构建提取所需的交线信息.
		construct_polylines(pm);

		// 转入算法的第二部分, 交给visitor, 由其负责对网格进行细化
		visitor.finalize(nodes, sm, pm, vpm_sm, vpm_pm);
	}

};

}	// namespace MCAL

#endif