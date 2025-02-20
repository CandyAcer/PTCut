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


namespace MCAL    // Mesh Cut Algorithm Libaray, ��Ŷ����������и��㷨�����д���
{

// �Ķ�������Ҫע���һ���ط�:
// filter_intersections(), compute_intersection_points()�Ⱥ��������������汾, 
// ��������polyhedral mesh��surface mesh����Ϊ��ͬ���µ�, 
// Ϊ������, �����ں���������˺�׺
// 
// sepf: sm edge to pm face
// pesf: pm edge to sm face
//


/*
* Intersection���и��㷨�ĵ�һ����, ������������:
*   1) ����������ı�����Щ�������������ཻ, �Կ����ཻ�Ķ��������������ǻ�,
*      ���㼸�μ���; ������������ı�����Щ��������������ཻ.
*   2) ����Щɸѡ������Ϣ, ������������εĽ���, ͬʱ��¼�����Ϣ, ��Ҫ�����¼��ࣺ
*        2.1) �����������Ϣ, ����IntersectionNodesͳһ����, ��NodeId���б�ʶ.
*        2.2) �������ĸ���, �����߻��ĸ�����, ��¼�����Ӧ��ϵ.
*        2.3) ������Ϣ, ���������ν����ϵ�����������.
*   3) ��cellΪ��λ��������.
*
* 2��¼����Ϣ��ϸ������, �ʼ�¼��ϸ�����data member��.
* 3��¼����Ϣ���������������ս��, �ʼ�¼��OutputBuilder��data member��.
*
* @param Corefine: ������������ϸ������, ��Ϊ�����Visitor.
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
				   "��������ĵ����Ͳ�һ��");

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

	// ���surface mesh��edge����Щpolyhedral mesh��face�ཻ(ָ��Χ���ཻ, ��ʵ���ཻ).
	// ���˵�surface mesh�а�Χ�в������edge.
	void filter_intersections_sepf(SurfaceMesh& sm,
								   PolyhedralMesh& pm,
								   VPMSM& vpm_sm,
								   VPMPM& vpm_pm)
	{
		std::vector<Box> sm_edge_boxes;
		std::vector<Box> pm_face_boxes;

		// ΪSurface mesh���е�edge��һ��AABB��Χ��.
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

		// Ϊpolyhedral mesh������face��һ��AABB��Χ��.
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

		// callback�ڼ�⵽��Χ���ཻʱ���м�¼.
		typedef CollectSMEdgeToPMFaces<SurfaceMesh, PolyhedralMesh, SMEdgeToPMFaces>
			Callback;
		Callback callback(sm, pm, sm_edge_to_pm_faces);

		// �����ཻ����㷨.
		CGAL::box_intersection_d(sm_edge_boxes.begin(), sm_edge_boxes.end(),
								 pm_face_boxes.begin(), pm_face_boxes.end(),
								 callback);

		// ���ڶ�������������ǿռ�����, ����������θ�Ϊ����, ��Ҫ�ڴ˽���⵽����������ǻ�.
		std::unordered_map<PMFace, std::vector<PMFace>> polygon_to_triangles;
		for (typename SMEdgeToPMFaces::iterator it = sm_edge_to_pm_faces.begin();
												it != sm_edge_to_pm_faces.end();
												++it)
		{
			SMEdge e_sm = it->first;
			// �˴�����������, ���ں�������set�����Ԫ��, ���ûᵼ��set�ж�̬����, 
			// ����Ҫ�����ľ���ԭ������Щ��.
			PMFaceSet polygon_set = it->second;

			for (PMFace polygon : polygon_set)
			{
				// �Ը���Ķ�����������������ǻ�, 
				// ���ǻ�������¼�ڴ����������, ����Ԫ��һ����polygon.
				std::vector<PMFace> refined_triangles;
				if (!polygon_to_triangles.count(polygon))
				{
					pm.triangulate_a_polygon(polygon, std::back_inserter(refined_triangles));
					assert(polygon == refined_triangles.front());
					polygon_to_triangles.insert(std::make_pair(polygon, refined_triangles));
				}
				else
					refined_triangles = polygon_to_triangles[polygon];

				// ��¼�����ǻ����������.
				sm_edge_to_pm_faces[e_sm].insert(refined_triangles.begin(), refined_triangles.end());
			}
		}
	}

	// ���polyhedral mesh��edge����ЩSurface mesh��face�ཻ(ָ��Χ���ཻ, ��ʵ���ཻ).
	// ���˵�polyhedral mesh�а�Χ�в������edge.
	void filter_intersections_pesf(PolyhedralMesh& pm,
								   SurfaceMesh& sm,
								   VPMPM& vpm_pm,
								   VPMSM& vpm_sm)
	{
		std::vector<Box> pm_cell_boxes;
		std::vector<Box> sm_face_boxes;

		// ����Ϊedge�װ�Χ��, ��Ϊ������������ܻ���ڷ���������, �����жϹ���ͻ���鷳.
		// ��Ȼ���, �Ƿ���Ե���box_intersection_dɸ�������ཻ��edge, Ȼ�����ǻ����������?
		// ����. ����ᵼ�¶�һ��filter_pm_intersections(), Ч�ʺܵ�.
		// 
		// Note: 
		// �����˼·���ǻ����²�����edge��Χ�п��ܸ���(�����б�ߵĻ�), �������д���, 
		// ��Ϊ�²�����edgeһ���ڶ���ε��ڲ�, ���������涼����������, ����������ѭ��.

		// Ϊpolyhedral mesh������cell����AABB��Χ��.
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

		// Ϊsurface mesh���е�face����AABB��Χ��.
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

		// callback�����⵽��Χ���ཻʱ���м�¼.
		typedef std::unordered_map<PMCell, std::unordered_set<SMFace>> PMCellToSMFaces;
		PMCellToSMFaces pm_cell_to_sm_faces;

		typedef CollectPMCellToSMFaces<PolyhedralMesh, SurfaceMesh, PMCellToSMFaces>
			Callback;
		Callback callback(pm, sm, pm_cell_to_sm_faces);

		// �����ཻ����㷨.
		CGAL::box_intersection_d(pm_cell_boxes.begin(), pm_cell_boxes.end(),
								 sm_face_boxes.begin(), sm_face_boxes.end(),
								 callback);

		// ֮����ʹ�ö���������ĵ�Ԫ�����Ǳ߽����ཻ���, ����Ϊ���İ�Χ�бȰ��������б߶���,
		for (typename PMCellToSMFaces::iterator it = pm_cell_to_sm_faces.begin();
												it != pm_cell_to_sm_faces.end();
												++it)
		{
			PMCell c = it->first;
			SMFaceSet& fset_sm = it->second;  // �˴�����������, ��Ϊ����insert, size�̶�����.

			std::vector<PMFace> cell_faces;
			pm.incident_faces(c, std::back_inserter(cell_faces));

			for (PMFace polygon : cell_faces)
				pm.triangulate_a_polygon(polygon, CGAL::Emptyset_iterator());

			// �ռ������polygon���ǻ��������edge.
			std::vector<PMEdge> all_edges;
			pm.incident_edges(c, std::back_inserter(all_edges));
			// �������ǻ�������б�, ���ڱߵİ�Χ��һ���cell��С,
			// ���Բ�һ����fset_sm�е����б������������, ��Ҫ�޳�.
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

					// ��pm edge�İ�Χ�к�sm face�İ�Χ��ȷʵ����, ��¼.
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
			// ���ձ�̹淶��˵, using����ǽ��õ�, ���ǽ���switch����ʹ��, ��֤Ӱ�췶Χ������С.
			// Ŀ����Ϊ��ʹ�ô��뾡���ܼ��.
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

	// ����sm edge��pm face�Ľ�������, �Ծ�ȷ��ʽ��¼��nodes��, ��node idΨһ��ʶ.
	void add_new_node_sepf(SMHalfedge h_sm,
						   PMFace f_pm,
						   SurfaceMesh& sm,
						   PolyhedralMesh& pm,
						   VPMSM& vpm_sm,
						   VPMPM& vpm_pm,
						   std::tuple<IntersectionType, PMHalfedge, bool, bool> res)
	{
		if (std::get<3>(res))  // ������h_sm���յ�
			nodes.add_new_node(get(vpm_sm, target(h_sm, sm)));
		else if (std::get<2>(res))  // ������h_sm�����
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

	// �����������ı���������������Ľ��㲢��¼��ص���Ϣ.
	void compute_intersection_points_sepf(SMEdgeToPMFaces& sm_edge_to_pm_faces,
										  SurfaceMesh& sm,
										  PolyhedralMesh& pm,
										  VPMSM& vpm_sm,
										  VPMPM& vpm_pm,
										  NodeId& current_node)
	{
		/*
		* �ڴ˽���һ��std::tuple<IntersectionType, Halfedge, bool, bool>�ĺ���.
		* tuple[0]: ON_FACE/ON_EDGE/ON_VERTEX, ˵���ཻ������.
		* tuple[1]: ���������ϵ�һ��halfedge. ע�������߶α��޹�.
		*			���ཻ����ΪON_FACE, halfedge��face�������м�¼�İ��;
		*			���ཻ����ΪON_EDGE, halfedge���ཻ�߶�Ӧ�İ��;
		*			���ཻ����ΪON_VERTEX, halfedgeָ���Ǹ�vertex;
		* tuple[2]: �߶αߵ���ʼ���Ƿ�����������.
		* tuple[3]: �߶αߵ��յ��Ƿ�����������.
		*
		* ֮�������������, ����Ϊ����Ҫ�ռ�ÿ��face/edge�ϵĽ���, ������tuple�з���1~3����Ϣ, ����ͨ�����˲�ѯ��¼�����Ϣ.
		*/
		typedef std::tuple<IntersectionType, PMHalfedge, bool, bool>  InterInfo;

		// �����ཻ����������Ԫ��.
		for (typename SMEdgeToPMFaces::iterator it = sm_edge_to_pm_faces.begin();
												it != sm_edge_to_pm_faces.end();
												++it)
		{
			SMEdge e_sm = it->first;
			SMHalfedge h_sm = halfedge(e_sm, sm);
			PMFaceSet& fset_pm = it->second;

			// ����ÿ��e_sm�������.
			while (!fset_pm.empty())
			{
				PMFace f_pm = *fset_pm.begin();
				// Note: e_sm��f_pmͨ��halfedge()����ȡ���İ�߾��в�ȷ����, 
				// ��h_sm�ķ���������, f_pmҲ�������෴������, �⵼�¼����ཻ���͵Ĳ�ȷ����, 
				// ������������������/��Ϊkey��¼��Ϣ��, ֻ�豣֤����res��¼����Ϣ�����к����Ĳ�������.
				InterInfo res = compute_intersection_type_sepf(h_sm, f_pm, sm, pm, vpm_sm, vpm_pm);
				IntersectionType type = std::get<0>(res);

				// ע��, �������ڱ���edge_to_faces�Ĺ�����, �ɱ��������μ���ó���,
				// ���Ǹ�ÿ���������һ��NodeId, Ψһ�ر�ʶ��, ���������������. 
				// ��������Ǳߵ�һ���˵�, ������������ܶ�����, �����������εĽ��㶼���������.
				// �����ͻ��ظ�����, ��������һ������, ȴ�����˶��NodeId, �������ԸΥ��.
				// ������Ҫ�ӹ����ı߼�¼��FaceSet��, ɾ�������, �����μ���.
				// 
				// all_halfedges���������ռ���Щ�����ıߵ�.
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
					// ���ǲ������������, ֮ǰ�Ĳ����Ѿ����˵���.
					/*std::cout << "Coplanar triangles: this point should never be reached!";*/
					break;

				case EMPTY:  // ���ཻ, ɾ�����Ŵ����¸��漴��.
					fset_pm.erase(f_pm);
					break;

					// Case when the edge pierces the face in its interior.
				case ON_FACE: // ����������ڲ�.
				{
					assert(f_pm == pm.face(std::get<1>(res)));

					NodeId node_id = ++current_node;  // ����node id
					// ��һ����Ϣ: �������������, �����nodes��.
					add_new_node_sepf(h_sm, f_pm, sm, pm, vpm_sm, vpm_pm, res);
					PMHalfedge h_pm = std::get<1>(res);

					// �ڶ�����Ϣ: ϸ���������Ϣ, ��visitor����.
					visitor.new_node_added_sepf(node_id, ON_FACE, h_sm, h_pm, sm, pm, std::get<2>(res), std::get<3>(res));

					// ��������Ϣ: polyline�����Ϣ, ���ڹ���output_bulider��Ҫ������.
					for (; it_halfedge != all_halfedges.end(); ++it_halfedge)
					{
						fill_face_pair_to_node_case_face_sepf(*it_halfedge, f_pm, sm, pm, node_id);

						if (it_halfedge == all_halfedges.begin())
							fset_pm.erase(f_pm);
						else
						{
							// �ӹ����ı߼�¼��FaceSet��, ɾ�������, �����μ���.
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
					// ��һ����Ϣ: �������������, �����nodes��.
					add_new_node_sepf(h_sm, f_pm, sm, pm, vpm_sm, vpm_pm, res);
					PMHalfedge h_pm = std::get<1>(res);

					// �ڶ�����Ϣ: ϸ���������Ϣ, ��visitor����.
					visitor.new_node_added_sepf(node_id, ON_EDGE, h_sm, h_pm, sm, pm, std::get<2>(res), std::get<3>(res));

					// ��������Ϣ: polyline�����Ϣ, ���ڹ���output_bulider��Ҫ������.
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

					// ��һ����Ϣ: �������������, �����nodes��.
					nodes.add_new_node(vpm_pm[pm.target(h_pm)]);

					// �ڶ�����Ϣ: ϸ���������Ϣ, ��visitor����.
					visitor.new_node_added_sepf(node_id, ON_VERTEX, h_sm, h_pm, sm, pm, std::get<2>(res), std::get<3>(res));

					// ��������Ϣ: polyline�����Ϣ, ���ڹ���output_bulider��Ҫ������.
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

	//================ ���´�����������첻��, �ɲο������ע�� =======================

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

				// ��������û��ʵ����halfedges_around_target�����ĵ�����, ����ʹ�����нӿ�.
				// ���Ҷ������incident halfedges of v��û�й��ɵ�, ��������������ת����, ������edge������.
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

		// degree����һ.
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
	
	// �������߻���������ʵ: 
	// 1) ����������������ζ�Ӧ�����������.
	// 2) surface meshһ����������������, ����һ��������ھ��������.
	// 
	// �㷨����������һ������i, �ҵ�i��һ���ھ�j, ��j��¼�����ڽ�����ɾ��i, ����j����¼�����ھ�,
	// ʣ�µ��Ǹ���Ȼ��polyline����һ����. ���ѭ����ȥ, �Ϳ����ҵ�polyline�����н���.
	// 
	// ����һ�������Ϊ�����������, ��ʵ2���ٳ���.
	// ��һ��pm edge ����sm face��ʱ, ����pm edge�����˶����, �ھӲ�ֹ����.
	// 
	// ����ϸ����˵, ������Ҫ֪��һ������������ھ�, �Է������Լ����/����. ����ܼ�, ȫ����¼�ͺ�.
	// ��������ȡ��˵, ��ʵ2���ٳ�����ɺܶ��鷳. 
	// 
	// һ�����������, ��cellΪ��λconstruct_polylines, ��ʱcell�ɿ���һ��surface mesh, 1��2������.
	//
	void construct_polylines(PolyhedralMesh& pm)
	{
		// ϸ������Ϣ.
		std::vector<boost::container::flat_set<NodeId>> node_neighbors(nodes.size());
		// ��ȡ����Ϣ.
		// Note: �˴���¼cell��polyline��Ϣ��Ӧʹ��std::map<cell, std::vector<GraphNode>>.
		//       ��Ϊֻʹ��nodes�еĺ���һ���ֵ�, ���node_idƫ��, ʹ��vector��ռ�ô����ڴ�.
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

			// ��¼ÿ��cell�Ͻ����������Ϣ.
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

		// ��cellΪ��λ, ����ÿ��cell��polyline��Ϣ
		for (typename GraphPerCell::iterator it = cell_graph.begin();
			it != cell_graph.end();
			++it)
		{
			PMCell c = it->first;
			std::map<NodeId, GraphNode>& graph_nodes = it->second;

			std::size_t nb_nodes = graph_nodes.size();

			// dynamic_bitsetҪ��������������, ��graph_nodes�е�node id��������.
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

			// �ȴ����պϵĽ���.
			while (terminal_nodes.any())
			{
				// ����ҵ����պϽ��ߵ�һ���˵�, ������ʼ����һ������.
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

			// ����պϵĽ���.
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

	// Intersection functor. �ɷº�����֯�ཻ���ֵ��㷨.
	void operator()()
	{
		SurfaceMesh& sm = nodes.sm;
		PolyhedralMesh& pm = nodes.pm;
		VPMSM& vpm_sm = nodes.vpm_sm;
		VPMPM& vpm_pm = nodes.vpm_pm;

		// Step1 �ཻ���.
		filter_intersections_sepf(sm, pm, vpm_sm, vpm_pm);
		filter_intersections_pesf(pm, sm, vpm_pm, vpm_sm);

		// Step2 ���㽻�㲢��¼�����Ϣ.
		NodeId current_node(std::numeric_limits<NodeId>::max());
		assert(current_node + 1 == 0);

		compute_intersection_points_sepf(sm_edge_to_pm_faces, sm, pm, vpm_sm, vpm_pm, current_node);
		compute_intersection_points_pesf(pm_edge_to_sm_faces, pm, sm, vpm_pm, vpm_sm, current_node);

		// Step3 ������ȡ����Ľ�����Ϣ.
		construct_polylines(pm);

		// ת���㷨�ĵڶ�����, ����visitor, ���为����������ϸ��
		visitor.finalize(nodes, sm, pm, vpm_sm, vpm_pm);
	}

};

}	// namespace MCAL

#endif