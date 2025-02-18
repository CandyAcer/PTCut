// $Intro: ���ļ����и��㷨�ĵ�һ����, ������㽻��, ����ϸ������ȡ�������Ϣ.
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

namespace MCAL    // Mesh Cut Algorithm Libaray, ��Ŷ����������и��㷨�����д���
{

	// �Ķ�������Ҫע���һ���ط�:
	// filter_intersections(), compute_intersection_points()�Ⱥ��������������汾, 
	// ��������polyhedron grid��surface mesh����Ϊ��ͬ���µ�, Ϊ������, �����ں���������˺�׺
	// 
	// sepf: sm edge to pg face
	// pesf: pg edge to sm face
	//


	/*
	* IntersectionImpl���и��㷨�ĵ�һ����, ������������:
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
	* @param IntersectionVisitor: ������������ϸ������, ��Ϊ�����Visitor.
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
			"��������ĵ����Ͳ�һ��");

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

		// ���surface mesh��edge����Щpolyhedron grid��face�ཻ(ָ��Χ���ཻ, ��ʵ���ཻ).
		// ���˵�surface mesh�а�Χ�в������edge.
		void filter_intersections_sepf(SurfaceMesh& sm,
			PolyhedronGrid& pg,
			VPMSM& vpm_sm,
			VPMPG& vpm_pg)
		{
			std::vector<Box> sm_edge_boxes;
			std::vector<Box> pg_face_boxes;

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

			// Ϊpolyhedron grid������face��һ��AABB��Χ��.
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

			// callback�ڼ�⵽��Χ���ཻʱ���м�¼.
			typedef CollectSMEdgeToPGFaces<SurfaceMesh, PolyhedronGrid, SMEdgeToPGFaces>
				Callback;
			Callback callback(sm, pg, sm_edge_to_pg_faces);

			// �����ཻ����㷨.
			CGAL::box_intersection_d(sm_edge_boxes.begin(), sm_edge_boxes.end(),
				pg_face_boxes.begin(), pg_face_boxes.end(),
				callback);

			// ���ڶ�������������ǿռ�����, ����������θ�Ϊ����, ��Ҫ�ڴ˽���⵽����������ǻ�.
			std::unordered_map<PGFace, std::vector<PGFace>> polygon_to_triangles;
			for (typename SMEdgeToPGFaces::iterator it = sm_edge_to_pg_faces.begin();
				it != sm_edge_to_pg_faces.end();
				++it)
			{
				SMEdge e_sm = it->first;
				// �˴�����������, ���ں�������set�����Ԫ��, ���ûᵼ��set�ж�̬����, 
				// ����Ҫ�����ľ���ԭ������Щ��.
				PGFaceSet polygon_set = it->second;

				for (PGFace polygon : polygon_set)
				{
					// �Ը���Ķ�����������������ǻ�, 
					// ���ǻ�������¼�ڴ����������, ����Ԫ��һ����polygon.
					std::vector<PGFace> refined_triangles;
					if (!polygon_to_triangles.count(polygon))
					{
						pg.triangulate_a_polygon(polygon, std::back_inserter(refined_triangles));
						assert(polygon == refined_triangles.front());
						polygon_to_triangles.insert(std::make_pair(polygon, refined_triangles));
					}
					else
						refined_triangles = polygon_to_triangles[polygon];

					// ��¼�����ǻ����������.
					sm_edge_to_pg_faces[e_sm].insert(refined_triangles.begin(), refined_triangles.end());
				}
			}
		}

		// ���polyhedron grid��edge����ЩSurface mesh��face�ཻ(ָ��Χ���ཻ, ��ʵ���ཻ).
		// ���˵�polyhedron grid�а�Χ�в������edge.
		void filter_intersections_pesf(PolyhedronGrid& pg,
			SurfaceMesh& sm,
			VPMPG& vpm_pg,
			VPMSM& vpm_sm)
		{
			std::vector<Box> pg_cell_boxes;
			std::vector<Box> sm_face_boxes;

			// ����Ϊedge�װ�Χ��, ��Ϊ������������ܻ���ڷ���������, �����жϹ���ͻ���鷳.
			// ��Ȼ���, �Ƿ���Ե���box_intersection_dɸ�������ཻ��edge, Ȼ�����ǻ����������?
			// ����. ����ᵼ�¶�һ��filter_pg_intersections(), Ч�ʺܵ�.
			// 
			// Note: 
			// �����˼·���ǻ����²�����edge��Χ�п��ܸ���(�����б�ߵĻ�), �������д���, 
			// ��Ϊ�²�����edgeһ���ڶ���ε��ڲ�, ���������涼����������, ����������ѭ��.

			// Ϊpolyhedron grid������cell����AABB��Χ��.
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
			typedef std::unordered_map<PGCell, std::unordered_set<SMFace>> PGCellToSMFaces;
			PGCellToSMFaces pg_cell_to_sm_faces;

			typedef CollectPGCellToSMFaces<PolyhedronGrid, SurfaceMesh, PGCellToSMFaces>
				Callback;
			Callback callback(pg, sm, pg_cell_to_sm_faces);

			// �����ཻ����㷨.
			CGAL::box_intersection_d(pg_cell_boxes.begin(), pg_cell_boxes.end(),
				sm_face_boxes.begin(), sm_face_boxes.end(),
				callback);

			// ֮����ʹ�ö���������ĵ�Ԫ�����Ǳ߽����ཻ���, ����Ϊ���İ�Χ�бȰ��������б߶���,
			for (typename PGCellToSMFaces::iterator it = pg_cell_to_sm_faces.begin();
				it != pg_cell_to_sm_faces.end();
				++it)
			{
				PGCell c = it->first;
				SMFaceSet& fset_sm = it->second;  // �˴�����������, ��Ϊ����insert, size�̶�����.

				std::vector<PGFace> cell_faces;
				pg.incident_faces(c, std::back_inserter(cell_faces));

				for (PGFace polygon : cell_faces)
					pg.triangulate_a_polygon(polygon, CGAL::Emptyset_iterator());

				// �ռ������polygon���ǻ��������edge.
				std::vector<PGEdge> all_edges;
				pg.incident_edges(c, std::back_inserter(all_edges));
				// �������ǻ�������б�, ���ڱߵİ�Χ��һ���cell��С,
				// ���Բ�һ����fset_sm�е����б������������, ��Ҫ�޳�.
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

						// ��pg edge�İ�Χ�к�sm face�İ�Χ��ȷʵ����, ��¼.
						if (CGAL::do_overlap(edge_bbox, face_bbox))
							pg_edge_to_sm_faces[e_pg].insert(f_sm);
					}
				}
			}
		}

		// �ж�sm face��pg face�Ƿ���.
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

			// �˴�Ӧ�����޳�, �ݲ�����.
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
			// ����sm_edge_to_pg_faces, ��⹲������.
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

			// ����pg_edge_to_sm_faces, ��⹲������.
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

		// ����sm edge��pg face�Ľ�������, �Ծ�ȷ��ʽ��¼��nodes��, ��node idΨһ��ʶ.
		void add_new_node_sepf(SMHalfedge h_sm,
			PGFace f_pg,
			SurfaceMesh& sm,
			PolyhedronGrid& pg,
			VPMSM& vpm_sm,
			VPMPG& vpm_pg,
			std::tuple<IntersectionType, PGHalfedge, bool, bool> res)
		{
			if (std::get<3>(res))  // ������h_sm���յ�
				nodes.add_new_node(get(vpm_sm, target(h_sm, sm)));
			else if (std::get<2>(res))  // ������h_sm�����
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

		// �����������ı���������������Ľ��㲢��¼��ص���Ϣ.
		void compute_intersection_points_sepf(SMEdgeToPGFaces& sm_edge_to_pg_faces,
			SurfaceMesh& sm,
			PolyhedronGrid& pg,
			VPMSM& vpm_sm,
			VPMPG& vpm_pg,
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
			typedef std::tuple<IntersectionType, PGHalfedge, bool, bool>  InterInfo;

			// �����ཻ����������Ԫ��.
			for (typename SMEdgeToPGFaces::iterator it = sm_edge_to_pg_faces.begin();
				it != sm_edge_to_pg_faces.end();
				++it)
			{
				SMEdge e_sm = it->first;
				SMHalfedge h_sm = halfedge(e_sm, sm);
				PGFaceSet& fset_pg = it->second;

				// ����ÿ��e_sm�������.
				while (!fset_pg.empty())
				{
					PGFace f_pg = *fset_pg.begin();
					// Note: e_sm��f_pgͨ��halfedge()����ȡ���İ�߾��в�ȷ����, 
					// ��h_sm�ķ���������, f_pgҲ�������෴������, �⵼�¼����ཻ���͵Ĳ�ȷ����, 
					// ������������������/��Ϊkey��¼��Ϣ��, ֻ�豣֤����res��¼����Ϣ�����к����Ĳ�������.
					InterInfo res = compute_intersection_type_sepf(h_sm, f_pg, sm, pg, vpm_sm, vpm_pg);
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
						fset_pg.erase(f_pg);
						break;

						// Case when the edge pierces the face in its interior.
					case ON_FACE: // ����������ڲ�.
					{
						assert(f_pg == pg.face(std::get<1>(res)));

						NodeId node_id = ++current_node;  // ����node id
						// ��һ����Ϣ: �������������, �����nodes��.
						add_new_node_sepf(h_sm, f_pg, sm, pg, vpm_sm, vpm_pg, res);
						PGHalfedge h_pg = std::get<1>(res);

						// �ڶ�����Ϣ: ϸ���������Ϣ, ��visitor����.
						visitor.new_node_added_sepf(node_id, ON_FACE, h_sm, h_pg, sm, pg, std::get<2>(res), std::get<3>(res));

						// ��������Ϣ: polyline�����Ϣ, ���ڹ���output_bulider��Ҫ������.
						for (; it_halfedge != all_halfedges.end(); ++it_halfedge)
						{
							fill_face_pair_to_node_case_face_sepf(*it_halfedge, f_pg, sm, pg, node_id);

							if (it_halfedge == all_halfedges.begin())
								fset_pg.erase(f_pg);
							else
							{
								// �ӹ����ı߼�¼��FaceSet��, ɾ�������, �����μ���.
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
						// ��һ����Ϣ: �������������, �����nodes��.
						add_new_node_sepf(h_sm, f_pg, sm, pg, vpm_sm, vpm_pg, res);
						PGHalfedge h_pg = std::get<1>(res);

						// �ڶ�����Ϣ: ϸ���������Ϣ, ��visitor����.
						visitor.new_node_added_sepf(node_id, ON_EDGE, h_sm, h_pg, sm, pg, std::get<2>(res), std::get<3>(res));

						// ��������Ϣ: polyline�����Ϣ, ���ڹ���output_bulider��Ҫ������.
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

						// ��һ����Ϣ: �������������, �����nodes��.
						nodes.add_new_node(vpm_pg[pg.target(h_pg)]);

						// �ڶ�����Ϣ: ϸ���������Ϣ, ��visitor����.
						visitor.new_node_added_sepf(node_id, ON_VERTEX, h_sm, h_pg, sm, pg, std::get<2>(res), std::get<3>(res));

						// ��������Ϣ: polyline�����Ϣ, ���ڹ���output_bulider��Ҫ������.
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

					// ��������û��ʵ����halfedges_around_target�����ĵ�����, ����ʹ�����нӿ�.
					// ���Ҷ������incident halfedges of v��û�й��ɵ�, ��������������ת����, ������edge������.
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

		// CGAL�������߻���������ʵ: 
		// 1) ����������������ζ�Ӧ�����������.
		// 2) surface meshһ����������������, ����һ��������ھ��������.
		// 
		// �㷨����������һ������i, �ҵ�i��һ���ھ�j, ��j��¼�����ڽ�����ɾ��i, ����j����¼�����ھ�,
		// ʣ�µ��Ǹ���Ȼ��polyline����һ����. ���ѭ����ȥ, �Ϳ����ҵ�polyline�����н���.
		// 
		// ����һ�������Ϊ�����������, ��ʵ2���ٳ���.
		// ��һ��pg edge ����sm face��ʱ, ����pg edge�����˶����, �ھӲ�ֹ����.
		// 
		// ����ϸ����˵, ������Ҫ֪��һ������������ھ�, �Է������Լ����/����. ����ܼ�, ȫ����¼�ͺ�.
		// ��������ȡ��˵, ��ʵ2���ٳ�����ɺܶ��鷳. 
		// 
		// һ�����������, ��cellΪ��λconstruct_polylines, ��ʱcell�ɿ���һ��surface mesh, 1��2������.
		//
		void construct_polylines(PolyhedronGrid& pg)
		{
			// ϸ������Ϣ.
			std::vector<boost::container::flat_set<NodeId>> node_neighbors(nodes.size());
			// ��ȡ����Ϣ.
			// Note: �˴���¼cell��polyline��Ϣ��Ӧʹ��std::map<cell, std::vector<GraphNode>>.
			//       ��Ϊֻʹ��nodes�еĺ���һ���ֵ�, ���node_idƫ��, ʹ��vector��ռ�ô����ڴ�.
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

				// ��¼ÿ��cell�Ͻ����������Ϣ.
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

			// ��cellΪ��λ, ����ÿ��cell��polyline��Ϣ
			for (typename GraphPerCell::iterator it = cell_graph.begin();
				it != cell_graph.end();
				++it)
			{
				PGCell c = it->first;
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
		IntersectionImpl(SurfaceMesh& sm,
			PolyhedronGrid& pg,
			VPMSM& vpm_sm,
			VPMPG& vpm_pg,
			const IntersectionVisitor& v = IntersectionVisitor())
			: nodes(sm, pg, vpm_sm, vpm_pg)
			, visitor(v)
		{
		}

		// Intersection functor. �ɷº�����֯�ཻ���ֵ��㷨.
		void operator()()
		{
			SurfaceMesh& sm = nodes.sm;
			PolyhedronGrid& pg = nodes.pg;
			VPMSM& vpm_sm = nodes.vpm_sm;
			VPMPG& vpm_pg = nodes.vpm_pg;

			// Step1 �ཻ���.
			filter_intersections_sepf(sm, pg, vpm_sm, vpm_pg);
			filter_intersections_pesf(pg, sm, vpm_pg, vpm_sm);
			//filter_coplanar_case(sm, pg, vpm_sm, vpm_pg);

			std::cout << "�ཻ������-----------------------------" << std::endl;
			std::cout << "��ȷ����Ĵ�����" << iExactComputeCount << std::endl;
			std::cout << "��������Ĵ�����" << iIntervalComputeCount << std::endl;
			std::cout << "----------------------------------------\n\n";


			std::cout << "ɸѡ����sm_edge������" << sm_edge_to_pg_faces.size() << std::endl;
			std::unordered_set<PGFace> pg_face_set;
			std::size_t cnt = 0;
			for (auto it = sm_edge_to_pg_faces.begin(); it != sm_edge_to_pg_faces.end(); ++it)
			{
				pg_face_set.insert(it->second.begin(), it->second.end());
				cnt += it->second.size();
			}
			std::cout << "ɸѡ����pg_face������û��ȥ�أ���" << cnt << std::endl;
			std::cout << "ɸѡ����pg_face������ȥ�غ󣩣�" << pg_face_set.size() << "\n\n";


			std::cout << "ɸѡ����pg_edge������" << pg_edge_to_sm_faces.size() << std::endl;
			std::unordered_set<SMFace> sm_face_set;
			cnt = 0;
			for (auto it = pg_edge_to_sm_faces.begin(); it != pg_edge_to_sm_faces.end(); ++it)
			{
				sm_face_set.insert(it->second.begin(), it->second.end());
				cnt += it->second.size();
			}
			std::cout << "ɸѡ����sm_face������û��ȥ�أ���" << cnt << std::endl;
			std::cout << "ɸѡ����sm_face������ȥ�أ���" << sm_face_set.size() << "\n\n";

			// Step2 ���㽻�㲢��¼�����Ϣ.
			NodeId current_node(std::numeric_limits<NodeId>::max());
			assert(current_node + 1 == 0);

			int cur_exact_cnt = iExactComputeCount;
			int cur_interval_cnt = iIntervalComputeCount;

			compute_intersection_points_sepf(sm_edge_to_pg_faces, sm, pg, vpm_sm, vpm_pg, current_node);
			compute_intersection_points_pesf(pg_edge_to_sm_faces, pg, sm, vpm_pg, vpm_sm, current_node);

			std::cout << "���㽻�����-----------------------------" << std::endl;
			std::cout << "��ȷ����Ĵ�����" << iExactComputeCount - cur_exact_cnt << std::endl;
			std::cout << "��������Ĵ�����" << iIntervalComputeCount - cur_interval_cnt << std::endl;
			std::cout << "----------------------------------------\n\n";

			// Step3 ������ȡ����Ľ�����Ϣ.
			construct_polylines(pg);

			// ת���㷨�ĵڶ�����, ����visitor, ���为����������ϸ��
			visitor.finalize(nodes, sm, pg, vpm_sm, vpm_pg);
		}

	};

}	// namespace MCAL


#endif