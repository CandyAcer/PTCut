// $Intro: 本文件是切割算法的第三部分, 构建最终的多面体网格.
// $License: GPL-3.0-or-later
// 
// 
// $Author: Liu Zhiqiang




#ifndef MCAL_ALGORITHM_OUTPUT_BUILDER_H
#define MCAL_ALGORITHM_OUTPUT_BUILDER_H

#include <CGAL/utility.h>

#include <boost/dynamic_bitset/dynamic_bitset.hpp>

#include <queue>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <bitset>

#include "utils/sorted_pair.h"
#include "predicates.h"
#include <ctime>

extern std::clock_t start;


namespace MCAL   // Mesh Cut Algorithm Libaray, 存放多面体网格切割算法的实现
{

	enum BO_Type { UNION = 0, INTERSECTION, SM_MINUS_PG, PG_MINUS_SM, NONE };

	/*
	* OutputBuilder类负责构建最终的多面体网格.
	*
	*/
	template <typename SurfaceMesh,
		typename PolyhedronGrid,
		typename VPMSM,
		typename VPMPG>
	class OutputBuilder
	{
		// typedefs
		typedef boost::graph_traits<SurfaceMesh>              GT;
		typedef typename GT::vertex_descriptor                SMVertex;
		typedef typename GT::edge_descriptor                  SMEdge;
		typedef typename GT::halfedge_descriptor              SMHalfedge;
		typedef typename GT::face_descriptor                  SMFace;

		typedef typename PolyhedronGrid::VertexIndex          PGVertex;
		typedef typename PolyhedronGrid::EdgeIndex            PGEdge;
		typedef typename PolyhedronGrid::FaceIndex            PGFace;
		typedef typename PolyhedronGrid::CellIndex            PGCell;
		typedef typename PolyhedronGrid::Halfedge             PGHalfedge;

		//============================ 交线信息的记录 ===================================
			// 在此解释一下polyline是如何记录的.
			// 首先, 为了尽可能复用CGAL的构建方法, 我们记录每个cell与surface mesh的polyline.
			// 其次, 对于每条polyline, 虽然它包含很多交边, 但我们只记录一条.
			// 
			// Question: 只记录一条, 那其他的交边从何得知呢?
			// Answer: 我们还会记录每个cell的IntersectSMEdge和IntersectPGEdge, polyline中的
			//         所有交边都在里面, 通过拓扑查询总可以搜索到.
			//
		typedef std::size_t                                   NodeId;
		typedef std::pair<NodeId, NodeId>                     NodeIdPair;
		// 一条polyline记录的信息.
		struct PolylineInfo
		{
			SMHalfedge sm_hedge;    // polyline某条交边对应的sm halfedge
			PGHalfedge pg_hedge;    // polyline某条交边对应的pg halfedge
			bool is_reverse;        // The key (pair<NodeId, NodeId>) was reversed?
			std::size_t node_num;   // polyline包含的交点数目

			PolylineInfo(SMHalfedge h_sm = SMHalfedge(),
				PGHalfedge h_pg = PGHalfedge(),
				bool b = false,
				std::size_t n = 0)
				:sm_hedge(h_sm), pg_hedge(h_pg), is_reverse(b), node_num(n)
			{
			}
		};
		// AnEdgePerPolyline + IntersectSMEdge + IntersectPGEdge = 所有的相交信息.
		// 我们为每个cell都记录这三类信息.
		typedef std::map<NodeIdPair, PolylineInfo>            AnEdgePerPolyline;
		typedef std::map<PGCell, AnEdgePerPolyline>           CellPolylinesMap;

		typedef std::unordered_set<SMEdge>                    IntersectSMEdge;
		typedef std::unordered_map<PGCell, IntersectSMEdge>   CellToIntersectSMEdge;
		typedef std::unordered_set<PGEdge>                    IntersectPGEdge;
		typedef std::unordered_map<PGCell, IntersectPGEdge>   CellToIntersectPGEdge;
		//==============================================================================

		typedef std::unordered_map<SMVertex, NodeId>          SMVertexToNode;
		typedef std::unordered_map<PGVertex, NodeId>          PGVertexToNode;

		typedef std::map<SortedPair<NodeId>, boost::container::flat_set<PGCell>>  IdCellMap;

		// data members
	private:
		SurfaceMesh& sm;
		PolyhedronGrid& pg;
		VPMSM& vpm_sm;
		VPMPG& vpm_pg;

		CellPolylinesMap cell_polylines;
		typename AnEdgePerPolyline::iterator last_polyline;
		CellToIntersectSMEdge cell_intersect_sm_edges;   // 每个cell与sm的交线的所有sm edge
		CellToIntersectPGEdge cell_intersect_pg_edges;   // 每个cell与sm的交线的所有pg edge	
		SMVertexToNode sm_vertex_to_node_id;
		PGVertexToNode pg_vertex_to_node_id;

		IdCellMap id_cell;

		const NodeId NID;

	public:
		OutputBuilder(SurfaceMesh& sm,
			PolyhedronGrid& pg,
			VPMSM& vpm_sm,
			VPMPG& vpm_pg)
			: sm(sm), pg(pg), vpm_sm(vpm_sm), vpm_pg(vpm_pg)
			, NID((std::numeric_limits<NodeId>::max)())
		{
		}

		// 设置sm vertex与node id的对应关系.
		void set_sm_vertex_id(SMVertex v, NodeId node_id)
		{
			sm_vertex_to_node_id.insert(std::make_pair(v, node_id));
		}

		// 设置pg vertex与node id的对应关系.
		void set_pg_vertex_id(PGVertex v, NodeId node_id)
		{
			pg_vertex_to_node_id.insert(std::make_pair(v, node_id));
		}

		// 记录交边与cell的对应关系.
		void record_segment_cell(SortedPair<NodeId>& segment, PGCell c)
		{
			id_cell[segment].insert(c);
		}

		void start_new_polyline(PGCell c, NodeId i, NodeId j)
		{
			std::pair<typename AnEdgePerPolyline::iterator, bool> res =
				cell_polylines[c].insert(std::make_pair(CGAL::make_sorted_pair(i, j), PolylineInfo()));

			assert(res.second);
			last_polyline = res.first;
			if (i != last_polyline->first.first)
			{
				last_polyline->second.is_reverse = true;
			}
		}

		void add_node_to_polyline(PGCell c, NodeId i)
		{
			++(last_polyline->second.node_num);
		}

		void set_sm_edge_per_polyline(SurfaceMesh& sm,
			NodeIdPair id_pair,
			SMHalfedge h_sm)
		{
			SortedPair<NodeId> segment(id_pair.first, id_pair.second);
			typename IdCellMap::iterator it = id_cell.find(segment);
			assert(it != id_cell.end());

			boost::container::flat_set<PGCell>& cset = it->second;
			assert(cset.size() <= 2);

			for (PGCell c : cset)
			{
				cell_intersect_sm_edges[c].insert(edge(h_sm, sm));

				if (id_pair.first > id_pair.second)
				{
					std::swap(id_pair.first, id_pair.second);
					h_sm = opposite(h_sm, sm);
				}

				assert(cell_polylines.count(c));
				AnEdgePerPolyline& an_edge_per_polyline = cell_polylines[c];
				typename AnEdgePerPolyline::iterator it_poly = an_edge_per_polyline.find(id_pair);

				if (it_poly != an_edge_per_polyline.end())
				{
					it_poly->second.sm_hedge = h_sm;
				}
			}
		}

		// 对应关系: v_pg对应id_pair.first, e_pg是id_pair所在的物理边.
		void set_pg_edge_per_polyline(NodeIdPair id_pair,
			PGVertex srcv,
			PGEdge target_edge,
			PolyhedronGrid& pg)
		{
			assert(id_pair.first < id_pair.second);
			assert(pg_vertex_to_node_id[srcv] == id_pair.first);

			SortedPair<NodeId> segment(id_pair.first, id_pair.second);
			typename IdCellMap::iterator it = id_cell.find(segment);
			assert(it != id_cell.end());

			boost::container::flat_set<PGCell>& cset = it->second;
			assert(cset.size() <= 2);

			for (PGCell c : cset)
			{
				cell_intersect_pg_edges[c].insert(target_edge);

				// 找到id_pair对应的半边, 这个半边满足:
				//   1) h_pg是c当中的半边
				//   2) h_pg的source是v_pg, 对应着id_pair.first
				PGHalfedge h_pg = pg.halfedge(target_edge, srcv, c);
				assert(h_pg != pg.null_halfedge());

				AnEdgePerPolyline& an_edge_per_polyline = cell_polylines[c];
				typename AnEdgePerPolyline::iterator it_poly = an_edge_per_polyline.find(id_pair);

				if (it_poly != an_edge_per_polyline.end())
					it_poly->second.pg_hedge = h_pg;
			}
		}

		template <typename Vertex, typename VertexToNodeId>
		NodeId get_node_id(Vertex v, const VertexToNodeId& node_ids)
		{
			typename VertexToNodeId::const_iterator it = node_ids.find(v);
			if (it == node_ids.end())
				return NID;
			return it->second;
		}

		int mark_sm_faces(IntersectSMEdge& intersect_sm_edge,
			std::vector<std::size_t>& sm_patch_ids,
			SurfaceMesh& sm)
		{
			int current_patch_id = 0;

			std::vector<bool> marked(num_faces(sm), false);
			for (SMFace f : faces(sm))
			{
				if (marked[f])
					continue;

				std::queue<SMFace> Queue;
				Queue.push(f);
				while (!Queue.empty())
				{
					SMFace current_face = Queue.front();
					Queue.pop();

					if (marked[current_face])
						continue;

					marked[current_face] = true;
					sm_patch_ids[current_face] = current_patch_id;

					for (SMHalfedge h : halfedges_around_face(halfedge(current_face, sm), sm))
					{
						if (is_border(edge(h, sm), sm) || intersect_sm_edge.count(edge(h, sm)))
							continue;

						SMFace f_opp = face(opposite(h, sm), sm);
						if (f_opp != GT::null_face() && !marked[f_opp])
							Queue.push(f_opp);
					}
				}
				++current_patch_id;
			}
			return current_patch_id;
		}

		int mark_pg_faces(PGCell c,
			IntersectPGEdge& intersect_pg_edge,
			std::map<PGFace, std::size_t>& pg_patch_ids,
			PolyhedronGrid& pg)
		{
			int current_patch_id = 0;

			std::vector<PGFace> c_faces;
			pg.incident_faces(c, std::back_inserter(c_faces));

			std::unordered_map<PGFace, bool> marked;

			for (std::size_t i = 0; i < c_faces.size(); ++i)
			{
				marked.insert(std::make_pair(c_faces[i], false));
				pg_patch_ids.insert(std::make_pair(c_faces[i], NID));
			}

			for (PGFace f : c_faces)
			{
				if (marked[f])
					continue;

				std::queue<PGFace> Queue;
				Queue.push(f);
				while (!Queue.empty())
				{
					PGFace current_face = Queue.front();
					Queue.pop();

					if (marked[current_face])
						continue;

					marked[current_face] = true;
					pg_patch_ids[current_face] = current_patch_id;

					std::vector<PGHalfedge> f_hedges;
					pg.halfedges_around_face(current_face, c, std::back_inserter(f_hedges));

					for (PGHalfedge h : f_hedges)
					{
						if (intersect_pg_edge.count(pg.edge(h)))
							continue;

						PGFace f_opp = pg.face(pg.polygon_mate(h));
						if (f_opp != pg.null_face() && !marked[f_opp])
							Queue.push(f_opp);
					}
				}
				++current_patch_id;
			}
			return current_patch_id;
		}

#ifdef MCAL_DEBUG
		struct IntersectPolylines
		{
			const std::vector<SMHalfedge>& sm_polylines;
			const std::vector<PGHalfedge>& pg_polylines;
			const std::vector<std::size_t>& poly_lengths;
			boost::dynamic_bitset<> to_skip;
			boost::dynamic_bitset<> to_skip_in_sm;
			boost::dynamic_bitset<> to_skip_in_pg;
			std::size_t nb_poly;

			IntersectPolylines(const std::vector<SMHalfedge>& sm_polylines_,
				const std::vector<PGHalfedge>& pg_polylines_,
				const std::vector<std::size_t>& poly_lengths_,
				std::size_t nb_poly_)
				: sm_polylines(sm_polylines_)
				, pg_polylines(pg_polylines_)
				, poly_lengths(poly_lengths_)
				, nb_poly(nb_poly_)
				, to_skip(nb_poly_, false)
				, to_skip_in_sm(nb_poly_, false)
				, to_skip_in_pg(nb_poly_, false)
			{
			}
		};

		// 如果交线两边的patch在当前的布尔操作中都不需要, 那么这条交线需要跳过.
		void fill_polyline_to_skip(IntersectPolylines& polylines,
			const std::vector<std::size_t>& sm_patch_ids,
			const std::map<PGFace, std::size_t>& pg_patch_ids,
			const boost::dynamic_bitset<> used_sm_patches,
			const boost::dynamic_bitset<> used_pg_patches,
			const SurfaceMesh& sm,
			const PolyhedronGrid& pg)
		{
			std::size_t nb_polylines = polylines.nb_poly;

			for (std::size_t i = 0; i < nb_polylines; ++i)
			{
				SMHalfedge h_sm = polylines.sm_polylines[i];
				PGHalfedge h_pg = polylines.pg_polylines[i];

				bool need_skip_in_sm = true;
				if (!is_border(h_sm, sm))
				{
					std::size_t patch_id = sm_patch_ids[face(h_sm, sm)];
					if (used_sm_patches.test(patch_id))
						need_skip_in_sm = false;
				}
				if (need_skip_in_sm && !is_border(opposite(h_sm, sm), sm))
				{
					std::size_t patch_id = sm_patch_ids[face(opposite(h_sm, sm), sm)];
					if (used_sm_patches.test(patch_id))
						need_skip_in_sm = false;
				}

				bool need_skip_in_pg = true;
				{
					typename std::map<PGFace, std::size_t>::const_iterator
						it = pg_patch_ids.find(pg.face(h_pg));
					assert(it != pg_patch_ids.end());
					std::size_t patch_id = it->second;
					if (used_pg_patches.test(patch_id))
						need_skip_in_pg = false;
				}
				if (need_skip_in_pg)
				{
					typename std::map<PGFace, std::size_t>::const_iterator
						it = pg_patch_ids.find(pg.face(pg.polygon_mate(h_pg)));
					assert(it != pg_patch_ids.end());
					std::size_t patch_id = it->second;
					if (used_pg_patches.test(patch_id))
						need_skip_in_pg = false;
				}

				if (need_skip_in_sm)
					polylines.to_skip_in_sm.set(i);
				if (need_skip_in_pg)
					polylines.to_skip_in_pg.set(i);
				if (need_skip_in_sm && need_skip_in_pg)
					polylines.to_skip.set(i);
			}
		}
#endif

		// surface mesh某个patch包含的所有几何元素.
		struct SMPatch
		{
			std::vector<SMFace> faces;         // patch包含的所有面.
			std::vector<SMVertex> interior_vertices; // patch不在交线上的点.
			std::vector<SMHalfedge> interior_edges;  // patch不再交线上的边(记录索引小的那条半边).
			std::vector<SMHalfedge> shared_edges;    // patch交线上的边.
			bool is_initialized;

			SMPatch() :is_initialized(false) {}
		};

		static
			void extract_patch_simplices(std::vector<SMFace>& faces,
				std::vector<SMVertex>& interior_vertices,
				std::vector<SMHalfedge>& interior_edges,
				std::vector<SMHalfedge>& shared_edges,
				const IntersectSMEdge& intersect_sm_edges,
				SurfaceMesh& sm)
		{
			for (SMFace f : faces)
			{
				for (SMHalfedge h : halfedges_around_face(halfedge(f, sm), sm))
				{
					if (!intersect_sm_edges.count(edge(h, sm)))
					{
						if (h < opposite(h, sm) || is_border(opposite(h, sm), sm))
							interior_edges.push_back(h);
					}
					else
						shared_edges.push_back(h);
				}
			}

			std::set<SMVertex> border_vertices;
			for (SMHalfedge h : shared_edges)
			{
				border_vertices.insert(source(h, sm));
				border_vertices.insert(target(h, sm));
			}

			for (SMHalfedge h : interior_edges)
			{
				if (!border_vertices.count(source(h, sm)))
					interior_vertices.push_back(source(h, sm));
				if (!border_vertices.count(target(h, sm)))
					interior_vertices.push_back(target(h, sm));
			}
		}

		// surface mesh若干patch, 注意, 由于算法只需要内部的patch, 外部的不再记录.
		struct SMPatchContainer
		{
			typedef std::size_t PatchId;
			std::map<PatchId, SMPatch> patches;

			SMPatchContainer(const std::vector<PatchId>& sm_patch_ids,
				const boost::dynamic_bitset<>& is_patch_inside_pg,
				const IntersectSMEdge& intersect_sm_edges,
				SurfaceMesh& sm)
			{
				for (SMFace f : faces(sm))
				{
					PatchId patch_id = sm_patch_ids[f];
					if (is_patch_inside_pg.test(patch_id))
						patches[patch_id].faces.push_back(f);
				}

				for (typename std::map<PatchId, SMPatch>::iterator
					it = patches.begin(); it != patches.end(); ++it)
				{
					if (it->second.is_initialized == false)
					{
						extract_patch_simplices(it->second.faces,
							it->second.interior_vertices,
							it->second.interior_edges,
							it->second.shared_edges,
							intersect_sm_edges, sm);
						it->second.is_initialized = true;
					}
				}
			}
		};

		SMHalfedge
			next_intersect_sm_halfedge(SMHalfedge h,
				const IntersectSMEdge& intersect_edge,
				const SurfaceMesh& sm)
		{
			assert(intersect_edge.count(edge(h, sm)));
			SMHalfedge nxt = next(h, sm);
			while (!intersect_edge.count(edge(nxt, sm)))
				nxt = next(opposite(nxt, sm), sm);
			assert(nxt != h);
			return nxt;
		}

		PGHalfedge
			next_intersect_pg_halfedge(PGHalfedge h,
				const IntersectPGEdge& intersect_edge,
				const PolyhedronGrid& pg)
		{
			assert(intersect_edge.count(pg.edge(h)));
			PGHalfedge nxt = pg.next(h);
			while (!intersect_edge.count(pg.edge(nxt)))
				nxt = pg.next(pg.polygon_mate(nxt));
			assert(nxt != h);
			return nxt;
		}

		// OutputBuilder functor. 由仿函数组织提取切割后网格的算法.
		template <typename NodeVector, typename NodeToSMVertex, typename NodeToPGVertex>
		void operator()(const NodeVector& nodes,
			const NodeToSMVertex& node_id_to_sm_vertex,
			const NodeToPGVertex& node_id_to_pg_vertex)
		{
			assert(node_id_to_sm_vertex.size() <= nodes.size());
			assert(node_id_to_pg_vertex.size() <= nodes.size());
			assert(sm_vertex_to_node_id.size() == pg_vertex_to_node_id.size());
			assert(sm_vertex_to_node_id.size() == nodes.size());

			int cur_exact_cnt = iExactComputeCount;
			int cur_interval_cnt = iIntervalComputeCount;

			std::clock_t after_refine = clock();
			std::cout << "This program cost " << (double)(after_refine - start) / CLOCKS_PER_SEC << "seconds" << std::endl;

			// part 1: 切开每个相交的cell.

			// 在之前的步骤中, 我们记录了每个cell的交线信息, 它由三部分组成:
			//   1) AnEdgePerPolyline: 交线的元信息(交点数等).
			//   2) IntersectSMEdge: 这个cell与表面网格交线包含的所有交边(sm中的).
			//   3) IntersectPGEdge: 这个cell与表面网格交线包含的所有交边(pg中的).
			// 
			// 依据这些信息, 将每个cell切开, 分为内和外两个部分.
			// 

#ifdef MCAL_DEBUG
			std::size_t sm_cnt = 0;
#endif

			for (typename CellPolylinesMap::iterator it = cell_polylines.begin();
				it != cell_polylines.end();
				++it)
			{
				PGCell c = it->first;
				assert(cell_intersect_sm_edges.count(c) && cell_intersect_pg_edges.count(c));

				AnEdgePerPolyline& polylines = it->second;
				IntersectSMEdge& intersect_sm_edge = cell_intersect_sm_edges[c];
				IntersectPGEdge& intersect_pg_edge = cell_intersect_pg_edges[c];

				// 交线将表面网格(一个cell也可视为表面网格)的面分成几个区域(patch),
				// 为每个区域的面赋一个patch_id.
				// (1) Assign a patch id to each facet indicating in which connected
				// component limited by intersection edges of the surface they are.

				// ... for sm
				std::vector<std::size_t> sm_patch_ids(num_faces(sm), NID);
				int nb_patches_sm = mark_sm_faces(intersect_sm_edge, sm_patch_ids, sm);

				std::vector<std::size_t> sm_patch_sizes(nb_patches_sm, 0);
				for (std::size_t i : sm_patch_ids)
				{
					if (i != NID)
						++sm_patch_sizes[i];
				}

				// ... for pg
				std::map<PGFace, std::size_t> pg_patch_ids;
				int nb_patches_pg = mark_pg_faces(c, intersect_pg_edge, pg_patch_ids, pg);

				std::vector<std::size_t> pg_patch_sizes(nb_patches_pg, 0);
				for (typename std::map<PGFace, std::size_t>::iterator
					it = pg_patch_ids.begin(); it != pg_patch_ids.end(); ++it)
				{
					if (it->second != NID)
						++pg_patch_sizes[it->second];
				}

				// 两个表面网格相交, 一条交边周围有4个面, 区分出内与外.
				// (2) Use the orientation around an edge to classify a patch

				boost::dynamic_bitset<> is_patch_inside_pg(nb_patches_sm, false);
				boost::dynamic_bitset<> is_patch_inside_sm(nb_patches_pg, false);
				boost::dynamic_bitset<> patch_not_set_sm(nb_patches_sm);
				boost::dynamic_bitset<> patch_not_set_pg(nb_patches_pg);
				patch_not_set_sm.set();
				patch_not_set_pg.set();

				for (typename AnEdgePerPolyline::iterator it_poly = polylines.begin();
					it_poly != polylines.end();
					++it_poly)
				{
					const std::pair<NodeId, NodeId>& ids = it_poly->first;

					SMHalfedge h_sm = it_poly->second.sm_hedge;
					PGHalfedge h_pg = it_poly->second.pg_hedge;

					assert(ids.first == sm_vertex_to_node_id[source(h_sm, sm)]);
					assert(ids.second == sm_vertex_to_node_id[target(h_sm, sm)]);
					assert(ids.first == pg_vertex_to_node_id[pg.source(h_pg)]);
					assert(ids.second == pg_vertex_to_node_id[pg.target(h_pg)]);

					assert(!is_border(h_sm, sm));

					// Sort the four triangle faces around their common edge.
					// we assume that the exterior of the volume is indicated by
					// counterclockwise oriented faces.

					// when looking from the side of ids.second,
					// the interior of the surface mesh is described 
					// by turning counterclockwise from p1 to p2.
					SMVertex p1 = target(next(opposite(h_sm, sm), sm), sm);
					SMVertex p2 = target(next(h_sm, sm), sm);
					// when looking from the side of ids.second,
					// the interior of the polyhedron grid is described
					// by turning from q1 to q2.
					PGVertex q1 = pg.target(pg.next(pg.polygon_mate(h_pg)));
					PGVertex q2 = pg.target(pg.next(h_pg));

					NodeId index_p1 = get_node_id(p1, sm_vertex_to_node_id);
					NodeId index_p2 = get_node_id(p2, sm_vertex_to_node_id);
					NodeId index_q1 = get_node_id(q1, pg_vertex_to_node_id);
					NodeId index_q2 = get_node_id(q2, pg_vertex_to_node_id);

					std::size_t patch_id_p1 = sm_patch_ids[face(opposite(h_sm, sm), sm)];
					std::size_t patch_id_p2 = sm_patch_ids[face(h_sm, sm)];
					std::size_t patch_id_q1 = pg_patch_ids[pg.face(pg.polygon_mate(h_pg))];
					std::size_t patch_id_q2 = pg_patch_ids[pg.face(h_pg)];

					// info on whether the patches were already classified.
					std::bitset<4> patch_status_was_not_already_set;
					std::bitset<4> prev_bitvalue;
					// for sm
					patch_status_was_not_already_set[0] = patch_not_set_sm.test(patch_id_p1);
					patch_status_was_not_already_set[1] = patch_not_set_sm.test(patch_id_p2);
					prev_bitvalue[0] = is_patch_inside_pg.test(patch_id_p1);
					prev_bitvalue[1] = is_patch_inside_pg.test(patch_id_p2);
					// for pg
					patch_status_was_not_already_set[2] = patch_not_set_pg.test(patch_id_q1);
					patch_status_was_not_already_set[3] = patch_not_set_pg.test(patch_id_q2);
					prev_bitvalue[2] = is_patch_inside_sm.test(patch_id_q1);
					prev_bitvalue[3] = is_patch_inside_sm.test(patch_id_q2);

					// 如果全部patch都已经标记, 无需再做一遍, 跳过即可.
					if (!patch_status_was_not_already_set[0] &&
						!patch_status_was_not_already_set[1] &&
						!patch_status_was_not_already_set[2] &&
						!patch_status_was_not_already_set[3])
						continue;

					// 更新相关的patch status, 因为接下来要设置它.
					patch_not_set_sm.reset(patch_id_p1);
					patch_not_set_sm.reset(patch_id_p2);
					patch_not_set_pg.reset(patch_id_q1);
					patch_not_set_pg.reset(patch_id_q2);

					// 保证这四个点各不相同.
					assert(
						(index_p1 == NID ? nodes.to_exact(get(vpm_sm, p1)) : nodes.exact_node(index_p1)) !=
						(index_q1 == NID ? nodes.to_exact(vpm_pg[q1]) : nodes.exact_node(index_q1))
						&&
						(index_p2 == NID ? nodes.to_exact(get(vpm_sm, p2)) : nodes.exact_node(index_p2)) !=
						(index_q1 == NID ? nodes.to_exact(vpm_pg[q1]) : nodes.exact_node(index_q1))
						&&
						(index_p1 == NID ? nodes.to_exact(get(vpm_sm, p1)) : nodes.exact_node(index_p1)) !=
						(index_q2 == NID ? nodes.to_exact(vpm_pg[q2]) : nodes.exact_node(index_q2))
						&&
						(index_p2 == NID ? nodes.to_exact(get(vpm_sm, p2)) : nodes.exact_node(index_p2)) !=
						(index_q2 == NID ? nodes.to_exact(vpm_pg[q2]) : nodes.exact_node(index_q2))
					);

					bool q1_is_inside_p1p2 = decide_q_inside_p(
						ids.first, ids.second,
						index_p1, index_p2, index_q1,
						p1, p2, q1,
						vpm_sm, vpm_pg,
						nodes);

					bool q2_is_inside_p1p2 = decide_q_inside_p(
						ids.first, ids.second,
						index_p1, index_p2, index_q2,
						p1, p2, q2,
						vpm_sm, vpm_pg,
						nodes);

					// 一般情况下, q1和q2应该是一内一外, 否则就需要特殊情形的处理.
					assert(q1_is_inside_p1p2 != q2_is_inside_p1p2);

					if (q1_is_inside_p1p2)
					{
						// 若q1在内而q2在外, 则p2在多面体cell内部.
						is_patch_inside_sm.set(patch_id_q1);
						is_patch_inside_pg.set(patch_id_p2);
					}
					else
					{
						// 若q1在外而q2在内, 则p1在多面体网格cell内部.
						is_patch_inside_sm.set(patch_id_q2);
						is_patch_inside_pg.set(patch_id_p1);
					}

					assert(patch_status_was_not_already_set[0] || prev_bitvalue[0] == is_patch_inside_pg[patch_id_p1]);
					assert(patch_status_was_not_already_set[1] || prev_bitvalue[1] == is_patch_inside_pg[patch_id_p2]);
					assert(patch_status_was_not_already_set[2] || prev_bitvalue[2] == is_patch_inside_sm[patch_id_q1]);
					assert(patch_status_was_not_already_set[3] || prev_bitvalue[3] == is_patch_inside_sm[patch_id_q2]);
				}

				// 至此, 所有的patch都必须标记好, 否则为错误.
				assert(patch_not_set_sm.none());
				assert(patch_not_set_pg.none());

				// info 1
				std::map<BO_Type, boost::dynamic_bitset<>> used_sm_patches;
				std::map<BO_Type, boost::dynamic_bitset<>> used_pg_patches;
				used_sm_patches[INTERSECTION] = is_patch_inside_pg;
				used_pg_patches[INTERSECTION] = is_patch_inside_sm;
				used_sm_patches[PG_MINUS_SM] = is_patch_inside_pg;
				used_pg_patches[PG_MINUS_SM] = ~is_patch_inside_sm;

				// info 2
				SMPatchContainer sm_patches(sm_patch_ids, is_patch_inside_pg, intersect_sm_edge, sm);

				// info 3
				std::vector<SMHalfedge> sm_polylines;
				std::vector<PGHalfedge> pg_polylines;
				std::vector<std::size_t> poly_lengths;
				for (typename AnEdgePerPolyline::iterator it_poly = polylines.begin();
					it_poly != polylines.end();
					++it_poly)
				{
					const PolylineInfo& polyline_info = it_poly->second;

					SMHalfedge h_sm = polyline_info.sm_hedge;
					PGHalfedge h_pg = polyline_info.pg_hedge;

					if (polyline_info.is_reverse)
					{
						h_sm = opposite(h_sm, sm);
						h_pg = pg.polygon_mate(h_pg);
					}

					sm_polylines.push_back(h_sm);
					pg_polylines.push_back(h_pg);
					poly_lengths.push_back(polyline_info.node_num + 1);
				}

#ifdef MCAL_DEBUG
				// 检测一下在提取cell内外的时候, 是否有不需要的交线.(检测完是没有)
				std::size_t nb_poly = sm_polylines.size();
				IntersectPolylines in_polys(sm_polylines, pg_polylines, poly_lengths, nb_poly);
				fill_polyline_to_skip(
					in_polys, sm_patch_ids, pg_patch_ids,
					used_sm_patches[INTERSECTION], used_pg_patches[INTERSECTION],
					sm, pg);

				IntersectPolylines out_polys(sm_polylines, pg_polylines, poly_lengths, nb_poly);
				fill_polyline_to_skip(
					out_polys, sm_patch_ids, pg_patch_ids,
					used_sm_patches[PG_MINUS_SM], used_pg_patches[PG_MINUS_SM],
					sm, pg);

				assert(in_polys.to_skip.none() && out_polys.to_skip.none());
#endif

				// (3) 收集sm与pg的几何元素的对应关系.
				std::unordered_map<SMVertex, PGVertex> sm_vertex_to_pg_vertex;
				std::unordered_map<PGFace, std::vector<PGVertex>> f_vseq;
				std::map<SortedPair<PGVertex>, PGEdge> endpoint_to_edge;

				std::size_t nb_poly = sm_polylines.size();
				for (std::size_t i = 0; i < nb_poly; ++i)
				{
					SMHalfedge h_sm = sm_polylines[i];
					PGHalfedge h_pg = pg_polylines[i];
					std::size_t nb_segment = poly_lengths[i];

					for (std::size_t k = 0;;)
					{
						SMVertex srcv_sm = source(h_sm, sm);
						PGVertex srcv_pg = pg.source(h_pg);
						if (!sm_vertex_to_pg_vertex.count(srcv_sm))
							sm_vertex_to_pg_vertex.insert(std::make_pair(srcv_sm, srcv_pg));

						SMVertex tgtv_sm = target(h_sm, sm);
						PGVertex tgtv_pg = pg.target(h_pg);
						if (!sm_vertex_to_pg_vertex.count(tgtv_sm))
							sm_vertex_to_pg_vertex.insert(std::make_pair(tgtv_sm, tgtv_pg));

						assert(get(vpm_sm, srcv_sm) == vpm_pg[srcv_pg]);
						assert(get(vpm_sm, tgtv_sm) == vpm_pg[tgtv_pg]);

						SortedPair<PGVertex> vpair(srcv_pg, tgtv_pg);
						assert(!endpoint_to_edge.count(vpair));
						endpoint_to_edge.insert(std::make_pair(vpair, pg.edge(h_pg)));

						if (++k == nb_segment)
							break;
						h_sm = next_intersect_sm_halfedge(h_sm, intersect_sm_edge, sm);
						h_pg = next_intersect_pg_halfedge(h_pg, intersect_pg_edge, pg);
					}
				}

				for (const std::pair<std::size_t, SMPatch>& kv : sm_patches.patches)
				{
					const SMPatch& sm_patch = kv.second;
					for (SMVertex v_sm : sm_patch.interior_vertices)
					{
						if (!sm_vertex_to_pg_vertex.count(v_sm))
						{
							PGVertex v_pg = pg.add_vertex(get(vpm_sm, v_sm));
							sm_vertex_to_pg_vertex.insert(std::make_pair(v_sm, v_pg));
						}
					}

					for (SMFace f_sm : sm_patch.faces)
					{
						PGFace f_pg = pg.add_face();
						std::vector<PGVertex> vseq;
						for (SMHalfedge h_sm : halfedges_around_face(halfedge(f_sm, sm), sm))
						{
							assert(sm_vertex_to_pg_vertex.count(source(h_sm, sm)));
							assert(sm_vertex_to_pg_vertex.count(target(h_sm, sm)));

							PGVertex srcv = sm_vertex_to_pg_vertex[source(h_sm, sm)];
							PGVertex tgtv = sm_vertex_to_pg_vertex[target(h_sm, sm)];
							vseq.push_back(srcv);

							SortedPair<PGVertex> vpair(srcv, tgtv);
							if (!endpoint_to_edge.count(vpair))
								endpoint_to_edge.insert(std::make_pair(vpair, pg.add_edge()));
						}
						// 确保没有重复的元素.
						assert(std::set<PGVertex>(vseq.begin(), vseq.end()).size() == vseq.size());
						f_vseq.insert(std::make_pair(f_pg, vseq));
					}
				}

				// (4) 添加一个cell, 我们规定原来的cell一分为二, 原来cell保存外部的拓扑, 新添加的保存内部拓扑.
				PGCell new_cell = pg.add_cell();
				pg.copy_cinfo(c, new_cell);

				for (typename std::unordered_map<PGFace, std::vector<PGVertex>>::iterator
					fit = f_vseq.begin(); fit != f_vseq.end(); ++fit)
				{
					pg.add_triangle_to_cell(fit->first, new_cell, fit->second, endpoint_to_edge, false);
					pg.add_triangle_to_cell(fit->first, c, fit->second, endpoint_to_edge, true);
				}

				for (typename std::map<PGFace, std::size_t>::iterator
					fit = pg_patch_ids.begin(); fit != pg_patch_ids.end(); ++fit)
				{
					PGFace f_pg = fit->first;
					std::size_t patch_id = fit->second;
					// 需要从原来cell中删除内部的面, 将其移动到新的cell上.
					if (is_patch_inside_sm.test(patch_id))
						pg.move_face(f_pg, c, new_cell);
				}

#ifdef MCAL_DEBUG   // 将两个cell转换为表面网格写出, 看一下效果.
				SurfaceMesh mesh_in, mesh_out;
				pg.convert_cell_to_surface_mesh(new_cell, mesh_in);
				pg.convert_cell_to_surface_mesh(c, mesh_out);

				std::string out_path = "E:\\tmp";
				CGAL::IO::write_polygon_mesh(out_path + "\\" + std::to_string(sm_cnt) +
					"_in.off", mesh_in, CGAL::parameters::stream_precision(17));
				CGAL::IO::write_polygon_mesh(out_path + "\\" + std::to_string(sm_cnt) +
					"_out.off", mesh_out, CGAL::parameters::stream_precision(17));
				++sm_cnt;
#endif
			}

			pg.mark_all_cells();

			std::cout << "构建过程----------------------------" << std::endl;
			std::cout << "精确计算的次数：" << iExactComputeCount - cur_exact_cnt << std::endl;
			std::cout << "区间代数的次数：" << iIntervalComputeCount - cur_interval_cnt << std::endl;
			std::cout << "----------------------------------------\n\n";

		}

	};

} // namespace MCAL

#endif