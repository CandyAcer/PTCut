// $Intro: 本文件主要实现算法需要的几何断言.
// $License: GPL-3.0-or-later
// 
// 
// $Author: Liu Zhiqiang


#ifndef MCAL_ALGORITHM_PREDICATES_H
#define MCAL_ALGORITHM_PREDICATES_H


#include <CGAL/Polygon_mesh_processing/internal/Corefinement/predicates.h>


namespace MCAL   // Mesh Cut Algorithm Libaray, 存放多面体网格切割算法的实现
{

	template <typename SurfaceMesh, typename VPMSM, typename NodeVector>
	struct Less_along_a_sm_halfedge
	{
		// typedefs
		typedef boost::graph_traits<SurfaceMesh>        GT;
		typedef typename GT::halfedge_descriptor        SMHalfedge;

		// data members
		SMHalfedge h_sm;
		SurfaceMesh& sm;
		VPMSM& vpm_sm;
		NodeVector& nodes;

		Less_along_a_sm_halfedge(SMHalfedge h_sm_,
			SurfaceMesh& sm_,
			VPMSM& vpm_sm_,
			NodeVector& nodes_)
			: h_sm(h_sm_), sm(sm_), vpm_sm(vpm_sm_), nodes(nodes_)
		{
		}

		bool operator()(std::size_t i, std::size_t j) const
		{
			// 当且仅当j在i和h_sm的target之间时, 返回true.
			return CGAL::collinear_are_strictly_ordered_along_line(
				nodes.to_exact(get(vpm_sm, target(h_sm, sm))),
				nodes.exact_node(j),
				nodes.exact_node(i));
		}
	};

	template <typename PolyhedronGrid, typename VPMPG, typename NodeVector>
	struct Less_along_a_pg_halfedge
	{
		// typedefs
		typedef typename PolyhedronGrid::Halfedge    PGHalfedge;

		// data members
		PGHalfedge h_pg;
		PolyhedronGrid& pg;
		VPMPG& vpm_pg;
		NodeVector& nodes;

		Less_along_a_pg_halfedge(PGHalfedge h_pg_,
			PolyhedronGrid& pg_,
			VPMPG& vpm_pg_,
			NodeVector& nodes_)
			:h_pg(h_pg_), pg(pg_), vpm_pg(vpm_pg_), nodes(nodes_)
		{
		}

		bool operator()(std::size_t i, std::size_t j) const
		{
			// 当且仅当j在i和h_pg的target之间时, 返回true.
			return CGAL::collinear_are_strictly_ordered_along_line(
				nodes.to_exact(vpm_pg[pg.target(h_pg)]),
				nodes.exact_node(j),
				nodes.exact_node(i));
		}
	};

	template <typename NodeId,
		typename NodeVector,
		typename SMVertex,
		typename PGVertex,
		typename VPMSM,
		typename VPMPG>
	bool decide_q_inside_p(NodeId o_prime_index,
		NodeId o_index,
		NodeId p1_index,
		NodeId p2_index,
		NodeId q_index,
		SMVertex p1,
		SMVertex p2,
		PGVertex q,
		const VPMSM& vpm_sm,
		const VPMPG& vpm_pg,
		const NodeVector& nodes)
	{
		const NodeId NID((std::numeric_limits<NodeId>::max)());
		using namespace CGAL::Polygon_mesh_processing::Corefinement;

		return sorted_around_edge<typename NodeVector::ExactKernel>(
			nodes.exact_node(o_prime_index),
			nodes.exact_node(o_index),
			p1_index == NID ? nodes.to_exact(get(vpm_sm, p1))
			: nodes.exact_node(p1_index),
			p2_index == NID ? nodes.to_exact(get(vpm_sm, p2))
			: nodes.exact_node(p2_index),
			q_index == NID ? nodes.to_exact(vpm_pg[q])
			: nodes.exact_node(q_index));
	}

}	// namespace MCAL

#endif