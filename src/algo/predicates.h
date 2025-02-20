// $License: GPL-3.0
// $Author: Liu Zhiqiang(from Beihang University)
// $Email: lzhq0930@gmail.com(liuzhiqiang@buaa.edu.cn)
// 


#ifndef MCAL_ALGO_PREDICATES_H
#define MCAL_ALGO_PREDICATES_H

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
	{}

	bool operator()(std::size_t i, std::size_t j) const
	{
		// 当且仅当j在i和h_sm的target之间时, 返回true.
		return CGAL::collinear_are_strictly_ordered_along_line(
			nodes.to_exact(get(vpm_sm, target(h_sm, sm))),
			nodes.exact_node(j),
			nodes.exact_node(i));
	}
};

template <typename PolyhedralMesh, typename VPMPM, typename NodeVector>
struct Less_along_a_pm_halfedge
{
	// typedefs
	typedef typename PolyhedralMesh::Halfedge    PMHalfedge;

	// data members
	PMHalfedge h_pm;
	PolyhedralMesh& pm;
	VPMPM& vpm_pm;
	NodeVector& nodes;

	Less_along_a_pm_halfedge(PMHalfedge h_pm_,
							 PolyhedralMesh& pm_,
							 VPMPM& vpm_pm_,
							 NodeVector& nodes_)
		:h_pm(h_pm_), pm(pm_), vpm_pm(vpm_pm_), nodes(nodes_)
	{}

	bool operator()(std::size_t i, std::size_t j) const
	{
		// 当且仅当j在i和h_pm的target之间时, 返回true.
		return CGAL::collinear_are_strictly_ordered_along_line(
			nodes.to_exact(vpm_pm[pm.target(h_pm)]),
			nodes.exact_node(j),
			nodes.exact_node(i));
	}
};

template <typename NodeId,
		  typename NodeVector,
		  typename SMVertex,
		  typename PMVertex,
		  typename VPMSM,
		  typename VPMPM>
bool decide_q_inside_p(NodeId o_prime_index,
					   NodeId o_index,
					   NodeId p1_index,
					   NodeId p2_index,
					   NodeId q_index,
					   SMVertex p1,
					   SMVertex p2,
					   PMVertex q,
					   const VPMSM& vpm_sm,
					   const VPMPM& vpm_pm,
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
				q_index == NID ? nodes.to_exact(vpm_pm[q])
								: nodes.exact_node(q_index));
}

}	// namespace MCAL

#endif