// $License: GPL-3.0
// $Author: Liu Zhiqiang(from Beihang University)
// $Email: lzhq0930@gmail.com(liuzhiqiang@buaa.edu.cn)
// 


#ifndef MCAL_ALGO_INTERSECTION_NODES_H
#define MCAL_ALGO_INTERSECTION_NODES_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Kernel_traits.h>

#include <boost/property_map/property_map.hpp>

#include <assert.h>

namespace MCAL    // Mesh Cut Algorithm Libaray, 存放多面体网格切割算法的实现
{
// 管理所有交点.
// Store an exact copy of the points so that we can answer exactly predicates.
template <typename SurfaceMesh,
		  typename PolyhedralMesh,
		  typename VPMSM,
		  typename VPMPM>
class IntersectionNodes
{
public:
	typedef CGAL::Exact_predicates_exact_constructions_kernel   ExactKernel;

private:
	typedef typename boost::property_traits<VPMSM>::value_type  Point3;
	static_assert((std::is_same<Point3, typename PolyhedralMesh::Point>::value),
				  "输入网格的点类型不一致");

	typedef typename CGAL::Kernel_traits<Point3>::Kernel        InputKernel;

	typedef CGAL::Cartesian_converter<InputKernel, ExactKernel>  DoubleToExact;
	typedef CGAL::Cartesian_converter<ExactKernel, InputKernel>  ExactToDouble;

	typedef boost::graph_traits<SurfaceMesh>                             GT;
	typedef typename GT::halfedge_descriptor                     SMHalfedge;
	typedef typename GT::face_descriptor                             SMFace;
	typedef typename GT::vertex_descriptor                         SMVertex;

	typedef typename PolyhedralMesh::Halfedge                    PMHalfedge;
	typedef typename PolyhedralMesh::VertexIndex                   PMVertex;
	typedef typename PolyhedralMesh::FaceIndex                       PMFace;

	typedef std::vector<ExactKernel::Point_3>                    ExactNodes;

	// data members
	ExactNodes              enodes;              // 以精确形式存储所有的交点
	// 提供精确/非精确的相互转化
	DoubleToExact           double_to_exact;
	ExactToDouble           exact_to_double;

public:
	SurfaceMesh& sm;
	PolyhedralMesh& pm;
	VPMSM& vpm_sm;
	VPMPM& vpm_pm;

	IntersectionNodes(SurfaceMesh& sm_,
					  PolyhedralMesh& pm_,
					  VPMSM& vpm_sm_,
					  VPMPM& vpm_pm_)
		:sm(sm_), pm(pm_), vpm_sm(vpm_sm_), vpm_pm(vpm_pm_)
	{}

	Point3 operator[](std::size_t i) const
	{
		return exact_to_double(enodes[i]);
	}

	const ExactKernel::Point_3
	exact_node(std::size_t i) const
	{
		return enodes[i];
	}

	static ExactKernel::Point_3
	to_exact(const Point3& p)
	{
		return ExactKernel::Point_3(p.x(), p.y(), p.z());
	}

	std::size_t size() const { return enodes.size(); }

	// Case 1: 交点是网格已有的点, 只需转换为精确表示进行存储即可
	void add_new_node(const Point3& p)
	{
		enodes.push_back(to_exact(p));
	}

	void add_new_node(const ExactKernel::Point_3& p)
	{
		const ExactKernel::Approximate_kernel::Point_3 p_approx = p.approx();
		const double precision =
			CGAL::Lazy_exact_nt<ExactKernel::FT>::get_relative_precision_of_to_double();

		if (!CGAL::has_smaller_relative_precision(p_approx.x(), precision) ||
			!CGAL::has_smaller_relative_precision(p_approx.y(), precision) ||
			!CGAL::has_smaller_relative_precision(p_approx.z(), precision))
		{
			p.exact();
		}
		enodes.push_back(p);
	}

	// Case 2: 计算SurfaceMesh的edge与PolyhedralMesh的face的交点
	void add_new_node(SMHalfedge h_sm,
					  PMFace f_pm,
					  SurfaceMesh& sm,
					  PolyhedralMesh& pm,
					  VPMSM& vpm_sm,
					  VPMPM& vpm_pm)
	{
		assert(pm.is_triangle(f_pm));
		PMHalfedge h_pm = pm.halfedge(f_pm);
		// 此处是精确构造, 避免因浮点误差出现错误.
		add_new_node(
			ExactKernel::Construct_plane_line_intersection_point_3()(
				to_exact(vpm_pm[pm.source(h_pm)]),
				to_exact(vpm_pm[pm.target(h_pm)]),
				to_exact(vpm_pm[pm.target(pm.next(h_pm))]),
				to_exact(get(vpm_sm, source(h_sm, sm))),
				to_exact(get(vpm_sm, target(h_sm, sm))))
		);
	}

	// Case 3: 计算PolyhedralMesh的edge与SurfaceMesh的face的交点
	void add_new_node(PMHalfedge h_pm,
					  SMFace f_sm,
					  PolyhedralMesh& pm,
					  SurfaceMesh& sm,
					  VPMPM& vpm_pm,
					  VPMSM& vpm_sm)
	{
		SMHalfedge h_sm = halfedge(f_sm, sm);
		add_new_node(
			ExactKernel::Construct_plane_line_intersection_point_3()(
				to_exact(get(vpm_sm, source(h_sm, sm))),
				to_exact(get(vpm_sm, target(h_sm, sm))),
				to_exact(get(vpm_sm, target(next(h_sm, sm), sm))),
				to_exact(vpm_pm[pm.source(h_pm)]),
				to_exact(vpm_pm[pm.target(h_pm)]))
		);
	}

	// 精确转换为非精确, 写出的点可能不对.
	void call_put(VPMSM& vpm_sm, SMVertex v_sm, std::size_t i, SurfaceMesh& sm)
	{
		put(vpm_sm, v_sm, exact_to_double(enodes[i]));
	}

	void call_put(VPMPM& vpm_pm, PMVertex v_pm, std::size_t i, PolyhedralMesh& pm)
	{
		vpm_pm[v_pm] = exact_to_double(enodes[i]);
	}

	template <typename NodeToSMVertex, typename NodeToPMVertex>
	void finalize(NodeToSMVertex& node_id_to_sm_vertex, NodeToPMVertex& node_id_to_pm_vertex)
	{
		typedef std::size_t NodeId;
		for (NodeId node_id = 0; node_id < enodes.size(); ++node_id)
		{
			Point3 p = exact_to_double(enodes[node_id]);

			SMVertex v_sm = node_id_to_sm_vertex[node_id];
			if (v_sm != GT::null_vertex())
				put(vpm_sm, v_sm, p);

			PMVertex v_pm = node_id_to_pm_vertex[node_id];
			if (v_pm != pm.null_vertex())
				vpm_pm[v_pm] = p;
		}
	}

	void check_no_duplicates() const
	{
		std::set<typename ExactKernel::Point_3> node_set(enodes.begin(), enodes.end());
		assert(enodes.size() == node_set.size());
	}
};

}	// namespace MCAL

#endif