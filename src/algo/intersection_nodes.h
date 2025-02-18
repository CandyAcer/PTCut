// $Intro: 本文件实现对交点的统一管理, 以精确形式存储, 提供精确的几何断言.
// $License: GPL-3.0-or-later
// 
// 
// $Author: Liu Zhiqiang



#ifndef MCAL_ALGORITHM_INTERSECTION_NODES_H
#define MCAL_ALGORITHM_INTERSECTION_NODES_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Kernel_traits.h>

#include <boost/property_map/property_map.hpp>

#include <assert.h>


namespace MCAL    // Mesh Cut Algorithm Libaray, 存放多面体网格切割算法的实现
{
	// 管理所有交点.
	// Store an exact copy of the points so that we can answer exactly predicates.
	template <typename SurfaceMesh,
		typename PolyhedronGrid,
		typename VPMSM,
		typename VPMPG>
	class IntersectionNodes
	{
	public:
		typedef CGAL::Exact_predicates_exact_constructions_kernel   ExactKernel;

	private:
		typedef typename boost::property_traits<VPMSM>::value_type  Point3;
		static_assert((std::is_same<Point3, typename PolyhedronGrid::Point>::value)
			, "输入网格的点类型不一致");

		typedef typename CGAL::Kernel_traits<Point3>::Kernel        InputKernel;

		typedef CGAL::Cartesian_converter<InputKernel, ExactKernel>  DoubleToExact;
		typedef CGAL::Cartesian_converter<ExactKernel, InputKernel>  ExactToDouble;

		typedef boost::graph_traits<SurfaceMesh>                             GT;
		typedef typename GT::halfedge_descriptor                     SMHalfedge;
		typedef typename GT::face_descriptor                             SMFace;
		typedef typename GT::vertex_descriptor                         SMVertex;

		typedef typename PolyhedronGrid::Halfedge                    PGHalfedge;
		typedef typename PolyhedronGrid::VertexIndex                   PGVertex;
		typedef typename PolyhedronGrid::FaceIndex                       PGFace;

		typedef std::vector<ExactKernel::Point_3>                    ExactNodes;

		// data members
		ExactNodes              enodes;              // 以精确形式存储所有的交点
		// 提供精确/非精确的相互转化
		DoubleToExact           double_to_exact;
		ExactToDouble           exact_to_double;

	public:
		SurfaceMesh& sm;
		PolyhedronGrid& pg;
		VPMSM& vpm_sm;
		VPMPG& vpm_pg;

		IntersectionNodes(SurfaceMesh& sm_,
			PolyhedronGrid& pg_,
			VPMSM& vpm_sm_,
			VPMPG& vpm_pg_)
			:sm(sm_), pg(pg_), vpm_sm(vpm_sm_), vpm_pg(vpm_pg_)
		{
		}

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

		// Case 2: 计算SurfaceMesh的edge与PolyhedronGrid的face的交点
		void add_new_node(SMHalfedge h_sm,
			PGFace f_pg,
			SurfaceMesh& sm,
			PolyhedronGrid& pg,
			VPMSM& vpm_sm,
			VPMPG& vpm_pg)
		{
			assert(pg.is_triangle(f_pg));
			PGHalfedge h_pg = pg.halfedge(f_pg);
			// 此处是精确构造, 避免因浮点误差出现错误.
			add_new_node(
				ExactKernel::Construct_plane_line_intersection_point_3()(
					to_exact(vpm_pg[pg.source(h_pg)]),
					to_exact(vpm_pg[pg.target(h_pg)]),
					to_exact(vpm_pg[pg.target(pg.next(h_pg))]),
					to_exact(get(vpm_sm, source(h_sm, sm))),
					to_exact(get(vpm_sm, target(h_sm, sm))))
			);
		}

		// Case 3: 计算PolyhedronGrid的edge与SurfaceMesh的face的交点
		void add_new_node(PGHalfedge h_pg,
			SMFace f_sm,
			PolyhedronGrid& pg,
			SurfaceMesh& sm,
			VPMPG& vpm_pg,
			VPMSM& vpm_sm)
		{
			SMHalfedge h_sm = halfedge(f_sm, sm);
			add_new_node(
				ExactKernel::Construct_plane_line_intersection_point_3()(
					to_exact(get(vpm_sm, source(h_sm, sm))),
					to_exact(get(vpm_sm, target(h_sm, sm))),
					to_exact(get(vpm_sm, target(next(h_sm, sm), sm))),
					to_exact(vpm_pg[pg.source(h_pg)]),
					to_exact(vpm_pg[pg.target(h_pg)]))
			);
		}

		// 精确转换为非精确, 写出的点可能不对.
		void call_put(VPMSM& vpm_sm, SMVertex v_sm, std::size_t i, SurfaceMesh& sm)
		{
			put(vpm_sm, v_sm, exact_to_double(enodes[i]));
		}

		void call_put(VPMPG& vpm_pg, PGVertex v_pg, std::size_t i, PolyhedronGrid& pg)
		{
			vpm_pg[v_pg] = exact_to_double(enodes[i]);
		}

		template <typename NodeToSMVertex, typename NodeToPGVertex>
		void finalize(NodeToSMVertex& node_id_to_sm_vertex, NodeToPGVertex& node_id_to_pg_vertex)
		{
			typedef std::size_t NodeId;
			for (NodeId node_id = 0; node_id < enodes.size(); ++node_id)
			{
				Point3 p = exact_to_double(enodes[node_id]);

				SMVertex v_sm = node_id_to_sm_vertex[node_id];
				if (v_sm != GT::null_vertex())
					put(vpm_sm, v_sm, p);

				PGVertex v_pg = node_id_to_pg_vertex[node_id];
				if (v_pg != pg.null_vertex())
					vpm_pg[v_pg] = p;
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