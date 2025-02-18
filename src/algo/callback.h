// $Intro: ���ļ�ʵ���ཻ���Ļص�����, ��¼�����ཻ��pair.
// $License: GPL-3.0-or-later
// 
// 
// $Author: Liu Zhiqiang



#ifndef MCAL_ALGORITHM_INTERSECTION_CALLBACK_H
#define MCAL_ALGORITHM_INTERSECTION_CALLBACK_H


#include <CGAL/Box_intersection_d/Box_with_info_d.h>
#include <boost/graph/graph_traits.hpp>


namespace MCAL   // Mesh Cut Algorithm Libaray, ��Ŷ����������и��㷨��ʵ��
{

	// helper class
	// CGAL��box_intersection_d�㷨ʹ�õ��Ǵ�info��AABB��Χ��, �и��㷨��info������İ��.
	// ��������, SurfaceMesh��PolyhedronGrid�İ���ǲ�һ����.
	// ����ʹ��BoxInfo�Զ��߽���ͳһ, �Ӷ����̶ȵظ���CGAL::box_intersection_d�㷨.
	//
	template <typename SurfaceMesh, typename PolyhedronGrid>
	struct BoxInfo
	{
		typedef typename SurfaceMesh::size_type  size_type;
		static_assert(std::is_same<size_type, typename PolyhedronGrid::size_type>::value,
			"��������������������size_type���Ͳ�һ�£�");

		typedef typename PolyhedronGrid::Halfedge PGHalfedge;

		size_type primitive_index;
		size_type array_offset;

		BoxInfo(size_type idx = std::numeric_limits<size_type>::max(),
			size_type off = std::numeric_limits<size_type>::max())
			:primitive_index(idx), array_offset(off)
		{
		}

		BoxInfo(const PGHalfedge& h)
			:primitive_index(h.he_cell), array_offset(h.he_off)
		{
		}

		BoxInfo(const BoxInfo& rhs) { operator=(rhs); }

		BoxInfo& operator=(const BoxInfo& rhs)
		{
			primitive_index = rhs.primitive_index;
			array_offset = rhs.array_offset;
			return *this;
		}
	};

	// ��⵽SurfaceMesh��edge��PolyhedronGrid��face��Χ���ཻʱ, ���ô�callback.
	template <typename SurfaceMesh,
		typename PolyhedronGrid,
		typename EdgeToFaces>
	class CollectSMEdgeToPGFaces
	{
		// typedefs
		typedef boost::graph_traits<SurfaceMesh>             GT;
		typedef typename GT::halfedge_descriptor             SMHalfedge;

		typedef typename PolyhedronGrid::FaceIndex           PGFace;
		typedef typename PolyhedronGrid::CellIndex           PGCell;
		typedef typename PolyhedronGrid::Halfedge            PGHalfedge;

		typedef BoxInfo<SurfaceMesh, PolyhedronGrid>                       Info;
		typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, Info>  Box;

		// data members
		SurfaceMesh& sm;
		PolyhedronGrid& pg;
		EdgeToFaces& sm_edge_to_pg_faces;

	public:
		CollectSMEdgeToPGFaces(SurfaceMesh& sm_,
			PolyhedronGrid& pg_,
			EdgeToFaces& edge_to_faces_)
			:sm(sm_), pg(pg_), sm_edge_to_pg_faces(edge_to_faces_)
		{
		}

		void operator()(const Box& sm_edge_box, const Box& pg_face_box)
		{
			Info sm_info = sm_edge_box.info();
			SMHalfedge h_sm(sm_info.primitive_index);

			Info pg_info = pg_face_box.info();
			PGHalfedge h_pg(PGCell(pg_info.primitive_index), pg_info.array_offset);
			PGFace f_pg = pg.face(h_pg);

			sm_edge_to_pg_faces[edge(h_sm, sm)].insert(f_pg);
		}
	};

	// ��⵽PolyhedronGrid��cell��SurfaceMesh��face��Χ���ཻʱ, ���ô�callback.
	template <typename PolyhedronGrid,
		typename SurfaceMesh,
		typename CellToFaces>
	class CollectPGCellToSMFaces
	{
	private:
		PolyhedronGrid& pg;
		SurfaceMesh& sm;
		CellToFaces& pg_cell_to_sm_faces;

		typedef boost::graph_traits<SurfaceMesh>             GT;
		typedef typename GT::halfedge_descriptor             SMHalfedge;

		typedef typename PolyhedronGrid::Halfedge            PGHalfedge;
		typedef typename PolyhedronGrid::CellIndex           PGCell;

		typedef BoxInfo<SurfaceMesh, PolyhedronGrid>         Info;
		typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, Info> Box;

	public:
		CollectPGCellToSMFaces(PolyhedronGrid& pg_,
			SurfaceMesh& sm_,
			CellToFaces& cell_to_faces)
			:pg(pg_), sm(sm_), pg_cell_to_sm_faces(cell_to_faces)
		{
		}

		void operator()(const Box& pg_cell_box, const Box& sm_face_box)
		{
			Info pg_info = pg_cell_box.info();
			PGHalfedge h_pg(PGCell(pg_info.primitive_index), pg_info.array_offset);

			Info sm_info = sm_face_box.info();
			SMHalfedge h_sm(sm_info.primitive_index);

			pg_cell_to_sm_faces[pg.cell(h_pg)].insert(face(h_sm, sm));
		}
	};

}	// namespace MCAL

#endif