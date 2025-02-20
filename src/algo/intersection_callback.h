// $License: GPL-3.0
// $Author: Liu Zhiqiang(from Beihang University)
// $Email: lzhq0930@gmail.com(liuzhiqiang@buaa.edu.cn)
// 


#ifndef MCAL_ALGO_INTERSECTION_CALLBACK_H
#define MCAL_ALGO_INTERSECTION_CALLBACK_H

#include <CGAL/Box_intersection_d/Box_with_info_d.h>
#include <boost/graph/graph_traits.hpp>

namespace MCAL   // Mesh Cut Algorithm Libaray, ��Ŷ����������и��㷨��ʵ��
{

// helper class
// CGAL��box_intersection_d�㷨ʹ�õ��Ǵ�info��AABB��Χ��, �и��㷨��info������İ��.
// ��������, SurfaceMesh��PolyhedralMesh�İ���ǲ�һ����.
// ����ʹ��BoxInfo�Զ��߽���ͳһ, �Ӷ����̶ȵظ���CGAL::box_intersection_d�㷨.
//
template <typename SurfaceMesh, typename PolyhedralMesh>
struct BoxInfo
{
	typedef typename SurfaceMesh::size_type  size_type;
	static_assert(std::is_same<size_type, typename PolyhedralMesh::size_type>::value,
				  "��������������������size_type���Ͳ�һ�£�");

	typedef typename PolyhedralMesh::Halfedge PMHalfedge;

	size_type primitive_index;
	size_type array_offset;

	BoxInfo(size_type idx = std::numeric_limits<size_type>::max(),
		size_type off = std::numeric_limits<size_type>::max())
		:primitive_index(idx), array_offset(off)
	{}

	BoxInfo(const PMHalfedge& h)
		:primitive_index(h.he_cell), array_offset(h.he_off)
	{}

	BoxInfo(const BoxInfo& rhs) { operator=(rhs); }

	BoxInfo& operator=(const BoxInfo& rhs)
	{
		primitive_index = rhs.primitive_index;
		array_offset = rhs.array_offset;
		return *this;
	}
};

// ��⵽SurfaceMesh��edge��PolyhedralMesh��face��Χ���ཻʱ, ���ô�callback.
template <typename SurfaceMesh,
		  typename PolyhedralMesh,
		  typename EdgeToFaces>
class CollectSMEdgeToPMFaces
{
	// typedefs
	typedef boost::graph_traits<SurfaceMesh>             GT;
	typedef typename GT::halfedge_descriptor             SMHalfedge;

	typedef typename PolyhedralMesh::FaceIndex           PMFace;
	typedef typename PolyhedralMesh::CellIndex           PMCell;
	typedef typename PolyhedralMesh::Halfedge            PMHalfedge;

	typedef BoxInfo<SurfaceMesh, PolyhedralMesh>                       Info;
	typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, Info>  Box;

	// data members
	SurfaceMesh& sm;
	PolyhedralMesh& pm;
	EdgeToFaces& sm_edge_to_pm_faces;

public:
	CollectSMEdgeToPMFaces(SurfaceMesh& sm_,
						   PolyhedralMesh& pm_,
						   EdgeToFaces& edge_to_faces_)
		:sm(sm_), pm(pm_), sm_edge_to_pm_faces(edge_to_faces_)
	{}

	void operator()(const Box& sm_edge_box, const Box& pm_face_box)
	{
		Info sm_info = sm_edge_box.info();
		SMHalfedge h_sm(sm_info.primitive_index);

		Info pm_info = pm_face_box.info();
		PMHalfedge h_pm(PMCell(pm_info.primitive_index), pm_info.array_offset);
		PMFace f_pm = pm.face(h_pm);

		sm_edge_to_pm_faces[edge(h_sm, sm)].insert(f_pm);
	}
};

// ��⵽PolyhedralMesh��cell��SurfaceMesh��face��Χ���ཻʱ, ���ô�callback.
template <typename PolyhedralMesh,
		  typename SurfaceMesh,
		  typename CellToFaces>
class CollectPMCellToSMFaces
{
private:
	PolyhedralMesh& pm;
	SurfaceMesh& sm;
	CellToFaces& pm_cell_to_sm_faces;

	typedef boost::graph_traits<SurfaceMesh>             GT;
	typedef typename GT::halfedge_descriptor             SMHalfedge;

	typedef typename PolyhedralMesh::Halfedge            PMHalfedge;
	typedef typename PolyhedralMesh::CellIndex           PMCell;

	typedef BoxInfo<SurfaceMesh, PolyhedralMesh>         Info;
	typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, Info> Box;

public:
	CollectPMCellToSMFaces(PolyhedralMesh& pm_,
						   SurfaceMesh& sm_,
						   CellToFaces& cell_to_faces)
		:pm(pm_), sm(sm_), pm_cell_to_sm_faces(cell_to_faces)
	{}

	void operator()(const Box& pm_cell_box, const Box& sm_face_box)
	{
		Info pm_info = pm_cell_box.info();
		PMHalfedge h_pm(PMCell(pm_info.primitive_index), pm_info.array_offset);

		Info sm_info = sm_face_box.info();
		SMHalfedge h_sm(sm_info.primitive_index);

		pm_cell_to_sm_faces[pm.cell(h_pm)].insert(face(h_sm, sm));
	}
};

}	// namespace MCAL

#endif