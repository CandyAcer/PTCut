// $License: GPL-3.0
// $Author: Liu Zhiqiang(from Beihang University)
// $Email: lzhq0930@gmail.com(liuzhiqiang@buaa.edu.cn)
// 


#ifndef MCAL_TRUNCATED_GRID_IO_H
#define MCAL_TRUNCATED_GRID_IO_H

#include <CGAL/Kernel_traits.h>

#include <vector>
#include <string>

#include "dat_file_io.h"

namespace MCAL {   // Mesh Cut Algorithm Libaray, ��Ŷ����������и��㷨�����д���
namespace IO {     // I/O

// helper function
template <typename T>
void clear_vec(std::vector<T>& v)
{
	std::vector<T>().swap(v);
}

/*
* �ض�����Ҳ��һ�ֶ���������.
* 
* @param P���������, ������double��ʾ�ĵ�, Ҳ������CGAL��Exact_kernel�ĵ�
*/
template <typename P>
class TruncatedGrid
{
	// typedefs
public:
	typedef P												Point;
	typedef typename CGAL::Kernel_traits<Point>::Kernel		Kernel;
	typedef typename Kernel::Vector_3						Normal;
	typedef int												index_type;
	typedef int												offset_type;

	// member functions
public:
	TruncatedGrid() = default;
	~TruncatedGrid() { release_geom_data(); }

	void set_root_path(std::string path) { root_path = path; }

	bool exist_geom_data() const; // �ض������Ƿ���ڼ�������
	bool load_geom_data();        // ��ȡdat�ļ�, �������񼸺�����
	bool unload_geom_data(const std::string out_path);	// ������ļ�������д��dat�ļ�
	void release_geom_data();     // �ͷŽض����񼸺�������ռ�õ��ڴ�

	bool exist_ijk_info() const;
	bool load_ijk_info();
	bool unload_ijk_info(const std::string out_path);
	void release_ijk_info();

	bool unload_planar_info(const std::string out_path);

	// data members
public:
	std::string root_path;  // �ض����������ļ��ĸ�Ŀ¼

	// geom data
	std::vector<Point> vertices;     // �ض���������ж���
	std::vector<index_type> face_vertices_indices;           // ������Ķ�������
	std::vector<offset_type> face_vertices_indices_offset;   // �涥����������ʼλ��
	std::vector<index_type> cell_faces_indices;              // ��Ԫ������������
	std::vector<offset_type> cell_faces_indices_offset;      // ��Ԫ����������ʼλ��
	std::vector<int> faces_directions_per_cell;             // ��Ԫ������������ʹ�ñ�ʶ
	std::vector<Normal> normal_per_face;                     // ÿ����Ԫ��ķ�����
	// ����2��������Ҫ���н���
	std::vector<index_type> face_edges_indices;              // ������ı�����
	std::vector<offset_type> face_edges_indices_offset;      // �����������ʼλ��

	std::vector<index_type> iidx_per_cell;
	std::vector<index_type> jidx_per_cell;
	std::vector<index_type> kidx_per_cell;

	std::vector<index_type> startend_per_planarcell;
};

template <typename P>
bool
TruncatedGrid<P>::
exist_geom_data() const
{
	// ����vectorΪ��, �򲻴��ڼ�������
	if (vertices.empty() || \
		face_vertices_indices.empty() || \
		face_vertices_indices_offset.empty() || \
		cell_faces_indices.empty() || \
		cell_faces_indices_offset.empty())
	{
		return false;
	}
	else
		return true;
}

template<typename P>
bool
TruncatedGrid<P>::
load_geom_data()
{
	if (root_path.empty())
		return false;

	// ��ȡ��������ǰ�����ԭ������.
	release_geom_data();

	// 1.��ȡvertices.dat, ��ȡ�ض���������е�
	std::string path = root_path + "\\geom\\vertices.dat";
	std::vector<double> vertex_coords;
	if (!read_dat_file(path, vertex_coords))
		return false;

	vertices.resize(vertex_coords.size() / 3);
	for (std::size_t i = 0; i < vertex_coords.size() / 3; ++i)
	{
		double x = vertex_coords[3 * i];
		double y = vertex_coords[3 * i + 1];
		double z = vertex_coords[3 * i + 2];
		vertices[i] = Point(x, y, z);
	}

	// 2.��ȡface_vertices_indices.dat
	path = root_path + "\\geom\\face_vertices_indices.dat";
	if (!read_dat_file(path, face_vertices_indices))
		return false;

	// 3.��ȡface_vertices_indices_offset.dat
	path = root_path + "\\geom\\face_vertices_indices_offset.dat";
	if (!read_dat_file(path, face_vertices_indices_offset))
		return false;

	// 4.��ȡcell_faces_indices.dat
	path = root_path + "\\geom\\cell_faces_indices.dat";
	if (!read_dat_file(path, cell_faces_indices))
		return false;

	// 5.��ȡcell_faces_indices_offset.dat
	path = root_path + "\\geom\\cell_faces_indices_offset.dat";
	if (!read_dat_file(path, cell_faces_indices_offset))
		return false;

	// 6.��ȡfaces_directions_per_cell.dat
	path = root_path + "\\geom\\faces_directions_per_cell.dat";
	if (!read_dat_file(path, faces_directions_per_cell))
		return false;

	// 7.��ȡnormal_per_face.dat
	path = root_path + "\\face\\normal_per_face.dat";
	std::vector<float> normal_coords;
	if (!read_dat_file(path, normal_coords))
		return false;

	normal_per_face.resize(normal_coords.size() / 3);
	for (std::size_t i = 0; i < normal_coords.size() / 3; ++i)
	{
		float x = normal_coords[3 * i];
		float y = normal_coords[3 * i + 1];
		float z = normal_coords[3 * i + 2];
		normal_per_face[i] = Normal(x, y, z);
	}

	return true;
}

template<typename P>
bool
TruncatedGrid<P>::
unload_geom_data(std::string out_path)
{
	std::string fpath;

	// 1.д��vertices.dat
	std::vector<double> vertex_coords;
	if (!vertices.empty())
	{
		vertex_coords.resize(3 * vertices.size());
		for (std::size_t i = 0; i < vertices.size(); ++i)
		{
			vertex_coords[3 * i] = vertices[i].x();
			vertex_coords[3 * i + 1] = vertices[i].y();
			vertex_coords[3 * i + 2] = vertices[i].z();
		}
	}

	fpath = out_path + "\\geom\\vertices.dat";
	if (!write_dat_file(fpath, vertex_coords, "m", 3, false))
		return false;

	// 2.д��face_vertices_indices.dat
	fpath = out_path + "\\geom\\face_vertices_indices.dat";
	if (!write_dat_file(fpath, face_vertices_indices, "int", 1, false))
		return false;

	// 3.д��face_vertices_indices_offset.dat
	fpath = out_path + "\\geom\\face_vertices_indices_offset.dat";
	if (!write_dat_file(fpath, face_vertices_indices_offset, "int", 1, false))
		return false;

	// 4.д��cell_faces_indices.dat
	fpath = out_path + "\\geom\\cell_faces_indices.dat";
	if (!write_dat_file(fpath, cell_faces_indices, "int", 1, false))
		return false;

	// 5.д��cell_faces_indices_offset.dat
	fpath = out_path + "\\geom\\cell_faces_indices_offset.dat";
	if (!write_dat_file(fpath, cell_faces_indices_offset, "int", 1, false))
		return false;

	// 6.д��faces_directions_per_cell.dat
	fpath = out_path + "\\geom\\faces_directions_per_cell.dat";
	if (!write_dat_file(fpath, faces_directions_per_cell, "int", 1, false))
		return false;

	// 7.д��normal_per_face.dat
	std::vector<float> normal_coords;
	if (!normal_per_face.empty())
	{
		normal_coords.resize(3 * normal_per_face.size());
		for (std::size_t i = 0; i < normal_per_face.size(); ++i)
		{
			normal_coords[3 * i] = normal_per_face[i].x();
			normal_coords[3 * i + 1] = normal_per_face[i].y();
			normal_coords[3 * i + 2] = normal_per_face[i].z();
		}
	}

	fpath = out_path + "\\face\\normal_per_face.dat";
	if (!write_dat_file(fpath, normal_coords, "m", 3, false))
		return false;

	return true;
}

template <typename P>
void TruncatedGrid<P>::release_geom_data()
{
	clear_vec(vertices);
	clear_vec(face_vertices_indices);
	clear_vec(face_vertices_indices_offset);
	clear_vec(cell_faces_indices);
	clear_vec(cell_faces_indices_offset);
	clear_vec(faces_directions_per_cell);
	clear_vec(normal_per_face);

	clear_vec(face_edges_indices);
	clear_vec(face_edges_indices_offset);
}

template<typename P>
bool
TruncatedGrid<P>::
exist_ijk_info() const
{
	if (iidx_per_cell.empty() || jidx_per_cell.empty() || kidx_per_cell.empty())
		return false;
	else
		return true;
}

template<typename P>
bool
TruncatedGrid<P>::
load_ijk_info()
{
	if (root_path.empty())
		return false;

	release_ijk_info();

	// 1.��ȡI.dat
	std::string path = root_path + "\\cell\\I.dat";
	if (!read_dat_file(path, iidx_per_cell))
		return false;

	// 2.��ȡJ.dat
	path = root_path + "\\cell\\J.dat";
	if (!read_dat_file(path, jidx_per_cell))
		return false;

	// 3.��ȡK.dat
	path = root_path + "\\cell\\K.dat";
	if (!read_dat_file(path, kidx_per_cell))
		return false;

	return true;
}

template<typename P>
bool
TruncatedGrid<P>::
unload_ijk_info(const std::string out_path)
{
	std::string fpath;

	// 1.д��I.dat
	fpath = out_path + "\\cell\\I.dat";
	if (!write_dat_file(fpath, iidx_per_cell, "int", 1, false))
		return false;

	// 2.д��J.dat
	fpath = out_path + "\\cell\\J.dat";
	if (!write_dat_file(fpath, jidx_per_cell, "int", 1, false))
		return false;

	// 3.д��K.dat
	fpath = out_path + "\\cell\\K.dat";
	if (!write_dat_file(fpath, kidx_per_cell, "int", 1, false))
		return false;

	return true;
}


template<typename P>
void
TruncatedGrid<P>::
release_ijk_info()
{
	clear_vec(iidx_per_cell);
	clear_vec(jidx_per_cell);
	clear_vec(kidx_per_cell);
}

template<typename P>
bool
TruncatedGrid<P>::
unload_planar_info(const std::string out_path)
{
	std::string fpath = out_path + "\\cell\\layercolumngrid_startend_per_planarcell.dat";
	if (!write_dat_file(fpath, startend_per_planarcell, "int", 2, false))
		return false;
	return true;
}

}	// namespace MCAL::IO
}	// namespace MCAL

#endif