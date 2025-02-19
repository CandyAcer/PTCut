// $License: GPL-3.0
// $Author: Liu Zhiqiang(from Beihang University)
// $Email: lzhq0930@gmail.com
// 


#ifndef MACL_EMBEDDED_BODY_IO_H
#define MACL_EMBEDDED_BODY_IO_H

#include <CGAL/Kernel_traits.h>
#include <CGAL/Kernel/global_functions_3.h>

#include <boost/graph/graph_traits.hpp>

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>

namespace MCAL {   // Mesh Cut Algorithm Libaray, ��Ŷ����������и��㷨�����д���
namespace IO {     // I/O

/*
* �и��㷨��������һ�������������һ���и���.
* ������������Ҫ�Լ�����, ��ʵ�ֵ�PolyhedralMesh,
* �и���(Embeddedbody)��������һ����������, ʹ����������, ˼·����ͬ�ض�����.
* ���ڰ�ߵ�����������Ŀǰ������ʵ�ַ�ʽ,
* ���ǽ��и���תΪCGAL::Surface_mesh, ����ʹ��CGAL���еĴ���, ��֤�ɿ���.
*
* ע�⣺���ǽ�ʵ�ֶ��Ľӿ�, д��ֱ�ӵ������еĽӿ�, д����off�ļ�����.
*/

struct Point3d
{
public:
	Point3d(double xcoord = 0, double ycoord = 0, double zcoord = 0)
		:x(xcoord), y(ycoord), z(zcoord)
	{}

	Point3d& operator=(const Point3d& rhs)
	{
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;
		return *this;
	}

public:
	double x;
	double y;
	double z;

};

struct Normal3d
{
public:
	Normal3d(float normx = .0f, float normy = .0f, float normz = .0f)
		:x(normx), y(normy), z(normz)
	{}

public:
	float x;
	float y;
	float z;
};


void construct_origin_and_matrix(const std::string& path,
								 std::array<double, 3>& origin,
								 std::array<std::array<double, 3>, 3>& matrix)
{
	std::ifstream input(path);

	{
		std::string line;
		std::getline(input, line);
		std::stringstream ss(line);
		std::vector<std::string> coords;
		std::string str;
		while (ss >> str)
			coords.push_back(str);

		for (int i = 1; i < 4; ++i)
		{
			std::stringstream get_coord(coords[i]);
			get_coord >> origin[i - 1];
		}
	}

	{
		std::string line;
		std::getline(input, line);
		std::stringstream ss(line);
		std::vector<std::string> coords;
		std::string str;
		while (ss >> str)
			coords.push_back(str);

		for (int i = 1; i < 4; ++i)
		{
			std::stringstream get_coord(coords[i]);
			get_coord >> matrix[0][i - 1];
		}
	}

	{
		std::string line;
		std::getline(input, line);
		std::stringstream ss(line);
		std::vector<std::string> coords;
		std::string str;
		while (ss >> str)
			coords.push_back(str);

		for (int i = 1; i < 4; ++i)
		{
			std::stringstream get_coord(coords[i]);
			get_coord >> matrix[1][i - 1];
		}
	}

	{
		std::string line;
		std::getline(input, line);
		std::stringstream ss(line);
		std::vector<std::string> coords;
		std::string str;
		while (ss >> str)
			coords.push_back(str);

		for (int i = 1; i < 4; ++i)
		{
			std::stringstream get_coord(coords[i]);
			get_coord >> matrix[2][i - 1];
		}
	}
}

// face_grid.dat���ļ���ʽ: 
// INT_MAX(4Byte)+�ļ��汾��(4�ֽ�����)+�и�����Ŀ(4�ֽ�����).
// ����ÿ���и���, �������֣�
//     1) �����Ŀ(4�ֽ�����)+�������(��¼xyz, ע����double����).
//     2) ���������Ŀ(4�ֽ�����)+ÿ��������ķ�����(��¼xyz, ע����float����).
//
// ע�⣺��¼�ĵ�������������ÿ�����������, point_num/tri_num == 3
//

// ���и���dat�ļ�, ��������������
// @param file_path : �и��������ȫ�������face_grid.dat��, ���ļ���·��
// @param sm_set : ת�����ݵ���������, �����ж��
template <typename SurfaceMesh>
bool load_embedded_body(const std::string& file_path,
						std::vector<SurfaceMesh>& sm_set,
						const std::string& local_coord_path)
{
	// ��local_coordinate.txt�ļ�����ȡorigin����ת����
	std::array<double, 3> origin;
	std::array<std::array<double, 3>, 3> matrix;
	construct_origin_and_matrix(local_coord_path, origin, matrix);

	typedef typename SurfaceMesh::Point         Point;
	typedef boost::graph_traits<SurfaceMesh>    GT;
	typedef typename GT::vertex_descriptor      SMVertex;

	FILE* fp;
	if (fopen_s(&fp, file_path.c_str(), "rb") == 0)
	{
		int first_num(0);   // ���а汾���Ƶ��ļ��е�һ���������INT_MAX
		int version(0);     // �ļ��汾��, ���ڰ汾����
		int body_num(0);    // �и�������, face_grid.dat��¼���Ǻܶ���и������һ��

		fread(&first_num, sizeof(int), 1, fp);
		if (first_num != INT_MAX)   // �ļ�û�а汾����, ��������
			return false;

		fread(&version, sizeof(int), 1, fp);
		fread(&body_num, sizeof(int), 1, fp);

		for (int i = 0; i < body_num; ++i)
		{
			int point_num(0);
			fread(&point_num, sizeof(int), 1, fp);

			std::vector<Point3d> points(point_num);  // �и����Լ�����ĵ�����
			std::vector<Point> cgal_points;   // CGAL�ĵ�����
			fread(&points[0], sizeof(Point3d), point_num, fp);

			// ���������ǲ�����Ҫ, ��ֻ��Ϊ�������ⲿ������
			int tri_num(0);
			fread(&tri_num, sizeof(int), 1, fp);

			std::vector<Normal3d> normals(tri_num);
			fread(&normals[0], sizeof(Normal3d), tri_num, fp);

			SurfaceMesh sm;
			std::map<Point, SMVertex> record;  // ȥ��
			for (int i = 0; i < point_num; ++i)
			{
				Point3d p = points[i];
				p.x -= origin[0];
				p.y -= origin[1];
				p.z -= origin[2];
				double t[3];
				for (int k = 0; k < 3; ++k)
					t[k] = matrix[k][0] * p.x + matrix[k][1] * p.y + matrix[k][2] * p.z;

				Point cgal_point(t[0], t[1], t[2]);
				cgal_points.push_back(cgal_point);
				if (!record.count(cgal_point))
				{
					SMVertex v = sm.add_vertex(cgal_point);
					record.insert(std::make_pair(cgal_point, v));
				}
			}

			for (int i = 0; i < tri_num; ++i)
			{
				std::vector<SMVertex> vertices;
				for (int k = 0; k < 3; ++k)
				{
					vertices.push_back(record[cgal_points[3 * i + k]]);
				}
				sm.add_face(vertices);
			}
			sm_set.push_back(sm);
		}
		fclose(fp);
		return true;
	}
	else
		return false;
}


template <typename SurfaceMesh>
bool unload_embedded_body(SurfaceMesh& sm, const std::string& file_path, double xm, double ym, double zm)
{
	typedef typename SurfaceMesh::Point         Point;
	typedef boost::graph_traits<SurfaceMesh>    GT;
	typedef typename GT::vertex_descriptor      SMVertex;
	typedef typename GT::face_descriptor        SMFace;
	typedef typename GT::halfedge_descriptor    SMHalfedge;

	typedef typename CGAL::Kernel_traits<Point>::Kernel      Kernel;
	typedef typename Kernel::Vector_3                        Normal;

	int nBody = 1;

	std::vector<Point3d> points;
	std::vector<Normal3d>  fnorms;
	for (SMFace f : sm.faces())
	{
		SMHalfedge h = sm.halfedge(f);
		Point p0 = sm.point(sm.target(h));
		Point p1 = sm.point(sm.target(sm.next(h)));
		Point p2 = sm.point(sm.source(h));

		points.push_back(Point3d(p0.x() - xm, p0.y() - ym, p0.z() - zm));
		points.push_back(Point3d(p1.x() - xm, p1.y() - ym, p1.z() - zm));
		points.push_back(Point3d(p2.x() - xm, p2.y() - ym, p2.z() - zm));

		Normal f_norm = Kernel::Construct_normal_3()(p0, p1, p2);
		f_norm /= CGAL::sqrt(f_norm.squared_length()); // ��λ��.
		fnorms.push_back(Normal3d(f_norm.x(), f_norm.y(), f_norm.z()));
	}

	FILE* fp;
	if (fopen_s(&fp, file_path.c_str(), "wb") == 0)
	{
		int nMagicNum = INT_MAX;
		fwrite(&nMagicNum, sizeof(int), 1, fp);
		int version = 1;
		fwrite(&version, sizeof(int), 1, fp);

		fwrite(&nBody, sizeof(int), 1, fp);

		{
			int nTri = fnorms.size();
			int nPt = points.size();

			fwrite(&nPt, sizeof(int), 1, fp);
			fwrite(&points[0], sizeof(Point3d), nPt, fp);

			fwrite(&nTri, sizeof(int), 1, fp);
			fwrite(&fnorms[0], sizeof(Normal3d), nTri, fp);
		}

		fflush(fp);
		fclose(fp);

		return true;
	}
	return false;
}

}	// namespace MCAL::IO
}	// namespace MCAL

#endif