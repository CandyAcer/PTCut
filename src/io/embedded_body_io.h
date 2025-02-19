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

namespace MCAL {   // Mesh Cut Algorithm Libaray, 存放多面体网格切割算法的所有代码
namespace IO {     // I/O

/*
* 切割算法的输入是一个多面体网格和一个切割体.
* 多面体网格需要自己定义, 即实现的PolyhedralMesh,
* 切割体(Embeddedbody)本质上是一个表面网格, 使用数组描述, 思路基本同截断网格.
* 基于半边的三角网格是目前最成熟的实现方式,
* 我们将切割体转为CGAL::Surface_mesh, 可以使用CGAL已有的代码, 保证可靠性.
*
* 注意：我们仅实现读的接口, 写出直接调用已有的接口, 写出到off文件即可.
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

// face_grid.dat的文件格式: 
// INT_MAX(4Byte)+文件版本号(4字节整数)+切割体数目(4字节整数).
// 对于每个切割体, 分两部分：
//     1) 点的数目(4字节整数)+点的坐标(记录xyz, 注意是double类型).
//     2) 三角面的数目(4字节整数)+每个三角面的法向量(记录xyz, 注意是float类型).
//
// 注意：记录的点是三角网格中每个面的三个点, point_num/tri_num == 3
//

// 读切割体dat文件, 并建立表面网格
// @param file_path : 切割体的数据全部存放在face_grid.dat中, 该文件的路径
// @param sm_set : 转存数据的三角网格, 可能有多个
template <typename SurfaceMesh>
bool load_embedded_body(const std::string& file_path,
						std::vector<SurfaceMesh>& sm_set,
						const std::string& local_coord_path)
{
	// 读local_coordinate.txt文件，获取origin和旋转矩阵
	std::array<double, 3> origin;
	std::array<std::array<double, 3>, 3> matrix;
	construct_origin_and_matrix(local_coord_path, origin, matrix);

	typedef typename SurfaceMesh::Point         Point;
	typedef boost::graph_traits<SurfaceMesh>    GT;
	typedef typename GT::vertex_descriptor      SMVertex;

	FILE* fp;
	if (fopen_s(&fp, file_path.c_str(), "rb") == 0)
	{
		int first_num(0);   // 在有版本控制的文件中第一个输出的是INT_MAX
		int version(0);     // 文件版本号, 用于版本控制
		int body_num(0);    // 切割体数量, face_grid.dat记录的是很多个切割体而非一个

		fread(&first_num, sizeof(int), 1, fp);
		if (first_num != INT_MAX)   // 文件没有版本控制, 不作处理
			return false;

		fread(&version, sizeof(int), 1, fp);
		fread(&body_num, sizeof(int), 1, fp);

		for (int i = 0; i < body_num; ++i)
		{
			int point_num(0);
			fread(&point_num, sizeof(int), 1, fp);

			std::vector<Point3d> points(point_num);  // 切割体自己定义的点类型
			std::vector<Point> cgal_points;   // CGAL的点类型
			fread(&points[0], sizeof(Point3d), point_num, fp);

			// 法向量我们并不需要, 读只是为了跳过这部分数据
			int tri_num(0);
			fread(&tri_num, sizeof(int), 1, fp);

			std::vector<Normal3d> normals(tri_num);
			fread(&normals[0], sizeof(Normal3d), tri_num, fp);

			SurfaceMesh sm;
			std::map<Point, SMVertex> record;  // 去重
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
		f_norm /= CGAL::sqrt(f_norm.squared_length()); // 单位化.
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