// $Intro: 多面体网格的I/O功能.
// $License: GPL-3.0-or-later
// 
// 
// $Author: Liu Zhiqiang



#ifndef MCAL_POLYHEDRON_GRID_IO_H
#define MCAL_POLYHEDRON_GRID_IO_H

#include <string>

#include "truncated_grid_io.h"


namespace MCAL {   // Mesh Cut Algorithm Libaray, 存放多面体网格切割算法的所有代码
	namespace IO {     // I/O


		// helper function, 拷贝一些必须的文件, 创建目录.
		void initialize_io(std::string input_path, std::string out_path)
		{
			assert(is_path_exist(input_path));

			// 创建必要的目录.
			if (!is_path_exist(out_path))
				make_directory(out_path);

			std::string path = out_path + "\\IN\\geom";
			if (!is_path_exist(path))
				make_directory(path);

			path = out_path + "\\IN\\face";
			if (!is_path_exist(path))
				make_directory(path);

			path = out_path + "\\IN\\cell";
			if (!is_path_exist(path))
				make_directory(path);

			// local_coordinate.txt是显示必需的文件, 拷贝即可.
			copy_file(input_path + "\\local_coordinate.txt", out_path + "\\IN\\local_coordinate.txt");

			path = out_path + "\\OUT\\geom";
			if (!is_path_exist(path))
				make_directory(path);

			path = out_path + "\\OUT\\face";
			if (!is_path_exist(path))
				make_directory(path);

			path = out_path + "\\OUT\\cell";
			if (!is_path_exist(path))
				make_directory(path);

			// local_coordinate.txt是显示必需的文件, 拷贝即可.
			copy_file(input_path + "\\local_coordinate.txt", out_path + "\\OUT\\local_coordinate.txt");
		}

		// PolyhedronGrid的Input函数, 读取dat文件, 建立多面体网格.
		// root_dir: 截断网格的根目录.
		// 
		template <typename PolyhedronGrid>
		bool read_polyhedron_grid(std::string root_dir, PolyhedronGrid& pg)
		{
			typedef typename PolyhedronGrid::Point       P;

			TruncatedGrid<P> tg;
			tg.set_root_path(root_dir);
			if (!tg.load_geom_data() || !tg.load_ijk_info())
			{
				std::cerr << "截断网格读取失败!" << std::endl;
				return false;
			}

			std::cout << "单元包含面的平均数" << (double)(tg.cell_faces_indices.size() / tg.cell_faces_indices_offset.size());

			if (!pg.load_from_truncated_grid(tg))
			{
				std::cerr << "截断网格转换为多面体网格失败!" << std::endl;
				return false;
			}

			return true;
		}

		// PolyhedronGrid的Output函数, 将多面体网格的数据按截断网格的文件格式写出.
		// root_dir 截断网格的根目录
		//
		template <typename PolyhedronGrid>
		bool write_polyhedron_grid(std::string out_path, PolyhedronGrid& pg, bool need_out)
		{
			typedef typename PolyhedronGrid::Point       P;

			TruncatedGrid<P> tg;

			if (!pg.transform_to_truncated_grid(tg, need_out))
			{
				std::cerr << "多面体网格转换为截断网格失败!" << std::endl;
				return false;
			}

			if (!tg.unload_geom_data(out_path) || !tg.unload_ijk_info(out_path) || !tg.unload_planar_info(out_path))
			{
				std::cerr << "截断网格写出失败!" << std::endl;
				return false;
			}

			return true;
		}


	}	// namespace MCAL::IO
}	// namespace MCAL


#endif