// $Intro: �����������I/O����.
// $License: GPL-3.0-or-later
// 
// 
// $Author: Liu Zhiqiang



#ifndef MCAL_POLYHEDRON_GRID_IO_H
#define MCAL_POLYHEDRON_GRID_IO_H

#include <string>

#include "truncated_grid_io.h"


namespace MCAL {   // Mesh Cut Algorithm Libaray, ��Ŷ����������и��㷨�����д���
	namespace IO {     // I/O


		// helper function, ����һЩ������ļ�, ����Ŀ¼.
		void initialize_io(std::string input_path, std::string out_path)
		{
			assert(is_path_exist(input_path));

			// ������Ҫ��Ŀ¼.
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

			// local_coordinate.txt����ʾ������ļ�, ��������.
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

			// local_coordinate.txt����ʾ������ļ�, ��������.
			copy_file(input_path + "\\local_coordinate.txt", out_path + "\\OUT\\local_coordinate.txt");
		}

		// PolyhedronGrid��Input����, ��ȡdat�ļ�, ��������������.
		// root_dir: �ض�����ĸ�Ŀ¼.
		// 
		template <typename PolyhedronGrid>
		bool read_polyhedron_grid(std::string root_dir, PolyhedronGrid& pg)
		{
			typedef typename PolyhedronGrid::Point       P;

			TruncatedGrid<P> tg;
			tg.set_root_path(root_dir);
			if (!tg.load_geom_data() || !tg.load_ijk_info())
			{
				std::cerr << "�ض������ȡʧ��!" << std::endl;
				return false;
			}

			std::cout << "��Ԫ�������ƽ����" << (double)(tg.cell_faces_indices.size() / tg.cell_faces_indices_offset.size());

			if (!pg.load_from_truncated_grid(tg))
			{
				std::cerr << "�ض�����ת��Ϊ����������ʧ��!" << std::endl;
				return false;
			}

			return true;
		}

		// PolyhedronGrid��Output����, ����������������ݰ��ض�������ļ���ʽд��.
		// root_dir �ض�����ĸ�Ŀ¼
		//
		template <typename PolyhedronGrid>
		bool write_polyhedron_grid(std::string out_path, PolyhedronGrid& pg, bool need_out)
		{
			typedef typename PolyhedronGrid::Point       P;

			TruncatedGrid<P> tg;

			if (!pg.transform_to_truncated_grid(tg, need_out))
			{
				std::cerr << "����������ת��Ϊ�ض�����ʧ��!" << std::endl;
				return false;
			}

			if (!tg.unload_geom_data(out_path) || !tg.unload_ijk_info(out_path) || !tg.unload_planar_info(out_path))
			{
				std::cerr << "�ض�����д��ʧ��!" << std::endl;
				return false;
			}

			return true;
		}


	}	// namespace MCAL::IO
}	// namespace MCAL


#endif