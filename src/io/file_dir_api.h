// $License: GPL-3.0
// $Author: Liu Zhiqiang(from Beihang University)
// $Email: lzhq0930@gmail.com
// 
// Basic file and directory I/O interface.
// 


#ifndef MCAL_FILE_DIRECTORY_API_H
#define MCAL_FILE_DIRECTORY_API_H

#include <fstream>
#include <iostream>
#include <string>
#include <stack>
#include <io.h>
#include <direct.h>
#include <assert.h>

namespace MCAL {   // Mesh Cut Algorithm Libaray, ��Ŷ����������и��㷨�����д���
namespace IO {     // I/O

// �ж��ļ���Ŀ¼�Ƿ����
bool is_path_exist(std::string path)
{
	return _access(path.c_str(), 0) == 0;
}

// ����Ŀ¼(�ɴ����༶)
bool make_directory(std::string path)
{
	std::stack<std::string> dirs;  // Ҫ�����Ķ༶Ŀ¼
	dirs.push(path);
	std::string str = path;

	// �ָ��ַ���path
	while (str.rfind('\\') != str.npos || str.rfind('/') != str.npos)
	{
		std::size_t pos = str.rfind('\\');
		if (pos != str.npos)
		{
			str = str.substr(0, pos);
			dirs.push(str);
		}
		else
		{
			assert(str.rfind('/') != str.npos);
			pos = str.rfind('/');
			str = str.substr(0, pos);
			dirs.push(str);
		}
	}

	dirs.pop();  // ��Ŀ¼(E:)����Ҫ����, ��pop
	bool is_sucess = true;
	while (!dirs.empty())
	{
		str = dirs.top();
		dirs.pop();
		if (!is_path_exist(str) && (_mkdir(str.c_str()) != 0))
			is_sucess = false;
	}
	return is_sucess;
}

// �����ļ�, ע��·�������ļ�·��, ������Ŀ¼.
void copy_file(std::string src_fpath, std::string dst_fpath)
{
	std::ifstream src(src_fpath, std::ios::binary);
	std::ofstream dst(dst_fpath, std::ios::binary);
	dst << src.rdbuf();
}

}	// namespace MCAL::IO
}	// namespace MCAL


#endif