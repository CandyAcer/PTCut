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

namespace MCAL {   // Mesh Cut Algorithm Libaray, 存放多面体网格切割算法的所有代码
namespace IO {     // I/O

// 判断文件或目录是否存在
bool is_path_exist(std::string path)
{
	return _access(path.c_str(), 0) == 0;
}

// 创建目录(可创建多级)
bool make_directory(std::string path)
{
	std::stack<std::string> dirs;  // 要创建的多级目录
	dirs.push(path);
	std::string str = path;

	// 分割字符串path
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

	dirs.pop();  // 根目录(E:)不需要创建, 故pop
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

// 复制文件, 注意路径都是文件路径, 不能是目录.
void copy_file(std::string src_fpath, std::string dst_fpath)
{
	std::ifstream src(src_fpath, std::ios::binary);
	std::ofstream dst(dst_fpath, std::ios::binary);
	dst << src.rdbuf();
}

}	// namespace MCAL::IO
}	// namespace MCAL


#endif