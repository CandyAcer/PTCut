// $License: GPL-3.0
// $Author: Liu Zhiqiang(from Beihang University)
// $Email: lzhq0930@gmail.com(liuzhiqiang@buaa.edu.cn)
// 


#ifndef MCAL_DAT_FILE_IO_H
#define MCAL_DAT_FILE_IO_H

#include <vector>
#include <string>

#include "file_dir_api.h"

namespace MCAL {   // Mesh Cut Algorithm Libaray, 存放多面体网格切割算法的所有代码
namespace IO {     // I/O

/*
* 截断网格的二进制文件名为 *.dat(例如vertices.dat)
* 一个dat文件对应截断网格的一个数组(例如vertices.dat对应TruncatedGrid的vertices)
*
* dat文件由文件头和二进制数据组成
* 1) 文件头
*    记录文件的元数据, 用于说明文件存储的数据内容
*    长度可变, 但始终是64字节的倍数, 不足的字节会进行填充
*
*
*    文件头中的元数据信息示例(以vertices.dat为例)：
*    offset = 128       // 二进制数据从第128字节开始, 0~127为文件头
*    units = m          // 元素值的计量单位, m代表meter(米), 索引的计量单位是int
*    ncols = 3          // 一个几何元素(点, 面, 单元)对应几个分量
*                       // 此处一个点要依次记录xyz三个值, 所以为3
*
*    dim = 2928         // 数组包含的几何元素个数, 此处说明记录了2928个点
*    type = double      // 数组的数据类型, 此处说明点坐标xyz都是double类型
*    nbytes = 8         // 数组元素的类型占多少字节(double为8字节)
*    big = 0            // true为大端格式, false为小端格式
*
*
* 2) 二进制数据
*    对应某个数组
*
* 截断网格文件的读写由两个类共同完成
*     1) TGFileBasicIO类负责操作dat文件
*     2) GridDataStore类负责将从dat文件读出的数组数据存储在vector中
*
*/


/*
* GridDataStore类负责将从dat文件读出的数组数据存储在vector中
*/
class GridDataStore
{
	friend class TGFileBasicIO;

public:
	GridDataStore()
		: unit(std::string()), col(0), dim(0)
		, type(std::string()), nbyte(0)
		, data_len(0), p_datablock(NULL)
	{}

	~GridDataStore() { release_data(); }

	// 将p_datablock指向的二进制数据存放在数组中
	template <typename T>
	bool get_data(std::vector<T>& data);

	void release_data();

public:
	// 与TGFileBasicIO的数据成员对应, 不再解释含义
	std::string unit;
	int col;
	int dim;
	std::string type;
	int nbyte;
	__int64 data_len;    // 数据的字节数

private:
	// 指向读出数据的指针, 之后做类型转换即可
	// 该指针不允许被外界访问, 获取数据请调用get_data()
	void* p_datablock;
};

template <typename T>
bool
GridDataStore::
get_data(std::vector<T>& data)
{
	if (dim <= 0 || p_datablock == NULL)
		return false;

	// 若获取数据类型与当前数据类型一致, 直接拷贝即可
	if (typeid(T).name() == type)
	{
		std::vector<T>* pvec = (std::vector<T>*)p_datablock;
		if (pvec)
		{
			pvec->swap(data);
			data_len = 0;
		}
	}
	else
	{
		int vec_size = col * dim;
		data.resize(vec_size);

		if (type == "double" || type == "d")
		{
			std::vector<double>& pdata = *(std::vector<double>*)p_datablock;
			for (int i = 0; i < vec_size; i++)
			{
				data[i] = pdata[i];
			}
		}
		else if (type == "int" || type == "i")
		{
			std::vector<int>& pdata = *(std::vector<int>*)p_datablock;
			for (int i = 0; i < vec_size; i++)
			{
				data[i] = pdata[i];
			}
		}
		else
			return false;
	}
	return true;
}

void GridDataStore::release_data()
{
	if (p_datablock != NULL)
	{
		if (type == "double" || type == "d")
		{
			std::vector<double>* pvec = (std::vector<double>*)p_datablock;
			if (pvec)
			{
				std::vector<double> v;
				pvec->swap(v);
				delete pvec;
				pvec = NULL;
				p_datablock = NULL;
			}
		}
		else if (type == "int" || type == "i")
		{
			std::vector<int>* pvec = (std::vector<int>*)p_datablock;
			if (pvec)
			{
				std::vector<int> v;
				pvec->swap(v);
				delete pvec;
				pvec = NULL;
				p_datablock = NULL;
			}
		}
		else
		{
			delete[] p_datablock;
			p_datablock = NULL;
		}
	}
	col = dim = nbyte = data_len = 0;
}

/*
* TGFileBasicIO类负责操作dat文件
*/
class TGFileBasicIO
{
public:
	// 默认构造函数会设置元数据的缺省值
	// 元数据作为一个整体出现, 没有必要单个数据成员的getter和setter
	TGFileBasicIO()
		: offset(64), unit(std::string()), col(0), dim(0)
		, type(std::string()), nbyte(0), big_endian(false) /*默认为小端格式*/
		, fp(NULL), pbuf(NULL)
	{}

	TGFileBasicIO(std::string str_unit, int n_col, int n_dim,
		std::string str_type, int n_byte, bool b_endian, FILE* p);

	~TGFileBasicIO();

	void set_file_path(std::string path) { file_path = path; }

	int estimate_meta_length();  // 计算元数据的字节数
	bool read_meta_data();   // 读取元数据
	bool write_meta_data(bool b);  // 写出元数据

	bool read_array_data(GridDataStore& data);
	bool write_array_data(int n_dim, int n_col, int n_byte, void* pdata, __int64 data_len);

private:
	bool inverse_data_byte_order(void* pdata, int data_unit, int ele_byte);

	template <typename T>
	bool read_big_data_wrapper(FILE* fp,  /*流关联的fp*/
		int data_start_byte,  /*二进制数据开始的字节数*/
		std::vector<T>& data,  /*存储数据的数组*/
		__int64 request_bytes,  /*读取的字节数*/
		bool is_big  /*是否为大端*/);

	void write_big_data_wrapper(const void* pdata, __int64 data_len, FILE* fp_write);

private: // data members
	std::string file_path;   // 读写的文件路径

	// 文件头的元数据
	int offset;              // 二进制数据的起始偏移字节数, 第一个关键字
	std::string unit;        // 元素值对应的计量单位, 如m, ft, int, 第二个关键字
	int col;                 // 一个几何元素(点, 面, 单元)对应几个分量, 
	// 若对应分量数不固定, 则为1; 固定则为固定数值

	int dim;                 // 包含的几何元素(点, 面, 单元)个数
	// 若对应分量数不固定, dim为所有几何元素的分量数之和
	// 若固定, 则dim为几何元素的个数
	// 总之, dim * col = 数组的size

	std::string type;        // 数组元素的数据类型, 如double, int
	int nbyte;               // 数组元素的字节长度, dim * col * nbyte = 二进制数组的字节数
	bool big_endian;         // true为大端格式, false为小端格式

	FILE* fp;     // 文件指针
	char* pbuf;   // 设置文件系统缓冲区
};

TGFileBasicIO::TGFileBasicIO(std::string str_unit,
							 int n_col,
							 int n_dim,
							 std::string str_type,
							 int n_byte,
							 bool b_endian,
							 FILE* p)
	:unit(str_unit), col(n_col), dim(n_dim)
	, type(str_type), nbyte(n_byte), big_endian(b_endian)
{
	offset = estimate_meta_length();

	assert(p != NULL);
	fp = p;

	if (pbuf)
	{
		delete[] pbuf;
		pbuf = NULL;
	}
	int buffer_size = 8 * 1024 * 1024;
	pbuf = new char[buffer_size];

	if (fp && pbuf)
		setvbuf(fp, pbuf, _IOFBF, buffer_size);
}

TGFileBasicIO::~TGFileBasicIO()
{
	if (fp != NULL)
	{
		fflush(fp);
		fclose(fp);
		fp = NULL;
	}

	if (pbuf)
	{
		delete[] pbuf;
		pbuf = NULL;
	}
}

int TGFileBasicIO::estimate_meta_length()
{
	// 元数据的格式为key=value
	// 计算元数据的方法为：key所占的字符数+value所占字节数+1(换行符'\n')
	int len = 8 + sizeof(int) +  /*"offset=?\n", offset字段的长度*/
		7 + 8 +  /*"units=?\n", ?为string类型, 预留8字节*/
		7 + sizeof(int) +  /*"ncols=?\n"*/
		5 + sizeof(int) +  /*"dim=?\n"*/
		6 + 8 +  /*"type=?\n", ?为string类型, 预留8字节*/
		8 + sizeof(int) +  /*"nbytes=?\n", ?代表对应的type类型数据需要的字节数 */
		5 + sizeof(int)  /*"big=0 or 1\n"*/;

	// 64字节对齐
	len = (len / 64 + 1) * 64;
	return len;
}

bool TGFileBasicIO::read_meta_data()
{
	if (fp == NULL)
		return false;

	// 将文件指针移动到文件的开始位置
	_fseeki64(fp, 0, SEEK_SET);

	// offset为64的倍数, 故64字节是读的基本单位
	char buffer[64];
	int buf_size = 64;
	memset(buffer, '\0', buf_size);

	int fp_location = 0;           // 文件指针的位置, 以字节为单位
	while (fp_location < offset)
	{
		// 从fp关联的流中读取一行到buffer中
		fgets(buffer, buf_size, fp);
		// 修改fp_location为文件指针的当前位置
		fp_location = ftell(fp);

		// 分割buffer中的字符串, 分割后的字符串在cutted_strs中
		const char* delimiters = "= \t\n";
		char* uncut_str = NULL;
		char* pstr = strtok_s(buffer, delimiters, &uncut_str);
		std::vector<std::string> cutted_strs;

		while (pstr != NULL)
		{
			cutted_strs.push_back(pstr);
			pstr = strtok_s(NULL, delimiters, &uncut_str);
		}

		if (cutted_strs.empty())
			break;

		if (cutted_strs[0] == "offset")
		{
			offset = atoi(cutted_strs[1].c_str());
			assert((offset > 0) && (offset % 64 == 0));
			continue;
		}

		if (cutted_strs[0] == "dim" && (cutted_strs.size() == 2))
		{
			dim = atoi(cutted_strs[1].c_str());
			assert(dim > 0);
			continue;
		}

		if (cutted_strs[0] == "units" && (cutted_strs.size() == 2))
		{
			unit = cutted_strs[1];
			continue;
		}

		if (cutted_strs[0] == "type" && (cutted_strs.size() >= 2))
		{
			type = "";
			for (int i = 1; i < cutted_strs.size(); ++i)
			{
				type += cutted_strs[i];
				if (i < cutted_strs.size() - 1)
					type += " ";
			}
		}

		if (cutted_strs[0] == "nbytes" && (cutted_strs.size() >= 2))
		{
			nbyte = atoi(cutted_strs[1].c_str());
			assert(nbyte > 0);
			continue;
		}

		if (cutted_strs[0] == "ncols" && (cutted_strs.size() >= 2))
		{
			col = atoi(cutted_strs[1].c_str());
			continue;
		}

		if (cutted_strs[0] == "big" && (cutted_strs.size() >= 2))
		{
			big_endian = atoi(cutted_strs[1].c_str());
			continue;
		}

		if (cutted_strs[0] == "<<")
			break;
	}

	// 校验元数据正确性
	if (offset <= 0 || (offset % 64 != 0) ||
		nbyte <= 0 || type.empty() || col < 1)
		return false;

	__int64 data_len = nbyte * col * dim + offset;

	_fseeki64(fp, 0L, SEEK_END);
	__int64 file_len = _ftelli64(fp);

	if (data_len != file_len)
		return false;

	return true;
}

bool TGFileBasicIO::write_meta_data(bool b)
{
	if (fp == NULL)
		return false;

	_fseeki64(fp, 0L, SEEK_SET);

	// 将元数据全部填充为'\0'
	char* pmeta = new char[offset];
	memset(pmeta, '\0', offset);
	fwrite(pmeta, offset, 1, fp);

	_fseeki64(fp, 0L, SEEK_SET);

	// 填写meta各字段

	// offset
	std::string key = "offset";
	sprintf_s(pmeta, offset, "%s=%d\n", key.c_str(), offset);
	std::size_t str_len = strlen(pmeta);
	fwrite(pmeta, str_len, 1, fp);

	// units
	std::string str_unit = "units=" + unit + "\n";
	str_len = str_unit.length();
	fwrite(str_unit.c_str(), str_len, 1, fp);

	// ncols
	key = "ncols";
	sprintf_s(pmeta, offset, "%s=%d\n", key.c_str(), col);
	str_len = strlen(pmeta);
	fwrite(pmeta, str_len, 1, fp);

	// dim
	key = "dim";
	sprintf_s(pmeta, offset, "%s=%d\n", key.c_str(), dim);
	str_len = strlen(pmeta);
	fwrite(pmeta, str_len, 1, fp);

	// type
	std::string str_type = "type=" + type + "\n";
	str_len = str_type.length();
	fwrite(str_type.c_str(), str_len, 1, fp);

	// nbytes
	key = "nbytes";
	sprintf_s(pmeta, offset, "%s=%d\n", key.c_str(), nbyte);
	str_len = strlen(pmeta);
	fwrite(pmeta, str_len, 1, fp);

	// big
	key = "big";
	sprintf_s(pmeta, offset, "%s=%d\n", key.c_str(), big_endian);
	str_len = strlen(pmeta);
	fwrite(pmeta, str_len, 1, fp);

	// free memory
	delete[] pmeta;
	pmeta = NULL;

	__int64 meta_len = _ftelli64(fp);
	assert(meta_len <= offset);

	if (b)
		_fseeki64(fp, offset, SEEK_SET);

	return meta_len <= offset;
}

bool TGFileBasicIO::read_array_data(GridDataStore& ob_store)
{
	if (file_path.empty())
		return false;

	if (fopen_s(&fp, file_path.c_str(), "rb") != 0)
		return false;

	// 读取元数据
	if (!read_meta_data())
		return false;

	ob_store.unit = unit;
	ob_store.col = col;
	ob_store.dim = dim;
	ob_store.type = type;
	ob_store.nbyte = nbyte;
	ob_store.data_len = (__int64)nbyte * col * dim;
	ob_store.p_datablock = NULL;

	// 文件指针移动到二进制数据的起始位置
	_fseeki64(fp, offset, SEEK_SET);

	if (ob_store.type == "double" || ob_store.type == "d")
	{
		std::vector<double>* pvec = new std::vector<double>;
		if (!read_big_data_wrapper(fp, offset, *pvec, ob_store.data_len, big_endian))
			return false;
		ob_store.p_datablock = pvec;
	}
	else if (ob_store.type == "int" || ob_store.type == "i")
	{
		std::vector<int>* pvec = new std::vector<int>;
		if (!read_big_data_wrapper(fp, offset, *pvec, ob_store.data_len, big_endian))
			return false;
		ob_store.p_datablock = pvec;
	}
	else if (ob_store.type == "float" || ob_store.type == "f")
	{
		std::vector<float>* pvec = new std::vector<float>;
		if (!read_big_data_wrapper(fp, offset, *pvec, ob_store.data_len, big_endian))
			return false;
		ob_store.p_datablock = pvec;
	}

	if (!ob_store.p_datablock)
		return false;

	fclose(fp);
	fp = NULL;

	return true;
}

bool
TGFileBasicIO::
write_array_data(int n_dim, int n_col, int n_byte, void* pdata, __int64 data_len)
{
	// 参数是否有错误？
	__int64 len = (__int64)n_dim * n_col * n_byte;
	if (data_len != len)
		return false;

	if (dim != n_dim || col != n_col || nbyte != n_byte)
		return false;

	_fseeki64(fp, 0, SEEK_END);
	__int64 file_len = _ftelli64(fp);

	if (file_len != offset)
		return false;

	if (big_endian && pdata)
		inverse_data_byte_order(pdata, n_dim * n_col, n_byte);

	if (pdata != NULL)
		write_big_data_wrapper(pdata, data_len, fp);

	return true;
}

static void inverse_byte_order(char* pdata, int count)
{
	assert(count > 0);
	int half = count / 2;
	for (int i = 0; i < half; ++i)
	{
		char c = pdata[i];
		pdata[i] = pdata[count - 1 - i];
		pdata[count - 1 - i] = c;
	}
}

bool
TGFileBasicIO::
inverse_data_byte_order(void* pdata, int data_unit, int ele_byte)
{
	assert(big_endian == true);

	if (pdata == NULL)
		return false;

	if (ele_byte != 8 && ele_byte != 4 && ele_byte != 2)
		return false;

	for (int i = 0; i < data_unit; ++i)
		inverse_byte_order((char*)((char*)pdata + (__int64)i * ele_byte), ele_byte);

	return true;
}


template<typename T>
bool
TGFileBasicIO::
read_big_data_wrapper(FILE* fp,  /*流关联的fp*/
					  int data_start_byte,  /*二进制数据开始的字节数*/
					  std::vector<T>& data,  /*存储数据的数组*/
					  __int64 request_bytes,  /*读取的字节数*/
					  bool is_big  /*是否为大端*/)
{
	if (fp == NULL || request_bytes <= 0)
		return false;

	int type_len = sizeof(T);
	data.resize(request_bytes / type_len);

	char* p = (char*)&data[0];

	_fseeki64(fp, data_start_byte, SEEK_SET);

	// 大文件需要分次读
	__int64 have_read_bytes = 0;                      // 已读取的字节数
	__int64 max_bytes_per_read = 256 * 1024 * 1024;   // 一次最多读取的字节数
	__int64 bytes_per_read = 0;                       // 一次需要读取的字节数

	do
	{
		if (request_bytes - have_read_bytes > max_bytes_per_read)
			bytes_per_read = max_bytes_per_read;
		else
			bytes_per_read = request_bytes - have_read_bytes;

		if (fread((char*)p + have_read_bytes, bytes_per_read, 1, fp) != 1)
			return false;

		have_read_bytes += bytes_per_read;

	} while (request_bytes > have_read_bytes);

	return true;
}


void
TGFileBasicIO::
write_big_data_wrapper(const void* pdata, __int64 data_len, FILE* fp_write)
{
	__int64 have_write_bytes = 0;
	__int64 max_bytes_per_write = 256 * 1024 * 1024;
	__int64 bytes_per_write = 0;

	do
	{
		if (data_len - have_write_bytes > max_bytes_per_write)
			bytes_per_write = max_bytes_per_write;
		else
			bytes_per_write = data_len - have_write_bytes;

		fwrite((char*)pdata + have_write_bytes, bytes_per_write, 1, fp);
		have_write_bytes += bytes_per_write;

	} while (data_len > have_write_bytes);
}


// global io function, 将读写dat文件封装为统一的接口, 之后只调用下面两个函数即可.

template <typename T>
bool read_dat_file(const std::string dat_file_path, std::vector<T>& data)
{
	if (dat_file_path.empty())
		return false;

	data.clear();

	TGFileBasicIO ob_file;
	GridDataStore ob_store;
	ob_file.set_file_path(dat_file_path);
	if (!ob_file.read_array_data(ob_store) || ob_store.data_len == 0)
		return false;

	if (!ob_store.get_data(data))
		return false;

	return true;
}

template <typename T>
bool write_dat_file(const std::string dat_file_path,  /*数据文件路径*/
					std::vector<T>& data,  /*要写入的数据*/
					std::string str_unit,
					int n_col,
					bool is_big)
{
	if (dat_file_path.empty())
		return false;

	if (data.empty() || n_col < 1)
		return false;

	FILE* fp_write = NULL;
	fopen_s(&fp_write, dat_file_path.c_str(), "wb");
	if (fp_write == NULL)
		return false;

	int n_byte = sizeof(T);
	std::string str_type = typeid(T).name();
	int n_dim = data.size() / n_col;
	TGFileBasicIO ob_file(str_unit, n_col, n_dim, str_type, n_byte, is_big, fp_write);

	if (!ob_file.write_meta_data(true))
		return false;

	__int64 data_len = data.size() * sizeof(T);
	void* pdata = &data[0];
	if (data_len > 0 && !ob_file.write_array_data(n_dim, n_col, n_byte, pdata, data_len))
		return false;

	return true;
}

}	// namespace MCAL::IO
}	// namespace MCAL

#endif