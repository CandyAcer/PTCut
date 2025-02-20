// $License: GPL-3.0
// $Author: Liu Zhiqiang(from Beihang University)
// $Email: lzhq0930@gmail.com(liuzhiqiang@buaa.edu.cn)
// 


#ifndef MCAL_DAT_FILE_IO_H
#define MCAL_DAT_FILE_IO_H

#include <vector>
#include <string>

#include "file_dir_api.h"

namespace MCAL {   // Mesh Cut Algorithm Libaray, ��Ŷ����������и��㷨�����д���
namespace IO {     // I/O

/*
* �ض�����Ķ������ļ���Ϊ *.dat(����vertices.dat)
* һ��dat�ļ���Ӧ�ض������һ������(����vertices.dat��ӦTruncatedGrid��vertices)
*
* dat�ļ����ļ�ͷ�Ͷ������������
* 1) �ļ�ͷ
*    ��¼�ļ���Ԫ����, ����˵���ļ��洢����������
*    ���ȿɱ�, ��ʼ����64�ֽڵı���, ������ֽڻ�������
*
*
*    �ļ�ͷ�е�Ԫ������Ϣʾ��(��vertices.datΪ��)��
*    offset = 128       // ���������ݴӵ�128�ֽڿ�ʼ, 0~127Ϊ�ļ�ͷ
*    units = m          // Ԫ��ֵ�ļ�����λ, m����meter(��), �����ļ�����λ��int
*    ncols = 3          // һ������Ԫ��(��, ��, ��Ԫ)��Ӧ��������
*                       // �˴�һ����Ҫ���μ�¼xyz����ֵ, ����Ϊ3
*
*    dim = 2928         // ��������ļ���Ԫ�ظ���, �˴�˵����¼��2928����
*    type = double      // �������������, �˴�˵��������xyz����double����
*    nbytes = 8         // ����Ԫ�ص�����ռ�����ֽ�(doubleΪ8�ֽ�)
*    big = 0            // trueΪ��˸�ʽ, falseΪС�˸�ʽ
*
*
* 2) ����������
*    ��Ӧĳ������
*
* �ض������ļ��Ķ�д�������๲ͬ���
*     1) TGFileBasicIO�ฺ�����dat�ļ�
*     2) GridDataStore�ฺ�𽫴�dat�ļ��������������ݴ洢��vector��
*
*/


/*
* GridDataStore�ฺ�𽫴�dat�ļ��������������ݴ洢��vector��
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

	// ��p_datablockָ��Ķ��������ݴ����������
	template <typename T>
	bool get_data(std::vector<T>& data);

	void release_data();

public:
	// ��TGFileBasicIO�����ݳ�Ա��Ӧ, ���ٽ��ͺ���
	std::string unit;
	int col;
	int dim;
	std::string type;
	int nbyte;
	__int64 data_len;    // ���ݵ��ֽ���

private:
	// ָ��������ݵ�ָ��, ֮��������ת������
	// ��ָ�벻����������, ��ȡ���������get_data()
	void* p_datablock;
};

template <typename T>
bool
GridDataStore::
get_data(std::vector<T>& data)
{
	if (dim <= 0 || p_datablock == NULL)
		return false;

	// ����ȡ���������뵱ǰ��������һ��, ֱ�ӿ�������
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
* TGFileBasicIO�ฺ�����dat�ļ�
*/
class TGFileBasicIO
{
public:
	// Ĭ�Ϲ��캯��������Ԫ���ݵ�ȱʡֵ
	// Ԫ������Ϊһ���������, û�б�Ҫ�������ݳ�Ա��getter��setter
	TGFileBasicIO()
		: offset(64), unit(std::string()), col(0), dim(0)
		, type(std::string()), nbyte(0), big_endian(false) /*Ĭ��ΪС�˸�ʽ*/
		, fp(NULL), pbuf(NULL)
	{}

	TGFileBasicIO(std::string str_unit, int n_col, int n_dim,
		std::string str_type, int n_byte, bool b_endian, FILE* p);

	~TGFileBasicIO();

	void set_file_path(std::string path) { file_path = path; }

	int estimate_meta_length();  // ����Ԫ���ݵ��ֽ���
	bool read_meta_data();   // ��ȡԪ����
	bool write_meta_data(bool b);  // д��Ԫ����

	bool read_array_data(GridDataStore& data);
	bool write_array_data(int n_dim, int n_col, int n_byte, void* pdata, __int64 data_len);

private:
	bool inverse_data_byte_order(void* pdata, int data_unit, int ele_byte);

	template <typename T>
	bool read_big_data_wrapper(FILE* fp,  /*��������fp*/
		int data_start_byte,  /*���������ݿ�ʼ���ֽ���*/
		std::vector<T>& data,  /*�洢���ݵ�����*/
		__int64 request_bytes,  /*��ȡ���ֽ���*/
		bool is_big  /*�Ƿ�Ϊ���*/);

	void write_big_data_wrapper(const void* pdata, __int64 data_len, FILE* fp_write);

private: // data members
	std::string file_path;   // ��д���ļ�·��

	// �ļ�ͷ��Ԫ����
	int offset;              // ���������ݵ���ʼƫ���ֽ���, ��һ���ؼ���
	std::string unit;        // Ԫ��ֵ��Ӧ�ļ�����λ, ��m, ft, int, �ڶ����ؼ���
	int col;                 // һ������Ԫ��(��, ��, ��Ԫ)��Ӧ��������, 
	// ����Ӧ���������̶�, ��Ϊ1; �̶���Ϊ�̶���ֵ

	int dim;                 // �����ļ���Ԫ��(��, ��, ��Ԫ)����
	// ����Ӧ���������̶�, dimΪ���м���Ԫ�صķ�����֮��
	// ���̶�, ��dimΪ����Ԫ�صĸ���
	// ��֮, dim * col = �����size

	std::string type;        // ����Ԫ�ص���������, ��double, int
	int nbyte;               // ����Ԫ�ص��ֽڳ���, dim * col * nbyte = ������������ֽ���
	bool big_endian;         // trueΪ��˸�ʽ, falseΪС�˸�ʽ

	FILE* fp;     // �ļ�ָ��
	char* pbuf;   // �����ļ�ϵͳ������
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
	// Ԫ���ݵĸ�ʽΪkey=value
	// ����Ԫ���ݵķ���Ϊ��key��ռ���ַ���+value��ռ�ֽ���+1(���з�'\n')
	int len = 8 + sizeof(int) +  /*"offset=?\n", offset�ֶεĳ���*/
		7 + 8 +  /*"units=?\n", ?Ϊstring����, Ԥ��8�ֽ�*/
		7 + sizeof(int) +  /*"ncols=?\n"*/
		5 + sizeof(int) +  /*"dim=?\n"*/
		6 + 8 +  /*"type=?\n", ?Ϊstring����, Ԥ��8�ֽ�*/
		8 + sizeof(int) +  /*"nbytes=?\n", ?�����Ӧ��type����������Ҫ���ֽ��� */
		5 + sizeof(int)  /*"big=0 or 1\n"*/;

	// 64�ֽڶ���
	len = (len / 64 + 1) * 64;
	return len;
}

bool TGFileBasicIO::read_meta_data()
{
	if (fp == NULL)
		return false;

	// ���ļ�ָ���ƶ����ļ��Ŀ�ʼλ��
	_fseeki64(fp, 0, SEEK_SET);

	// offsetΪ64�ı���, ��64�ֽ��Ƕ��Ļ�����λ
	char buffer[64];
	int buf_size = 64;
	memset(buffer, '\0', buf_size);

	int fp_location = 0;           // �ļ�ָ���λ��, ���ֽ�Ϊ��λ
	while (fp_location < offset)
	{
		// ��fp���������ж�ȡһ�е�buffer��
		fgets(buffer, buf_size, fp);
		// �޸�fp_locationΪ�ļ�ָ��ĵ�ǰλ��
		fp_location = ftell(fp);

		// �ָ�buffer�е��ַ���, �ָ����ַ�����cutted_strs��
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

	// У��Ԫ������ȷ��
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

	// ��Ԫ����ȫ�����Ϊ'\0'
	char* pmeta = new char[offset];
	memset(pmeta, '\0', offset);
	fwrite(pmeta, offset, 1, fp);

	_fseeki64(fp, 0L, SEEK_SET);

	// ��дmeta���ֶ�

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

	// ��ȡԪ����
	if (!read_meta_data())
		return false;

	ob_store.unit = unit;
	ob_store.col = col;
	ob_store.dim = dim;
	ob_store.type = type;
	ob_store.nbyte = nbyte;
	ob_store.data_len = (__int64)nbyte * col * dim;
	ob_store.p_datablock = NULL;

	// �ļ�ָ���ƶ������������ݵ���ʼλ��
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
	// �����Ƿ��д���
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
read_big_data_wrapper(FILE* fp,  /*��������fp*/
					  int data_start_byte,  /*���������ݿ�ʼ���ֽ���*/
					  std::vector<T>& data,  /*�洢���ݵ�����*/
					  __int64 request_bytes,  /*��ȡ���ֽ���*/
					  bool is_big  /*�Ƿ�Ϊ���*/)
{
	if (fp == NULL || request_bytes <= 0)
		return false;

	int type_len = sizeof(T);
	data.resize(request_bytes / type_len);

	char* p = (char*)&data[0];

	_fseeki64(fp, data_start_byte, SEEK_SET);

	// ���ļ���Ҫ�ִζ�
	__int64 have_read_bytes = 0;                      // �Ѷ�ȡ���ֽ���
	__int64 max_bytes_per_read = 256 * 1024 * 1024;   // һ������ȡ���ֽ���
	__int64 bytes_per_read = 0;                       // һ����Ҫ��ȡ���ֽ���

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


// global io function, ����дdat�ļ���װΪͳһ�Ľӿ�, ֮��ֻ��������������������.

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
bool write_dat_file(const std::string dat_file_path,  /*�����ļ�·��*/
					std::vector<T>& data,  /*Ҫд�������*/
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