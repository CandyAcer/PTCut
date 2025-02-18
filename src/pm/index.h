// $License: GPL-3.0
// $Author: Liu Zhiqiang (from Beihang University)
// $Email: lzhq0930@gmail.com
// $Description: 
// 


#ifndef MCAL_POLYHEDRAL_MESH_INDEX_H
#define MCAL_POLYHEDRAL_MESH_INDEX_H

#include <iostream>
#include <limits>
#include <cstddef>
#include <functional>


namespace MCAL {   // Mesh Cut Algorithm Libaray, ��Ŷ����������и��㷨�����д���
namespace MESH {   // ��Ŷ���������ʵ�ֵ���ش���

/*
* Base class for vertex, edge, face and cell index.
* �����һ�����������з�װ
* ��������: ����, �����Լ�, ����ת��
*/
class BaseIndex
{
public:
    typedef std::uint32_t size_type;

public:
    // ���캯��Ĭ�Ϲ���Ƿ�(invalid)ֵ, �Ƿ�ֵֻ��һ��.
    // index��32λ�޷�������, ��������������ļ���Ԫ��һ�㲻�ᳬ�������Χ.
    // �ʷǷ�ֵѡ�����ֵ0xffffffff.
    explicit BaseIndex(size_type idx = std::numeric_limits<size_type>::max())
        : index(idx)
    {}

    BaseIndex(const BaseIndex& rhs) :index(rhs.index) {}

    BaseIndex& operator=(const BaseIndex& rhs)
    {
        index = rhs.index;
        return *this;
    }

    // ����ת��
    operator size_type() const { return index; }

    // ��index����Ϊ�Ƿ�ֵ
    void reset() { index = std::numeric_limits<size_type>::max(); }

    // �ж�index�Ƿ�Ϸ�
    bool is_valid() const
    {
        size_type inf = std::numeric_limits<size_type>::max();
        return index != inf;
    }

    // �����Լ�(ǰ׺��׺), �ޱ߽���
    BaseIndex& operator++() { ++index; return *this; }
    BaseIndex& operator--() { --index; return *this; }
    BaseIndex operator++(int) { BaseIndex tmp(*this); ++index; return tmp; }
    BaseIndex operator--(int) { BaseIndex tmp(*this); --index; return tmp; }

    // +n or -n, �ޱ߽���
    BaseIndex operator+=(std::ptrdiff_t n)
    {
        index = size_type(std::ptrdiff_t(index) + n);
        return *this;
    }

protected:
    size_type index;       // index��������һ��������
};

std::size_t hash_value(const BaseIndex& i)
{
    std::size_t ret = i;
    return ret;
}

/*
* Implementation for PolyhedralMesh::VertexIndex.
*/
class PMVertexIndex : public BaseIndex
{
public:
    explicit PMVertexIndex(size_type idx = std::numeric_limits<size_type>::max())
        : BaseIndex(idx)
    {}

    // ���´����Ŀ���ǽ�ֹPMVertexIndex�뱾����������������Ƚ�, ���ⲻ��Ҫ�Ĵ���
    template<typename T> bool operator==(const T&) const = delete;
    template<typename T> bool operator!=(const T&) const = delete;
    template<typename T> bool operator<(const T&) const = delete;

    bool operator==(const PMVertexIndex& rhs) const
    {
        return this->index == rhs.index;
    }

    bool operator!=(const PMVertexIndex& rhs) const
    {
        return this->index != rhs.index;
    }

    bool operator<(const PMVertexIndex& rhs) const
    {
        return this->index < rhs.index;
    }

    // �������
    friend std::ostream& operator<<(std::ostream& os, PMVertexIndex const& v)
    {
        return (os << 'V' << (size_type)v);
    }
};

/*
* Implementation for PolyhedralMesh::EdgeIndex.
* Note: CGAL��ʵ����HalfedgeIndex, �������ڱ��������Edge��Ӧ����Halfedge,
*       �������������edge���پ�����������, �ʲ�������.
*/
class PGEdgeIndex : public BaseIndex
{
public:
    explicit PGEdgeIndex(size_type idx = std::numeric_limits<size_type>::max())
        : BaseIndex(idx)
    {
    }

    // ���´����Ŀ���ǽ�ֹPGEdgeIndex�뱾����������������Ƚ�, ���ⲻ��Ҫ�Ĵ���
    template<typename T> bool operator==(const T&) const = delete;
    template<typename T> bool operator!=(const T&) const = delete;
    template<typename T> bool operator<(const T&) const = delete;

    bool operator==(const PGEdgeIndex& rhs) const
    {
        return this->index == rhs.index;
    }

    bool operator!=(const PGEdgeIndex& rhs) const
    {
        return this->index != rhs.index;
    }

    bool operator<(const PGEdgeIndex& rhs) const
    {
        return this->index < rhs.index;
    }

    // �������
    friend std::ostream& operator<<(std::ostream& os, PGEdgeIndex const& e)
    {
        return (os << 'E' << (size_type)e);
    }
};

/*
* Implementation for PolyhedralMesh::FaceIndex.
*/
class PGFaceIndex : public BaseIndex
{
public:
    explicit PGFaceIndex(size_type idx = std::numeric_limits<size_type>::max())
        : BaseIndex(idx)
    {
    }

    // ���´����Ŀ���ǽ�ֹPGFaceIndex�뱾����������������Ƚ�, ���ⲻ��Ҫ�Ĵ���
    template<typename T> bool operator==(const T&)const = delete;
    template<typename T> bool operator!=(const T&)const = delete;
    template<typename T> bool operator<(const T&)const = delete;

    bool operator==(const PGFaceIndex& rhs) const
    {
        return this->index == rhs.index;
    }

    bool operator!=(const PGFaceIndex& rhs) const
    {
        return this->index != rhs.index;
    }

    bool operator<(const PGFaceIndex& rhs) const
    {
        return this->index < rhs.index;
    }

    // �������
    friend std::ostream& operator<<(std::ostream& os, PGFaceIndex const& f)
    {
        return (os << 'F' << (size_type)f);
    }
};

/*
* Implementation for PolyhedralMesh::CellIndex.
*/
class PGCellIndex : public BaseIndex
{
public:
    explicit PGCellIndex(size_type idx = std::numeric_limits<size_type>::max())
        : BaseIndex(idx)
    {
    }

    // ���´����Ŀ���ǽ�ֹPGFaceIndex�뱾����������������Ƚ�, ���ⲻ��Ҫ�Ĵ���
    template<typename T> bool operator==(const T&)const = delete;
    template<typename T> bool operator!=(const T&)const = delete;
    template<typename T> bool operator<(const T&)const = delete;

    bool operator==(const PGCellIndex& rhs)const
    {
        return this->index == rhs.index;
    }

    bool operator!=(const PGCellIndex& rhs)const
    {
        return this->index != rhs.index;
    }

    bool operator<(const PGCellIndex& rhs)const
    {
        return this->index < rhs.index;
    }

    // �������
    friend std::ostream& operator<<(std::ostream& os, PGCellIndex const& c)
    {
        return (os << 'C' << (size_type)c);
    }
};

}	// namespace MCAL::GRID
}	// namespace MCAL



namespace std
{
// �㷨Ϊ�˼�¼ĳЩ�м�����, ��Ҫ��һЩ����
// ��ʹ��unordered_map��unordered_setʱ, ��Ҫ��key����hash�Լ�ȷ����ֵ���ظ�
// �������£�
// template<class Key,     // 1.keyֵ
//          class Ty,      // 2.valueֵ
//          class Hash = std::hash<Key>,       // 3.��ζ�key����hash
//          class Pred = std::equal_to<Key>,   // 4.�ж�key�Ƿ����
//          class Alloc = std::allocator<std::pair<const Key, Ty>>  // 5.�ڴ����
//         >
// class unordered_map;
// 
// �����Զ��������, �������Ϊunordered_map��key, ����ʵ��3��4���ֹ���
// ʵ��4��������������operator==
// ʵ��3�ķ����ܶ�, ����һ���Ƕ����Լ���ģ����std::hash<Key>
// ���ַ�ʽ�Ϻ�, ������ֻ��������ļ�, �Ϳ���ֱ��ʹ��, ������⹤��
//

// �ṩstd::hash<Key>���ػ��汾

template <>
struct hash<MCAL::GRID::PMVertexIndex>
    :public std::unary_function< MCAL::GRID::PMVertexIndex, std::size_t>
{
    std::size_t operator()(const MCAL::GRID::PMVertexIndex& i) const
    {
        return i;
    }
};

template <>
struct hash<MCAL::GRID::PGEdgeIndex>
    :public std::unary_function< MCAL::GRID::PGEdgeIndex, std::size_t>
{
    std::size_t operator()(const MCAL::GRID::PGEdgeIndex& i) const
    {
        return i;
    }
};

template <>
struct hash<MCAL::GRID::PGFaceIndex>
    :public std::unary_function< MCAL::GRID::PGFaceIndex, std::size_t>
{
    std::size_t operator()(const MCAL::GRID::PGFaceIndex& i) const
    {
        return i;
    }
};

template <>
struct hash<MCAL::GRID::PGCellIndex>
    :public std::unary_function< MCAL::GRID::PGCellIndex, std::size_t>
{
    std::size_t operator()(const MCAL::GRID::PGCellIndex& i) const
    {
        return i;
    }
};

}   // namespace std

#endif