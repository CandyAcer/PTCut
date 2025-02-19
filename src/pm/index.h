// $License: GPL-3.0
// $Author: Liu Zhiqiang(from Beihang University)
// $Email: lzhq0930@gmail.com
// 
// Geometric elements(vertex, edge, face, cell) in polyhedral mesh are represented
// by index class, which encapsulate an integer and define its behavior.
// 


#ifndef MCAL_POLYHEDRAL_MESH_INDEX_H
#define MCAL_POLYHEDRAL_MESH_INDEX_H

#include <iostream>
#include <limits>
#include <cstddef>
#include <functional>

namespace MCAL {   // Mesh Cut Algorithm Libaray, 存放多面体网格切割算法的所有代码
namespace MESH {   // 存放多面体网格实现的相关代码

/*
* Base class for vertex, edge, face and cell index.
* 负责对一个正整数进行封装
* 基本功能: 构造, 自增自减, 类型转换
*/
class BaseIndex
{
public:
    typedef std::uint32_t size_type;

public:
    // 构造函数默认构造非法(invalid)值, 非法值只有一个.
    // index是32位无符号整数, 多面体网格包含的几何元素一般不会超过这个范围.
    // 故非法值选择最大值0xffffffff.
    explicit BaseIndex(size_type idx = std::numeric_limits<size_type>::max())
        : index(idx)
    {}

    BaseIndex(const BaseIndex& rhs) :index(rhs.index) {}

    BaseIndex& operator=(const BaseIndex& rhs)
    {
        index = rhs.index;
        return *this;
    }

    // 类型转换
    operator size_type() const { return index; }

    void reset() { index = std::numeric_limits<size_type>::max(); }

    bool is_valid() const
    {
        size_type inf = std::numeric_limits<size_type>::max();
        return index != inf;
    }

    // 自增自减(前缀后缀), 无边界检查
    BaseIndex& operator++() { ++index; return *this; }
    BaseIndex& operator--() { --index; return *this; }
    BaseIndex operator++(int) { BaseIndex tmp(*this); ++index; return tmp; }
    BaseIndex operator--(int) { BaseIndex tmp(*this); --index; return tmp; }

    // +n or -n, 无边界检查
    BaseIndex operator+=(std::ptrdiff_t n)
    {
        index = size_type(std::ptrdiff_t(index) + n);
        return *this;
    }

protected:
    size_type index;       // index本质上是一个正整数
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

    // 以下代码的目的是禁止PMVertexIndex与本类型以外的数据做比较, 避免不必要的错误
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

    // 测试输出
    friend std::ostream& operator<<(std::ostream& os, PMVertexIndex const& v)
    {
        return (os << 'V' << (size_type)v);
    }
};

/*
* Implementation for PolyhedralMesh::EdgeIndex.
*/
class PMEdgeIndex : public BaseIndex
{
public:
    explicit PMEdgeIndex(size_type idx = std::numeric_limits<size_type>::max())
        : BaseIndex(idx)
    {}

    // 以下代码的目的是禁止PMEdgeIndex与本类型以外的数据做比较, 避免不必要的错误
    template<typename T> bool operator==(const T&) const = delete;
    template<typename T> bool operator!=(const T&) const = delete;
    template<typename T> bool operator<(const T&) const = delete;

    bool operator==(const PMEdgeIndex& rhs) const
    {
        return this->index == rhs.index;
    }

    bool operator!=(const PMEdgeIndex& rhs) const
    {
        return this->index != rhs.index;
    }

    bool operator<(const PMEdgeIndex& rhs) const
    {
        return this->index < rhs.index;
    }

    // 测试输出
    friend std::ostream& operator<<(std::ostream& os, PMEdgeIndex const& e)
    {
        return (os << 'E' << (size_type)e);
    }
};

/*
* Implementation for PolyhedralMesh::FaceIndex.
*/
class PMFaceIndex : public BaseIndex
{
public:
    explicit PMFaceIndex(size_type idx = std::numeric_limits<size_type>::max())
        : BaseIndex(idx)
    {}

    // 以下代码的目的是禁止PMFaceIndex与本类型以外的数据做比较, 避免不必要的错误
    template<typename T> bool operator==(const T&)const = delete;
    template<typename T> bool operator!=(const T&)const = delete;
    template<typename T> bool operator<(const T&)const = delete;

    bool operator==(const PMFaceIndex& rhs) const
    {
        return this->index == rhs.index;
    }

    bool operator!=(const PMFaceIndex& rhs) const
    {
        return this->index != rhs.index;
    }

    bool operator<(const PMFaceIndex& rhs) const
    {
        return this->index < rhs.index;
    }

    // 测试输出
    friend std::ostream& operator<<(std::ostream& os, PMFaceIndex const& f)
    {
        return (os << 'F' << (size_type)f);
    }
};

/*
* Implementation for PolyhedralMesh::CellIndex.
*/
class PMCellIndex : public BaseIndex
{
public:
    explicit PMCellIndex(size_type idx = std::numeric_limits<size_type>::max())
        : BaseIndex(idx)
    {}

    // 以下代码的目的是禁止PMFaceIndex与本类型以外的数据做比较, 避免不必要的错误
    template<typename T> bool operator==(const T&)const = delete;
    template<typename T> bool operator!=(const T&)const = delete;
    template<typename T> bool operator<(const T&)const = delete;

    bool operator==(const PMCellIndex& rhs)const
    {
        return this->index == rhs.index;
    }

    bool operator!=(const PMCellIndex& rhs)const
    {
        return this->index != rhs.index;
    }

    bool operator<(const PMCellIndex& rhs)const
    {
        return this->index < rhs.index;
    }

    // 测试输出
    friend std::ostream& operator<<(std::ostream& os, PMCellIndex const& c)
    {
        return (os << 'C' << (size_type)c);
    }
};

}	// namespace MCAL::MESH
}	// namespace MCAL



namespace std
{
// 算法为了记录某些中间数据, 需要用一些容器
// 当使用unordered_map或unordered_set时, 需要对key进行hash以及确保键值不重复
// 定义如下：
// template<class Key,     // 1.key值
//          class Ty,      // 2.value值
//          class Hash = std::hash<Key>,       // 3.如何对key进行hash
//          class Pred = std::equal_to<Key>,   // 4.判断key是否相等
//          class Alloc = std::allocator<std::pair<const Key, Ty>>  // 5.内存分配
//         >
// class unordered_map;
// 
// 对于自定义的类型, 如果想作为unordered_map的key, 必须实现3和4两种功能
// 实现4可以在类中重载operator==
// 实现3的方法很多, 其中一种是定制自己的模板类std::hash<Key>
// 这种方式较好, 其他人只需包含本文件, 就可以直接使用, 无需额外工作
//

// 提供std::hash<Key>的特化版本

template <>
struct hash<MCAL::MESH::PMVertexIndex>
    :public std::unary_function< MCAL::MESH::PMVertexIndex, std::size_t>
{
    std::size_t operator()(const MCAL::MESH::PMVertexIndex& i) const
    {
        return i;
    }
};

template <>
struct hash<MCAL::MESH::PMEdgeIndex>
    :public std::unary_function< MCAL::MESH::PMEdgeIndex, std::size_t>
{
    std::size_t operator()(const MCAL::MESH::PMEdgeIndex& i) const
    {
        return i;
    }
};

template <>
struct hash<MCAL::MESH::PMFaceIndex>
    :public std::unary_function< MCAL::MESH::PMFaceIndex, std::size_t>
{
    std::size_t operator()(const MCAL::MESH::PMFaceIndex& i) const
    {
        return i;
    }
};

template <>
struct hash<MCAL::MESH::PMCellIndex>
    :public std::unary_function< MCAL::MESH::PMCellIndex, std::size_t>
{
    std::size_t operator()(const MCAL::MESH::PMCellIndex& i) const
    {
        return i;
    }
};

}   // namespace std

#endif