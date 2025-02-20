// $License: GPL-3.0
// $Author: Liu Zhiqiang(from Beihang University)
// $Email: lzhq0930@gmail.com(liuzhiqiang@buaa.edu.cn)
// 


#ifndef MCAL_POLYHEDRAL_MESH_ITERATOR_H
#define MCAL_POLYHEDRAL_MESH_ITERATOR_H

#include <utility>
#include <tuple>

namespace MCAL {   // Mesh Cut Algorithm Libaray, 存放多面体网格切割算法的所有代码
namespace MESH {   // 存放多面体网格实现的相关代码


/*
* 使用一个pair指出几何元素迭代的头和尾, 方便迭代
*/
template <typename I>
class IteratorRange : public std::pair<I, I>
{
    typedef std::pair<I, I> Base;

public:
    typedef I iterator;
    typedef I const_iterator;

    IteratorRange(I b, I e)
        :Base(b, e)
    {}

    IteratorRange(const std::pair<I, I>& rhs) :Base(rhs) {}

    I begin() const
    {
        return this->first;
    }

    I end() const
    {
        return this->second;
    }

    std::size_t size() const
    {
        return static_cast<std::size_t>(std::distance(begin(), end()));
    }

    bool empty() const
    {
        return begin() == end();
    }

    operator std::tuple<I&, I&>()
    {
        return std::tuple<I&, I&>{this->first, this->second};
    }

    operator std::tuple<const I&, const I&>() const
    {
        return std::tuple<const I&, const I&>{this->first, this->second};
    }
};

template <typename T>
IteratorRange<T>
make_range(const T& b, const T& e)
{
    return IteratorRange<T>(b, e);
}

template <typename T>
IteratorRange<T>
make_range(const std::pair<T, T>& rhs)
{
    return IteratorRange<T>(rhs.first, rhs.second);
}

}	// namespace MCAL::MESH
}	// namespace MCAL

#endif