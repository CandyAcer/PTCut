// $Intro: 排序pair, first < second.
// $License: GPL-3.0-or-later
// 
// 
// $Author: Liu Zhiqiang



#ifndef MCAL_SORTED_PAIR_H
#define MCAL_SORTED_PAIR_H

#include <assert.h>

namespace MCAL   // Mesh Cut Algorithm Libaray, 存放多面体网格切割算法的实现
{

	// helper class, 目的是保证二元组的唯一性, 即(a, b)=(b, a).
	// 具体做法是在插入时进行调整, 保证first < second.
	template <typename T>
	class SortedPair
	{
	public:
		T first;
		T second;
		unsigned dim; // dim的取值为0,1,2.

		SortedPair() :dim(0) {}
		SortedPair(T a, T b) :dim(0)
		{
			assert(a != b);
			insert(a);
			insert(b);
		}

		void insert(T i)
		{
			if (dim == 0)
			{
				first = i;
				++dim;
			}
			else if (dim == 1)
			{
				assert(i != first);
				if (i < first)
				{
					second = first;
					first = i;
				}
				else
					second = i;
				++dim;
			}
		}

		std::size_t dimension() { return dim; }

		// 作为map的key, 必须实现严格弱序, 即小于关系.
		bool operator<(const SortedPair<T>& rhs) const
		{
			return first < rhs.first || (first == rhs.first && second < rhs.second);
		}
	};


}	// namespace MCAL

#endif