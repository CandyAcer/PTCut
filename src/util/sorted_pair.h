// $Intro: ����pair, first < second.
// $License: GPL-3.0-or-later
// 
// 
// $Author: Liu Zhiqiang



#ifndef MCAL_SORTED_PAIR_H
#define MCAL_SORTED_PAIR_H

#include <assert.h>

namespace MCAL   // Mesh Cut Algorithm Libaray, ��Ŷ����������и��㷨��ʵ��
{

	// helper class, Ŀ���Ǳ�֤��Ԫ���Ψһ��, ��(a, b)=(b, a).
	// �����������ڲ���ʱ���е���, ��֤first < second.
	template <typename T>
	class SortedPair
	{
	public:
		T first;
		T second;
		unsigned dim; // dim��ȡֵΪ0,1,2.

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

		// ��Ϊmap��key, ����ʵ���ϸ�����, ��С�ڹ�ϵ.
		bool operator<(const SortedPair<T>& rhs) const
		{
			return first < rhs.first || (first == rhs.first && second < rhs.second);
		}
	};


}	// namespace MCAL

#endif