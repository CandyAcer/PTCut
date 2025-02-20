// $License: GPL-3.0
// $Author: Liu Zhiqiang(from Beihang University)
// $Email: lzhq0930@gmail.com(liuzhiqiang@buaa.edu.cn)
// 


#ifndef MCAL_ALGO_INTERSECTION_TYPE_H
#define MCAL_ALGO_INTERSECTION_TYPE_H

namespace MCAL    // Mesh Cut Algorithm Libaray, ��Ŷ����������и��㷨��ʵ��
{

enum IntersectionType 
{ 
	ON_VERTEX, ON_EDGE, ON_FACE, EMPTY, COPLANAR_TRIANGLES 
};

}	// namespace MCAL

#endif