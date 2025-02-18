// $Intro: 本文件定义线段与三角形的相交类型.
// $License: GPL-3.0-or-later
// 
// 
// $Author: Liu Zhiqiang



#ifndef MCAL_ALGORITHM_INTERSECTION_TYPE_H
#define MCAL_ALGORITHM_INTERSECTION_TYPE_H


namespace MCAL    // Mesh Cut Algorithm Libaray, 存放多面体网格切割算法的实现
{

	enum IntersectionType { ON_VERTEX, ON_EDGE, ON_FACE, EMPTY, COPLANAR_TRIANGLES };

}	// namespace MCAL

#endif