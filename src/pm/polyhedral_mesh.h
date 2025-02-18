// $Intro: 多面体网格的代码实现, 借鉴了Mallison和CGAL::Surface_mesh的设计思路.
// $License: GPL-3.0-or-later
// 
// 
// $Author: Liu Zhiqiang



#ifndef MCAL_POLYHEDRON_GRID_IMPL_H
#define MCAL_POLYHEDRON_GRID_IMPL_H

#include <CGAL/Kernel_traits.h>
#include <CGAL/Kernel/global_functions_3.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <CGAL/boost/graph/iterator.h>

#include <array>
#include <queue>
#include <set>
#include <map>
#include <unordered_map>

#include "pg_index.h"
#include "pg_properties.h"
#include "pg_iterator.h"
#include "utils/sorted_pair.h"



namespace MCAL {   // Mesh Cut Algorithm Libaray, 存放多面体网格切割算法的所有代码
    namespace GRID {   // 存放多面体网格实现的相关代码


        /*
        * 多面体网格的实现代码, 包含属性管理, 拓扑查询, 拓扑修改等功能.
        *
        * @param P：点的类型, 可以是double三元组或cgal的精确类型.
        * @param T：CellTopo中的halfedges使用局部索引指示引用的vertex与edge,
        *           局部索引的数据类型, 其数据宽度决定一个单元包含多少几何元素.
        *
        */
        template <typename P, typename T>
        class PolyhedronGrid
        {
            typedef PolyhedronGrid<P, T>   Self;

        public:

            template <typename I, typename T>
            struct PropertyMap : public PropertyMapBase<I, T>
            {
                typedef PropertyMapBase<I, T> Base;
                typedef typename Base::reference reference;

                PropertyMap() = default;
                PropertyMap(const Base& pm) : Base(pm) {}
            };

            template <typename Key, typename T>
            struct GetPropertyMap
            {
                typedef PropertyMap<Key, T>  type;
            };

            enum CellType { IN = 0, OUT, UNDEFINED };

            struct CellInfo
            {
                int iidx, jidx, kidx; // Cell的IJK属性
                CellType ctype;

                CellInfo(int i = std::numeric_limits<int>::max(),
                    int j = std::numeric_limits<int>::max(),
                    int k = std::numeric_limits<int>::max(),
                    CellType type = UNDEFINED)
                    :iidx(i), jidx(j), kidx(k), ctype(type)
                {
                }
            };

            typedef P                                                Point;
            typedef typename CGAL::Kernel_traits<Point>::Kernel      Kernel;
            typedef typename Kernel::Vector_3                        Normal;
            typedef T                                                bias_type;
            typedef std::uint32_t                                    size_type;

            typedef PGVertexIndex          VertexIndex;
            typedef PGEdgeIndex            EdgeIndex;
            typedef PGFaceIndex            FaceIndex;
            typedef PGCellIndex            CellIndex;

        public:// 拓扑定义

            // 多面体网格的halfedge用一个二元组表示
            // 指示其所在的cell以及CellTopo的halfedges数组中的下标
            struct Halfedge
            {
                CellIndex he_cell;  // 所在cell
                size_type he_off;   // CellTopo的halfedges数组中的下标

                Halfedge(CellIndex c = CellIndex(),
                    size_type i = std::numeric_limits<size_type>::max())
                    :he_cell(c), he_off(i)
                {
                }

                Halfedge(const Halfedge& rhs) { operator=(rhs); }

                Halfedge& operator=(const Halfedge& rhs)
                {
                    he_cell = rhs.he_cell;
                    he_off = rhs.he_off;

                    return *this;
                }

                bool operator!=(const Halfedge& rhs)
                {
                    return !(he_cell == rhs.he_cell && he_off == rhs.he_off);
                }

                bool operator==(const Halfedge& rhs)
                {
                    return he_cell == rhs.he_cell && he_off == rhs.he_off;
                }

                // 作为map的key, 必须实现严格弱序, 即小于关系.
                bool operator<(const Halfedge& rhs) const
                {
                    return he_cell < rhs.he_cell || (he_cell == rhs.he_cell && he_off < rhs.he_off);
                }

                // for debug
                friend std::ostream& operator<<(std::ostream& os, Halfedge const& he)
                {
                    return os << "he(" << (size_type)he.he_cell
                        << ", " << he.he_off << ")";
                }
            };

            struct VertexTopo
            {
                Halfedge halfedge;  // vertex为halfedge的source.
            };

            struct EdgeTopo
            {
                Halfedge halfedge;  // edge对应的某条halfedge(无固定数量关系, 表面网格为1:2)
            };

            struct FaceTopo
            {
                Halfedge halfedge;
                std::array<CellIndex, 2> incident_cells;  // face关联的cell
            };

            struct CellTopo
            {
                // 初始状态除了offset的容器都为空, offset有一个元素
                CellTopo()
                {
                    offset.push_back(std::numeric_limits<size_type>::max());
                }

                std::vector<VertexIndex>      vertices;    // cell包含的所有vertex
                std::vector<EdgeIndex>        edges;       // cell包含的所有edge
                std::vector<FaceIndex>        faces;       // cell包含的所有face
                std::vector<size_type>        offset;      // face在halfedges中的开始位置, 比faces多1
                std::vector<std::pair<bias_type, bias_type>> halfedges;
            };

        public:
            // Constructor.
            PolyhedronGrid();

            // Copy constructor: copies `rhs` to `*this`. Performs a deep copy of all properties.
            PolyhedronGrid(const PolyhedronGrid& rhs) { *this = rhs; }

            // assigns `rhs` to `*this`. Performs a deep copy of all properties.
            PolyhedronGrid& operator=(const PolyhedronGrid& rhs);

        private: // 几何元素的迭代

            template <typename Index>
            class IndexIterator
                :public boost::iterator_facade<IndexIterator<Index>,/*Derived Class*/
                Index,/*迭代器值的类型*/
                std::random_access_iterator_tag/*迭代器类型*/>
            {
            public:
                typedef boost::iterator_facade<IndexIterator<Index>,
                    Index,
                    std::random_access_iterator_tag
                > Facade;

                IndexIterator() :idx(), grid(nullptr) {}
                IndexIterator(const Index& i, const PolyhedronGrid* pg)
                    :idx(i), grid(pg)
                {
                    if (grid && grid->has_garbage())
                    {
                        // 当idx被标记为已删除时, 应当+1, 直至找到不越界且未删除的索引
                        while (grid->has_valid_index(idx) && grid->is_removed(idx))
                            ++idx;
                    }
                }

            private:
                friend class boost::iterator_core_access;

                void increment()
                {
                    ++idx;
                    assert(grid != nullptr);

                    if (grid->has_garbage())
                    {
                        while (grid->has_valid_index(idx) && grid->is_removed(idx))
                            ++idx;
                    }
                }

                void decrement()
                {
                    --idx;
                    assert(grid != nullptr);

                    if (grid->has_garbage())
                    {
                        while (grid->has_valid_index(idx) && grid->is_removed(idx))
                            --idx;
                    }
                }

                void advance(std::ptrdiff_t n)
                {
                    assert(grid != nullptr);

                    if (grid->has_garbage())
                    {
                        if (n > 0)
                            for (std::ptrdiff_t i = 0; i < n; ++i)
                                increment();
                        else
                            for (std::ptrdiff_t i = 0; i < -n; ++i)
                                decrement();
                    }
                    else
                        idx += n;
                }

                std::ptrdiff_t distance_to(const IndexIterator& other) const
                {
                    if (grid->has_garbage())
                    {
                        bool forward = (other.idx > this->idx);

                        std::ptrdiff_t ans = 0;
                        IndexIterator it = *this;
                        while (!it.equal(other))
                        {
                            if (forward)
                            {
                                ++it;
                                ++ans;
                            }
                            else
                            {
                                --it;
                                --ans;
                            }
                        }
                        return ans;
                    }
                    else
                        return std::ptrdiff_t(other.idx) - std::ptrdiff_t(this->idx);
                }

                bool equal(const IndexIterator& other) const
                {
                    return this->idx == other.idx;
                }

                Index& dereference() const { return const_cast<Index&>(idx); }

            private:
                Index idx;                    // 迭代的几何元素index
                const PolyhedronGrid* grid;   // index所在的网格, 便于访问index是否被标记删除
            };

        public:
            typedef IndexIterator<VertexIndex>        VertexIterator;
            typedef IteratorRange<VertexIterator>     VertexRange;

            typedef IndexIterator<EdgeIndex>          EdgeIterator;
            typedef IteratorRange<EdgeIterator>       EdgeRange;

            typedef IndexIterator<FaceIndex>          FaceIterator;
            typedef IteratorRange<FaceIterator>       FaceRange;

            typedef IndexIterator<CellIndex>          CellIterator;
            typedef IteratorRange<CellIterator>       CellRange;

            // Start iterator for vertices.
            VertexIterator vertices_begin() const
            {
                return VertexIterator(VertexIndex(0), this);
            }

            // End iterator for vertices.
            VertexIterator vertices_end() const
            {
                return VertexIterator(VertexIndex(num_vertices()), this);
            }

            // The iterator range of the vertices of the PolyhedronGrid.
            VertexRange vertices() const
            {
                return make_range(vertices_begin(), vertices_end());
            }

            // Start iterator for edges.
            EdgeIterator edges_begin() const
            {
                return EdgeIterator(EdgeIndex(0), this);
            }

            // End iterator for edges.
            EdgeIterator edges_end() const
            {
                return EdgeIterator(EdgeIndex(num_edges()), this);
            }

            // The iterator range of the edges of the PolyhedronGrid.
            EdgeRange edges() const
            {
                return make_range(edges_begin(), edges_end());
            }

            // Start iterator for faces.
            FaceIterator faces_begin() const
            {
                return FaceIterator(FaceIndex(0), this);
            }

            // End iterator for faces.
            FaceIterator faces_end() const
            {
                return FaceIterator(FaceIndex(num_faces()), this);
            }

            // The iterator range of the faces of the PolyhedronGrid.
            FaceRange faces() const
            {
                return make_range(faces_begin(), faces_end());
            }

            // Start iterator for cells.
            CellIterator cells_begin() const
            {
                return CellIterator(CellIndex(0), this);
            }

            // End iterator for cells.
            CellIterator cells_end() const
            {
                return CellIterator(CellIndex(num_cells()), this);
            }

            // The iterator range of the cells of the PolyhedronGrid.
            CellRange cells() const
            {
                return make_range(cells_begin(), cells_end());
            }

        public:
            /*-------------------- 内存管理模块 -------------------------------
            *
            * 该模块的主要功能是管理几何元素的增减, 并在此过程中管理内存
            *
            * 在算法执行过程中, 会删除和添加几何元素
            * 几何元素的的属性统一放在一个vector中, 在数组的随机位置删除的性能消耗较大
            * 故采用以下方法：
            *     1) 删除某个几何元素时, 并不实际删除, 而是将其标记为removed
            *     2) 被标记为删除的index在添加时可以重新分配(循环利用)
            *
            */


            // 此处解释一下vertices_freelist的作用
            // 在算法增减元素的过程中, 可以按删除的先后顺序将VertexIndex进行排列
            // 新增元素时, 可以将最近被删除的分配给新元素, 此时这个线性表的行为类似于栈
            // 逻辑上看, 这个排列是一个线性表, 似乎可以用链表或数组进行记录, 但消耗较大
            //
            //
            // 我们不用数组或链表实现, 一个无符号整数vertices_freelist足矣
            // 它类似于栈顶指针, 始终记录最近被删除的元素的索引
            // 之所以命名为freelist, 说明该变量在模拟一个内存释放(free)的链表(list)
            //
            // 但这个排列是一种线性关系, 节点间的前驱后继关系怎么表示呢？
            // 我们注意到, 标记为删除的VertexIndex对应的所有属性值都是不可使用的
            // 那么可以借VertexTopo来记录下一个结点
            //
            // 这种方式很巧妙, 通过借用废弃的内存空间, 只用一个无符号整数便实现出了一个线性表
            // 极大地减少了内存占用
            //
            //
            // CGAL::Surface_mesh中的这种方式可以从两个方面对内存进行优化：
            // 1)删除时仅作标记, 避免了数组在删除时移动内存的消耗
            // 2)vertices_freelist和废弃内存的结合, 将一个数组变为一个整数, 类似于滚动DP
            //
            // edges_freelist,faces_freelist,cells_freelist同理
            // freelist的初始值应为非法值(即最大值), 可视为链表的空结点(空结点可简化操作)
            //


            VertexIndex add_vertex()
            {
                size_type inf = std::numeric_limits<size_type>::max();

                // recycle为true表示网格会回收利用已标记为删除的Index
                // vertices_freelist != inf 表示删除过元素, 可以进行分配
                if (recycle && (vertices_freelist != inf))
                {
                    VertexIndex v(vertices_freelist);   // 分配最近已删除的索引
                    vertices_freelist = vtopo[v].halfedge.he_off;   //更新freelist
                    --removed_vertices;
                    vremoved[v] = false;   // 重置removed标记
                    vprops.reset(v);       // 重置v对应的所有属性值, 后续可修改
                    return v;
                }
                else    // 不可分配时push_back即可(index+1), 较简单
                {
                    vprops.push_back();
                    return VertexIndex(num_vertices() - 1);
                }
            }

            VertexIndex add_vertex(const Point& p)
            {
                VertexIndex v = add_vertex();
                vpoint[v] = p;
                return v;
            }

            EdgeIndex add_edge()
            {
                size_type inf = std::numeric_limits<size_type>::max();

                if (recycle && (edges_freelist != inf))
                {
                    EdgeIndex e(edges_freelist);
                    edges_freelist = etopo[e].halfedge.he_off;
                    --removed_edges;
                    eremoved[e] = false;
                    eprops.reset(e);
                    return e;
                }
                else
                {
                    eprops.push_back();
                    return EdgeIndex(num_edges() - 1);
                }
            }

            FaceIndex add_face()
            {
                size_type inf = std::numeric_limits<size_type>::max();

                if (recycle && (faces_freelist != inf))
                {
                    FaceIndex f(faces_freelist);
                    faces_freelist = ftopo[f].halfedge.he_off;
                    --removed_faces;
                    fremoved[f] = false;
                    fprops.reset(f);
                    return f;
                }
                else
                {
                    fprops.push_back();
                    return FaceIndex(num_faces() - 1);
                }
            }

            CellIndex add_cell()
            {
                size_type inf = std::numeric_limits<size_type>::max();

                if (recycle && (cells_freelist != inf))
                {
                    CellIndex c(cells_freelist);
                    cells_freelist = ctopo[c].offset[0];
                    --removed_cells;
                    cremoved[c] = false;
                    cprops.reset(c);
                    return c;
                }
                else
                {
                    cprops.push_back();
                    return CellIndex(num_cells() - 1);
                }
            }

            // 删除Vertex, 仅标记为删除, 其关联的所有属性仍在内存中
            void remove_vertex(VertexIndex v)
            {
                vremoved[v] = true; ++removed_vertices; garbage = true;

                // VertexTopo记录了一条halfedge, Halfedge有两个数据成员
                // 我们用he_off来记录vertices_freelist(因为类型都是size_type)
                vtopo[v].halfedge.he_off = vertices_freelist;
                vertices_freelist = (size_type)v;
            }

            void remove_edge(EdgeIndex e)
            {
                eremoved[e] = true; ++removed_edges; garbage = true;
                etopo[e].halfedge.he_off = edges_freelist;
                edges_freelist = (size_type)e;
            }

            void remove_face(FaceIndex f)
            {
                fremoved[f] = true; ++removed_faces; garbage = true;
                ftopo[f].halfedge.he_off = faces_freelist;
                faces_freelist = (size_type)f;
            }

            void remove_cell(CellIndex c)
            {
                cremoved[c] = true; ++removed_cells; garbage = true;
                // CellTopo中不记录halfedge, 所以不能像前三个函数一样
                // CellTopo的默认构造函数构造出的对象中offset有一个元素, 其余为空
                // 可以使用这个元素记录
                ctopo[c].offset[0] = cells_freelist;
                cells_freelist = (size_type)c;
            }

            // The number of used and removed vertices in the gird.
            // 多面体网格中几何元素(vertex, edge, face, cell)的总数(包括已标记删除的)
            // 实质是已经分配的index总数
            size_type num_vertices() const { return (size_type)vprops.size(); }
            size_type num_edges() const { return (size_type)eprops.size(); }
            size_type num_faces() const { return (size_type)fprops.size(); }
            size_type num_cells() const { return (size_type)cprops.size(); }

            // 多面体网格中被标记为已删除的几何元素(vertex, edge, face, cell)总数
            size_type number_of_removed_vertices() const { return removed_vertices; }
            size_type number_of_removed_edges() const { return removed_edges; }
            size_type number_of_removed_faces() const { return removed_faces; }
            size_type number_of_removed_cells() const { return removed_cells; }

            // 多面体网格包含的几何元素(vertex, edge, face, cell)总数
            size_type number_of_used_vertices() const
            {
                return num_vertices() - number_of_removed_vertices();
            }

            size_type number_of_used_edges() const
            {
                return num_edges() - number_of_removed_edges();
            }

            size_type number_of_used_faces() const
            {
                return num_faces() - number_of_removed_faces();
            }

            size_type number_of_used_cells() const
            {
                return num_cells() - number_of_removed_cells();
            }

            bool is_removed(VertexIndex v) const { return vremoved[v]; }
            bool is_removed(EdgeIndex e) const { return eremoved[e]; }
            bool is_removed(FaceIndex f) const { return fremoved[f]; }
            bool is_removed(CellIndex c) const { return cremoved[c]; }

            void clear();

            // 是否有元素被标记为删除
            bool has_garbage() const { return garbage; }

        public: /*-------------------- 合法性检查 -------------------------*/

            bool has_valid_index(VertexIndex v) const
            {
                return ((size_type)v < num_vertices());
            }

            bool has_valid_index(EdgeIndex e) const
            {
                return ((size_type)e < num_edges());
            }

            bool has_valid_index(FaceIndex f) const
            {
                return ((size_type)f < num_faces());
            }

            bool has_valid_index(CellIndex c) const
            {
                return ((size_type)c < num_cells());
            }

            // is_valid还需要扩充检查的内容

            bool is_valid(VertexIndex v) const
            {
                if (!has_valid_index(v))
                    return false;

                return true;
            }

            bool is_valid(EdgeIndex e) const
            {
                if (!has_valid_index(e))
                    return false;

                return true;
            }

            bool is_valid(FaceIndex f) const
            {
                if (!has_valid_index(f))
                    return false;

                return true;
            }

            bool is_valid(CellIndex c) const
            {
                if (!has_valid_index(c))
                    return false;

                return true;
            }

            static VertexIndex null_vertex() { return VertexIndex(); }
            static EdgeIndex null_edge() { return EdgeIndex(); }
            static Halfedge null_halfedge() { return Halfedge(); }
            static FaceIndex null_face() { return FaceIndex(); }
            static CellIndex null_cell() { return CellIndex(); }

        private:
            /* helper functions used in topological query.
            * 将拓扑查询相关函数中公用的一些代码提取出来, 写成helper function, 简化代码.
            */

            // halfedge所在的face在cell拓扑的faces数组中的索引.
            size_type halfedge_related_face_index(Halfedge h, CellIndex c) const
            {
                assert(h.he_cell == c);
                const CellTopo& cell_topo = ctopo[c];
                size_type off = h.he_off;
                assert(off < cell_topo.offset.back());  // off不应大于offset数组的最后一个元素.

                size_type fi = 0;
                while (off >= cell_topo.offset[fi])
                    ++fi;
                return fi - 1;
            }

            // face在cell拓扑的faces数组中是否存在, 若存在返回其索引, 不存在返回invalid值.
            size_type face_location_in_cell(FaceIndex f, CellIndex c) const
            {
                const CellTopo& cell_topo = ctopo[c];
                for (size_type i = 0; i < cell_topo.faces.size(); ++i)
                {
                    if (cell_topo.faces[i] == f)
                        return i;
                }
                // 不存在返回invalid值.
                return std::numeric_limits<size_type>::max();
            }



        public: /*-------------------- 拓扑查询相关实现 ------------------------*/

            // halfedge的起始点.
            // halfedges数组中记录的<V,E>, V是对应半边的起始点而非终点.
            VertexIndex source(Halfedge h) const
            {
                CellIndex& c = h.he_cell;
                const CellTopo& cell_topo = ctopo[c];
                size_type off = h.he_off;

                bias_type vbias = cell_topo.halfedges[off].first;
                return cell_topo.vertices[vbias];
            }

            // halfedge指向的点
            VertexIndex target(Halfedge h) const
            {
                CellIndex& c = h.he_cell;
                const CellTopo& cell_topo = ctopo[c];
                size_type off = h.he_off;

                size_type fi = halfedge_related_face_index(h, c);

                // h对应的target是off所指位置的下一个, 遇到边界时为face点序列的第一个
                bias_type vbias;
                if (off + 1 == cell_topo.offset[fi + 1]) // off为face点序列的最后一个
                {
                    size_type fstart = cell_topo.offset[fi];
                    vbias = cell_topo.halfedges[fstart].first;
                }
                else
                    vbias = cell_topo.halfedges[off + 1].first;

                return cell_topo.vertices[vbias];
            }

            // halfedge所在的edge
            EdgeIndex edge(Halfedge h) const
            {
                CellIndex& c = h.he_cell;
                const CellTopo& cell_topo = ctopo[c];
                size_type off = h.he_off;

                bias_type ebias = cell_topo.halfedges[off].second;
                return cell_topo.edges[ebias];
            }

            // halfedge所在的face
            FaceIndex face(Halfedge h) const
            {
                CellIndex c = h.he_cell;
                const CellTopo& cell_topo = ctopo[c];
                size_type fi = halfedge_related_face_index(h, c);

                return cell_topo.faces[fi];
            }

            // halfedge所在的cell
            CellIndex cell(Halfedge h) const
            {
                return h.he_cell;
            }

            // halfedge的后继
            Halfedge next(Halfedge h) const
            {
                CellIndex& c = h.he_cell;
                const CellTopo& cell_topo = ctopo[c];
                size_type off = h.he_off;

                size_type fi = halfedge_related_face_index(h, c);

                if (off + 1 == cell_topo.offset[fi + 1])
                {
                    size_type fnext = cell_topo.offset[fi];
                    return Halfedge(c, fnext);
                }
                else
                    return Halfedge(c, off + 1);
            }

            // halfedge的前驱
            Halfedge prev(Halfedge h) const
            {
                CellIndex c = h.he_cell;
                const CellTopo& cell_topo = ctopo[c];
                size_type off = h.he_off;

                size_type fi = halfedge_related_face_index(h, c);

                if (off == cell_topo.offset[fi])
                {
                    size_type fprev = cell_topo.offset[fi + 1] - 1;
                    return Halfedge(c, fprev);
                }
                else
                    return Halfedge(c, off - 1);
            }

            // polygon_mate: 一个cell中共边的兄弟半边
            Halfedge polygon_mate(Halfedge h) const
            {
                CellIndex c = h.he_cell;
                const CellTopo& cell_topo = ctopo[c];
                size_type off = h.he_off;

                bias_type ebias = cell_topo.halfedges[off].second;

                for (size_type i = 0; i < cell_topo.halfedges.size(); ++i)
                {
                    if (i != off && (cell_topo.halfedges[i].second == ebias))
                        return Halfedge(c, i);
                }
                assert(false); // 执行到此处说明出现错误
            }

            // polyhedron_mate: 共面的cell中共边的兄弟半边.
            // Note: polygon_mate一定存在, 但polyhedron_mate不一定存在.
            Halfedge polyhedron_mate(Halfedge h) const
            {
                // step 1: 找到halfedge所在的面 
                FaceIndex f = face(h);

                // step 2: 进入另一个cell
                const FaceTopo& face_topo = ftopo[f];

                CellIndex c = h.he_cell;
                assert(c == face_topo.incident_cells[0] || c == face_topo.incident_cells[1]);

                CellIndex opp_c = (c == face_topo.incident_cells[0]) ?
                    face_topo.incident_cells[1] : face_topo.incident_cells[0];

                // 边界面, 只记录一个cell, 另一个为空.
                if (opp_c == null_cell())
                    return null_halfedge();

                // 找到f在opp_c的faces对应的索引.
                size_type opp_fi = face_location_in_cell(f, opp_c);
                assert(opp_fi != std::numeric_limits<size_type>::max());

                // step 3: opp_c.faces[opp_fi]包含的边中, 记录同一edge的即为所要找的结果
                const CellTopo& opp_cell_topo = ctopo[opp_c];
                size_type fstart = opp_cell_topo.offset[opp_fi];
                size_type fend = opp_cell_topo.offset[opp_fi + 1];
                EdgeIndex target_edge = edge(h);
                for (size_type i = fstart; i < fend; ++i)
                {
                    bias_type ebias = opp_cell_topo.halfedges[i].second;
                    EdgeIndex e = opp_cell_topo.edges[ebias];

                    if (e == target_edge)
                        return Halfedge(opp_c, i);
                }
                assert(false); // 执行到此处说明出现错误.
            }

            // 与点, 边, 面, 单元incident的halfedge.
            // 'source(halfedge(v)) == v', not target.
            Halfedge halfedge(VertexIndex v) const { return vtopo[v].halfedge; }
            Halfedge halfedge(EdgeIndex e) const { return etopo[e].halfedge; }
            Halfedge halfedge(FaceIndex f) const { return ftopo[f].halfedge; }
            Halfedge halfedge(CellIndex c) const { return Halfedge(c, 0); }

            Halfedge halfedge(EdgeIndex target_edge, VertexIndex srcv, CellIndex c) const
            {
                const CellTopo& cell_topo = ctopo[c];

                for (size_type i = 0; i < cell_topo.halfedges.size(); ++i)
                {
                    bias_type vbias = cell_topo.halfedges[i].first;
                    VertexIndex v = cell_topo.vertices[vbias];
                    bias_type ebias = cell_topo.halfedges[i].second;
                    EdgeIndex e = cell_topo.edges[ebias];

                    if (e == target_edge && v == srcv)
                        return Halfedge(c, i);
                }
                return null_halfedge();
            }

            CellIndex cell(FaceIndex f) const
            {
                Halfedge h = ftopo[f].halfedge;
                return h.he_cell;
            }

            // 01查询: 和一个点incident的所有边
            // 基本思想是, 对初始cell, 找到以v为顶点的那些面, 继续搜索这些面的相邻cell, 直至没有新的cell
            //
            template <typename OutputIterator>
            void incident_edges(VertexIndex v, OutputIterator out) const
            {
                if (!is_valid(v))
                    return;

                std::vector<CellIndex> crecord; // 去重
                crecord.push_back(cell(halfedge(v)));
                boost::container::flat_set<EdgeIndex> record;

                for (int k = 0; k < crecord.size(); k++) //注意, 循环中crecord的size可能会增加
                {
                    CellIndex c = crecord[k];
                    const CellTopo& cell_topo = ctopo[c];

                    // cell是watertight 流形, 以v为端点的(物理)边, 有且仅有一条以v为source的半边;
                    // 搜出cell中以v为source的所有半边 (i.e., incident faces of v)
                    for (size_type i = 0; i < cell_topo.halfedges.size(); ++i)
                    {
                        bias_type vbias = cell_topo.halfedges[i].first;
                        VertexIndex vi = cell_topo.vertices[vbias];
                        if (vi != v)
                            continue;

                        bias_type ebias = cell_topo.halfedges[i].second;
                        EdgeIndex ei = cell_topo.edges[ebias];
                        //只要cell是watertight 可定向的流形, 就不会遗漏与v incident的边;
                        record.insert(ei);

                        //当前cell的某个face的半边 source为v (该半边的"前序"半边以v为target)
                        FaceIndex f = face(Halfedge(c, i));
                        const FaceTopo& face_topo = ftopo[f];
                        assert(c == face_topo.incident_cells[0] || c == face_topo.incident_cells[1]);
                        CellIndex opp_cell = (c == face_topo.incident_cells[0]) ?
                            face_topo.incident_cells[1] : face_topo.incident_cells[0];

                        // 若opp_cell合法且之前没有遍历过, 插入.
                        if (is_valid(opp_cell) && std::find(crecord.begin(), crecord.end(), opp_cell) == crecord.end())
                            crecord.push_back(opp_cell);
                    }
                }

                for (EdgeIndex e : record)
                    *out++ = e;
            }

            // 02查询: 和一个点incident的所有面
            // 与01查询的思路类似
            template <typename OutputIterator>
            void incident_faces(VertexIndex v, OutputIterator out) const
            {
                if (!is_valid(v))
                    return;

                std::vector<CellIndex> crecord;
                crecord.push_back(cell(halfedge(v)));
                boost::container::flat_set<FaceIndex> frecord;

                for (int k = 0; k < crecord.size(); k++) //注意,循环中crecord的size可能会增加
                {
                    CellIndex c = crecord[k];
                    const CellTopo& cell_topo = ctopo[c];

                    //cell是watertight 流形, 以v为端点的(物理)边,有且仅有一条以v为source的半边;
                    //搜出cell中以v为source的所有半边 (i.e., incident faces of v)
                    for (size_type i = 0; i < cell_topo.halfedges.size(); ++i)
                    {
                        bias_type vbias = cell_topo.halfedges[i].first;
                        VertexIndex vi = cell_topo.vertices[vbias];
                        if (vi != v)
                            continue;

                        FaceIndex f = face(Halfedge(c, i));
                        frecord.insert(f);

                        const FaceTopo& face_topo = ftopo[f];
                        assert(c == face_topo.incident_cells[0] || c == face_topo.incident_cells[1]);
                        CellIndex opp_cell = (c == face_topo.incident_cells[0]) ?
                            face_topo.incident_cells[1] : face_topo.incident_cells[0];

                        // 若opp_cell合法且之前没有遍历过, 插入.
                        if (is_valid(opp_cell) && std::find(crecord.begin(), crecord.end(), opp_cell) == crecord.end())
                            crecord.push_back(opp_cell);
                    }
                }

                for (FaceIndex face : frecord)
                    *out++ = face;
            }

            // 03查询: 和一个点incident的所有单元
            template <typename OutputIterator>
            void incident_cells(VertexIndex v, OutputIterator out)
            {
                if (!is_valid(v))
                    return;

                std::vector<CellIndex> crecord;
                crecord.push_back(cell(halfedge(v)));

                for (int k = 0; k < crecord.size(); k++) //注意,循环中crecord的size可能会增加的
                {
                    CellIndex c = crecord[k];
                    const CellTopo& cell_topo = ctopo[c];

                    //cell是watertight 流形, 以v为端点的(物理)边,有且仅有一条以v为source的半边;
                    //搜出cell中以v为source的所有半边 (i.e., incident faces of v)
                    for (size_type i = 0; i < cell_topo.halfedges.size(); ++i)
                    {
                        bias_type vbias = cell_topo.halfedges[i].first;
                        VertexIndex vi = cell_topo.vertices[vbias];
                        if (vi != v)
                            continue;

                        FaceIndex f = face(Halfedge(c, i));
                        FaceTopo& face_topo = ftopo[f];
                        CellIndex opp_cell = (c == face_topo.incident_cells[0]) ?
                            face_topo.incident_cells[1] : face_topo.incident_cells[0];

                        // 若opp_cell合法且之前没有遍历过, 插入.
                        if (is_valid(opp_cell) && std::find(crecord.begin(), crecord.end(), opp_cell) == crecord.end())
                            crecord.push_back(opp_cell);
                    }
                }

                for (CellIndex c : crecord)
                    *out++ = c;
            }

            // 10查询: 和一条边incident的所有点
            template <typename OutputIterator>
            void incident_vertices(EdgeIndex e, OutputIterator out) const
            {
                if (!is_valid(e))
                    return;

                Halfedge h = halfedge(e);
                *out++ = source(h);
                *out++ = target(h);
            }

            // 12查询: 和一条边incident的所有面
            // 十字形的特殊情况目前无法处理
            template <typename OutputIterator>
            void incident_faces(EdgeIndex e, OutputIterator out) const
            {
                if (!is_valid(e))
                    return;

                Halfedge h = halfedge(e);
                Halfedge start(h);

                // 为了应对cell没有围成一圈的情况(扇形), 需要两个方向都搜一次.
                std::vector<FaceIndex> record;
                record.push_back(face(h));
                do {
                    h = polygon_mate(h);
                    record.push_back(face(h));
                    h = polyhedron_mate(h);
                } while (h != start && h != null_halfedge());

                if (h == start) //回到start, 说明不是扇形的情况, 此时record.back()==record.front(), 需要删掉一个.
                    record.pop_back();
                else // 说明搜索到达边界, 需要换个方向
                {
                    h = polyhedron_mate(start);
                    while (h != null_halfedge())
                    {
                        h = polygon_mate(h);
                        record.push_back(face(h));
                        h = polyhedron_mate(h);
                    };
                }

                for (FaceIndex f : record)
                    *out++ = f;
            }

            // 13查询: 和一条边incident的所有单元
            template <typename OutputIterator>
            void incident_cells(EdgeIndex e, OutputIterator out) const
            {
                if (!is_valid(e))
                    return;

                Halfedge h = halfedge(e);
                Halfedge start(h);

                std::vector<CellIndex> record;
                do {
                    record.push_back(cell(h));
                    h = polyhedron_mate(polygon_mate(h));
                } while (h != start && h != null_halfedge());

                if (h != start) // 说明搜索到达边界, 需要换个方向
                {
                    h = polyhedron_mate(start);
                    while (h != null_halfedge())
                    {
                        record.push_back(cell(h));
                        h = polyhedron_mate(polygon_mate(h));
                    };
                }

                for (CellIndex c : record)
                    *out++ = c;
            }

            // 20查询: 和一个面incident的所有点
            template <typename OutputIterator>
            void incident_vertices(FaceIndex f, OutputIterator out) const
            {
                if (!is_valid(f))
                    return;

                Halfedge h = halfedge(f);
                CellIndex c = cell(h);
                size_type fi = face_location_in_cell(f, c);
                assert(fi != std::numeric_limits<size_type>::max());

                const CellTopo& cell_topo = ctopo[c];

                size_type fstart = cell_topo.offset[fi];
                size_type fend = cell_topo.offset[fi + 1];

                for (size_type i = fstart; i < fend; ++i)
                {
                    bias_type vbias = cell_topo.halfedges[i].first;
                    VertexIndex v = cell_topo.vertices[vbias];
                    *out++ = v;
                }
            }

            // 21查询: 和一个面incident的所有边
            template <typename OutputIterator>
            void incident_edges(FaceIndex f, OutputIterator out) const
            {
                if (!is_valid(f))
                    return;

                Halfedge h = halfedge(f);
                CellIndex c = cell(h);
                size_type fi = face_location_in_cell(f, c);
                assert(fi != std::numeric_limits<size_type>::max());

                const CellTopo& cell_topo = ctopo[c];
                size_type fstart = cell_topo.offset[fi];
                size_type fend = cell_topo.offset[fi + 1];

                for (size_type i = fstart; i < fend; ++i)
                {
                    bias_type ebias = cell_topo.halfedges[i].second;
                    EdgeIndex e = cell_topo.edges[ebias];
                    *out++ = e;
                }
            }

            // 23查询: 和一个面incident的所有单元
            template <typename OutputIterator>
            void incident_cells(FaceIndex f, OutputIterator out) const
            {
                if (!is_valid(f))
                    return;

                for (int i = 0; i < 2; ++i)
                {
                    CellIndex c = ftopo[f].incident_cells[i];
                    if (is_valid(c))
                        *out++ = c;
                }
            }

            // 30查询: 和一个单元incident的所有点
            template <typename OutputIterator>
            void incident_vertices(CellIndex c, OutputIterator out) const
            {
                if (!is_valid(c))
                    return;

                for (VertexIndex v : ctopo[c].vertices)
                    *out++ = v;
            }

            // 31查询: 和一个单元incident的所有边
            template <typename OutputIterator>
            void incident_edges(CellIndex c, OutputIterator out) const
            {
                if (!is_valid(c))
                    return;

                for (EdgeIndex e : ctopo[c].edges)
                    *out++ = e;
            }

            // 32查询: 和一个单元incident的所有面
            template <typename OutputIterator>
            void incident_faces(CellIndex c, OutputIterator out) const
            {
                if (!is_valid(c))
                    return;

                for (FaceIndex f : ctopo[c].faces)
                    *out++ = f;
            }

            // 一个面上的一圈半边. 注意, cell不同, f对应的半边也不同.
            template <typename OutputIterator>
            void halfedges_around_face(FaceIndex f, CellIndex c, OutputIterator out) const
            {
                size_type fi = face_location_in_cell(f, c);
                assert(fi != std::numeric_limits<size_type>::max());

                const CellTopo& cell_topo = ctopo[c];
                size_type fbegin = cell_topo.offset[fi];
                size_type fend = cell_topo.offset[fi + 1];

                for (size_type i = fbegin; i < fend; ++i)
                    *out++ = Halfedge(c, i);
            }


        public: /*------------------- 拓扑修改的相关函数 ---------------------------*/

            // 判断一个面是否为三角形面.
            bool is_triangle(FaceIndex f)
            {
                CellIndex c = cell(f);
                CellTopo& cell_topo = ctopo[c];

                size_type fi = face_location_in_cell(f, c);
                assert(fi != std::numeric_limits<size_type>::max());

                size_type point_num = cell_topo.offset[fi + 1] - cell_topo.offset[fi];
                return point_num == 3;
            }

            // TODO: 删除face的时候如果相关的点和边不再使用, 应该删除, 性能更优.
            // 从cell的拓扑中删除一个面, 并维护拓扑正确.
            void remove_face_from_cell(FaceIndex f, CellIndex c)
            {
                CellTopo& cell_topo = ctopo[c];

                // 找到f的位置.
                size_type fi = face_location_in_cell(f, c);
                assert(fi != std::numeric_limits<size_type>::max());

                size_type fbegin = cell_topo.offset[fi];
                size_type fend = cell_topo.offset[fi + 1];
                size_type vertex_num = fend - fbegin;

                // faces, offset, halfedges删除相关元素, vertices和edges无需修改, 因为不会少任何点和边.
                // 1) faces数组
                cell_topo.faces.erase(cell_topo.faces.begin() + fi);
                // 2) offset数组
                // 2.1) 删除元素
                cell_topo.offset.erase(cell_topo.offset.begin() + fi);
                for (size_type i = fi; i < cell_topo.offset.size(); ++i)
                    cell_topo.offset[i] -= vertex_num;

                // 2.2) 修改face的拓扑(要删除的f不改, f之后的要改).
                for (size_type i = fi; i < cell_topo.faces.size(); ++i)
                {
                    FaceIndex f_ = cell_topo.faces[i];
                    if (cell(f_) == c)
                        ftopo[f_].halfedge = Halfedge(c, cell_topo.offset[i]);
                }

                // 3) halfedges数组
                // 3.1) 删除元素
                cell_topo.halfedges.erase(cell_topo.halfedges.begin() + fbegin, cell_topo.halfedges.begin() + fend);
                // 3.2) 修改拓扑
                // vertex, edge的拓扑中会记录halfedge, 若halfedge是本cell中的,
                // 由于halfedges删除了元素, 偏移值已经不再正确, 需要更新它们的拓扑.
                boost::container::flat_set<VertexIndex> vmodified; // 去重
                boost::container::flat_set<EdgeIndex> emodified;

                for (size_type i = 0; i < cell_topo.halfedges.size(); ++i)
                {
                    bias_type vbias = cell_topo.halfedges[i].first;
                    VertexIndex v = cell_topo.vertices[vbias];
                    Halfedge vh = halfedge(v);
                    // 仅当vh指向本cell, 索引位置在fbegin之后, 
                    // 且是第一次修改(v在halfedegs中不止一个), 才会修改.
                    if (cell(vh) == c && !vmodified.count(v))
                    {
                        vtopo[v].halfedge = Halfedge(c, i);
                        vmodified.insert(v);
                    }

                    bias_type ebias = cell_topo.halfedges[i].second;
                    EdgeIndex e = cell_topo.edges[ebias];
                    Halfedge eh = halfedge(e);
                    // 同理
                    if (cell(eh) == c && !emodified.count(e))
                    {
                        etopo[e].halfedge = Halfedge(c, i);
                        emodified.insert(e);
                    }
                }
            }

            // 将一个face加入cell拓扑中并维护拓扑正确.
            // bool变量标明是在本cell还是oppo_cell添加, 由于f_vseq对应的是本cell,
            // 所以在oppo_cell中添加应该将旋向反一下.
            void add_triangle_to_cell(FaceIndex f,
                CellIndex c,
                std::vector<VertexIndex>& f_vseq,
                std::map<SortedPair<VertexIndex>, EdgeIndex>& endpoint_to_edge,
                bool use_reverse_order)
            {
                // 要实现的效果: [0,1,2,3]->[0,3,2,1]
                if (use_reverse_order)
                    std::reverse(f_vseq.begin() + 1, f_vseq.end());

                CellTopo& cell_topo = ctopo[c];
                if (cell_topo.faces.size() == 0)
                    cell_topo.offset[0] = 0;

                // 1) faces数组
                cell_topo.faces.push_back(f);
                // 2) offset数组
                size_type fstart = cell_topo.offset.back();
                cell_topo.offset.push_back(fstart + 3); // 三角形, 固定加3即可.
                // 3) halfedges数组
                // global index与local index的映射关系, 去重.
                std::map<VertexIndex, bias_type> vlocation;
                std::map<EdgeIndex, bias_type> elocation;

                // Note: 请注意, 由于多面体网格一个cell包含的点和边的数目由bias_type的数据宽度决定,
                // 因此在添加时必须加入此assert, 以便及时报告.
                for (bias_type i = 0; i < cell_topo.vertices.size(); ++i)
                {
                    vlocation.insert(std::make_pair(cell_topo.vertices[i], i));
                    assert(cell_topo.vertices.size() < std::numeric_limits<bias_type>::max());
                }
                for (bias_type i = 0; i < cell_topo.edges.size(); ++i)
                {
                    elocation.insert(std::make_pair(cell_topo.edges[i], i));
                    assert(cell_topo.edges.size() < std::numeric_limits<bias_type>::max());
                }

                // 已有的vertex和edge的拓扑不需要调整, 新增的需要调整.
                for (std::size_t k = 0; k < 3; ++k)
                {
                    VertexIndex v = f_vseq[k];
                    bias_type vbias;
                    if (!vlocation.count(v))
                    {
                        vbias = cell_topo.vertices.size();
                        cell_topo.vertices.push_back(v);
                        vlocation.insert(std::make_pair(v, vbias));

                        // 只在本cell添加face时调整拓扑, oppo_cell添加时不调整.
                        if (!use_reverse_order)
                            vtopo[v].halfedge = Halfedge(c, fstart + k);
                    }
                    else
                        vbias = vlocation[v];

                    SortedPair<VertexIndex> vpair(f_vseq[k], f_vseq[(k + 1) % 3]);
                    assert(endpoint_to_edge.count(vpair));
                    EdgeIndex e = endpoint_to_edge[vpair];
                    bias_type ebias;
                    if (!elocation.count(e))    // 说明是新的边
                    {
                        ebias = cell_topo.edges.size();
                        cell_topo.edges.push_back(e);
                        elocation.insert(std::make_pair(e, ebias));

                        // 只在本cell添加face时调整拓扑, oppo_cell添加时不调整.
                        if (!use_reverse_order)
                            etopo[e].halfedge = Halfedge(c, fstart + k);
                    }
                    else
                        ebias = elocation[e];

                    cell_topo.halfedges.push_back(std::make_pair(vbias, ebias));
                }

                // 调整face的拓扑.
                // 只在本cell添加face时设置halfedge, oppo_cell添加时不调整.
                FaceTopo& face_topo = ftopo[f];
                if (!use_reverse_order)
                {
                    face_topo.halfedge = Halfedge(c, fstart);
                    face_topo.incident_cells[0] = c;
                }
                else
                    face_topo.incident_cells[1] = c;
            }

            // 使用分治策略将一个空间多边形面三角化.
            // Input: 面的点序列polygon, 按逆时针顺序排列.
            // Output: 三角化后的所有点序列triangles, 依然按逆时针排列.
            // 
            // 大体思路是, 对于点序列(V0 V1 ...Vi... Vn), V0和V1唯一确定一条边, 
            // 在V2~Vn中找出相对于这条边张角最大的点Vi, 
            // 此时将点序列分成了三部分：left_polygon, right_polygon, 三角形(V0, V1, Vi)
            // 三角形(V0, V1, Vi)插入triangles中, left_polygon和right_polygon继续三角化.
            //
            bool triangulation(std::vector<VertexIndex>& polygon,
                std::vector<std::vector<VertexIndex>>& triangle_faces)
            {
                std::size_t point_num = polygon.size();
                assert(point_num > 2);

                // 递归边界
                if (point_num == 3)
                {
                    triangle_faces.push_back(polygon);
                    return true;
                }

                typedef typename Kernel::FT  FT;

                // 最大张角和其对应的vertex在polygon中的索引
                FT max_angle = std::numeric_limits<FT>::min();
                std::size_t max_idx = std::numeric_limits<std::size_t>::max();

                // 求最大张角
                Point p0 = point(polygon[0]);
                Point p1 = point(polygon[1]);
                for (std::size_t i = 2; i < point_num; ++i)
                {
                    Point pi = point(polygon[i]);
                    FT angle = CGAL::approximate_angle(p0, pi, p1);  // 计算p0-pi-p1的角度, 返回的是角度值不是弧度值.
                    if (angle > max_angle)
                    {
                        max_angle = angle;
                        max_idx = i;
                    }
                }

                assert(max_idx != std::numeric_limits<std::size_t>::max() && max_idx < point_num);
                bool is_sucess = true;

                // 最大张角对应一个三角形
                triangle_faces.push_back({ polygon[0],polygon[1],polygon[max_idx] });

                // left
                if (max_idx < point_num - 1)
                {
                    std::vector<VertexIndex> left_polygon;
                    left_polygon.push_back(polygon[0]);
                    for (std::size_t i = max_idx; i < polygon.size(); ++i)
                        left_polygon.push_back(polygon[i]);

                    if (!triangulation(left_polygon, triangle_faces))
                        is_sucess = false;
                }

                // right
                if (max_idx > 2)
                {
                    std::vector<VertexIndex> right_polygon;
                    for (int i = 1; i <= max_idx; ++i)
                        right_polygon.push_back(polygon[i]);

                    if (!triangulation(right_polygon, triangle_faces))
                        is_sucess = false;
                }
                return is_sucess;
            }

            // 将一个面三角化后, 将三角化的结果导入多面体网格, 同时修改拓扑.
            // 
            // @param original_face: 被三角化的面, 需要从关联cell中删除.
            // @param f_to_vseq: 三角化后面与点序列的对应关系.
            // @param endpoint_to_edge: 物理边与端点的对应关系.
            //
            void import_triangulation_info(FaceIndex original_face,
                std::unordered_map<FaceIndex, std::vector<VertexIndex>>& f_to_vseq,
                std::map<SortedPair<VertexIndex>, EdgeIndex>& endpoint_to_edge)
            {
                CellIndex c = cell(original_face);
                FaceTopo& face_topo = ftopo[original_face];
                assert(c == face_topo.incident_cells[0] || c == face_topo.incident_cells[1]);
                CellIndex opp_c = (c == face_topo.incident_cells[0]) ?
                    face_topo.incident_cells[1] : face_topo.incident_cells[0];

                remove_face_from_cell(original_face, c);
                if (is_valid(opp_c))
                    remove_face_from_cell(original_face, opp_c);

                // Note:
                // 删除face后, vertex和edge的拓扑都已经调整正确, 但original_face的拓扑仍保留原来的值, 
                // 这些值是错误的, 但由于face没有加入cell, 所以暂时无法调整拓扑, 
                // 必须在加入面之后更改, 否则会引起拓扑错误.
                fprops.reset(static_cast<std::size_t>(original_face));

                // 将三角化后的若干三角形插入到原始面关联的两个cell的拓扑中并维护拓扑正确.
                // Note: original_face在c中的旋向保证面朝外, 在opp_c中应当把旋向反一下.
                typename std::unordered_map<FaceIndex, std::vector<VertexIndex>>::iterator it = f_to_vseq.begin();
                for (; it != f_to_vseq.end(); ++it)
                {
                    add_triangle_to_cell(it->first, c, it->second, endpoint_to_edge, false);
                    if (is_valid(opp_c))
                        add_triangle_to_cell(it->first, opp_c, it->second, endpoint_to_edge, true);
                }
            }

            // 将一个空间多边形面三角化, 在此过程中维护拓扑正确, 返回多边形三角化后的所有面.
            template <typename OutputIterator>
            void triangulate_a_polygon(FaceIndex f_split, OutputIterator out)
            {
                if (is_triangle(f_split))
                {
                    *out++ = f_split;
                    return;
                }

                std::vector<VertexIndex> f_vertices;
                incident_vertices(f_split, std::back_inserter(f_vertices));

                std::vector<std::vector<VertexIndex>> triangle_faces;

                // Step1 通过三角化, 将原始的多边形点序列分裂成若干三角形的点序列(都是逆时针)
                if (!triangulation(f_vertices, triangle_faces))
                    assert(false);

                // Step2 triangle_faces中的点序列对应一个face, 点序列相邻的两个点对应一条edge, 
                // 我们需要建立起这种对应关系, 必要时add_face, add_edge.
                std::vector<EdgeIndex> f_edges;
                incident_edges(f_split, std::back_inserter(f_edges));

                assert(f_vertices.size() == f_edges.size());

                // face与点序列的对应关系
                std::unordered_map<FaceIndex, std::vector<VertexIndex>> f_to_vseq;
                // 三角化后edge与端点的对应关系, 使用SortedPair的目的是为了让vpair唯一, 
                // 即(v0,v1)==(v1,v0), 否则add_edge()会被多调用从而造成错误.
                std::map<SortedPair<VertexIndex>, EdgeIndex> endpoint_to_edge;

                // 收集原始多边形面中边与端点的对应关系.
                std::size_t sz = f_vertices.size();
                for (std::size_t i = 0; i < sz; ++i)
                {
                    SortedPair<VertexIndex> vpair(f_vertices[i], f_vertices[(i + 1) % sz]);
                    endpoint_to_edge.insert(std::make_pair(vpair, f_edges[i]));
                }

                std::size_t f_num = triangle_faces.size();

                bool is_first = true;
                for (std::size_t i = 0; i < f_num; ++i)
                {
                    // 多边形三角化后, 要复用原有的面索引, 不会删除, 固定将triangle_faces[0]分配给f_split.
                    if (is_first)
                    {
                        is_first = false;
                        f_to_vseq.insert(std::make_pair(f_split, triangle_faces[i]));
                        *out++ = f_split; // 三角化后的面记录在传入的容器中, 方便函数后续使用.
                    }
                    else // 需要新加面.
                    {
                        FaceIndex new_face = add_face();
                        f_to_vseq.insert(std::make_pair(new_face, triangle_faces[i]));
                        *out++ = new_face;
                    }

                    std::vector<VertexIndex>& vset = triangle_faces[i];
                    assert(vset.size() == 3);
                    for (int k = 0; k < 3; ++k)
                    {
                        SortedPair<VertexIndex> vpair(vset[k], vset[(k + 1) % 3]);
                        if (!endpoint_to_edge.count(vpair)) // 非原始边, 需要新加边.
                        {
                            EdgeIndex new_edge = add_edge();
                            endpoint_to_edge.insert(std::make_pair(vpair, new_edge));
                        }
                    }
                }

                // Step3 将三角化的结果修改到多面体网格上.
                import_triangulation_info(f_split, f_to_vseq, endpoint_to_edge);
            }

            // 判断一个面记录的半边和半边h是否在一个cell中.
            bool in_same_cell(FaceIndex f, Halfedge h)
            {
                return cell(f) == cell(h);
            }

            void collect_split_faces(Halfedge h_split,
                std::vector<std::pair<FaceIndex, bool>>& record)
            {
                FaceIndex f = face(h_split);
                bool is_same_direction = in_same_cell(f, h_split);
                record.push_back(std::make_pair(f, is_same_direction));

                Halfedge anchor = h_split;
                Halfedge pos = h_split;
                do {
                    pos = polygon_mate(pos);
                    f = face(pos);
                    is_same_direction = !in_same_cell(f, pos);
                    record.push_back(std::make_pair(f, is_same_direction));
                    pos = polyhedron_mate(pos);
                } while (pos != anchor && pos != null_halfedge());

                if (pos == anchor) // 回到anchor, 说明不是扇形的情况, 此时record.back()==record.front(), 需要删掉一个.
                    record.pop_back();
                else // 说明搜索到达边界, 需要换个方向
                {
                    pos = polyhedron_mate(anchor);
                    while (pos != null_halfedge())
                    {
                        pos = polygon_mate(pos);
                        f = face(pos);
                        is_same_direction = in_same_cell(f, pos);
                        record.push_back(std::make_pair(f, is_same_direction));
                        pos = polyhedron_mate(pos);
                    }
                }
            }

            void move_face(FaceIndex f, CellIndex old_cell, CellIndex new_cell)
            {
                std::vector<VertexIndex> f_vertices;
                incident_vertices(f, std::back_inserter(f_vertices));

                std::vector<EdgeIndex> f_edges;
                incident_edges(f, std::back_inserter(f_edges));

                assert(f_vertices.size() == f_edges.size());

                std::map<SortedPair<VertexIndex>, EdgeIndex> endpoint_to_edge;
                std::size_t sz = f_vertices.size();
                for (std::size_t i = 0; i < sz; ++i)
                {
                    SortedPair<VertexIndex> vpair(f_vertices[i], f_vertices[(i + 1) % sz]);
                    endpoint_to_edge.insert(std::make_pair(vpair, f_edges[i]));
                }

                bool need_reverse = (cell(f) != old_cell);
                remove_face_from_cell(f, old_cell);
                add_triangle_to_cell(f, new_cell, f_vertices, endpoint_to_edge, need_reverse);
            }

            template <typename OutputIterator>
            void adjacent_cells(CellIndex c, OutputIterator out)
            {
                std::vector<FaceIndex> fset;
                incident_faces(c, std::back_inserter(fset));

                std::set<CellIndex> cset;
                for (FaceIndex f : fset)
                {
                    FaceTopo face_topo = ftopo[f];
                    for (int i = 0; i < 2; ++i)
                    {
                        if (is_valid(face_topo.incident_cells[i]) && face_topo.incident_cells[i] != c)
                            cset.insert(face_topo.incident_cells[i]);
                    }
                }

                for (CellIndex c : cset)
                    *out++ = c;
            }

            // 标记所有单元的内外属性
            void mark_all_cells()
            {
                std::vector<bool> marked(num_cells(), false);

                for (CellIndex c : cells())
                    if (cinfo[c].ctype != UNDEFINED)
                        marked[c] = true;

                for (CellIndex c : cells())
                {
                    if (!marked[c] || cinfo[c].ctype != IN)
                        continue;

                    std::queue<CellIndex> cell_queue;
                    cell_queue.push(c);
                    while (!cell_queue.empty())
                    {
                        CellIndex cur_cell = cell_queue.front();
                        cell_queue.pop();

                        std::vector<CellIndex> cset;
                        adjacent_cells(cur_cell, std::back_inserter(cset));

                        for (CellIndex ci : cset)
                        {
                            if (!marked[ci])
                            {
                                marked[ci] = true;
                                cinfo[ci].ctype = IN;
                                cell_queue.push(ci);
                            }
                        }
                    }
                }

                for (CellIndex c : cells())
                {
                    if (!marked[c] || cinfo[c].ctype != OUT)
                        continue;

                    std::queue<CellIndex> cell_queue;
                    cell_queue.push(c);
                    while (!cell_queue.empty())
                    {
                        CellIndex cur_cell = cell_queue.front();
                        cell_queue.pop();

                        std::vector<CellIndex> cset;
                        adjacent_cells(cur_cell, std::back_inserter(cset));

                        for (CellIndex ci : cset)
                        {
                            if (!marked[ci])
                            {
                                marked[ci] = true;
                                cinfo[ci].ctype = OUT;
                                cell_queue.push(ci);
                            }
                        }
                    }
                }
            }

        public: /*-------------------- I/O功能 -------------------------*/

            // 将截断网格的数据转化为多面体网格的形式.
            template <typename TruncatedGrid>
            bool load_from_truncated_grid(TruncatedGrid& tg)
            {
                if (!tg.exist_geom_data() || !tg.exist_ijk_info())
                    return false;

                typedef typename TruncatedGrid::index_type  index_type;
                // step 1: 创建几何元素的全局索引, 此时所有的拓扑属性都为默认值.

                // 1.1) add_vertex.
                size_type vertex_num = tg.vertices.size();
                for (size_type i = 0; i < vertex_num; ++i)
                    add_vertex(tg.vertices[i]);

                // 1.2) add_face的过程中add_edge
                // 截断网格没有边的概念, 但多面体网格的拓扑需要边的全局索引, 怎么办呢?
                // 截断网格的face_vertices_indices数组按逆时针方向记录着网格中所有面的点.
                // 由此可知, 一个面记录着多少点, 就对应多少边.
                // 但是, 多面体网格的一条边被多个面共享, 需要去重.

                typedef SortedPair<VertexIndex> Vpair;
                std::map<Vpair, EdgeIndex> endpoint_to_edge;

                size_type face_num = tg.face_vertices_indices_offset.size() - 1;
                for (size_type i = 0; i < face_num; ++i)
                {
                    // 添加面的过程中, f与i在数值上相等.
                    FaceIndex f = add_face();

                    size_type vbegin = tg.face_vertices_indices_offset[f];
                    size_type vend = tg.face_vertices_indices_offset[f + 1];

                    // 按逆时针遍历一个面的所有点时, 建立edge的global index
                    // Note: 此处的处理决定了边对应的点是起始点而非终点(source, not target)
                    for (size_type k = vbegin; k < vend; ++k)
                    {
                        Vpair vpair;

                        vpair.insert(VertexIndex(tg.face_vertices_indices[k]));
                        if (k + 1 == vend)
                            vpair.insert(VertexIndex(tg.face_vertices_indices[vbegin]));
                        else
                            vpair.insert(VertexIndex(tg.face_vertices_indices[k + 1]));

                        if (!endpoint_to_edge.count(vpair))
                        {
                            EdgeIndex e = add_edge();
                            endpoint_to_edge.insert(std::make_pair(vpair, e));
                            tg.face_edges_indices.push_back(e);
                        }
                        else
                            tg.face_edges_indices.push_back(endpoint_to_edge[vpair]);
                    }
                }
                // 二者的值完全一样, 直接拷贝即可.
                tg.face_edges_indices_offset = tg.face_vertices_indices_offset;

                // 1.3) add_cell的过程中记录面关联的cell.
                size_type cell_num = tg.cell_faces_indices_offset.size() - 1;
                int imax = -1, jmax = -1;
                for (size_type i = 0; i < cell_num; ++i)
                {
                    CellIndex c = add_cell();

                    index_type iidx = tg.iidx_per_cell[c];
                    index_type jidx = tg.jidx_per_cell[c];
                    index_type kidx = tg.kidx_per_cell[c];
                    cinfo[c] = CellInfo(iidx, jidx, kidx);
                    imax = std::max(iidx, imax);
                    jmax = std::max(jidx, jmax);

                    size_type fbegin = tg.cell_faces_indices_offset[c];
                    size_type fend = tg.cell_faces_indices_offset[c + 1];

                    for (size_type k = fbegin; k < fend; ++k)
                    {
                        FaceIndex f(tg.cell_faces_indices[k]);
                        FaceTopo& face_topo = ftopo[f];

                        if (face_topo.incident_cells[0] == null_cell())
                            face_topo.incident_cells[0] = c;
                        else
                            face_topo.incident_cells[1] = c;
                    }
                }
                IMax = imax + 1;
                JMax = jmax + 1;
                // step 2: 设置几何元素的拓扑.
                // 上面的步骤创建了所有几何元素, 相关属性只有拓扑属性仍然为空, 需要设置.

                // vertex, edge, face的拓扑中只记录着一条halfedge, 由于他们会被很多个cell共享, 
                // 所以会被多次修改, 增加该属性可以保证只修改一次拓扑
                PropertyMap<VertexIndex, bool> vmodified =
                    add_property_map<VertexIndex, bool>("v:first", false).first;
                PropertyMap<EdgeIndex, bool> emodified =
                    add_property_map<EdgeIndex, bool>("e:first", false).first;
                PropertyMap<FaceIndex, bool> fmodified =
                    add_property_map<FaceIndex, bool>("f:first", false).first;

                for (CellIndex c : cells())
                {
                    CellTopo& cell_topo = ctopo[c];

                    size_type fbegin = tg.cell_faces_indices_offset[c];
                    size_type fend = tg.cell_faces_indices_offset[c + 1];
                    size_type fnum = fend - fbegin;

                    // global index与local index的映射关系, 去重
                    std::map<VertexIndex, bias_type> vlocal;
                    std::map<EdgeIndex, bias_type> elocal;

                    size_type off(0);
                    cell_topo.offset[0] = off;
                    for (size_type i = fbegin; i < fend; ++i)
                    {
                        FaceIndex f(tg.cell_faces_indices[i]);

                        // 填充faces数组
                        cell_topo.faces.push_back(f);

                        // 填充FaceTopo
                        FaceTopo& face_topo = ftopo[f];
                        if (fmodified[f] == false)
                        {
                            Halfedge h(c, static_cast<size_type>(off));
                            face_topo.halfedge = h;
                            fmodified[f] = true;
                        }

                        // face对应的vertex与edge的范围, 二者范围一致.
                        size_type start = tg.face_vertices_indices_offset[f];
                        size_type finish = tg.face_vertices_indices_offset[f + 1];

                        std::vector<index_type> f_vertices;
                        std::vector<index_type> f_edges;

                        for (size_type k = start; k < finish; ++k)
                        {
                            f_vertices.push_back(tg.face_vertices_indices[k]);
                            f_edges.push_back(tg.face_edges_indices[k]);
                        }

                        // 要考虑cell引用face时的旋向
                        bool use_born_order = tg.faces_directions_per_cell[i];
                        if (!use_born_order)
                        {
                            // Note: 简单的反向会使得各点是对应有向边的末点, 而我们期望是起点, 
                            // 所以反转的是[begin + 1, end)而非[begin, end).
                            // 期望的效果: [0,1,2,3]-->[0,3,2,1]
                            std::reverse(f_vertices.begin() + 1, f_vertices.end());
                            std::reverse(f_edges.begin(), f_edges.end());
                        }

                        for (size_type k = 0; k < f_vertices.size(); ++k)
                        {
                            bias_type vbias;
                            VertexIndex vi(f_vertices[k]);
                            if (!vlocal.count(vi))
                            {
                                vbias = cell_topo.vertices.size();
                                cell_topo.vertices.push_back(vi);  // 填充vertices数组
                                vlocal.insert(std::make_pair(vi, vbias));
                            }
                            else
                                vbias = vlocal[vi];

                            bias_type ebias;
                            EdgeIndex ei(f_edges[k]);
                            if (!elocal.count(ei))
                            {
                                ebias = cell_topo.edges.size();
                                cell_topo.edges.push_back(ei);     // 填充edges数组
                                elocal.insert(std::make_pair(ei, ebias));
                            }
                            else
                                ebias = elocal[ei];

                            // 填充halfedges
                            cell_topo.halfedges.push_back(std::make_pair(vbias, ebias));

                            // 填充VertexTopo
                            if (vmodified[vi] == false)
                            {
                                VertexTopo& vertex_topo = vtopo[vi];

                                Halfedge h(c, static_cast<size_type>(off + k));
                                vertex_topo.halfedge = h;
                                vmodified[vi] = true;
                            }

                            // 填充EdgeTopo
                            if (emodified[ei] == false)
                            {
                                EdgeTopo& edge_topo = etopo[ei];

                                Halfedge h(c, static_cast<size_type>(off + k));
                                edge_topo.halfedge = h;
                                emodified[ei] = true;
                            }
                        }
                        off += (finish - start);
                        cell_topo.offset.push_back(off);  // 填充offset数组
                    }
                }
                // 删除附加属性.
                remove_all_attached_property_maps();
                return true;
            }


            // 将多面体网格的数据转化为截断网格的形式.
            template <typename TruncatedGrid>
            bool transform_to_truncated_grid(TruncatedGrid& tg, bool need_out)
            {
                CellType in_out_type = (need_out == true) ? OUT : IN;

                tg.startend_per_planarcell.resize(2 * IMax * JMax, 0);
                std::map<int, std::vector<CellIndex>> planar_cells;

                int in_cnt = 0, out_cnt = 0, undef_cnt = 0;
                for (CellIndex c : cells())
                {
                    if (cinfo[c].ctype == IN)
                        in_cnt++;
                    else if (cinfo[c].ctype == OUT)
                        ++out_cnt;
                    else
                        ++undef_cnt;
                }

                std::cout << "总单元数：" << num_cells() << '\n'
                    << "内部单元数" << in_cnt << '\n'
                    << "外部单元数" << out_cnt << '\n'
                    << "未定义单元数" << undef_cnt << "\n\n";

                for (CellIndex c : cells())
                {
                    CellInfo& cell_info = cinfo[c];


                    if (cell_info.ctype != in_out_type)
                        continue;

                    int i = cell_info.iidx;
                    int j = cell_info.jidx;
                    planar_cells[j * IMax + i].push_back(c);
                }

                std::vector<CellIndex> cell_set;
                for (auto it = planar_cells.begin(); it != planar_cells.end(); ++it)
                {
                    int loc = it->first;
                    tg.startend_per_planarcell[2 * loc] = cell_set.size();
                    std::vector<CellIndex>& cset = it->second;
                    for (CellIndex c : cset)
                        cell_set.push_back(c);
                    tg.startend_per_planarcell[2 * loc + 1] = cell_set.size();
                }

                tg.cell_faces_indices_offset.push_back(0);

                std::unordered_map<FaceIndex, int> face_int;
                std::vector<FaceIndex> face_set;

                for (CellIndex c : cell_set)
                {
                    CellTopo& cell_topo = ctopo[c];

                    for (FaceIndex f : cell_topo.faces)
                    {
                        if (!face_int.count(f))
                        {
                            int fi = face_set.size();
                            face_set.push_back(f);
                            face_int.insert(std::make_pair(f, fi));
                        }
                        tg.cell_faces_indices.push_back(face_int[f]);

                        if (cell(f) == c)
                            tg.faces_directions_per_cell.push_back(1);
                        else
                            tg.faces_directions_per_cell.push_back(0);
                    }

                    tg.cell_faces_indices_offset.push_back(tg.cell_faces_indices.size());
                    CellInfo& cell_info = cinfo[c];
                    tg.iidx_per_cell.push_back(cell_info.iidx);
                    tg.jidx_per_cell.push_back(cell_info.jidx);
                    tg.kidx_per_cell.push_back(cell_info.kidx);
                }

                tg.face_vertices_indices_offset.push_back(0);
                std::map<VertexIndex, std::size_t> vertex_int;
                std::vector<VertexIndex> vertex_set;

                for (FaceIndex f : face_set)
                {
                    std::vector<VertexIndex> f_vertices;
                    incident_vertices(f, std::back_inserter(f_vertices));

                    for (VertexIndex v : f_vertices)
                    {
                        if (!vertex_int.count(v))
                        {
                            int vi = vertex_set.size();
                            vertex_int.insert(std::make_pair(v, vi));
                            vertex_set.push_back(v);
                        }
                        tg.face_vertices_indices.push_back(vertex_int[v]);
                    }

                    tg.face_vertices_indices_offset.push_back(tg.face_vertices_indices.size());

                    // 面的法向量需要重新计算一下, 重新计算会导致很少几个法向量反向, 基本不影响显示.
                    Point p0 = point(f_vertices[0]);
                    Point p1 = point(f_vertices[1]);
                    for (std::size_t i = 2; i < f_vertices.size(); ++i)
                    {
                        Point pi = point(f_vertices[i]);
                        if (!Kernel::Collinear_3()(p0, p1, pi)) // 注意跳过共线的情况.
                        {
                            Normal f_norm = Kernel::Construct_normal_3()(p0, p1, pi);
                            f_norm /= CGAL::sqrt(f_norm.squared_length()); // 单位化.
                            tg.normal_per_face.push_back(f_norm);
                            break;
                        }
                    }
                }

                // 遍历点, 填充vertices数组.
                for (VertexIndex v : vertex_set)
                    tg.vertices.push_back(point(v));

                return true;
            }

            template <typename SurfaceMesh>
            void convert_cell_to_surface_mesh(CellIndex c, SurfaceMesh& sm)
            {
                assert(sm.is_empty());

                typedef typename SurfaceMesh::Vertex_index     SMVertex;
                typedef VertexIndex                            PGVertex;
                std::unordered_map<PGVertex, SMVertex> pg_vertex_to_sm_vertex;

                CellTopo& cell_topo = ctopo[c];

                for (VertexIndex v_pg : cell_topo.vertices)
                {
                    SMVertex v_sm = sm.add_vertex(point(v_pg));
                    pg_vertex_to_sm_vertex.insert(std::make_pair(v_pg, v_sm));
                }

                for (size_type i = 0; i < cell_topo.faces.size(); ++i)
                {
                    size_type fbegin = cell_topo.offset[i];
                    size_type fend = cell_topo.offset[i + 1];
                    assert(fend - fbegin == 3);

                    std::vector<SMVertex> f_vertices;
                    for (size_type k = fbegin; k < fend; ++k)
                    {
                        bias_type vbias = cell_topo.halfedges[k].first;
                        PGVertex v = cell_topo.vertices[vbias];
                        assert(pg_vertex_to_sm_vertex.count(v));
                        f_vertices.push_back(pg_vertex_to_sm_vertex[v]);
                    }
                    assert(std::set<SMVertex>(f_vertices.begin(), f_vertices.end()).size() == f_vertices.size());
                    sm.add_face(f_vertices);
                }
            }


            template <typename SurfaceMesh>
            void load_from_surface_mesh(SurfaceMesh& sm)
            {
                typedef boost::graph_traits<SurfaceMesh>              GT;
                typedef typename GT::vertex_descriptor                SMVertex;
                typedef typename GT::edge_descriptor                  SMEdge;
                typedef typename GT::halfedge_descriptor              SMHalfedge;
                typedef typename GT::face_descriptor                  SMFace;

                std::unordered_map<SMVertex, VertexIndex> sm_vertex_to_pg_vertex;

                // step 1: 创建几何元素的全局索引, 此时所有的拓扑属性都为默认值.

                // 1.1) 添加点
                for (SMVertex v_sm : sm.vertices())
                {
                    VertexIndex v_pg = add_vertex(sm.point(v_sm));
                    sm_vertex_to_pg_vertex.insert(std::make_pair(v_sm, v_pg));
                }


                // 1.2) add_face的过程中add_edge
                typedef SortedPair<VertexIndex> Vpair;
                std::map<Vpair, EdgeIndex> endpoint_to_edge;

                std::unordered_map<FaceIndex, std::vector<VertexIndex>> f_to_vseq;

                for (SMFace f_sm : sm.faces())
                {
                    FaceIndex f_pg = add_face();

                    for (SMHalfedge h_sm : sm.halfedges_around_face(sm.halfedge(f_sm)))
                    {
                        VertexIndex vsrc = sm_vertex_to_pg_vertex[sm.source(h_sm)];
                        VertexIndex vtgt = sm_vertex_to_pg_vertex[sm.target(h_sm)];

                        f_to_vseq[f_pg].push_back(vtgt);

                        Vpair vpair(vsrc, vtgt);
                        if (!endpoint_to_edge.count(vpair))
                        {
                            EdgeIndex e_pg = add_edge();
                            endpoint_to_edge.insert(std::make_pair(vpair, e_pg));
                        }
                    }
                }

                // 1.3) add_cell
                CellIndex c = add_cell();

                // step 2: 设置几何元素的拓扑.
                // 上面的步骤创建了所有几何元素, 相关属性只有拓扑属性仍然为空, 需要设置.

                {
                    CellTopo& cell_topo = ctopo[c];

                    for (auto it = f_to_vseq.begin(); it != f_to_vseq.end(); ++it)
                        add_triangle_to_cell(it->first, c, it->second, endpoint_to_edge, false);
                }
            }

            template<typename SurfaceMesh>
            void transform_to_surface_mesh(SurfaceMesh& sm)
            {
                assert(sm.is_empty());

                typedef typename SurfaceMesh::Vertex_index     SMVertex;
                typedef VertexIndex                            PGVertex;
                std::unordered_map<PGVertex, SMVertex> pg_vertex_to_sm_vertex;

                CellTopo& cell_topo = ctopo[CellIndex(0)];

                for (VertexIndex v_pg : cell_topo.vertices)
                {
                    SMVertex v_sm = sm.add_vertex(point(v_pg));
                    pg_vertex_to_sm_vertex.insert(std::make_pair(v_pg, v_sm));
                }

                for (size_type i = 0; i < cell_topo.faces.size(); ++i)
                {
                    size_type fbegin = cell_topo.offset[i];
                    size_type fend = cell_topo.offset[i + 1];
                    assert(fend - fbegin == 3);

                    std::vector<SMVertex> f_vertices;
                    for (size_type k = fbegin; k < fend; ++k)
                    {
                        bias_type vbias = cell_topo.halfedges[k].first;
                        PGVertex v = cell_topo.vertices[vbias];
                        assert(pg_vertex_to_sm_vertex.count(v));
                        f_vertices.push_back(pg_vertex_to_sm_vertex[v]);
                    }
                    assert(std::set<SMVertex>(f_vertices.begin(), f_vertices.end()).size() == f_vertices.size());
                    sm.add_face(f_vertices);
                }
            }

        private:
            /*----------------------- 属性处理模块 ----------------------------------
            *
            * 基于pg_properties.h提供的四个类, 可以很方便地管理几何元素的属性
            * 例如, vertex的一种属性由PropertyArray<T>存储, vertex的所有属性都交出一个
            * 指针, 统一放在PropertyContainer中进行管理, 如果想修改某个属性值, 则从
            * PropertyContainer中取出对应的指针, 初始化PropertyMap, 由PropertyMap进行修改
            *
            *
            * 属性分为固有属性和附加属性两类
            *     1) 固有属性: 指类自带的属性, 如v:point, v:removed, v:topology
            *     2) 附加属性: 根据实际需要添加的属性
            *
            * 本模块主要实现对属性的处理, 即属性的添加, 获取, 删除
            *
            */

            // 主要功能：选取某个几何元素的属性, PropertySelector的可见性应为private
            template <typename, bool = true>
            struct PropertySelector {};

            // 偏特化1：获取vertex的PropertyContainer
            template <bool dummy>
            struct PropertySelector<typename PolyhedronGrid<P, T>::VertexIndex, dummy>
            {
            public:
                PropertySelector(PolyhedronGrid<P, T>* pg) : p_grid(pg) {}

                // 该仿函数用于获取PropertyContainer
                PropertyContainer<Self, typename PolyhedronGrid<P, T>::VertexIndex>& // 返回值
                    operator()()
                {
                    return p_grid->vprops;
                }

                // 删除附加属性, 只保留固有属性
                // vertex的固有属性有point, topology, removed三个属性, 故为3
                void resize_property_array() { p_grid->vprops.resize_property_array(3); }

            private:
                PolyhedronGrid<P, T>* p_grid;        // 指向多面体网格的指针
            };

            // 偏特化2：获取edge的PropertyContainer
            template <bool dummy>
            struct PropertySelector<typename PolyhedronGrid<P, T>::EdgeIndex, dummy>
            {
            public:
                PropertySelector(PolyhedronGrid<P, T>* pg) : p_grid(pg) {}

                // 该仿函数用于获取PropertyContainer
                PropertyContainer<Self, typename PolyhedronGrid<P, T>::EdgeIndex>&    // 返回值
                    operator()()
                {
                    return p_grid->eprops;
                }

                // 删除附加属性, 只保留固有属性
                // edge的固有属性有topology, removed两个属性, 故为2
                void resize_property_array() { p_grid->eprops.resize_property_array(2); }

            private:
                PolyhedronGrid<P, T>* p_grid;        // 指向多面体网格的指针
            };

            // 偏特化3：获取face的PropertyContainer
            template <bool dummy>
            struct PropertySelector<typename PolyhedronGrid<P, T>::FaceIndex, dummy>
            {
            public:
                PropertySelector(PolyhedronGrid<P, T>* pg) :p_grid(pg) {}

                // 该仿函数用于获取PropertyContainer
                PropertyContainer<Self, typename PolyhedronGrid<P, T>::FaceIndex>&    // 返回值
                    operator()()
                {
                    return p_grid->fprops;
                }

                // 删除附加属性, 只保留固有属性
                void resize_property_array() { p_grid->fprops.resize_property_array(2); }

            private:
                PolyhedronGrid<P, T>* p_grid;        // 指向多面体网格的指针
            };

            // 偏特化4：获取cell的PropertyContainer
            template <bool dummy>
            struct PropertySelector<typename PolyhedronGrid<P, T>::CellIndex, dummy>
            {
            public:
                PropertySelector(PolyhedronGrid<P, T>* pg) :p_grid(pg) {}

                // 该仿函数用于获取PropertyContainer
                PropertyContainer<Self, typename PolyhedronGrid<P, T>::CellIndex>&    // 返回值
                    operator()()
                {
                    return p_grid->cprops;
                }

                // 删除附加属性, 只保留固有属性
                void resize_property_array() { p_grid->cprops.resize_property_array(3); }

            private:
                PolyhedronGrid<P, T>* p_grid;        // 指向多面体网格的指针
            };

        public:

            // 向PropertyContainer中添加属性
            // name为属性的名称, t为该属性的默认值
            // 返回值是一个pair, bool为true时表明没有该属性, 成功添加, 返回PropertyMap指向该属性数组
            //                      为false时表明该属性已存在, 无需添加
            template <typename I, typename T>
            std::pair<PropertyMap<I, T>, bool>
                add_property_map(std::string name = std::string(), const T t = T())
            {
                return PropertySelector<I>(this)().template add<T>(name, t);
            }

            // 在PropertyContainer中搜索属性
            template <typename I, typename T>
            std::pair<PropertyMap<I, T>, bool>
                get_property_map(const std::string& name) const
            {
                return PropertySelector<I>(const_cast<PolyhedronGrid*>(this))().template get<T>(name);
            }

            // 从PropertyContainer中删除某属性
            template <typename I, typename T>
            void remove_property_map(PropertyMap<I, T>& pm)
            {
                PropertySelector<I>(this)().template remove<T>(pm);
            }

            // 删除附加属性, 即由add_property_map函数添加的属性
            // 如vertex只保留'v:point','v:topology','v:removed'这三个, 其余全部删除
            template <typename I>
            void remove_attached_property_maps()
            {
                PropertySelector<I>(this).resize_property_array();
            }

            // 删除所有几何元素的附加属性
            // 注意：删除会释放属性所在的内存空间
            void remove_all_attached_property_maps()
            {
                remove_attached_property_maps<VertexIndex>();
                remove_attached_property_maps<EdgeIndex>();
                remove_attached_property_maps<FaceIndex>();
                remove_attached_property_maps<CellIndex>();
            }

            // 获取某种属性的值的类型, 即PropertyMap的value_type
            // 依靠属性名进行查询, 若属性不存在, 则会得到typeid(void)
            template <typename I>
            const std::type_info& property_type(const std::string& name)
            {
                return PropertySelector<I>(this)().get_type(name);
            }

            PropertyMap<VertexIndex, Point>&
                points()
            {
                return vpoint;
            }

            PropertyMap<VertexIndex, Point>&
                points() const
            {
                return vpoint;
            }

            Point& point(VertexIndex v) { return vpoint[v]; }
            const Point& point(VertexIndex v) const { return vpoint[v]; }

            void copy_cinfo(CellIndex old_cell, CellIndex new_cell)
            {
                cinfo[new_cell] = cinfo[old_cell];
                cinfo[new_cell].ctype = IN;
                cinfo[old_cell].ctype = OUT;
            }

            // for debug, 返回索引类型为I的所有属性名称
            template <typename I>
            std::vector<std::string> properties() const
            {
                return PropertySelector<I>(const_cast<Self*>(this))().properties();
            }

            // for debug
            void propertie_stats() const;

            // data member
        private:
            PropertyContainer<Self, VertexIndex>   vprops;    // vertex的所有属性
            PropertyContainer<Self, EdgeIndex>     eprops;    // edge的所有属性
            PropertyContainer<Self, FaceIndex>     fprops;    // face的所有属性
            PropertyContainer<Self, CellIndex>     cprops;    // cell的所有属性

            PropertyMap<VertexIndex, VertexTopo>   vtopo;     // vertex的拓扑属性
            PropertyMap<EdgeIndex, EdgeTopo>       etopo;     // edge的拓扑属性
            PropertyMap<FaceIndex, FaceTopo>       ftopo;     // face的拓扑属性
            PropertyMap<CellIndex, CellTopo>       ctopo;     // cell的拓扑属性

            PropertyMap<VertexIndex, Point>        vpoint;    // vertex的点属性
            PropertyMap<CellIndex, CellInfo>       cinfo;     // cell的IJK属性以及内外属性
            int IMax;
            int JMax;

            PropertyMap<VertexIndex, bool>      vremoved;     // vertex是否被标记为删除
            PropertyMap<EdgeIndex, bool>        eremoved;     // edge是否被标记为删除
            PropertyMap<FaceIndex, bool>        fremoved;     // face是否被标记为删除
            PropertyMap<CellIndex, bool>        cremoved;     // cell是否被标记为删除

            size_type removed_vertices;            // 标记为删除的vertex数量
            size_type removed_edges;               // 标记为删除的edge数量
            size_type removed_faces;               // 标记为删除的face数量
            size_type removed_cells;               // 标记为删除的cell数量

            size_type vertices_freelist;           // 目前可重新分配的废弃VertexIndex
            size_type edges_freelist;              // 目前可重新分配的废弃EdgeIndex
            size_type faces_freelist;              // 目前可重新分配的废弃FaceIndex
            size_type cells_freelist;              // 目前可重新分配的废弃CellIndex

            bool garbage;   // 是否有垃圾(即元素被标记为删除)
            bool recycle;   // 是否循环利用已标记为删除的index
        };


        template <typename P, typename T>
        PolyhedronGrid<P, T>::
            PolyhedronGrid()
        {
            vtopo = add_property_map<VertexIndex, VertexTopo>("v:topology").first;
            etopo = add_property_map<EdgeIndex, EdgeTopo>("e:topology").first;
            ftopo = add_property_map<FaceIndex, FaceTopo>("f:topology").first;
            ctopo = add_property_map<CellIndex, CellTopo>("c:topology").first;
            vpoint = add_property_map<VertexIndex, Point>("v:point").first;
            cinfo = add_property_map<CellIndex, CellInfo>("c:info").first;
            vremoved = add_property_map<VertexIndex, bool>("v:removed", false).first;
            eremoved = add_property_map<EdgeIndex, bool>("e:removed", false).first;
            fremoved = add_property_map<FaceIndex, bool>("f:removed", false).first;
            cremoved = add_property_map<CellIndex, bool>("c:removed", false).first;

            IMax = JMax = -1;

            removed_vertices = removed_edges = removed_faces = removed_cells = 0;
            vertices_freelist = edges_freelist = faces_freelist = cells_freelist
                = std::numeric_limits<size_type>::max();

            garbage = false;
            recycle = true;
        }

        template <typename P, typename T>
        PolyhedronGrid<P, T>&
            PolyhedronGrid<P, T>::
            operator=(const PolyhedronGrid<P, T>& rhs)
        {
            if (this != &rhs)
            {
                // 属性容器的深拷贝
                vprops = rhs.vprops;
                eprops = rhs.eprops;
                fprops = rhs.fprops;
                cprops = rhs.cprops;

                // PropertyMap存着指向属性的指针, 必须重新赋值
                // 由于PropertyContainer已更新, 所以得到的是新的pmap
                vtopo = get_property_map<VertexIndex, VertexTopo>("v:topology").first;
                etopo = get_property_map<EdgeIndex, EdgeTopo>("e:topology").first;
                ftopo = get_property_map<FaceIndex, FaceTopo>("f:topology").first;
                ctopo = get_property_map<CellIndex, CellTopo>("c:topology").first;
                vpoint = get_property_map<VertexIndex, Point>("v:point").first;
                cinfo = get_property_map<CellIndex, CellInfo>("c:info").first;
                vremoved = get_property_map<VertexIndex, bool>("v:removed", false).first;
                eremoved = get_property_map<EdgeIndex, bool>("e:removed", false).first;
                fremoved = get_property_map<FaceIndex, bool>("f:removed", false).first;
                cremoved = get_property_map<CellIndex, bool>("c:removed", false).first;

                IMax = rhs.IMax;
                JMax = rhs.JMax;

                removed_vertices = rhs.removed_vertices;
                removed_edges = rhs.removed_edges;
                removed_faces = rhs.removed_faces;
                removed_cells = rhs.removed_cells;
                vertices_freelist = rhs.vertices_freelist;
                edges_freelist = rhs.edges_freelist;
                faces_freelist = rhs.faces_freelist;
                cells_freelist = rhs.cells_freelist;
                garbage = rhs.garbage;
                recycle = rhs.recycle;
            }
            return *this;
        }

        template<typename P, typename T>
        void PolyhedronGrid<P, T>::
            clear()
        {
            // step1 清除所有属性占用的内存空间

            // 此时内存空间仍未释放
            vprops.resize(0);
            eprops.resize(0);
            fprops.resize(0);
            cprops.resize(0);

            // 释放内存空间
            vprops.shrink_to_fit();
            eprops.shrink_to_fit();
            fprops.shrink_to_fit();
            cprops.shrink_to_fit();

            //TODO 还没写完呢喂(#`O′)
            removed_vertices = removed_edges = removed_faces = removed_cells = 0;
            vertices_freelist = edges_freelist = faces_freelist = cells_freelist
                = std::numeric_limits<size_type>::max();

            garbage = false;
            recycle = true;

            // step2 清除附加的属性
            remove_all_attached_property_maps();
        }

        template <typename P, typename T>
        void PolyhedronGrid<P, T>::
            propertie_stats() const
        {
            std::vector<std::string> props;

            props = properties<VertexIndex>();
            std::cout << "Vertex has " << props.size() << " properties:\n";
            for (std::size_t i = 0; i < props.size(); ++i)
                std::cout << "\t" << props[i] << std::endl;

            props = properties<EdgeIndex>();
            std::cout << "Edge has " << props.size() << " properties:\n";
            for (std::size_t i = 0; i < props.size(); ++i)
                std::cout << "\t" << props[i] << std::endl;

            props = properties<FaceIndex>();
            std::cout << "Face has " << props.size() << " properties:\n";
            for (std::size_t i = 0; i < props.size(); ++i)
                std::cout << "\t" << props[i] << std::endl;

            props = properties<CellIndex>();
            std::cout << "Cell has " << props.size() << " properties:\n";
            for (std::size_t i = 0; i < props.size(); ++i)
                std::cout << "\t" << props[i] << std::endl;
        }


    }	// namespace MCAL::GRID
}	// namespace MCAL

#endif