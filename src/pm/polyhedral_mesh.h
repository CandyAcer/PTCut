// $Intro: ����������Ĵ���ʵ��, �����Mallison��CGAL::Surface_mesh�����˼·.
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



namespace MCAL {   // Mesh Cut Algorithm Libaray, ��Ŷ����������и��㷨�����д���
    namespace GRID {   // ��Ŷ���������ʵ�ֵ���ش���


        /*
        * �����������ʵ�ִ���, �������Թ���, ���˲�ѯ, �����޸ĵȹ���.
        *
        * @param P���������, ������double��Ԫ���cgal�ľ�ȷ����.
        * @param T��CellTopo�е�halfedgesʹ�þֲ�����ָʾ���õ�vertex��edge,
        *           �ֲ���������������, �����ݿ�Ⱦ���һ����Ԫ�������ټ���Ԫ��.
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
                int iidx, jidx, kidx; // Cell��IJK����
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

        public:// ���˶���

            // �����������halfedge��һ����Ԫ���ʾ
            // ָʾ�����ڵ�cell�Լ�CellTopo��halfedges�����е��±�
            struct Halfedge
            {
                CellIndex he_cell;  // ����cell
                size_type he_off;   // CellTopo��halfedges�����е��±�

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

                // ��Ϊmap��key, ����ʵ���ϸ�����, ��С�ڹ�ϵ.
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
                Halfedge halfedge;  // vertexΪhalfedge��source.
            };

            struct EdgeTopo
            {
                Halfedge halfedge;  // edge��Ӧ��ĳ��halfedge(�޹̶�������ϵ, ��������Ϊ1:2)
            };

            struct FaceTopo
            {
                Halfedge halfedge;
                std::array<CellIndex, 2> incident_cells;  // face������cell
            };

            struct CellTopo
            {
                // ��ʼ״̬����offset��������Ϊ��, offset��һ��Ԫ��
                CellTopo()
                {
                    offset.push_back(std::numeric_limits<size_type>::max());
                }

                std::vector<VertexIndex>      vertices;    // cell����������vertex
                std::vector<EdgeIndex>        edges;       // cell����������edge
                std::vector<FaceIndex>        faces;       // cell����������face
                std::vector<size_type>        offset;      // face��halfedges�еĿ�ʼλ��, ��faces��1
                std::vector<std::pair<bias_type, bias_type>> halfedges;
            };

        public:
            // Constructor.
            PolyhedronGrid();

            // Copy constructor: copies `rhs` to `*this`. Performs a deep copy of all properties.
            PolyhedronGrid(const PolyhedronGrid& rhs) { *this = rhs; }

            // assigns `rhs` to `*this`. Performs a deep copy of all properties.
            PolyhedronGrid& operator=(const PolyhedronGrid& rhs);

        private: // ����Ԫ�صĵ���

            template <typename Index>
            class IndexIterator
                :public boost::iterator_facade<IndexIterator<Index>,/*Derived Class*/
                Index,/*������ֵ������*/
                std::random_access_iterator_tag/*����������*/>
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
                        // ��idx�����Ϊ��ɾ��ʱ, Ӧ��+1, ֱ���ҵ���Խ����δɾ��������
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
                Index idx;                    // �����ļ���Ԫ��index
                const PolyhedronGrid* grid;   // index���ڵ�����, ���ڷ���index�Ƿ񱻱��ɾ��
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
            /*-------------------- �ڴ����ģ�� -------------------------------
            *
            * ��ģ�����Ҫ�����ǹ�����Ԫ�ص�����, ���ڴ˹����й����ڴ�
            *
            * ���㷨ִ�й�����, ��ɾ������Ӽ���Ԫ��
            * ����Ԫ�صĵ�����ͳһ����һ��vector��, ����������λ��ɾ�����������Ľϴ�
            * �ʲ������·�����
            *     1) ɾ��ĳ������Ԫ��ʱ, ����ʵ��ɾ��, ���ǽ�����Ϊremoved
            *     2) �����Ϊɾ����index�����ʱ�������·���(ѭ������)
            *
            */


            // �˴�����һ��vertices_freelist������
            // ���㷨����Ԫ�صĹ�����, ���԰�ɾ�����Ⱥ�˳��VertexIndex��������
            // ����Ԫ��ʱ, ���Խ������ɾ���ķ������Ԫ��, ��ʱ������Ա����Ϊ������ջ
            // �߼��Ͽ�, ���������һ�����Ա�, �ƺ������������������м�¼, �����Ľϴ�
            //
            //
            // ���ǲ������������ʵ��, һ���޷�������vertices_freelist����
            // ��������ջ��ָ��, ʼ�ռ�¼�����ɾ����Ԫ�ص�����
            // ֮��������Ϊfreelist, ˵���ñ�����ģ��һ���ڴ��ͷ�(free)������(list)
            //
            // �����������һ�����Թ�ϵ, �ڵ���ǰ����̹�ϵ��ô��ʾ�أ�
            // ����ע�⵽, ���Ϊɾ����VertexIndex��Ӧ����������ֵ���ǲ���ʹ�õ�
            // ��ô���Խ�VertexTopo����¼��һ�����
            //
            // ���ַ�ʽ������, ͨ�����÷������ڴ�ռ�, ֻ��һ���޷���������ʵ�ֳ���һ�����Ա�
            // ����ؼ������ڴ�ռ��
            //
            //
            // CGAL::Surface_mesh�е����ַ�ʽ���Դ�����������ڴ�����Ż���
            // 1)ɾ��ʱ�������, ������������ɾ��ʱ�ƶ��ڴ������
            // 2)vertices_freelist�ͷ����ڴ�Ľ��, ��һ�������Ϊһ������, �����ڹ���DP
            //
            // edges_freelist,faces_freelist,cells_freelistͬ��
            // freelist�ĳ�ʼֵӦΪ�Ƿ�ֵ(�����ֵ), ����Ϊ����Ŀս��(�ս��ɼ򻯲���)
            //


            VertexIndex add_vertex()
            {
                size_type inf = std::numeric_limits<size_type>::max();

                // recycleΪtrue��ʾ�������������ѱ��Ϊɾ����Index
                // vertices_freelist != inf ��ʾɾ����Ԫ��, ���Խ��з���
                if (recycle && (vertices_freelist != inf))
                {
                    VertexIndex v(vertices_freelist);   // ���������ɾ��������
                    vertices_freelist = vtopo[v].halfedge.he_off;   //����freelist
                    --removed_vertices;
                    vremoved[v] = false;   // ����removed���
                    vprops.reset(v);       // ����v��Ӧ����������ֵ, �������޸�
                    return v;
                }
                else    // ���ɷ���ʱpush_back����(index+1), �ϼ�
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

            // ɾ��Vertex, �����Ϊɾ��, ��������������������ڴ���
            void remove_vertex(VertexIndex v)
            {
                vremoved[v] = true; ++removed_vertices; garbage = true;

                // VertexTopo��¼��һ��halfedge, Halfedge���������ݳ�Ա
                // ������he_off����¼vertices_freelist(��Ϊ���Ͷ���size_type)
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
                // CellTopo�в���¼halfedge, ���Բ�����ǰ��������һ��
                // CellTopo��Ĭ�Ϲ��캯��������Ķ�����offset��һ��Ԫ��, ����Ϊ��
                // ����ʹ�����Ԫ�ؼ�¼
                ctopo[c].offset[0] = cells_freelist;
                cells_freelist = (size_type)c;
            }

            // The number of used and removed vertices in the gird.
            // �����������м���Ԫ��(vertex, edge, face, cell)������(�����ѱ��ɾ����)
            // ʵ�����Ѿ������index����
            size_type num_vertices() const { return (size_type)vprops.size(); }
            size_type num_edges() const { return (size_type)eprops.size(); }
            size_type num_faces() const { return (size_type)fprops.size(); }
            size_type num_cells() const { return (size_type)cprops.size(); }

            // �����������б����Ϊ��ɾ���ļ���Ԫ��(vertex, edge, face, cell)����
            size_type number_of_removed_vertices() const { return removed_vertices; }
            size_type number_of_removed_edges() const { return removed_edges; }
            size_type number_of_removed_faces() const { return removed_faces; }
            size_type number_of_removed_cells() const { return removed_cells; }

            // ��������������ļ���Ԫ��(vertex, edge, face, cell)����
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

            // �Ƿ���Ԫ�ر����Ϊɾ��
            bool has_garbage() const { return garbage; }

        public: /*-------------------- �Ϸ��Լ�� -------------------------*/

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

            // is_valid����Ҫ�����������

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
            * �����˲�ѯ��غ����й��õ�һЩ������ȡ����, д��helper function, �򻯴���.
            */

            // halfedge���ڵ�face��cell���˵�faces�����е�����.
            size_type halfedge_related_face_index(Halfedge h, CellIndex c) const
            {
                assert(h.he_cell == c);
                const CellTopo& cell_topo = ctopo[c];
                size_type off = h.he_off;
                assert(off < cell_topo.offset.back());  // off��Ӧ����offset��������һ��Ԫ��.

                size_type fi = 0;
                while (off >= cell_topo.offset[fi])
                    ++fi;
                return fi - 1;
            }

            // face��cell���˵�faces�������Ƿ����, �����ڷ���������, �����ڷ���invalidֵ.
            size_type face_location_in_cell(FaceIndex f, CellIndex c) const
            {
                const CellTopo& cell_topo = ctopo[c];
                for (size_type i = 0; i < cell_topo.faces.size(); ++i)
                {
                    if (cell_topo.faces[i] == f)
                        return i;
                }
                // �����ڷ���invalidֵ.
                return std::numeric_limits<size_type>::max();
            }



        public: /*-------------------- ���˲�ѯ���ʵ�� ------------------------*/

            // halfedge����ʼ��.
            // halfedges�����м�¼��<V,E>, V�Ƕ�Ӧ��ߵ���ʼ������յ�.
            VertexIndex source(Halfedge h) const
            {
                CellIndex& c = h.he_cell;
                const CellTopo& cell_topo = ctopo[c];
                size_type off = h.he_off;

                bias_type vbias = cell_topo.halfedges[off].first;
                return cell_topo.vertices[vbias];
            }

            // halfedgeָ��ĵ�
            VertexIndex target(Halfedge h) const
            {
                CellIndex& c = h.he_cell;
                const CellTopo& cell_topo = ctopo[c];
                size_type off = h.he_off;

                size_type fi = halfedge_related_face_index(h, c);

                // h��Ӧ��target��off��ָλ�õ���һ��, �����߽�ʱΪface�����еĵ�һ��
                bias_type vbias;
                if (off + 1 == cell_topo.offset[fi + 1]) // offΪface�����е����һ��
                {
                    size_type fstart = cell_topo.offset[fi];
                    vbias = cell_topo.halfedges[fstart].first;
                }
                else
                    vbias = cell_topo.halfedges[off + 1].first;

                return cell_topo.vertices[vbias];
            }

            // halfedge���ڵ�edge
            EdgeIndex edge(Halfedge h) const
            {
                CellIndex& c = h.he_cell;
                const CellTopo& cell_topo = ctopo[c];
                size_type off = h.he_off;

                bias_type ebias = cell_topo.halfedges[off].second;
                return cell_topo.edges[ebias];
            }

            // halfedge���ڵ�face
            FaceIndex face(Halfedge h) const
            {
                CellIndex c = h.he_cell;
                const CellTopo& cell_topo = ctopo[c];
                size_type fi = halfedge_related_face_index(h, c);

                return cell_topo.faces[fi];
            }

            // halfedge���ڵ�cell
            CellIndex cell(Halfedge h) const
            {
                return h.he_cell;
            }

            // halfedge�ĺ��
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

            // halfedge��ǰ��
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

            // polygon_mate: һ��cell�й��ߵ��ֵܰ��
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
                assert(false); // ִ�е��˴�˵�����ִ���
            }

            // polyhedron_mate: �����cell�й��ߵ��ֵܰ��.
            // Note: polygon_mateһ������, ��polyhedron_mate��һ������.
            Halfedge polyhedron_mate(Halfedge h) const
            {
                // step 1: �ҵ�halfedge���ڵ��� 
                FaceIndex f = face(h);

                // step 2: ������һ��cell
                const FaceTopo& face_topo = ftopo[f];

                CellIndex c = h.he_cell;
                assert(c == face_topo.incident_cells[0] || c == face_topo.incident_cells[1]);

                CellIndex opp_c = (c == face_topo.incident_cells[0]) ?
                    face_topo.incident_cells[1] : face_topo.incident_cells[0];

                // �߽���, ֻ��¼һ��cell, ��һ��Ϊ��.
                if (opp_c == null_cell())
                    return null_halfedge();

                // �ҵ�f��opp_c��faces��Ӧ������.
                size_type opp_fi = face_location_in_cell(f, opp_c);
                assert(opp_fi != std::numeric_limits<size_type>::max());

                // step 3: opp_c.faces[opp_fi]�����ı���, ��¼ͬһedge�ļ�Ϊ��Ҫ�ҵĽ��
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
                assert(false); // ִ�е��˴�˵�����ִ���.
            }

            // ���, ��, ��, ��Ԫincident��halfedge.
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

            // 01��ѯ: ��һ����incident�����б�
            // ����˼����, �Գ�ʼcell, �ҵ���vΪ�������Щ��, ����������Щ�������cell, ֱ��û���µ�cell
            //
            template <typename OutputIterator>
            void incident_edges(VertexIndex v, OutputIterator out) const
            {
                if (!is_valid(v))
                    return;

                std::vector<CellIndex> crecord; // ȥ��
                crecord.push_back(cell(halfedge(v)));
                boost::container::flat_set<EdgeIndex> record;

                for (int k = 0; k < crecord.size(); k++) //ע��, ѭ����crecord��size���ܻ�����
                {
                    CellIndex c = crecord[k];
                    const CellTopo& cell_topo = ctopo[c];

                    // cell��watertight ����, ��vΪ�˵��(����)��, ���ҽ���һ����vΪsource�İ��;
                    // �ѳ�cell����vΪsource�����а�� (i.e., incident faces of v)
                    for (size_type i = 0; i < cell_topo.halfedges.size(); ++i)
                    {
                        bias_type vbias = cell_topo.halfedges[i].first;
                        VertexIndex vi = cell_topo.vertices[vbias];
                        if (vi != v)
                            continue;

                        bias_type ebias = cell_topo.halfedges[i].second;
                        EdgeIndex ei = cell_topo.edges[ebias];
                        //ֻҪcell��watertight �ɶ��������, �Ͳ�����©��v incident�ı�;
                        record.insert(ei);

                        //��ǰcell��ĳ��face�İ�� sourceΪv (�ð�ߵ�"ǰ��"�����vΪtarget)
                        FaceIndex f = face(Halfedge(c, i));
                        const FaceTopo& face_topo = ftopo[f];
                        assert(c == face_topo.incident_cells[0] || c == face_topo.incident_cells[1]);
                        CellIndex opp_cell = (c == face_topo.incident_cells[0]) ?
                            face_topo.incident_cells[1] : face_topo.incident_cells[0];

                        // ��opp_cell�Ϸ���֮ǰû�б�����, ����.
                        if (is_valid(opp_cell) && std::find(crecord.begin(), crecord.end(), opp_cell) == crecord.end())
                            crecord.push_back(opp_cell);
                    }
                }

                for (EdgeIndex e : record)
                    *out++ = e;
            }

            // 02��ѯ: ��һ����incident��������
            // ��01��ѯ��˼·����
            template <typename OutputIterator>
            void incident_faces(VertexIndex v, OutputIterator out) const
            {
                if (!is_valid(v))
                    return;

                std::vector<CellIndex> crecord;
                crecord.push_back(cell(halfedge(v)));
                boost::container::flat_set<FaceIndex> frecord;

                for (int k = 0; k < crecord.size(); k++) //ע��,ѭ����crecord��size���ܻ�����
                {
                    CellIndex c = crecord[k];
                    const CellTopo& cell_topo = ctopo[c];

                    //cell��watertight ����, ��vΪ�˵��(����)��,���ҽ���һ����vΪsource�İ��;
                    //�ѳ�cell����vΪsource�����а�� (i.e., incident faces of v)
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

                        // ��opp_cell�Ϸ���֮ǰû�б�����, ����.
                        if (is_valid(opp_cell) && std::find(crecord.begin(), crecord.end(), opp_cell) == crecord.end())
                            crecord.push_back(opp_cell);
                    }
                }

                for (FaceIndex face : frecord)
                    *out++ = face;
            }

            // 03��ѯ: ��һ����incident�����е�Ԫ
            template <typename OutputIterator>
            void incident_cells(VertexIndex v, OutputIterator out)
            {
                if (!is_valid(v))
                    return;

                std::vector<CellIndex> crecord;
                crecord.push_back(cell(halfedge(v)));

                for (int k = 0; k < crecord.size(); k++) //ע��,ѭ����crecord��size���ܻ����ӵ�
                {
                    CellIndex c = crecord[k];
                    const CellTopo& cell_topo = ctopo[c];

                    //cell��watertight ����, ��vΪ�˵��(����)��,���ҽ���һ����vΪsource�İ��;
                    //�ѳ�cell����vΪsource�����а�� (i.e., incident faces of v)
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

                        // ��opp_cell�Ϸ���֮ǰû�б�����, ����.
                        if (is_valid(opp_cell) && std::find(crecord.begin(), crecord.end(), opp_cell) == crecord.end())
                            crecord.push_back(opp_cell);
                    }
                }

                for (CellIndex c : crecord)
                    *out++ = c;
            }

            // 10��ѯ: ��һ����incident�����е�
            template <typename OutputIterator>
            void incident_vertices(EdgeIndex e, OutputIterator out) const
            {
                if (!is_valid(e))
                    return;

                Halfedge h = halfedge(e);
                *out++ = source(h);
                *out++ = target(h);
            }

            // 12��ѯ: ��һ����incident��������
            // ʮ���ε��������Ŀǰ�޷�����
            template <typename OutputIterator>
            void incident_faces(EdgeIndex e, OutputIterator out) const
            {
                if (!is_valid(e))
                    return;

                Halfedge h = halfedge(e);
                Halfedge start(h);

                // Ϊ��Ӧ��cellû��Χ��һȦ�����(����), ��Ҫ����������һ��.
                std::vector<FaceIndex> record;
                record.push_back(face(h));
                do {
                    h = polygon_mate(h);
                    record.push_back(face(h));
                    h = polyhedron_mate(h);
                } while (h != start && h != null_halfedge());

                if (h == start) //�ص�start, ˵���������ε����, ��ʱrecord.back()==record.front(), ��Ҫɾ��һ��.
                    record.pop_back();
                else // ˵����������߽�, ��Ҫ��������
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

            // 13��ѯ: ��һ����incident�����е�Ԫ
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

                if (h != start) // ˵����������߽�, ��Ҫ��������
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

            // 20��ѯ: ��һ����incident�����е�
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

            // 21��ѯ: ��һ����incident�����б�
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

            // 23��ѯ: ��һ����incident�����е�Ԫ
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

            // 30��ѯ: ��һ����Ԫincident�����е�
            template <typename OutputIterator>
            void incident_vertices(CellIndex c, OutputIterator out) const
            {
                if (!is_valid(c))
                    return;

                for (VertexIndex v : ctopo[c].vertices)
                    *out++ = v;
            }

            // 31��ѯ: ��һ����Ԫincident�����б�
            template <typename OutputIterator>
            void incident_edges(CellIndex c, OutputIterator out) const
            {
                if (!is_valid(c))
                    return;

                for (EdgeIndex e : ctopo[c].edges)
                    *out++ = e;
            }

            // 32��ѯ: ��һ����Ԫincident��������
            template <typename OutputIterator>
            void incident_faces(CellIndex c, OutputIterator out) const
            {
                if (!is_valid(c))
                    return;

                for (FaceIndex f : ctopo[c].faces)
                    *out++ = f;
            }

            // һ�����ϵ�һȦ���. ע��, cell��ͬ, f��Ӧ�İ��Ҳ��ͬ.
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


        public: /*------------------- �����޸ĵ���غ��� ---------------------------*/

            // �ж�һ�����Ƿ�Ϊ��������.
            bool is_triangle(FaceIndex f)
            {
                CellIndex c = cell(f);
                CellTopo& cell_topo = ctopo[c];

                size_type fi = face_location_in_cell(f, c);
                assert(fi != std::numeric_limits<size_type>::max());

                size_type point_num = cell_topo.offset[fi + 1] - cell_topo.offset[fi];
                return point_num == 3;
            }

            // TODO: ɾ��face��ʱ�������صĵ�ͱ߲���ʹ��, Ӧ��ɾ��, ���ܸ���.
            // ��cell��������ɾ��һ����, ��ά��������ȷ.
            void remove_face_from_cell(FaceIndex f, CellIndex c)
            {
                CellTopo& cell_topo = ctopo[c];

                // �ҵ�f��λ��.
                size_type fi = face_location_in_cell(f, c);
                assert(fi != std::numeric_limits<size_type>::max());

                size_type fbegin = cell_topo.offset[fi];
                size_type fend = cell_topo.offset[fi + 1];
                size_type vertex_num = fend - fbegin;

                // faces, offset, halfedgesɾ�����Ԫ��, vertices��edges�����޸�, ��Ϊ�������κε�ͱ�.
                // 1) faces����
                cell_topo.faces.erase(cell_topo.faces.begin() + fi);
                // 2) offset����
                // 2.1) ɾ��Ԫ��
                cell_topo.offset.erase(cell_topo.offset.begin() + fi);
                for (size_type i = fi; i < cell_topo.offset.size(); ++i)
                    cell_topo.offset[i] -= vertex_num;

                // 2.2) �޸�face������(Ҫɾ����f����, f֮���Ҫ��).
                for (size_type i = fi; i < cell_topo.faces.size(); ++i)
                {
                    FaceIndex f_ = cell_topo.faces[i];
                    if (cell(f_) == c)
                        ftopo[f_].halfedge = Halfedge(c, cell_topo.offset[i]);
                }

                // 3) halfedges����
                // 3.1) ɾ��Ԫ��
                cell_topo.halfedges.erase(cell_topo.halfedges.begin() + fbegin, cell_topo.halfedges.begin() + fend);
                // 3.2) �޸�����
                // vertex, edge�������л��¼halfedge, ��halfedge�Ǳ�cell�е�,
                // ����halfedgesɾ����Ԫ��, ƫ��ֵ�Ѿ�������ȷ, ��Ҫ�������ǵ�����.
                boost::container::flat_set<VertexIndex> vmodified; // ȥ��
                boost::container::flat_set<EdgeIndex> emodified;

                for (size_type i = 0; i < cell_topo.halfedges.size(); ++i)
                {
                    bias_type vbias = cell_topo.halfedges[i].first;
                    VertexIndex v = cell_topo.vertices[vbias];
                    Halfedge vh = halfedge(v);
                    // ����vhָ��cell, ����λ����fbegin֮��, 
                    // ���ǵ�һ���޸�(v��halfedegs�в�ֹһ��), �Ż��޸�.
                    if (cell(vh) == c && !vmodified.count(v))
                    {
                        vtopo[v].halfedge = Halfedge(c, i);
                        vmodified.insert(v);
                    }

                    bias_type ebias = cell_topo.halfedges[i].second;
                    EdgeIndex e = cell_topo.edges[ebias];
                    Halfedge eh = halfedge(e);
                    // ͬ��
                    if (cell(eh) == c && !emodified.count(e))
                    {
                        etopo[e].halfedge = Halfedge(c, i);
                        emodified.insert(e);
                    }
                }
            }

            // ��һ��face����cell�����в�ά��������ȷ.
            // bool�����������ڱ�cell����oppo_cell���, ����f_vseq��Ӧ���Ǳ�cell,
            // ������oppo_cell�����Ӧ�ý�����һ��.
            void add_triangle_to_cell(FaceIndex f,
                CellIndex c,
                std::vector<VertexIndex>& f_vseq,
                std::map<SortedPair<VertexIndex>, EdgeIndex>& endpoint_to_edge,
                bool use_reverse_order)
            {
                // Ҫʵ�ֵ�Ч��: [0,1,2,3]->[0,3,2,1]
                if (use_reverse_order)
                    std::reverse(f_vseq.begin() + 1, f_vseq.end());

                CellTopo& cell_topo = ctopo[c];
                if (cell_topo.faces.size() == 0)
                    cell_topo.offset[0] = 0;

                // 1) faces����
                cell_topo.faces.push_back(f);
                // 2) offset����
                size_type fstart = cell_topo.offset.back();
                cell_topo.offset.push_back(fstart + 3); // ������, �̶���3����.
                // 3) halfedges����
                // global index��local index��ӳ���ϵ, ȥ��.
                std::map<VertexIndex, bias_type> vlocation;
                std::map<EdgeIndex, bias_type> elocation;

                // Note: ��ע��, ���ڶ���������һ��cell�����ĵ�ͱߵ���Ŀ��bias_type�����ݿ�Ⱦ���,
                // ��������ʱ��������assert, �Ա㼰ʱ����.
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

                // ���е�vertex��edge�����˲���Ҫ����, ��������Ҫ����.
                for (std::size_t k = 0; k < 3; ++k)
                {
                    VertexIndex v = f_vseq[k];
                    bias_type vbias;
                    if (!vlocation.count(v))
                    {
                        vbias = cell_topo.vertices.size();
                        cell_topo.vertices.push_back(v);
                        vlocation.insert(std::make_pair(v, vbias));

                        // ֻ�ڱ�cell���faceʱ��������, oppo_cell���ʱ������.
                        if (!use_reverse_order)
                            vtopo[v].halfedge = Halfedge(c, fstart + k);
                    }
                    else
                        vbias = vlocation[v];

                    SortedPair<VertexIndex> vpair(f_vseq[k], f_vseq[(k + 1) % 3]);
                    assert(endpoint_to_edge.count(vpair));
                    EdgeIndex e = endpoint_to_edge[vpair];
                    bias_type ebias;
                    if (!elocation.count(e))    // ˵�����µı�
                    {
                        ebias = cell_topo.edges.size();
                        cell_topo.edges.push_back(e);
                        elocation.insert(std::make_pair(e, ebias));

                        // ֻ�ڱ�cell���faceʱ��������, oppo_cell���ʱ������.
                        if (!use_reverse_order)
                            etopo[e].halfedge = Halfedge(c, fstart + k);
                    }
                    else
                        ebias = elocation[e];

                    cell_topo.halfedges.push_back(std::make_pair(vbias, ebias));
                }

                // ����face������.
                // ֻ�ڱ�cell���faceʱ����halfedge, oppo_cell���ʱ������.
                FaceTopo& face_topo = ftopo[f];
                if (!use_reverse_order)
                {
                    face_topo.halfedge = Halfedge(c, fstart);
                    face_topo.incident_cells[0] = c;
                }
                else
                    face_topo.incident_cells[1] = c;
            }

            // ʹ�÷��β��Խ�һ���ռ����������ǻ�.
            // Input: ��ĵ�����polygon, ����ʱ��˳������.
            // Output: ���ǻ�������е�����triangles, ��Ȼ����ʱ������.
            // 
            // ����˼·��, ���ڵ�����(V0 V1 ...Vi... Vn), V0��V1Ψһȷ��һ����, 
            // ��V2~Vn���ҳ�������������Ž����ĵ�Vi, 
            // ��ʱ�������зֳ��������֣�left_polygon, right_polygon, ������(V0, V1, Vi)
            // ������(V0, V1, Vi)����triangles��, left_polygon��right_polygon�������ǻ�.
            //
            bool triangulation(std::vector<VertexIndex>& polygon,
                std::vector<std::vector<VertexIndex>>& triangle_faces)
            {
                std::size_t point_num = polygon.size();
                assert(point_num > 2);

                // �ݹ�߽�
                if (point_num == 3)
                {
                    triangle_faces.push_back(polygon);
                    return true;
                }

                typedef typename Kernel::FT  FT;

                // ����ŽǺ����Ӧ��vertex��polygon�е�����
                FT max_angle = std::numeric_limits<FT>::min();
                std::size_t max_idx = std::numeric_limits<std::size_t>::max();

                // ������Ž�
                Point p0 = point(polygon[0]);
                Point p1 = point(polygon[1]);
                for (std::size_t i = 2; i < point_num; ++i)
                {
                    Point pi = point(polygon[i]);
                    FT angle = CGAL::approximate_angle(p0, pi, p1);  // ����p0-pi-p1�ĽǶ�, ���ص��ǽǶ�ֵ���ǻ���ֵ.
                    if (angle > max_angle)
                    {
                        max_angle = angle;
                        max_idx = i;
                    }
                }

                assert(max_idx != std::numeric_limits<std::size_t>::max() && max_idx < point_num);
                bool is_sucess = true;

                // ����ŽǶ�Ӧһ��������
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

            // ��һ�������ǻ���, �����ǻ��Ľ���������������, ͬʱ�޸�����.
            // 
            // @param original_face: �����ǻ�����, ��Ҫ�ӹ���cell��ɾ��.
            // @param f_to_vseq: ���ǻ�����������еĶ�Ӧ��ϵ.
            // @param endpoint_to_edge: �������˵�Ķ�Ӧ��ϵ.
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
                // ɾ��face��, vertex��edge�����˶��Ѿ�������ȷ, ��original_face�������Ա���ԭ����ֵ, 
                // ��Щֵ�Ǵ����, ������faceû�м���cell, ������ʱ�޷���������, 
                // �����ڼ�����֮�����, ������������˴���.
                fprops.reset(static_cast<std::size_t>(original_face));

                // �����ǻ�������������β��뵽ԭʼ�����������cell�������в�ά��������ȷ.
                // Note: original_face��c�е�����֤�泯��, ��opp_c��Ӧ��������һ��.
                typename std::unordered_map<FaceIndex, std::vector<VertexIndex>>::iterator it = f_to_vseq.begin();
                for (; it != f_to_vseq.end(); ++it)
                {
                    add_triangle_to_cell(it->first, c, it->second, endpoint_to_edge, false);
                    if (is_valid(opp_c))
                        add_triangle_to_cell(it->first, opp_c, it->second, endpoint_to_edge, true);
                }
            }

            // ��һ���ռ����������ǻ�, �ڴ˹�����ά��������ȷ, ���ض�������ǻ����������.
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

                // Step1 ͨ�����ǻ�, ��ԭʼ�Ķ���ε����з��ѳ����������εĵ�����(������ʱ��)
                if (!triangulation(f_vertices, triangle_faces))
                    assert(false);

                // Step2 triangle_faces�еĵ����ж�Ӧһ��face, ���������ڵ��������Ӧһ��edge, 
                // ������Ҫ���������ֶ�Ӧ��ϵ, ��Ҫʱadd_face, add_edge.
                std::vector<EdgeIndex> f_edges;
                incident_edges(f_split, std::back_inserter(f_edges));

                assert(f_vertices.size() == f_edges.size());

                // face������еĶ�Ӧ��ϵ
                std::unordered_map<FaceIndex, std::vector<VertexIndex>> f_to_vseq;
                // ���ǻ���edge��˵�Ķ�Ӧ��ϵ, ʹ��SortedPair��Ŀ����Ϊ����vpairΨһ, 
                // ��(v0,v1)==(v1,v0), ����add_edge()�ᱻ����ôӶ���ɴ���.
                std::map<SortedPair<VertexIndex>, EdgeIndex> endpoint_to_edge;

                // �ռ�ԭʼ��������б���˵�Ķ�Ӧ��ϵ.
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
                    // ��������ǻ���, Ҫ����ԭ�е�������, ����ɾ��, �̶���triangle_faces[0]�����f_split.
                    if (is_first)
                    {
                        is_first = false;
                        f_to_vseq.insert(std::make_pair(f_split, triangle_faces[i]));
                        *out++ = f_split; // ���ǻ�������¼�ڴ����������, ���㺯������ʹ��.
                    }
                    else // ��Ҫ�¼���.
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
                        if (!endpoint_to_edge.count(vpair)) // ��ԭʼ��, ��Ҫ�¼ӱ�.
                        {
                            EdgeIndex new_edge = add_edge();
                            endpoint_to_edge.insert(std::make_pair(vpair, new_edge));
                        }
                    }
                }

                // Step3 �����ǻ��Ľ���޸ĵ�������������.
                import_triangulation_info(f_split, f_to_vseq, endpoint_to_edge);
            }

            // �ж�һ�����¼�İ�ߺͰ��h�Ƿ���һ��cell��.
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

                if (pos == anchor) // �ص�anchor, ˵���������ε����, ��ʱrecord.back()==record.front(), ��Ҫɾ��һ��.
                    record.pop_back();
                else // ˵����������߽�, ��Ҫ��������
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

            // ������е�Ԫ����������
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

        public: /*-------------------- I/O���� -------------------------*/

            // ���ض����������ת��Ϊ�������������ʽ.
            template <typename TruncatedGrid>
            bool load_from_truncated_grid(TruncatedGrid& tg)
            {
                if (!tg.exist_geom_data() || !tg.exist_ijk_info())
                    return false;

                typedef typename TruncatedGrid::index_type  index_type;
                // step 1: ��������Ԫ�ص�ȫ������, ��ʱ���е��������Զ�ΪĬ��ֵ.

                // 1.1) add_vertex.
                size_type vertex_num = tg.vertices.size();
                for (size_type i = 0; i < vertex_num; ++i)
                    add_vertex(tg.vertices[i]);

                // 1.2) add_face�Ĺ�����add_edge
                // �ض�����û�бߵĸ���, �������������������Ҫ�ߵ�ȫ������, ��ô����?
                // �ض������face_vertices_indices���鰴��ʱ�뷽���¼��������������ĵ�.
                // �ɴ˿�֪, һ�����¼�Ŷ��ٵ�, �Ͷ�Ӧ���ٱ�.
                // ����, �����������һ���߱�����湲��, ��Ҫȥ��.

                typedef SortedPair<VertexIndex> Vpair;
                std::map<Vpair, EdgeIndex> endpoint_to_edge;

                size_type face_num = tg.face_vertices_indices_offset.size() - 1;
                for (size_type i = 0; i < face_num; ++i)
                {
                    // �����Ĺ�����, f��i����ֵ�����.
                    FaceIndex f = add_face();

                    size_type vbegin = tg.face_vertices_indices_offset[f];
                    size_type vend = tg.face_vertices_indices_offset[f + 1];

                    // ����ʱ�����һ��������е�ʱ, ����edge��global index
                    // Note: �˴��Ĵ�������˱߶�Ӧ�ĵ�����ʼ������յ�(source, not target)
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
                // ���ߵ�ֵ��ȫһ��, ֱ�ӿ�������.
                tg.face_edges_indices_offset = tg.face_vertices_indices_offset;

                // 1.3) add_cell�Ĺ����м�¼�������cell.
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
                // step 2: ���ü���Ԫ�ص�����.
                // ����Ĳ��贴�������м���Ԫ��, �������ֻ������������ȻΪ��, ��Ҫ����.

                // vertex, edge, face��������ֻ��¼��һ��halfedge, �������ǻᱻ�ܶ��cell����, 
                // ���Իᱻ����޸�, ���Ӹ����Կ��Ա�ֻ֤�޸�һ������
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

                    // global index��local index��ӳ���ϵ, ȥ��
                    std::map<VertexIndex, bias_type> vlocal;
                    std::map<EdgeIndex, bias_type> elocal;

                    size_type off(0);
                    cell_topo.offset[0] = off;
                    for (size_type i = fbegin; i < fend; ++i)
                    {
                        FaceIndex f(tg.cell_faces_indices[i]);

                        // ���faces����
                        cell_topo.faces.push_back(f);

                        // ���FaceTopo
                        FaceTopo& face_topo = ftopo[f];
                        if (fmodified[f] == false)
                        {
                            Halfedge h(c, static_cast<size_type>(off));
                            face_topo.halfedge = h;
                            fmodified[f] = true;
                        }

                        // face��Ӧ��vertex��edge�ķ�Χ, ���߷�Χһ��.
                        size_type start = tg.face_vertices_indices_offset[f];
                        size_type finish = tg.face_vertices_indices_offset[f + 1];

                        std::vector<index_type> f_vertices;
                        std::vector<index_type> f_edges;

                        for (size_type k = start; k < finish; ++k)
                        {
                            f_vertices.push_back(tg.face_vertices_indices[k]);
                            f_edges.push_back(tg.face_edges_indices[k]);
                        }

                        // Ҫ����cell����faceʱ������
                        bool use_born_order = tg.faces_directions_per_cell[i];
                        if (!use_born_order)
                        {
                            // Note: �򵥵ķ����ʹ�ø����Ƕ�Ӧ����ߵ�ĩ��, ���������������, 
                            // ���Է�ת����[begin + 1, end)����[begin, end).
                            // ������Ч��: [0,1,2,3]-->[0,3,2,1]
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
                                cell_topo.vertices.push_back(vi);  // ���vertices����
                                vlocal.insert(std::make_pair(vi, vbias));
                            }
                            else
                                vbias = vlocal[vi];

                            bias_type ebias;
                            EdgeIndex ei(f_edges[k]);
                            if (!elocal.count(ei))
                            {
                                ebias = cell_topo.edges.size();
                                cell_topo.edges.push_back(ei);     // ���edges����
                                elocal.insert(std::make_pair(ei, ebias));
                            }
                            else
                                ebias = elocal[ei];

                            // ���halfedges
                            cell_topo.halfedges.push_back(std::make_pair(vbias, ebias));

                            // ���VertexTopo
                            if (vmodified[vi] == false)
                            {
                                VertexTopo& vertex_topo = vtopo[vi];

                                Halfedge h(c, static_cast<size_type>(off + k));
                                vertex_topo.halfedge = h;
                                vmodified[vi] = true;
                            }

                            // ���EdgeTopo
                            if (emodified[ei] == false)
                            {
                                EdgeTopo& edge_topo = etopo[ei];

                                Halfedge h(c, static_cast<size_type>(off + k));
                                edge_topo.halfedge = h;
                                emodified[ei] = true;
                            }
                        }
                        off += (finish - start);
                        cell_topo.offset.push_back(off);  // ���offset����
                    }
                }
                // ɾ����������.
                remove_all_attached_property_maps();
                return true;
            }


            // �����������������ת��Ϊ�ض��������ʽ.
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

                std::cout << "�ܵ�Ԫ����" << num_cells() << '\n'
                    << "�ڲ���Ԫ��" << in_cnt << '\n'
                    << "�ⲿ��Ԫ��" << out_cnt << '\n'
                    << "δ���嵥Ԫ��" << undef_cnt << "\n\n";

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

                    // ��ķ�������Ҫ���¼���һ��, ���¼���ᵼ�º��ټ�������������, ������Ӱ����ʾ.
                    Point p0 = point(f_vertices[0]);
                    Point p1 = point(f_vertices[1]);
                    for (std::size_t i = 2; i < f_vertices.size(); ++i)
                    {
                        Point pi = point(f_vertices[i]);
                        if (!Kernel::Collinear_3()(p0, p1, pi)) // ע���������ߵ����.
                        {
                            Normal f_norm = Kernel::Construct_normal_3()(p0, p1, pi);
                            f_norm /= CGAL::sqrt(f_norm.squared_length()); // ��λ��.
                            tg.normal_per_face.push_back(f_norm);
                            break;
                        }
                    }
                }

                // ������, ���vertices����.
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

                // step 1: ��������Ԫ�ص�ȫ������, ��ʱ���е��������Զ�ΪĬ��ֵ.

                // 1.1) ��ӵ�
                for (SMVertex v_sm : sm.vertices())
                {
                    VertexIndex v_pg = add_vertex(sm.point(v_sm));
                    sm_vertex_to_pg_vertex.insert(std::make_pair(v_sm, v_pg));
                }


                // 1.2) add_face�Ĺ�����add_edge
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

                // step 2: ���ü���Ԫ�ص�����.
                // ����Ĳ��贴�������м���Ԫ��, �������ֻ������������ȻΪ��, ��Ҫ����.

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
            /*----------------------- ���Դ���ģ�� ----------------------------------
            *
            * ����pg_properties.h�ṩ���ĸ���, ���Ժܷ���ع�����Ԫ�ص�����
            * ����, vertex��һ��������PropertyArray<T>�洢, vertex���������Զ�����һ��
            * ָ��, ͳһ����PropertyContainer�н��й���, ������޸�ĳ������ֵ, ���
            * PropertyContainer��ȡ����Ӧ��ָ��, ��ʼ��PropertyMap, ��PropertyMap�����޸�
            *
            *
            * ���Է�Ϊ�������Ժ͸�����������
            *     1) ��������: ָ���Դ�������, ��v:point, v:removed, v:topology
            *     2) ��������: ����ʵ����Ҫ��ӵ�����
            *
            * ��ģ����Ҫʵ�ֶ����ԵĴ���, �����Ե����, ��ȡ, ɾ��
            *
            */

            // ��Ҫ���ܣ�ѡȡĳ������Ԫ�ص�����, PropertySelector�Ŀɼ���ӦΪprivate
            template <typename, bool = true>
            struct PropertySelector {};

            // ƫ�ػ�1����ȡvertex��PropertyContainer
            template <bool dummy>
            struct PropertySelector<typename PolyhedronGrid<P, T>::VertexIndex, dummy>
            {
            public:
                PropertySelector(PolyhedronGrid<P, T>* pg) : p_grid(pg) {}

                // �÷º������ڻ�ȡPropertyContainer
                PropertyContainer<Self, typename PolyhedronGrid<P, T>::VertexIndex>& // ����ֵ
                    operator()()
                {
                    return p_grid->vprops;
                }

                // ɾ����������, ֻ������������
                // vertex�Ĺ���������point, topology, removed��������, ��Ϊ3
                void resize_property_array() { p_grid->vprops.resize_property_array(3); }

            private:
                PolyhedronGrid<P, T>* p_grid;        // ָ������������ָ��
            };

            // ƫ�ػ�2����ȡedge��PropertyContainer
            template <bool dummy>
            struct PropertySelector<typename PolyhedronGrid<P, T>::EdgeIndex, dummy>
            {
            public:
                PropertySelector(PolyhedronGrid<P, T>* pg) : p_grid(pg) {}

                // �÷º������ڻ�ȡPropertyContainer
                PropertyContainer<Self, typename PolyhedronGrid<P, T>::EdgeIndex>&    // ����ֵ
                    operator()()
                {
                    return p_grid->eprops;
                }

                // ɾ����������, ֻ������������
                // edge�Ĺ���������topology, removed��������, ��Ϊ2
                void resize_property_array() { p_grid->eprops.resize_property_array(2); }

            private:
                PolyhedronGrid<P, T>* p_grid;        // ָ������������ָ��
            };

            // ƫ�ػ�3����ȡface��PropertyContainer
            template <bool dummy>
            struct PropertySelector<typename PolyhedronGrid<P, T>::FaceIndex, dummy>
            {
            public:
                PropertySelector(PolyhedronGrid<P, T>* pg) :p_grid(pg) {}

                // �÷º������ڻ�ȡPropertyContainer
                PropertyContainer<Self, typename PolyhedronGrid<P, T>::FaceIndex>&    // ����ֵ
                    operator()()
                {
                    return p_grid->fprops;
                }

                // ɾ����������, ֻ������������
                void resize_property_array() { p_grid->fprops.resize_property_array(2); }

            private:
                PolyhedronGrid<P, T>* p_grid;        // ָ������������ָ��
            };

            // ƫ�ػ�4����ȡcell��PropertyContainer
            template <bool dummy>
            struct PropertySelector<typename PolyhedronGrid<P, T>::CellIndex, dummy>
            {
            public:
                PropertySelector(PolyhedronGrid<P, T>* pg) :p_grid(pg) {}

                // �÷º������ڻ�ȡPropertyContainer
                PropertyContainer<Self, typename PolyhedronGrid<P, T>::CellIndex>&    // ����ֵ
                    operator()()
                {
                    return p_grid->cprops;
                }

                // ɾ����������, ֻ������������
                void resize_property_array() { p_grid->cprops.resize_property_array(3); }

            private:
                PolyhedronGrid<P, T>* p_grid;        // ָ������������ָ��
            };

        public:

            // ��PropertyContainer���������
            // nameΪ���Ե�����, tΪ�����Ե�Ĭ��ֵ
            // ����ֵ��һ��pair, boolΪtrueʱ����û�и�����, �ɹ����, ����PropertyMapָ�����������
            //                      Ϊfalseʱ�����������Ѵ���, �������
            template <typename I, typename T>
            std::pair<PropertyMap<I, T>, bool>
                add_property_map(std::string name = std::string(), const T t = T())
            {
                return PropertySelector<I>(this)().template add<T>(name, t);
            }

            // ��PropertyContainer����������
            template <typename I, typename T>
            std::pair<PropertyMap<I, T>, bool>
                get_property_map(const std::string& name) const
            {
                return PropertySelector<I>(const_cast<PolyhedronGrid*>(this))().template get<T>(name);
            }

            // ��PropertyContainer��ɾ��ĳ����
            template <typename I, typename T>
            void remove_property_map(PropertyMap<I, T>& pm)
            {
                PropertySelector<I>(this)().template remove<T>(pm);
            }

            // ɾ����������, ����add_property_map������ӵ�����
            // ��vertexֻ����'v:point','v:topology','v:removed'������, ����ȫ��ɾ��
            template <typename I>
            void remove_attached_property_maps()
            {
                PropertySelector<I>(this).resize_property_array();
            }

            // ɾ�����м���Ԫ�صĸ�������
            // ע�⣺ɾ�����ͷ��������ڵ��ڴ�ռ�
            void remove_all_attached_property_maps()
            {
                remove_attached_property_maps<VertexIndex>();
                remove_attached_property_maps<EdgeIndex>();
                remove_attached_property_maps<FaceIndex>();
                remove_attached_property_maps<CellIndex>();
            }

            // ��ȡĳ�����Ե�ֵ������, ��PropertyMap��value_type
            // �������������в�ѯ, �����Բ�����, ���õ�typeid(void)
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

            // for debug, ������������ΪI��������������
            template <typename I>
            std::vector<std::string> properties() const
            {
                return PropertySelector<I>(const_cast<Self*>(this))().properties();
            }

            // for debug
            void propertie_stats() const;

            // data member
        private:
            PropertyContainer<Self, VertexIndex>   vprops;    // vertex����������
            PropertyContainer<Self, EdgeIndex>     eprops;    // edge����������
            PropertyContainer<Self, FaceIndex>     fprops;    // face����������
            PropertyContainer<Self, CellIndex>     cprops;    // cell����������

            PropertyMap<VertexIndex, VertexTopo>   vtopo;     // vertex����������
            PropertyMap<EdgeIndex, EdgeTopo>       etopo;     // edge����������
            PropertyMap<FaceIndex, FaceTopo>       ftopo;     // face����������
            PropertyMap<CellIndex, CellTopo>       ctopo;     // cell����������

            PropertyMap<VertexIndex, Point>        vpoint;    // vertex�ĵ�����
            PropertyMap<CellIndex, CellInfo>       cinfo;     // cell��IJK�����Լ���������
            int IMax;
            int JMax;

            PropertyMap<VertexIndex, bool>      vremoved;     // vertex�Ƿ񱻱��Ϊɾ��
            PropertyMap<EdgeIndex, bool>        eremoved;     // edge�Ƿ񱻱��Ϊɾ��
            PropertyMap<FaceIndex, bool>        fremoved;     // face�Ƿ񱻱��Ϊɾ��
            PropertyMap<CellIndex, bool>        cremoved;     // cell�Ƿ񱻱��Ϊɾ��

            size_type removed_vertices;            // ���Ϊɾ����vertex����
            size_type removed_edges;               // ���Ϊɾ����edge����
            size_type removed_faces;               // ���Ϊɾ����face����
            size_type removed_cells;               // ���Ϊɾ����cell����

            size_type vertices_freelist;           // Ŀǰ�����·���ķ���VertexIndex
            size_type edges_freelist;              // Ŀǰ�����·���ķ���EdgeIndex
            size_type faces_freelist;              // Ŀǰ�����·���ķ���FaceIndex
            size_type cells_freelist;              // Ŀǰ�����·���ķ���CellIndex

            bool garbage;   // �Ƿ�������(��Ԫ�ر����Ϊɾ��)
            bool recycle;   // �Ƿ�ѭ�������ѱ��Ϊɾ����index
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
                // �������������
                vprops = rhs.vprops;
                eprops = rhs.eprops;
                fprops = rhs.fprops;
                cprops = rhs.cprops;

                // PropertyMap����ָ�����Ե�ָ��, �������¸�ֵ
                // ����PropertyContainer�Ѹ���, ���Եõ������µ�pmap
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
            // step1 �����������ռ�õ��ڴ�ռ�

            // ��ʱ�ڴ�ռ���δ�ͷ�
            vprops.resize(0);
            eprops.resize(0);
            fprops.resize(0);
            cprops.resize(0);

            // �ͷ��ڴ�ռ�
            vprops.shrink_to_fit();
            eprops.shrink_to_fit();
            fprops.shrink_to_fit();
            cprops.shrink_to_fit();

            //TODO ��ûд����ι(#`O��)
            removed_vertices = removed_edges = removed_faces = removed_cells = 0;
            vertices_freelist = edges_freelist = faces_freelist = cells_freelist
                = std::numeric_limits<size_type>::max();

            garbage = false;
            recycle = true;

            // step2 ������ӵ�����
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