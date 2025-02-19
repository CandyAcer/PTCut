// $License: GPL-3.0
// $Author: Liu Zhiqiang(from Beihang University)
// $Email: lzhq0930@gmail.com
// 
// Part 2 of PTCut, add intersection points and segments into surface and 
// polyhedral mesh, triangulate mesh to form a well-defined topology.
// 


#ifndef MCAL_ALGO_COREFINE_H
#define MCAL_ALGO_COREFINE_H

#include <CGAL/utility.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <boost/dynamic_bitset/dynamic_bitset.hpp>

#include <vector>
#include <map>
#include <unordered_map>

#include "intersection_nodes.h"
#include "predicates.h"

namespace MCAL   // Mesh Cut Algorithm Libaray, ��Ŷ����������и��㷨��ʵ��
{

/*
* Corefine���𽫱�������Ͷ�����������ཻ�������ǻ�, �������ǻ����������������.
* ���ǻ�������Ϣ������Intersection, �����Ϊ�������ݳ�Ա(����˵Visitor), �����㷨�������Ѽ������Ϣ.
*/
template <typename SurfaceMesh,
		  typename PolyhedralMesh,
		  typename VPMSM,
		  typename VPMPM,
		  typename OutputBuilder>
class Corefine
{
	// typedefs
	typedef boost::graph_traits<SurfaceMesh>             GT;
	typedef typename GT::vertex_descriptor               SMVertex;
	typedef typename GT::edge_descriptor                 SMEdge;
	typedef typename GT::halfedge_descriptor             SMHalfedge;
	typedef typename GT::face_descriptor                 SMFace;

	typedef typename PolyhedralMesh::VertexIndex         PMVertex;
	typedef typename PolyhedralMesh::EdgeIndex           PMEdge;
	typedef typename PolyhedralMesh::FaceIndex           PMFace;
	typedef typename PolyhedralMesh::CellIndex           PMCell;
	typedef typename PolyhedralMesh::Halfedge            PMHalfedge;

	typedef std::size_t                                  NodeId;
	typedef std::pair<NodeId, NodeId>                    NodeIdPair;
	typedef std::vector<NodeId>                          NodeIdVec;

	typedef std::unordered_map<SMFace, NodeIdVec>        OnSMFace;
	typedef std::unordered_map<SMEdge, NodeIdVec>        OnSMEdge;
	typedef std::vector<SMVertex>                        NodeToSMVertex;
	typedef std::unordered_map<SMVertex, NodeId>         SMVertexToNode;
	typedef std::multimap<NodeId, SMHalfedge>            NodeToSMTarget;

	typedef std::unordered_map<PMFace, NodeIdVec>        OnPMFace;
	typedef std::unordered_map<PMEdge, NodeIdVec>        OnPMEdge;
	typedef std::vector<PMVertex>                        NodeToPMVertex;
	typedef std::unordered_map<PMVertex, NodeId>         PMVertexToNode;
	typedef std::multimap<NodeId, PMHalfedge>            NodeToPMTarget;

	typedef IntersectionNodes<SurfaceMesh, PolyhedralMesh, VPMSM, VPMPM>    INodes;
	typedef typename INodes::ExactKernel                                    EK;

	typedef CGAL::Projection_traits_3<EK>                                   CDTtraits;
	typedef CGAL::Triangulation_vertex_base_with_info_2<NodeId, CDTtraits>  Vb;
	typedef CGAL::Constrained_triangulation_face_base_2<CDTtraits>          Fb;
	typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                    TDS2;
	typedef CGAL::Constrained_Delaunay_triangulation_2<CDTtraits, TDS2>     CDT;
	typedef typename CDT::Vertex_handle                                     CDT_Vertex_handle;

	// data members
private:
	// ������ϸ����Ҫ��������Ϣ:

	// 1) corresponding information. ���㵽������ʲôλ��, ע��sm��pm��������һ��.

	OnSMFace on_sm_face;  // ����sm���ڲ��Ľ���
	OnSMEdge on_sm_edge;  // ����sm���ϵĽ���
	NodeToSMVertex node_id_to_sm_vertex;   // node id��sm vertex�Ķ�Ӧ��ϵ
	SMVertexToNode sm_vertex_to_node_id;   // sm vertex��node id�Ķ�Ӧ��ϵ
	NodeToSMTarget node_to_sm_target_halfedge;  // ������ָ�����İ��

	OnPMFace on_pm_face;  // ����pm���ڲ��Ľ���
	OnPMEdge on_pm_edge;  // ����pm���ϵĽ���
	NodeToPMVertex node_id_to_pm_vertex;   // node id��pm vertex�Ķ�Ӧ��ϵ
	PMVertexToNode pm_vertex_to_node_id;   // pm vertex��node id�Ķ�Ӧ��ϵ
	NodeToPMTarget node_to_pm_target_halfedge;  // ������ָ�����İ��

	// 2) adjacent information. ����֮����ڽӹ�ϵ.

	// �洢һ����������ڽ���, һ���������ڽ����Ӧһ������, ���������߹��ɽ���.
	std::vector<std::vector<NodeId>> graph_of_constraints;
	// ��ʶһ�������Ƿ�ֻ��һ���ھ�(���߲��պ�ʱ�Ķ˵�).
	boost::dynamic_bitset<> is_node_of_degree_one;

	OutputBuilder& output_builder;

public:
	Corefine(OutputBuilder& ob) :output_builder(ob) {}

	// ��surface mesh edge �� polyhedral mesh face�ཻʱ, �����ཻ����, ��¼��Ӧ��ϵ.
	void new_node_added_sepf(NodeId node_id,
							 IntersectionType type,
							 SMHalfedge h_sm,
							 PMHalfedge h_pm,
							 SurfaceMesh& sm,
							 PolyhedralMesh& pm,
							 bool is_source_coplanar,
							 bool is_target_coplanar)
	{
		// Step1 ��¼node id�����������face, edge, vertex�Ķ�Ӧ��ϵ.
		// Note: ����������/�����Ϊkey���м�¼��.
		switch (type)
		{
		case ON_FACE: // ��sm edge����pm face���ڲ�ʱ, Ҫ��¼��on_pm_face��.
		{
			on_pm_face[pm.face(h_pm)].push_back(node_id);
		}
		break;

		case ON_EDGE: // ��sm edge����pm face��ĳ������ʱ, Ҫ��¼��on_pm_edge��.
		{
			on_pm_edge[pm.edge(h_pm)].push_back(node_id);
		}
		break;

		case ON_VERTEX: // ��sm edge����pm face��ĳ��������ʱ, Ҫ��¼vertex��node id�Ķ�Ӧ��ϵ.
		{
			node_to_pm_target_halfedge.insert(std::make_pair(node_id, h_pm));
			if (node_id_to_pm_vertex.size() <= node_id)
				node_id_to_pm_vertex.resize(node_id + 1, pm.null_vertex());
			node_id_to_pm_vertex[node_id] = pm.target(h_pm);
			bool insert_ok =
				pm_vertex_to_node_id.insert(std::make_pair(pm.target(h_pm), node_id)).second;
			assert(insert_ok || pm_vertex_to_node_id[pm.target(h_pm)] == node_id);

			// ע��, ��ʱ�����Ӧ��vertex�ڶ����������л�û��, Ҫ�����ǻ��Ĺ������𲽴������ɸú�����¼.
			output_builder.set_pm_vertex_id(pm.target(h_pm), node_id);
		}
		break;

		default:
			return;
		}

		// Step2 ��¼node id���������face, edge, vertex�Ķ�Ӧ��ϵ.
		assert(!is_source_coplanar || !is_target_coplanar);

		if (is_target_coplanar)
		{
			node_to_sm_target_halfedge.insert(std::make_pair(node_id, h_sm));
			if (node_id_to_sm_vertex.size() <= node_id)
				node_id_to_sm_vertex.resize(node_id + 1, GT::null_vertex());
			node_id_to_sm_vertex[node_id] = target(h_sm, sm);
			bool insert_ok =
				sm_vertex_to_node_id.insert(std::make_pair(target(h_sm, sm), node_id)).second;
			assert(insert_ok || sm_vertex_to_node_id[target(h_sm, sm)] == node_id);

			output_builder.set_sm_vertex_id(target(h_sm, sm), node_id);
		}
		else if (is_source_coplanar)
		{
			SMHalfedge opp_h_sm = opposite(h_sm, sm);
			node_to_sm_target_halfedge.insert(std::make_pair(node_id, opp_h_sm));
			if (node_id_to_sm_vertex.size() <= node_id)
				node_id_to_sm_vertex.resize(node_id + 1, GT::null_vertex());
			node_id_to_sm_vertex[node_id] = target(opp_h_sm, sm);
			bool insert_ok =
				sm_vertex_to_node_id.insert(std::make_pair(target(opp_h_sm, sm), node_id)).second;
			assert(insert_ok || sm_vertex_to_node_id[target(opp_h_sm, sm)] == node_id);

			output_builder.set_sm_vertex_id(source(h_sm, sm), node_id);
		}
		else // ��������sm edge��.
			on_sm_edge[edge(h_sm, sm)].push_back(node_id);
	}

	// ��polyhedral mesh edge �� surface mesh face�ཻʱ, �����ཻ����, ��¼��Ӧ��ϵ.
	void new_node_added_pesf(NodeId node_id,
							 IntersectionType type,
							 PMHalfedge h_pm,
							 SMHalfedge h_sm,
							 PolyhedralMesh& pm,
							 SurfaceMesh& sm,
							 bool is_source_coplanar,
							 bool is_target_coplanar)
	{
		switch (type)
		{
		case ON_FACE:
		{
			on_sm_face[face(h_sm, sm)].push_back(node_id);
		}
		break;

		case ON_EDGE:
		{
			on_sm_edge[edge(h_sm, sm)].push_back(node_id);
		}
		break;

		case ON_VERTEX:
		{
			node_to_sm_target_halfedge.insert(std::make_pair(node_id, h_sm));
			if (node_id_to_sm_vertex.size() <= node_id)
				node_id_to_sm_vertex.resize(node_id + 1, GT::null_vertex());
			node_id_to_sm_vertex[node_id] = target(h_sm, sm);
			bool insert_ok =
				sm_vertex_to_node_id.insert(std::make_pair(target(h_sm, sm), node_id)).second;
			assert(insert_ok || sm_vertex_to_node_id[target(h_sm, sm)] == node_id);

			output_builder.set_sm_vertex_id(target(h_sm, sm), node_id);
		}
		break;

		default:
			return;
		}

		assert(!is_source_coplanar || !is_target_coplanar);

		if (is_target_coplanar)
		{
			node_to_pm_target_halfedge.insert(std::make_pair(node_id, h_pm));
			if (node_id_to_pm_vertex.size() <= node_id)
				node_id_to_pm_vertex.resize(node_id + 1, pm.null_vertex());
			node_id_to_pm_vertex[node_id] = pm.target(h_pm);
			bool insert_ok =
				pm_vertex_to_node_id.insert(std::make_pair(pm.target(h_pm), node_id)).second;
			assert(insert_ok || pm_vertex_to_node_id[pm.target(h_pm)] == node_id);

			output_builder.set_pm_vertex_id(pm.target(h_pm), node_id);
		}
		else if (is_source_coplanar)
		{
			PMHalfedge opp_h_pm = pm.polygon_mate(h_pm);
			node_to_pm_target_halfedge.insert(std::make_pair(node_id, opp_h_pm));
			if (node_id_to_pm_vertex.size() <= node_id)
				node_id_to_pm_vertex.resize(node_id + 1, pm.null_vertex());
			node_id_to_pm_vertex[node_id] = pm.target(opp_h_pm);
			bool insert_ok =
				pm_vertex_to_node_id.insert(std::make_pair(pm.target(opp_h_pm), node_id)).second;
			assert(insert_ok || pm_vertex_to_node_id[pm.target(opp_h_pm)] == node_id);

			output_builder.set_pm_vertex_id(pm.source(h_pm), node_id);
		}
		else
			on_pm_edge[pm.edge(h_pm)].push_back(node_id);
	}

	// �ռ��ڽ���Ϣ, �򵥸��Ƽ���.
	template <typename Neighbors>
	void annotate_graph(std::vector<Neighbors>& node_neighbors)
	{
		std::size_t nb_nodes = node_neighbors.size();
		graph_of_constraints.resize(nb_nodes);
		is_node_of_degree_one.resize(nb_nodes);

		for (NodeId node_id = 0; node_id < nb_nodes; ++node_id)
		{
			graph_of_constraints[node_id].assign(
				node_neighbors[node_id].begin(),
				node_neighbors[node_id].end());

			if (graph_of_constraints[node_id].size() == 1)
				is_node_of_degree_one.set(node_id);
		}
	}

	void record_segment_cell(SortedPair<NodeId>& segment, PMCell cell)
	{
		output_builder.record_segment_cell(segment, cell);
	}

	void start_new_polyline(PMCell c, NodeId i, NodeId j)
	{
		assert(i != j);
		output_builder.start_new_polyline(c, i, j);
	}

	void add_node_to_polyline(PMCell c, NodeId i)
	{
		output_builder.add_node_to_polyline(c, i);
	}

	// һ��sm face�ı߽���Ϣ, ������������, ��������Լ�ÿ������ϵĽ���.
	// split_halfedges()��ʱ����Щ��Ϣ�ռ�����, ��on_sm_faceһ������CDT, ʵ�ֶ�һ��������ǻ�.
	// SMFaceBoundary�Ǳ߽���Ϣ, on_sm_face�����ڲ�����Ϣ, graph_of_constraints�ǽ�����Ϣ,
	// ���߼�¼�����ǻ�һ���������ȫ����Ϣ.
	struct SMFaceBoundary
	{
		SMVertex                    vertices[3]; // the three vertices of the original face.
		SMHalfedge                  halfedges[3]; //the three halfedges of the original face.
		std::map<SMHalfedge, int>   hedge_id_map; // map a halfedge to an index.
		std::vector<NodeId>	        ids_on_edge[3]; // the node_ids on each halfedge.

		// ��Ӧ��ϵ
		// (splitǰ):
		// vertices[0,1]-->halfedges[0]-->ids_on_edge[0]
		// vertices[1,2]-->halfedges[1]-->ids_on_edge[1]
		// vertices[2,0]-->halfedges[2]-->ids_on_edge[2]
		// 
		// (split��):
		// halfedges[0]��target��vertices[1], ��source����vertices[0]
		// halfedges[1]��target��vertices[2], ��source����vertices[1]
		// halfedges[2]��target��vertices[0], ��source����vertices[2]
		//

		SMFaceBoundary(SMHalfedge first, SurfaceMesh& sm)
		{
			assert(is_triangle(first, sm));
			halfedges[0] = first;
			halfedges[1] = next(first, sm);
			halfedges[2] = next(halfedges[1], sm);

			vertices[0] = source(halfedges[0], sm);
			vertices[1] = source(halfedges[1], sm);
			vertices[2] = source(halfedges[2], sm);

			hedge_id_map.insert(std::make_pair(halfedges[0], 0));
			hedge_id_map.insert(std::make_pair(halfedges[1], 1));
			hedge_id_map.insert(std::make_pair(halfedges[2], 2));
		}

		// ��halfedge��Ӧ��SMFaceBoundary������halfedge���ֵܰ�߱�split��ʱ��, ��Ҫ���°��.
		// ���º�halfedge�����, ��ids_on_edge�ϵĽ�����Ȼ��δ�з�ǰ�ĳ���, ��һ�ֲ���Ӧ�ĸо�.
		void update_original_halfedge(SMHalfedge original_halfedge,
									  SMHalfedge new_halfedge,
									  SurfaceMesh& sm)
		{
			typename std::map<SMHalfedge, int>::iterator it
				= hedge_id_map.find(original_halfedge);
			assert(it != hedge_id_map.end());
			int idx = it->second;
			assert(halfedges[idx] == original_halfedge);
			hedge_id_map.erase(it);
			hedge_id_map.insert(std::make_pair(new_halfedge, idx));
			halfedges[idx] = new_halfedge;
		}

		// ��on_sm_edge�Ľ��㿽������Ӧ�İ����.
		template <typename Iterator>
		void copy_ids_on_edge(SMHalfedge h, Iterator begin, Iterator end)
		{
			typename std::map<SMHalfedge, int>::iterator it
				= hedge_id_map.find(h);
			assert(it != hedge_id_map.end());
			std::copy(begin, end, std::back_inserter(ids_on_edge[it->second]));
		}
	};
	typedef std::unordered_map<SMFace, SMFaceBoundary> SMFaceBoundaries;

	// ��h_split�ϵĵ��������, ����Զ��������.
	// �������ids_on_edge������ԽС�ĵ���target(h_split, sm)ԽԶ.
	void sort_vertices_along_sm_halfedge(NodeIdVec& ids_on_edge,
										 SMHalfedge h_split,
										 SurfaceMesh& sm,
										 VPMSM& vpm_sm,
										 INodes& nodes)
	{
		std::sort(ids_on_edge.begin(),
				  ids_on_edge.end(),
				  Less_along_a_sm_halfedge<SurfaceMesh, VPMSM, INodes>
					(h_split, sm, vpm_sm, nodes));
	}

	// split every intersected edge in surface mesh.
	void split_sm_halfedges(OnSMEdge& on_sm_edge,
							SurfaceMesh& sm,
							VPMSM& vpm_sm,
							INodes& nodes,
							SMFaceBoundaries& sm_face_boundaries)
	{
		// ���������ཻ��sm edge, ������, �����µ�, ά��������.
		// Note: split����Щ�������Ƕ������Χ�ɵ�, ������������.
		for (typename OnSMEdge::iterator it = on_sm_edge.begin();
										 it != on_sm_edge.end();
										 ++it)
		{
			SMHalfedge h_split = halfedge(it->first, sm); // the edge to be splited.
			NodeIdVec& ids_on_edge = it->second; // node ids to be inserted.

			// ȷ��ids_on_edge��û���ظ���node_id
			assert(std::set<NodeId>(ids_on_edge.begin(), ids_on_edge.end()).size() == ids_on_edge.size());

			sort_vertices_along_sm_halfedge(ids_on_edge, h_split, sm, vpm_sm, nodes);

			// ����h_split����faceԭʼ����Ϣ, ��ids_on_edge��������Ӧ��face boundary��.
			if (!is_border(h_split, sm))
			{
				SMFace f = face(h_split, sm);
				typename SMFaceBoundaries::iterator it_face = sm_face_boundaries.find(f);

				if (it_face == sm_face_boundaries.end())
				{
					SMFaceBoundary face_bound(h_split, sm);
					it_face = sm_face_boundaries.insert(std::make_pair(f, face_bound)).first;
				}
				it_face->second.copy_ids_on_edge(h_split, ids_on_edge.begin(), ids_on_edge.end());
			}

			// ����h_split�ֵܰ������faceԭʼ����Ϣ, ����.
			SMHalfedge opp_h_split = opposite(h_split, sm);
			typename SMFaceBoundaries::iterator it_opp_face = sm_face_boundaries.end();
			if (!is_border(opp_h_split, sm))
			{
				SMFace opp_face = face(opp_h_split, sm);
				it_opp_face = sm_face_boundaries.find(opp_face);
				if (it_opp_face == sm_face_boundaries.end())
				{
					SMFaceBoundary opp_face_bound(opp_h_split, sm);
					it_opp_face = sm_face_boundaries.insert(std::make_pair(opp_face, opp_face_bound)).first;
				}
				// ע��, �ֵܰ����Ҫ����.
				it_opp_face->second.copy_ids_on_edge(opp_h_split, ids_on_edge.rbegin(), ids_on_edge.rend());
			}

			SMVertex original_srcv = source(h_split, sm);

			bool first = true;
			SMHalfedge h_incident_to_src = GT::null_halfedge();
			SMVertex expected_src = source(h_split, sm);

			// insert intersection nodes consecutively.
			for (NodeId node_id : ids_on_edge)
			{
				// �����±ߺ��µ�, ������������, ��ʱ�Ѿ���������surface mesh��.
				SMHalfedge hnew = CGAL::Euler::split_edge(h_split, sm);
				assert(expected_src == source(hnew, sm));
				SMVertex vnew = target(hnew, sm);
				nodes.call_put(vpm_sm, vnew, node_id, sm);
				output_builder.set_sm_vertex_id(vnew, node_id);
				node_id_to_sm_vertex[node_id] = vnew;

				if (first)
				{
					first = false;
					h_incident_to_src = next(opposite(h_split, sm), sm);
				}
				expected_src = vnew;
			}

			assert(target(h_incident_to_src, sm) == original_srcv);
			assert(face(h_incident_to_src, sm) == face(opp_h_split, sm));

			// ����opp_h_split��Ӧface boundary�е�vertices[0,1]
			// ��h_split�ϲ������ɽ����, opp_h_split��source��Ȼ��vertices[0], ��target����vertices[1],
			// ���и���Ӧ������: opp_h_split��source����vertices[0], ��target��Ȼ��vertices[1], 
			// ������split��ʱ����h_incident_to_src��������Ҫ���������, 
			// ����opp_h_split��Ӧ��face boundary�Ѿ�����ʱ, Ӧ������Ϊh_incident_to_src, �����Ӱ����һ�������ǻ�.
			//
			if (!is_border(opp_h_split, sm))
			{
				assert(it_opp_face != sm_face_boundaries.end());
				it_opp_face->second.update_original_halfedge(
					opp_h_split, h_incident_to_src, sm);
			}

			// ֮ǰ�Ĳ�����, on_sm_face��¼����ֻ�����ڲ��Ľ���, ��Ҳ��Щ��ֻ�б����н���,
			// ��û���ڲ�����������on_sm_face, ��֤�����ཻ���涼�����ǻ�.
			if (!is_border(h_split, sm))
				on_sm_face[face(h_split, sm)];
			if (!is_border(opp_h_split, sm))
				on_sm_face[face(opp_h_split, sm)];
		}
	}

	struct PMFaceBoundary
	{
		PMVertex                     vertices[3];
		PMHalfedge                   halfedges[3];
		std::map<PMHalfedge, int>    hedge_id_map;
		std::vector<NodeId>          ids_on_edge[3];
		std::map<SortedPair<PMVertex>, PMEdge>      endpoint_to_edge;

		PMFaceBoundary(PMHalfedge first, PolyhedralMesh& pm)
		{
			halfedges[0] = first;
			halfedges[1] = pm.next(first);
			halfedges[2] = pm.next(halfedges[1]);

			vertices[0] = pm.source(halfedges[0]);
			vertices[1] = pm.source(halfedges[1]);
			vertices[2] = pm.source(halfedges[2]);

			hedge_id_map.insert(std::make_pair(halfedges[0], 0));
			hedge_id_map.insert(std::make_pair(halfedges[1], 1));
			hedge_id_map.insert(std::make_pair(halfedges[2], 2));

			for (int i = 0; i < 3; ++i)
			{
				SortedPair<PMVertex> vpair(vertices[i], vertices[(i + 1) % 3]);
				endpoint_to_edge.insert(std::make_pair(vpair, pm.edge(halfedges[i])));
			}
		}

		template <typename Iterator>
		void copy_ids_on_edge(PMHalfedge h, Iterator begin, Iterator end)
		{
			typename std::map<PMHalfedge, int>::iterator it
				= hedge_id_map.find(h);
			assert(it != hedge_id_map.end());
			std::copy(begin, end, std::back_inserter(ids_on_edge[it->second]));
		}

		void record_endpoint_to_edge(PMVertex v1, PMVertex v2, PMEdge e)
		{
			SortedPair<PMVertex> vpair(v1, v2);
			endpoint_to_edge.insert(std::make_pair(vpair, e));
		}

		void remove_endpoint_to_edge(PMVertex v1, PMVertex v2)
		{
			SortedPair<PMVertex> vpair(v1, v2);
			assert(endpoint_to_edge.count(vpair));
			endpoint_to_edge.erase(vpair);
		}
	};
	typedef std::unordered_map<PMFace, PMFaceBoundary> PMFaceBoundaries;

	// ��h_split�ϵĵ��������, ����Զ��������.
	// �������ids_on_edge������ԽС�ĵ���pm.target(h_split)ԽԶ.
	void sort_vertices_along_pm_halfedge(NodeIdVec& ids_on_edge,
										 PMHalfedge h_split,
										 PolyhedralMesh& pm,
										 VPMPM& vpm_pm,
										 INodes& nodes)
	{
		std::sort(ids_on_edge.begin(),
				  ids_on_edge.end(),
				  Less_along_a_pm_halfedge<PolyhedralMesh, VPMPM, INodes>
					(h_split, pm, vpm_pm, nodes));
	}

	// split every intersected edge in polyhedral mesh.
	// 
	// Note: ��split_sm_halfedges()��ͬ, split_pm_halfedges()������ϱ�, Ҳ����������.
	// ��������������������������û���κθı�, ����ֻ���ռ��±�����˵����Ϣ.
	//
	void split_pm_halfedges(OnPMEdge& on_pm_edge,
							PolyhedralMesh& pm,
							VPMPM& vpm_pm,
							INodes& nodes,
							PMFaceBoundaries& pm_face_boundaries)
	{
		for (typename OnPMEdge::iterator it = on_pm_edge.begin();
										 it != on_pm_edge.end();
										 ++it)
		{
			PMEdge e_split = it->first;
			NodeIdVec& ids_on_edge = it->second;

			// ȷ��ids_on_edge��û���ظ���node_id.
			assert(std::set<NodeId>(ids_on_edge.begin(), ids_on_edge.end()).size() == ids_on_edge.size());

			PMHalfedge h_split = pm.halfedge(e_split);
			sort_vertices_along_pm_halfedge(ids_on_edge, h_split, pm, vpm_pm, nodes);

			std::vector<std::pair<PMFace, bool>> split_faces;
			pm.collect_split_faces(h_split, split_faces);

			for (std::size_t i = 0; i < split_faces.size(); ++i)
			{
				PMFace f = split_faces[i].first;
				assert(pm.is_triangle(f));
				bool is_same_direction = split_faces[i].second;

				PMHalfedge h = pm.halfedge(f);
				while (pm.edge(h) != e_split)
					h = pm.next(h);

				if (is_same_direction)
					assert(pm.target(h) == pm.target(h_split) && pm.source(h) == pm.source(h_split));
				else
					assert(pm.target(h) == pm.source(h_split) && pm.source(h) == pm.target(h_split));

				typename PMFaceBoundaries::iterator it_fb = pm_face_boundaries.find(f);
				if (it_fb == pm_face_boundaries.end())
				{
					PMFaceBoundary face_bound(h, pm);
					it_fb = pm_face_boundaries.insert(std::make_pair(f, face_bound)).first;
				}

				if (is_same_direction)
					it_fb->second.copy_ids_on_edge(h, ids_on_edge.begin(), ids_on_edge.end());
				else
					it_fb->second.copy_ids_on_edge(h, ids_on_edge.rbegin(), ids_on_edge.rend());
			}

			PMVertex original_srcv = pm.source(h_split);
			PMVertex original_tgtv = pm.target(h_split);

			PMVertex pre_vertex = pm.source(h_split);
			for (NodeId node_id : ids_on_edge)
			{
				PMEdge enew = pm.add_edge();
				PMVertex vnew = pm.add_vertex();
				nodes.call_put(vpm_pm, vnew, node_id, pm);
				output_builder.set_pm_vertex_id(vnew, node_id);
				node_id_to_pm_vertex[node_id] = vnew;

				for (std::size_t i = 0; i < split_faces.size(); ++i)
				{
					typename PMFaceBoundaries::iterator it_fb =
						pm_face_boundaries.find(split_faces[i].first);
					assert(it_fb != pm_face_boundaries.end());
					it_fb->second.record_endpoint_to_edge(pre_vertex, vnew, enew);
				}
				pre_vertex = vnew;
			}

			// ȷ��û���޸�ԭ��������.
			assert(pm.source(h_split) == original_srcv);
			assert(pm.target(h_split) == original_tgtv);

			for (std::size_t i = 0; i < split_faces.size(); ++i)
			{
				PMFace f = split_faces[i].first;
				typename PMFaceBoundaries::iterator it_fb = pm_face_boundaries.find(f);
				assert(it_fb != pm_face_boundaries.end());
				it_fb->second.record_endpoint_to_edge(pre_vertex, original_tgtv, e_split);
				if (!ids_on_edge.empty())
					it_fb->second.remove_endpoint_to_edge(original_srcv, original_tgtv);

				on_pm_face[f];
			}
		}
	}

	typename CDT::Vertex_handle
	insert_point_on_edge(CDT& cdt,
						 typename CDT::Face_handle& fh,
						 const typename CDT::Point& p)
	{
		assert(cdt.is_infinite(fh));
		int fi = fh->index(cdt.infinite_vertex());
		typename CDT::Vertex_handle vh = cdt.insert(p, CDT::EDGE, fh, fi);
		// �ҵ�infinite_edge, �������˵�ֱ�Ϊ������infinite_vertex�ı�
		typename CDT::Edge_circulator ec = cdt.incident_edges(vh);
		while (ec->first->vertex(CDT::ccw(ec->second)) != cdt.infinite_vertex())
			++ec;

		// ����fh
		fh = ec->first->neighbor(ec->second);
		assert(cdt.is_valid());
		return vh;
	}

	void insert_constrained_edges(NodeIdVec& node_ids,
								  CDT& cdt,
								  std::map<NodeId, CDT_Vertex_handle>& id_to_CDT_vh,
								  std::vector<std::pair<NodeId, NodeId>>& constrained_edges)
	{
		// ����node_ids�е�id, �ҵ����������ھ�, ������������һ������(Լ����), ������ֻҪ��Щ��ϸ�����ϵ��ھ�,
		// ��ô�ж���? 
		// ǰ��Ĳ�����, id_to_CDT_vh�Ѿ���ϸ�����ϵ����е㶼��¼�ڰ���, ֻҪ�ܲ鵽, ��������.
		// �˴�Ҳ������Ϊʲô���ǲ�����id���ھӴ�������, ��Ϊ��Ӱ�������߼�.
		for (NodeId id : node_ids)
		{
			assert(id < graph_of_constraints.size());
			NodeIdVec& neighbors = graph_of_constraints[id];
			if (!neighbors.empty())
			{
				CDT_Vertex_handle vh = id_to_CDT_vh.find(id)->second;
				for (NodeId id_neighbor : neighbors)
				{
					typename std::map<NodeId, CDT_Vertex_handle>::iterator it_vh
						= id_to_CDT_vh.find(id_neighbor);

					// ��֤��ͬһ����������, �������ڵ�CDT_Vertex_handle����id_to_CDT_vh��, �������Ҳ���.
					if (it_vh != id_to_CDT_vh.end())
					{
						cdt.insert_constraint(vh, it_vh->second); // ����Լ����.
						constrained_edges.push_back(std::make_pair(id, id_neighbor));
						constrained_edges.push_back(std::make_pair(id_neighbor, id));
					}
				}
			}
		}
	}

	void update_sm_face_indices(std::array<SMVertex, 3>& f_vertices,
								std::array<NodeId, 3>& f_indices)
	{
		for (int i = 0; i < 3; ++i)
		{
			typename SMVertexToNode::iterator it =
				sm_vertex_to_node_id.find(f_vertices[i]);
			if (it != sm_vertex_to_node_id.end())
				f_indices[i] = it->second;
		}
	}

	// ��CDTϸ���Ľ�����뵽�������������.
	void import_cdt_to_sm_face(SMFace current_face,
							   CDT& cdt,
							   std::vector<NodeId>& ids_in_face,
							   std::map<std::pair<NodeId, NodeId>, SMHalfedge>& edge_to_halfedge,
							   SurfaceMesh& sm,
							   VPMSM& vpm_sm,
							   INodes& nodes)
	{
		assert(ids_in_face.size() == std::set<NodeId>(ids_in_face.begin(), ids_in_face.end()).size());

		// 1) �����ڲ���.
		// 
		// �����α��ϵĵ��Ѿ���ӵ�������������, �����ڲ��ĵ㻹û�в���, ������ӵ�����������.
		// insert the intersection point interior to the face inside the mesh
		// and save their vertex_descriptor.
		for (NodeId node_id : ids_in_face)
		{
			SMVertex v = add_vertex(sm);
			nodes.call_put(vpm_sm, v, node_id, sm);

			output_builder.set_sm_vertex_id(v, node_id);
			assert(node_id < node_id_to_sm_vertex.size());
			node_id_to_sm_vertex[node_id] = v;
		}

		// 2) ����߲���������(��û��ȫ��������, ��Ϊ��ûadd_face, ֻ������ָ���vertex).
		// 
		// insert the new halfedge and set their incident vertex.
		for (typename CDT::Finite_edges_iterator it = cdt.finite_edges_begin();
												 it != cdt.finite_edges_end();
												 ++it)
		{
			typename CDT::Vertex_handle cdt_v0 = it->first->vertex(cdt.ccw(it->second));
			typename CDT::Vertex_handle cdt_v1 = it->first->vertex(cdt.cw(it->second));

			// �����εı��Ѿ���split_halfedges()��ʱ�򶼵�������, ����.
			// grab edges that are not on the convex hull (these have already been created).
			if (!cdt.is_infinite(it->first->vertex(it->second)) &&
				!cdt.is_infinite(cdt.mirror_vertex(it->first, it->second)))
			{
				SMEdge e = add_edge(sm);
				SMHalfedge h = halfedge(e, sm), h_opp = opposite(h, sm);

				NodeId i0 = cdt_v0->info(), i1 = cdt_v1->info();
				assert(node_id_to_sm_vertex[i0] != GT::null_vertex());
				assert(node_id_to_sm_vertex[i1] != GT::null_vertex());

				SMVertex v0 = node_id_to_sm_vertex[i0];
				SMVertex v1 = node_id_to_sm_vertex[i1];

				set_target(h, v0, sm);
				set_target(h_opp, v1, sm);
				set_halfedge(v0, h, sm);
				set_halfedge(v1, h_opp, sm);

				edge_to_halfedge[std::make_pair(i0, i1)] = h_opp;
				edge_to_halfedge[std::make_pair(i1, i0)] = h;
			}
		}

		// 3) ����沢����������. ע��ԭ������Ҫ����, ����ɾ��.
		for (typename CDT::Finite_faces_iterator it = cdt.finite_faces_begin(),
				it_end = cdt.finite_faces_end();;)
		{
			// CDT����������������Ϊ012, ����ʱ��˳��.
			typename CDT::Vertex_handle cdt_v0 = it->vertex(0);
			typename CDT::Vertex_handle cdt_v1 = it->vertex(1);
			typename CDT::Vertex_handle cdt_v2 = it->vertex(2);

			NodeId i0 = cdt_v0->info(), i1 = cdt_v1->info(), i2 = cdt_v2->info();

			assert(edge_to_halfedge.count(std::make_pair(i0, i1)) != 0);
			assert(edge_to_halfedge.count(std::make_pair(i1, i2)) != 0);
			assert(edge_to_halfedge.count(std::make_pair(i2, i0)) != 0);

			SMHalfedge h01 = edge_to_halfedge[std::make_pair(i0, i1)];
			SMHalfedge h12 = edge_to_halfedge[std::make_pair(i1, i2)];
			SMHalfedge h20 = edge_to_halfedge[std::make_pair(i2, i0)];

			assert(target(h01, sm) == node_id_to_sm_vertex[i1]);
			assert(target(h12, sm) == node_id_to_sm_vertex[i2]);
			assert(target(h20, sm) == node_id_to_sm_vertex[i0]);

			set_next(h01, h12, sm);
			set_next(h12, h20, sm);
			set_next(h20, h01, sm);

			//update face halfedge
			set_halfedge(current_face, h01, sm);

			//update face of halfedges
			set_face(h01, current_face, sm);
			set_face(h12, current_face, sm);
			set_face(h20, current_face, sm);

			if (++it != it_end)
			{
				current_face = add_face(sm);
			}
			else
				break;
		}
	}

	// �������ཻ�������ǻ�.
	void triangulate_sm_intersected_faces(OnSMFace& on_sm_face,
										  SurfaceMesh& sm,
										  VPMSM& vpm_sm,
										  INodes& nodes,
										  SMFaceBoundaries& sm_face_boundaries)
	{
		const std::size_t nb_nodes = nodes.size();

		// �������е�sm face, �������ǻ�.
		for (typename OnSMFace::iterator it = on_sm_face.begin();
										 it != on_sm_face.end();
										 ++it)
		{
			SMFace f = it->first; // the face to be triangulated.
			NodeIdVec& ids_in_face = it->second;  // node ids in the interior of f.

			std::map<NodeId, CDT_Vertex_handle> id_to_CDT_vh;
			std::map<std::pair<NodeId, NodeId>, SMHalfedge> edge_to_halfedge;
			std::array<SMVertex, 3> f_vertices; // f����������.
			std::array<NodeId, 3> f_indices = { nb_nodes, nb_nodes + 1, nb_nodes + 2 };

			// Step1 CDT���ǻ�.
			// Note: 
			// CDT������/��/�������Լ���һ�ױ�ʾ, ��Face_handle/Vertex_handle.
			// Vertex_handle��ģ�������һ��Info, �˴�����ʹ��NodeId.
			// 
			// CDT�Ĵ�������Ϊ: �Ȳ�����������������, �ٲ�����Ϻ����ڲ��ĵ�, ���������Delaunay���ǻ�, 
			//                 ������Լ����, �㷨������Լ�����ཻ�ı�, �γ����������, �ٸ������ǻ�.
			// 
			// ���ǻ����������CDT��iteratorȡ�����ǻ�������бߺ���.
			// 
			// ���Խ�CDT��Ϊһ���ں�, ������Ҫ�����ǽ�������node id��CDT_Vertex_handle�Ķ�Ӧ��ϵ, 

			// 1.1) ����ԭʼ����������.
			typename SMFaceBoundaries::iterator it_fb = sm_face_boundaries.find(f);
			if (it_fb != sm_face_boundaries.end())
			{
				f_vertices[0] = it_fb->second.vertices[0];
				f_vertices[1] = it_fb->second.vertices[1];
				f_vertices[2] = it_fb->second.vertices[2];

				update_sm_face_indices(f_vertices, f_indices);
			}
			else
			{
				assert(is_triangle(halfedge(f, sm), sm));
				SMHalfedge h0 = halfedge(f, sm), h1 = next(h0, sm), h2 = next(h1, sm);
				f_vertices[0] = target(h0, sm); //nb_nodes
				f_vertices[1] = target(h1, sm); //nb_nodes+1
				f_vertices[2] = target(h2, sm); //nb_nodes+2

				update_sm_face_indices(f_vertices, f_indices);
				edge_to_halfedge[std::make_pair(f_indices[2], f_indices[0])] = h0;
				edge_to_halfedge[std::make_pair(f_indices[0], f_indices[1])] = h1;
				edge_to_halfedge[std::make_pair(f_indices[1], f_indices[2])] = h2;
			}

			if (CGAL::collinear(get(vpm_sm, f_vertices[0]), get(vpm_sm, f_vertices[1]), get(vpm_sm, f_vertices[2])))
			{
				std::cout << "�˻����";
			}

			typename EK::Point_3 p = nodes.to_exact(get(vpm_sm, f_vertices[0]));
			typename EK::Point_3 q = nodes.to_exact(get(vpm_sm, f_vertices[1]));
			typename EK::Point_3 r = nodes.to_exact(get(vpm_sm, f_vertices[2]));

			CDTtraits traits(typename EK::Construct_normal_3()(p, q, r));
			CDT cdt(traits);

			std::array<CDT_Vertex_handle, 3> triangle_vertices;
			triangle_vertices[0] = cdt.insert_outside_affine_hull(p);
			triangle_vertices[1] = cdt.insert_outside_affine_hull(q);
			triangle_vertices[2] = cdt.tds().insert_dim_up(cdt.infinite_vertex(), false);
			triangle_vertices[2]->set_point(r);

			triangle_vertices[0]->info() = f_indices[0];
			triangle_vertices[1]->info() = f_indices[1];
			triangle_vertices[2]->info() = f_indices[2];

			// Note: node_id_to_sm_vertex���������Ԫ�ض�Ӧ��vertex��仯, ��Ϊnb_nodes~nb_nodes+2��fake node_id.
			node_id_to_sm_vertex[nb_nodes] = f_vertices[0];
			node_id_to_sm_vertex[nb_nodes + 1] = f_vertices[1];
			node_id_to_sm_vertex[nb_nodes + 2] = f_vertices[2];

			for (int i = 0; i < 3; ++i)
			{
				if (f_indices[i] < nb_nodes)
				{
					id_to_CDT_vh.insert(std::make_pair(f_indices[i], triangle_vertices[i]));
				}
			}

			// 1.2) ���������θ����ϵĽ���.
			if (it_fb != sm_face_boundaries.end())
			{
				// collect infinite faces incident to the initial triangle
				typename CDT::Face_handle infinite_faces[3];
				for (int i = 0; i < 3; ++i)
				{
					int oi = -1;
					bool is_edge =
						cdt.is_edge(triangle_vertices[i], triangle_vertices[(i + 1) % 3], infinite_faces[i], oi);
					assert(is_edge);
					assert(cdt.is_infinite(infinite_faces[i]->vertex(oi)));
				}

				SMFaceBoundary& f_bound = it_fb->second;
				for (int i = 0; i < 3; ++i)
				{
					NodeIdVec& ids_on_edge = f_bound.ids_on_edge[i];
					CDT_Vertex_handle prev_CDT_vh = triangle_vertices[i];
					NodeId prev_node_id = f_indices[i];
					// �˴���h_sm�ļ��㷽ʽҪ�ο�SMFaceBoundary�ж�Ӧ��ϵ.
					SMHalfedge h_sm = next(f_bound.halfedges[(i + 2) % 3], sm);
					assert(source(h_sm, sm) == f_bound.vertices[i]);
					if (!ids_on_edge.empty()) // ������������һ�������н���.
					{
						for (NodeId id : ids_on_edge)
						{
							CDT_Vertex_handle vh =
								insert_point_on_edge(cdt, infinite_faces[i], nodes.exact_node(id));
							vh->info() = id;
							id_to_CDT_vh.insert(std::make_pair(id, vh));
							edge_to_halfedge[std::make_pair(prev_node_id, id)] = h_sm;
							prev_CDT_vh = vh;
							prev_node_id = id;
							h_sm = next(h_sm, sm);
						}
					}
					else
					{
						SMHalfedge h = f_bound.halfedges[i];
						assert(target(h, sm) == f_bound.vertices[(i + 1) % 3]);
						assert(source(h, sm) == f_bound.vertices[i]);
					}
					assert(h_sm == f_bound.halfedges[i]);
					edge_to_halfedge[std::make_pair(prev_node_id, f_indices[(i + 1) % 3])] =
						it_fb->second.halfedges[i];
				}
			}

			// 1.3) �����������ڲ��Ľ���.
			for (NodeId node_id : ids_in_face)
			{
				CDT_Vertex_handle vh = cdt.insert(nodes.exact_node(node_id));
				vh->info() = node_id;
				id_to_CDT_vh.insert(std::make_pair(node_id, vh));
			}

			// 1.4) ����Լ����(����), �����Ǳ��뱣����.
			std::vector<std::pair<NodeId, NodeId>> constrained_edges;

			// insert constraints that are interior to the triangle.
			insert_constrained_edges(ids_in_face, cdt, id_to_CDT_vh, constrained_edges);

			// insert constraints between points that are on the boundary.
			if (it_fb != sm_face_boundaries.end())
			{
				for (int i = 0; i < 3; ++i)
				{
					NodeIdVec& ids_on_edge = it_fb->second.ids_on_edge[i];
					insert_constrained_edges(ids_on_edge, cdt, id_to_CDT_vh, constrained_edges);
				}
			}

			// Step2 ��CDTϸ���Ľ�����뵽��������, �޸�����.
			// 
			// import the triangle in `cdt` to the face `f` of `sm`.
			// 
			// Note: ��Step1��, �Ѿ���������ǻ�, �����ǻ��Ľ��ȫ����CDT������,
			//       ��Ҫ�����ǻ��Ľ��ȡ��, ���õ���f��, ����޸�sm.
			//
			import_cdt_to_sm_face(f, cdt, ids_in_face, edge_to_halfedge, sm, vpm_sm, nodes);

			// Step3 �ռ�������Ϣ����¼.
			// ��������һ����ȡ����Ҫ��, ֮ǰ�Ĳ�������δ����, ��ʱ�Ѿ�����, ��Ҫ��¼��OutputBuilder�����ݳ�Ա��.
			for (const NodeIdPair& id_pair : constrained_edges)
			{
				typename std::map<NodeIdPair, SMHalfedge>::iterator
					it_poly_hedge = edge_to_halfedge.find(id_pair);

				assert(it_poly_hedge != edge_to_halfedge.end());
				if (it_poly_hedge != edge_to_halfedge.end())
				{
					output_builder.set_sm_edge_per_polyline(sm, id_pair, it_poly_hedge->second);
				}
				else
				{
					NodeIdPair opp_id_pair(id_pair.second, id_pair.first);
					it_poly_hedge = edge_to_halfedge.find(opp_id_pair);
					assert(it_poly_hedge != edge_to_halfedge.end());

					output_builder.set_sm_edge_per_polyline(sm, opp_id_pair, it_poly_hedge->second);
				}
			}
		}
	}

	void update_pm_face_indices(std::array<PMVertex, 3>& f_vertices,
								std::array<NodeId, 3>& f_indices)
	{
		for (int i = 0; i < 3; ++i)
		{
			typename PMVertexToNode::iterator it =
				pm_vertex_to_node_id.find(f_vertices[i]);
			if (it != pm_vertex_to_node_id.end())
				f_indices[i] = it->second;
		}
	}

	void collect_triangulation_info(PMFace current_face,
									std::vector<NodeId>& ids_in_face,
									CDT& cdt,
									PolyhedralMesh& pm,
									VPMPM& vpm_pm,
									INodes& nodes,
									std::unordered_map<PMFace, std::vector<PMVertex>>& f_to_vseq,
									std::map<SortedPair<PMVertex>, PMEdge>& endpoint_to_edge)
	{
		assert(std::set<NodeId>(ids_in_face.begin(), ids_in_face.end()).size() == ids_in_face.size());

		for (NodeId node_id : ids_in_face)
		{
			PMVertex v = pm.add_vertex();
			nodes.call_put(vpm_pm, v, node_id, pm);
			output_builder.set_pm_vertex_id(v, node_id);

			assert(node_id < node_id_to_pm_vertex.size());
			node_id_to_pm_vertex[node_id] = v;
		}

		for (typename CDT::Finite_edges_iterator it = cdt.finite_edges_begin();
			it != cdt.finite_edges_end();
			++it)
		{
			typename CDT::Vertex_handle cdt_v0 = it->first->vertex(cdt.ccw(it->second));
			typename CDT::Vertex_handle cdt_v1 = it->first->vertex(cdt.cw(it->second));

			if (!cdt.is_infinite(it->first->vertex(it->second)) &&
				!cdt.is_infinite(cdt.mirror_vertex(it->first, it->second)))
			{
				PMEdge e = pm.add_edge();

				NodeId i0 = cdt_v0->info(), i1 = cdt_v1->info();
				assert(node_id_to_pm_vertex[i0] != pm.null_vertex());
				assert(node_id_to_pm_vertex[i1] != pm.null_vertex());

				PMVertex v0 = node_id_to_pm_vertex[i0];
				PMVertex v1 = node_id_to_pm_vertex[i1];

				endpoint_to_edge.insert(std::make_pair(SortedPair<PMVertex>(v0, v1), e));
			}
		}

		for (typename CDT::Finite_faces_iterator it = cdt.finite_faces_begin(),
			it_end = cdt.finite_faces_end();;)
		{
			typename CDT::Vertex_handle cdt_v0 = it->vertex(0);
			typename CDT::Vertex_handle cdt_v1 = it->vertex(1);
			typename CDT::Vertex_handle cdt_v2 = it->vertex(2);

			NodeId i0 = cdt_v0->info(), i1 = cdt_v1->info(), i2 = cdt_v2->info();

			assert(node_id_to_pm_vertex[i0] != pm.null_vertex());
			assert(node_id_to_pm_vertex[i1] != pm.null_vertex());
			assert(node_id_to_pm_vertex[i2] != pm.null_vertex());

			PMVertex v0 = node_id_to_pm_vertex[i0];
			PMVertex v1 = node_id_to_pm_vertex[i1];
			PMVertex v2 = node_id_to_pm_vertex[i2];

			std::vector<PMVertex> vseq{ v0,v1,v2 };
			f_to_vseq.insert(std::make_pair(current_face, vseq));

			if (++it != it_end)
			{
				current_face = pm.add_face();
			}
			else
				break;
		}
	}

	void triangulate_pm_intersected_faces(OnPMFace& on_pm_face,
										  PolyhedralMesh& pm,
										  VPMPM& vpm_pm,
										  INodes& nodes,
										  PMFaceBoundaries& pm_face_boundaries)
	{
		const std::size_t nb_nodes = nodes.size();

		std::map<NodeIdPair, std::pair<PMVertex, PMEdge>> build_info;

		for (typename OnPMFace::iterator it = on_pm_face.begin();
										 it != on_pm_face.end();
										 ++it)
		{
			PMFace f = it->first;
			assert(pm.is_triangle(f));
			NodeIdVec& ids_in_face = it->second;

			std::map<NodeId, CDT_Vertex_handle> id_to_CDT_vh;
			std::map<SortedPair<PMVertex>, PMEdge> endpoint_to_edge;
			std::array<PMVertex, 3> f_vertices;
			std::array<NodeId, 3> f_indices = { nb_nodes, nb_nodes + 1, nb_nodes + 2 };

			typename PMFaceBoundaries::iterator it_fb = pm_face_boundaries.find(f);
			if (it_fb != pm_face_boundaries.end())
			{
				f_vertices[0] = it_fb->second.vertices[0];
				f_vertices[1] = it_fb->second.vertices[1];
				f_vertices[2] = it_fb->second.vertices[2];

				update_pm_face_indices(f_vertices, f_indices);
				endpoint_to_edge = it_fb->second.endpoint_to_edge;
			}
			else
			{
				PMHalfedge h0 = pm.halfedge(f), h1 = pm.next(h0), h2 = pm.next(h1);
				f_vertices[0] = pm.target(h0);
				f_vertices[1] = pm.target(h1);
				f_vertices[2] = pm.target(h2);

				update_pm_face_indices(f_vertices, f_indices);

				endpoint_to_edge.insert(std::make_pair(SortedPair<PMVertex>(f_vertices[2], f_vertices[0]), pm.edge(h0)));
				endpoint_to_edge.insert(std::make_pair(SortedPair<PMVertex>(f_vertices[0], f_vertices[1]), pm.edge(h1)));
				endpoint_to_edge.insert(std::make_pair(SortedPair<PMVertex>(f_vertices[1], f_vertices[2]), pm.edge(h2)));
			}

			typename EK::Point_3 p = nodes.to_exact(vpm_pm[f_vertices[0]]);
			typename EK::Point_3 q = nodes.to_exact(vpm_pm[f_vertices[1]]);
			typename EK::Point_3 r = nodes.to_exact(vpm_pm[f_vertices[2]]);

			CDTtraits traits(typename EK::Construct_normal_3()(p, q, r));
			CDT cdt(traits);

			std::array<CDT_Vertex_handle, 3> triangle_vertices;
			triangle_vertices[0] = cdt.insert_outside_affine_hull(p);
			triangle_vertices[1] = cdt.insert_outside_affine_hull(q);
			triangle_vertices[2] = cdt.tds().insert_dim_up(cdt.infinite_vertex(), false);
			triangle_vertices[2]->set_point(r);

			triangle_vertices[0]->info() = f_indices[0];
			triangle_vertices[1]->info() = f_indices[1];
			triangle_vertices[2]->info() = f_indices[2];

			node_id_to_pm_vertex[nb_nodes] = f_vertices[0];
			node_id_to_pm_vertex[nb_nodes + 1] = f_vertices[1];
			node_id_to_pm_vertex[nb_nodes + 2] = f_vertices[2];

			for (int i = 0; i < 3; ++i)
			{
				if (f_indices[i] < nb_nodes)
				{
					id_to_CDT_vh.insert(std::make_pair(f_indices[i], triangle_vertices[i]));
				}
			}

			if (it_fb != pm_face_boundaries.end())
			{
				typename CDT::Face_handle infinite_faces[3];
				for (int i = 0; i < 3; ++i)
				{
					int oi = -1;
					bool is_edge =
						cdt.is_edge(triangle_vertices[i], triangle_vertices[(i + 1) % 3], infinite_faces[i], oi);

					assert(is_edge);
					assert(cdt.is_infinite(infinite_faces[i]->vertex(oi)));
				}

				PMFaceBoundary& f_bound = it_fb->second;
				for (int i = 0; i < 3; ++i)
				{
					NodeIdVec& ids_on_edge = f_bound.ids_on_edge[i];
					CDT_Vertex_handle prev_CDT_vh = triangle_vertices[i];
					NodeId prev_node_id = f_indices[i];

					if (!ids_on_edge.empty())
					{
						for (NodeId id : ids_on_edge)
						{
							CDT_Vertex_handle vh =
								insert_point_on_edge(cdt, infinite_faces[i], nodes.exact_node(id));
							vh->info() = id;
							id_to_CDT_vh.insert(std::make_pair(id, vh));
							prev_CDT_vh = vh;
							prev_node_id = id;
						}
					}
					//else
					//{
					//	PMHalfedge h = f_bound.halfedges[i];
					//	assert(pm.target(h) == f_bound.vertices[(i+1)%3]);
					//	assert(pm.source(h) == f_bound.vertices[i]);
					//}
				}
			}

			for (NodeId node_id : ids_in_face)
			{
				CDT_Vertex_handle vh = cdt.insert(nodes.exact_node(node_id));
				vh->info() = node_id;
				id_to_CDT_vh.insert(std::make_pair(node_id, vh));
			}

			std::vector<std::pair<NodeId, NodeId>> constrained_edges;

			// insert constraints that are interior to the triangle.
			insert_constrained_edges(ids_in_face, cdt, id_to_CDT_vh, constrained_edges);

			// insert constraints between points that are on the boundary.
			if (it_fb != pm_face_boundaries.end())
			{
				for (int i = 0; i < 3; ++i)
				{
					NodeIdVec& ids_on_edge = it_fb->second.ids_on_edge[i];
					insert_constrained_edges(ids_on_edge, cdt, id_to_CDT_vh, constrained_edges);
				}
			}

			// Step2 ��CDTϸ���Ľ�����뵽����������, �޸�����.
			std::unordered_map<PMFace, std::vector<PMVertex>> f_to_vseq;
			collect_triangulation_info(f, ids_in_face, cdt, pm, vpm_pm, nodes, f_to_vseq, endpoint_to_edge);

			pm.import_triangulation_info(f, f_to_vseq, endpoint_to_edge);

			// Step3 �ռ�������Ϣ����¼.
			// ��������һ����ȡ����Ҫ��, ֮ǰ�Ĳ�������δ����, ��ʱ�Ѿ�����, ��Ҫ��¼��OutputBuilder�����ݳ�Ա��.
			for (const NodeIdPair& id_pair : constrained_edges)
			{
				NodeId i0 = id_pair.first;
				NodeId i1 = id_pair.second;

				if (i0 > i1)
					continue;

				assert(node_id_to_pm_vertex[i0] != pm.null_vertex());
				assert(node_id_to_pm_vertex[i1] != pm.null_vertex());

				PMVertex v0 = node_id_to_pm_vertex[i0];
				PMVertex v1 = node_id_to_pm_vertex[i1];
				SortedPair<PMVertex> vpair(v0, v1);

				typename std::map<SortedPair<PMVertex>, PMEdge>::iterator
					it_poly_edge = endpoint_to_edge.find(vpair);
				assert(it_poly_edge != endpoint_to_edge.end());

				build_info.insert(std::make_pair(id_pair, std::make_pair(v0, it_poly_edge->second)));
			}
		}

		for (typename std::map<NodeIdPair, std::pair<PMVertex, PMEdge>>::iterator
			it = build_info.begin(); it != build_info.end(); ++it)
		{
			NodeIdPair id_pair = it->first;
			PMVertex srcv = it->second.first;
			PMEdge target_edge = it->second.second;
			output_builder.set_pm_edge_per_polyline(id_pair, srcv, target_edge, pm);
		}
	}

	// finalize: ����Ϊ����ȷ��, �ڴ˵���˼�����ǻ��ཻ�������õ�������.
	void finalize(INodes& nodes,
				  SurfaceMesh& sm,
				  PolyhedralMesh& pm,
				  VPMSM& vpm_sm,
				  VPMPM& vpm_pm)
	{
		const std::size_t nb_nodes = nodes.size();

		node_id_to_sm_vertex.resize(nb_nodes + 3, GT::null_vertex());
		node_id_to_pm_vertex.resize(nb_nodes + 3, pm.null_vertex());

		SMFaceBoundaries sm_face_boundaries;
		PMFaceBoundaries pm_face_boundaries;

		split_sm_halfedges(on_sm_edge, sm, vpm_sm, nodes, sm_face_boundaries);
		split_pm_halfedges(on_pm_edge, pm, vpm_pm, nodes, pm_face_boundaries);
		triangulate_sm_intersected_faces(on_sm_face, sm, vpm_sm, nodes, sm_face_boundaries);
		triangulate_pm_intersected_faces(on_pm_face, pm, vpm_pm, nodes, pm_face_boundaries);

		nodes.finalize(node_id_to_sm_vertex, node_id_to_pm_vertex);

		node_id_to_sm_vertex.resize(nb_nodes);
		node_id_to_pm_vertex.resize(nb_nodes);
		output_builder(nodes, node_id_to_sm_vertex, node_id_to_pm_vertex);
	}

};

} // namespace MCAL


#endif