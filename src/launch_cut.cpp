#define CGAL_NO_STATIC_FILTERS 

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>


#include "grid/polyhedron_grid.h"
#include "io/polyhedron_grid_io.h"
#include "io/embedded_body_io.h"

#include "cut_algo/intersection_impl.h"
#include "cut_algo/intersection_visitor.h"
#include "cut_algo/output_builder.h"

#include "ctime"

typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
typedef CGAL::Surface_mesh<K::Point_3>                         SM;
// Note: �и�������Ԫ�ػᵼ�³���256, ��16λ��ʾbias_type������.
typedef MCAL::GRID::PolyhedronGrid<K::Point_3, std::uint16_t>   PG;

namespace PGIO = MCAL::IO;
namespace SMIO = CGAL::Polygon_mesh_processing::IO;

int iExactComputeCount = 0;
int iIntervalComputeCount = 0;
std::clock_t start = std::clock();

// �����̻�����exe, ����̽����, ������Ҫ��������: �ض������ļ���Ŀ¼, �и���·��, ����ĸ�Ŀ¼.
int main(int argc, char* argv[])
{
	std::string tg_root_path = (argc > 1) ? argv[1] :
		"D:\\B-���ɽ��ģ\\SimModel\\3814578003\\";
	std::string eb_file_name = (argc > 2) ? argv[2] :
		"D:\\B-���ɽ��ģ\\SimModel\\embedded\\1\\face_grid.dat";
	std::string out_path = (argc > 3) ? argv[3] : "E:/tmp";

	PG pg;
	std::vector<SM> sm_set;
	// ���ļ�, ��������������ͱ�������.
	std::string coord_txt_path = tg_root_path + "\\local_coordinate.txt";
	if (!PGIO::read_polyhedron_grid(tg_root_path, pg) ||
		!PGIO::load_embedded_body(eb_file_name, sm_set, coord_txt_path))
	{
		std::cerr << "Invalid input." << std::endl;
		return 1;
	}
	SM sm = sm_set[0];

	//CGAL::IO::write_polygon_mesh("E:\\stone.obj", sm);
	//if (CGAL::is_closed(sm) && !CGAL::Polygon_mesh_processing::does_self_intersect(sm) && CGAL::Polygon_mesh_processing::does_bound_a_volume(sm))
	//{
	//	std::cout << "sm is valid" << std::endl;
	//

	//	CGAL::IO::write_polygon_mesh("E:/stone.off", sm, CGAL::parameters::stream_precision(17));
	//}
	//else
	//	std::cerr << "sm is invalid" << std::endl;

	std::cout << "����������------------------" << std::endl;
	std::cout << "��������" << pg.num_vertices()
		<< "\n������" << pg.num_edges()
		<< "\n������" << pg.num_faces()
		<< "\n��Ԫ����" << pg.num_cells() << "\n\n";

	std::cout << "��������-----------------" << std::endl;
	std::cout << "��������" << sm.num_vertices()
		<< "\n������" << sm.num_edges()
		<< "\n������" << sm.num_faces() << "\n\n";

	// ׼����ز���, ʵ���������.
	typedef typename SM::template Property_map<SM::Vertex_index, K::Point_3> VPMSM;
	typedef typename PG::template PropertyMap<PG::VertexIndex, K::Point_3>   VPMPG;

	VPMSM vpm_sm = sm.points();    // ���������vertex point map
	VPMPG vpm_pg = pg.points();

	typedef MCAL::OutputBuilder<SM, PG, VPMSM, VPMPG> OB;
	typedef MCAL::IntersectionVisitor<SM, PG, VPMSM, VPMPG, OB> AlgoVisitor;
	typedef MCAL::IntersectionImpl<SM, PG, VPMSM, VPMPG, AlgoVisitor> InterFunctor;

	OB ob(sm, pg, vpm_sm, vpm_pg);
	InterFunctor functor(sm, pg, vpm_sm, vpm_pg, AlgoVisitor(ob));

	iExactComputeCount = 0;
	iIntervalComputeCount = 0;

	start = std::clock();

	// ��ʼ�и�
	functor();

	std::cout << "��ȷ������ܴ�����" << iExactComputeCount << std::endl;
	std::cout << "����������ܴ�����" << iIntervalComputeCount << std::endl;
	std::cout << "----------------------------------------\n\n";

	std::clock_t finish = std::clock();
	std::cout << "This program cost " << (double)(finish - start) / CLOCKS_PER_SEC << "seconds" << std::endl;

	PGIO::initialize_io(tg_root_path, out_path);
	std::string outer_path = out_path + "\\OUT";
	std::string inner_path = out_path + "\\IN";
	PGIO::write_polyhedron_grid(outer_path, pg, true);
	PGIO::write_polyhedron_grid(inner_path, pg, false);


	CGAL::IO::write_polygon_mesh(out_path + "\\out.off", sm, CGAL::parameters::stream_precision(17));

	return 0;
}