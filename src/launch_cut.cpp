#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>

#include "pm/polyhedral_mesh.h"
#include "io/polyhedral_mesh_io.h"
#include "io/embedded_body_io.h"

#include "algo/intersection.h"
#include "algo/corefine.h"
#include "algo/output_builder.h"


// 调用需要三个参数: 截断网格文件根目录, 切割体路径, 输出的根目录.
int main(int argc, char* argv[])
{
	typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
	typedef CGAL::Surface_mesh<K::Point_3>                         SM;
	typedef MCAL::MESH::PolyhedralMesh<K::Point_3, std::uint16_t>  PM;

	namespace PMIO = MCAL::IO;
	namespace SMIO = CGAL::Polygon_mesh_processing::IO;
	
	
	std::string tg_root_path = argv[1];
	std::string eb_file_name = argv[2];
	std::string out_path = argv[3];

	PM pm;
	std::vector<SM> sm_set;
	// 读文件, 建立多面体网格和表面网格.
	std::string coord_txt_path = tg_root_path + "\\local_coordinate.txt";
	if (!PMIO::read_polyhedron_grid(tg_root_path, pm) ||
		!PMIO::load_embedded_body(eb_file_name, sm_set, coord_txt_path))
	{
		std::cerr << "Invalid input." << std::endl;
		return 1;
	}
	SM sm = sm_set[0];

	// 准备相关参数, 实例化相关类.
	typedef typename SM::template Property_map<SM::Vertex_index, K::Point_3> VPMSM;
	typedef typename PM::template PropertyMap<PM::VertexIndex, K::Point_3>   VPMPM;

	VPMSM vpm_sm = sm.points();    // 表面网格的vertex point map
	VPMPM vpm_pm = pm.points();

	typedef MCAL::OutputBuilder<SM, PM, VPMSM, VPMPM> OB;
	typedef MCAL::Corefine<SM, PM, VPMSM, VPMPM, OB> AlgoVisitor;
	typedef MCAL::Intersection<SM, PM, VPMSM, VPMPM, AlgoVisitor> InterFunctor;

	OB ob(sm, pm, vpm_sm, vpm_pm);
	InterFunctor functor(sm, pm, vpm_sm, vpm_pm, AlgoVisitor(ob));

	// 开始切割
	functor();

	PMIO::initialize_io(tg_root_path, out_path);
	std::string outer_path = out_path + "\\OUT";
	std::string inner_path = out_path + "\\IN";
	PMIO::write_polyhedron_grid(outer_path, pm, true);
	PMIO::write_polyhedron_grid(inner_path, pm, false);

	CGAL::IO::write_polygon_mesh(out_path + "\\out.off", sm, CGAL::parameters::stream_precision(17));

	return 0;
}