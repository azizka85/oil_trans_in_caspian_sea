#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_VISCOSITY_IGENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_VISCOSITY_IGENERATOR_H

#include <vector>

#include <filesystem>

using namespace std;
using namespace std::filesystem;

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Viscosity {
	class IGenerator {
	public:
		virtual path createDirectory(path outDir) = 0;

		virtual vector<float> generateNU(
			int nx, int ny, int nz,
			vector<float>& dz, vector<float>& h,
			vector<float>& u10, vector<float>& v10,
			vector<float> &qx, vector<float>& qy,
			vector<float> &ua, vector<float> &va
		) = 0;
	};
}

#endif