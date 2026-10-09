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

		virtual vector<double> generateNU(
			int nx, int ny, int nz,
			vector<double>& dz, vector<double>& h,
			vector<double>& u10, vector<double>& v10,
			vector<double> &qx, vector<double>& qy,
			vector<double> &ua, vector<double> &va
		) = 0;
	};
}

#endif