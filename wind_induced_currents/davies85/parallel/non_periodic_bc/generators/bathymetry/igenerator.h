#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_BATHYMETRY_IGENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_BATHYMETRY_IGENERATOR_H

#include <vector>

#include <filesystem>

using namespace std;
using namespace std::filesystem;

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Bathymetry {
	class IGenerator {
		public:
			virtual path createDirectory(path outDir) = 0;

			virtual vector<float> generateH(int nx, int ny) = 0;
	};
}

#endif