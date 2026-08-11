#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_DZ_IGENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_DZ_IGENERATOR_H

#include <vector>

#include <filesystem>

using namespace std;
using namespace std::filesystem;

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::DZ {
	class IGenerator {
		public:			
			virtual path createDirectory(path outDir) = 0;

			virtual vector<float> generateDZ() = 0;
	};
}

#endif