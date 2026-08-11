#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_AREA_IGENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_AREA_IGENERATOR_H

#include <filesystem>

using namespace std::filesystem;

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Area {
	struct Geometry {
		float l;
		float w;
	};

	class IGenerator {
		public:
			virtual path createDirectory(path outDir) = 0;

			virtual Geometry generate() = 0;
	};
}

#endif