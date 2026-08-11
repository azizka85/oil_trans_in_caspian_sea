#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_WIND_IGENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_WIND_IGENERATOR_H

#include <vector>

#include <filesystem>

using namespace std;
using namespace std::filesystem;

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Wind {
	struct Data {
		int64_t time;

		vector<float> u10;
		vector<float> v10;

		vector<float> qx;
		vector<float> qy;
	};

	class IGenerator {
		public:
			virtual path createDirectory(path outDir) = 0;

			virtual vector<Data> generate(int nx, int ny) = 0;
	};
}

#endif