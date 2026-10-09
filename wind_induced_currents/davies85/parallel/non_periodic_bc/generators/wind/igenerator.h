#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_WIND_IGENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_WIND_IGENERATOR_H

#include <vector>

#include <filesystem>

using namespace std;
using namespace std::filesystem;

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Wind {
	struct Data {
		int64_t time;

		vector<double> u10;
		vector<double> v10;

		vector<double> qx;
		vector<double> qy;
	};

	class IGenerator {
		public:
			virtual path createDirectory(path outDir) = 0;

			virtual vector<Data> generate(int nx, int ny) = 0;
	};
}

#endif