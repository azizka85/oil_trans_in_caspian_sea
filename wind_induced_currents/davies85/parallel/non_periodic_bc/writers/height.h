#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_WRITERS_HEIGHT_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_WRITERS_HEIGHT_H

#include <vector>

#include <filesystem>

using namespace std;
using namespace std::filesystem;

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Writers::Height {
	void write(vector<float> &h, int nx, int ny, float dx, float dy, path outDir);
}

#endif