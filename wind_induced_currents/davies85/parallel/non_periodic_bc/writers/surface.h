#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_WRITERS_SURFACE_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_WRITERS_SURFACE_H

#include <vector>

#include <filesystem>

using namespace std;
using namespace std::filesystem;

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Writers::Surface {
	void write(
        float t, int m,
        int nx, int ny,
        float dx, float dy,
        vector<float>& ua, vector<float>& va,
        vector<float>& qx, vector<float>& qy,
        vector<float>& z, path outDir
    );
}

#endif