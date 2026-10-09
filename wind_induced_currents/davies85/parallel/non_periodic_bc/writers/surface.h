#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_WRITERS_SURFACE_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_WRITERS_SURFACE_H

#include <vector>

#include <filesystem>

using namespace std;
using namespace std::filesystem;

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Writers::Surface {
	void write(
        double t, int m,
        int nx, int ny,
        double dx, double dy,
        vector<double>& ua, vector<double>& va,
        vector<double>& qx, vector<double>& qy,
        vector<double>& z, path outDir
    );
}

#endif