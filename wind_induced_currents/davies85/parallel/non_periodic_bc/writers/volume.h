#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_WRITERS_VOLUME_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_WRITERS_VOLUME_H

#include <vector>

#include <filesystem>

using namespace std;
using namespace std::filesystem;

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Writers::Volume {
    void write(
        double t, int m,
        int nx, int ny, int nz,
        double dx, double dy,
        vector<double>& dz, vector<double>& h,
        vector<double>& ua, vector<double>& va,
        vector<double>& uf, vector<double>& vf,
        path outDir
    );

    void writeViscosity(
        double t, int m,
        int nx, int ny, int nz,
        double dx, double dy,
        vector<double>& dz, vector<double>& h,
        vector<double>& nu,
        path outDir
    );
}

#endif