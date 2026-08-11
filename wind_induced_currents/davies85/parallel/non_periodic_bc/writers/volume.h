#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_WRITERS_VOLUME_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_WRITERS_VOLUME_H

#include <vector>

#include <filesystem>

using namespace std;
using namespace std::filesystem;

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Writers::Volume {
    void write(
        float t, int m,
        int nx, int ny, int nz,
        float dx, float dy,
        vector<float>& dz, vector<float>& h,
        vector<float>& ua, vector<float>& va,
        vector<float>& uf, vector<float>& vf,
        path outDir
    );

    void writeViscosity(
        float t, int m,
        int nx, int ny, int nz,
        float dx, float dy,
        vector<float>& dz, vector<float>& h,
        vector<float>& nu,
        path outDir
    );
}

#endif