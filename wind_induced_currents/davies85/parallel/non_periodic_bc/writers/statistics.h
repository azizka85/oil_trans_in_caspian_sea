#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_WRITERS_STATISTICS_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_WRITERS_STATISTICS_H

#include <tuple>
#include <vector>

#include <filesystem>

using namespace std;
using namespace std::filesystem;

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Writers::Statistics {
    void write(
        vector<tuple<int, float, long long, float, float, float>>& statistics,
        path outDir
    );
}

#endif