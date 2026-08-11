#include <utils/fs.h>

#include "statistics.h"

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Writers;

void Statistics::write(
    vector<tuple<int, float, long long, float, float, float>>& statistics,
    path outDir
) {
    auto file = FS::createFileByPath(
        outDir / path("convergence.csv")
    );

    file << "n,  t, calc_time,  u,  v,  z" << endl;

    for (auto t : statistics) {
        file << format(
            "{},  {:.3f}, {},  {:.5f},  {:.5f},  {:.5f}",
            get<0>(t),
            get<1>(t),
            get<2>(t),
            get<3>(t),
            get<4>(t),
            get<5>(t)
        ) << endl;
    }
}