#include <stdexcept>

#include <utils/fs.h>

#include "uniform_generator.h"

using namespace std;

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Viscosity;

UniformGenerator::UniformGenerator(double num) {
	setNUM(num);
}

double UniformGenerator::getNUM() {
    return num;
}

void UniformGenerator::setNUM(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("NUM should be > 0, but it is {}", val)
        );
    }

    num = val;
}

path UniformGenerator::createDirectory(path outDir) {
    return FS::createDirByPath(
        outDir / path(
            format("UNU, nu={}", num)
        )
    );
}

vector<double> UniformGenerator::generateNU(
    int nx, int ny, int nz,
    vector<double>& dz, vector<double>& h,
    vector<double>& u10, vector<double>& v10,
    vector<double>& qx, vector<double>& qy,
    vector<double>& ua, vector<double>& va
) {
    vector<double> nu(nx * ny * nz);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i * ny;

            for (int k = 0; k < nz; k++) {
                int id = k + p * nz;

                nu[id] = num;
            }
        }
    }

    return nu;
}
