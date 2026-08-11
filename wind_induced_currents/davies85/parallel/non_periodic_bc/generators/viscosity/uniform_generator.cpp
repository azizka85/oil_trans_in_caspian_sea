#include <stdexcept>

#include <utils/fs.h>

#include "uniform_generator.h"

using namespace std;

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Viscosity;

UniformGenerator::UniformGenerator(float num) {
	setNUM(num);
}

float UniformGenerator::getNUM() {
    return num;
}

void UniformGenerator::setNUM(float val) {
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

vector<float> UniformGenerator::generateNU(
    int nx, int ny, int nz,
    vector<float>& dz, vector<float>& h,
    vector<float>& u10, vector<float>& v10,
    vector<float>& qx, vector<float>& qy,
    vector<float>& ua, vector<float>& va
) {
    vector<float> nu(nx * ny * nz);

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
