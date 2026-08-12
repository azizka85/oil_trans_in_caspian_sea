#include <stdexcept>

#include <utils/fs.h>

#include "linear_generator.h"

using namespace std;

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Viscosity;

LinearGenerator::LinearGenerator(
    float ht,
    float nus, float nut
) {
    setHT(ht);

    setNUS(nus);
    setNUT(nut);    
}

float LinearGenerator::getHT() {
    return ht;
}

void LinearGenerator::setHT(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("HT should be > 0, but it is {}", val)
        );
    }

    ht = val;
}

float LinearGenerator::getNUS() {
    return nus;
}

void LinearGenerator::setNUS(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("NUS should be > 0, but it is {}", val)
        );
    }

    nus = val;
}

float LinearGenerator::getNUT() {
    return nut;
}

void LinearGenerator::setNUT(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("NUT should be > 0, but it is {}", val)
        );
    }

    nut = val;
}

path LinearGenerator::createDirectory(path outDir) {
    return FS::createDirByPath(
        outDir / path(
            format("LNU, ht={}, nus={}, nut={}", ht, nus, nut)
        )
    );
}

vector<float> LinearGenerator::generateNU(
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

            float z = 0;

            for (int k = 0; k < nz; k++) {
                int id = k + p * nz;

                if (z < ht) {
                    nu[id] = nus - (nus - nut) * z / ht;
                }
                else {
                    nu[id] = nut;
                }

                if (k < nz - 1) {
                    z += h[p] * dz[k];
                }
            }
        }
    }

    return nu;
}
