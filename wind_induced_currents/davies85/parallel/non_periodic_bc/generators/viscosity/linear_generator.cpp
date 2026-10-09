#include <stdexcept>

#include <utils/fs.h>

#include "linear_generator.h"

using namespace std;

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Viscosity;

LinearGenerator::LinearGenerator(
    double ht,
    double nus, double nut
) {
    setHT(ht);

    setNUS(nus);
    setNUT(nut);    
}

double LinearGenerator::getHT() {
    return ht;
}

void LinearGenerator::setHT(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("HT should be > 0, but it is {}", val)
        );
    }

    ht = val;
}

double LinearGenerator::getNUS() {
    return nus;
}

void LinearGenerator::setNUS(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("NUS should be > 0, but it is {}", val)
        );
    }

    nus = val;
}

double LinearGenerator::getNUT() {
    return nut;
}

void LinearGenerator::setNUT(double val) {
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

vector<double> LinearGenerator::generateNU(
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

            double z = 0;

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
