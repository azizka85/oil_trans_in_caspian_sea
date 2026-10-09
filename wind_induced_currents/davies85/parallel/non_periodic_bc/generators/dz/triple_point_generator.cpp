#include <stdexcept>

#include <utils/fs.h>

#include "triple_point_generator.h"

using namespace std;

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::DZ;

TriplePointGenerator::TriplePointGenerator(double dzMin, double dzMax, double zm) {
    setDZMinMax(dzMin, dzMax);

    setZM(zm);
}

double TriplePointGenerator::getDZMin() {
    return dzMin;
}

double TriplePointGenerator::getDZMax() {
    return dzMax;
}

void TriplePointGenerator::setDZMinMax(double dzMin, double dzMax) {
    if (dzMin <= 0) {
        throw runtime_error(
            format("DZMin should be > 0, but it is {}", dzMin)
        );
    }

    if (dzMax <= 0) {
        throw runtime_error(
            format("DZMax should be > 0, but it is {}", dzMax)
        );
    }

    if (dzMin >= dzMax) {
        throw runtime_error(
            format("DZMin should be < DZMax, but they are {} and {}", dzMin, dzMax)
        );
    }

    this->dzMin = dzMin;
    this->dzMax = dzMax;
}

double TriplePointGenerator::getZM() {
    return zm;
}

void TriplePointGenerator::setZM(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("ZM should be > 0, but it is {}", val)
        );
    }

    zm = val;
}

path TriplePointGenerator::createDirectory(path outDir) {
    return FS::createDirByPath(
        outDir / path(
            format("TPDZ, dzMin={}, dzMax={}, zm={}", dzMin, dzMax, zm)
        )
    );
}

vector<double> TriplePointGenerator::generateDZ() {
    vector<double> dz;

    double z = 0;
    double vdz = calcDZ(z);

    dz.push_back(vdz);

    while (z < 1) {
        z += vdz;

        if (z >= 1) {
            if (z > 1) {
                z -= vdz;

                int n = dz.size();

                double dzd = 1 - z;

                if (n > 1) {
                    double dzp = dz[n - 2];

                    dz[n - 2] = (dzd + dzp) / 2;
                    dz[n - 1] = (dzd + dzp) / 2;
                }
                else {
                    dz[n - 1] = dzd;
                }
            }

            break;
        }

        vdz = calcDZ(z);
        dz.push_back(vdz);
    }

    return dz;
}

double TriplePointGenerator::calcDZ(double z) {
    double k = dzMax / (zm * zm * zm * zm / 4 - (zm + 1) * zm * zm * zm / 3 + zm * zm * zm / 2);

    return k * (z * z * z * z / 4 - (zm + 1) * z * z * z / 3 + zm * z * z / 2) + dzMin;
}
