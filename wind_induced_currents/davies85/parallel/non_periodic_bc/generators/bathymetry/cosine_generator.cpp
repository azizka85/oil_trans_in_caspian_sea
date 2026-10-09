#define _USE_MATH_DEFINES

#include <cmath>

#include <stdexcept>

#include <utils/fs.h>

#include "cosine_generator.h"

using namespace std;

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Bathymetry;

CosineGenerator::CosineGenerator(double hm) {
    setHM(hm);
}

double CosineGenerator::getHM() {
    return hm;
}

void CosineGenerator::setHM(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("HM should be > 0, but it is {}", val)
        );
    }

    hm = val;
}

path CosineGenerator::createDirectory(path outDir) {
    return FS::createDirByPath(
        outDir / path(
            format("UH, h={}", hm)
        )
    );
}

vector<double> CosineGenerator::generateH(int nx, int ny) {
    vector<double> h(nx * ny);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int id = j + i * ny;

            h[id] = hm * (1 - 0.8 * cos(2 * M_PI * (i / (nx - 1.) + j / (ny - 1.)))) / 2;
        }
    }

    return h;
}
