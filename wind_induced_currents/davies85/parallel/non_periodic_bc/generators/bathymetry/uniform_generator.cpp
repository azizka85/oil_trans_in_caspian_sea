#include <stdexcept>

#include <utils/fs.h>

#include "uniform_generator.h"

using namespace std;

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Bathymetry;

UniformGenerator::UniformGenerator(float hm) {
	setHM(hm);
}

float UniformGenerator::getHM() {
    return hm;
}

void UniformGenerator::setHM(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("HM should be > 0, but it is {}", val)
        );
    }

    hm = val;
}

path UniformGenerator::createDirectory(path outDir) {
    return FS::createDirByPath(
        outDir / path(
            format("UH, h={}", hm)
        )
    );
}

vector<float> UniformGenerator::generateH(int nx, int ny) {
    vector<float> h(nx * ny);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int id = j + i * ny;

            h[id] = hm;
        }
    }

    return h;
}
