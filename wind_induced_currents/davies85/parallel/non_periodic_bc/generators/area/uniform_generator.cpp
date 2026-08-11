#include <stdexcept>

#include <utils/fs.h>

#include "uniform_generator.h"

using namespace std;

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Area;

UniformGenerator::UniformGenerator(float w, float l) {
    setW(w);
    setL(l);
}

float UniformGenerator::getW() {
    return w;
}

void UniformGenerator::setW(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("W should be > 0, but it is {}", val)
        );
    }

    w = val;
}

float UniformGenerator::getL() {
    return l;
}

void UniformGenerator::setL(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("L should be > 0, but it is {}", val)
        );
    }

    l = val;
}

path UniformGenerator::createDirectory(path outDir) {
    return FS::createDirByPath(
        outDir / path(
            format("w={}, l={}", w, l)
        )
    );
}

Geometry UniformGenerator::generate() {
    return Geometry{
        .l = l,
        .w = w
    };
}
