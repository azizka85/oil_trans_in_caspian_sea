#include <stdexcept>

#include <utils/fs.h>

#include "uniform_generator.h"

using namespace std;

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::DZ;

UniformGenerator::UniformGenerator(float dzm) {
	setDZM(dzm);
}

float UniformGenerator::getDZM() {
    return dzm;
}

void UniformGenerator::setDZM(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("DZM should be > 0, but it is {}", val)
        );
    }

    dzm = val;
}

path UniformGenerator::createDirectory(path outDir) {        
    return FS::createDirByPath(
        outDir / path(
            format("UDZ, dz={}", dzm)
        )
    );
}

vector<float> UniformGenerator::generateDZ(){
    int nz = static_cast<int>(
        ceil(1. / dzm)
    );

    vector<float> dz(nz);

    for (int k = 0; k < nz; k++) {
        dz[k] = dzm;
    }

    return dz;
}
