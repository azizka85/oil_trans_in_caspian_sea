#include <stdexcept>

#include <utils/fs.h>

#include "uniform_generator.h"

using namespace std;

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Wind;

UniformGenerator::UniformGenerator(
    float u10m, float v10m,
    float qxm, float qym
) {
    setU10M(u10m);
    setV10M(v10m);

	setQXM(qxm);
	setQYM(qym);
}

float UniformGenerator::getU10M() {
    return u10m;
}

void UniformGenerator::setU10M(float val) {
    u10m = val;
}

float UniformGenerator::getV10M() {
    return v10m;
}

void UniformGenerator::setV10M(float val) {
    v10m = val;
}

float UniformGenerator::getQXM() {
    return qxm;
}

void UniformGenerator::setQXM(float val) {
    qxm = val;
}

float UniformGenerator::getQYM() {
    return qym;
}

void UniformGenerator::setQYM(float val) {
    qym = val;
}

path UniformGenerator::createDirectory(path outDir) {
    return FS::createDirByPath(
        outDir / path(
            format("UQ, qx={}, qy={}", qxm, qym)
        )
    );
}

vector<Data> UniformGenerator::generate(int nx, int ny) {
    vector<float> u10(nx * ny);
    vector<float> v10(nx * ny);

    vector<float> qx(nx * ny);
    vector<float> qy(nx * ny);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int id = j + i * ny;

            u10[id] = u10m;
            v10[id] = v10m;

            qx[id] = qxm;
            qy[id] = qym;
        }
    }

    return {
        Data {
            .time = 0,

            .u10 = u10,
            .v10 = v10,

            .qx = qx,
            .qy = qy
        }
    };
}