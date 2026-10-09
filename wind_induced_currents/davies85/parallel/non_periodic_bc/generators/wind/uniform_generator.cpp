#include <stdexcept>

#include <utils/fs.h>

#include "uniform_generator.h"

using namespace std;

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Wind;

UniformGenerator::UniformGenerator(
    double u10m, double v10m,
    double qxm, double qym
) {
    setU10M(u10m);
    setV10M(v10m);

	setQXM(qxm);
	setQYM(qym);
}

double UniformGenerator::getU10M() {
    return u10m;
}

void UniformGenerator::setU10M(double val) {
    u10m = val;
}

double UniformGenerator::getV10M() {
    return v10m;
}

void UniformGenerator::setV10M(double val) {
    v10m = val;
}

double UniformGenerator::getQXM() {
    return qxm;
}

void UniformGenerator::setQXM(double val) {
    qxm = val;
}

double UniformGenerator::getQYM() {
    return qym;
}

void UniformGenerator::setQYM(double val) {
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
    vector<double> u10(nx * ny);
    vector<double> v10(nx * ny);

    vector<double> qx(nx * ny);
    vector<double> qy(nx * ny);

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