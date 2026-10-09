#include <stdexcept>

#include <utils/fs.h>
#include <utils/calc.h>

#include "proj_generator.h"

using namespace std;

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Area;

ProjGenerator::ProjGenerator(
	double latMin, double latMax,
	double lonMin, double lonMax
) {
	setLatMinMax(latMin, latMax);
	setLonMinMax(lonMin, lonMax);
}

double ProjGenerator::getLatMin() {
    return latMin;
}

double ProjGenerator::getLatMax() {
    return latMax;
}

void ProjGenerator::setLatMinMax(double latMin, double latMax) {
    if (latMin >= latMax) {
        throw runtime_error(
            format("LatMin should be < LatMax, but LatMin={} and LatMax={}", latMin, latMax)
        );
    }

    this->latMin = latMin;
    this->latMax = latMax;
}

double ProjGenerator::getLonMin() {
    return lonMin;
}

double ProjGenerator::getLonMax() {
    return lonMax;
}

void ProjGenerator::setLonMinMax(double lonMin, double lonMax) {
    if (lonMin >= lonMax) {
        throw runtime_error(
            format("LonMin should be < LonMax, but LonMin={} and LonMax={}", lonMin, lonMax)
        );
    }

    this->lonMin = lonMin;
    this->lonMax = lonMax;
}

path ProjGenerator::createDirectory(path outDir) {
    auto geo = generate();

    return FS::createDirByPath(
        outDir / path(
            format("w={}, l={}", geo.w, geo.l)
        )
    );
}

Geometry ProjGenerator::generate() {
    auto [minX, maxX, minY, maxY] = Calc::project(latMin, latMax, lonMin, lonMax);

    return Geometry{
        .l = maxX - minX,
        .w = maxY - minY
    };
}
