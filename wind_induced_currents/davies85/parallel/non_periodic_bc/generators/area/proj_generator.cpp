#include <stdexcept>

#include <utils/fs.h>
#include <utils/calc.h>

#include "proj_generator.h"

using namespace std;

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Area;

ProjGenerator::ProjGenerator(
	float latMin, float latMax,
	float lonMin, float lonMax
) {
	setLatMinMax(latMin, latMax);
	setLonMinMax(lonMin, lonMax);
}

float ProjGenerator::getLatMin() {
    return latMin;
}

float ProjGenerator::getLatMax() {
    return latMax;
}

void ProjGenerator::setLatMinMax(float latMin, float latMax) {
    if (latMin >= latMax) {
        throw runtime_error(
            format("LatMin should be < LatMax, but LatMin={} and LatMax={}", latMin, latMax)
        );
    }

    this->latMin = latMin;
    this->latMax = latMax;
}

float ProjGenerator::getLonMin() {
    return lonMin;
}

float ProjGenerator::getLonMax() {
    return lonMax;
}

void ProjGenerator::setLonMinMax(float lonMin, float lonMax) {
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
