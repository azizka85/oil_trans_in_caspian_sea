#include <stdexcept>

#include <netcdf>

#include <proj.h>

#include <boost/math/interpolators/bilinear_uniform.hpp>

#include <utils/fs.h>
#include <utils/calc.h>

#include "gebco_generator.h"

using namespace std;

using namespace netCDF;

using namespace boost::math::interpolators;

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Bathymetry;

GEBCOGenerator::GEBCOGenerator(
    float latMin, float latMax,
    float lonMin, float lonMax,
    float refDepth, float minDepth,
    string filePath
) {
    setLatMinMax(latMin, latMax);
    setLonMinMax(lonMin, lonMax);

    setRefDepth(refDepth);
    setMinDepth(minDepth);

    setFilePath(filePath);
}

float GEBCOGenerator::getLatMin() {
    return latMin;
}

float GEBCOGenerator::getLatMax() {
    return latMax;
}

void GEBCOGenerator::setLatMinMax(float latMin, float latMax) {
    if (latMin >= latMax) {
        throw runtime_error(
            format("LatMin should be < LatMax, but LatMin={} and LatMax={}", latMin, latMax)
        );
    }

    this->latMin = latMin;
    this->latMax = latMax;
}

float GEBCOGenerator::getLonMin() {
    return lonMin;
}

float GEBCOGenerator::getLonMax() {
    return lonMax;
}

void GEBCOGenerator::setLonMinMax(float lonMin, float lonMax) {
    if (lonMin >= lonMax) {
        throw runtime_error(
            format("LonMin should be < LonMax, but LonMin={} and LonMax={}", lonMin, lonMax)
        );
    }

    this->lonMin = lonMin;
    this->lonMax = lonMax;
}

float GEBCOGenerator::getRefDepth() {
    return refDepth;
}

void GEBCOGenerator::setRefDepth(float val) {
    refDepth = val;
}

float GEBCOGenerator::getMinDepth() {
    return minDepth;
}

void GEBCOGenerator::setMinDepth(float val) {
    if (val < 0) {
        throw runtime_error(
            format("MinDepth should be >= 0, but it is {}", val)
        );
    }

    minDepth = val;
}

string GEBCOGenerator::getFilePath() {
    return filePath;
}

void GEBCOGenerator::setFilePath(string val) {
    if (val.empty()) {
        throw runtime_error(
            format("filePath should not be empty, but it is {}", val)
        );
    }

    filePath = val;
}

path GEBCOGenerator::createDirectory(path outDir) {
    return FS::createDirByPath(
        outDir / path(
            format("Bathymetry, GEBCO, lat={}-{}, lon={}-{}", latMin, latMax, lonMin, lonMax)
        )
    );
}

vector<float> GEBCOGenerator::generateH(int nx, int ny) {
    if (nx <= 1 && ny <= 1) {
        return vector<float>();
    }

    NcFile bathymetryFile(filePath, NcFile::read);

    NcVar latVar = bathymetryFile.getVar("lat");
    size_t latSize = latVar.getDim(0).getSize();
    vector<float> lats(latSize);

    latVar.getVar(lats.data());

    float latStep = lats.size() > 1 ? lats[1] - lats[0] : 0;

    NcVar lonVar = bathymetryFile.getVar("lon");
    size_t lonSize = lonVar.getDim(0).getSize();
    vector<float> lons(lonSize);

    lonVar.getVar(lons.data());

    float lonStep = lons.size() > 1 ? lons[1] - lons[0] : 0;

    NcVar elevVar = bathymetryFile.getVar("elevation");
    vector<float> elevations(latSize * lonSize);

    elevVar.getVar(elevations.data());

    auto [minX, maxX, minY, maxY] = Calc::project(latMin, latMax, lonMin, lonMax);

    float dx = (maxX - minX) / (nx - 1);
    float dy = (maxY - minY) / (ny - 1);

    unique_ptr<PJ_CONTEXT, PJ_CONTEXT* (*)(PJ_CONTEXT*)> ctx(
        proj_context_create(),
        proj_context_destroy
    );

    unique_ptr<PJ, PJ* (*)(PJ*)> transformer(
        proj_create_crs_to_crs(ctx.get(), "EPSG:4326", "EPSG:32639", nullptr),
        proj_destroy
    );

    auto interpFunc = bilinear_uniform(
        move(elevations),
        latSize, lonSize,
        lonStep, latStep,
        lons.front(), lats.front()
    );

    vector<float> depths(nx * ny);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i * ny;

            float x = minX + dx * i;
            float y = minY + dy * j;

            PJ_COORD coord = proj_coord(x, y, 0, 0);
            PJ_COORD geo = proj_trans(transformer.get(), PJ_INV, coord);

            if (
                geo.lp.lam >= lats.front() && geo.lp.lam <= lats.back() &&
                geo.lp.phi >= lons.front() && geo.lp.phi <= lons.back()
            ) {
                float elevation = interpFunc(geo.lp.phi, geo.lp.lam);

                if (elevation < refDepth - minDepth) {
                    depths[p] = abs(elevation - refDepth);
                }
            }
        }
    }

    return depths;
}
