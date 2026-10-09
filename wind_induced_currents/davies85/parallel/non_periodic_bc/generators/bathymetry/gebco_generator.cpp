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
    double latMin, double latMax,
    double lonMin, double lonMax,
    double refDepth, double minDepth,
    string filePath
) {
    setLatMinMax(latMin, latMax);
    setLonMinMax(lonMin, lonMax);

    setRefDepth(refDepth);
    setMinDepth(minDepth);

    setFilePath(filePath);
}

double GEBCOGenerator::getLatMin() {
    return latMin;
}

double GEBCOGenerator::getLatMax() {
    return latMax;
}

void GEBCOGenerator::setLatMinMax(double latMin, double latMax) {
    if (latMin >= latMax) {
        throw runtime_error(
            format("LatMin should be < LatMax, but LatMin={} and LatMax={}", latMin, latMax)
        );
    }

    this->latMin = latMin;
    this->latMax = latMax;
}

double GEBCOGenerator::getLonMin() {
    return lonMin;
}

double GEBCOGenerator::getLonMax() {
    return lonMax;
}

void GEBCOGenerator::setLonMinMax(double lonMin, double lonMax) {
    if (lonMin >= lonMax) {
        throw runtime_error(
            format("LonMin should be < LonMax, but LonMin={} and LonMax={}", lonMin, lonMax)
        );
    }

    this->lonMin = lonMin;
    this->lonMax = lonMax;
}

double GEBCOGenerator::getRefDepth() {
    return refDepth;
}

void GEBCOGenerator::setRefDepth(double val) {
    refDepth = val;
}

double GEBCOGenerator::getMinDepth() {
    return minDepth;
}

void GEBCOGenerator::setMinDepth(double val) {
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

vector<double> GEBCOGenerator::generateH(int nx, int ny) {
    if (nx <= 1 && ny <= 1) {
        return vector<double>();
    }

    NcFile bathymetryFile(filePath, NcFile::read);

    NcVar latVar = bathymetryFile.getVar("lat");
    size_t latSize = latVar.getDim(0).getSize();
    vector<double> lats(latSize);

    latVar.getVar(lats.data());

    double latStep = lats.size() > 1 ? lats[1] - lats[0] : 0;

    NcVar lonVar = bathymetryFile.getVar("lon");
    size_t lonSize = lonVar.getDim(0).getSize();
    vector<double> lons(lonSize);

    lonVar.getVar(lons.data());

    double lonStep = lons.size() > 1 ? lons[1] - lons[0] : 0;

    NcVar elevVar = bathymetryFile.getVar("elevation");
    vector<double> elevations(latSize * lonSize);

    elevVar.getVar(elevations.data());

    auto [minX, maxX, minY, maxY] = Calc::project(latMin, latMax, lonMin, lonMax);

    double dx = (maxX - minX) / (nx - 1);
    double dy = (maxY - minY) / (ny - 1);

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

    vector<double> depths(nx * ny);

    double epsilon = 0.003;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i * ny;

            double x = minX + dx * i;
            double y = minY + dy * j;

            PJ_COORD coord = proj_coord(x, y, 0, 0);
            PJ_COORD geo = proj_trans(transformer.get(), PJ_INV, coord);

            if (
                geo.lp.phi >= 55.9964 - epsilon && geo.lp.phi <= 55.9964 + epsilon && 
                geo.lp.lam >= 46.0642 - epsilon && geo.lp.lam <= 46.0642 + epsilon
            ) {
                bool ok = true;
            }

            if (
                geo.lp.lam >= lats.front() + epsilon && geo.lp.lam <= lats.back() - epsilon &&
                geo.lp.phi >= lons.front() + epsilon && geo.lp.phi <= lons.back() - epsilon
            ) {
                double elevation = interpFunc(geo.lp.phi, geo.lp.lam);

                if (elevation < refDepth - minDepth) {
                    depths[p] = abs(elevation - refDepth);
                }
            }
        }
    }

    return depths;
}
