#include <span>

#include <stdexcept>

#include <netcdf>

#include <proj.h>

#include <boost/math/interpolators/bilinear_uniform.hpp>

#include <utils/fs.h>
#include <utils/calc.h>

#include "ecmwf_generator.h"

using namespace std;

using namespace netCDF;

using namespace boost::math::interpolators;

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Wind;

ECMWFGenerator::ECMWFGenerator(
    float latMin, float latMax,
    float lonMin, float lonMax,
    float rhoAir, float Cd,
    string filePath
) {
    setLatMinMax(latMin, latMax);
    setLonMinMax(lonMin, lonMax);

    setRhoAir(rhoAir);
    setCd(Cd);

    setFilePath(filePath);
}

float ECMWFGenerator::getLatMin() {
    return latMin;
}

float ECMWFGenerator::getLatMax() {
    return latMax;
}

void ECMWFGenerator::setLatMinMax(float latMin, float latMax) {
    if (latMin >= latMax) {
        throw runtime_error(
            format("LatMin should be < LatMax, but LatMin={} and LatMax={}", latMin, latMax)
        );
    }

    this->latMin = latMin;
    this->latMax = latMax;
}

float ECMWFGenerator::getLonMin() {
    return lonMin;
}

float ECMWFGenerator::getLonMax() {
    return lonMax;
}

void ECMWFGenerator::setLonMinMax(float lonMin, float lonMax) {
    if (lonMin >= lonMax) {
        throw runtime_error(
            format("LonMin should be < LonMax, but LonMin={} and LonMax={}", lonMin, lonMax)
        );
    }

    this->lonMin = lonMin;
    this->lonMax = lonMax;
}

float ECMWFGenerator::getRhoAir() {
    return rhoAir;
}

void ECMWFGenerator::setRhoAir(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("RhoAir should be > 0, but it is {}", val)
        );
    }

    rhoAir = val;
}

float ECMWFGenerator::getCd() {
    return Cd;
}

void ECMWFGenerator::setCd(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("Cd should be > 0, but it is {}", val)
        );
    }

    Cd = val;
}

string ECMWFGenerator::getFilePath() {
    return filePath;
}

void ECMWFGenerator::setFilePath(string val) {
    if (val.empty()) {
        throw runtime_error(
            format("filePath should not be empty, but it is {}", val)
        );
    }

    filePath = val;
}

path ECMWFGenerator::createDirectory(path outDir) {
    return FS::createDirByPath(
        outDir / path(
            format("Wind, ECMWF, lat={}-{}, lon={}-{}", latMin, latMax, lonMin, lonMax)
        )
    );
}

vector<Data> ECMWFGenerator::generate(int nx, int ny) {
    if (nx <= 1 && ny <= 1) {
        return vector<Data>();
    }

    NcFile windFile(filePath, NcFile::read);

    NcVar timeVar = windFile.getVar("valid_time");
    size_t timeSize = timeVar.getDim(0).getSize();

    vector<int64_t> times(timeSize);

    timeVar.getVar(times.data());

    NcVar latVar = windFile.getVar("latitude");
    size_t latSize = latVar.getDim(0).getSize();
    vector<double> lats(latSize);

    latVar.getVar(lats.data());

    float latStep = lats.size() > 1 ? lats[0] - lats[1] : 0;

    NcVar lonVar = windFile.getVar("longitude");
    size_t lonSize = lonVar.getDim(0).getSize();
    vector<double> lons(lonSize);

    lonVar.getVar(lons.data());

    float lonStep = lons.size() > 1 ? lons[1] - lons[0] : 0;

    NcVar u10Var = windFile.getVar("u10");
    vector<double> u10Arr(timeSize * latSize * lonSize);

    u10Var.getVar(u10Arr.data());

    NcVar v10Var = windFile.getVar("v10");
    vector<double> v10Arr(timeSize * latSize * lonSize);

    v10Var.getVar(v10Arr.data());

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

    vector<Data> res(timeSize);

    for (int k = 0; k < timeSize; k++) {
        double* u10SlicePtr = u10Arr.data() + k * latSize * lonSize;
        span<double> u10SliceSpan(u10SlicePtr, latSize * lonSize);

        auto u10Func = bilinear_uniform(
            move(u10SliceSpan),
            latSize, lonSize,
            lonStep, latStep,
            lons.front(), lats.back()
        );

        double* v10SlicePtr = v10Arr.data() + k * latSize * lonSize;
        span<double> v10SliceSpan(v10SlicePtr, latSize * lonSize);

        auto v10Func = bilinear_uniform(
            move(v10SliceSpan),
            latSize, lonSize,
            lonStep, latStep,
            lons.front(), lats.back()
        );

        vector<float> u10(nx * ny);
        vector<float> v10(nx * ny);

        vector<float> qx(nx * ny);
        vector<float> qy(nx * ny);

        int64_t refTime = timeSize > 0 ? times[0] : 0;

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                int id = j + i * ny;

                float x = minX + dx * i;
                float y = minY + dy * j;

                PJ_COORD coord = proj_coord(x, y, 0, 0);
                PJ_COORD geo = proj_trans(transformer.get(), PJ_INV, coord);

                float vu10 = 0;
                float vv10 = 0;

                float vqx = 0;
                float vqy = 0;

                if (
                    geo.lp.lam >= lats.back() && geo.lp.lam <= lats.front() &&
                    geo.lp.phi >= lons.front() && geo.lp.phi <= lons.back()
                    
                ) {
                    vu10 = u10Func(geo.lp.phi, geo.lp.lam);
                    vv10 = v10Func(geo.lp.phi, geo.lp.lam);

                    vqx = rhoAir * Cd * vu10 * sqrt(vu10 * vu10 + vv10 * vv10);
                    vqy = rhoAir * Cd * vv10 * sqrt(vu10 * vu10 + vv10 * vv10);
                }

                u10[id] = vu10;
                v10[id] = vv10;

                qx[id] = vqx;
                qy[id] = vqy;
            }
        }

        res[k].time = times[k] - refTime;

        res[k].u10 = u10;
        res[k].v10 = v10;

        res[k].qx = qx;
        res[k].qy = qy;
    }

    return res;
}
