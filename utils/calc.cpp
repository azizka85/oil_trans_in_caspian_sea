#include <memory>

#include <proj.h>

#include "calc.h"

using namespace Utils;

float Calc::adjustTimeStep(float b, float t, float dt, float tMax, float dtMax, bool mult) {
    if (dt >= dtMax) {
        return dtMax;
    }

    if (t < tMax && t + dt >= tMax) {
        return tMax - t;
    }

    if (mult) {
        float dtp = b * dt;

        if (dtp > dtMax) {
            return dtMax;
        }
        else {
            return dtp;
        }
    }

    return dt;
}

float Calc::maxAbsDifference(int nx, int ny, vector<float>& u, vector<float>& u1) {
    float maxDiff = 0.;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i * ny;

            maxDiff = max(
                maxDiff,
                abs(u[p] - u1[p])
            );
        }
    }

    return maxDiff;
}

float Calc::maxAbsDifference(int nx, int ny, int nz, vector<float>& u, vector<float>& u1) {
    float maxDiff = 0.;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i * ny;

            for (int k = 0; k < nz; k++) {
                int id = k + p * nz;

                maxDiff = max(
                    maxDiff,
                    abs(u[id] - u1[id])
                );
            }
        }
    }

    return maxDiff;
}

void Calc::updateData(int nx, int ny, vector<float>& u, vector<float>& u1) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i * ny;

            u[p] = u1[p];
        }
    }
}

void Calc::updateData(int nx, int ny, int nz, vector<float>& u, vector<float>& u1) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i * ny;

            for (int k = 0; k < nz; k++) {
                int id = k + p * nz;

                u[id] = u1[id];
            }
        }
    }
}

tuple<float, float, float, float> Calc::project(
    float latMin, float latMax,
    float lonMin, float lonMax
) {
    unique_ptr<PJ_CONTEXT, PJ_CONTEXT* (*)(PJ_CONTEXT*)> ctx(
        proj_context_create(),
        proj_context_destroy
    );

    unique_ptr<PJ, PJ* (*)(PJ*)> transformer(
        proj_create_crs_to_crs(ctx.get(), "EPSG:4326", "EPSG:32639", nullptr),
        proj_destroy
    );

    PJ_COORD coord = proj_coord(latMin, lonMin, 0, 0);
    PJ_COORD utm = proj_trans(transformer.get(), PJ_FWD, coord);

    float minX = utm.xy.x;
    float maxX = utm.xy.x;

    float minY = utm.xy.y;
    float maxY = utm.xy.y;

    coord = proj_coord(latMin, lonMax, 0, 0);
    utm = proj_trans(transformer.get(), PJ_FWD, coord);

    minX = min(minX, (float)utm.xy.x);
    maxX = max(maxX, (float)utm.xy.x);

    minY = min(minY, (float)utm.xy.y);
    maxY = max(maxY, (float)utm.xy.y);

    coord = proj_coord(latMax, lonMin, 0, 0);
    utm = proj_trans(transformer.get(), PJ_FWD, coord);

    minX = min(minX, (float)utm.xy.x);
    maxX = max(maxX, (float)utm.xy.x);

    minY = min(minY, (float)utm.xy.y);
    maxY = max(maxY, (float)utm.xy.y);

    coord = proj_coord(latMax, lonMax, 0, 0);
    utm = proj_trans(transformer.get(), PJ_FWD, coord);

    minX = min(minX, (float)utm.xy.x);
    maxX = max(maxX, (float)utm.xy.x);

    minY = min(minY, (float)utm.xy.y);
    maxY = max(maxY, (float)utm.xy.y);

    return {
        minX, maxX,
        minY, maxY
    };
}
