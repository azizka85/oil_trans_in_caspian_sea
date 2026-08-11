#ifndef UTILS_CALC_H
#define UTILS_CALC_H

#include <tuple>
#include <vector>

using namespace std;

namespace Utils::Calc {
	float adjustTimeStep(float b, float t, float dt, float tMax, float dtMax, bool mult);

    float maxAbsDifference(
        int nx, int ny,
        vector<float>& u,
        vector<float>& u1
    );

    float maxAbsDifference(
        int nx, int ny, int nz,
        vector<float>& u,
        vector<float>& u1
    );

    void updateData(
        int nx, int ny,
        vector<float>& u,
        vector<float>& u1
    );

    void updateData(
        int nx, int ny, int nz,
        vector<float>& u,
        vector<float>& u1
    );

    tuple<float, float, float, float> project(
        float latMin, float latMax,
        float lonMin, float lonMax
    );
}

#endif