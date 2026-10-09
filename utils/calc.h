#ifndef UTILS_CALC_H
#define UTILS_CALC_H

#include <tuple>
#include <vector>

using namespace std;

namespace Utils::Calc {
	double adjustTimeStep(double b, double t, double dt, double tMax, double dtMax, bool mult);

    double maxAbsDifference(
        int nx, int ny,
        vector<double>& u,
        vector<double>& u1
    );

    double maxAbsDifference(
        int nx, int ny, int nz,
        vector<double>& u,
        vector<double>& u1
    );

    void updateData(
        int nx, int ny,
        vector<double>& u,
        vector<double>& u1
    );

    void updateData(
        int nx, int ny, int nz,
        vector<double>& u,
        vector<double>& u1
    );

    tuple<double, double, double, double> project(
        double latMin, double latMax,
        double lonMin, double lonMax
    );
}

#endif