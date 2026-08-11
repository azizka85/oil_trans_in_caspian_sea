#ifndef UTILS_BATHYMETRY_H
#define UTILS_BATHYMETRY_H

#include <tuple>
#include <vector>

using namespace std;

namespace Utils::Bathymetry {
	tuple<float, float> minMaxH(vector<float> &h, int nx, int ny);
}

#endif