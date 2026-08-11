#ifndef UTILS_VISCOSITY_H
#define UTILS_VISCOSITY_H

#include <vector>

using namespace std;

namespace Utils::Viscosity {
	float maxNU(vector<float> &nu, int nx, int ny, int nz);
}

#endif