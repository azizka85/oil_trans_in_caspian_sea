#include "viscosity.h"

using namespace Utils;

float Viscosity::maxNU(vector<float>& nu, int nx, int ny, int nz) {
    float nuMax = nu[0];

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i * ny;

            for (int k = 0; k < nz; k++) {
                int id = k + p * nz;

                nuMax = max(nuMax, nu[id]);
            }
        }
    }

    return nuMax;
}