#include "bathymetry.h"

using namespace Utils;

tuple<float, float> Bathymetry::minMaxH(vector<float>& h, int nx, int ny) {
    float minH = h[0];
    float maxH = h[0];

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int id = j + i * ny;

            if (h[id] > 0) {
                minH = min(minH, h[id]);
                maxH = max(maxH, h[id]);
            }
        }
    }

    return { minH, maxH };
}