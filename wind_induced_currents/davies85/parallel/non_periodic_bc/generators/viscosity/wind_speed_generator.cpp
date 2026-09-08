#include <stdexcept>

#include <utils/fs.h>

#include "wind_speed_generator.h"

using namespace std;

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Viscosity;

WindSpeedGenerator::WindSpeedGenerator(
    float k0, float k,
    float f, float sigma,
    float rho, float nu0
) {
    setK0(k0);
    setK(k);

    setF(f);
    setSigma(sigma);

    setRho(rho);
    setNU0(nu0);
}

float WindSpeedGenerator::getK0() {
    return k0;
}

void WindSpeedGenerator::setK0(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("K0 should be > 0, but it is {}", val)
        );
    }

    k0 = val;
}

float WindSpeedGenerator::getK() {
    return k;
}

void WindSpeedGenerator::setK(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("K should be > 0, but it is {}", val)
        );
    }

    k = val;
}

float WindSpeedGenerator::getF() {
    return f;
}

void WindSpeedGenerator::setF(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("F should be > 0, but it is {}", val)
        );
    }

    f = val;
}

float WindSpeedGenerator::getSigma() {
    return sigma;
}

void WindSpeedGenerator::setSigma(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("Sigma should be > 0, but it is {}", val)
        );
    }

    sigma = val;
}

float WindSpeedGenerator::getRho() {
    return rho;
}

void WindSpeedGenerator::setRho(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("Rho should be > 0, but it is {}", val)
        );
    }

    rho = val;
}

float WindSpeedGenerator::getNU0() {
    return nu0;
}

void WindSpeedGenerator::setNU0(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("NU0 should be > 0, but it is {}", val)
        );
    }

    nu0 = val;
}

path WindSpeedGenerator::createDirectory(path outDir) {
    return FS::createDirByPath(
        outDir / path(
            format("NU from Wind Speed, k0={}, k={}, sigma={}, nu0={}", k0, k, sigma, nu0)
        )
    );
}

vector<float> WindSpeedGenerator::generateNU(
    int nx, int ny, int nz,
    vector<float>& dz, vector<float>& h,
    vector<float>& u10, vector<float>& v10,
    vector<float>& qx, vector<float>& qy,
    vector<float>& ua, vector<float>& va
) {
    vector<float> nu(nx * ny * nz);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i * ny;

            float vqx = qx[p];
            float vqy = qy[p];

            float vq = sqrt(vqx*vqx + vqy*vqy);

            float ut = sqrt(vq / rho);
            float ht = k0 * ut / f;
            
            float vua = ua[p];
            float vva = va[p];

            float vu10 = u10[p];
            float vv10 = v10[p];
            
            float nut = k * (vua * vua + vva * vva) / sigma + nu0;
            float nus = 0.1825e-3 * pow(vu10 * vu10 + vv10 * vv10, 1.25) + nu0;

            float z = 0;

            for (int k = 0; k < nz; k++) {
                int id = k + p * nz;

                if (z < ht) {
                    nu[id] = nus - (nus - nut) * z / ht;
                }
                else {
                    nu[id] = nut;
                }

                if (k < nz - 1) {
                    z += h[p] * dz[k];
                }
            }
        }
    }

    return nu;
}
