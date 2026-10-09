#include <stdexcept>

#include <utils/fs.h>

#include "wind_speed_generator.h"

using namespace std;

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Viscosity;

WindSpeedGenerator::WindSpeedGenerator(
    double k0, double k,
    double f, double sigma,
    double rho, double nu0
) {
    setK0(k0);
    setK(k);

    setF(f);
    setSigma(sigma);

    setRho(rho);
    setNU0(nu0);
}

double WindSpeedGenerator::getK0() {
    return k0;
}

void WindSpeedGenerator::setK0(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("K0 should be > 0, but it is {}", val)
        );
    }

    k0 = val;
}

double WindSpeedGenerator::getK() {
    return k;
}

void WindSpeedGenerator::setK(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("K should be > 0, but it is {}", val)
        );
    }

    k = val;
}

double WindSpeedGenerator::getF() {
    return f;
}

void WindSpeedGenerator::setF(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("F should be > 0, but it is {}", val)
        );
    }

    f = val;
}

double WindSpeedGenerator::getSigma() {
    return sigma;
}

void WindSpeedGenerator::setSigma(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("Sigma should be > 0, but it is {}", val)
        );
    }

    sigma = val;
}

double WindSpeedGenerator::getRho() {
    return rho;
}

void WindSpeedGenerator::setRho(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("Rho should be > 0, but it is {}", val)
        );
    }

    rho = val;
}

double WindSpeedGenerator::getNU0() {
    return nu0;
}

void WindSpeedGenerator::setNU0(double val) {
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

vector<double> WindSpeedGenerator::generateNU(
    int nx, int ny, int nz,
    vector<double>& dz, vector<double>& h,
    vector<double>& u10, vector<double>& v10,
    vector<double>& qx, vector<double>& qy,
    vector<double>& ua, vector<double>& va
) {
    vector<double> nu(nx * ny * nz);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i * ny;

            double vqx = qx[p];
            double vqy = qy[p];

            double vq = sqrt(vqx*vqx + vqy*vqy);

            double ut = sqrt(vq / rho);
            double ht = k0 * ut / f;
            
            double vua = ua[p];
            double vva = va[p];

            double vu10 = u10[p];
            double vv10 = v10[p];
            
            double nut = k * (vua * vua + vva * vva) / sigma + nu0;
            double nus = 0.1825e-3 * pow(vu10 * vu10 + vv10 * vv10, 1.25) + nu0;

            double z = 0;

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
