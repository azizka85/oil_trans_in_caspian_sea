#define _USE_MATH_DEFINES

#include <format>

#include <cmath>

#include <fstream>
#include <iostream>

#include <chrono>

#include <stdexcept>

#include <slae/direct/tridiagonal.h>

#include "solver.h"

using namespace std::chrono;

using namespace SLAE::Direct;

using namespace WindInducedCurrents::Davies85::VariableParameters;

Solver::Solver(
    double f, double b,
    double hm, GenerateH hg,
    double w, double l,
    double g, double rho, double kb,
    double num, GenerateNU nug,
    double qxm, double qym, GenerateQ qg,
    double dx, double dy,
    double dzm, GenerateDZ dzg,
    double endTime, double outputTimeStep, string dir)
{
    setF(f);
    setB(b);

    setHM(hm);
    setHG(hg);

    setW(w);
    setL(l);

    setG(g);
    setRHO(rho);
    setKB(kb);

    setNUM(num);
    setNUG(nug);
    
    setQXM(qxm);
    setQYM(qym);
    setQG(qg);

    setDX(dx);
    setDY(dy);

    setDZM(dzm);
    setDZG(dzg);

    setEndTime(endTime);
    setOutputTimeStep(outputTimeStep);
    setDir(dir);
}

double Solver::getF() {
    return f;
}

void Solver::setF(double val) {
    f = val;
}

double Solver::getB() {
    return b;
}

void Solver::setB(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("B should be > 0, but it is {}", val)
        );
    }

    b = val;
}

double Solver::getHM() {
    return hm;
}

void Solver::setHM(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("HM should be > 0, but it is {}", val)
        );
    }

    hm = val;
}

GenerateH Solver::getHG() {
    return hg;
}

void Solver::setHG(GenerateH val) {
    hg = val;
}

vector<vector<double>> Solver::getH() {
    return h;
}

double Solver::getW() {
    return w;
}

void Solver::setW(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("W should be > 0, but it is {}", val)
        );
    }

    w = val;
}

double Solver::getL() {
    return l;
}

void Solver::setL(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("L should be > 0, but it is {}", val)
        );
    }

    l = val;
}

double Solver::getG() {
    return g;
}

void Solver::setG(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("G should be > 0, but it is {}", val)
        );
    }

    g = val;
}

double Solver::getRHO() {
    return rho;
}

void Solver::setRHO(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("RHO should be > 0, but it is {}", val)
        );
    }

    rho = val;
}

double Solver::getKB() {
    return kb;
}

void Solver::setKB(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("KB should be > 0, but it is {}", val)
        );
    }

    kb = val;
}

double Solver::getNUM() {
    return num;
}

void Solver::setNUM(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("NUM should be > 0, but it is {}", val)
        );
    }

    num = val;
}

GenerateNU Solver::getNUG() {
    return nug;
}

void Solver::setNUG(GenerateNU val) {
    nug = val;
}

vector<vector<vector<double>>> Solver::getNU() {
    return nu;
}

double Solver::getQXM() {
    return qxm;
}

void Solver::setQXM(double val) {
    qxm = val;
}

double Solver::getQYM() {
    return qym;
}

void Solver::setQYM(double val) {
    qym = val;
}

GenerateQ Solver::getQG() {
    return qg;
}

void Solver::setQG(GenerateQ val) {
    qg = val;
}

vector<vector<double>> Solver::getQX() {
    return qx;
}

vector<vector<double>> Solver::getQY() {
    return qy;
}

double Solver::getDX() {
    return dx;
}

void Solver::setDX(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("DX should be > 0, but it is {}", val)
        );
    }

    dx = val;
}

double Solver::getDY() {
    return dy;
}

void Solver::setDY(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("DY should be > 0, but it is {}", val)
        );
    }

    dy = val;
}

double Solver::getDZM() {
    return dzm;
}

void Solver::setDZM(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("DZM should be > 0, but it is {}", val)
        );
    }

    dzm = val;
}

GenerateDZ Solver::getDZG() {
    return dzg;
}

void Solver::setDZG(GenerateDZ val) {
    dzg = val;
}

vector<double> Solver::getDZ() {
    return dz;
}

double Solver::getEndTime() {
    return endTime;
}

void Solver::setEndTime(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("endTime should be > 0, but it is {}", val)
        );
    }

    endTime = val;
}

double Solver::getOutputTimeStep() {
    return outputTimeStep;
}

void Solver::setOutputTimeStep(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("outputTimeStep should be > 0, but it is {}", val)
        );
    }

    outputTimeStep = val;
}

string Solver::getDir() {
    return dir;
}

void Solver::setDir(string val) {
    if (val.empty()) {
        throw runtime_error(
            format("dir should not be empty, but it is {}", val)
        );
    }

    dir = val;
}

void Solver::generateH(int nx, int ny) {
    switch(hg) {
        case GenerateH::UniformH:
            generateUniformH(nx, ny);
            break;
        case GenerateH::CosineH:
            generateCosineH(nx, ny);
            break;
    }
}

void Solver::generateUniformH(int nx, int ny) {
    h = vector<vector<double>>(nx, vector<double>(ny));

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            h[i][j] = hm;
        }
    }
}

void Solver::generateCosineH(int nx, int ny) {
    h = vector<vector<double>>(nx, vector<double>(ny));

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            h[i][j] = hm*(1 - 0.8*cos(2*M_PI*(i/(nx-1.) + j/(ny-1.))))/2;
        }
    }
}

tuple<double, double> Solver::minMaxH(int nx, int ny) {
    double minH = h[0][0];
    double maxH = h[0][0];

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            minH = min(minH, h[i][j]);
            maxH = max(maxH, h[i][j]);
        }
    }

    return {minH, maxH};
}

void Solver::generateNU(
    int nx, int ny, int nz, 
    vector<vector<double>> &ua, 
    vector<vector<double>> &va, 
    vector<vector<double>> &z, 
    vector<vector<vector<double>>> &uf, 
    vector<vector<vector<double>>> &vf
) {
    switch(nug) {
        case GenerateNU::UniformNU:
            generateUniformNU(nx, ny, nz);
            break;
        case GenerateNU::FromWindSpeedNU:
            generateNUFromWindSpeed(nx, ny, nz);
            break;
    }
}

void Solver::updateNU(
    double t, int n, int nx, int ny, int nz, 
    vector<vector<double>> &ua, 
    vector<vector<double>> &va, 
    vector<vector<double>> &z, 
    vector<vector<vector<double>>> &uf, 
    vector<vector<vector<double>>> &vf
) {
    switch(nug) {
        case GenerateNU::UniformNU:
            updateUniformNU(nx, ny, nz);
            break;
        case GenerateNU::FromWindSpeedNU:
            updateNUFromWindSpeed(nx, ny, nz);
            break;
    }
}

void Solver::generateUniformNU(int nx, int ny, int nz)
{
    nu = vector<vector<vector<double>>>(
        nx, 
        vector<vector<double>>(
            ny, 
            vector<double>(nz)
        )
    );

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                nu[i][j][k] = num;
            }
        }
    }
}

void Solver::updateUniformNU(int nx, int ny, int nz) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                nu[i][j][k] = num;
            }
        }
    }
}

void Solver::generateNUFromWindSpeed(int nx, int ny, int nz) {
    generateUniformNU(nx, ny, nz);
}

void Solver::updateNUFromWindSpeed(int nx, int ny, int nz) {
    updateUniformNU(nx, ny, nz);
}

double Solver::maxNU(int nx, int ny, int nz) {
    double nuMax = nu[0][0][0];

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                nuMax = max(nuMax, nu[i][j][k]);
            }
        }
    }

    return nuMax;
}

void Solver::generateQX(
    int nx, int ny, int nz, 
    vector<vector<double>> &ua, 
    vector<vector<double>> &va, 
    vector<vector<double>> &z, 
    vector<vector<vector<double>>> &uf, 
    vector<vector<vector<double>>> &vf
) {
    switch (qg) {
        case GenerateQ::UniformQ:
            generateUniformQX(nx, ny);
            break;
        case GenerateQ::FromWindSpeedQ:
            generateQXFromWindSpeed(nx, ny);
            break;
    }
}

void Solver::updateQX(
    double t, int n, int nx, int ny, int nz, 
    vector<vector<double>> &ua, 
    vector<vector<double>> &va, 
    vector<vector<double>> &z, 
    vector<vector<vector<double>>> &uf, 
    vector<vector<vector<double>>> &vf
) {
    switch (qg) {
        case GenerateQ::UniformQ:
            updateUniformQX(nx, ny);
            break;
        case GenerateQ::FromWindSpeedQ:
            updateQXFromWindSpeed(nx, ny);
            break;
    }
}

void Solver::generateUniformQX(int nx, int ny) {
    qx = vector<vector<double>>(nx, vector<double>(ny));

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            qx[i][j] = qxm;
        }
    }
}

void Solver::updateUniformQX(int nx, int ny) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            qx[i][j] = qxm;
        }
    }
}

void Solver::generateQXFromWindSpeed(int nx, int ny) {
    generateUniformQX(nx, ny);
}

void Solver::updateQXFromWindSpeed(int nx, int ny) {
    updateUniformQX(nx, ny);
}

void Solver::generateQY(
    int nx, int ny, int nz, 
    vector<vector<double>> &ua, 
    vector<vector<double>> &va, 
    vector<vector<double>> &z, 
    vector<vector<vector<double>>> &uf, 
    vector<vector<vector<double>>> &vf
) {
    switch (qg) {
        case GenerateQ::UniformQ:
            generateUniformQY(nx, ny);
            break;
        case GenerateQ::FromWindSpeedQ:
            generateQYFromWindSpeed(nx, ny);
            break;
    }
}

void Solver::updateQY(
    double t, int n, int nx, int ny, int nz, 
    vector<vector<double>> &ua, 
    vector<vector<double>> &va, 
    vector<vector<double>> &z, 
    vector<vector<vector<double>>> &uf, 
    vector<vector<vector<double>>> &vf
) {
    switch (qg) {
        case GenerateQ::UniformQ:
            updateUniformQY(nx, ny);
            break;
        case GenerateQ::FromWindSpeedQ:
            updateQYFromWindSpeed(nx, ny);
            break;
    }
}

void Solver::generateUniformQY(int nx, int ny) {
    qy = vector<vector<double>>(nx, vector<double>(ny));

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            qy[i][j] = qym;
        }
    }
}

void Solver::updateUniformQY(int nx, int ny) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            qy[i][j] = qym;
        }
    }
}

void Solver::generateQYFromWindSpeed(int nx, int ny) {
    generateUniformQY(nx, ny);
}

void Solver::updateQYFromWindSpeed(int nx, int ny) {
    updateUniformQY(nx, ny);
}

void Solver::generateDZ() {
    switch (dzg) {
        case GenerateDZ::UniformDZ:
            generateUniformDZ();
            break;
        case GenerateDZ::ParabolicDZ:
            generateParabolicDZ();
            break;
    }
}

void Solver::generateUniformDZ() {
    int nz = static_cast<int>(
        ceil(1./dzm)
    );

    dz = vector<double>(nz);

    for (int k = 0; k < nz; k++) {
        dz[k] = dzm;
    }
}

void Solver::generateParabolicDZ() {
    double z = 0;
    double dzf = calcParabolicDZFactor(z);

    dz.push_back(dzm*dzf);    

    while (z < 1) {
        z += dzm * dzf;

        if (z >= 1) {
            if (z > 1) {
                z -= dzm * dzf;

                int n = dz.size();

                double dzd = 1 - z;

                if (n > 1) {
                    double dzp = dz[n-2];

                    dz[n-2] = (dzd + dzp) / 2;
                    dz[n-1] = (dzd + dzp) / 2;
                } else {
                    dz[n-1] = dzd;
                }
            }

            break;
        }

        dzf = calcParabolicDZFactor(z);
        dz.push_back(dzm*dzf);
    }
}

double Solver::calcParabolicDZFactor(double z) {
    return -2.96*z*z + 3.44*z + 0.02;
}

path Solver::createDirectory() {
    auto hStr = format("UH, h={}", hm);

    if (hg == GenerateH::CosineH) {
        hStr = format("CH, h={}", hm);
    }

    auto nuStr = format("UNU, nu={}", num);    
    auto qStr = format("UQ, qx={}, qy={}", qxm, qym);
    
    auto dzStr = format("UDZ, dz={}", dzm);

    if (dzg == GenerateDZ::ParabolicDZ) {
        dzStr = format("PDZ, dz={}", dzm);
    }

    auto dirPath = path(
        format("{}/f={}, g={}, rho={}, kb={}/{}/{}/{}/dx={}, dy={}/{}", dir, f, g, rho, kb, hStr, nuStr, qStr, dx, dy, dzStr)
    );

    create_directories(dirPath);    

    auto surfacePath = dirPath / path("surface");

    create_directories(surfacePath);

    auto volumePath = dirPath / path("volume");

    create_directories(volumePath);

    auto viscosityPath = dirPath / path("viscosity");

    create_directories(viscosityPath);

    return dirPath;
}

void Solver::setInitialCondition(
    int nx, int ny, int nz, 
    vector<vector<vector<double>>>& uf, 
    vector<vector<vector<double>>>& vf, 
    vector<vector<double>>& ua,
    vector<vector<double>>& va,
    vector<vector<double>>& z
) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                uf[i][j][k] = 0;
                vf[i][j][k] = 0;
            }

            ua[i][j] = 0;
            va[i][j] = 0;
            z[i][j] = 0;
        }
    }
}

void Solver::calculateResidualElements(
    double dt, int nx, int ny, int nz,    
    vector<vector<double>>& ua,
    vector<vector<double>>& u1a,
    vector<vector<double>>& va,    
    vector<vector<double>>& v1a,
    vector<vector<vector<double>>>& uf, 
    vector<vector<vector<double>>>& vf,
    vector<vector<vector<double>>>& ud, 
    vector<vector<vector<double>>>& vd
) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 1; k < nz-1; k++) {
                double u = uf[i][j][k] + ua[i][j];
                double v = vf[i][j][k] + va[i][j];

                ud[i][j][k] = uf[i][j][k] + f*vf[i][j][k]*dt + dt*(
                    nu[i][j][k+1]*(uf[i][j][k+1] - uf[i][j][k])/(dz[k+1] + dz[k]) - 
                    nu[i][j][k]*(uf[i][j][k] - uf[i][j][k-1])/(dz[k] + dz[k-1])
                )/h[i][j]/h[i][j]/dz[k] + kb*u*dt/h[i][j] - qx[i][j]*dt/rho/h[i][j];   

                vd[i][j][k] = vf[i][j][k] - f*uf[i][j][k]*dt + dt*(
                    nu[i][j][k+1]*(vf[i][j][k+1] - vf[i][j][k])/(dz[k+1] + dz[k]) - 
                    nu[i][j][k]*(vf[i][j][k] - vf[i][j][k-1])/(dz[k] + dz[k-1])
                )/h[i][j]/h[i][j]/dz[k] + kb*v*dt/h[i][j] - qy[i][j]*dt/rho/h[i][j];   
            }

            double um1 = uf[i][j][0] + h[i][j]*qx[i][j]*dz[0]/rho/nu[i][j][0];
            double un = (
                uf[i][j][nz-1] - kb*h[i][j]*dz[nz-1]*uf[i][j][nz-1]/2/nu[i][j][nz] 
                - kb*h[i][j]*dz[nz-1]*ua[i][j]/nu[i][j][nz]
            )/(1 + kb*h[i][j]*dz[nz-1]/2/nu[i][j][nz]);

            double vm1 = vf[i][j][0] + h[i][j]*qy[i][j]*dz[0]/rho/nu[i][j][0];
            double vn = (
                vf[i][j][nz-1] - kb*h[i][j]*dz[nz-1]*vf[i][j][nz-1]/2/nu[i][j][nz] 
                - kb*h[i][j]*dz[nz-1]*va[i][j]/nu[i][j][nz]
            )/(1 + kb*h[i][j]*dz[nz-1]/2/nu[i][j][nz]);

            ud[i][j][0] = uf[i][j][0] + f*vf[i][j][0]*dt + dt*(
                nu[i][j][1]*(uf[i][j][1] - uf[i][j][0])/(dz[1] + dz[0]) - 
                nu[i][j][0]*(uf[i][j][0] - um1)/2/dz[0]
            )/h[i][j]/h[i][j]/dz[0] + kb*(uf[i][j][0] + ua[i][j])*dt/h[i][j] - qx[i][j]*dt/rho/h[i][j]
            + dt*qx[i][j]/rho/h[i][j]/2/dz[0];

            ud[i][j][nz-1] = uf[i][j][nz-1] + f*vf[i][j][nz-1]*dt + dt*(
                nu[i][j][nz]*(un - uf[i][j][nz-1])/2/dz[nz-1] - 
                nu[i][j][nz-1]*(uf[i][j][nz-1] - uf[i][j][nz-2])/(dz[nz-1] + dz[nz-2])
            )/h[i][j]/h[i][j]/dz[nz-1] + kb*(uf[i][j][nz-1] + ua[i][j])*dt/h[i][j] - qx[i][j]*dt/rho/h[i][j]
            - dt*kb*ua[i][j]/h[i][j]/2/dz[nz-1]/(1 + kb*h[i][j]*dz[nz-1]/2/nu[i][j][nz]);

            vd[i][j][0] = vf[i][j][0] - f*uf[i][j][0]*dt + dt*(
                nu[i][j][1]*(vf[i][j][1] - vf[i][j][0])/(dz[1] + dz[0]) - 
                nu[i][j][0]*(vf[i][j][0] - vm1)/2/dz[0]
            )/h[i][j]/h[i][j]/dz[0] + kb*(vf[i][j][0] + va[i][j])*dt/h[i][j] - qy[i][j]*dt/rho/h[i][j]
            + dt*qy[i][j]/rho/h[i][j]/2/dz[0];

            vd[i][j][nz-1] = vf[i][j][nz-1] - f*uf[i][j][nz-1]*dt + dt*(
                nu[i][j][nz]*(vn - vf[i][j][nz-1])/2/dz[nz-1] - 
                nu[i][j][nz-1]*(vf[i][j][nz-1] - vf[i][j][nz-2])/(dz[nz-1] + dz[nz-2])
            )/h[i][j]/h[i][j]/dz[nz-1] + kb*(vf[i][j][nz-1] + va[i][j])*dt/h[i][j] - qy[i][j]*dt/rho/h[i][j]
            - dt*kb*va[i][j]/h[i][j]/2/dz[nz-1]/(1 + kb*h[i][j]*dz[nz-1]/2/nu[i][j][nz]);
        }
    }
}

void Solver::createTridiagonalMatrix(
    int i, int j, 
    int nz, double dt, 
    vector<double> &al, 
    vector<double> &ac, 
    vector<double> &ar
) {
    ac[0] = 1 + dt*nu[i][j][1]/h[i][j]/h[i][j]/dz[0]/(dz[1] + dz[0]);
    ar[0] = -dt*nu[i][j][1]/h[i][j]/h[i][j]/dz[0]/(dz[1] + dz[0]);

    for (int k = 1; k < nz-1; k++) {
        al[k] = -dt*nu[i][j][k]/h[i][j]/h[i][j]/dz[k]/(dz[k] + dz[k-1]);
        ac[k] = 1 + dt*(
            nu[i][j][k+1]/(dz[k+1] + dz[k]) +
            nu[i][j][k]/(dz[k] + dz[k-1])
        )/h[i][j]/h[i][j]/dz[k];
        ar[k] = -dt*nu[i][j][k+1]/h[i][j]/h[i][j]/dz[k]/(dz[k+1] + dz[k]);
    }

    al[nz-1] = -dt*nu[i][j][nz-1]/h[i][j]/h[i][j]/dz[nz-1]/(dz[nz-1] + dz[nz-2]);
    ac[nz-1] = 1 + dt*(
        kb*h[i][j]*dz[nz-1]/(1 + kb*h[i][j]*dz[nz-1]/2/nu[i][j][nz])/2/dz[nz-1] +
        nu[i][j][nz-1]/(dz[nz-1] + dz[nz-2])
    )/h[i][j]/h[i][j]/dz[nz-1];
}

double Solver::adjustTimeStep(double t, double dt, double tMax, double dtMax, bool mult) {
    if (dt >= dtMax) {
        return dtMax;
    }

    if (t < tMax && t + dt >= tMax) {
        return tMax - t;
    }

    if (mult) {
        return b * dt;
    }

    return dt;
}

double Solver::maxAbsDifference(int nx, int ny, vector<vector<double>> &u, vector<vector<double>> &u1) {
    double maxDiff = 0.;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            maxDiff = max(
                maxDiff,
                abs(u[i][j] - u1[i][j])
            );
        }
    }

    return maxDiff;
}

double Solver::maxAbsDifference(int nx, int ny, int nz, vector<vector<vector<double>>> &u, vector<vector<vector<double>>> &u1) {
    double maxDiff = 0.;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                maxDiff = max(
                    maxDiff,
                    abs(u[i][j][k] - u1[i][j][k])
                );
            }
        }
    }

    return maxDiff;
}

void Solver::updateData(int nx, int ny, vector<vector<double>> &u, vector<vector<double>> &u1) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            u[i][j] = u1[i][j];
        }
    }
}

void Solver::updateData(int nx, int ny, int nz, vector<vector<vector<double>>> &u, vector<vector<vector<double>>> &u1) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                u[i][j][k] = u1[i][j][k];
            }
        }
    }
}

void Solver::writeHeights(int nx, int ny, path outDir) {
    auto filePath = path("heights.vtk");

    filePath = outDir / filePath;

    ofstream file(filePath);

    if (file.bad()) {
        throw runtime_error(
            format("Failed to open file at: {}", filePath.string())
        );
    }

    file << "# vtk DataFile Version 3.0" << endl;
    file << format("TIME {:.3f}", 0.) << endl;
    file << "ASCII" << endl;
    file << "DATASET STRUCTURED_GRID" << endl;
    file << format("DIMENSIONS {} {} 1", nx, ny) << endl;    
    file << format("POINTS {} float", nx*ny) << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double x = dx * i;
            double y = dy * j;

            file << format("{:.3f} {:.3f} 0.0", x, y) << endl;
        }
    }

    file << "FIELD FieldData 1" << endl;
    file << "Time 1 1 float" << endl;
    file << format("{:.3f}", 0.) << endl;
    file << format("POINT_DATA {}", nx*ny) << endl;

    file << "SCALARS h float" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << format("{:.3f}", h[i][j]) << endl;
        }
    }
}

void Solver::writeSurfaceData(
    double t, int m, 
    int nx, int ny, 
    vector<vector<double>> &ua, 
    vector<vector<double>> &va, 
    vector<vector<double>> &z,
    path outDir
) {
    auto filePath = path(
        format("data.{:03}.vtk", m)
    );

    filePath = outDir / path("surface") / filePath;

    ofstream file(filePath);

    if (file.bad()) {
        throw runtime_error(
            format("Failed to open file at: {}", filePath.string())
        );
    }

    file << "# vtk DataFile Version 3.0" << endl;
    file << format("TIME {:.3f}", t) << endl;
    file << "ASCII" << endl;
    file << "DATASET STRUCTURED_GRID" << endl;
    file << format("DIMENSIONS {} {} 1", nx, ny) << endl;    
    file << format("POINTS {} float", nx*ny) << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double x = dx * i;
            double y = dy * j;

            file << format("{:.3f} {:.3f} 0.0", x, y) << endl;
        }
    }

    file << "FIELD FieldData 1" << endl;
    file << "Time 1 1 float" << endl;
    file << format("{:.3f}", t) << endl;
    file << format("POINT_DATA {}", nx*ny) << endl;

    file << "VECTORS VA float" << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << format("{:.3f} {:.3f} 0", ua[i][j], va[i][j]) << endl;
        }
    }

    file << "VECTORS Q float" << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << format("{:.3f} {:.3f} 0", qx[i][j], qy[i][j]) << endl;
        }
    }

    file << "SCALARS z float" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << format("{:.7f}", z[i][j]) << endl;
        }
    }
}

void Solver::writeVolumeData(
    double t, int m, 
    int nx, int ny, int nz, 
    vector<vector<double>> &ua, 
    vector<vector<double>> &va, 
    vector<vector<vector<double>>> &uf, 
    vector<vector<vector<double>>> &vf, 
    path outDir
) {
    auto filePath = path(
        format("data.{:03}.vtk", m)
    );

    filePath = outDir / path("volume") / filePath;

    ofstream file(filePath);

    if (file.bad()) {
        throw runtime_error(
            format("Failed to open file at: {}", filePath.string())
        );
    }

    file << "# vtk DataFile Version 3.0" << endl;
    file << format("TIME {:.3f}", t) << endl;
    file << "ASCII" << endl;
    file << "DATASET STRUCTURED_GRID" << endl;
    file << format("DIMENSIONS {} {} {}", nx, ny, nz) << endl;    
    file << format("POINTS {} float", nx*ny*nz) << endl;    

    vector<vector<double>> z(nx, vector<double>(ny, 0));

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                double x = dx * i;
                double y = dy * j;

                z[i][j] += h[i][j]*dz[k] / 2;
                
                file << format("{:.3f} {:.3f} {:.3f}", x, y, z[i][j]) << endl;

                z[i][j] += h[i][j]*dz[k] / 2;
            }
        }        
    }

    file << "FIELD FieldData 1" << endl;
    file << "Time 1 1 float" << endl;
    file << format("{:.3f}", t) << endl;
    file << format("POINT_DATA {}", nx*ny*nz) << endl;

    file << "VECTORS V float" << endl;

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                file << format("{:.3f} {:.3f} 0", uf[i][j][k] + ua[i][j], vf[i][j][k] + va[i][j]) << endl;
            }
        } 
    }    
}

void Solver::writeViscosity(
    double t, int m, 
    int nx, int ny, int nz, 
    path outDir
) {
    auto filePath = path(
        format("data.{:03}.vtk", m)
    );

    filePath = outDir / path("viscosity") / filePath;

    ofstream file(filePath);

    if (file.bad()) {
        throw runtime_error(
            format("Failed to open file at: {}", filePath.string())
        );
    }

    file << "# vtk DataFile Version 3.0" << endl;
    file << format("TIME {:.3f}", t) << endl;
    file << "ASCII" << endl;
    file << "DATASET STRUCTURED_GRID" << endl;
    file << format("DIMENSIONS {} {} {}", nx, ny, nz) << endl;    
    file << format("POINTS {} float", nx*ny*nz) << endl;    

    vector<vector<double>> z(nx, vector<double>(ny, 0));

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                double x = dx * i;
                double y = dy * j;

                double dzc = dz[nz-2];

                if (k < nz-1) {
                    dzc = dz[k];
                }

                z[i][j] += h[i][j]*dzc / 2;
                
                file << format("{:.3f} {:.3f} {:.3f}", x, y, z[i][j]) << endl;

                z[i][j] += h[i][j]*dzc / 2;
            }
        }        
    }

    file << "FIELD FieldData 1" << endl;
    file << "Time 1 1 float" << endl;
    file << format("{:.3f}", t) << endl;
    file << format("POINT_DATA {}", nx*ny*nz) << endl;

    file << "SCALARS nu float" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                file << format("{:.3f}", nu[i][j][k]) << endl;
            }
        }
    }
}

void Solver::writeData(
    double t, int m, 
    int nx, int ny, int nz, 
    vector<vector<double>> &ua, 
    vector<vector<double>> &va, 
    vector<vector<double>> &z,
    vector<vector<vector<double>>> &uf, 
    vector<vector<vector<double>>> &vf, 
    path outDir
) {
    writeSurfaceData(t, m, nx, ny, ua, va, z, outDir);
    writeVolumeData(t, m, nx, ny, nz, ua, va, uf, vf, outDir);
    writeViscosity(t, m, nx, ny, nz + 1, outDir);
}

void Solver::writeStatistics(vector<tuple<int, double, long long, double, double, double>> &statistics, path outDir) {    
    auto filePath = outDir / path("convergence.csv");

    ofstream file(filePath);

    if (file.bad()) {
        throw runtime_error(
            format("Failed to open file at: {}", filePath.string())
        );
    }

    file << "n,  t, calc_time,  u,  v,  z" << endl;

    for (auto t: statistics) {
        file << format(
            "{},  {:.3f}, {},  {:.5f},  {:.5f},  {:.5f}", 
            get<0>(t),
            get<1>(t),
            get<2>(t),
            get<3>(t),
            get<4>(t),
            get<5>(t)
        ) << endl;
    }
}

void Solver::solve() {
    auto outDir = createDirectory();

    int nx = static_cast<int>(
        ceil(l/dx)
    ) + 1;

    int ny = static_cast<int>(
        ceil(w/dy)
    ) + 1;

    generateDZ();

    int nz = dz.size();

    generateH(nx, ny);    

    auto [hMin, hMax] = minMaxH(nx, ny);
    double dzMin = *min_element(dz.begin(), dz.end());

    const double dtMax = min(dx, dy)/sqrt(2*g*hMax)/1.5;
    
    vector<vector<double>> ua(nx, vector<double>(ny));
    vector<vector<double>> u1a(nx, vector<double>(ny));

    vector<vector<double>> va(nx, vector<double>(ny));
    vector<vector<double>> v1a(nx, vector<double>(ny));

    vector<vector<double>> z(nx, vector<double>(ny));

    vector<vector<vector<double>>> uf(nx, vector<vector<double>>(ny, vector<double>(nz)));
    vector<vector<vector<double>>> ud(nx, vector<vector<double>>(ny, vector<double>(nz)));

    vector<vector<vector<double>>> vf(nx, vector<vector<double>>(ny, vector<double>(nz)));
    vector<vector<vector<double>>> vd(nx, vector<vector<double>>(ny, vector<double>(nz)));

    setInitialCondition(nx, ny, nz, uf, vf, ua, va, z);    
    
    generateQX(nx, ny, nz, ua, va, z, uf, vf);
    generateQY(nx, ny, nz, ua, va, z, uf, vf);

    generateNU(nx, ny, nz + 1, ua, va, z, uf, vf);

    double nuMax = maxNU(nx, ny, nz + 1);

    double dt = min(dtMax, hMin*dzMin*hMin*dzMin/2/nuMax);
    double dtp = dt;

    vector<vector<double>> zp(nx, vector<double>(ny));

    vector<vector<vector<double>>> up(nx, vector<vector<double>>(ny, vector<double>(nz)));
    vector<vector<vector<double>>> vp(nx, vector<vector<double>>(ny, vector<double>(nz)));

    updateData(nx, ny, zp, z);
    updateData(nx, ny, nz, up, uf);
    updateData(nx, ny, nz, vp, vf);

    double t = 0;
    double tn = outputTimeStep;

    int n = 1;
    int m = 0;

    vector<tuple<int, double, long long, double, double, double>> statistics;

    writeHeights(nx, ny, outDir);
    writeData(t, m, nx, ny, nz, ua, va, z, uf, vf, outDir);

    m += 1;

    vector<double> al(nz);
    vector<double> ac(nz);
    vector<double> ar(nz); 

    auto start = high_resolution_clock::now();    

    long long calcTime = 0;    

    while (t <= endTime) {
        for (int j = 0; j < ny; j++) {
            for (int i = 1; i < nx-1; i++) {
                u1a[i][j] = ua[i][j] + f*va[i][j]*dt - g*dt*(z[i+1][j] - z[i-1][j])/2/dx - kb*dt*(uf[i][j][nz-1] + ua[i][j])/h[i][j] + qx[i][j]*dt/rho/h[i][j];
            }

            u1a[0][j] = ua[0][j] + f*va[0][j]*dt - g*dt*(z[1][j] - z[nx-1][j])/2/dx - kb*dt*(uf[0][j][nz-1] + ua[0][j])/h[0][j] + qx[0][j]*dt/rho/h[0][j];
            u1a[nx-1][j] = ua[nx-1][j] + f*va[nx-1][j]*dt - g*dt*(z[0][j] - z[nx-2][j])/2/dx - kb*dt*(uf[nx-1][j][nz-1] + ua[nx-1][j])/h[nx-1][j] + qx[nx-1][j]*dt/rho/h[nx-1][j];
        }

        for (int i = 0; i < nx; i++) {
            for (int j = 1; j < ny-1; j++) {
                v1a[i][j] = va[i][j] - f*ua[i][j]*dt - g*dt*(z[i][j+1] - z[i][j-1])/2/dy - kb*dt*(vf[i][j][nz-1] + va[i][j])/h[i][j] + qy[i][j]*dt/rho/h[i][j];
            }

            v1a[i][0] = va[i][0] - f*ua[i][0]*dt - g*dt*(z[i][1] - z[i][ny-1])/2/dy - kb*dt*(vf[i][0][nz-1] + va[i][0])/h[i][0] + qy[i][0]*dt/rho/h[i][0];
            v1a[i][ny-1] = va[i][ny-1] - f*ua[i][ny-1]*dt - g*dt*(z[i][0] - z[i][ny-2])/2/dy - kb*dt*(vf[i][ny-1][nz-1] + va[i][ny-1])/h[i][ny-1] + qy[i][ny-1]*dt/rho/h[i][ny-1];
        }   
        
        calculateResidualElements(dt, nx, ny, nz, ua, u1a, va, v1a, uf, vf, ud, vd);        

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                createTridiagonalMatrix(i, j, nz, dt, al, ac, ar);

                Tridiagonal::solve(al, ac, ar, ud[i][j], uf[i][j]);
                Tridiagonal::solve(al, ac, ar, vd[i][j], vf[i][j]); 
            }
        }

        t += dt;
        n += 1;

        updateQX(t, n, nx, ny, nz, ua, va, z, uf, vf);
        updateQY(t, n, nx, ny, nz, ua, va, z, uf, vf);
        updateNU(t, n, nx, ny, nz + 1, ua, va, z, uf, vf);        

        double nuMax = maxNU(nx, ny, nz + 1);

        double dt1 = min(dtMax, hMin*dzMin*hMin*dzMin/2/nuMax);

        if (dt1 < dtp) {
            dtp = dt1;
            dt = adjustTimeStep(t, dt1, outputTimeStep,  dtMax, false);
        } else {
            dt = adjustTimeStep(t, dt, outputTimeStep,  dtMax, true);
        }

        if (t >= tn) {
            auto end = high_resolution_clock::now();

            auto duration = duration_cast<seconds>(end - start).count();
            
            calcTime += duration;

            writeData(t, m, nx, ny, nz, ua, va, z, uf, vf, outDir);

            auto umd = maxAbsDifference(nx, ny, nz, up, uf);
            auto vmd = maxAbsDifference(nx, ny, nz, vp, vf);
            auto zmd = maxAbsDifference(nx, ny, zp, z);

            cout << format(
                "Write data in file t={:.3f}, convergence of u={:.5f}, v={:.5f}, z={:.7f} with dt={:.5}, calc time={}", 
                t, umd, vmd, zmd, dt, calcTime
            ) << endl;

            statistics.push_back({n, tn, calcTime, umd, vmd, zmd});

            updateData(nx, ny, zp, z);
            updateData(nx, ny, nz, up, uf);
            updateData(nx, ny, nz, vp, vf);

            m += 1;
            tn = t + outputTimeStep;

            start = high_resolution_clock::now();
        }

        updateData(nx, ny, ua, u1a);
        updateData(nx, ny, va, v1a);

        for (int i = 1; i < nx-1; i++) {
            for (int j = 1; j < ny-1; j++) {
                z[i][j] += -(h[i+1][j]*ua[i+1][j] - h[i-1][j]*ua[i-1][j])/2/dx - (h[i][j+1]*va[i][j+1] - h[i][j-1]*va[i][j-1])/2/dy;
            }
        }

        for (int j = 1; j < ny-1; j++) {
            z[0][j] += -(h[1][j]*ua[1][j] - h[nx-1][j]*ua[nx-1][j])/2/dx - (h[0][j+1]*va[0][j+1] - h[0][j-1]*va[0][j-1])/2/dy;
            z[nx-1][j] += -(h[0][j]*ua[0][j] - h[nx-2][j]*ua[nx-2][j])/2/dx - (h[nx-1][j+1]*va[nx-1][j+1] - h[nx-1][j-1]*va[nx-1][j-1])/2/dy;
        }

        for (int i = 1; i < nx-1; i++) {
            z[i][0] += -(h[i+1][0]*ua[i+1][0] - h[i-1][0]*ua[i-1][0])/2/dx - (h[i][1]*va[i][1] - h[i][ny-1]*va[i][ny-1])/2/dy;
            z[i][ny-1] += -(h[i+1][ny-1]*ua[i+1][ny-1] - h[i-1][ny-1]*ua[i-1][ny-1])/2/dx - (h[i][0]*va[i][0] - h[i][ny-2]*va[i][ny-2])/2/dy;
        }

        z[0][0] += -(h[1][0]*ua[1][0] - h[nx-1][0]*ua[nx-1][0])/2/dx - (h[0][1]*va[0][1] - h[0][ny-1]*va[0][ny-1])/2/dy;
        z[nx-1][0] += -(h[0][0]*ua[0][0] - h[nx-2][0]*ua[nx-2][0])/2/dx - (h[nx-1][1]*va[nx-1][1] - h[nx-1][ny-1]*va[nx-1][ny-1])/2/dy;

        z[0][ny-1] += -(h[1][ny-1]*ua[1][ny-1] - h[nx-1][ny-1]*ua[nx-1][ny-1])/2/dx - (h[0][0]*va[0][0] - h[0][ny-2]*va[0][ny-2])/2/dy;
        z[nx-1][ny-1] += -(h[0][ny-1]*ua[0][ny-1] - h[nx-2][ny-1]*ua[nx-2][ny-1])/2/dx - (h[nx-1][0]*va[nx-1][0] - h[nx-1][ny-2]*va[nx-1][ny-2])/2/dy;
    }

    writeStatistics(statistics, outDir);

    cout << "Total number of iterations: " << n;
}
