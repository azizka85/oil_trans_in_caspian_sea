#define _USE_MATH_DEFINES

#include <format>

#include <cmath>

#include <algorithm>

#include <fstream>
#include <iostream>

#include <stdexcept>

#include <slae/direct/tridiagonal.h>

#include "solver.h"

using namespace SLAE::Direct;

using namespace WindInducedCurrents::Davies85::VariableStepSize::VariableDepth;

Solver::Solver(
    double f, double b, 
    double h, double w, double l,
    double g, double rho, double nu, 
    double kb, double qx, double qy, 
    double dx, double dy, double dz, 
    double endTime, double outputTimeStep, string dir
) {
    setF(f);
    setB(b);
    
    setW(w);
    setL(l);

    setG(g);
    setRHO(rho);
    setNU(nu);

    setKB(kb);
    setQX(qx);
    setQY(qy);

    setDX(dx);
    setDY(dy);
    generateDZ(dz);

    generateH(h);

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

vector<vector<double>> Solver::getH() {
    return h;
}

void Solver::generateH(double hMax) {
    if (hMax <= 0) {
        throw runtime_error(
            format("H should be > 0, but it is {}", hMax)
        );
    }

    int nx = static_cast<int>(
        ceil(l/dx)
    ) + 1;

    int ny = static_cast<int>(
        ceil(w/dy)
    ) + 1;

    h = vector<vector<double>>(nx, vector<double>(ny));

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            h[i][j] = hMax*(1 - 0.8*cos(2*M_PI*(i/(nx-1.) + j/(ny-1.))))/2;
        }
    }
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

double Solver::getNU() {
    return nu;
}

void Solver::setNU(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("NU should be > 0, but it is {}", val)
        );
    }

    nu = val;
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

double Solver::getQX() {
    return qx;
}

void Solver::setQX(double val) {
    qx = val;
}

double Solver::getQY() {
    return qy;
}

void Solver::setQY(double val) {
    qy = val;
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

vector<double> Solver::getDZ() {
    return dz;
}

void Solver::generateDZ(double dzMax) {
    if (dzMax <= 0) {
        throw runtime_error(
            format("DZ should be > 0, but it is {}", dzMax)
        );
    }

    double z = 0;
    double dzf = calculateDZFactor(z);

    dz.push_back(dzMax*dzf);

    while (z < 1) {
        z += dzMax*dzf/2;

        if (z >= 1) {
            z -= dzMax*dzf/2;

            int n = dz.size();

            dz[n-1] /= 2;
            dz.push_back(2*(1 - z));

            break;
        }

        dzf = calculateDZFactor(z);
        dz.push_back(dzMax*dzf);

        z += dzMax*dzf/2;

        if (z >= 1) {
            z -= dzMax*dzf/2;

            int n = dz.size();

            dz[n-1] = 2*(1 - z);

            break;
        }
    }
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

double Solver::calculateDZFactor(double z) {
    return -2.96*z*z + 3.44*z + 0.02;
}

tuple<double, double> Solver::minMaxH() {    
    int nx = h.size();
    
    if (nx > 0) {
        int ny = h[0].size();

        if (nx > 0 && ny > 0) {
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
    }

    return {0, 0};
}

path Solver::createDirectory(double dzm, double hMin, double hMax)
{
    auto dirPath = path(
        format("{}/f={}, H={:.3f}-{:.3f}/nu={}, g={}, rho={}/kb={}, qx={}, qy={}/dx={}, dy={}, dz={}", dir, f, hMin, hMax, nu, g, rho, kb, qx, qy, dx, dy, dzm)
    );

    create_directories(dirPath);    

    auto surfacePath = dirPath / path("surface");

    create_directories(surfacePath);

    auto volumePath = dirPath / path("volume");

    create_directories(volumePath);

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

void Solver::calculate_residual_elements(
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
                ud[i][j][k] = uf[i][j][k] + f*dt*vf[i][j][k] + nu*dt*(uf[i][j][k+1] - 2*uf[i][j][k] + uf[i][j][k-1])/2/h[i][j]/h[i][j]/dz[k]/dz[k] 
                                + kb*dt*(uf[i][j][nz-1] + ua[i][j])/h[i][j] - qx*dt/rho/h[i][j];

                vd[i][j][k] = vf[i][j][k] - f*dt*uf[i][j][k] + nu*dt*(vf[i][j][k+1] - 2*vf[i][j][k] + vf[i][j][k-1])/2/h[i][j]/h[i][j]/dz[k]/dz[k] 
                                + kb*dt*(vf[i][j][nz-1] + va[i][j])/h[i][j] - qy*dt/rho/h[i][j];
            }

            ud[i][j][0] = uf[i][j][0] + f*dt*vf[i][j][0] + nu*dt*(uf[i][j][1] - uf[i][j][0])/h[i][j]/h[i][j]/dz[0]/dz[0] + kb*dt*(uf[i][j][nz-1] + ua[i][j])/h[i][j]
                            - qx*dt/rho/h[i][j] + 2*qx*dt/rho/h[i][j]/dz[0];
            ud[i][j][nz-1] = uf[i][j][nz-1] + f*dt*vf[i][j][nz-1] + nu*dt*(-(1 + kb*h[i][j]*dz[nz-1]/nu)*uf[i][j][nz-1] + uf[i][j][nz-2])/h[i][j]/h[i][j]/dz[nz-1]/dz[nz-1]
                            + kb*dt*(uf[i][j][nz-1] + ua[i][j])/h[i][j] - qx*dt/rho/h[i][j] - kb*dt*(u1a[i][j] + ua[i][j])/h[i][j]/dz[nz-1];

            vd[i][j][0] = vf[i][j][0] - f*dt*uf[i][j][0] + nu*dt*(vf[i][j][1] - vf[i][j][0])/h[i][j]/h[i][j]/dz[0]/dz[0] + kb*dt*(vf[i][j][nz-1] + va[i][j])/h[i][j]
                            - qy*dt/rho/h[i][j] + 2*qy*dt/rho/h[i][j]/dz[0];
            vd[i][j][nz-1] = vf[i][j][nz-1] - f*dt*uf[i][j][nz-1] + nu*dt*(-(1 + kb*h[i][j]*dz[nz-1]/nu)*vf[i][j][nz-1] + vf[i][j][nz-2])/h[i][j]/h[i][j]/dz[nz-1]/dz[nz-1]
                            + kb*dt*(vf[i][j][nz-1] + va[i][j])/h[i][j] - qy*dt/rho/h[i][j] - kb*dt*(v1a[i][j] + va[i][j])/h[i][j]/dz[nz-1];
        }
    }
}

double Solver::adjustTimeStep(double t, double dt, double dtMax) {
    if (t >= dtMax) {
        return dtMax;
    }

    if (t < dtMax && t + dt >= dtMax) {
        return dtMax - t;
    }

    return b * dt;
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

    file << "SCALARS z float" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << format("{:.3f}", z[i][j]) << endl;
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

    double z = 0;

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                if (k > 0) {
                    z += h[i][j]*dz[k]/2;
                }

                double x = dx * i;
                double y = dy * j;
                
                file << format("{:.3f} {:.3f} {:.3f}", x, y, z) << endl;

                z += h[i][j]*dz[k]/2;
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
}

void Solver::writeStatistics(vector<tuple<int, double, double, double, double>> &statistics, path outDir) {    
    auto filePath = outDir / path("convergence.csv");

    ofstream file(filePath);

    if (file.bad()) {
        throw runtime_error(
            format("Failed to open file at: {}", filePath.string())
        );
    }

    file << "n,  t,  u,  v,  z" << endl;

    for (auto t: statistics) {
        file << format(
            "{},  {:.3f},  {:.5f},  {:.5f},  {:.5f}", 
            get<0>(t),
            get<1>(t),
            get<2>(t),
            get<3>(t),
            get<4>(t)
        ) << endl;
    }
}

void Solver::solve() {
    auto [hMin, hMax] = minMaxH();

    const double dtMax = min(dx, dy)/sqrt(2*g*hMax);
    double dzMin = *min_element(dz.begin(), dz.end());
    double dt = min(dtMax, hMin*dzMin*hMin*dzMin/2/nu);

    auto outDir = createDirectory(dzMin, hMin, hMax);

    int nx = static_cast<int>(
        ceil(l/dx)
    ) + 1;

    int ny = static_cast<int>(
        ceil(w/dy)
    ) + 1;

    int nz = dz.size();

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

    vector<tuple<int, double, double, double, double>> statistics;

    writeHeights(nx, ny, outDir);
    writeData(t, m, nx, ny, nz, ua, va, z, uf, vf, outDir);

    m += 1;

    vector<double> al(nz);
    vector<double> ac(nz);
    vector<double> ar(nz);    

    while (t <= endTime) {
        for (int j = 0; j < ny; j++) {
            for (int i = 1; i < nx-1; i++) {
                u1a[i][j] = ua[i][j] + f*va[i][j]*dt - g*dt*(z[i+1][j] - z[i-1][j])/2/dx - kb*dt*(uf[i][j][nz-1] + ua[i][j])/h[i][j] + qx*dt/rho/h[i][j];                
            }

            u1a[0][j] = ua[0][j] + f*va[0][j]*dt - g*dt*(z[1][j] - z[nx-1][j])/2/dx - kb*dt*(uf[0][j][nz-1] + ua[0][j])/h[0][j] + qx*dt/rho/h[0][j];
            u1a[nx-1][j] = ua[nx-1][j] + f*va[nx-1][j]*dt - g*dt*(z[0][j] - z[nx-2][j])/2/dx - kb*dt*(uf[nx-1][j][nz-1] + ua[nx-1][j])/h[nx-1][j] + qx*dt/rho/h[nx-1][j];
        }

        for (int i = 0; i < nx; i++) {
            for (int j = 1; j < ny-1; j++) {
                v1a[i][j] = va[i][j] - f*ua[i][j]*dt - g*dt*(z[i][j+1] - z[i][j-1])/2/dy - kb*dt*(vf[i][j][nz-1] + va[i][j])/h[i][j] + qy*dt/rho/h[i][j];
            }

            v1a[i][0] = va[i][0] - f*ua[i][0]*dt - g*dt*(z[i][1] - z[i][ny-1])/2/dy - kb*dt*(vf[i][0][nz-1] + va[i][0])/h[i][0] + qy*dt/rho/h[i][0];
            v1a[i][ny-1] = va[i][ny-1] - f*ua[i][ny-1]*dt - g*dt*(z[i][0] - z[i][ny-2])/2/dy - kb*dt*(vf[i][ny-1][nz-1] + va[i][ny-1])/h[i][ny-1] + qy*dt/rho/h[i][ny-1];
        }                

        calculate_residual_elements(dt, nx, ny, nz, ua, u1a, va, v1a, uf, vf, ud, vd);        

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                ac[0] = 1 + dt*nu/h[i][j]/h[i][j]/dz[0]/dz[0];
                ar[0] = -dt*nu/h[i][j]/h[i][j]/dz[0]/dz[0];

                for (int k = 1; k < nz-1; k++) {
                    al[k] = -dt*nu/2/h[i][j]/h[i][j]/dz[k]/dz[k];
                    ac[k] = 1 + dt*nu/h[i][j]/h[i][j]/dz[k]/dz[k];
                    ar[k] = -dt*nu/2/h[i][j]/h[i][j]/dz[k]/dz[k];
                }

                al[nz-1] = -dt*nu/h[i][j]/h[i][j]/dz[nz-1]/dz[nz-1];
                ac[nz-1] = 1 + dt*nu/h[i][j]/h[i][j]/dz[nz-1]/dz[nz-1] + dt*kb/h[i][j]/dz[nz-1];

                Tridiagonal::solve(al, ac, ar, ud[i][j], uf[i][j]);
                Tridiagonal::solve(al, ac, ar, vd[i][j], vf[i][j]);                
            }
        }

        t += dt;

        dt = adjustTimeStep(t, dt, dtMax);

        if (t >= tn) {
            writeData(t, m, nx, ny, nz, ua, va, z, uf, vf, outDir);

            auto umd = maxAbsDifference(nx, ny, nz, up, uf);
            auto vmd = maxAbsDifference(nx, ny, nz, vp, vf);
            auto zmd = maxAbsDifference(nx, ny, zp, z);

            cout << format("Write data in file t={:.3f}, convergence of u={:.5f}, v={:.5f}, z={:.5f}", t, umd, vmd, zmd) << endl;

            statistics.push_back({n, tn, umd, vmd, zmd});

            updateData(nx, ny, zp, z);
            updateData(nx, ny, nz, up, uf);
            updateData(nx, ny, nz, vp, vf);

            m += 1;
            tn += outputTimeStep;
        }

        updateData(nx, ny, ua, u1a);
        updateData(nx, ny, va, v1a);

        for (int i = 1; i < nx-1; i++) {
            for (int j = 1; j < ny-1; j++) {
                z[i][j] += -h[i][j]*(ua[i+1][j] - ua[i-1][j])/2/dx - h[i][j]*(va[i][j+1] - va[i][j-1])/2/dy;
            }
        }

        for (int j = 1; j < ny-1; j++) {
            z[0][j] += -h[0][j]*(ua[1][j] - ua[nx-1][j])/2/dx - h[0][j]*(va[0][j+1] - va[0][j-1])/2/dy;
            z[nx-1][j] += -h[nx-1][j]*(ua[0][j] - ua[nx-2][j])/2/dx - h[nx-1][j]*(va[nx-1][j+1] - va[nx-1][j-1])/2/dy;
        }

        for (int i = 1; i < nx-1; i++) {
            z[i][0] += -h[i][0]*(ua[i+1][0] - ua[i-1][0])/2/dx - h[i][0]*(va[i][1] - va[i][ny-1])/2/dy;
            z[i][ny-1] += -h[i][ny-1]*(ua[i+1][ny-1] - ua[i-1][ny-1])/2/dx - h[i][ny-1]*(va[i][0] - va[i][ny-2])/2/dy;
        }

        z[0][0] += -h[0][0]*(ua[1][0] - ua[nx-1][0])/2/dx - h[0][0]*(va[0][1] - va[0][ny-1])/2/dy;
        z[nx-1][0] += -h[nx-1][0]*(ua[0][0] - ua[nx-2][0])/2/dx - h[nx-1][0]*(va[nx-1][1] - va[nx-1][ny-1])/2/dy;

        z[0][ny-1] += -h[0][ny-1]*(ua[1][ny-1] - ua[nx-1][ny-1])/2/dx - h[0][ny-1]*(va[0][0] - va[0][ny-2])/2/dy;
        z[nx-1][ny-1] += -h[nx-1][ny-1]*(ua[0][ny-1] - ua[nx-2][ny-1])/2/dx - h[nx-1][ny-1]*(va[nx-1][0] - va[nx-1][ny-2])/2/dy;

        n += 1;
    }

    writeStatistics(statistics, outDir);

    cout << "Total number of iterations: " << n;
}
