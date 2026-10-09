#define _USE_MATH_DEFINES

#include <format>

#include <cmath>

#include <fstream>
#include <sstream>
#include <iostream>

#include <chrono>

#include <stdexcept>

#include "solver.h"

using namespace std::chrono;

using namespace WindInducedCurrents::Davies85::Parallel::VariableParameters;

Solver::Solver(
    double f, double b,
    double hm, GenerateH hg,
    double w, double l,
    double g, double rho, double kb,
    double num, double nub, double hp, GenerateNU nug,
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
	setNUB(nub);
    setHP(hp);
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

vector<double> Solver::getH() {
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

double Solver::getNUB() {
    return nub;
}

void Solver::setNUB(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("NUB should be > 0, but it is {}", val)
        );
    }

    nub = val;
}

double Solver::getHP() {
    return hp;
}

void Solver::setHP(double val) {
    if (val < 0) {
        throw runtime_error(
            format("HP should be >= 0, but it is {}", val)
        );
    }

    hp = val;
}

GenerateNU Solver::getNUG() {
    return nug;
}

void Solver::setNUG(GenerateNU val) {
    nug = val;
}

vector<double> Solver::getNU() {
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

vector<double> Solver::getQX() {
    return qx;
}

vector<double> Solver::getQY() {
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
    h = vector<double>(nx*ny);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int id = j + i*ny;

            h[id] = hm;
        }
    }
}

void Solver::generateCosineH(int nx, int ny) {
    h = vector<double>(nx*ny);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int id = j + i*ny;

            h[id] = hm*(1 - 0.8*cos(2*M_PI*(i/(nx-1.) + j/(ny-1.))))/2;
        }
    }
}

tuple<double, double> Solver::minMaxH(int nx, int ny) {
    double minH = h[0];
    double maxH = h[0];

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int id = j + i*ny;

            minH = min(minH, h[id]);
            maxH = max(maxH, h[id]);
        }
    }

    return {minH, maxH};
}

void Solver::generateNU(
    int nx, int ny, int nz, 
    vector<double> &ua, 
    vector<double> &va, 
    vector<double> &z, 
    vector<double> &uf, 
    vector<double> &vf
) {
    switch(nug) {
        case GenerateNU::UniformNU:
            generateUniformNU(nx, ny, nz);
            break;
        case GenerateNU::FromWindSpeedNU:
            generateNUFromWindSpeed(nx, ny, nz);
            break;
        case GenerateNU::LinearNU:
            generateLinearNU(nx, ny, nz);
			break;
    }
}

void Solver::generateUniformNU(int nx, int ny, int nz) {
    nu = vector<double>(nx*ny*nz);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i*ny;

            for (int k = 0; k < nz; k++) {
                int id = k + p*nz;

                nu[id] = num;
            }
        }
    }
}

void Solver::generateNUFromWindSpeed(int nx, int ny, int nz) {
    generateUniformNU(nx, ny, nz);
}

void Solver::generateLinearNU(int nx, int ny, int nz) {
    nu = vector<double>(nx * ny * nz);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i * ny;

            double z = 0;

            for (int k = 0; k < nz; k++) {
                int id = k + p * nz;

                if (z < hp) {
                    nu[id] = num - (num - nub) * z / hp;
                }
                else {
                    nu[id] = nub;
                }

                if (k < nz - 1) {
                    z += h[p] * dz[k];
                }
            }
        }
    }
}

cl::Kernel Solver::createUpdateNUKernel(
    int ny, int nz, 
    cl::Buffer& bufferNU,
    cl::Program& program
) {
    switch(nug) {
        case GenerateNU::UniformNU:
            return createUpdateUniformNUKernel(ny, nz, bufferNU, program);
        case GenerateNU::FromWindSpeedNU:
            return createUpdateNUFromWindSpeed(ny, nz, bufferNU, program);
		case GenerateNU::LinearNU:
			return createUpdateLinearNUKernel(program);
    }
}

cl::Kernel Solver::createUpdateUniformNUKernel(int ny, int nz, cl::Buffer& bufferNU, cl::Program& program) {
    cl::Kernel kernel(program, "wind_induced_currents_davies85_variable_parameters_update_uniform_nu");

    kernel.setArg(0, num);
    kernel.setArg(1, ny);
    kernel.setArg(2, nz);
    kernel.setArg(3, bufferNU);

    return kernel;
}

cl::Kernel Solver::createUpdateNUFromWindSpeed(int ny, int nz, cl::Buffer& bufferNU, cl::Program& program) {
    return createUpdateUniformNUKernel(ny, nz, bufferNU, program);
}

cl::Kernel Solver::createUpdateLinearNUKernel(cl::Program& program) {
    cl::Kernel kernel(program, "wind_induced_currents_davies85_variable_parameters_update_linear_nu");

    return kernel;
}

double Solver::maxNU(int nx, int ny, int nz)
{
    double nuMax = nu[0];

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i*ny;

            for (int k = 0; k < nz; k++) {
                int id = k + p*nz;

                nuMax = max(nuMax, nu[id]);
            }
        }
    }

    return nuMax;
}

void Solver::generateQX(
    int nx, int ny, int nz, 
    vector<double> &ua, 
    vector<double> &va, 
    vector<double> &z, 
    vector<double> &uf, 
    vector<double> &vf
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

void Solver::generateUniformQX(int nx, int ny) {
    qx = vector<double>(nx*ny);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int id = j + i*ny;

            qx[id] = qxm;
        }
    }
}

void Solver::generateQXFromWindSpeed(int nx, int ny) {
    generateUniformQX(nx, ny);
}

void Solver::generateQY(
    int nx, int ny, int nz, 
    vector<double> &ua, 
    vector<double> &va, 
    vector<double> &z, 
    vector<double> &uf, 
    vector<double> &vf
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

void Solver::generateUniformQY(int nx, int ny) {
    qy = vector<double>(nx*ny);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int id = j + i*ny;

            qy[id] = qym;
        }
    }
}

void Solver::generateQYFromWindSpeed(int nx, int ny) {
    generateUniformQY(nx, ny);
}

cl::Kernel Solver::createUpdateQKernel(double qm, int ny, cl::Buffer &bufferQ, cl::Program& program) {
    switch(qg) {
        case GenerateQ::UniformQ:
            return createUpdateUniformQKernel(qm, ny, bufferQ, program);
        case GenerateQ::FromWindSpeedQ:
            return createUpdateQFromWindSpeedKernel(qm, ny, bufferQ, program);
    }
}

cl::Kernel Solver::createUpdateUniformQKernel(double qm, int ny, cl::Buffer &bufferQ, cl::Program &program) {
    cl::Kernel kernel(program, "wind_induced_currents_davies85_variable_parameters_update_uniform_q");

    kernel.setArg(0, qm);
    kernel.setArg(1, ny);
    kernel.setArg(2, bufferQ);

    return kernel;
}

cl::Kernel Solver::createUpdateQFromWindSpeedKernel(double qm, int ny, cl::Buffer &bufferQ, cl::Program &program) {
    return createUpdateUniformQKernel(qm, ny, bufferQ, program);
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

    dz.push_back(dzm * dzf);

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
        dz.push_back(dzm * dzf);
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

    if (nug == GenerateNU::LinearNU) {
		nuStr = format("LNU, nus={}, num={}, hp={}", num, nub, hp);
    }

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

string Solver::loadKernelSource(path filePath) {
    ifstream file(filePath);    

    if (file.bad()) throw runtime_error("Cannot open file: " + filePath.string());

    return std::string(istreambuf_iterator<char>(file),
                        istreambuf_iterator<char>());
}

void Solver::loadKernelSources(cl::Program::Sources &sources) {
    string solverSrc = loadKernelSource("opencl/wind_induced_currents/davies85/variable_parameters/solver.cl");

    sources.push_back(solverSrc);

    string maxSrc = loadKernelSource("opencl/utils/max.cl");

    sources.push_back(maxSrc);

    switch(qg) {
        case GenerateQ::UniformQ:
            sources.push_back(
                loadKernelSource("opencl/wind_induced_currents/davies85/variable_parameters/q/update_uniform_q.cl")
            );
            break;
        case GenerateQ::FromWindSpeedQ:
            sources.push_back(
                loadKernelSource("opencl/wind_induced_currents/davies85/variable_parameters/q/update_uniform_q.cl")
            );
            break;
    }

    switch(nug) {
        case GenerateNU::UniformNU:
            sources.push_back(
                loadKernelSource("opencl/wind_induced_currents/davies85/variable_parameters/nu/update_uniform_nu.cl")
            );
            break;
        case GenerateNU::FromWindSpeedNU:
            sources.push_back(
                loadKernelSource("opencl/wind_induced_currents/davies85/variable_parameters/nu/update_uniform_nu.cl")
            );
            break;
    }
}

void Solver::setInitialCondition(
    int nx, int ny, int nz, 
    vector<double>& uf, 
    vector<double>& vf, 
    vector<double>& ua,
    vector<double>& va,
    vector<double>& z
) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i*ny;

            for (int k = 0; k < nz; k++) {                
                int id = k + p*nz;

                uf[id] = 0;
                vf[id] = 0;
            }

            ua[p] = 0;
            va[p] = 0;
            z[p] = 0;
        }
    }
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

double Solver::maxAbsDifference(int nx, int ny, vector<double> &u, vector<double> &u1) {
    double maxDiff = 0.;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i*ny;

            maxDiff = max(
                maxDiff,
                abs(u[p] - u1[p])
            );
        }
    }

    return maxDiff;
}

double Solver::maxAbsDifference(int nx, int ny, int nz, vector<double> &u, vector<double> &u1) {
    double maxDiff = 0.;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i*ny;

            for (int k = 0; k < nz; k++) {
                int id = k + p*nz;

                maxDiff = max(
                    maxDiff,
                    abs(u[id] - u1[id])
                );
            }
        }
    }

    return maxDiff;
}

void Solver::updateData(int nx, int ny, vector<double> &u, vector<double> &u1) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i*ny;

            u[p] = u1[p];
        }
    }
}

void Solver::updateData(int nx, int ny, int nz, vector<double> &u, vector<double> &u1) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i*ny;

            for (int k = 0; k < nz; k++) {
                int id = k + p*nz;

                u[id] = u1[id];
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
    file << format("POINTS {} double", nx*ny) << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double x = dx * i;
            double y = dy * j;

            file << format("{:.3f} {:.3f} 0.0", x, y) << endl;
        }
    }

    file << "FIELD FieldData 1" << endl;
    file << "Time 1 1 double" << endl;
    file << format("{:.3f}", 0.) << endl;
    file << format("POINT_DATA {}", nx*ny) << endl;

    file << "SCALARS h double" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int p = j + i*ny;

            file << format("{:.3f}", h[p]) << endl;
        }
    }
}

void Solver::writeSurfaceData(
    double t, int m, 
    int nx, int ny, 
    vector<double> &ua, 
    vector<double> &va, 
    vector<double> &z,
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
    file << format("POINTS {} double", nx*ny) << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double x = dx * i;
            double y = dy * j;

            file << format("{:.3f} {:.3f} 0.0", x, y) << endl;
        }
    }

    file << "FIELD FieldData 1" << endl;
    file << "Time 1 1 double" << endl;
    file << format("{:.3f}", t) << endl;
    file << format("POINT_DATA {}", nx*ny) << endl;

    file << "VECTORS VA double" << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int p = j + i*ny;

            file << format("{:.3f} {:.3f} 0", ua[p], va[p]) << endl;
        }
    }

    file << "VECTORS Q double" << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int p = j + i*ny;

            file << format("{:.3f} {:.3f} 0", qx[p], qy[p]) << endl;
        }
    }

    file << "SCALARS z double" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int p = j + i*ny;

            file << format("{:.7f}", z[p]) << endl;
        }
    }
}

void Solver::writeVolumeData(
    double t, int m, 
    int nx, int ny, int nz, 
    vector<double> &ua, 
    vector<double> &va, 
    vector<double> &uf, 
    vector<double> &vf, 
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
    file << format("POINTS {} double", nx*ny*nz) << endl;    

    vector<vector<double>> z(nx, vector<double>(ny, 0));

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                double x = dx * i;
                double y = dy * j;

                int p = j + i*ny;

                z[i][j] += h[p]*dz[k] / 2;
                
                file << format("{:.3f} {:.3f} {:.3f}", x, y, z[i][j]) << endl;

                z[i][j] += h[p]*dz[k] / 2;
            }
        }        
    }

    file << "FIELD FieldData 1" << endl;
    file << "Time 1 1 double" << endl;
    file << format("{:.3f}", t) << endl;
    file << format("POINT_DATA {}", nx*ny*nz) << endl;

    file << "VECTORS V double" << endl;

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int p = j + i*ny;
                int id = k + p*nz;

                file << format("{:.3f} {:.3f} 0", uf[id] + ua[p], vf[id] + va[p]) << endl;
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
    file << format("POINTS {} double", nx*ny*nz) << endl;    

    vector<vector<double>> z(nx, vector<double>(ny, 0));

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                double x = dx * i;
                double y = dy * j;

                file << format("{:.3f} {:.3f} {:.3f}", x, y, z[i][j]) << endl;

                int p = j + i * ny;

                if (k < nz-1) {
                    z[i][j] += h[p] * dz[k];
                }                     
            }
        }        
    }

    file << "FIELD FieldData 1" << endl;
    file << "Time 1 1 double" << endl;
    file << format("{:.3f}", t) << endl;
    file << format("POINT_DATA {}", nx*ny*nz) << endl;

    file << "SCALARS nu double" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int p = j + i*ny;
                int id = k + p*nz;

                file << format("{:.3f}", nu[id]) << endl;
            }
        }
    }
}

void Solver::writeData(
    double t, int m, 
    int nx, int ny, int nz, 
    vector<double> &ua, 
    vector<double> &va, 
    vector<double> &z,
    vector<double> &uf, 
    vector<double> &vf, 
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
    
    vector<double> ua(nx*ny);
    vector<double> u1a(nx*ny);

    vector<double> va(nx*ny);
    vector<double> v1a(nx*ny);

    vector<double> z(nx*ny);

    vector<double> uf(nx*ny*nz);
    vector<double> ud(nx*ny*nz);

    vector<double> vf(nx*ny*nz);
    vector<double> vd(nx*ny*nz);

    setInitialCondition(nx, ny, nz, uf, vf, ua, va, z);    
    
    generateQX(nx, ny, nz, ua, va, z, uf, vf);
    generateQY(nx, ny, nz, ua, va, z, uf, vf);

    generateNU(nx, ny, nz + 1, ua, va, z, uf, vf);

    double nuMax = maxNU(nx, ny, nz + 1);

    double dt = min(dtMax, hMin*dzMin*hMin*dzMin/2/nuMax);
    double dtp = dt;

    vector<double> zp(nx*ny);

    vector<double> up(nx*ny*nz);
    vector<double> vp(nx*ny*nz);

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

    vector<cl::Platform> platforms;

    cl::Platform::get(&platforms);

    cl::Platform platform = platforms.front();

    vector<cl::Device> devices;

    platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);

    cl::Device device = devices.front();

    cl::Context context(device);

    cl::Program::Sources sources;

    loadKernelSources(sources);

    cl::Program program(context, sources);

    auto err = program.build(device);

    cl::Buffer bufferV(context, CL_MEM_READ_WRITE, sizeof(double));
    cl::Buffer bufferR(context, CL_MEM_READ_WRITE, sizeof(double) * nx);
    cl::Buffer bufferP(context, CL_MEM_READ_WRITE, sizeof(double) * nx * ny);

    cl::Buffer bufferDZ(context, CL_MEM_READ_WRITE, sizeof(double) * nz);

    cl::Buffer bufferH(context, CL_MEM_READ_WRITE, sizeof(double) * nx * ny);
    cl::Buffer bufferQX(context, CL_MEM_READ_WRITE, sizeof(double) * nx * ny);
    cl::Buffer bufferQY(context, CL_MEM_READ_WRITE, sizeof(double) * nx * ny);
    cl::Buffer bufferZ(context, CL_MEM_READ_WRITE, sizeof(double) * nx * ny);

    cl::Buffer bufferUA(context, CL_MEM_READ_WRITE, sizeof(double) * nx * ny);
    cl::Buffer bufferVA(context, CL_MEM_READ_WRITE, sizeof(double) * nx * ny);
    cl::Buffer bufferU1A(context, CL_MEM_READ_WRITE, sizeof(double) * nx * ny);
    cl::Buffer bufferV1A(context, CL_MEM_READ_WRITE, sizeof(double) * nx * ny);

    cl::Buffer bufferNU(context, CL_MEM_READ_WRITE, sizeof(double) * nx * ny * (nz + 1));

    cl::Buffer bufferUF(context, CL_MEM_READ_WRITE, sizeof(double) * nx * ny * nz);
    cl::Buffer bufferVF(context, CL_MEM_READ_WRITE, sizeof(double) * nx * ny * nz);

    cl::Buffer bufferUD(context, CL_MEM_READ_WRITE, sizeof(double) * nx * ny * nz);
    cl::Buffer bufferVD(context, CL_MEM_READ_WRITE, sizeof(double) * nx * ny * nz);

    cl::Buffer bufferAL(context, CL_MEM_READ_WRITE, sizeof(double) * nx * ny * nz);
    cl::Buffer bufferAC(context, CL_MEM_READ_WRITE, sizeof(double) * nx * ny * nz);
    cl::Buffer bufferAR(context, CL_MEM_READ_WRITE, sizeof(double) * nx * ny * nz);

    cl::CommandQueue queue(context, device);

    err = queue.enqueueWriteBuffer(bufferDZ, CL_TRUE, 0, sizeof(double) * nz, dz.data());

    err = queue.enqueueWriteBuffer(bufferH, CL_TRUE, 0, sizeof(double) * nx * ny, h.data());
    err = queue.enqueueWriteBuffer(bufferQX, CL_TRUE, 0, sizeof(double) * nx * ny, qx.data());
    err = queue.enqueueWriteBuffer(bufferQY, CL_TRUE, 0, sizeof(double) * nx * ny, qy.data());
    err = queue.enqueueWriteBuffer(bufferZ, CL_TRUE, 0, sizeof(double) * nx * ny, z.data());

    err = queue.enqueueWriteBuffer(bufferUA, CL_TRUE, 0, sizeof(double) * nx * ny, ua.data());
    err = queue.enqueueWriteBuffer(bufferVA, CL_TRUE, 0, sizeof(double) * nx * ny, va.data());

    err = queue.enqueueWriteBuffer(bufferNU, CL_TRUE, 0, sizeof(double) * nx * ny * (nz + 1), nu.data());

    err = queue.enqueueWriteBuffer(bufferUF, CL_TRUE, 0, sizeof(double) * nx * ny * nz, uf.data());
    err = queue.enqueueWriteBuffer(bufferVF, CL_TRUE, 0, sizeof(double) * nx * ny * nz, vf.data());

    cl::Kernel updateUAKernel(program, "wind_induced_currents_davies85_variable_parameters_calc_ua");

    updateUAKernel.setArg(0, f);
    updateUAKernel.setArg(1, kb);
    updateUAKernel.setArg(2, g);    
    updateUAKernel.setArg(3, rho);    
    updateUAKernel.setArg(4, dx);    
    updateUAKernel.setArg(5, dt);        
    updateUAKernel.setArg(6, nx);        
    updateUAKernel.setArg(7, ny);        
    updateUAKernel.setArg(8, nz);        
    updateUAKernel.setArg(9, bufferH);        
    updateUAKernel.setArg(10, bufferQX);        
    updateUAKernel.setArg(11, bufferUA);        
    updateUAKernel.setArg(12, bufferVA);        
    updateUAKernel.setArg(13, bufferZ);        
    updateUAKernel.setArg(14, bufferUF);        
    updateUAKernel.setArg(15, bufferU1A);        

    cl::Kernel updateVAKernel(program, "wind_induced_currents_davies85_variable_parameters_calc_va");

    updateVAKernel.setArg(0, f);
    updateVAKernel.setArg(1, kb);
    updateVAKernel.setArg(2, g);
    updateVAKernel.setArg(3, rho);
    updateVAKernel.setArg(4, dy);
    updateVAKernel.setArg(5, dt);
    updateVAKernel.setArg(6, nx);
    updateVAKernel.setArg(7, ny);
    updateVAKernel.setArg(8, nz);
    updateVAKernel.setArg(9, bufferH);
    updateVAKernel.setArg(10, bufferQY);
    updateVAKernel.setArg(11, bufferUA);
    updateVAKernel.setArg(12, bufferVA);
    updateVAKernel.setArg(13, bufferZ);
    updateVAKernel.setArg(14, bufferVF);
    updateVAKernel.setArg(15, bufferV1A);

    cl::Kernel calcRHSKernel(program, "wind_induced_currents_davies85_variable_parameters_calc_rhs");

    calcRHSKernel.setArg(0, f);
    calcRHSKernel.setArg(1, kb);
    calcRHSKernel.setArg(2, rho);
    calcRHSKernel.setArg(3, dx);
    calcRHSKernel.setArg(4, dy);
    calcRHSKernel.setArg(5, dt);
    calcRHSKernel.setArg(6, ny);
    calcRHSKernel.setArg(7, nz);
    calcRHSKernel.setArg(8, bufferNU);
    calcRHSKernel.setArg(9, bufferDZ);
    calcRHSKernel.setArg(10, bufferH);
    calcRHSKernel.setArg(11, bufferQX);
    calcRHSKernel.setArg(12, bufferQY);
    calcRHSKernel.setArg(13, bufferUA);
    calcRHSKernel.setArg(14, bufferU1A);
    calcRHSKernel.setArg(15, bufferVA);
    calcRHSKernel.setArg(16, bufferV1A);
    calcRHSKernel.setArg(17, bufferUF);
    calcRHSKernel.setArg(18, bufferVF);
    calcRHSKernel.setArg(19, bufferUD);
    calcRHSKernel.setArg(20, bufferVD);

    cl::Kernel createTridiagonalMatrixKernel(program, "wind_induced_currents_davies85_variable_parameters_create_tridiagonal_matrix");

    createTridiagonalMatrixKernel.setArg(0, kb);
    createTridiagonalMatrixKernel.setArg(1, dt);
    createTridiagonalMatrixKernel.setArg(2, ny);
    createTridiagonalMatrixKernel.setArg(3, nz);
    createTridiagonalMatrixKernel.setArg(4, bufferNU);
    createTridiagonalMatrixKernel.setArg(5, bufferDZ);
    createTridiagonalMatrixKernel.setArg(6, bufferH);
    createTridiagonalMatrixKernel.setArg(7, bufferAL);
    createTridiagonalMatrixKernel.setArg(8, bufferAC);
    createTridiagonalMatrixKernel.setArg(9, bufferAR);

    cl::Kernel updateUVFKernel(program, "wind_induced_currents_davies85_variable_parameters_calc_uvf");

    updateUVFKernel.setArg(0, ny);
    updateUVFKernel.setArg(1, nz);
    updateUVFKernel.setArg(2, bufferAL);
    updateUVFKernel.setArg(3, bufferAC);
    updateUVFKernel.setArg(4, bufferAR);
    updateUVFKernel.setArg(5, bufferUF);
    updateUVFKernel.setArg(6, bufferVF);
    updateUVFKernel.setArg(7, bufferUD);
    updateUVFKernel.setArg(8, bufferVD);

    cl::Kernel updateZKernel(program, "wind_induced_currents_davies85_variable_parameters_calc_z");

    updateZKernel.setArg(0, dx);
    updateZKernel.setArg(1, dy);
    updateZKernel.setArg(2, nx);
    updateZKernel.setArg(3, ny);
    updateZKernel.setArg(4, nz);
    updateZKernel.setArg(5, bufferH);
    updateZKernel.setArg(6, bufferUA);
    updateZKernel.setArg(7, bufferVA);
    updateZKernel.setArg(8, bufferZ);

    cl::Kernel calcMaxPlaneKernel(program, "utils_max_plane");

    calcMaxPlaneKernel.setArg(0, ny);
    calcMaxPlaneKernel.setArg(1, nz + 1);
    calcMaxPlaneKernel.setArg(2, bufferNU);
    calcMaxPlaneKernel.setArg(3, bufferP);

    cl::Kernel calcMaxRowKernel(program, "utils_max_row");

    calcMaxRowKernel.setArg(0, ny);
    calcMaxRowKernel.setArg(1, bufferP);
    calcMaxRowKernel.setArg(2, bufferR);

    cl::Kernel calcMaxKernel(program, "utils_max");

    calcMaxKernel.setArg(0, nx);
    calcMaxKernel.setArg(1, bufferR);
    calcMaxKernel.setArg(2, bufferV);

    // cl::Kernel updateNUKernel = createUpdateNUKernel(ny, nz, bufferNU, program);
    // cl::Kernel updateQKernel = createUpdateQKernel(qxm, ny, bufferQX, program);

    cl::NDRange valueRange(1);
    cl::NDRange rowRange(nx);
    cl::NDRange surfaceRange(nx, ny);
    cl::NDRange volumeRange(nx, ny, nz);
    cl::NDRange nuRange(nx, ny, nz + 1);

    long long calcTime = 0;

    auto start = high_resolution_clock::now();            

    while (t <= endTime) {
        updateUAKernel.setArg(5, dt);
        updateUAKernel.setArg(11, bufferUA);
        updateUAKernel.setArg(15, bufferU1A);

        err = queue.enqueueNDRangeKernel(updateUAKernel, cl::NullRange, surfaceRange, cl::NullRange);

        updateVAKernel.setArg(5, dt);
        updateVAKernel.setArg(12, bufferVA);
        updateVAKernel.setArg(15, bufferV1A);

        err = queue.enqueueNDRangeKernel(updateVAKernel, cl::NullRange, surfaceRange, cl::NullRange);
        
        calcRHSKernel.setArg(5, dt);
        calcRHSKernel.setArg(13, bufferUA);
        calcRHSKernel.setArg(14, bufferU1A);
        calcRHSKernel.setArg(15, bufferVA);
        calcRHSKernel.setArg(16, bufferV1A);

        err = queue.enqueueNDRangeKernel(calcRHSKernel, cl::NullRange, volumeRange, cl::NullRange);      
        
        createTridiagonalMatrixKernel.setArg(1, dt);

        err = queue.enqueueNDRangeKernel(createTridiagonalMatrixKernel, cl::NullRange, volumeRange, cl::NullRange);

        err = queue.enqueueNDRangeKernel(updateUVFKernel, cl::NullRange, surfaceRange, cl::NullRange);

        t += dt;
        n += 1;

        /* updateQKernel.setArg(0, qxm);
        updateQKernel.setArg(2, bufferQX);

        err = queue.enqueueNDRangeKernel(updateQKernel, cl::NullRange, surfaceRange, cl::NullRange);

        updateQKernel.setArg(0, qym);
        updateQKernel.setArg(2, bufferQY);

        err = queue.enqueueNDRangeKernel(updateQKernel, cl::NullRange, surfaceRange, cl::NullRange);

        err = queue.enqueueNDRangeKernel(updateNUKernel, cl::NullRange, nuRange, cl::NullRange); */       

        err = queue.enqueueNDRangeKernel(calcMaxPlaneKernel, cl::NullRange, surfaceRange, cl::NullRange);
        err = queue.enqueueNDRangeKernel(calcMaxRowKernel, cl::NullRange, rowRange, cl::NullRange);
        err = queue.enqueueNDRangeKernel(calcMaxKernel, cl::NullRange, valueRange, cl::NullRange);

        err = queue.enqueueReadBuffer(bufferV, CL_TRUE, 0, sizeof(double), &nuMax);

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

            err = queue.enqueueReadBuffer(bufferUA, CL_TRUE, 0, sizeof(double) * nx * ny, ua.data());
            err = queue.enqueueReadBuffer(bufferVA, CL_TRUE, 0, sizeof(double) * nx * ny, va.data());
            err = queue.enqueueReadBuffer(bufferZ, CL_TRUE, 0, sizeof(double) * nx * ny, z.data());

            err = queue.enqueueReadBuffer(bufferQX, CL_TRUE, 0, sizeof(double) * nx * ny, qx.data());
            err = queue.enqueueReadBuffer(bufferQY, CL_TRUE, 0, sizeof(double) * nx * ny, qy.data());

            err = queue.enqueueReadBuffer(bufferUF, CL_TRUE, 0, sizeof(double) * nx * ny * nz, uf.data());
            err = queue.enqueueReadBuffer(bufferVF, CL_TRUE, 0, sizeof(double) * nx * ny * nz, vf.data());

            err = queue.enqueueReadBuffer(bufferNU, CL_TRUE, 0, sizeof(double) * nx * ny * (nz + 1), nu.data());

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

        swap(bufferU1A, bufferUA);
        swap(bufferV1A, bufferVA);

        updateZKernel.setArg(6, bufferUA);
        updateZKernel.setArg(7, bufferVA);

        queue.enqueueNDRangeKernel(updateZKernel, cl::NullRange, surfaceRange, cl::NullRange);
    }

    queue.finish();

    writeStatistics(statistics, outDir);

    cout << "Total number of iterations: " << n;
}
