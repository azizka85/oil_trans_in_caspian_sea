#include <iostream>

#include <format>

#include <chrono>

#include <stdexcept>

#include <utils/fs.h>
#include <utils/bathymetry.h>
#include <utils/viscosity.h>
#include <utils/calc.h>
#include <utils/opencl.h>

#include "writers/height.h"
#include "writers/surface.h"
#include "writers/volume.h"
#include "writers/statistics.h"

#include "solver.h"

using namespace std::chrono;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC;

Solver::Solver(
    float b,
    float f, float g, float rho, float kb,    
    float dx, float dy,
    float endTime, float outputTimeStep, string outDir,
    unique_ptr<Generators::Area::IGenerator> areaGenerator,
    unique_ptr<Generators::DZ::IGenerator> dzGenerator,
    unique_ptr<Generators::Bathymetry::IGenerator> hGenerator,
    unique_ptr<Generators::Wind::IGenerator> qGenerator,
    unique_ptr<Generators::Viscosity::IGenerator> nuGenerator
) {
    setB(b);

    setF(f);
    setG(g);
    setRHO(rho);
    setKB(kb);    

    setDX(dx);
    setDY(dy);

    setEndTime(endTime);
    setOutputTimeStep(outputTimeStep);
    setOutDir(outDir);

    setAreaGenerator(move(areaGenerator));
    setDZGenerator(move(dzGenerator));  
    setHGenerator(move(hGenerator));
    setQGenerator(move(qGenerator));
    setNUGenerator(move(nuGenerator));
}

float Solver::getB() {
    return b;
}

void Solver::setB(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("B should be > 0, but it is {}", val)
        );
    }

    b = val;
}

float Solver::getF() {
    return f;
}

void Solver::setF(float val) {
    f = val;
}

float Solver::getG() {
    return g;
}

void Solver::setG(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("G should be > 0, but it is {}", val)
        );
    }

    g = val;
}

float Solver::getRHO() {
    return rho;
}

void Solver::setRHO(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("RHO should be > 0, but it is {}", val)
        );
    }

    rho = val;
}

float Solver::getKB() {
    return kb;
}

void Solver::setKB(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("KB should be > 0, but it is {}", val)
        );
    }

    kb = val;
}

float Solver::getDX() {
    return dx;
}

void Solver::setDX(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("DX should be > 0, but it is {}", val)
        );
    }

    dx = val;
}

float Solver::getDY() {
    return dy;
}

void Solver::setDY(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("DY should be > 0, but it is {}", val)
        );
    }

    dy = val;
}

float Solver::getEndTime() {
    return endTime;
}

void Solver::setEndTime(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("endTime should be > 0, but it is {}", val)
        );
    }

    endTime = val;
}

float Solver::getOutputTimeStep() {
    return outputTimeStep;
}

void Solver::setOutputTimeStep(float val) {
    if (val <= 0) {
        throw runtime_error(
            format("outputTimeStep should be > 0, but it is {}", val)
        );
    }

    outputTimeStep = val;
}

string Solver::getOutDir() {
    return outDir;
}

void Solver::setOutDir(string val) {
    if (val.empty()) {
        throw runtime_error(
            format("outDir should not be empty, but it is {}", val)
        );
    }

    outDir = val;
}

void Solver::setAreaGenerator(unique_ptr<Generators::Area::IGenerator> val) {
    if (!val) {
        throw invalid_argument(
            "A null pointer (nullptr) was received instead of a valid Area generator."
        );
    }

    areaGenerator = move(val);
}

void Solver::setDZGenerator(unique_ptr<Generators::DZ::IGenerator> val) {
    if (!val) {
        throw invalid_argument(
            "A null pointer (nullptr) was received instead of a valid DZ generator."
        );
    }

    dzGenerator = move(val);
}

void Solver::setHGenerator(unique_ptr<Generators::Bathymetry::IGenerator> val) {
    if (!val) {
        throw invalid_argument(
            "A null pointer (nullptr) was received instead of a valid Bathymetry generator."
        );
    }

    hGenerator = move(val);
}

void Solver::setQGenerator(unique_ptr<Generators::Wind::IGenerator> val) {
    if (!val) {
        throw invalid_argument(
            "A null pointer (nullptr) was received instead of a valid Wind generator."
        );
    }

    qGenerator = move(val);
}

void Solver::setNUGenerator(unique_ptr<Generators::Viscosity::IGenerator> val) {
    if (!val) {
        throw invalid_argument(
            "A null pointer (nullptr) was received instead of a valid Viscosity generator."
        );
    }

    nuGenerator = move(val);
}

Directories Solver::createDirectory() {
    auto dirPath = Utils::FS::createDirByPath(
        path(
            format("{}/f={}, g={}, rho={}, kb={}", outDir, f, g, rho, kb)
        )
    );

    dirPath = hGenerator->createDirectory(dirPath);
    dirPath = nuGenerator->createDirectory(dirPath);
    dirPath = qGenerator->createDirectory(dirPath);
    dirPath = areaGenerator->createDirectory(dirPath);

    dirPath = Utils::FS::createDirByPath(
        dirPath / format("dx={}, dy={}", dx, dy)
    );

    dirPath = dzGenerator->createDirectory(dirPath);

    auto surfacePath = Utils::FS::createDirByPath(
        dirPath / path("surface")
    );

    auto volumePath = Utils::FS::createDirByPath(
        dirPath / path("volume")
    );

    auto viscosityPath = Utils::FS::createDirByPath(
        dirPath / path("viscosity")
    );

    return Directories{
        .root = dirPath, 
        .surface = surfacePath, 
        .volume = volumePath, 
        .viscosity = viscosityPath 
    };
}

void Solver::loadKernelSources(cl::Program::Sources& sources) {
    sources.push_back(
        Utils::OpenCL::loadKernelSource(
            "opencl/wind_induced_currents/davies85/non_periodic_bc/solver.cl"
        )
    );

    sources.push_back(
        Utils::OpenCL::loadKernelSource("opencl/utils/max.cl")
    );
}

void Solver::setInitialCondition(
    int nx, int ny, int nz,
    vector<float>& uf, vector<float>& vf,
    vector<float>& ua, vector<float>& va,
    vector<float>& z
) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int p = j + i * ny;

            for (int k = 0; k < nz; k++) {
                int id = k + p * nz;

                uf[id] = 0;
                vf[id] = 0;
            }

            ua[p] = 0;
            va[p] = 0;
            z[p] = 0;
        }
    }
}

void Solver::writeData(
    float t, int m,
    int nx, int ny, int nz,
    vector<float>& dz, 
    vector<float>& h, vector<float>& z,
    vector<float>& ua, vector<float>& va,
    vector<float>& uf, vector<float>& vf,
    vector<float>& qx, vector<float>& qy,
    vector<float>& nu, Directories &dirs
) {
    Writers::Surface::write(t, m, nx, ny, dx, dy, ua, va, qx, qy, z, dirs.surface);
    Writers::Volume::write(t, m, nx, ny, nz, dx, dy, dz, h, ua, va, uf, vf, dirs.volume);
    Writers::Volume::writeViscosity(t, m, nx, ny, nz + 1, dx, dy, dz, h, nu, dirs.viscosity);
}

void Solver::solve() {
    auto dirs = createDirectory();

    auto geo = areaGenerator->generate();

    int nx = static_cast<int>(
        ceil(geo.l / dx)
    ) + 1;

    int ny = static_cast<int>(
        ceil(geo.w / dy)
    ) + 1;

    vector<float> dz = dzGenerator->generateDZ();

    int nz = dz.size();

    vector<float> h = hGenerator->generateH(nx, ny);

    if (h.size() == 0) {
        throw runtime_error(
            format("The size of the heights array h must be > 0, but it is {}", h.size())
        );
    }

    auto [hMin, hMax] = Utils::Bathymetry::minMaxH(h, nx, ny);

    float dzMin = *min_element(dz.begin(), dz.end());

    const float dtMax = min(dx, dy) / sqrt(2 * g * hMax) / 1.5;

    vector<float> ua(nx * ny);
    vector<float> u1a(nx * ny);

    vector<float> va(nx * ny);
    vector<float> v1a(nx * ny);

    vector<float> z(nx * ny);

    vector<float> uf(nx * ny * nz);
    vector<float> ud(nx * ny * nz);

    vector<float> vf(nx * ny * nz);
    vector<float> vd(nx * ny * nz);

    setInitialCondition(
        nx, ny, nz, 
        uf, vf, 
        ua, va, 
        z
    );

    vector<Generators::Wind::Data> windData = qGenerator->generate(nx, ny);    

    if (windData.size() == 0) {
        throw runtime_error(
            format("The size of the wind data array windData must be > 0, but it is {}", windData.size())
        );
    }

    int currentWindIndex = 0;
    Generators::Wind::Data currentWindData = windData[currentWindIndex];    

    vector<float> nu = nuGenerator->generateNU(
        nx, ny, nz + 1,
        dz, h,
        currentWindData.u10, currentWindData.v10,
        currentWindData.qx, currentWindData.qy,
        ua, va
    );

    float nuMax = Utils::Viscosity::maxNU(nu, nx, ny, nz + 1);

    float dt = min(dtMax, hMin * dzMin * hMin * dzMin / 2 / nuMax);
    float dtp = dt;

    vector<float> zp(nx * ny);

    vector<float> up(nx * ny * nz);
    vector<float> vp(nx * ny * nz);

    Utils::Calc::updateData(nx, ny, zp, z);
    Utils::Calc::updateData(nx, ny, nz, up, uf);
    Utils::Calc::updateData(nx, ny, nz, vp, vf);

    float t = 0;
    float tn = outputTimeStep;

    int n = 1;
    int m = 0;

    vector<tuple<int, float, long long, float, float, float>> statistics;

    Writers::Height::write(h, nx, ny, dx, dy, dirs.root);

    writeData(
        t, m, nx, ny, nz, 
        dz, h, z, 
        ua, va, uf, vf, 
        currentWindData.qx, currentWindData.qy,
        nu, dirs
    );

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

    if (err == -11) {
        auto build_log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>();

        bool ok = true;
    }

    cl::Buffer bufferV(context, CL_MEM_READ_WRITE, sizeof(float));
    cl::Buffer bufferR(context, CL_MEM_READ_WRITE, sizeof(float) * nx);
    cl::Buffer bufferP(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);

    cl::Buffer bufferDZ(context, CL_MEM_READ_WRITE, sizeof(float) * nz);

    cl::Buffer bufferH(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);
    cl::Buffer bufferQX(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);
    cl::Buffer bufferQY(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);
    cl::Buffer bufferZ(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);

    cl::Buffer bufferUA(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);
    cl::Buffer bufferVA(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);
    cl::Buffer bufferU1A(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);
    cl::Buffer bufferV1A(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny);

    cl::Buffer bufferNU(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny * (nz + 1));

    cl::Buffer bufferUF(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny * nz);
    cl::Buffer bufferVF(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny * nz);

    cl::Buffer bufferUD(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny * nz);
    cl::Buffer bufferVD(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny * nz);

    cl::Buffer bufferAL(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny * nz);
    cl::Buffer bufferAC(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny * nz);
    cl::Buffer bufferAR(context, CL_MEM_READ_WRITE, sizeof(float) * nx * ny * nz);

    cl::CommandQueue queue(context, device);

    err = queue.enqueueWriteBuffer(bufferDZ, CL_TRUE, 0, sizeof(float) * nz, dz.data());

    err = queue.enqueueWriteBuffer(bufferH, CL_TRUE, 0, sizeof(float) * nx * ny, h.data());
    err = queue.enqueueWriteBuffer(bufferQX, CL_TRUE, 0, sizeof(float) * nx * ny, currentWindData.qx.data());
    err = queue.enqueueWriteBuffer(bufferQY, CL_TRUE, 0, sizeof(float) * nx * ny, currentWindData.qy.data());
    err = queue.enqueueWriteBuffer(bufferZ, CL_TRUE, 0, sizeof(float) * nx * ny, z.data());

    err = queue.enqueueWriteBuffer(bufferUA, CL_TRUE, 0, sizeof(float) * nx * ny, ua.data());
    err = queue.enqueueWriteBuffer(bufferVA, CL_TRUE, 0, sizeof(float) * nx * ny, va.data());

    err = queue.enqueueWriteBuffer(bufferNU, CL_TRUE, 0, sizeof(float) * nx * ny * (nz + 1), nu.data());

    err = queue.enqueueWriteBuffer(bufferUF, CL_TRUE, 0, sizeof(float) * nx * ny * nz, uf.data());
    err = queue.enqueueWriteBuffer(bufferVF, CL_TRUE, 0, sizeof(float) * nx * ny * nz, vf.data());

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
    updateUVFKernel.setArg(2, bufferH);
    updateUVFKernel.setArg(3, bufferAL);
    updateUVFKernel.setArg(4, bufferAC);
    updateUVFKernel.setArg(5, bufferAR);
    updateUVFKernel.setArg(6, bufferUF);
    updateUVFKernel.setArg(7, bufferVF);
    updateUVFKernel.setArg(8, bufferUD);
    updateUVFKernel.setArg(9, bufferVD);

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

        swap(bufferU1A, bufferUA);
        swap(bufferV1A, bufferVA);

        updateZKernel.setArg(6, bufferUA);
        updateZKernel.setArg(7, bufferVA);

        queue.enqueueNDRangeKernel(updateZKernel, cl::NullRange, surfaceRange, cl::NullRange);

        t += dt;
        n += 1;     

        err = queue.enqueueNDRangeKernel(calcMaxPlaneKernel, cl::NullRange, surfaceRange, cl::NullRange);
        err = queue.enqueueNDRangeKernel(calcMaxRowKernel, cl::NullRange, rowRange, cl::NullRange);
        err = queue.enqueueNDRangeKernel(calcMaxKernel, cl::NullRange, valueRange, cl::NullRange);

        err = queue.enqueueReadBuffer(bufferV, CL_TRUE, 0, sizeof(float), &nuMax);

        double dt1 = min(dtMax, hMin * dzMin * hMin * dzMin / 2 / nuMax);

        if (dt1 < dtp) {
            dtp = dt1;
            dt = Utils::Calc::adjustTimeStep(b, t, dt1, outputTimeStep, dtMax, false);
        }
        else {
            dt = Utils::Calc::adjustTimeStep(b, t, dt, outputTimeStep, dtMax, true);
        }        

        if (currentWindIndex < windData.size() - 1 && t >= windData[currentWindIndex + 1].time) {
            currentWindIndex += 1;
            currentWindData = windData[currentWindIndex];

            err = queue.enqueueReadBuffer(bufferUA, CL_TRUE, 0, sizeof(float) * nx * ny, ua.data());
            err = queue.enqueueReadBuffer(bufferVA, CL_TRUE, 0, sizeof(float) * nx * ny, va.data());

            nu = nuGenerator->generateNU(
                nx, ny, nz + 1, dz, h, 
                currentWindData.u10, currentWindData.v10,
                currentWindData.qx, currentWindData.qy,
                ua, va
            );

            err = queue.enqueueWriteBuffer(bufferNU, CL_TRUE, 0, sizeof(float) * nx * ny * (nz + 1), nu.data());
        }

        if (t >= tn) {
            auto end = high_resolution_clock::now();

            auto duration = duration_cast<milliseconds>(end - start).count();

            calcTime += duration;

            err = queue.enqueueReadBuffer(bufferUA, CL_TRUE, 0, sizeof(float) * nx * ny, ua.data());
            err = queue.enqueueReadBuffer(bufferVA, CL_TRUE, 0, sizeof(float) * nx * ny, va.data());
            err = queue.enqueueReadBuffer(bufferZ, CL_TRUE, 0, sizeof(float) * nx * ny, z.data());

            err = queue.enqueueReadBuffer(bufferUF, CL_TRUE, 0, sizeof(float) * nx * ny * nz, uf.data());
            err = queue.enqueueReadBuffer(bufferVF, CL_TRUE, 0, sizeof(float) * nx * ny * nz, vf.data());

            writeData(
                t, m, nx, ny, nz, 
                dz, h, z, 
                ua, va, uf, vf, 
                currentWindData.qx, currentWindData.qy, 
                nu, dirs
            );

            auto umd = Utils::Calc::maxAbsDifference(nx, ny, nz, up, uf);
            auto vmd = Utils::Calc::maxAbsDifference(nx, ny, nz, vp, vf);
            auto zmd = Utils::Calc::maxAbsDifference(nx, ny, zp, z);

            cout << format(
                "Write data in file t={:.3f}, convergence of u={:.5f}, v={:.5f}, z={:.7f} with dt={:.5}, calc time={}",
                t, umd, vmd, zmd, dt, calcTime / 1000
            ) << endl;

            statistics.push_back({ n, tn, calcTime / 1000, umd, vmd, zmd });

            Utils::Calc::updateData(nx, ny, zp, z);
            Utils::Calc::updateData(nx, ny, nz, up, uf);
            Utils::Calc::updateData(nx, ny, nz, vp, vf);

            m += 1;
            tn = t + outputTimeStep;

            start = high_resolution_clock::now();
        }                
    }

    queue.finish();

    Writers::Statistics::write(statistics, dirs.root);

    cout << "Total number of iterations: " << n;
}