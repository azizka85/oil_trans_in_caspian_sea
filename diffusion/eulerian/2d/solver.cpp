#include <format>

#include <cmath>

#include <fstream>
#include <iostream>

#include <stdexcept>

#include "solver.h"

using namespace Diffusion::Eulerian;

Solver::Solver(
	double l, double h,
	double a, double b,
    double D, double r,
	double dx, double dy,
    double endTime, double outputTimeStep, string dir
) {
    setL(l);
	setH(h);

	setA(a);
	setB(b);

    setD(D);
    setR(r);

    setDX(dx);
	setDY(dy);

    setEndTime(endTime);
    setOutputTimeStep(outputTimeStep);
	setDir(dir);
}

double Solver::getL() {
    return l;
}

void Solver::setL(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("l should be > 0, but it is {}", val)
        );
    }

    l = val;
}

double Solver::getH() {
    return h;
}

void Solver::setH(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("h should be > 0, but it is {}", val)
        );
    }

    h = val;
}

double Solver::getA() {
    return a;
}

void Solver::setA(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("a should be > 0, but it is {}", val)
        );
	}

    if (val > l) {
        throw runtime_error(
            format("a should be <= l, but it is {} and l={}", val, l)
        );
    }

    a = val;
}

double Solver::getB() {
    return b;
}

void Solver::setB(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("b should be > 0, but it is {}", val)
        );
    }

    if (val > h) {
        throw runtime_error(
            format("b should be <= h, but it is {} and h={}", val, h)
        );
    }

    b = val;
}

double Solver::getD() {
    return D;
}

void Solver::setD(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("D should be > 0, but it is {}", val)
        );
    }

    D = val;
}

double Solver::getR() {
    return r;
}

void Solver::setR(double val) {
    if (val <= 0) {
        throw runtime_error(
            format("r should be > 0, but it is {}", val)
        );
    }

    r = val;
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

    if (val > l) {
        throw runtime_error(
            format("DX should be <= l, but it is {} and l={}", val, l)
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

    if (val > h) {
        throw runtime_error(
            format("DY should be <= h, but it is {} and h={}", val, h)
        );
    }

    dy = val;
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

path Solver::createDirectory() {
    auto dirPath = path(
        format("{}/D={}, L={}, H={}/dx={}, dy={}, r={}", dir, D, l, h, dx, dy, r)
    );

    create_directories(dirPath);

    return dirPath;
}

void Solver::setInitialCondition(int nx, int ny, vector<vector<double>>& C) {
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double x = dx * i;
            double y = dy * j;

            C[i][j] = (x <= a && y <= b) ? 1. : 0.;
		}
	}
}

double Solver::maxAbsDifference(int nx, int ny, vector<vector<double>>& C, vector<vector<double>>& C1) {
    double maxDiff = 0.;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double diff = abs(C[i][j] - C1[i][j]);

            if (diff > maxDiff) {
                maxDiff = diff;
            }
		}
    }

    return maxDiff;
}

void Solver::updateData(int nx, int ny, vector<vector<double>>& C, vector<vector<double>>& C1) {
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            C[i][j] = C1[i][j];
		}
    }
}

void Solver::writeData(vector<vector<double>>& C, double t, int nx, int ny, int m, path outDir) {
    auto filePath = path(
        format("data.{:03}.vtk", m)
    );

    filePath = outDir / filePath;

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
    file << format("POINTS {} float", nx * ny) << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double x = dx * i;
            double y = dy * j;

            file << format("{:.3f} {:.3f} 0.000", x, y) << endl;
		}
    }

    file << "FIELD FieldData 1" << endl;
    file << "Time 1 1 float" << endl;
    file << format("{:.3f}", t) << endl;
    file << format("POINT_DATA {}", nx * ny) << endl;
    file << "SCALARS C float" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << format("{:.3f}", C[i][j]) << endl;
        }
    }
}

void Solver::writeStatistics(vector<tuple<int, double, double>>& statistics, path outDir) {
    auto dirPath = outDir / path("statistics");

    create_directory(dirPath);

    auto filePath = dirPath / path("convergence.csv");

    ofstream file(filePath);

    if (file.bad()) {
        throw runtime_error(
            format("Failed to open file at: {}", filePath.string())
        );
    }

    file << "n,  t,   max_diff" << endl;

    for (auto t : statistics) {
        file << format(
            "{},  {:.3f},  {:.5f}",
            get<0>(t),
            get<1>(t),
            get<2>(t)
        ) << endl;
    }
}

void Solver::solve() {
    auto outDir = createDirectory();

    double dt = r * dx * dx * dy * dy / D / (dx * dx + dy * dy);

    double rx = D * dt / dx / dx;
    double ry = D * dt / dy / dy;

    int nx = static_cast<int>(
        ceil(l / dx)
    ) + 1;

    int ny = static_cast<int>(
        ceil(h / dy)
    ) + 1;

    vector<vector<double>> C(nx, vector<double>(ny));
    vector<vector<double>> C1(nx, vector<double>(ny));
    vector<vector<double>> Cp(nx, vector<double>(ny));

    setInitialCondition(nx, ny, C);

    updateData(nx, ny, Cp, C);

    double t = 0;
    double tn = outputTimeStep;

    int n = 1;
    int m = 0;

    vector<tuple<int, double, double>> statistics;

    writeData(C, t, nx, ny, m, outDir);

    m += 1;

    while (t <= endTime) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                double Cl = i == 0 ? C[1][j] : C[i - 1][j];
                double Cr = i == nx - 1 ? C[nx - 2][j] : C[i + 1][j];

				double Cd = j == 0 ? C[i][1] : C[i][j - 1];
				double Cu = j == ny - 1 ? C[i][ny - 2] : C[i][j + 1];

                C1[i][j] = C[i][j] + rx * (Cr - 2 * C[i][j] + Cl) + ry * (Cu - 2 * C[i][j] + Cd);
            }
        }

        updateData(nx, ny, C, C1);

        t += dt;

        if (t >= tn) {
            writeData(C, t, nx, ny, m, outDir);

            auto maxDiff = maxAbsDifference(nx, ny, C, Cp);

            cout << format("Write data in file t={:.3f}, convergence of C={:.5f}", t, maxDiff) << endl;

            statistics.push_back({ n, tn, maxDiff });

            updateData(nx, ny, Cp, C);

            m += 1;

            tn += outputTimeStep;
        }

        n += 1;
    }

    writeStatistics(statistics, outDir);
}