#include <format>

#include <cmath>

#include <fstream>
#include <iostream>

#include <random>

#include <stdexcept>

#include "solver.h"

using namespace Diffusion::Lagrangian;

Solver::Solver(
    double l, double a,
    int Np,
    double D, double r,
    double dx,
    double endTime, double outputTimeStep, string dir
) {
    setL(l);
	setA(a);

	setNp(Np);

    setD(D);
    setR(r);

    setDX(dx);

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

int Solver::getNp() {
    return Np;
}

void Solver::setNp(int val) {
    if (val <= 0) {
        throw runtime_error(
            format("Np should be > 0, but it is {}", val)
        );
    }

    Np = val;
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
        format("{}/D={}, L={}, Np={}/dx={}, r={}", dir, D, l, Np, dx, r)
    );

    create_directories(dirPath);

    return dirPath;
}

tuple<int, vector<double>> Solver::setInitialCondition(int nx, vector<int>& C) {
    int N = 0;    

	vector<double> xp;

    for (int i = 0; i < nx; i++) {
        double x = dx * i;

        if (x <= a) {
            if (i == 0 || i == nx - 1) {
                C[i] = Np / 2;
				N += Np / 2;

                if (i == nx - 1) {
					x -= dx / 2;
                }

                for (int k = 0; k < Np / 2; k++) {
                    xp.push_back(x + k * dx / Np);
                }
            }
            else {
				C[i] = Np;
				N += Np;

                for (int k = 0; k < Np; k++) {
                    xp.push_back(x - 0.5 * dx + k * dx / Np);                    
                }
            }
        } else {
            C[i] = 0.;
        }
	}

    return {N, xp};
}

int Solver::maxAbsDifference(int nx, vector<int>& C, vector<int>& C1) {
    int maxDiff = 0.;

    for (int i = 0; i < nx; i++) {
        maxDiff = max(
            maxDiff,
            abs(C[i] - C1[i])
        );
    }

    return maxDiff;
}

void Solver::updateData(int N, int nx, vector<double>& xp, vector<int>& C, vector<double>& xp1, vector<int>& C1) {
    for (int i = 0; i < nx; i++) {
        C[i] = C1[i];
    }

    for (int i = 0; i < N; i++) {
        xp[i] = xp1[i];
	}
}

void Solver::writeData(vector<int>& C, double t, int nx, int m, path outDir) {
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
    file << format("DIMENSIONS {} 1 1", nx) << endl;
    file << format("POINTS {} float", nx) << endl;

    for (int i = 0; i < nx; i++) {
        double x = dx * i;

        file << format("{:.3f} 0.0 0.0", x) << endl;
    }

    file << "FIELD FieldData 1" << endl;
    file << "Time 1 1 float" << endl;
    file << format("{:.3f}", t) << endl;
    file << format("POINT_DATA {}", nx) << endl;
    file << "SCALARS C float" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < nx; i++) {
		int N = i == 0 || i == nx - 1 ? Np / 2 : Np;

        file << format("{:.3f}", (C[i] * 1.0) / N) << endl;
    }
}

void Solver::writeStatistics(vector<tuple<int, double, int>>& statistics, path outDir) {
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

    for (auto& t : statistics) {
        file << format(
            "{},  {:.3f},  {}",
            get<0>(t),
            get<1>(t),
            get<2>(t)
        ) << endl;
    }
}

void Solver::solve() {
    auto outDir = createDirectory();

    double dt = r * dx * dx / D;

    int nx = static_cast<int>(
        floor(l / dx)
    ) + 1;

    vector<int> C(nx);
    vector<int> C1(nx);
    vector<int> Cp(nx);    	

    auto [N, xp] = setInitialCondition(nx, C);

    vector<double> xp1(N);

    updateData(N, nx, xp, Cp, xp, C);

    const double mean = 0.0;
    const double std_dev = 1.0;

	mt19937 rd = default_random_engine();
	normal_distribution<double> dist(mean, std_dev);

    double t = 0;
    double tn = outputTimeStep;

    int n = 1;
    int m = 0;

    vector<tuple<int, double, int>> statistics;

    writeData(C, t, nx, m, outDir);

    m += 1;

    while (t <= endTime) {
        for (int i = 0; i < nx; i++) {
			C1[i] = C[i];
        }

        for (int i = 0; i < N; i++) {
			double dw = dist(rd);

            double x = xp[i] + dw*sqrt(2*D*dt);

            if (x < 0) {
                x = -x;
            }
            else if (x > l) {
                x = 2 * l - x;
            }

            int j = static_cast<int>(floor((x + dx / 2) / dx));
            int k = static_cast<int>(floor((xp[i] + dx / 2) / dx));

            if (j != k && (((j == 0 || j == nx - 1) && C1[j] == Np/2) || C1[j] == Np)) {
                xp1[i] = xp[i];
            }
            else {
                xp1[i] = x;

                if (j != k) {
                    C1[j] += 1;
                    C1[k] -= 1;

                    if (C1[j] < 0 || C1[k] < 0 || C1[j] > Np || C1[k] > Np) {
                        cout << format("{:.3f}, {}: Error: C1[{}]={}, C1[{}]={}", t + dt, n + 1, j, C1[j], k, C1[k]) << endl;
                    }
                }
            }            
        }

        updateData(N, nx, xp, C, xp1, C1);

        t += dt;

        if (t >= tn) {
            writeData(C, t, nx, m, outDir);

            int maxDiff = maxAbsDifference(nx, C, Cp);

            cout << format("Write data in file t={:.3f}, convergence of C={}", t, maxDiff) << endl;

            statistics.push_back({ n, tn, maxDiff });

            updateData(N, nx, xp, Cp, xp, C);

            m += 1;

            tn += outputTimeStep;
        }

        n += 1;
    }

    writeStatistics(statistics, outDir);
}