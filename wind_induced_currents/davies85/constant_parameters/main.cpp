#include <iostream>

#include "solver.h"

using namespace std;

using namespace WindInducedCurrents::Davies85::ConstantParameters;

int main() {
    const double f = 1.2e-4;
    const double b = 1.1;

    const double g = 9.81;
    const double rho = 1025;
    const double nu = 0.4;

    const double h = 260;
    const double w = 260;
    const double l = 260;

    const double kb = 0.002;
    const double qx = 1.5;
    const double qy = 1.5;

    const double dx = 13;
    const double dy = 13;
    const double dz = 0.05;

    const double endTime = 60000;
    const double outputTimeStep = 600;

    const string dir = "data";


    try {
        Solver solver(f, b, h, w, l, g, rho, nu, kb, qx, qy, dx, dy, dz, endTime, outputTimeStep, dir);

        solver.solve();
    } catch (const exception& e) {
        cout << "Caught exception: " << e.what() << std::endl;
    }

    return 0;
}
