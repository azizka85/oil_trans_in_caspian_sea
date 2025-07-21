#include <iostream>

#include "solver.h"

using namespace std;

using namespace WindInducedCurrents::Davies85::NonPeriodicBC;

int main() {
    const double f = 1.2e-4;
    const double b = 1.1;

    const double g = 9.81;
    const double rho = 1025;
    const double kb = 0.002;

    const double num = 0.4;
    const GenerateNU nug = GenerateNU::UniformNU;

    const double hm = 260;
    const GenerateH hg = GenerateH::UniformH;

    const double w = 260;
    const double l = 260;

    const double qxm = 1.5;
    const double qym = 1.5;
    const GenerateQ qg = GenerateQ::UniformQ;

    const double dx = 13;
    const double dy = 13;

    const double dzm = 0.1;
    const GenerateDZ dzg = GenerateDZ::ParabolicDZ;

    const double endTime = 60000;
    const double outputTimeStep = 600;

    const string dir = "data";

    try {
        Solver solver(f, b, hm, hg, w, l, g, rho, kb, num, nug, qxm, qym, qg, dx, dy, dzm, dzg, endTime, outputTimeStep, dir);

        solver.solve();
    } catch (const exception& e) {
        cout << "Caught exception: " << e.what() << std::endl;
    }

    return 0;
}
