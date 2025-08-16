#include <iostream>

#include "solver.h"

using namespace std;

using namespace WindInducedCurrents::Davies85::Parallel::VariableParameters;

int main() {
    const float f = 1.2e-4;
    const float b = 1.1;

    const float g = 9.81;
    const float rho = 1025;
    const float kb = 0.002;

    const float num = 0.4;
    const GenerateNU nug = GenerateNU::UniformNU;

    const float hm = 260;
    const GenerateH hg = GenerateH::UniformH;

    const float w = 260;
    const float l = 260;

    const float qxm = 1.5;
    const float qym = 1.5;
    const GenerateQ qg = GenerateQ::UniformQ;

    const float dx = 13;
    const float dy = 13;

    const float dzm = 0.1;
    const GenerateDZ dzg = GenerateDZ::ParabolicDZ;

    const float endTime = 60000;
    const float outputTimeStep = 600;

    const string dir = "data";

    try {
        Solver solver(f, b, hm, hg, w, l, g, rho, kb, num, nug, qxm, qym, qg, dx, dy, dzm, dzg, endTime, outputTimeStep, dir);

        solver.solve();
    } catch (const exception& e) {
        cout << "Caught exception: " << e.what() << std::endl;
    }

    return 0;
}
