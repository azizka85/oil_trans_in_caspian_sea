#include <iostream>

#include "solver.h"

using namespace std;

using namespace Diffusion::Lagrangian;

int main() {
    const double l = 1;
	const double a = 0.1;

	const int Np = 10;

    const double D = 1;

    const double dx = 0.1;

	const double r = 0.0003;

    const double endTime = 0.03;
    const double outputTimeStep = 0.0003;

    const string dir = "data";


    try {
        Solver solver(l, a, Np, D, r, dx, endTime, outputTimeStep, dir);

        solver.solve();
    } catch (const exception& e) {
        cout << "Caught exception: " << e.what() << std::endl;
    }

    return 0;
}
