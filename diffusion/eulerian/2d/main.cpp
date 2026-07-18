#include <iostream>

#include "solver.h"

using namespace std;

using namespace Diffusion::Eulerian;

int main() {
    const double l = 1;
	const double h = 1;

	const double a = 0.1;
	const double b = 0.1;

    const double D = 1;

    const double dx = 0.01;
	const double dy = 0.01;

	const double r = 0.06;

    const double endTime = 0.03;
    const double outputTimeStep = 0.0003;

    const string dir = "data";

    try {
        Solver solver(l, h, a, b, D, r, dx, dy, endTime, outputTimeStep, dir);

        solver.solve();
    } catch (const exception& e) {
        cout << "Caught exception: " << e.what() << std::endl;
    }

    return 0;
}
