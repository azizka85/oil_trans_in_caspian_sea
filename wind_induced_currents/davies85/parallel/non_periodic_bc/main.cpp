#include <string>

#include <memory>

#include <iostream>

#include <stdexcept>

#include "generators/area/proj_generator.h"
#include "generators/dz/triple_point_generator.h"
#include "generators/bathymetry/gebco_generator.h"
#include "generators/wind/ecmwf_generator.h"
#include "generators/viscosity/wind_speed_generator.h"

#include "solver.h"

using namespace std;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC;
using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators;

int main() {
	const float b = 1.1;

	const float f = 1.2e-4;
	const float g = 9.81;
	const float rho = 1025;
	const float kb = 0.002;

	const float k0 = 0.4;
	const float k = 2e-5;
	const float sigma = 1.2e-4;

	const float rhoAir = 1.225;
	const float Cd = 2.5e-3;

	const float latMin = 36;
	const float latMax = 47;

	const float lonMin = 46;
	const float lonMax = 56;

	double refDepth = -25;
	double minDepth = 1;

	const float dx = 5000;
	const float dy = 5000;

	const float dzMin = 0.002;
	const float dzMax = 0.1;

	const float zm = 0.51;

	const float num = 0.4;
	const float nu0 = 1.15e-6;

	const float w = 260;
	const float l = 260;

	const float hm = 260;

	const float u10m = 22;
	const float v10m = 22;

	const float qxm = 1.5;
	const float qym = 1.5;

	const float endTime = 60000;
	const float outputTimeStep = 600;

	const string gebcoFilePath = "data/bathymetry/gebco_2026_n47.0_s36.0_w46.0_e56.0.nc";
	const string ecmwfFilePath = "data/wind/ecmwf_2026_07_31_n47_e56_s36_w46.nc";

	const string outDir = "out";

	try {
		auto areaGenerator = make_unique<Area::ProjGenerator>(latMin, latMax, lonMin, lonMax);
		auto dzGenerator = make_unique<DZ::TriplePointGenerator>(dzMin, dzMax, zm);

		auto hGenerator = make_unique<Bathymetry::GEBCOGenerator>(
			latMin, latMax, 
			lonMin, lonMax, 
			refDepth, minDepth, 
			gebcoFilePath
		);

		auto qGenerator = make_unique<Wind::ECMWFGenerator>(
			latMin, latMax, 
			lonMin, lonMax, 
			rhoAir, Cd, 
			ecmwfFilePath
		);

		auto nuGenerator = make_unique<Viscosity::WindSpeedGenerator>(k0, k, f, sigma, rho, nu0);

		Solver solver(
			b, f, g, rho, kb, 			
			dx, dy, 
			endTime, outputTimeStep, outDir, 
			move(areaGenerator),
			move(dzGenerator),			
			move(hGenerator),
			move(qGenerator),
			move(nuGenerator)
		);
	
		solver.solve();
	}
	catch (const exception& e) {
		cout << "Caught exception: " << e.what() << std::endl;
	}

	return 0;
}