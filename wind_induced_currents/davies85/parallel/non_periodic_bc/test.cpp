#include <string>
#include <vector>

#include <ranges>

#include <chrono>

#include <format>

#include <memory>

#include <iostream>
#include <fstream>

#include <filesystem>

#include <netcdf>

#include <proj.h>

#include <boost/math/interpolators/bilinear_uniform.hpp>

#include "generators/default_generator.h"

using namespace std;
using namespace std::filesystem;

using namespace netCDF;

using namespace boost::math::interpolators;

int main() {
	double latMin = 36;
	double latMax = 47;

	double lonMin = 46;
	double lonMax = 56;

	double dx = 5000;
	double dy = 5000;

	double refDepth = -25;
	double minDepth = -1;

	string filePath = "data/bathymetry/gebco_2026_n47.0_s36.0_w46.0_e56.0.nc";

	NcFile bathymetryFile(filePath, NcFile::read);

	NcVar latVar = bathymetryFile.getVar("lat");
	size_t latSize = latVar.getDim(0).getSize();
	vector<double> lats(latSize);

	latVar.getVar(lats.data());	

	double latStep = lats.size() > 1 ? lats[1] - lats[0] : 0;

	NcVar lonVar = bathymetryFile.getVar("lon");
	size_t lonSize = lonVar.getDim(0).getSize();
	vector<double> lons(lonSize);

	lonVar.getVar(lons.data());

	double lonStep = lons.size() > 1 ? lons[1] - lons[0] : 0;

	NcVar elevVar = bathymetryFile.getVar("elevation");
	vector<double> elevations(latSize * lonSize);

	elevVar.getVar(elevations.data());	

	unique_ptr<PJ_CONTEXT, PJ_CONTEXT* (*)(PJ_CONTEXT*)> ctx(
		proj_context_create(),
		proj_context_destroy
	);

	unique_ptr<PJ, PJ* (*)(PJ*)> transformer(
		proj_create_crs_to_crs(ctx.get(), "EPSG:4326", "EPSG:32639", nullptr),
		proj_destroy
	);	

	PJ_COORD coord = proj_coord(latMin, lonMin, 0, 0);	

	PJ_COORD utm = proj_trans(transformer.get(), PJ_FWD, coord);
	
	double minX = utm.xy.x;
	double minY = utm.xy.y;

	coord = proj_coord(latMin, lonMax, 0, 0);
	utm = proj_trans(transformer.get(), PJ_FWD, coord);

	double maxX = utm.xy.x;

	coord = proj_coord(latMax, lonMax, 0, 0);
	utm = proj_trans(transformer.get(), PJ_FWD, coord);

	double maxY = utm.xy.y;

	cout << "Min: " << minX << ", " << minY << endl;
	cout << "Max: " << maxX << ", " << maxY << endl;
	
	auto interpFunc = bilinear_uniform(
		move(elevations),
		latSize, lonSize,
		lonStep, latStep,
		lons.front(), lats.front()
	);

	int nx = static_cast<int>(
			ceil((maxX - minX) / dx)
	) + 1;

	int ny = static_cast<int>(
		ceil((maxY - minY) / dy)
	) + 1;

	cout << "nx: " << nx << ", ny: " << ny << endl;

	vector<double> depths(nx * ny);
	
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			int p = i + j * nx;

			double x = minX + dx * i;
			double y = minY + dy * j;
			
			PJ_COORD coord = proj_coord(x, y, 0, 0);
			PJ_COORD geo = proj_trans(transformer.get(), PJ_INV, coord);

			if (
				geo.lp.lam >= lats.front() && geo.lp.lam <= lats.back() &&
				geo.lp.phi >= lons.front() && geo.lp.phi <= lons.back()
			) {
				double elevation = interpFunc(geo.lp.phi, geo.lp.lam);

				if (elevation < refDepth + minDepth) {
					depths[p] = abs(elevation) + refDepth;
				}
			}			
		}
	}

	auto dirPath = path("data");

	create_directories(dirPath);

	auto outPath = dirPath / path("depths.vtk");

	ofstream file(outPath);

	if (file.bad()) {
		throw runtime_error(
			format("Failed to open file at: {}", outPath.string())
		);
	}

	file << "# vtk DataFile Version 3.0" << endl;
	file << format("TIME {:.3f}", 0.) << endl;
	file << "ASCII" << endl;
	file << "DATASET STRUCTURED_GRID" << endl;
	file << format("DIMENSIONS {} {} 1", nx, ny) << endl;
	file << format("POINTS {} float", nx * ny) << endl;

	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			float x = dx * i;
			float y = dy * j;

			file << format("{:.3f} {:.3f} 0.0", floor(x / 1000), floor(y / 1000)) << endl;
		}
	}

	file << "FIELD FieldData 1" << endl;
	file << "Time 1 1 float" << endl;
	file << format("{:.3f}", 0.) << endl;
	file << format("POINT_DATA {}", nx * ny) << endl;

	file << "SCALARS h float" << endl;
	file << "LOOKUP_TABLE default" << endl;

	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			int p = i + j * nx;

			file << format("{:.3f}", depths[p]) << endl;
		}
	}

	filePath = "data/wind/ecmwf_2026_07_31_n47_e56_s36_w46.nc";

	NcFile windFile(filePath, NcFile::read);

	NcVar timeVar = windFile.getVar("valid_time");
	size_t timeSize = timeVar.getDim(0).getSize();	

	vector<int64_t> times(timeSize);

	timeVar.getVar(times.data());

	std::string timeUnits;

	timeVar.getAtt("units").getValues(timeUnits);

	cout << "time units: " << timeUnits << endl;

	auto pipeline = times |
		views::transform([](int64_t seconds) {
				auto duration = chrono::seconds(seconds);

				chrono::sys_time<chrono::seconds> timePoint(duration);

				std::time_t t = chrono::system_clock::to_time_t(timePoint);
				std::tm tm = *std::gmtime(&t);
				std::ostringstream oss;
				
				oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
				
				return oss.str();
			}
		);

	vector<string> utcTimes(pipeline.begin(), pipeline.end());

	latVar = windFile.getVar("latitude");
	latSize = latVar.getDim(0).getSize();
	lats = vector<double>(latSize);

	latVar.getVar(lats.data());

	latStep = lats.size() > 1 ? lats[0] - lats[1] : 0;

	lonVar = windFile.getVar("longitude");
	lonSize = lonVar.getDim(0).getSize();
	lons = vector<double>(lonSize);

	lonVar.getVar(lons.data());

	lonStep = lons.size() > 1 ? lons[1] - lons[0] : 0;	

	NcVar u10Var = windFile.getVar("u10");
	vector<double> u10Arr(timeSize* latSize* lonSize);

	u10Var.getVar(u10Arr.data());

	NcVar v10Var = windFile.getVar("v10");
	vector<double> v10Arr(timeSize * latSize * lonSize);

	v10Var.getVar(v10Arr.data());

	auto surfacePath = dirPath / path("surface");

	create_directories(surfacePath);

	for (int k = 0; k < timeSize; k++) {
		outPath = surfacePath / path(
			format("data.{:03}.vtk", k)
		);

		file = ofstream(outPath);

		if (file.bad()) {
			throw runtime_error(
				format("Failed to open file at: {}", outPath.string())
			);
		}

		auto t = times[k];

		file << "# vtk DataFile Version 3.0" << endl;
		file << format("TIME {:03}", t) << endl;
		file << "ASCII" << endl;
		file << "DATASET STRUCTURED_GRID" << endl;
		file << format("DIMENSIONS {} {} 1", nx, ny) << endl;
		file << format("POINTS {} float", nx * ny) << endl;

		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				float x = dx * i;
				float y = dy * j;

				file << format("{:.3f} {:.3f} 0.0", floor(x / 1000), floor(y / 1000)) << endl;
			}
		}

		file << "FIELD FieldData 1" << endl;
		file << "Time 1 1 string" << endl;
		file << format("{}", utcTimes[k]) << endl;
		file << format("POINT_DATA {}", nx * ny) << endl;

		file << "VECTORS VA float" << endl;

		double* u10SlicePtr = u10Arr.data() + k * latSize * lonSize;
		span<double> u10SliceSpan(u10SlicePtr, latSize * lonSize);

		auto u10Func = bilinear_uniform(
			move(u10SliceSpan),
			latSize, lonSize,
			lonStep, latStep,
			lons.front(), lats.back()
		);

		double* v10SlicePtr = v10Arr.data() + k * latSize * lonSize;	
		span<double> v10SliceSpan(v10SlicePtr, latSize * lonSize);

		auto v10Func = bilinear_uniform(
			move(v10SliceSpan),
			latSize, lonSize,
			lonStep, latStep,
			lons.front(), lats.back()
		);

		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				int p = i + j * nx;

				double x = minX + dx * i;
				double y = minY + dy * j;

				PJ_COORD coord = proj_coord(x, y, 0, 0);
				PJ_COORD geo = proj_trans(transformer.get(), PJ_INV, coord);

				double u10 = 0;
				double v10 = 0;

				if (
					geo.lp.lam >= lats.back() && geo.lp.lam <= lats.front() &&
					geo.lp.phi >= lons.front() && geo.lp.phi <= lons.back() &&
					depths[p] > 0
				) {
					u10 = u10Func(geo.lp.phi, geo.lp.lam);
					v10 = v10Func(geo.lp.phi, geo.lp.lam);					
				}

				file << format("{:.3f} {:.3f} 0", u10, v10) << endl;
			}
		}
	}

	cout << "Finish" << endl;

	return 0;
}