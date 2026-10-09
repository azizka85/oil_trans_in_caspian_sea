#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_BATHYMETRY_GEBCO_GENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_BATHYMETRY_GEBCO_GENERATOR_H

#include <string>

#include "igenerator.h"

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Bathymetry {
	class GEBCOGenerator : public IGenerator {
	private:
		double latMin;
		double latMax;

		double lonMin;
		double lonMax;

		double refDepth;
		double minDepth;

		string filePath;

	public:
		GEBCOGenerator(
			double latMin, double latMax,
			double lonMin, double lonMax,
			double refDepth, double minDepth,
			string filePath
		);

		double getLatMin();
		double getLatMax();
		void setLatMinMax(double latMin, double latMax);

		double getLonMin();
		double getLonMax();
		void setLonMinMax(double lonMin, double lonMax);

		double getRefDepth();
		void setRefDepth(double val);

		double getMinDepth();
		void setMinDepth(double val);

		string getFilePath();
		void setFilePath(string val);

		path createDirectory(path outDir) override;

		vector<double> generateH(int nx, int ny) override;
	};
}

#endif