#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_BATHYMETRY_GEBCO_GENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_BATHYMETRY_GEBCO_GENERATOR_H

#include <string>

#include "igenerator.h"

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Bathymetry {
	class GEBCOGenerator : public IGenerator {
	private:
		float latMin;
		float latMax;

		float lonMin;
		float lonMax;

		float refDepth;
		float minDepth;

		string filePath;

	public:
		GEBCOGenerator(
			float latMin, float latMax,
			float lonMin, float lonMax,
			float refDepth, float minDepth,
			string filePath
		);

		float getLatMin();
		float getLatMax();
		void setLatMinMax(float latMin, float latMax);

		float getLonMin();
		float getLonMax();
		void setLonMinMax(float lonMin, float lonMax);

		float getRefDepth();
		void setRefDepth(float val);

		float getMinDepth();
		void setMinDepth(float val);

		string getFilePath();
		void setFilePath(string val);

		path createDirectory(path outDir) override;

		vector<float> generateH(int nx, int ny) override;
	};
}

#endif