#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_WIND_ECMWF_GENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_WIND_ECMWF_GENERATOR_H

#include <string>

#include "igenerator.h"

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Wind {
	class ECMWFGenerator : public IGenerator {
		private:
			float latMin;
			float latMax;

			float lonMin;
			float lonMax;

			float rhoAir;
			float Cd;

			string filePath;

		public:
			ECMWFGenerator(
				float latMin, float latMax,
				float lonMin, float lonMax,
				float rhoAir, float Cd,
				string filePath
			);

			float getLatMin();
			float getLatMax();
			void setLatMinMax(float latMin, float latMax);

			float getLonMin();
			float getLonMax();
			void setLonMinMax(float lonMin, float lonMax);

			float getRhoAir();
			void setRhoAir(float val);

			float getCd();
			void setCd(float val);

			string getFilePath();
			void setFilePath(string val);

			path createDirectory(path outDir) override;

			vector<Data> generate(int nx, int ny) override;
	};
}

#endif