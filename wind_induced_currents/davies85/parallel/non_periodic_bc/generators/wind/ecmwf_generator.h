#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_WIND_ECMWF_GENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_WIND_ECMWF_GENERATOR_H

#include <string>

#include "igenerator.h"

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Wind {
	class ECMWFGenerator : public IGenerator {
		private:
			double latMin;
			double latMax;

			double lonMin;
			double lonMax;

			double rhoAir;
			double Cd;

			string filePath;

		public:
			ECMWFGenerator(
				double latMin, double latMax,
				double lonMin, double lonMax,
				double rhoAir, double Cd,
				string filePath
			);

			double getLatMin();
			double getLatMax();
			void setLatMinMax(double latMin, double latMax);

			double getLonMin();
			double getLonMax();
			void setLonMinMax(double lonMin, double lonMax);

			double getRhoAir();
			void setRhoAir(double val);

			double getCd();
			void setCd(double val);

			string getFilePath();
			void setFilePath(string val);

			path createDirectory(path outDir) override;

			vector<Data> generate(int nx, int ny) override;
	};
}

#endif