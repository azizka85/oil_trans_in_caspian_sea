#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_AREA_PROJ_GENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_AREA_PROJ_GENERATOR_H

#include "igenerator.h"

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Area {
	class ProjGenerator : public IGenerator {
		private:
			double latMin;
			double latMax;

			double lonMin;
			double lonMax;

		public:
			ProjGenerator(
				double latMin, double latMax,
				double lonMin, double lonMax
			);

			double getLatMin();
			double getLatMax();
			void setLatMinMax(double latMin, double latMax);

			double getLonMin();
			double getLonMax();
			void setLonMinMax(double lonMin, double lonMax);

			path createDirectory(path outDir) override;

			virtual Geometry generate() override;
	};
}

#endif