#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_AREA_PROJ_GENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_AREA_PROJ_GENERATOR_H

#include "igenerator.h"

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Area {
	class ProjGenerator : public IGenerator {
		private:
			float latMin;
			float latMax;

			float lonMin;
			float lonMax;

		public:
			ProjGenerator(
				float latMin, float latMax,
				float lonMin, float lonMax
			);

			float getLatMin();
			float getLatMax();
			void setLatMinMax(float latMin, float latMax);

			float getLonMin();
			float getLonMax();
			void setLonMinMax(float lonMin, float lonMax);

			path createDirectory(path outDir) override;

			virtual Geometry generate() override;
	};
}

#endif