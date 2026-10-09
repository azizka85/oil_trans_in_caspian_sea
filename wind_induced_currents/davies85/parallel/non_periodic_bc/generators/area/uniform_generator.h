#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_AREA_UNIFORM_GENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_AREA_UNIFORM_GENERATOR_H

#include "igenerator.h"

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Area {
	class UniformGenerator : public IGenerator {
		private:
			double w;
			double l;

		public:
			UniformGenerator(double w, double l);

			double getW();
			void setW(double val);

			double getL();
			void setL(double val);

			path createDirectory(path outDir) override;

			virtual Geometry generate() override;
	};
}

#endif