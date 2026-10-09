#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_BATHYMETRY_UNIFORM_GENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_BATHYMETRY_UNIFORM_GENERATOR_H

#include "igenerator.h"

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Bathymetry {
	class UniformGenerator : public IGenerator {
		private:
			double hm;

		public:
			UniformGenerator(double hm);

			double getHM();
			void setHM(double val);

			path createDirectory(path outDir) override;

			vector<double> generateH(int nx, int ny) override;
	};
}

#endif