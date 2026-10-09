#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_DZ_UNIFORM_GENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_DZ_UNIFORM_GENERATOR_H

#include "igenerator.h"

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::DZ {
	class UniformGenerator : public IGenerator {
		private:
			double dzm;

		public:
			UniformGenerator(double dzm);

			double getDZM();
			void setDZM(double val);

			path createDirectory(path outDir) override;

			vector<double> generateDZ() override;
	};
}

#endif