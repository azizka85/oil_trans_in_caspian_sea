#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_DZ_TRIPLE_POINT_GENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_DZ_TRIPLE_POINT_GENERATOR_H

#include "igenerator.h"

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::DZ {
	class TriplePointGenerator : public IGenerator {
		private:
			double dzMax;
			double dzMin;
			double zm;

		public:
			TriplePointGenerator(double dzMin, double dzMax, double zm);

			double getDZMin();
			double getDZMax();
			void setDZMinMax(double dzMin, double dzMax);

			double getZM();
			void setZM(double val);

			path createDirectory(path outDir) override;

			vector<double> generateDZ() override;

			double calcDZ(double z);
	};
}

#endif