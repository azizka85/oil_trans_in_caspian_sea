#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_DZ_TRIPLE_POINT_GENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_DZ_TRIPLE_POINT_GENERATOR_H

#include "igenerator.h"

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::DZ {
	class TriplePointGenerator : public IGenerator {
		private:
			float dzMax;
			float dzMin;
			float zm;

		public:
			TriplePointGenerator(float dzMin, float dzMax, float zm);

			float getDZMin();
			float getDZMax();
			void setDZMinMax(float dzMin, float dzMax);

			float getZM();
			void setZM(float val);

			path createDirectory(path outDir) override;

			vector<float> generateDZ() override;

			float calcDZ(float z);
	};
}

#endif