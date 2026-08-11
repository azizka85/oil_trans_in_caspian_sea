#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_WIND_UNIFORM_GENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_WIND_UNIFORM_GENERATOR_H

#include "igenerator.h"

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Wind {
	class UniformGenerator : public IGenerator {
	private:
		float u10m;
		float v10m;

		float qxm;
		float qym;

	public:
		UniformGenerator(
			float u10m, float v10m, 
			float qxm, float qym
		);

		float getU10M();
		void setU10M(float val);

		float getV10M();
		void setV10M(float val);

		float getQXM();
		void setQXM(float val);

		float getQYM();
		void setQYM(float val);

		path createDirectory(path outDir) override;

		vector<Data> generate(int nx, int ny) override;
	};
}

#endif