#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_WIND_UNIFORM_GENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_WIND_UNIFORM_GENERATOR_H

#include "igenerator.h"

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Wind {
	class UniformGenerator : public IGenerator {
	private:
		double u10m;
		double v10m;

		double qxm;
		double qym;

	public:
		UniformGenerator(
			double u10m, double v10m, 
			double qxm, double qym
		);

		double getU10M();
		void setU10M(double val);

		double getV10M();
		void setV10M(double val);

		double getQXM();
		void setQXM(double val);

		double getQYM();
		void setQYM(double val);

		path createDirectory(path outDir) override;

		vector<Data> generate(int nx, int ny) override;
	};
}

#endif