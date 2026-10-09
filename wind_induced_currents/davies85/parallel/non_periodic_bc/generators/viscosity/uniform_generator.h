#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_VISCOSITY_UNIFORM_GENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_VISCOSITY_UNIFORM_GENERATOR_H

#include "igenerator.h"

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Viscosity {
	class UniformGenerator : public IGenerator {
		private:
			double num;

		public:
			UniformGenerator(double num);

			double getNUM();
			void setNUM(double val);

			path createDirectory(path outDir) override;

			vector<double> generateNU(
				int nx, int ny, int nz,
				vector<double>& dz, vector<double>& h,
				vector<double>& u10, vector<double>& v10,
				vector<double>& qx, vector<double>& qy,
				vector<double>& ua, vector<double>& va
			) override;
	};
}

#endif