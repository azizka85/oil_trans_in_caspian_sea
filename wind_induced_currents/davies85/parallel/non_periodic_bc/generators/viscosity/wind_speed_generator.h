#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_VISCOSITY_WIND_SPEED_GENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_VISCOSITY_WIND_SPEED_GENERATOR_H

#include "igenerator.h"

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Viscosity {
	class WindSpeedGenerator : public IGenerator {
		private:
			double k0;
			double k;

			double f;
			double sigma;

			double rho;
			double nu0;

		public:
			WindSpeedGenerator(
				double k0, double k, 
				double f, double sigma,
				double rho, double nu0
			);

			double getK0();
			void setK0(double val);

			double getK();
			void setK(double val);

			double getF();
			void setF(double val);

			double getSigma();
			void setSigma(double val);

			double getRho();
			void setRho(double val);

			double getNU0();
			void setNU0(double val);

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