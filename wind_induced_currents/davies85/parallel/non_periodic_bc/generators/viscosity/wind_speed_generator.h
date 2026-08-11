#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_VISCOSITY_WIND_SPEED_GENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_VISCOSITY_WIND_SPEED_GENERATOR_H

#include "igenerator.h"

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Viscosity {
	class WindSpeedGenerator : public IGenerator {
		private:
			float k0;
			float k;

			float f;
			float sigma;

			float rho;
			float nu0;

		public:
			WindSpeedGenerator(
				float k0, float k, 
				float f, float sigma,
				float rho, float nu0
			);

			float getK0();
			void setK0(float val);

			float getK();
			void setK(float val);

			float getF();
			void setF(float val);

			float getSigma();
			void setSigma(float val);

			float getRho();
			void setRho(float val);

			float getNU0();
			void setNU0(float val);

			path createDirectory(path outDir) override;

			vector<float> generateNU(
				int nx, int ny, int nz,
				vector<float>& dz, vector<float>& h,
				vector<float>& u10, vector<float>& v10,
				vector<float>& qx, vector<float>& qy,
				vector<float>& ua, vector<float>& va
			) override;
	};
}

#endif