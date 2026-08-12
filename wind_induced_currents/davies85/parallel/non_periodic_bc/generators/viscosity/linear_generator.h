#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_VISCOSITY_LINEAR_GENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_VISCOSITY_LINEAR_GENERATOR_H

#include "igenerator.h"

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Viscosity {
	class LinearGenerator : public IGenerator {
	private:
		float ht;

		float nus;
		float nut;

	public:
		LinearGenerator(
			float ht, 
			float nus, float nut
		);

		float getHT();
		void setHT(float val);

		float getNUS();
		void setNUS(float val);

		float getNUT();
		void setNUT(float val);

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