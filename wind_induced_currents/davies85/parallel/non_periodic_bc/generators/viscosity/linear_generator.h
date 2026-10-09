#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_VISCOSITY_LINEAR_GENERATOR_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_GENERATORS_VISCOSITY_LINEAR_GENERATOR_H

#include "igenerator.h"

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Generators::Viscosity {
	class LinearGenerator : public IGenerator {
	private:
		double ht;

		double nus;
		double nut;

	public:
		LinearGenerator(
			double ht, 
			double nus, double nut
		);

		double getHT();
		void setHT(double val);

		double getNUS();
		void setNUS(double val);

		double getNUT();
		void setNUT(double val);

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