#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_SOLVER_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_NON_PERIODIC_BC_SOLVER_H

#include <string>

#include <memory>

#include <CL/opencl.hpp>

#include "generators/area/igenerator.h"
#include "generators/dz/igenerator.h"
#include "generators/bathymetry/igenerator.h"
#include "generators/wind/igenerator.h"
#include "generators/viscosity/igenerator.h"

using namespace std;

namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC {
    struct Directories {
        path root;
        path surface;
        path volume;
        path viscosity;
    };

    class Solver {
        private:
            double b;

            double f;
            double g;
            double rho;
            double kb;            

            double dx;
            double dy;

            double endTime;
            double outputTimeStep;
            string outDir;

            unique_ptr<Generators::Area::IGenerator> areaGenerator;
            unique_ptr<Generators::DZ::IGenerator> dzGenerator;           
            unique_ptr<Generators::Bathymetry::IGenerator> hGenerator;
            unique_ptr<Generators::Wind::IGenerator> qGenerator;
            unique_ptr<Generators::Viscosity::IGenerator> nuGenerator;

            Directories createDirectory();

            void loadKernelSources(cl::Program::Sources& sources);

            void setInitialCondition(
                int nx, int ny, int nz,
                vector<double>& uf,
                vector<double>& vf,
                vector<double>& ua,
                vector<double>& va,
                vector<double>& z
            );

            void writeData(
                double t, int m,
                int nx, int ny, int nz,
                vector<double>& dz, 
                vector<double>& h, vector<double>& z,
                vector<double>& ua, vector<double>& va,                 
                vector<double>& uf, vector<double>& vf,
                vector<double>& qx, vector<double>& qy,
                vector<double>& nu,
                Directories &dirs
            );

        public:
            Solver(
                double b,
                double f, double g, double rho, double kb,                
                double dx, double dy,
                double endTime, double outputTimeStep, string outDir,
                unique_ptr<Generators::Area::IGenerator> areaGenerator,
                unique_ptr<Generators::DZ::IGenerator> dzGenerator,                
                unique_ptr<Generators::Bathymetry::IGenerator> hGenerator,
                unique_ptr<Generators::Wind::IGenerator> qGenerator,
                unique_ptr<Generators::Viscosity::IGenerator> nuGenerator
            );

            double getB();
            void setB(double val);

            double getF();
            void setF(double val);

            double getG();
            void setG(double val);

            double getRHO();
            void setRHO(double val);

            double getKB();
            void setKB(double val);            

            double getDX();
            void setDX(double val);

            double getDY();
            void setDY(double val);

            double getEndTime();
            void setEndTime(double val);

            double getOutputTimeStep();
            void setOutputTimeStep(double val);

            string getOutDir();
            void setOutDir(string val);

            void setAreaGenerator(unique_ptr<Generators::Area::IGenerator> val);
            void setDZGenerator(unique_ptr<Generators::DZ::IGenerator> val);            
            void setHGenerator(unique_ptr<Generators::Bathymetry::IGenerator> val);
            void setQGenerator(unique_ptr<Generators::Wind::IGenerator> val);
            void setNUGenerator(unique_ptr<Generators::Viscosity::IGenerator> val);

            void solve();
    };
}

#endif