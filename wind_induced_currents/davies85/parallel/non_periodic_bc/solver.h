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
            float b;

            float f;
            float g;
            float rho;
            float kb;            

            float dx;
            float dy;

            float endTime;
            float outputTimeStep;
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
                vector<float>& uf,
                vector<float>& vf,
                vector<float>& ua,
                vector<float>& va,
                vector<float>& z
            );

            void writeData(
                float t, int m,
                int nx, int ny, int nz,
                vector<float>& dz, 
                vector<float>& h, vector<float>& z,
                vector<float>& ua, vector<float>& va,                 
                vector<float>& uf, vector<float>& vf,
                vector<float>& qx, vector<float>& qy,
                vector<float>& nu,
                Directories &dirs
            );

        public:
            Solver(
                float b,
                float f, float g, float rho, float kb,                
                float dx, float dy,
                float endTime, float outputTimeStep, string outDir,
                unique_ptr<Generators::Area::IGenerator> areaGenerator,
                unique_ptr<Generators::DZ::IGenerator> dzGenerator,                
                unique_ptr<Generators::Bathymetry::IGenerator> hGenerator,
                unique_ptr<Generators::Wind::IGenerator> qGenerator,
                unique_ptr<Generators::Viscosity::IGenerator> nuGenerator
            );

            float getB();
            void setB(float val);

            float getF();
            void setF(float val);

            float getG();
            void setG(float val);

            float getRHO();
            void setRHO(float val);

            float getKB();
            void setKB(float val);            

            float getDX();
            void setDX(float val);

            float getDY();
            void setDY(float val);

            float getEndTime();
            void setEndTime(float val);

            float getOutputTimeStep();
            void setOutputTimeStep(float val);

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