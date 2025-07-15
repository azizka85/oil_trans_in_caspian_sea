#ifndef WIND_INDUCED_CURRENTS_DAVIES85_VARIABLE_STEP_SIZE_SOLVER_H
#define WIND_INDUCED_CURRENTS_DAVIES85_VARIABLE_STEP_SIZE_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace WindInducedCurrents::Davies85::VariableStepSize {
    class Solver {
        private:
            double f;
            double b;

            double h;
            double w;
            double l;
            
            double g;
            double rho;
            double nu;
            
            double kb;
            double qx;
            double qy;

            double dx;
            double dy;
            vector<double> dz;

            double endTime;
            double outputTimeStep;
            string dir;
            
            double calculateDZFactor(double z);

            path createDirectory(double dzm);

            void setInitialCondition(
                int nx, int ny, int nz, 
                vector<vector<vector<double>>>& uf, 
                vector<vector<vector<double>>>& vf, 
                vector<vector<double>>& ua,
                vector<vector<double>>& va,
                vector<vector<double>>& z
            );                   
            
            void calculate_residual_elements(
                double dt, int nx, int ny, int nz,    
                vector<vector<double>>& ua,
                vector<vector<double>>& u1a,
                vector<vector<double>>& va,    
                vector<vector<double>>& v1a,
                vector<vector<vector<double>>>& uf, 
                vector<vector<vector<double>>>& vf,
                vector<vector<vector<double>>>& ud, 
                vector<vector<vector<double>>>& vd
            );

            double adjustTimeStep(double t, double dt, double dtMax);
            
            double maxAbsDifference(
                int nx, int ny, 
                vector<vector<double>>& u,
                vector<vector<double>>& u1
            );

            double maxAbsDifference(
                int nx, int ny, int nz, 
                vector<vector<vector<double>>> &u, 
                vector<vector<vector<double>>> &u1
            );

            void updateData(
                int nx, int ny, 
                vector<vector<double>>& u,
                vector<vector<double>>& u1
            );

            void updateData(
                int nx, int ny, int nz,
                vector<vector<vector<double>>>& u, 
                vector<vector<vector<double>>>& u1
            );

            void writeSurfaceData(
                double t, int m, 
                int nx, int ny, 
                vector<vector<double>> &ua, 
                vector<vector<double>> &va, 
                vector<vector<double>> &z, 
                path outDir
            );

            void writeVolumeData(
                double t, int m, 
                int nx, int ny, int nz, 
                vector<vector<double>> &ua, 
                vector<vector<double>> &va, 
                vector<vector<vector<double>>> &uf, 
                vector<vector<vector<double>>> &vf, 
                path outDir
            );

            void writeData(
                double t, int m,
                int nx, int ny, int nz,
                vector<vector<double>>& ua,
                vector<vector<double>>& va,
                vector<vector<double>>& z,
                vector<vector<vector<double>>>& uf, 
                vector<vector<vector<double>>>& vf,
                path outDir
            );

            void writeStatistics(
                vector<tuple<int, double, double, double, double>>& statistics, 
                path outDir
            );

        public:
            Solver(
                double f, double b,
                double h, double w, double l,
                double g, double rho, double nu,
                double kb, double qx, double qy,
                double dx, double dy, double dz,
                double endTime, double outputTimeStep, string dir
            );

            double getF();
            void setF(double val);

            double getB();
            void setB(double val);

            double getH();
            void setH(double val);

            double getW();
            void setW(double val);

            double getL();
            void setL(double val);

            double getG();
            void setG(double val);

            double getRHO();
            void setRHO(double val);

            double getNU();
            void setNU(double val);

            double getKB();
            void setKB(double val);

            double getQX();
            void setQX(double val);

            double getQY();
            void setQY(double val);

            double getDX();
            void setDX(double val);

            double getDY();
            void setDY(double val);

            vector<double> getDZ();
            void generateDZ(double dzMax);

            double getEndTime();
            void setEndTime(double val);

            double getOutputTimeStep();
            void setOutputTimeStep(double val);

            string getDir();                                                                        
            void setDir(string val);

            void solve();
    };
}

#endif
