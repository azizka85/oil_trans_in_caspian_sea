#ifndef WIND_INDUCED_CURRENTS_DAVIES85_VARIABLE_PARAMETERS_SOLVER_H
#define WIND_INDUCED_CURRENTS_DAVIES85_VARIABLE_PARAMETERS_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace WindInducedCurrents::Davies85::VariableParameters {
    enum GenerateH {
        UniformH = 0,
        CosineH = 1
    };

    enum GenerateNU {
        UniformNU = 0,
        FromWindSpeedNU = 1
    };

    enum GenerateQ {
        UniformQ = 0,
        FromWindSpeedQ = 1
    };

    enum GenerateDZ {
        UniformDZ = 0,
        ParabolicDZ = 1
    };

    class Solver {
        private:
            double f;
            double b;

            double hm;
            GenerateH hg;
            vector<vector<double>> h;
            
            double w;
            double l;
            
            double g;
            double rho;
            double kb;

            double num;
            GenerateNU nug;
            vector<vector<vector<double>>> nu;

            double qxm;
            double qym;
            GenerateQ qg;
            vector<vector<double>> qx;
            vector<vector<double>> qy;

            double dx;
            double dy;

            double dzm;
            GenerateDZ dzg;
            vector<double> dz;

            double endTime;
            double outputTimeStep;
            string dir;

            void generateH(int nx, int ny);
            void generateUniformH(int nx, int ny);
            void generateCosineH(int nx, int ny);

            tuple<double, double> minMaxH(int nx, int ny);

            void generateNU(
                int nx, int ny, int nz,
                vector<vector<double>>& ua,
                vector<vector<double>>& va,
                vector<vector<double>>& z,
                vector<vector<vector<double>>>& uf, 
                vector<vector<vector<double>>>& vf
            );
            
            void updateNU(
                double t, int n,
                int nx, int ny, int nz,
                vector<vector<double>>& ua,
                vector<vector<double>>& va,
                vector<vector<double>>& z,
                vector<vector<vector<double>>>& uf, 
                vector<vector<vector<double>>>& vf
            );

            void generateUniformNU(int nx, int ny, int nz);
            void updateUniformNU(int nx, int ny, int nz);

            void generateNUFromWindSpeed(int nx, int ny, int nz);
            void updateNUFromWindSpeed(int nx, int ny, int nz);

            double maxNU(int nx, int ny, int nz);

            void generateQX(
                int nx, int ny, int nz,
                vector<vector<double>>& ua,
                vector<vector<double>>& va,
                vector<vector<double>>& z,
                vector<vector<vector<double>>>& uf, 
                vector<vector<vector<double>>>& vf
            );

            void updateQX(
                double t, int n,
                int nx, int ny, int nz,
                vector<vector<double>>& ua,
                vector<vector<double>>& va,
                vector<vector<double>>& z,
                vector<vector<vector<double>>>& uf, 
                vector<vector<vector<double>>>& vf
            );

            void generateUniformQX(int nx, int ny);
            void updateUniformQX(int nx, int ny);

            void generateQXFromWindSpeed(int nx, int ny);
            void updateQXFromWindSpeed(int nx, int ny);

            void generateQY(                
                int nx, int ny, int nz,
                vector<vector<double>>& ua,
                vector<vector<double>>& va,
                vector<vector<double>>& z,
                vector<vector<vector<double>>>& uf, 
                vector<vector<vector<double>>>& vf
            );

            void updateQY(
                double t, int n,
                int nx, int ny, int nz,
                vector<vector<double>>& ua,
                vector<vector<double>>& va,
                vector<vector<double>>& z,
                vector<vector<vector<double>>>& uf, 
                vector<vector<vector<double>>>& vf
            );

            void generateUniformQY(int nx, int ny);
            void updateUniformQY(int nx, int ny);

            void generateQYFromWindSpeed(int nx, int ny);
            void updateQYFromWindSpeed(int nx, int ny);
            
            void generateDZ();
            void generateUniformDZ();
            void generateParabolicDZ();
            double calcParabolicDZFactor(double z);            

            path createDirectory();

            void setInitialCondition(
                int nx, int ny, int nz, 
                vector<vector<vector<double>>>& uf, 
                vector<vector<vector<double>>>& vf, 
                vector<vector<double>>& ua,
                vector<vector<double>>& va,
                vector<vector<double>>& z
            );                   
            
            void calculateResidualElements(
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

            void createTridiagonalMatrix(
                int i, int j, int nz,
                double dt,
                vector<double> &al,
                vector<double> &ac,
                vector<double> &ar                
            );

            double adjustTimeStep(double t, double dt, double tMax, double dtMax, bool mult);
            
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

            void writeHeights(int nx, int ny, path outDir);

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

            void writeViscosity(
                double t, int m, 
                int nx, int ny, int nz,
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
                double hm, GenerateH hg,
                double w, double l,
                double g, double rho, double kb, 
                double num, GenerateNU nug,
                double qxm, double qym, GenerateQ qg,
                double dx, double dy, 
                double dzm, GenerateDZ dzg,
                double endTime, double outputTimeStep, string dir
            );

            double getF();
            void setF(double val);

            double getB();
            void setB(double val);

            double getHM();
            void setHM(double val);

            GenerateH getHG();
            void setHG(GenerateH val);

            vector<vector<double>> getH();                        

            double getW();
            void setW(double val);

            double getL();
            void setL(double val);

            double getG();
            void setG(double val);

            double getRHO();
            void setRHO(double val);

            double getKB();
            void setKB(double val);

            double getNUM();
            void setNUM(double val);

            GenerateNU getNUG();
            void setNUG(GenerateNU val);

            vector<vector<vector<double>>> getNU();                        

            double getQXM();
            void setQXM(double val);

            double getQYM();
            void setQYM(double val);

            GenerateQ getQG();
            void setQG(GenerateQ val);

            vector<vector<double>> getQX();            
            vector<vector<double>> getQY();            

            double getDX();
            void setDX(double val);

            double getDY();
            void setDY(double val);

            double getDZM();
            void setDZM(double val);

            GenerateDZ getDZG();
            void setDZG(GenerateDZ dzg);

            vector<double> getDZ();            

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
