#ifndef WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_VARIABLE_PARAMETERS_SOLVER_H
#define WIND_INDUCED_CURRENTS_DAVIES85_PARALLEL_VARIABLE_PARAMETERS_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

#include <CL/opencl.hpp>

using namespace std;

using namespace std::filesystem;

namespace WindInducedCurrents::Davies85::Parallel::VariableParameters {
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
            float f;
            float b;

            float hm;
            GenerateH hg;
            vector<float> h;
            
            float w;
            float l;
            
            float g;
            float rho;
            float kb;

            float num;
            GenerateNU nug;
            vector<float> nu;

            float qxm;
            float qym;
            GenerateQ qg;
            vector<float> qx;
            vector<float> qy;

            float dx;
            float dy;

            float dzm;
            GenerateDZ dzg;
            vector<float> dz;

            float endTime;
            float outputTimeStep;
            string dir;

            void generateH(int nx, int ny);
            void generateUniformH(int nx, int ny);
            void generateCosineH(int nx, int ny);

            tuple<float, float> minMaxH(int nx, int ny);

            void generateNU(
                int nx, int ny, int nz,
                vector<float>& ua,
                vector<float>& va,
                vector<float>& z,
                vector<float>& uf, 
                vector<float>& vf
            );

            void generateUniformNU(int nx, int ny, int nz);
            void generateNUFromWindSpeed(int nx, int ny, int nz);

            cl::Kernel createUpdateNUKernel(                
                int ny, int nz,
                cl::Buffer& bufferNU,
                cl::Program& program
            );

            cl::Kernel createUpdateUniformNUKernel(
                int ny, int nz,
                cl::Buffer& bufferNU,
                cl::Program& program
            );

            cl::Kernel createUpdateNUFromWindSpeed(
                int ny, int nz,
                cl::Buffer& bufferNU,
                cl::Program& program
            );

            float maxNU(int nx, int ny, int nz);

            void generateQX(
                int nx, int ny, int nz,
                vector<float>& ua,
                vector<float>& va,
                vector<float>& z,
                vector<float>& uf, 
                vector<float>& vf
            );

            void generateUniformQX(int nx, int ny);
            void generateQXFromWindSpeed(int nx, int ny);

            void generateQY(                
                int nx, int ny, int nz,
                vector<float>& ua,
                vector<float>& va,
                vector<float>& z,
                vector<float>& uf, 
                vector<float>& vf
            );

            void generateUniformQY(int nx, int ny);
            void generateQYFromWindSpeed(int nx, int ny);

            cl::Kernel createUpdateQKernel(float qm, int ny, cl::Buffer& bufferQ, cl::Program& program);
            cl::Kernel createUpdateUniformQKernel(float qm, int ny, cl::Buffer& bufferQ, cl::Program& program);
            cl::Kernel createUpdateQFromWindSpeedKernel(float qm, int ny, cl::Buffer& bufferQ, cl::Program& program);

            void generateDZ();
            void generateUniformDZ();
            void generateParabolicDZ();
            float calcParabolicDZFactor(float z);            

            path createDirectory();

            string loadKernelSource(path filePath);
            
            void loadKernelSources(cl::Program::Sources &sources);            

            void setInitialCondition(
                int nx, int ny, int nz, 
                vector<float>& uf, 
                vector<float>& vf, 
                vector<float>& ua,
                vector<float>& va,
                vector<float>& z
            );                                       

            float adjustTimeStep(float t, float dt, float tMax, float dtMax, bool mult);
            
            float maxAbsDifference(
                int nx, int ny, 
                vector<float>& u,
                vector<float>& u1
            );

            float maxAbsDifference(
                int nx, int ny, int nz, 
                vector<float> &u, 
                vector<float> &u1
            );

            void updateData(
                int nx, int ny, 
                vector<float>& u,
                vector<float>& u1
            );

            void updateData(
                int nx, int ny, int nz,
                vector<float>& u, 
                vector<float>& u1
            );

            void writeHeights(int nx, int ny, path outDir);

            void writeSurfaceData(
                float t, int m, 
                int nx, int ny, 
                vector<float> &ua, 
                vector<float> &va, 
                vector<float> &z,
                path outDir
            );

            void writeVolumeData(
                float t, int m, 
                int nx, int ny, int nz, 
                vector<float> &ua, 
                vector<float> &va, 
                vector<float> &uf, 
                vector<float> &vf, 
                path outDir
            );

            void writeViscosity(
                float t, int m, 
                int nx, int ny, int nz,
                path outDir
            );

            void writeData(
                float t, int m,
                int nx, int ny, int nz,
                vector<float>& ua,
                vector<float>& va,
                vector<float>& z,
                vector<float>& uf, 
                vector<float>& vf,
                path outDir
            );

            void writeStatistics(
                vector<tuple<int, float, long long, float, float, float>>& statistics, 
                path outDir
            );

        public:
            Solver(
                float f, float b,
                float hm, GenerateH hg,
                float w, float l,
                float g, float rho, float kb, 
                float num, GenerateNU nug,
                float qxm, float qym, GenerateQ qg,
                float dx, float dy, 
                float dzm, GenerateDZ dzg,
                float endTime, float outputTimeStep, string dir
            );

            float getF();
            void setF(float val);

            float getB();
            void setB(float val);

            float getHM();
            void setHM(float val);

            GenerateH getHG();
            void setHG(GenerateH val);

            vector<float> getH();                        

            float getW();
            void setW(float val);

            float getL();
            void setL(float val);

            float getG();
            void setG(float val);

            float getRHO();
            void setRHO(float val);

            float getKB();
            void setKB(float val);

            float getNUM();
            void setNUM(float val);

            GenerateNU getNUG();
            void setNUG(GenerateNU val);

            vector<float> getNU();                        

            float getQXM();
            void setQXM(float val);

            float getQYM();
            void setQYM(float val);

            GenerateQ getQG();
            void setQG(GenerateQ val);

            vector<float> getQX();            
            vector<float> getQY();            

            float getDX();
            void setDX(float val);

            float getDY();
            void setDY(float val);

            float getDZM();
            void setDZM(float val);

            GenerateDZ getDZG();
            void setDZG(GenerateDZ dzg);

            vector<float> getDZ();            

            float getEndTime();
            void setEndTime(float val);

            float getOutputTimeStep();
            void setOutputTimeStep(float val);

            string getDir();                                                                        
            void setDir(string val);

            void solve();
    };
}

#endif
