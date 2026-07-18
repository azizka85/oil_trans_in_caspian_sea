#ifndef DIFFUSION_LAGRANGIAN_1D_SOLVER_H
#define DIFFUSION_LAGRANGIAN_1D_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

#include <random>

using namespace std;

using namespace std::filesystem;

namespace Diffusion::Lagrangian {
    class Solver {
        private:
            double l;
            double a;

            int Np;
            
            double D;

            double r;
            
            double dx;

            double endTime;
            double outputTimeStep;
            string dir;

            path createDirectory();

            tuple<int, vector<double>> setInitialCondition(int nx, vector<int>& C);
            
            int maxAbsDifference(int nx, vector<int>& C, vector<int>& C1);

            void updateData(int N, int nx, vector<double>& xp, vector<int>& C, vector<double>& xp1, vector<int>& C1);
            void writeData(vector<int>& C, double t, int nx, int m, path outDir);
            void writeStatistics(vector<tuple<int, double, int>>& statistics, path outDir);

        public:
            Solver(
                double l, double a,
				int Np,
                double D, double r,
                double dx,
                double endTime, double outputTimeStep, string dir
            );

            double getL();
            void setL(double val);

			double getA();
			void setA(double val);

			int getNp();
			void setNp(int val);

            double getD();
            void setD(double val);

            double getR();
            void setR(double val);

            double getDX();
            void setDX(double val);

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
