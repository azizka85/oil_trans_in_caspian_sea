#ifndef DIFFUSION_EULERIAN_2D_SOLVER_H
#define DIFFUSION_EULERIAN_2D_SOLVER_H

#include <tuple>
#include <string>
#include <vector>

#include <filesystem>

using namespace std;

using namespace std::filesystem;

namespace Diffusion::Eulerian {
    class Solver {
        private:
            double l;
			double h;

            double a;
            double b;
            
            double D;
            double r;
            
            double dx;
			double dy;

            double endTime;
            double outputTimeStep;
            string dir;

            path createDirectory();

            void setInitialCondition(int nx, int ny, vector<vector<double>>& C);
            
            double maxAbsDifference(int nx, int ny, vector<vector<double>>& C, vector<vector<double>>& C1);

            void updateData(int nx, int ny, vector<vector<double>>& C, vector<vector<double>>& C1);
            void writeData(vector<vector<double>>& C, double t, int nx, int ny, int m, path outDir);
            void writeStatistics(vector<tuple<int, double, double>>& statistics, path outDir);

        public:
            Solver(
                double l, double h,
				double a, double b,
                double D, double r,
				double dx, double dy,
                double endTime, double outputTimeStep, string dir
            );

            double getL();
            void setL(double val);

            double getH();
            void setH(double val);

			double getA();
			void setA(double val);

            double getB();
            void setB(double val);

            double getD();
            void setD(double val);

            double getR();
            void setR(double val);

            double getDX();
            void setDX(double val);

            double getDY();
            void setDY(double val);

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
