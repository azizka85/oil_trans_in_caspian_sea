#include <utils/fs.h>

#include "volume.h"

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Writers;

void Volume::write(
    double t, int m,
    int nx, int ny, int nz,
    double dx, double dy,
    vector<double>& dz, vector<double>& h,
    vector<double>& ua, vector<double>& va,
    vector<double>& uf, vector<double>& vf,
    path outDir
) {
    auto file = FS::createFileByPath(
        outDir /
        path(
            format("data.{:03}.vtk", m)
        )
    );

    file << "# vtk DataFile Version 3.0" << endl;
    file << format("TIME {:.3f}", t) << endl;
    file << "ASCII" << endl;
    file << "DATASET STRUCTURED_GRID" << endl;
    file << format("DIMENSIONS {} {} {}", nx, ny, nz) << endl;
    file << format("POINTS {} double", nx * ny * nz) << endl;

    vector<vector<double>> z(nx, vector<double>(ny, 0));

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                double x = dx * i;
                double y = dy * j;

                int p = j + i * ny;

                z[i][j] += h[p] * dz[k] / 2;

                file << format("{:.3f} {:.3f} {:.3f}", x, y, z[i][j]) << endl;

                z[i][j] += h[p] * dz[k] / 2;
            }
        }
    }

    file << "FIELD FieldData 1" << endl;
    file << "Time 1 1 double" << endl;
    file << format("{:.3f}", t) << endl;
    file << format("POINT_DATA {}", nx * ny * nz) << endl;

    file << "VECTORS V double" << endl;

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int p = j + i * ny;
                int id = k + p * nz;

                file << format("{:.3f} {:.3f} 0", uf[id] + ua[p], vf[id] + va[p]) << endl;
            }
        }
    }
}

void Volume::writeViscosity(
    double t, int m, 
    int nx, int ny, int nz, 
    double dx, double dy, 
    vector<double>& dz, vector<double>& h, 
    vector<double>& nu, 
    path outDir
) {
    auto file = FS::createFileByPath(
        outDir /
        path(
            format("data.{:03}.vtk", m)
        )
    );

    file << "# vtk DataFile Version 3.0" << endl;
    file << format("TIME {:.3f}", t) << endl;
    file << "ASCII" << endl;
    file << "DATASET STRUCTURED_GRID" << endl;
    file << format("DIMENSIONS {} {} {}", nx, ny, nz) << endl;
    file << format("POINTS {} double", nx * ny * nz) << endl;

    vector<vector<double>> z(nx, vector<double>(ny, 0));

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                double x = dx * i;
                double y = dy * j;

                file << format("{:.3f} {:.3f} {:.3f}", x, y, z[i][j]) << endl;

                int p = j + i * ny;

                if (k < nz - 1) {
                    z[i][j] += h[p] * dz[k];
                }                                
            }
        }
    }

    file << "FIELD FieldData 1" << endl;
    file << "Time 1 1 double" << endl;
    file << format("{:.3f}", t) << endl;
    file << format("POINT_DATA {}", nx * ny * nz) << endl;

    file << "SCALARS nu double" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int p = j + i * ny;
                int id = k + p * nz;

                file << format("{:.3f}", nu[id]) << endl;
            }
        }
    }
}
