#include <utils/fs.h>

#include "surface.h"

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Writers;

void Surface::write(
    float t, int m,
    int nx, int ny,
    float dx, float dy,
    vector<float>& ua, vector<float>& va,
    vector<float>& qx, vector<float>& qy,
    vector<float>& z, path outDir
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
    file << format("DIMENSIONS {} {} 1", nx, ny) << endl;
    file << format("POINTS {} float", nx * ny) << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            float x = dx * i;
            float y = dy * j;

            file << format("{:.3f} {:.3f} 0.0", x, y) << endl;
        }
    }

    file << "FIELD FieldData 1" << endl;
    file << "Time 1 1 float" << endl;
    file << format("{:.3f}", t) << endl;
    file << format("POINT_DATA {}", nx * ny) << endl;

    file << "VECTORS VA float" << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int p = j + i * ny;

            file << format("{:.3f} {:.3f} 0", ua[p], va[p]) << endl;
        }
    }

    file << "VECTORS Q float" << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int p = j + i * ny;

            file << format("{:.3f} {:.3f} 0", qx[p], qy[p]) << endl;
        }
    }

    file << "SCALARS z float" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int p = j + i * ny;

            file << format("{:.7f}", z[p]) << endl;
        }
    }
}