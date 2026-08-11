#include <utils/fs.h>

#include "height.h"

using namespace Utils;

using namespace WindInducedCurrents::Davies85::Parallel::NonPeriodicBC::Writers;

void Height::write(vector<float>& h, int nx, int ny, float dx, float dy, path outDir) {
    auto file = FS::createFileByPath(outDir / path("height.vtk"));

    file << "# vtk DataFile Version 3.0" << endl;
    file << format("TIME {:.3f}", 0.) << endl;
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
    file << format("{:.3f}", 0.) << endl;
    file << format("POINT_DATA {}", nx * ny) << endl;

    file << "SCALARS h float" << endl;
    file << "LOOKUP_TABLE default" << endl;

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int p = j + i * ny;

            file << format("{:.3f}", h[p]) << endl;
        }
    }
}