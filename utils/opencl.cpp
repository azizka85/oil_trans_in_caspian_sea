#include <fstream>

#include <stdexcept>

#include "opencl.h"

using namespace Utils;

string OpenCL::loadKernelSource(path filePath) {
    ifstream file(filePath);

    if (file.bad()) throw runtime_error("Cannot open file: " + filePath.string());

    return std::string(
        istreambuf_iterator<char>(file),
        istreambuf_iterator<char>()
    );
}