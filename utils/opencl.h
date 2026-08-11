#ifndef UTILS_OPENCL_H
#define UTILS_OPENCL_H

#include <string>

#include <filesystem>

using namespace std;
using namespace std::filesystem;

namespace Utils::OpenCL {
	string loadKernelSource(path filePath);
}

#endif