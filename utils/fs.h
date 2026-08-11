#ifndef UTILS_FS_H
#define UTILS_FS_H

#include <fstream>

#include <filesystem>

using namespace std;
using namespace std::filesystem;

namespace Utils::FS {
	path createDirByPath(path dirPath);
	ofstream createFileByPath(path filePath);
}

#endif