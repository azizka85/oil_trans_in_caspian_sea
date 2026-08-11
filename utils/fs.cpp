#include <stdexcept>

#include "fs.h"

using namespace Utils;

path FS::createDirByPath(path dirPath) {
    error_code ec;

    create_directories(dirPath, ec);

    if (ec || !filesystem::exists(dirPath)) {
        throw runtime_error(
            format("Could not create directory {}, message: {}", dirPath.string(), ec.message())
        );
    }

    return dirPath;
}

ofstream FS::createFileByPath(path filePath) {
    ofstream file(filePath);

    if (file.bad()) {
        throw runtime_error(
            format("Failed to open file at: {}", filePath.string())
        );
    }

    return file;
}