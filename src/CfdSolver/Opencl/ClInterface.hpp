/*
 * BSD 2-Clause License
 *
 * Copyright (c) 2019, Christoph Neuhauser, Stefan Haas, Paul Ng
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef CLINTERFACE_HPP_
#define CLINTERFACE_HPP_

#define CL_TARGET_OPENCL_VERSION 120
//#define CL_HPP_ENABLE_EXCEPTIONS
//#define __CL_ENABLE_EXCEPTIONS

#include <vector>
#include <string>
#include "CL/cl.hpp"
#include "Singleton.hpp"

/// Represents information on OpenCL context creation
struct CLContextInfo {
    /// The number of the OpenCL platform to be used (usually 0)
    size_t platformNum;
    /**
     * false: Use only one (default) device
     * true: Use all devices available (usually multiple GPUs), experimental
     */
    bool useAllDevices;

    CLContextInfo(size_t platformNum = 0, bool useAllDevices = false)
        : platformNum(platformNum), useAllDevices(useAllDevices)
    {
    }
};

/**
 * A wrapper class for accessing initialization & program loading functionality of OpenCL in an easy way.
 */
class ClInterface : public Singleton<ClInterface> {
public:
    virtual ~ClInterface();

    /// Initializes an OpenCL context for the default platform & default device.
	void initialize(CLContextInfo contextInfo = CLContextInfo());
    /// Prints information about the installed OpenCL platform/devices on the command line
	void printInfo();

    /// Loads and compiles a compute program from OpenCL C source files.
	cl::Program loadProgramFromSourceFile(const char *filename);
	cl::Program loadProgramFromSourceFiles(const std::vector<std::string> &filenames);
    /// Loads a compute program from a pre-compiled OpenCL C binary file (device specific!).
	cl::Program loadProgramFromBinaryFile(const char *filename);

	/// Sets a header which should be prepended to all OpenCL C source files.
	void setPrependHeader(const std::string &prependHeader);

    /// Utility functions for enqueueing kernels
    inline int iceil(int x, int y) { return (x - 1) / y + 1; }
    cl::NDRange rangePadding1D(int w, int local);
    cl::NDRange rangePadding2D(int w, int h, cl::NDRange local);
    cl::NDRange rangePadding3D(int w, int h, int d, cl::NDRange local);

    /// Getter functions
    inline cl::Context &getContext() { return context; }
    inline cl::Platform &getPlatform() { return platform; }
    inline std::vector<cl::Device> &getDevices() { return devices; }

private:
    /**
     * Compiles a source file and returns a cl::Program object.
     * @param filename: The filename of the source file
     */
    cl::Program compileSourceFile(const std::string &filename);

    std::vector<cl::Platform> allPlatforms;
	std::vector<cl::Device> allDevices;
    cl::Platform platform; ///< Used platform
	std::vector<cl::Device> devices; ///< Used devices
    cl::Context context;
    std::string prependHeader;
};

std::string loadTextFile(const std::string &filename);

#endif /* CLINTERFACE_HPP_ */
