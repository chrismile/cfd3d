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

#include <iostream>
#include <cassert>
#include <fstream>
#include "ClInterface.hpp"

ClInterface::~ClInterface()
{
}

void ClInterface::initialize(CLContextInfo contextInfo)
{
	// 1. Get right OpenCL platform
    cl::Platform::get(&allPlatforms);
    if (allPlatforms.size() <= contextInfo.platformNum) {
        std::cerr << "Fatal Error: No OpenCL platform with specified ID found!" << std::endl;
    	exit(1);
    }
    platform = allPlatforms.at(contextInfo.platformNum);

	// 2. Get right OpenCL device(s)
	platform.getDevices(CL_DEVICE_TYPE_ALL, &allDevices);
	if(allDevices.size() == 0){
        std::cout << "No OpenCL devices found. Please check your installation!" << std::endl;
		exit(1);
	}
	if (contextInfo.useAllDevices) {
		devices = allDevices;
	} else {
		devices = {allDevices.at(0)};
	}


	// 3. Create OpenCL context
	context = cl::Context(devices);
}

void ClInterface::printInfo()
{
    std::cout << "Number of OpenCL platforms: " << allPlatforms.size() << std::endl;
    std::cout << "Standard OpenCL platform:" << std::endl;
    std::cout << platform.getInfo<CL_PLATFORM_NAME>() << std::endl;
    std::cout << platform.getInfo<CL_PLATFORM_VERSION>() << std::endl << std::endl;
	int i = 0;
	for (cl::Device &device : devices) {
        std::cout << "Using device #" << i << ": " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
		i++;
	}
}


std::string loadTextFile(const std::string &filename) {
	std::ifstream ifile(filename.c_str());
	if (!ifile.is_open()) {
		std::cerr << "ERROR in loadTextFile: Couldn't open file \"" << filename <<  "\"." << std::endl;
		exit(1);
	}
	std::string kernelCode((std::istreambuf_iterator<char>(ifile)), std::istreambuf_iterator<char>());
	ifile.close();
	return kernelCode;
}

void ClInterface::setPrependHeader(const std::string &prependHeader)
{
    this->prependHeader = prependHeader;
}

cl::Program ClInterface::loadProgramFromSourceFile(const char *filename)
{
	return loadProgramFromSourceFiles({filename});
}

cl::Program ClInterface::compileSourceFile(const std::string &filename) {
	std::cout << "Compiling " << filename << "..." << std::endl;

	cl::Program::Sources sources;
	std::string kernelCode = prependHeader + "\n" + loadTextFile(filename.c_str());
	sources.push_back({kernelCode.c_str(), kernelCode.length()});
	cl::Program program(context, sources);
	std::string compileArgs = "-Werror";
	if (strcmp(platform.getInfo<CL_PLATFORM_NAME>().c_str(), "AMD Accelerated Parallel Processing") == 0) {
		// Odd bug on AMD APP SDK with "-Werror"
		compileArgs = "";
	}
	if (program.compile(compileArgs.c_str()) != CL_SUCCESS) {
        std::cerr << "Error while building " << filename << ":" << std::endl
                << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]) << std::endl;
		exit(1);
	}
	return program;
}

cl::Program ClInterface::loadProgramFromSourceFiles(const std::vector<std::string> &filenames)
{
	// 4. Load and compile OpenCL C source files
	std::vector<cl::Program> programs;
	for (const std::string& filename : filenames) {
		programs.push_back(compileSourceFile(filename));
	}

	// 5. Link all source files together into one executable
	cl_int err = CL_SUCCESS;
	cl::Program linkedProgram = cl::linkProgram(programs, NULL, NULL, NULL, &err);
	if (err != CL_SUCCESS) {
		std::cerr << "Error while linking programs:" << std::endl
		        << linkedProgram.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]) << std::endl;
		exit(1);
	}
	return linkedProgram;
}

cl::Program ClInterface::loadProgramFromBinaryFile(const char *filename)
{
	cl::Program::Binaries binaries;
	std::ifstream ifile(filename, std::ifstream::binary);
	ifile.seekg(0, std::ios::end);
	std::streampos size = ifile.tellg();
	ifile.seekg(0, std::ios::beg);
	std::vector<char> buffer(size);
	ifile.read(buffer.data(), size);
	ifile.close();
	binaries.push_back({std::make_pair(buffer.data(), size)});
	
	cl::Program program(context, devices, binaries);
	if (program.build(devices) != CL_SUCCESS){
        std::cerr << "Error while builing cl::Program: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0])
                << std::endl;
		exit(1);
	}
	return program;
}



cl::NDRange ClInterface::rangePadding1D(int w, int local) {
	return cl::NDRange(iceil(w, local)*local);
}

cl::NDRange ClInterface::rangePadding2D(int w, int h, cl::NDRange local) {
    return cl::NDRange(iceil(w, local[0])*local[0], iceil(h, local[1])*local[1]);
}

cl::NDRange ClInterface::rangePadding3D(int w, int h, int d, cl::NDRange local) {
    return cl::NDRange(iceil(w, local[0])*local[0], iceil(h, local[1])*local[1], iceil(d, local[2])*local[2]);
}

