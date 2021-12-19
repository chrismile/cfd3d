/*
 * BSD 2-Clause License
 *
 * Copyright (c) 2019, Christoph Neuhauser
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

#include "BinaryStream.hpp"
#include <algorithm>
#include <iostream>
#include <cstring>
#include <cmath>

BinaryWriteStream::BinaryWriteStream(size_t size /* = STD_BUFFER_SIZE */)
{
    bufferSize = 0;
    capacity = 0;
    buffer = nullptr;
    reserve(size);
}

BinaryWriteStream::~BinaryWriteStream()
{
    if (buffer) {
        delete[] buffer;
        buffer = nullptr;
        capacity = 0;
        bufferSize = 0;
    }
}

void BinaryWriteStream::reserve(size_t size /* = STD_BUFFER_SIZE */)
{
    size = std::max((size_t)4, size); // Minimum buffer size: 32 bits
    if (size > capacity) {
        auto *_buffer = new uint8_t[size];
        if (buffer) {
            memcpy(_buffer, buffer, bufferSize);
            delete[] buffer;
        }
        buffer = _buffer;
        capacity = size;
    }
}

void BinaryWriteStream::write(const void *data, size_t size)
{
    // Check if we need to increase the buffer size
    if (bufferSize + size > capacity) {
        reserve(std::max(bufferSize + size, bufferSize*2));
    }

    assert(bufferSize + size <= capacity);
    memcpy(buffer + bufferSize, data, size);
    bufferSize += size;
}

void BinaryWriteStream::write(const char *str)
{
    uint32_t strSize = strlen(str);
    write(strSize);
    write((void*)str, strSize);
}

void BinaryWriteStream::write(const std::string &str)
{
    uint32_t strSize = str.size();
    write(strSize);
    write((void*)str.c_str(), strSize);
}



BinaryReadStream::BinaryReadStream(BinaryWriteStream &stream)
{
    // Copy the buffer address to this stream
    buffer = stream.buffer;
    bufferSize = stream.bufferSize;
    bufferStart = 0;

    // Delete the buffer from the old stream
    stream.buffer = nullptr;
    stream.bufferSize = 0;
    stream.capacity = 0;
}

BinaryReadStream::BinaryReadStream(void *_buffer, size_t _bufferSize)
{
    buffer = (uint8_t*)_buffer;
    bufferSize = _bufferSize;
    bufferStart = 0;
}

BinaryReadStream::BinaryReadStream(const void *_buffer, size_t _bufferSize)
{
    buffer = new uint8_t[_bufferSize];
    memcpy(buffer, _buffer, _bufferSize);
    bufferSize = _bufferSize;
    bufferStart = 0;
}

BinaryReadStream::~BinaryReadStream()
{
    if (buffer) {
        delete[] buffer;
        buffer = nullptr;
        bufferStart = 0;
        bufferSize = 0;
    }
}

void BinaryReadStream::read(void *data, size_t size)
{
    if (bufferStart + size > bufferSize) {
        std::cerr << "FATAL ERROR: BinaryReadStream::read(void*, size_t)" << std::endl;
        return;
    }
    memcpy(data, buffer + bufferStart, size);
    bufferStart += size;
}

void BinaryReadStream::read(std::string &str)
{
    uint32_t strSize = 0;
    read(strSize);
    if (bufferStart + (size_t)strSize > bufferSize) {
        std::cerr << "FATAL ERROR: BinaryReadStream::read(string&)" << std::endl;
        return;
    }
    char *cstr = new char[strSize + 1];
    cstr[strSize] = '\0';
    read((void*)cstr, strSize);
    str = cstr;
    delete[] cstr;
}
