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

#ifndef BINARYSTREAM_HPP_
#define BINARYSTREAM_HPP_

#include <string>
#include <vector>
#include <cassert>
#include <cstdint>
#include <memory>

class BinaryReadStream;

/// Standard size: 256 bytes
const size_t STD_BUFFER_SIZE = 256;

class BinaryWriteStream
{
    friend class BinaryReadStream;
public:
    /// @param size: Standard buffer capacity (gets increased if not sufficient).
    explicit BinaryWriteStream(size_t size = STD_BUFFER_SIZE);
    ~BinaryWriteStream();
    /// @return Current size of the used buffer (not the capacity)
    inline size_t getSize() const { return bufferSize; }
    inline const uint8_t *getBuffer() const { return buffer; }
    /// Manually make sure buffer holds at least passed size bytes.
    void reserve(size_t size = STD_BUFFER_SIZE);

    /// Write "size"-bytes of array "data"
    void write(const void *data, size_t size);
    /// Write a typed primitive value to the file (i.e. no pointers in type T).
    template<typename T>
    void write(const T &val) { write((const void*)&val, sizeof(T)); }
    /// Write a string to the file
    void write(const char *str);
    void write(const std::string &str);

    /// Write an array of primitive values to the file (i.e. no pointers in type T).
    template<typename T>
    void writeArray(const std::vector<T> &v)
    {
        uint32_t size = v.size();
        write(size);
        if (size > 0) {
            write((const void*)&v.front(), sizeof(T)*v.size());
        }
    }

    /// Serialization with pipe operator
    template<typename T>
    BinaryWriteStream& operator<<(const T &val) { write(val); return *this; }
    BinaryWriteStream& operator<<(const char *str) { write(str); return *this; }
    BinaryWriteStream& operator<<(const std::string &str) { write(str); return *this; }

protected:
    /// The current buffer size (only the used part of the buffer counts)
    size_t bufferSize;
    /// The maximum buffer size before it needs to be increased/reallocated
    size_t capacity;
    uint8_t *buffer;
};

class BinaryReadStream
{
public:
    /// Read from passed input stream.
    explicit BinaryReadStream(BinaryWriteStream &stream);
    /// Read from passed input buffer.
    BinaryReadStream(void *_buffer, size_t _bufferSize);
    BinaryReadStream(const void *_buffer, size_t _bufferSize);
    ~BinaryReadStream();
    inline size_t getSize() const { return bufferSize; }

    /// Deserialization (see BinaryWriteStream for details).
    void read(void *data, size_t size);
    template<typename T>
    void read(T &val) { read((void*)&val, sizeof(T)); }
    void read(std::string &str);

    template<typename T>
    void readArray(std::vector<T> &v)
    {
        uint32_t size;
        read(size);
        if (size > 0) {
            v.resize(size);
            read((void*)&v.front(), sizeof(T)*v.size());
        }
    }

    /// Deserialization with pipe operator
    template<typename T>
    BinaryReadStream& operator>>(T &val) { read(val); return *this; }
    BinaryReadStream& operator>>(std::string &str) { read(str); return *this; }

protected:
    /// The total buffer size
    size_t bufferSize;
    /// The current point in the buffer where the code reads from
    size_t bufferStart;
    uint8_t *buffer;
};

/*! BINARYSTREAM_HPP_ */
#endif
