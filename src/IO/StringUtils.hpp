/*
 * BSD 2-Clause License
 *
 * Copyright (c) 2021, Christoph Neuhauser
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

#ifndef CFD3D_STRINGUTILS_HPP
#define CFD3D_STRINGUTILS_HPP

#include <string>

/**
 * Returns whether @str starts with @prefix.
 * @param str The full string.
 * @param prefix The prefix.
 * @return Returns whether the full string starts with the passed prefix.
 */
bool startsWith(const std::string& str, const std::string& prefix);

/**
 * Returns whether @str ends with @postfix.
 * @param str The full string.
 * @param postfix The postfix.
 * @return Returns whether the full string starts with the passed postfix.
 */
bool endsWith(const std::string& str, const std::string& postfix);

/**
 * Converts strings like "This is a test!" with separator ' ' to { "This", "is", "a", "test!" }.
 * @tparam InputIterator The list class to use.
 * @param stringObject The string to split.
 * @param separator The separator to use for splitting.
 * @param listObject The split parts.
 */
template<class InputIterator>
void splitString(const std::string &stringObject, char separator, InputIterator& listObject) {
    std::string buffer;
    for (char c : stringObject) {
        if (c != separator) {
            buffer += c;
        } else {
            if (buffer.length() > 0) {
                listObject.push_back(buffer);
                buffer = "";
            }
        }
    }
    if (buffer.length() > 0) {
        listObject.push_back(buffer);
        buffer = "";
    }
}

/**
 * Converts strings like "This is a test!" with separators ' ' and '\t' to { "This", "is", "a", "test!" }.
 * @tparam InputIterator The list class to use.
 * @param stringObject The string to split.
 * @param listObject The split parts.
 * NOTE: This is equivalent to
 * 'boost::algorithm::split(listObject, stringObject, boost::is_any_of("\t "), boost::token_compress_on);'.
 */
template<class InputIterator>
void splitStringWhitespace(const std::string &stringObject, InputIterator& listObject) {
    std::string buffer;
    for (char c : stringObject) {
        if (c != ' ' && c != '\t') {
            buffer += c;
        } else {
            if (buffer.length() > 0) {
                listObject.push_back(buffer);
                buffer = "";
            }
        }
    }
    if (buffer.length() > 0) {
        listObject.push_back(buffer);
        buffer = "";
    }
}

#endif //CFD3D_STRINGUTILS_HPP
