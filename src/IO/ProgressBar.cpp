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
#include "ProgressBar.hpp"

const char prefix[] = "Progress: [";
const char suffix[] = "]";
const size_t prefixLength = sizeof(prefix) - 1;
const size_t suffixLength = sizeof(suffix) - 1;

void ProgressBar::printOutput(int n, Real t, std::size_t maxVal) {
    char output[] = "\rWrote file at %5s %6i %5s %10f\n";
    size_t lengthToOverwrite = maxVal + prefixLength + suffixLength + 30;
    std::string buffer(lengthToOverwrite, ' ');

    std::cout << "\r" << buffer;
    std::cout.flush();

    printf(output, "n: ", n, "t:", t);
    std::cout.flush();
}

void ProgressBar::printProgress(Real t, Real tEnd, std::size_t maxVal) {
    int currentProgressPercent = (int)((t / tEnd) * 100);
    if (lastProgressPercent == currentProgressPercent) {
        // Don't print anything if the percentage hasn't changed.
        return;
    }
    lastProgressPercent = currentProgressPercent;

    //const size_t suffix_length = sizeof(suffix) - 1;
    auto progressChars = size_t(t / tEnd * Real(maxVal));
    std::string buffer =
            prefix + std::string(progressChars, '#') + std::string(maxVal - progressChars, ' ') + suffix;

    std::cout << "\r" << char(27) << "[2K\r" << buffer;
    std::cout << currentProgressPercent << "%\tt = " << t;
    std::cout.flush();
}
