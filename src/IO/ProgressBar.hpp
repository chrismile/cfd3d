//
// Created by christoph on 01.06.19.
//

#ifndef CFD3D_PROGRESSBAR_HPP
#define CFD3D_PROGRESSBAR_HPP

#include <cstdint>
#include <string>
#include "Defines.hpp"

class ProgressBar {
public:
    /**
     * Prints a progress bar in the terminal. Therefore it overwrites the current line on stdout.
     * @param t The current time of the simulation.
     * @param tEnd The end time of the simulation.
     * @param max The maximal amount of # to show.
     */
    void printProgress(Real t, Real tEnd, size_t max);

    /**
     * Prints a message that a file was written at step n and time t.
     * @param n The current step of the simulation.
     * @param t The current time of the simulation.
     * @param max The amount of # used in progress bar to overwrite it.
     */
    void printOutput(int n, Real t, size_t max);
};


#endif //CFD3D_PROGRESSBAR_HPP
