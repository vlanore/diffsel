/*Copyright or © or Copr. Centre National de la Recherche Scientifique (CNRS) (2017-06-14).
Contributors:
* Nicolas LARTILLOT - nicolas.lartillot@univ-lyon1.fr

This software is a computer program whose purpose is to detect convergent evolution using Bayesian
phylogenetic codon models.

This software is governed by the CeCILL-C license under French law and abiding by the rules of
distribution of free software. You can use, modify and/ or redistribute the software under the terms
of the CeCILL-C license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy, modify and redistribute
granted by the license, users are provided only with a limited warranty and the software's author,
the holder of the economic rights, and the successive licensors have only limited liability.

In this respect, the user's attention is drawn to the risks associated with loading, using,
modifying and/or developing or reproducing the software by the user in light of its specific status
of free software, that may mean that it is complicated to manipulate, and that also therefore means
that it is reserved for developers and experienced professionals having in-depth computer knowledge.
Users are therefore encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or data to be ensured and,
more generally, to use and operate it in the same conditions as regards security.

The fact that you are presently reading this means that you have had knowledge of the CeCILL-C
license and that you accept its terms.*/

#ifndef CHRONO_H
#define CHRONO_H

#include <sys/time.h>
#include <chrono>
#include <ctime>
#include <iostream>
#include <sstream>

class Chrono {
  public:
    void Reset();
    void Start();
    void Stop();

    inline int operator++() { return N++; }

    inline double GetTime() { return TotalTime; }

    inline double GetTimePerCount() { return TotalTime / N; }

    inline int GetCount() { return N; }

  private:
    // this is in milli seconds
    double sec1;
    double sec2;
    double milli1;
    double milli2;
    double TotalTime;
    int N;
};

class MeasureTime : public std::stringstream {
  public:
    MeasureTime() : stopped(false) { start(); }

    void start() {
        counter = std::chrono::high_resolution_clock::now();
        stopped = false;
    }

    void stop() {
        stopped = true;
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - counter);
    }

    template <int i>
    void print(std::string message) {
        if (!stopped) {
            stop();
        }
        std::string left(2 * i, ' ');
        std::cout << left << "* " << message << str() << "Time: " << duration.count() << "ms."
                  << std::endl;
        str("");
        start();
    }

    template <int i>
    void print() {
        print<i>("");
    }

  private:
    std::chrono::time_point<std::chrono::high_resolution_clock> counter;
    std::chrono::milliseconds duration;
    bool stopped;
    std::stringstream ss;
};

#endif  // CHRONO_H
