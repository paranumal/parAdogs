/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#ifndef PARADOGS_GRAPH_HPP
#define PARADOGS_GRAPH_HPP 1

#include "parAdogs.hpp"
#include "parAdogs/parAdogsDefines.h"
#include "parAdogs/parAdogsMatrix.hpp"
#include "parAdogs/parAdogsMultigrid.hpp"

namespace paradogs {

class graph_t {
public:
  dlong Nverts=0;
  hlong NVertsGlobal=0;
  
  /*Multilevel Laplacian*/
  int Nlevels=0;
  mgLevel_t L[PARADOGS_MAX_LEVELS];
  coarseSolver_t coarseSolver;

  graph_t() {};

  /*Build a graph from mesh connectivity info*/
  graph_t(const dlong Nelements,
          const int Nfaces,
          const dlong* EToE,
          const int* EToP,
          MPI_Comm comm);


  /*Divide graph into two pieces according to a bisection*/
  void Split(const int partition[],
             graph_t &graphL, dlong* &mapL,
             graph_t &graphR, dlong* &mapR);

  /*Compute Fiedler vector of graph */
  dfloat* FiedlerVector();

  /*Improve a Fiedler vector*/
  void Refine(const int level);

  /* Solve A_{l}*x = b*/
  int Solve(const int level,
            const dfloat TOL,
            dfloat r[],
            dfloat x[],
            dfloat scratch[]);

  /*Create multilevel heirarchy*/
  void MultigridSetup();

  void MultigridVcycle(const int l,
                      dfloat r[],
                      dfloat x[]);

  /*Clear multilevel heirarchy*/
  void MultigridDestroy();
};

}

#endif

