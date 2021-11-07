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

#include "parAdogs.hpp"
#include "parAdogs/parAdogsGraph.hpp"

namespace paradogs {

/*Build a graph from mesh connectivity info*/
graph_t::graph_t(const dlong Nelements, const int Nfaces,
                 const dlong* EToE, const int* EToP,
                 MPI_Comm comm) {

  Nverts = Nelements;

  /*Create graph Laplacian from mesh*/
  Nlevels=1;
  L[0].CreateLaplacian(Nelements, Nfaces, EToE, EToP, comm);

  /*Create multilevel heirarchy*/
  MultigridSetup();
}

/*Divide graph into two pieces according to a bipartition*/
void graph_t::Split(const int partition[],
                    graph_t &graphL, dlong* &mapL,
                    graph_t &graphR, dlong* &mapR) {

  graphL.Nverts=0;
  graphR.Nverts=0;
  for (dlong n=0;n<Nverts;++n) {
    if (partition[n]==0) graphL.Nverts++;
    else                 graphR.Nverts++;
  }

  graphL.Nlevels=1;
  graphR.Nlevels=1;
  graphL.L[0].A.Nrows=graphL.Nverts;
  graphR.L[0].A.Nrows=graphR.Nverts;
  mapL = new dlong[graphL.Nverts];
  mapR = new dlong[graphR.Nverts];
  L[0].SplitLaplacian(partition,
                      graphL.L[0], mapL,
                      graphR.L[0], mapR);

  graphL.MultigridSetup();
  graphR.MultigridSetup();
}

}
