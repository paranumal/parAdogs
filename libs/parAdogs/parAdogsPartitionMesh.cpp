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
#include "parAdogs/parAdogsPartition.hpp"

namespace paradogs {

std::mt19937 RNG; //Standard mersenne_twister_engine seeded with rank

void PartitionMesh(const dlong Nelements, const int Nfaces,
                   const dlong EToE[], const int EToP[],
                   MPI_Comm comm, int partition[]) {

  /* Create RNG*/
  paradogs::RNG = std::mt19937(0);

  /* Create graph from mesh info, and initialize multilevel Laplacian*/
  graph_t graph(Nelements, Nfaces, EToE, EToP, comm);

  const int Nparts=7;

  dfloat *targetFraction = new dfloat[Nparts];
  for (int n=0;n<Nparts;++n) targetFraction[n] = 1.0/Nparts;

  Partition(graph, Nparts, targetFraction, partition);

  delete[] targetFraction;

  /*Sanity check*/
  dfloat *partitionWeight = new dfloat[Nparts];
  for (int n=0;n<Nparts;++n) partitionWeight[n] = 0.0;

  dfloat testcut=0.0;
  parCSR &A = graph.L[0].A;

  for (dlong n=0;n<graph.Nverts;++n) {
    const dlong start = A.diag.rowStarts[n];
    const dlong end   = A.diag.rowStarts[n+1];

    const int partn = partition[n];
    partitionWeight[partn] += 1.0;

    for (dlong j=start;j<end;++j) {
      const int partj = partition[A.diag.cols[j]];
      if (partn != partj) {
        testcut += 1.0;
      }
    }
  }
  testcut /= 2.0;

  printf("Stats: cut = %f, Balance = { ", testcut);
  for (int n=0;n<Nparts;++n) printf(" %f, ", partitionWeight[n]/A.Nrows);
  printf("}\n");

  delete[] partitionWeight;
}

}

