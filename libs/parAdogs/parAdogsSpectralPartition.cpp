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

/*************************************************/
/* k-Way Recusive Spectral Partitioning          */
/*************************************************/

void graph_t::SpectralPartition() {

  const int Nparts=size;

  if (Nparts==1) return;

  /*Determine size of left and right partitions*/
  const int NpartsL = (Nparts+1)/2;
  const int NpartsR = Nparts-NpartsL;

  /*Set target */
  dfloat bipartitionFraction[2] = {0.0, 0.0};
  bipartitionFraction[0] = static_cast<dfloat>(NpartsL)/Nparts;
  bipartitionFraction[1] = 1.0 - bipartitionFraction[0];

  SpectralBipartition(bipartitionFraction);

  // if (Nparts>2) {
  //   /*Split graph according to partitioning*/
  //   graph_t graphL, graphR;
  //   dlong *mapL, *mapR;

  //   graph.Split(partition,
  //               graphL, mapL,
  //               graphR, mapR);

  //   /*Recursive calls*/
  //   if (NpartsL>1) {
  //     int* partitionL = new int[graphL.Nverts];
  //     dfloat *targetFractionL = new dfloat[NpartsL];

  //     for (int n=0;n<NpartsL;++n) {
  //       targetFractionL[n] = targetFraction[n]/bisectionFraction[0];
  //     }

  //     SpectralPartition(graphL, NpartsL, targetFractionL, partitionL);

  //     /*Inject partitioning*/
  //     for (dlong n=0;n<graphL.Nverts;++n){
  //       partition[mapL[n]] = partitionL[n];
  //     }

  //     delete[] targetFractionL;
  //     delete[] partitionL;
  //   }

  //   if (NpartsR>1) {
  //     int* partitionR = new int[graphR.Nverts];
  //     dfloat *targetFractionR = new dfloat[NpartsR];

  //     for (int n=0;n<NpartsR;++n) {
  //       targetFractionR[n] = targetFraction[n+NpartsL]/bisectionFraction[1];
  //     }

  //     SpectralPartition(graphR, NpartsR, targetFractionR, partitionR);

  //     /*Inject partitioning*/
  //     for (dlong n=0;n<graphR.Nverts;++n){
  //       partition[mapR[n]] = partitionR[n] + NpartsL;
  //     }

  //     delete[] targetFractionR;
  //     delete[] partitionR;
  //   } else if (NpartsL>1) {
  //     /*Inject partitioning*/
  //     for (dlong n=0;n<graphR.Nverts;++n){
  //       partition[mapR[n]] = 1 + NpartsL-1;
  //     }

  //   }
  //   delete[] mapL;
  //   delete[] mapR;
  // }
}

}

