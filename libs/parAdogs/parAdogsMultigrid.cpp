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
#include "parAdogs/parAdogsMultigrid.hpp"

namespace paradogs {

void mgLevel_t::Residual(dfloat r[], dfloat x[], dfloat res[]) {
  A.SpMV(-1.0, x, 1.0, r, res);
}

void mgLevel_t::Coarsen(dfloat x[], dfloat xC[]) {
  for (dlong n=0;n<Nc;n++) {
    xC[n]=0.0;
  }

  for (dlong n=0;n<A.Nrows;n++) {
    xC[FineToCoarse[n]] += P[n]*x[n];
  }  
}

void mgLevel_t::Prolongate(dfloat xC[], dfloat x[]) {
  for (dlong n=0;n<A.Nrows;n++) {
    x[n] += P[n]*xC[FineToCoarse[n]];
  }
}

void mgLevel_t::Smooth(dfloat r[], dfloat x[], const bool xIsZero) {
  const int ChebyshevIterations=2;
  A.smoothChebyshev(r, x, lambda0, lambda1,
                    xIsZero, scratch,
                    ChebyshevIterations);
}

void coarseSolver_t::Solve(dfloat rhs[], dfloat x[]) {

  //queue local part of gemv
  // const dfloat one=1.0;
  // const dfloat zero=0.0;
  // if (N)
  //   dGEMVKernel(N,N,one,o_diagInvAT,o_rhs, zero, o_x);

  // if(offdTotal) {
  //   //wait for data to arrive on host
  //   platform.device.setStream(ogs::dataStream);
  //   platform.device.finish();


  //   //gather the offd rhs entries
  //   MPI_Alltoallv(diagRhs,   sendCounts,   sendOffsets, MPI_DFLOAT,
  //                 offdRhs, coarseCounts, coarseOffsets, MPI_DFLOAT, comm);

  //   //queue transfering coarse vector to device
  //   o_offdRhs.copyFrom(offdRhs, offdTotal*sizeof(dfloat), 0, "async: true");
  //   platform.device.finish(); //wait for transfer to complete

  //   platform.device.setStream(currentStream);

  //   //queue offd part of gemv
  //   if (N)
  //     dGEMVKernel(N,offdTotal, one, o_offdInvAT,o_offdRhs, one, o_x);
  // }

  for (int n=0;n<N;++n) {
    dfloat xn=0.0;
    for (int m=0;m<N;++m) {
      xn += diagInvA[n*N + m]*rhs[m];
    }
    x[n] = xn;
  }
}

}
