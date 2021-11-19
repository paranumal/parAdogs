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
#include <algorithm>
#include <limits>

#ifdef GLIBCXX_PARALLEL
#include <parallel/algorithm>
using __gnu_parallel::partition;
#else
using std::partition;
#endif

namespace paradogs {

static dfloat Pivot(dfloat A[],
                    const dlong left,
                    const dlong right,
                    const hlong k,
                    const dfloat min,
                    const dfloat max,
                    MPI_Comm comm) {
  /*Start with guessing a pivot halfway between min and max*/
  const dfloat pivot = (min+max)/2.0;

  /*Bail out if we're looking at a tiny window*/
  if (max-min < 1.0E-13) return pivot;

  dfloat* Am = partition(A+left, A+right, [pivot](const dfloat& a){ return a <= pivot; });

  /*Get how many entries are globally <= pivot*/
  hlong localCnt = Am-A;
  hlong globalCnt = localCnt;
  MPI_Allreduce(&localCnt, &globalCnt, 1, MPI_HLONG, MPI_SUM, comm);

  if (globalCnt==k) return pivot;

  if (k<globalCnt) {
    return Pivot(A, left, localCnt, k, min, pivot, comm);
  } else {
    return Pivot(A, localCnt, right, k, pivot, max, comm);
  }
}

/* Given a distributed vector F in comm, find a pivot value,
   such that there are globally k entries of F which are <= pivot. */
dfloat ParallelPivot(const dlong N, dfloat F[], 
                     const hlong k, MPI_Comm comm) {

  /*Make a copy of input vector*/
  dfloat *A = new dfloat[N];
  
  #pragma omp parallel for
  for (dlong n=0;n<N;++n) {
    A[n] = F[n];
  }

  /*Find global minimum/maximum*/
  dfloat localMin=std::numeric_limits<dfloat>::max();
  dfloat localMax=std::numeric_limits<dfloat>::min();
  for (dlong n=0;n<N;++n) {
    localMax = std::max(A[n], localMax);
    localMin = std::min(A[n], localMin);
  }
  dfloat globalMin=localMin;
  dfloat globalMax=localMax;
  MPI_Allreduce(&localMin, &globalMin, 1, MPI_DFLOAT, MPI_MIN, comm);
  MPI_Allreduce(&localMax, &globalMax, 1, MPI_DFLOAT, MPI_MAX, comm);

  /*Find pivot point via binary search*/
  dfloat pivot = Pivot(A, 0, N, k, globalMin, globalMax, comm);

  delete[] A;

  return pivot;
}

}

