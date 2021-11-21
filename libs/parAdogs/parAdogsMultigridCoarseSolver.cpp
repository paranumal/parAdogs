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

void coarseSolver_t::Solve(dfloat rhs[], dfloat x[]) {

  //gather the global rhs
  MPI_Allgatherv(rhs, N, MPI_DFLOAT,
                 grhs, coarseCounts, coarseOffsets, MPI_DFLOAT, comm);

  #pragma omp parallel for
  for (int n=0;n<N;++n) {
    dfloat xn=0.0;
    for (int m=0;m<coarseTotal;++m) {
      xn += invA[n*coarseTotal + m]*grhs[m];
    }
    x[n] = xn;
  }
}

void coarseSolver_t::Setup(parCSR* A, dfloat null[]) {

  comm = A->comm;
  int size;
  MPI_Comm_size(comm, &size);

  N = static_cast<int>(A->Nrows);
  Nrows = A->Nrows;
  Ncols = A->Ncols;

  coarseCounts = new int[size];
  coarseOffsets = new int[size];

  //collect partitioning info
  MPI_Allgather(&N, 1, MPI_INT,
                coarseCounts, 1, MPI_INT, comm);

  coarseTotal=0;
  for (int r=0;r<size;++r) {
    coarseTotal+=coarseCounts[r];
  }
  coarseOffsets[0]=0;
  for (int r=1;r<size;++r) {
    coarseOffsets[r]= coarseOffsets[r-1] + coarseCounts[r-1];
  }

  //gather global null vector
  dfloat *gnull = new dfloat[coarseTotal];

  MPI_Allgatherv( null, N, MPI_DFLOAT,
                 gnull, coarseCounts, coarseOffsets, MPI_DFLOAT,
                 comm);

  //populate local dense matrix
  dfloat *localA = new dfloat[N*coarseTotal];

  /*Fill the matrix with the null boost*/
  #pragma omp parallel for
  for (int n=0;n<N;n++) {
    for (int m=0;m<coarseTotal;m++) {
      localA[n*coarseTotal+m] = null[n]*gnull[m];
    }
  }
  delete[] gnull;

  /*Add sparse entries*/
  #pragma omp parallel for
  for (int n=0;n<N;n++) {
    const int start = static_cast<int>(A->diag.rowStarts[n]);
    const int end   = static_cast<int>(A->diag.rowStarts[n+1]);
    for (int m=start;m<end;m++) {
      const int col = static_cast<int>(A->diag.cols[m] + A->colOffsetL);
      localA[n*coarseTotal+col] += A->diag.vals[m];
    }
  }
  #pragma omp parallel for
  for (int n=0;n<A->offd.nzRows;n++) {
    const int row   = static_cast<int>(A->offd.rows[n]);
    const int start = static_cast<int>(A->offd.mRowStarts[n]);
    const int end   = static_cast<int>(A->offd.mRowStarts[n+1]);
    for (int m=start;m<end;m++) {
      const int col = static_cast<int>(A->colMap[A->offd.cols[m]]);
      localA[row*coarseTotal+col] += A->offd.vals[m];
    }
  }

  //assemble the full matrix
  dfloat *gA = new dfloat[coarseTotal*coarseTotal];

  for (int r=0;r<size;++r) {
    coarseCounts[r] *= coarseTotal;
    coarseOffsets[r] *= coarseTotal;
  }

  MPI_Allgatherv(localA, N*coarseTotal, MPI_DFLOAT,
                 gA, coarseCounts, coarseOffsets, MPI_DFLOAT,
                 comm);

  MPI_Barrier(comm);
  delete[] localA;

  for (int r=0;r<size;++r) {
    coarseCounts[r]  /= coarseTotal;
    coarseOffsets[r] /= coarseTotal;
  }

  matrixInverse(coarseTotal, gA);

  //diag piece of invA
  invA = new dfloat[N*coarseTotal];

  #pragma omp parallel for
  for (int n=0;n<N;n++) {
    for (int m=0;m<coarseTotal;m++) {
      invA[n*coarseTotal+m] = gA[(n+A->rowOffsetL)*coarseTotal+m];
    }
  }

  delete[] gA;

  /*Space for global rhs*/
  grhs = new dfloat[coarseTotal];
}

void coarseSolver_t::Free() {
  N=0;
  Nrows=0;
  Ncols=0;
  coarseTotal=0;
  if(coarseOffsets) {delete[] coarseOffsets; coarseOffsets=nullptr;}
  if(coarseCounts)  {delete[] coarseCounts; coarseCounts=nullptr;}
  if(invA)          {delete[] invA; invA=nullptr; }
  if(grhs)          {delete[] grhs; grhs=nullptr; }
}

}


