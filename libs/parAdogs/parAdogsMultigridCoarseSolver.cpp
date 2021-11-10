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

  #pragma omp parallel for
  for (int n=0;n<N;++n) {
    dfloat xn=0.0;
    for (int m=0;m<N;++m) {
      xn += diagInvA[n*N + m]*rhs[m];
    }
    x[n] = xn;
  }
}

void coarseSolver_t::Setup(parCSR &A, dfloat null[]) {

  // comm = A->comm;
  // MPI_Comm_rank(comm, &rank);
  // MPI_Comm_size(comm, &size);

  //copy the global coarse partition as ints
  // coarseOffsets = (int* ) calloc(size+1,sizeof(int));
  // for (int r=0;r<size+1;r++) coarseOffsets[r] = (int) A->globalRowStarts[r];

  // coarseTotal   = coarseOffsets[size];
  // coarseOffset  = coarseOffsets[rank];

  N = (int) A.Nrows;
  Nrows = A.Nrows;
  Ncols = A.Ncols;

  coarseTotal=N;

  // coarseCounts = (int*) calloc(size,sizeof(int));

  // int sendNNZ = (int) (A->diag.nnz+A->offd.nnz);

  // if((rank==0)&&(settings.compareSetting("VERBOSE","TRUE")))
  //   {printf("Setting up coarse solver...");fflush(stdout);}

  // parCOO::nonZero_t *sendNonZeros = (parCOO::nonZero_t *) calloc(sendNNZ, sizeof(parCOO::nonZero_t));

  //populate matrix
  // int cnt = 0;
  // for (int n=0;n<N;n++) {
  //   const int start = (int) A->diag.rowStarts[n];
  //   const int end   = (int) A->diag.rowStarts[n+1];
  //   for (int m=start;m<end;m++) {
  //     sendNonZeros[cnt].row = n + coarseOffset;
  //     sendNonZeros[cnt].col = A->diag.cols[m] + coarseOffset;
  //     sendNonZeros[cnt].val = A->diag.vals[m];
  //     cnt++;
  //   }
  // }

  // for (int n=0;n<A->offd.nzRows;n++) {
  //   const int row   = (int) A->offd.rows[n];
  //   const int start = (int) A->offd.mRowStarts[n];
  //   const int end   = (int) A->offd.mRowStarts[n+1];
  //   for (int m=start;m<end;m++) {
  //     sendNonZeros[cnt].row = row + coarseOffset;
  //     sendNonZeros[cnt].col = A->colMap[A->offd.cols[m]];
  //     sendNonZeros[cnt].val = A->offd.vals[m];
  //     cnt++;
  //   }
  // }

  // //get the nonzero counts from all ranks
  // int *recvNNZ    = (int*) calloc(size,sizeof(int));
  // int *NNZoffsets = (int*) calloc(size+1,sizeof(int));
  // MPI_Allgather(&sendNNZ, 1, MPI_INT,
  //                recvNNZ, 1, MPI_INT, comm);

  // int totalNNZ = 0;
  // for (int r=0;r<size;r++) {
  //   totalNNZ += recvNNZ[r];
  //   NNZoffsets[r+1] = NNZoffsets[r] + recvNNZ[r];
  // }

  // parCOO::nonZero_t *recvNonZeros = (parCOO::nonZero_t *) calloc(totalNNZ, sizeof(parCOO::nonZero_t));

  // MPI_Allgatherv(sendNonZeros, sendNNZ,             MPI_NONZERO_T,
  //                recvNonZeros, recvNNZ, NNZoffsets, MPI_NONZERO_T, comm);

  // //gather null vector
  // dfloat *nullTotal = (dfloat*) calloc(coarseTotal,sizeof(dfloat));

  // for (int r=0;r<size;r++)
  //   coarseCounts[r] = coarseOffsets[r+1]-coarseOffsets[r];

  // MPI_Allgatherv(nullVector,            N,                MPI_DFLOAT,
  //                 nullTotal, coarseCounts, coarseOffsets, MPI_DFLOAT,
  //                comm);

  // //clean up
  // MPI_Barrier(comm);
  // free(sendNonZeros);
  // free(NNZoffsets);
  // free(recvNNZ);

  //assemble the full matrix
  dfloat *coarseA = new dfloat[coarseTotal*coarseTotal];
  // for (int i=0;i<totalNNZ;i++) {
  //   int n = recvNonZeros[i].row;
  //   int m = recvNonZeros[i].col;
  //   coarseA[n*coarseTotal+m] = recvNonZeros[i].val;
  // }

  // if (nullSpace) { //A is dense due to nullspace augmentation
  //   for (int n=0;n<coarseTotal;n++) {
  //     for (int m=0;m<coarseTotal;m++) {
  //       coarseA[n*coarseTotal+m] += nullSpacePenalty*nullTotal[n]*nullTotal[m];
  //     }
  //   }
  // }

  // free(recvNonZeros);
  // free(nullTotal);

  #pragma omp parallel for
  for (int n=0;n<coarseTotal;n++) {
    for (int m=0;m<coarseTotal;m++) {
      coarseA[n*coarseTotal+m] = 0.0;
    }
  }

  #pragma omp parallel for
  for (int n=0;n<A.Nrows;n++) {
    const dlong start = A.diag.rowStarts[n];
    const dlong end   = A.diag.rowStarts[n+1];
    for (int j=start;j<end;++j) {
      coarseA[n*N + A.diag.cols[j]] = A.diag.vals[j];
    }
  }

  /*Augment null space*/

  #pragma omp parallel for
  for (int n=0;n<coarseTotal;n++) {
    for (int m=0;m<coarseTotal;m++) {
      coarseA[n*coarseTotal+m] += null[n]*null[m];
    }
  }

  matrixInverse(coarseTotal, coarseA);

  //determine size of offd piece
  // offdTotal = coarseTotal - N;

  // //shift offsets for MPI_AllGatherv of offd pieces
  // for (int r=rank+1;r<=size;r++)
  //   coarseOffsets[r]-= N;

  // //dont copy the local piece in MPI_AllGatherv
  // coarseCounts[rank]=0;

  // //counts for all-to-all
  // sendCounts = (int* ) calloc(size,sizeof(int));
  // sendOffsets = (int* ) calloc(size,sizeof(int));
  // for (int r=0;r<size;r++) {
  //   sendCounts[r] = N;
  //   sendOffsets[r] = 0;
  // }
  // sendCounts[rank] = 0;

  // //diag piece of invA
  // diagInvAT = (dfloat *) calloc(N*N,sizeof(dfloat));
  // for (int n=0;n<N;n++) {
  //   for (int m=0;m<N;m++) {
  //     diagInvAT[n+m*N] = coarseA[(n+coarseOffset)*coarseTotal+(m+coarseOffset)];
  //   }
  // }

  // //offd piece of invA
  // offdInvAT = (dfloat *) calloc(N*offdTotal,sizeof(dfloat));
  // for (int n=0;n<N;n++) {
  //   for (int m=0;m<coarseOffset;m++) {
  //     offdInvAT[n+m*N] = coarseA[(n+coarseOffset)*coarseTotal+m];
  //   }
  //   for (int m=coarseOffset+N;m<coarseTotal;m++) {
  //     offdInvAT[n+(m-N)*N] = coarseA[(n+coarseOffset)*coarseTotal+m];
  //   }
  // }

  // o_diagInvAT = platform.malloc(N*N*sizeof(dfloat), diagInvAT);
  // o_offdInvAT = platform.malloc(N*offdTotal*sizeof(dfloat), offdInvAT);

  // diagRhs = (dfloat*) calloc(N,sizeof(dfloat));
  // offdRhs = (dfloat*) calloc(offdTotal,sizeof(dfloat));

  // o_offdRhs = platform.malloc(offdTotal*sizeof(dfloat));

  // if((rank==0)&&(settings.compareSetting("VERBOSE","TRUE"))) printf("done.\n");


  //diag piece of invA
  diagInvA = new dfloat[N*N];

  #pragma omp parallel for
  for (int n=0;n<N;n++) {
    for (int m=0;m<N;m++) {
      diagInvA[n*N+m] = coarseA[n*N+m];
    }
  }

  delete[] coarseA;
}

void coarseSolver_t::Free() {
  Nrows=0;
  Ncols=0;
  coarseTotal=0;
  coarseOffset=0;
  N=0;
  offdTotal=0;
  if(coarseOffsets) {delete[] coarseOffsets; coarseOffsets=nullptr;}
  if(coarseCounts)  {delete[] coarseCounts; coarseCounts=nullptr;}
  if(sendOffsets)   {delete[] sendOffsets; sendOffsets=nullptr;}
  if(sendCounts)    {delete[] sendCounts; sendCounts=nullptr;}
  if(diagInvA)      {delete[] diagInvA; diagInvA=nullptr; }
  if(offdInvA)      {delete[] offdInvA; offdInvA=nullptr;}
  if(diagRhs)       {delete[] diagRhs; diagRhs=nullptr; }
  if(offdRhs)       {delete[] offdRhs; offdRhs=nullptr;}
}

}


