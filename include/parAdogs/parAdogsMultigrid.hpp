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

#ifndef PARADOGS_MULTIGRID_HPP
#define PARADOGS_MULTIGRID_HPP 1

#include "parAdogs.hpp"
#include "parAdogs/parAdogsDefines.h"
#include "parAdogs/parAdogsMatrix.hpp"

namespace paradogs {

class mgLevel_t {
public:
  // parCSR A, P, R;
  parCSR A, R, P;

  /*null vector*/
  dfloat *null=nullptr;

  /*Fiedler vector*/
  dfloat *Fiedler=nullptr;

  /*Vcycle storage*/
  dfloat *RHS=nullptr;
  dfloat *X=nullptr;
  dfloat *RES=nullptr;
  dfloat *scratch=nullptr;

  dfloat lambda1, lambda0; //smoothing params

  ~mgLevel_t();

  /*Create graph Laplacian*/
  void CreateLaplacian(const dlong Nelements,
                       const int Nfaces,
                       const dlong* EToE,
                       const int* EToP,
                       MPI_Comm comm);

  /*Split a graph Laplacian in two based on a partitioning*/
  void SplitLaplacian(const int partition[],
                      mgLevel_t &L, dlong mapL[],
                      mgLevel_t &R, dlong mapR[]);

  /*Construct a coarse level*/
  void CoarsenLevel(mgLevel_t &Lf);

  void SetupSmoother();

  /*Compute Fiedler vector directly*/
  void FiedlerVector();

  /*Multigrid functions*/
  void Smooth(dfloat r[], dfloat x[], const bool xIsZero);
  void Residual(dfloat r[], dfloat x[], dfloat res[]);
  void Coarsen(dfloat x[], dfloat xC[]);
  void Prolongate(dfloat xC[], dfloat x[]);
};

parCSR& TentativeProlongator(const dlong Nf,
                             const dlong Nc,
                             dlong FineToCoarse[],
                             dfloat FineNull[],
                             dfloat CoarseNull[]);

parCSR& SmoothProlongator(const parCSR& A,
                          const parCSR& T);

parCSR& Transpose(const parCSR& A);

parCSR& SpMM(const parCSR& A, const parCSR& B);

class coarseSolver_t {

public:
  int Nrows=0;
  int Ncols=0;

  int coarseTotal;
  int coarseOffset;
  int *coarseOffsets=nullptr;
  int *coarseCounts=nullptr;
  int *sendOffsets=nullptr;
  int *sendCounts=nullptr;

  int N;
  int offdTotal=0;

  dfloat *diagInvA=nullptr, *offdInvA=nullptr;
  dfloat *diagRhs=nullptr, *offdRhs=nullptr;

  ~coarseSolver_t();

  void Setup(parCSR &A, dfloat null[]);
  void Solve(dfloat r[], dfloat x[]);
};


}

#endif

