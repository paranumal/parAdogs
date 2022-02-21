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
#include "parAdogs/parAdogsMatrix.hpp"

namespace libp {

namespace paradogs {

class mgLevel_t {
public:
  dlong Nrows=0, Ncols=0;
  hlong Nglobal=0;

  parCSR A, P, R;

  /*null vector*/
  libp::memory<dfloat> null;

  /*Fiedler vector*/
  libp::memory<dfloat> Fiedler;

  /*Vcycle storage*/
  libp::memory<dfloat> RHS;
  libp::memory<dfloat> X;
  libp::memory<dfloat> RES;
  libp::memory<dfloat> scratch;

  dfloat lambda1, lambda0; //smoothing params

  /*Create graph Laplacian*/
  void CreateLaplacian(const dlong Nelements,
                       const int Nfaces,
                       const libp::memory<hlong>& EToE,
                       MPI_Comm comm);

  /*Construct a coarse level*/
  void CoarsenLevel(mgLevel_t &Lf, const dfloat theta);

  void SetupSmoother();

  void AllocateScratch(const int l);

  /*Compute Fiedler vector directly*/
  void FiedlerVector();

  /*Multigrid functions*/
  void Smooth(libp::memory<dfloat>& r, libp::memory<dfloat>& x, const bool xIsZero);
  void Residual(libp::memory<dfloat>& r, libp::memory<dfloat>& x, libp::memory<dfloat>& res);
  void Coarsen(libp::memory<dfloat>& x, libp::memory<dfloat>& xC);
  void Prolongate(libp::memory<dfloat>& xC, libp::memory<dfloat>& x);
};

parCSR TentativeProlongator(const dlong Nf,
                            const dlong Nc,
                            platform_t& platform,
                            MPI_Comm comm,
                            libp::memory<hlong>& FineToCoarse,
                            libp::memory<dfloat>& FineNull,
                            libp::memory<dfloat>& CoarseNull);

parCSR SmoothProlongator(const parCSR& A,
                         const parCSR& T);

parCSR Transpose(const parCSR& A);

parCSR SpMM(const parCSR& A, const parCSR& B);

class coarseSolver_t {

public:
  MPI_Comm comm;

  int N=0;
  int Nrows=0;
  int Ncols=0;

  int coarseTotal=0;
  libp::memory<int> coarseCounts;
  libp::memory<int> coarseOffsets;

  libp::memory<dfloat> invA;
  libp::memory<dfloat> grhs;

  void Setup(parCSR& A, libp::memory<dfloat>& null);
  void Solve(libp::memory<dfloat>& r, libp::memory<dfloat>& x);
};

} //namespace paradogs

} //namespace libp

#endif

