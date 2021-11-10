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

#ifndef PARADOGS_MATRIX_HPP
#define PARADOGS_MATRIX_HPP 1

#include "parAdogs.hpp"
#include "parAdogs/parAdogsDefines.h"
#include "ogs.hpp"

namespace paradogs {

struct nonZero_t {
  hlong row;
  hlong col;
  dfloat val;
};


class parCSR {
public:
  dlong Nrows=0;
  dlong Ncols=0;

  //local sparse matrix
  struct CSR {
    dlong nnz=0;
    dlong  *rowStarts=nullptr;
    dlong  *cols=nullptr;
    pfloat *vals=nullptr;
  };
  CSR diag;

  //non-local sparse matrix
  struct MCSR {
    dlong nnz=0;
    dlong nzRows=0;

    dlong  *rowStarts=nullptr;
    dlong  *mRowStarts=nullptr; //compressed version of rowStarts
    dlong  *rows=nullptr;
    dlong  *cols=nullptr;
    pfloat *vals=nullptr;
  };
  MCSR offd;

  dfloat *diagA=nullptr;
  dfloat *diagInv=nullptr;

  //partition info
  hlong *globalRowStarts=nullptr;
  hlong *globalColStarts=nullptr;
  hlong *colMap=nullptr;

  ogs::halo_t *halo = nullptr;
  dlong NlocalCols = 0;

  //rho ~= cond(invD * A)
  dfloat rho=0.0;

  parCSR() {};
  parCSR(dlong N, dlong M): Nrows(N), Ncols(M) {}

  //build a parCSR matrix from a distributed COO matrix
  parCSR(dlong _Nrows, dlong _Ncols,
         const dlong NNZ,
         nonZero_t entries[]);

  /*Assignment - Copy & swap idiom*/
  parCSR& operator=(parCSR A);

  ~parCSR() {Free();}
  void Free();

  void haloSetup(hlong *colIds);

  // estimate rho(invD * A)
  dfloat rhoDinvA(dfloat null[]);

  /*Aggregate via distance-2 PMIS*/
  void Aggregate(dlong& cNverts,
                 dlong FineToCoarse[]);

  void GalerkinProduct(const parCSR &A, const parCSR &P);

  void SpMV(const dfloat alpha, dfloat x[],
            const dfloat beta, dfloat y[]);
  void SpMV(const dfloat alpha, dfloat x[],
            const dfloat beta, const dfloat y[], dfloat z[]);

  void SmoothChebyshev(dfloat b[], dfloat x[],
                       const dfloat lambda0, const dfloat lambda1,
                       const bool xIsZero, dfloat scratch[],
                       const int ChebyshevIterations);
};

}

#endif

