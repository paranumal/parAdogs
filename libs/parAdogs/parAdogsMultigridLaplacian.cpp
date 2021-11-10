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

/*Create graph Laplacian*/
void mgLevel_t::CreateLaplacian(const dlong Nelements, 
                                const int Nfaces,
                                const dlong* EToE, 
                                const int* EToP,
                                MPI_Comm comm) {

  /*Create a graph Laplacian from mesh info*/
  A.Nrows = Nelements;
  A.Ncols = Nelements;
  A.diag.rowStarts = new dlong[Nelements+1];

  #pragma omp parallel for
  for (dlong n=0;n<Nelements+1;++n) {
    A.diag.rowStarts[n] = 0;
  }

  A.diag.nnz=0;
  for (dlong e=0;e<Nelements;++e) {
    A.diag.rowStarts[e+1]++; /*Count diagonal*/
    A.diag.nnz++;

    for (int n=0;n<Nfaces;++n) {
      if (EToE[Nfaces*e + n]>-1) {
        A.diag.rowStarts[e+1]++; /*count connections per vert*/
        A.diag.nnz++;            /*count total connections*/
      }
    }
  }

  /*cumulative sum*/
  for (dlong e=0;e<Nelements;++e) {
    A.diag.rowStarts[e+1] += A.diag.rowStarts[e];
  }

  /*Build connectivity*/
  A.diagA = new dfloat[A.Nrows];
  A.diagInv = new dfloat[A.Nrows];
  A.diag.cols = new dlong[A.diag.nnz];
  A.diag.vals = new dfloat[A.diag.nnz];

  A.diag.nnz=0;
  for (dlong e=0;e<Nelements;++e) {
    A.diag.cols[A.diag.nnz] = e;
    dfloat& Ann = A.diag.vals[A.diag.nnz];
    A.diag.nnz++;

    Ann = 0.0;

    for (int n=0;n<Nfaces;++n) {
      if (EToE[Nfaces*e + n]>-1) {
        A.diag.cols[A.diag.nnz] = EToE[Nfaces*e + n];
        A.diag.vals[A.diag.nnz] = -1.0;
        A.diag.nnz++;
        Ann += 1.0;
      }
    }
    A.diagA[e] = Ann;
    A.diagInv[e] = 1.0/Ann;
  }

  /*Construct fine null vector*/
  null = new dfloat[Nelements];

  #pragma omp parallel for
  for (dlong n=0;n<Nelements;++n) {
    null[n] = 1.0/sqrt(Nelements);
  }

  /*Space for Fiedler*/
  Fiedler = new dfloat[Nelements];

  /*Multigrid buffers*/
  RES = new dfloat[Nelements];
}

/*Split a graph Laplacian in two based on a partitioning*/
void mgLevel_t::SplitLaplacian(const int partition[],
                                mgLevel_t &LL, dlong mapL[],
                                mgLevel_t &LR, dlong mapR[]) {

  dlong* map = new dlong[A.Nrows];

  parCSR &AL = LL.A;
  parCSR &AR = LR.A;

  AL.diag.rowStarts = new dlong[AL.Nrows+1];
  AR.diag.rowStarts = new dlong[AR.Nrows+1];

  #pragma omp parallel for
  for (dlong n=0;n<AL.Nrows+1;++n) {
    AL.diag.rowStarts[n] = 0;
  }
  #pragma omp parallel for
  for (dlong n=0;n<AR.Nrows+1;++n) {
    AR.diag.rowStarts[n] = 0;
  }

  AL.Nrows=0;
  AR.Nrows=0;
  for (dlong n=0;n<A.Nrows;++n) {
    if (partition[n]==0) {
      map[n] = AL.Nrows;
      mapL[AL.Nrows++] = n;
      AL.diag.rowStarts[AL.Nrows]++; //count diagonal

      const dlong start = A.diag.rowStarts[n];
      const dlong end   = A.diag.rowStarts[n+1];
      for (dlong j=start;j<end;++j) {
        const dlong k = A.diag.cols[j];
        if (partition[k]==0) AL.diag.rowStarts[AL.Nrows]++;
      }
    } else {
      map[n] = AR.Nrows;
      mapR[AR.Nrows++] = n;
      AR.diag.rowStarts[AR.Nrows]++; //count diagonal

      const dlong start = A.diag.rowStarts[n];
      const dlong end   = A.diag.rowStarts[n+1];
      for (dlong j=start;j<end;++j) {
        const dlong k = A.diag.cols[j];
        if (partition[k]==1) AR.diag.rowStarts[AR.Nrows]++;
      }
    }
  }

  for (dlong n=0;n<AL.Nrows;++n) {
    AL.diag.rowStarts[n+1] += AL.diag.rowStarts[n];
  }
  for (dlong n=0;n<AR.Nrows;++n) {
    AR.diag.rowStarts[n+1] += AR.diag.rowStarts[n];
  }

  AL.Ncols=AL.Nrows;
  AR.Ncols=AR.Nrows;

  AL.diag.nnz=AL.diag.rowStarts[AL.Nrows];
  AR.diag.nnz=AR.diag.rowStarts[AR.Nrows];

  AL.diag.cols = new dlong[AL.diag.nnz];
  AR.diag.cols = new dlong[AR.diag.nnz];

  AL.diag.vals = new dfloat[AL.diag.nnz];
  AR.diag.vals = new dfloat[AR.diag.nnz];

  AL.diagA = new dfloat[AL.Nrows];
  AR.diagA = new dfloat[AR.Nrows];
  AL.diagInv = new dfloat[AL.Nrows];
  AR.diagInv = new dfloat[AR.Nrows];

  AL.diag.nnz=0;
  for (dlong nn=0;nn<AL.Nrows;++nn) {
    const dlong n = mapL[nn];
    AL.diag.cols[AL.diag.nnz] = nn;
    dfloat& Ann = AL.diag.vals[AL.diag.nnz];
    AL.diag.nnz++;

    Ann = 0.0;

    const dlong start = A.diag.rowStarts[n];
    const dlong end   = A.diag.rowStarts[n+1];
    for (dlong j=start;j<end;++j) {
      const dlong k = A.diag.cols[j];
      if (partition[k]==0) {
        AL.diag.cols[AL.diag.nnz] = map[k];
        AL.diag.vals[AL.diag.nnz] = -1.0;
        Ann += 1.0;
        AL.diag.nnz++;
      }
    }
    AL.diagA[nn] = Ann;
    AL.diagInv[nn] = 1.0/Ann;
  }

  AR.diag.nnz=0;
  for (dlong nn=0;nn<AR.Nrows;++nn) {
    const dlong n = mapR[nn];
    AR.diag.cols[AR.diag.nnz] = nn;
    dfloat& Ann = AR.diag.vals[AR.diag.nnz];
    AR.diag.nnz++;

    Ann = 0.0;

    const dlong start = A.diag.rowStarts[n];
    const dlong end   = A.diag.rowStarts[n+1];
    for (dlong j=start;j<end;++j) {
      const dlong k = A.diag.cols[j];
      if (partition[k]==1) {
        AR.diag.cols[AR.diag.nnz] = map[k];
        AR.diag.vals[AR.diag.nnz] = -1.0;
        Ann += 1.0;
        AR.diag.nnz++;
      }
    }
    AR.diagA[nn] = Ann;
    AR.diagInv[nn] = 1.0/Ann;
  }

  /*Construct null vectors*/
  LL.null = new dfloat[AL.Nrows];
  #pragma omp parallel for
  for (dlong n=0;n<AL.Nrows;++n) {
    LL.null[n] = 1.0/sqrt(AL.Nrows);
  }
  LR.null = new dfloat[AR.Nrows];
  #pragma omp parallel for
  for (dlong n=0;n<AR.Nrows;++n) {
    LR.null[n] = 1.0/sqrt(AR.Nrows);
  }

  /*Space for Fiedler*/
  LL.Fiedler = new dfloat[AL.Nrows];
  LR.Fiedler = new dfloat[AR.Nrows];

  /*Multigrid buffers*/
  LL.RES = new dfloat[AL.Nrows];
  LR.RES = new dfloat[AR.Nrows];

  delete[] map;
}

void mgLevel_t::Free() {
  A.Free();
  P.Free();
  R.Free();
  if (null) {delete[] null; null=nullptr; }
  if (Fiedler) {delete[] Fiedler; Fiedler=nullptr; }
  if (X) {delete[] X; X=nullptr; }
  if (RHS) {delete[] RHS; RHS=nullptr; }
  if (RES) {delete[] RES; RES=nullptr; }
  if (scratch) {delete[] scratch; scratch=nullptr; }
}

}
