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
#include <limits>

extern "C" {
  void dsyev_ (char *JOBZ, char *UPLO, int *N, double *A, int *LDA, double *W, double *WORK, int *LWORK, int *INFO);
}

dfloat* gFiedler;

namespace paradogs {

/*Compute Fiedler vector of graph via multilevel heirarchy*/
dfloat* graph_t::FiedlerVector() {

  /*Fiedler vector on coarsest level*/
  L[Nlevels-1].FiedlerVector();

  /*Project and improve the Fiedler vector to the fine level*/
  for (int l=Nlevels-2;l>=0;--l) {
    Project(l);
  }

  // gFiedler = new dfloat[graph.Nverts];
  // for (dlong n=0;n<graph.Nverts;++n) gFiedler[n] = L[0].Fiedler[n];
  return L[0].Fiedler;
}

/*Compute Fiedler vector of graph Laplacian*/
void mgLevel_t::FiedlerVector() {

  int N = A.Nrows;

  /*Create full matrix from sparse csr*/
  double *M = new double[N*N];
  for (int n=0;n<N*N;++n) M[n] = 0.0;

  for (dlong n=0;n<N;++n) {
    const dlong start=A.diag.rowStarts[n];
    const dlong end=A.diag.rowStarts[n+1];
    for (dlong j=start;j<end;++j) {
      M[n*N+A.diag.cols[j]] = A.diag.vals[j];
    }
  }

  /*Call LaPack to find eigen pairs*/
  int INFO = -999;
  char JOBZ='V';
  char UPLO='L';
  int LWORK = -1;
  int LDA = N;
  double WORKSIZE=0.0;
  double *W= new double[N];
  dsyev_(&JOBZ, &UPLO, &N, M, &LDA, W, &WORKSIZE, &LWORK, &INFO); //Size query

  LWORK = int(WORKSIZE);
  double *WORK= new double[LWORK];
  dsyev_(&JOBZ, &UPLO, &N, M, &LDA, W, WORK, &LWORK, &INFO);
  delete[] WORK;

  if(INFO) {
    std::stringstream ss;
    ss << "Paradogs: dsyev_ reports info = " << INFO << " in FiedlerVector";
    LIBP_ABORT(ss.str());
  }

  /*Find the second smallest eigenvalue (the smallest is 0)*/
  double min0 = std::numeric_limits<double>::max();
  double min1 = std::numeric_limits<double>::max();
  int minloc0 = -1;
  int minloc1 = -1;
  for (int i=0;i<N;++i) {
    // printf("Eig[%d] = %f\n", i, W[i]);

    if (W[i]<min0) {
      min1 = min0;
      min0 = W[i];
      minloc1 = minloc0;
      minloc0 = i;
    } else if (W[i]<min1) {
      min1 = W[i];
      minloc1 = i;
    }
  }

  // printf("min1 = %f, minloc1 = %d \n", min1, minloc1);

  double* minV = M + minloc1*N;
  for (int i=0;i<N;++i) {
    Fiedler[i] = minV[i];
  }

  // Fiedler vector is already orthogonal to null

  /* Fiedler vector is probably already normalized, but just in case */
  dfloat norm = 0.0;
  for (dlong n=0;n<N;++n) norm += Fiedler[n]*Fiedler[n];
  norm = sqrt(norm);
  for (dlong n=0;n<N;++n) Fiedler[n] /= norm;

  delete[] W;
  delete[] M;
}

}

