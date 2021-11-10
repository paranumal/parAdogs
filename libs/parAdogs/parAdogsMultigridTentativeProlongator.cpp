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

parCSR TentativeProlongator(const dlong Nf,
                            const dlong Nc,
                            dlong FineToCoarse[],
                            dfloat FineNull[],
                            dfloat CoarseNull[]) {
  dlong nnz = Nf;
  nonZero_t *entries = new nonZero_t[nnz];

  /* Each entry is the CoarseNull vector entry*/
  #pragma omp parallel for
  for (dlong n=0;n<Nf;++n) {
    entries[n].row = n;
    entries[n].col = FineToCoarse[n];
    entries[n].val = FineNull[n];
  }

  parCSR T(Nf, Nc, nnz, entries);
  delete[] entries;

  /*Init coarse null*/
  #pragma omp parallel for
  for (dlong v=0;v<Nc;++v) CoarseNull[v] = 0.0;

  /*Sum columns of T*/
  //add local nonzeros
  for (dlong n=0;n<T.diag.nnz;++n)
    CoarseNull[T.diag.cols[n]] += T.diag.vals[n] * T.diag.vals[n];

  //add nonlocal nonzeros
  // for(dlong i=0; i<P->offd.nnz; i++)
  //   CoarseNull[P->offd.cols[i]] += P->offd.vals[i] * P->offd.vals[i];

  //add the halo values to their origins
  // P->halo->Combine(CoarseNull, 1, ogs::Dfloat);

  #pragma omp parallel for
  for (dlong n=0;n<Nc;++n)
    CoarseNull[n] = sqrt(CoarseNull[n]);

  //share the results
  // P->halo->Exchange(CoarseNull, 1, ogs::Dfloat);

  #pragma omp parallel for
  for (dlong n=0;n<T.diag.nnz;++n)
    T.diag.vals[n] /= CoarseNull[T.diag.cols[n]];
  // for(dlong v=0; v<Nf; v++)
  //   P[v] /= cgraph.CoarseNull[FtoC[v]];

  return T;
}

}


