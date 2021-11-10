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

#ifdef GLIBCXX_PARALLEL_SORT
#include <parallel/algorithm>
using __gnu_parallel::sort;
#else
using std::sort;
#endif

namespace paradogs {

parCSR Transpose(const parCSR& A) {

  // MPI info
  // int rank, size;
  // MPI_Comm_rank(A.comm, &rank);
  // MPI_Comm_size(A.comm, &size);

  // copy data from nonlocal entries into send buffer
  // parCOO::nonZero_t *sendNonZeros = (parCOO::nonZero_t *)
  //                                   calloc(A.offd.nnz, sizeof(parCOO::nonZero_t));
  // for(dlong i=0;i<A.offd.nzRows;++i){
  //   const hlong row = A.offd.rows[i] + A.globalRowStarts[rank]; //global ids
  //   for (dlong j=A.offd.mRowStarts[i];j<A.offd.mRowStarts[i+1];j++) {
  //     const hlong col =  A.colMap[A.offd.cols[j]]; //global ids
  //     sendNonZeros[j].row = col;
  //     sendNonZeros[j].col = row;
  //     sendNonZeros[j].val = A.offd.vals[j];
  //   }
  // }

  // //sort by destination row
  // std::sort(sendNonZeros, sendNonZeros+A.offd.nnz,
  //           [](const parCOO::nonZero_t& a, const parCOO::nonZero_t& b) {
  //             if (a.row < b.row) return true;
  //             if (a.row > b.row) return false;

  //             return a.col < b.col;
  //           });

  // //count number of non-zeros we're sending
  // int *sendCounts = (int*) calloc(size, sizeof(int));
  // int *recvCounts = (int*) calloc(size, sizeof(int));
  // int *sendOffsets = (int*) calloc(size+1, sizeof(int));
  // int *recvOffsets = (int*) calloc(size+1, sizeof(int));

  // int r=0;
  // for (dlong n=0;n<A.offd.nnz;n++) {
  //   dlong row = sendNonZeros[n].row;
  //   while(row>=A.globalColStarts[r+1]) r++;
  //   sendCounts[r]++;
  // }

  // MPI_Alltoall(sendCounts, 1, MPI_INT,
  //              recvCounts, 1, MPI_INT, A.comm);

  // for (r=0;r<size;r++) {
  //   sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
  //   recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  // }
  // dlong offdnnz = recvOffsets[size]; //total offd nonzeros


  // parCOO cooAt(A.platform, A.comm);

  // //copy global partition
  // cooAt.globalRowStarts = (hlong *) calloc(size+1,sizeof(hlong));
  // cooAt.globalColStarts = (hlong *) calloc(size+1,sizeof(hlong));
  // memcpy(cooAt.globalRowStarts, A.globalColStarts, (size+1)*sizeof(hlong));
  // memcpy(cooAt.globalColStarts, A.globalRowStarts, (size+1)*sizeof(hlong));

  // cooAt.nnz = A.diag.nnz+offdnnz;

  dlong NNZ = A.diag.nnz;

  nonZero_t *entries = new nonZero_t[NNZ];

  dlong *rowStarts = new dlong[A.Ncols+1];
  dlong *rowCounts = new dlong[A.Ncols];
  #pragma omp parallel for
  for(dlong i=0; i<A.Ncols+1; i++) rowStarts[i]=0;

  #pragma omp parallel for
  for(dlong i=0; i<A.Ncols; i++) rowCounts[i]=0;

  /*Count entries per row*/
  for(dlong i=0; i<A.Nrows; i++) {
    const dlong Jstart = A.diag.rowStarts[i];
    const dlong Jend   = A.diag.rowStarts[i+1];

    for(dlong jj=Jstart; jj<Jend; jj++){
      rowStarts[A.diag.cols[jj]+1]++;
    }
  }

  /*Cumulative sum*/
  for(dlong i=1; i<A.Ncols+1; i++) {
    rowStarts[i] += rowStarts[i-1];
  }

  //fill local nonzeros
  // #pragma omp parallel for
  for(dlong i=0; i<A.Nrows; i++){
    const dlong Jstart = A.diag.rowStarts[i];
    const dlong Jend   = A.diag.rowStarts[i+1];

    for(dlong jj=Jstart; jj<Jend; jj++){
      const dlong row = A.diag.cols[jj];
      const dlong c = rowStarts[row] + rowCounts[row];

      // cooAt.entries[jj].row = A.diag.cols[jj] + A.globalColStarts[rank];
      entries[c].row = row;
      // cooAt.entries[jj].col = i + A.globalRowStarts[rank];
      entries[c].col = i;
      entries[c].val = A.diag.vals[jj];
      rowCounts[row]++;
    }
  }

  // receive non-local nonzeros
  // MPI_Alltoallv(sendNonZeros,              sendCounts, sendOffsets, MPI_NONZERO_T,
  //               cooAt.entries+A.diag.nnz, recvCounts, recvOffsets, MPI_NONZERO_T,
  //               A.comm);

  // //clean up
  // MPI_Barrier(A.comm);
  // free(sendNonZeros);
  // free(sendCounts);
  // free(recvCounts);
  // free(sendOffsets);
  // free(recvOffsets);

  //sort each row by column id
  #pragma omp parallel for
  for(dlong i=0; i<A.Ncols; i++){
    const dlong Nentries = rowStarts[i+1]-rowStarts[i];
    const dlong c = rowStarts[i];
    std::sort(entries+c, entries+c+Nentries,
          [](const nonZero_t& a, const nonZero_t& b) {
            return a.col < b.col;
          });
  }

  delete[] rowStarts;
  delete[] rowCounts;

  parCSR AT(A.Ncols, A.Nrows, NNZ, entries);

  delete[] entries;

  return AT;
}

}


