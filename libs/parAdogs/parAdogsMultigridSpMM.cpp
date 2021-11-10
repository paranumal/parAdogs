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
#include "parAdogs/parAdogsMatrix.hpp"
#include "parAdogs/parAdogsPartition.hpp"

namespace paradogs {

struct entry_t {
  dlong col;
  dfloat val;

  entry_t() {}
  entry_t(const dlong _col, const dfloat _val):
    col{_col}, val{_val} {}

  friend bool operator<(const entry_t& a, const entry_t& b) {
    return a.col < b.col;
  }
};

parCSR SpMM(const parCSR& A, const parCSR& B){

  // MPI info
  // int rank, size;
  // MPI_Comm_rank(A.comm, &rank);
  // MPI_Comm_size(A.comm, &size);

  // To compute C = A*B we need all the rows B(j,:) for which
  // j is a column index for the nonzeros of A on this rank.
  // For all local column indices in A.diag, we will already
  // have the row of B on this rank, so we just need to gather
  // the offd colIds

  // hlong *recvRows = (hlong *) calloc(A.Ncols-A.NlocalCols, sizeof(hlong));
  // int *sendCounts = (int*) calloc(size, sizeof(int));
  // int *recvCounts = (int*) calloc(size, sizeof(int));
  // int *sendOffsets = (int*) calloc(size+1, sizeof(int));
  // int *recvOffsets = (int*) calloc(size+1, sizeof(int));

  // //use the colMap of A to list the needed rows of B
  // int r=0;
  // for (dlong n=A.NlocalCols;n<A.Ncols;n++) {
  //   const hlong id = A.colMap[n];
  //   while (id>=B.globalRowStarts[r+1]) r++; //assumes the halo is sorted
  //   recvCounts[r]++;
  //   recvRows[n-A.NlocalCols] = id; //record the row to recv
  // }

  // //share the counts
  // MPI_Alltoall(recvCounts, 1, MPI_INT,
  //              sendCounts, 1, MPI_INT, A.comm);

  // for (r=0;r<size;r++) {
  //   sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
  //   recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  // }

  // int sendTotal = sendOffsets[size];
  // hlong *sendRows = (hlong *) calloc(sendTotal, sizeof(hlong));

  // //share the rowIds
  // MPI_Alltoallv(recvRows, recvCounts, recvOffsets, MPI_HLONG,
  //               sendRows, sendCounts, sendOffsets, MPI_HLONG,
  //               B.comm);

  // //we now have a list of rows to send, count the nnz to send
  // dlong NNZ=0;
  // for (r=0;r<size;r++) {
  //   sendCounts[r] =0; //reset
  //   for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
  //     dlong i = (dlong) (sendRows[n]-B.globalRowStarts[rank]); //local row id
  //     sendCounts[r]+= B.diag.rowStarts[i+1]-B.diag.rowStarts[i]; //count entries in this row
  //     sendCounts[r]+= B.offd.rowStarts[i+1]-B.offd.rowStarts[i]; //count entries in this row
  //   }
  //   NNZ += sendCounts[r]; //tally the total
  // }

  // parCOO::nonZero_t *sendNonZeros = (parCOO::nonZero_t *) calloc(NNZ, sizeof(parCOO::nonZero_t));

  // NNZ=0; //reset
  // for (r=0;r<size;r++) {
  //   for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
  //     dlong i = (dlong) (sendRows[n] - B.globalRowStarts[rank]); //local row id
  //     for (dlong jj=B.diag.rowStarts[i]; jj<B.diag.rowStarts[i+1];jj++){
  //       sendNonZeros[NNZ].row = sendRows[n];
  //       sendNonZeros[NNZ].col = B.diag.cols[jj] + B.globalColStarts[rank];
  //       sendNonZeros[NNZ].val = B.diag.vals[jj];
  //       NNZ++;
  //     }
  //     for (dlong jj=B.offd.rowStarts[i]; jj<B.offd.rowStarts[i+1];jj++){
  //       sendNonZeros[NNZ].row = sendRows[n];
  //       sendNonZeros[NNZ].col = B.colMap[B.offd.cols[jj]];
  //       sendNonZeros[NNZ].val = B.offd.vals[jj];
  //       NNZ++;
  //     }
  //   }
  // }

  // MPI_Alltoall(sendCounts, 1, MPI_INT,
  //              recvCounts, 1, MPI_INT, A.comm);

  // for (r=0;r<size;r++) {
  //   sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
  //   recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  // }


  // dlong Boffdnnz = recvOffsets[size]; //total nonzeros
  // parCOO::nonZero_t *BoffdRows = (parCOO::nonZero_t *)
  //                                calloc(Boffdnnz, sizeof(parCOO::nonZero_t));

  // MPI_Alltoallv(sendNonZeros, sendCounts, sendOffsets, MPI_NONZERO_T,
  //               BoffdRows, recvCounts, recvOffsets, MPI_NONZERO_T,
  //               B.comm);

  // //clean up
  // MPI_Barrier(B.comm);
  // free(sendNonZeros);
  // free(sendCounts);
  // free(recvCounts);
  // free(sendOffsets);
  // free(recvOffsets);

  // //we now have all the needed nonlocal rows (should also be sorted by row then col)

  // //make an array of row offsets so we know how large each row is
  // dlong *BoffdRowOffsets = (dlong *) calloc(A.Ncols-A.NlocalCols+1, sizeof(dlong));

  // dlong id=0;
  // for (dlong n=0;n<Boffdnnz;n++) {
  //   hlong row = BoffdRows[n].row;

  //   while(A.colMap[id+A.NlocalCols]!=row) id++;

  //   BoffdRowOffsets[id+1]++; //count entry in row
  // }

  // //cumulative sum
  // for (dlong n=0;n<A.Ncols-A.NlocalCols;n++)
  //   BoffdRowOffsets[n+1] += BoffdRowOffsets[n];


  // The next step to compute C = A*B is to multiply each entry A(i,j) by the
  // row B(j,:), store the all the results, sort them by row+col, and compress
  // the entries

  // Find how big the intermediate form is
  dlong *rowStarts = new dlong[A.Nrows+1];
  dlong *rowCounts = new dlong[A.Nrows];

  #pragma omp parallel for
  for(dlong i=0; i<A.Nrows+1; i++) rowStarts[i]=0;

  #pragma omp parallel for
  for(dlong i=0; i<A.Nrows; i++) rowCounts[i]=0;

  #pragma omp parallel for
  for (dlong i=0;i<A.Nrows;i++) {
    //local entries
    dlong start = A.diag.rowStarts[i];
    dlong end   = A.diag.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A.diag.cols[j];
      const int nnzBj =  B.diag.rowStarts[col+1]-B.diag.rowStarts[col];
                        // +B.offd.rowStarts[col+1]-B.offd.rowStarts[col];
      rowStarts[i+1] += nnzBj;
    }
    //non-local entries
    // start = A.offd.rowStarts[i];
    // end   = A.offd.rowStarts[i+1];
    // for (dlong j=start;j<end;j++) {
    //   const dlong col = A.offd.cols[j]-A.NlocalCols;
    //   const int nnzBj = BoffdRowOffsets[col+1] - BoffdRowOffsets[col];
    //   NNZ += nnzBj;
    // }
  }

  /*Cumulative sum*/
  for(dlong i=1; i<A.Nrows+1; i++) {
    rowStarts[i] += rowStarts[i-1];
  }

  // NNZ = T.diag.nnz+T.offd.nnz; //start with T populated
  dlong NNZ = rowStarts[A.Nrows];

  nonZero_t *Ctmp = new nonZero_t[NNZ];

  //count total number of nonzeros;
  dlong nnz =0;

  // Fill the intermediate form of C
  // #pragma omp parallel for
  for (dlong i=0;i<A.Nrows;i++) {
    const dlong cStart = rowStarts[i];
    dlong& c = rowCounts[i];

    //local A entries
    dlong start = A.diag.rowStarts[i];
    dlong end   = A.diag.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A.diag.cols[j];
      const dfloat Aval = A.diag.vals[j];

      //local B entries
      dlong Bstart = B.diag.rowStarts[col];
      dlong Bend   = B.diag.rowStarts[col+1];
      for (dlong jj=Bstart;jj<Bend;jj++) {
        // Ctmp[cStart+c].row = i + A.globalRowStarts[rank];
        Ctmp[cStart+c].row = i;
        // Ctmp[cStart+c].col = B.diag.cols[jj]+B.globalColStarts[rank]; //global id
        Ctmp[cStart+c].col = B.diag.cols[jj]; //global id
        Ctmp[cStart+c].val = Aval*B.diag.vals[jj];
        c++;
      }
      //non-local B entries
      // Bstart = B.offd.rowStarts[col];
      // Bend   = B.offd.rowStarts[col+1];
      // for (dlong jj=Bstart;jj<Bend;jj++) {
      //   Ctmp[cStart+c].row = i + A.globalRowStarts[rank];
      //   Ctmp[cStart+c].col = B.colMap[B.offd.cols[jj]]; //global id
      //   Ctmp[cStart+c].val = Aval*B.offd.vals[jj];
      //   c++;
      // }
    }
    //non-local A entries
    // start = A.offd.rowStarts[i];
    // end   = A.offd.rowStarts[i+1];
    // for (dlong j=start;j<end;j++) {
    //   const dlong col = A.offd.cols[j]-A.NlocalCols;
    //   const dfloat Aval = A.offd.vals[j];

    //   // entries from recived rows of B
    //   dlong Bstart = BoffdRowOffsets[col];
    //   dlong Bend   = BoffdRowOffsets[col+1];
    //   for (dlong jj=Bstart;jj<Bend;jj++) {
    //     Ctmp[cStart+c].row = i + A.globalRowStarts[rank];
    //     Ctmp[cStart+c].col = BoffdRows[jj].col; //global id
    //     Ctmp[cStart+c].val = Aval*BoffdRows[jj].val;
    //     c++;
    //   }
    // }

    //sort entries in this row by col id
    std::sort(Ctmp+cStart, Ctmp+cStart+c,
              [](const nonZero_t& a, const nonZero_t& b) {
                return a.col < b.col;
              });

    /*Count how many actual nonzeros will be in this row*/
    dlong nnzRow=0;
    if (c>0) nnzRow++;
    for (dlong j=1;j<c;j++) {
      if ((Ctmp[cStart+j].col!=Ctmp[cStart+j-1].col)) nnzRow++;
    }

    nnz+=nnzRow; //Add to total
  }
  // free(BoffdRowOffsets);
  // free(BoffdRows);

  delete[] rowStarts;
  delete[] rowCounts;

  // parCOO cooC(A.platform, A.comm);

  // //copy global partition
  // cooC.globalRowStarts = (hlong *) calloc(size+1,sizeof(hlong));
  // cooC.globalColStarts = (hlong *) calloc(size+1,sizeof(hlong));
  // memcpy(cooC.globalRowStarts, A.globalRowStarts, (size+1)*sizeof(hlong));
  // memcpy(cooC.globalColStarts, B.globalColStarts, (size+1)*sizeof(hlong));

  // cooC.nnz = nnz;
  nonZero_t *entries = new nonZero_t[nnz];

  //compress nonzeros
  nnz = 0;
  if (NNZ) entries[nnz++] = Ctmp[0];
  for (dlong i=1;i<NNZ;i++) {
    if ((Ctmp[i].row!=Ctmp[i-1].row)||
        (Ctmp[i].col!=Ctmp[i-1].col)) {
      entries[nnz++] = Ctmp[i];
    } else {
      entries[nnz-1].val += Ctmp[i].val;
    }
  }
  //clean up
  delete [] Ctmp;

  //build C from coo matrix
  parCSR C(A.Nrows, B.Ncols, nnz, entries);

  delete[] entries;

  return C;
}

void parCSR::GalerkinProduct(const parCSR &A, const parCSR &P) {

  const dlong Nf = A.Nrows;
  const dlong Nc = P.Ncols;

  Nrows = P.Ncols;
  Ncols = P.Ncols;
  diag.rowStarts = new dlong[Nc+1];
  for (dlong n=0;n<Nc+1;++n) diag.rowStarts[n]=0;

  /*Make an array to hold all the uncompressed nonZeros*/
  entry_t* nonZeros = new entry_t[A.diag.nnz];
  dlong* nzStarts = new dlong[Nc+1];
  dlong* nzCounts = new dlong[Nc];
  #pragma omp parallel for
  for (dlong n=0;n<Nc+1;++n) nzStarts[n]=0;

  #pragma omp parallel for
  for (dlong n=0;n<Nc;++n) nzCounts[n]=0;

  /*Count how many coarse nonZeros will be in each list (upper bound)*/
  for (dlong v=0;v<Nf;++v) {
    const dlong c = P.diag.cols[v];
    nzStarts[c+1] += A.diag.rowStarts[v+1]-A.diag.rowStarts[v];
  }

  /*Cumulative sum*/
  for (dlong c=0;c<Nc;++c) {
    nzStarts[c+1] += nzStarts[c];
  }

  /*Fill the nonzero list*/
  for (dlong v=0;v<Nf;++v) {
    const dlong c = P.diag.cols[v];

    entry_t* nonZerosn = nonZeros + nzStarts[c] + nzCounts[c];

    /*nonZeros*/
    int cnt=0;
    const dlong start = A.diag.rowStarts[v];
    const dlong end   = A.diag.rowStarts[v+1];
    for (dlong j=start;j<end;++j) {
      const dlong k = A.diag.cols[j];
      const dlong cj = P.diag.cols[k];
      nonZerosn[cnt++] = entry_t(cj, A.diag.vals[j]*P.diag.vals[k]*P.diag.vals[v]);
    }

    nzCounts[c]+=cnt;
  }

  /*Sort each row's entries and count*/
  for (dlong c=0;c<Nc;++c) {
    entry_t* nonZerosn = nonZeros + nzStarts[c];
    const int Nentries = nzCounts[c];

    std::sort(nonZerosn, nonZerosn+Nentries);

    /*Count the real number of nonZeros for this coarse node*/
    int cnt=0;
    if (Nentries>0) cnt++;
    for (int j=1;j<Nentries;++j) {
      if (nonZerosn[j].col != nonZerosn[j-1].col) cnt++;
    }
    diag.rowStarts[c+1] = cnt;
  }

  /*Cumulative sum*/
  for (dlong c=0;c<Nc;++c) {
    diag.rowStarts[c+1] += diag.rowStarts[c];
  }

  diag.nnz = diag.rowStarts[Nc];
  diag.cols = new dlong[diag.nnz];
  diag.vals = new dfloat[diag.nnz];

  /*Write nonZeros into coarse graph and compress weights*/
  for (dlong c=0;c<Nc;++c) {
    entry_t* nonZerosn = nonZeros + nzStarts[c];
    const int Nentries = nzCounts[c];

    dlong*  colsc = diag.cols + diag.rowStarts[c];
    dfloat* valsc = diag.vals + diag.rowStarts[c];

    int cnt=0;
    if (Nentries>0) {
      colsc[0] = nonZerosn[0].col;
      valsc[0] = nonZerosn[0].val;
    }
    for (dlong j=1;j<Nentries;++j) {
      if (nonZerosn[j].col != nonZerosn[j-1].col) {
        cnt++;
        colsc[cnt] = nonZerosn[j].col;
        valsc[cnt] = nonZerosn[j].val;
      } else {
        valsc[cnt] += nonZerosn[j].val;
      }
    }
  }



  // for (dlong v=0;v<Nf;++v) {
  //   printf("FineToCoarse[%d] = %d, P[%d] = %f, null[%d] = %f\n", v, FtoC[v], v, cP[v], v, graph.null[v]);
  // }

  // for (dlong v=0;v<Nc;++v) {
  //   printf("cnull[%d] = %f\n", v, cgraph.null[v]);
  // }

  // for (dlong c=0;c<Nc;++c) {

  //   printf("Node %d, Dweight %f, Edges", c, cgraph.Dweight[c]);
  //   dlong start = cgraph.rowStarts[c];
  //   dlong end   = cgraph.rowStarts[c+1];
  //   for (dlong j=start;j<end;++j) {
  //     printf(" %d (%f),", cgraph.colIds[j], -cgraph.Eweight[j]);
  //   }
  //   printf("\n");
  // }

  // for (dlong c=0;c<Nc;++c) {
  //   dfloat a = cgraph.Dweight[c]*cgraph.null[c];
  //   dlong start = cgraph.rowStarts[c];
  //   dlong end   = cgraph.rowStarts[c+1];
  //   for (dlong j=start;j<end;++j) {
  //     a -= cgraph.Eweight[j]*cgraph.null[cgraph.colIds[j]];
  //   }
  //   printf("null test %d = %f \n", c, a);
  // }

  delete[] nonZeros;
  delete[] nzStarts;
  delete[] nzCounts;
}



}
