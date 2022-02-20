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

namespace libp {

namespace paradogs {

parCSR* SpMM(const parCSR* A, const parCSR* B){

  // MPI info
  int rank, size;
  MPI_Comm_rank(A->comm, &rank);
  MPI_Comm_size(A->comm, &size);

  // To compute C = A*B we need all the rows B(j,:) for which
  // j is a column index for the nonzeros of A on this rank.
  // For all local column indices in A->diag, we will already
  // have the row of B on this rank, so we just need to gather
  // the offd colIds

  hlong *recvRows = new hlong[A->Ncols-A->NlocalCols];
  int *sendCounts = new int[size];
  int *recvCounts = new int[size];
  int *sendOffsets = new int[size+1];
  int *recvOffsets = new int[size+1];

  hlong *globalRowStarts = new hlong[size+1];
  globalRowStarts[0]=0;
  MPI_Allgather(&(B->rowOffsetU), 1, MPI_HLONG,
                 globalRowStarts+1, 1, MPI_HLONG, B->comm);

  for (int r=0;r<size;++r) {
    recvCounts[r]=0;
  }
  //use the colMap of A to list the needed rows of B
  int r=0;
  for (dlong n=A->NlocalCols;n<A->Ncols;n++) {
    const hlong id = A->colMap[n];
    while (id>=globalRowStarts[r+1]) r++; //assumes the halo is sorted
    recvCounts[r]++;
    recvRows[n-A->NlocalCols] = id; //record the row to recv
  }
  delete[] globalRowStarts;

  //share the counts
  MPI_Alltoall(recvCounts, 1, MPI_INT,
               sendCounts, 1, MPI_INT, A->comm);

  sendOffsets[0]=0;
  recvOffsets[0]=0;
  for (r=0;r<size;r++) {
    sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  }

  int sendTotal = sendOffsets[size];
  hlong *sendRows = new hlong[sendTotal];

  //share the rowIds
  MPI_Alltoallv(recvRows, recvCounts, recvOffsets, MPI_HLONG,
                sendRows, sendCounts, sendOffsets, MPI_HLONG,
                B->comm);

  //we now have a list of rows to send, count the nnz to send
  dlong NNZ=0;
  for (r=0;r<size;r++) {
    sendCounts[r] =0; //reset
    for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      const dlong i = static_cast<dlong>(sendRows[n]-B->rowOffsetL); //local row id
      sendCounts[r]+= B->diag.rowStarts[i+1]-B->diag.rowStarts[i]; //count entries in this row
      sendCounts[r]+= B->offd.rowStarts[i+1]-B->offd.rowStarts[i]; //count entries in this row
    }
    NNZ += sendCounts[r]; //tally the total
  }

  nonZero_t *sendNonZeros = new nonZero_t[NNZ];

  NNZ=0; //reset
  for (r=0;r<size;r++) {
    for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      const dlong i = static_cast<dlong>(sendRows[n] - B->rowOffsetL); //local row id
      for (dlong jj=B->diag.rowStarts[i]; jj<B->diag.rowStarts[i+1];jj++){
        sendNonZeros[NNZ].row = sendRows[n];
        sendNonZeros[NNZ].col = B->diag.cols[jj] + B->colOffsetL;
        sendNonZeros[NNZ].val = B->diag.vals[jj];
        NNZ++;
      }
      for (dlong jj=B->offd.rowStarts[i]; jj<B->offd.rowStarts[i+1];jj++){
        sendNonZeros[NNZ].row = sendRows[n];
        sendNonZeros[NNZ].col = B->colMap[B->offd.cols[jj]];
        sendNonZeros[NNZ].val = B->offd.vals[jj];
        NNZ++;
      }
    }
  }

  MPI_Alltoall(sendCounts, 1, MPI_INT,
               recvCounts, 1, MPI_INT, A->comm);

  for (r=0;r<size;r++) {
    sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  }


  dlong Boffdnnz = recvOffsets[size]; //total nonzeros
  nonZero_t *BoffdRows = new nonZero_t[Boffdnnz];

  // Make the MPI_NONZERO_T data type
  MPI_Datatype MPI_NONZERO_T;
  MPI_Datatype dtype[3] = {MPI_HLONG, MPI_HLONG, MPI_DFLOAT};
  int blength[3] = {1, 1, 1};
  MPI_Aint addr[3], displ[3];
  MPI_Get_address ( &(BoffdRows[0]    ), addr+0);
  MPI_Get_address ( &(BoffdRows[0].col), addr+1);
  MPI_Get_address ( &(BoffdRows[0].val), addr+2);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  MPI_Type_create_struct (3, blength, displ, dtype, &MPI_NONZERO_T);
  MPI_Type_commit (&MPI_NONZERO_T);

  MPI_Alltoallv(sendNonZeros, sendCounts, sendOffsets, MPI_NONZERO_T,
                BoffdRows, recvCounts, recvOffsets, MPI_NONZERO_T,
                B->comm);

  //clean up
  MPI_Barrier(B->comm);
  MPI_Type_free(&MPI_NONZERO_T);
  delete[] sendNonZeros;
  delete[] sendRows;
  delete[] recvRows;
  delete[] sendCounts;
  delete[] recvCounts;
  delete[] sendOffsets;
  delete[] recvOffsets;

  //we now have all the needed nonlocal rows (should also be sorted by row then col)

  //make an array of row offsets so we know how large each row is
  dlong *BoffdRowOffsets = new dlong[A->Ncols-A->NlocalCols+1];

  for (dlong n=0;n<A->Ncols-A->NlocalCols+1;n++) {
    BoffdRowOffsets[n]=0;
  }

  dlong id=0;
  for (dlong n=0;n<Boffdnnz;n++) {
    hlong row = BoffdRows[n].row;

    while(A->colMap[id+A->NlocalCols]!=row) id++;

    BoffdRowOffsets[id+1]++; //count entry in row
  }

  //cumulative sum
  for (dlong n=0;n<A->Ncols-A->NlocalCols;n++)
    BoffdRowOffsets[n+1] += BoffdRowOffsets[n];


  // The next step to compute C = A*B is to multiply each entry A(i,j) by the
  // row B(j,:), store the all the results, sort them by row+col, and compress
  // the entries

  // Find how big the intermediate form is
  dlong *rowStarts = new dlong[A->Nrows+1];
  dlong *rowCounts = new dlong[A->Nrows];

  #pragma omp parallel for
  for(dlong i=0; i<A->Nrows+1; i++) rowStarts[i]=0;

  #pragma omp parallel for
  for(dlong i=0; i<A->Nrows; i++) rowCounts[i]=0;

  /*Count entries per row*/
  #pragma omp parallel for
  for (dlong i=0;i<A->Nrows;i++) {
    //local entries
    dlong start = A->diag.rowStarts[i];
    dlong end   = A->diag.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A->diag.cols[j];
      rowStarts[i+1] +=  B->diag.rowStarts[col+1]-B->diag.rowStarts[col]
                        +B->offd.rowStarts[col+1]-B->offd.rowStarts[col];
    }
    //non-local entries
    start = A->offd.rowStarts[i];
    end   = A->offd.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A->offd.cols[j]-A->NlocalCols;
      rowStarts[i+1] += BoffdRowOffsets[col+1] - BoffdRowOffsets[col];
    }
  }

  /*Cumulative sum*/
  for(dlong i=1; i<A->Nrows+1; i++) {
    rowStarts[i] += rowStarts[i-1];
  }

  NNZ = rowStarts[A->Nrows];

  nonZero_t *Ctmp = new nonZero_t[NNZ];

  //count total number of nonzeros;
  dlong nnz =0;

  // Fill the intermediate form of C
  // #pragma omp parallel for
  for (dlong i=0;i<A->Nrows;i++) {
    const dlong cStart = rowStarts[i];
    dlong& c = rowCounts[i];

    //local A entries
    dlong start = A->diag.rowStarts[i];
    dlong end   = A->diag.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A->diag.cols[j];
      const dfloat Aval = A->diag.vals[j];

      //local B entries
      dlong Bstart = B->diag.rowStarts[col];
      dlong Bend   = B->diag.rowStarts[col+1];
      for (dlong jj=Bstart;jj<Bend;jj++) {
        Ctmp[cStart+c].row = i + A->rowOffsetL;
        Ctmp[cStart+c].col = B->diag.cols[jj] + B->colOffsetL; //global id
        Ctmp[cStart+c].val = Aval*B->diag.vals[jj];
        c++;
      }
      //non-local B entries
      Bstart = B->offd.rowStarts[col];
      Bend   = B->offd.rowStarts[col+1];
      for (dlong jj=Bstart;jj<Bend;jj++) {
        Ctmp[cStart+c].row = i + A->rowOffsetL;
        Ctmp[cStart+c].col = B->colMap[B->offd.cols[jj]]; //global id
        Ctmp[cStart+c].val = Aval*B->offd.vals[jj];
        c++;
      }
    }
    //non-local A entries
    start = A->offd.rowStarts[i];
    end   = A->offd.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A->offd.cols[j]-A->NlocalCols;
      const dfloat Aval = A->offd.vals[j];

      // entries from recived rows of B
      dlong Bstart = BoffdRowOffsets[col];
      dlong Bend   = BoffdRowOffsets[col+1];
      for (dlong jj=Bstart;jj<Bend;jj++) {
        Ctmp[cStart+c].row = i + A->rowOffsetL;
        Ctmp[cStart+c].col = BoffdRows[jj].col; //global id
        Ctmp[cStart+c].val = Aval*BoffdRows[jj].val;
        c++;
      }
    }

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
  delete[] BoffdRowOffsets;
  delete[] BoffdRows;

  delete[] rowStarts;
  delete[] rowCounts;

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
  parCSR* C = new parCSR(A->Nrows, B->NlocalCols,
                         nnz, entries,
                         A->platform, A->comm);

  delete[] entries;

  return C;
}

} //namespace paradogs

} //namespace libp
