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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT-> IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "parAdogs.hpp"
#include "parAdogs/parAdogsGraph.hpp"

namespace libp {

namespace paradogs {

parCSR* SmoothProlongator(const parCSR* A, const parCSR* T) {

  // MPI info
  int rank, size;
  MPI_Comm_rank(A->comm, &rank);
  MPI_Comm_size(A->comm, &size);

  // This function computes a smoothed prologation operator
  // via a single weighted Jacobi iteration on the tentative
  // prologator, i.e.,
  //
  //   P = (I - omega*D^{-1}*A)*T
  //
  // To compute D^{-1}*A*T we need all the rows T(j,:) for which
  // j is a column index for the nonzeros of A on this rank.
  // For all local column indices in A->diag, we will already
  // have the row of T on this rank, so we just need to gather
  // the offd colIds

  //Jacobi weight
  const dfloat omega = (4./3.)/A->rho;

  hlong *recvRows = new hlong[A->Ncols-A->NlocalCols];
  int *sendCounts = new int[size];
  int *recvCounts = new int[size];
  int *sendOffsets = new int[size+1];
  int *recvOffsets = new int[size+1];

  hlong *globalRowStarts = new hlong[size+1];
  globalRowStarts[0]=0;
  MPI_Allgather(&(T->rowOffsetU), 1, MPI_HLONG,
                globalRowStarts+1, 1, MPI_HLONG, T->comm);

  for (int r=0;r<size;++r) {
    recvCounts[r]=0;
  }
  //use the colMap of A to list the needed rows of T
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
                T->comm);

  //we now have a list of rows to send, count the nnz to send
  dlong nnzTotal=0;
  for (r=0;r<size;r++) {
    sendCounts[r] =0; //reset
    for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      const dlong i = static_cast<dlong>(sendRows[n]-T->rowOffsetL); //local row id
      sendCounts[r]+= T->diag.rowStarts[i+1]-T->diag.rowStarts[i]; //count entries in this row
      sendCounts[r]+= T->offd.rowStarts[i+1]-T->offd.rowStarts[i]; //count entries in this row
    }
    nnzTotal += sendCounts[r]; //tally the total
  }

  nonZero_t *sendNonZeros = new nonZero_t[nnzTotal];

  nnzTotal=0; //reset
  for (r=0;r<size;r++) {
    for (int n=sendOffsets[r];n<sendOffsets[r+1];n++) {
      const dlong i = static_cast<dlong>(sendRows[n] - T->rowOffsetL); //local row id
      for (dlong jj=T->diag.rowStarts[i]; jj<T->diag.rowStarts[i+1];jj++){
        sendNonZeros[nnzTotal].row = sendRows[n];
        sendNonZeros[nnzTotal].col = T->diag.cols[jj] + T->colOffsetL;
        sendNonZeros[nnzTotal].val = T->diag.vals[jj];
        nnzTotal++;
      }
      for (dlong jj=T->offd.rowStarts[i]; jj<T->offd.rowStarts[i+1];jj++){
        sendNonZeros[nnzTotal].row = sendRows[n];
        sendNonZeros[nnzTotal].col = T->colMap[T->offd.cols[jj]];
        sendNonZeros[nnzTotal].val = T->offd.vals[jj];
        nnzTotal++;
      }
    }
  }

  MPI_Alltoall(sendCounts, 1, MPI_INT,
               recvCounts, 1, MPI_INT, A->comm);

  for (r=0;r<size;r++) {
    sendOffsets[r+1] = sendOffsets[r]+sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r]+recvCounts[r];
  }


  dlong Toffdnnz = recvOffsets[size]; //total nonzeros
  nonZero_t *ToffdRows = new nonZero_t[Toffdnnz];

  // Make the MPI_NONZERO_T data type
  MPI_Datatype MPI_NONZERO_T;
  MPI_Datatype dtype[3] = {MPI_HLONG, MPI_HLONG, MPI_DFLOAT};
  int blength[3] = {1, 1, 1};
  MPI_Aint addr[3], displ[3];
  MPI_Get_address ( &(ToffdRows[0]    ), addr+0);
  MPI_Get_address ( &(ToffdRows[0].col), addr+1);
  MPI_Get_address ( &(ToffdRows[0].val), addr+2);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  MPI_Type_create_struct (3, blength, displ, dtype, &MPI_NONZERO_T);
  MPI_Type_commit (&MPI_NONZERO_T);

  MPI_Alltoallv(sendNonZeros, sendCounts, sendOffsets, MPI_NONZERO_T,
                ToffdRows, recvCounts, recvOffsets, MPI_NONZERO_T,
                T->comm);

  //clean up
  MPI_Barrier(T->comm);
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
  dlong *ToffdRowOffsets = new dlong[A->Ncols-A->NlocalCols+1];

  for (dlong n=0;n<A->Ncols-A->NlocalCols+1;n++) {
    ToffdRowOffsets[n]=0;
  }

  dlong id=0;
  for (dlong n=0;n<Toffdnnz;n++) {
    hlong row = ToffdRows[n].row;

    while(A->colMap[id+A->NlocalCols]!=row) id++;

    ToffdRowOffsets[id+1]++; //count entry in row
  }

  //cumulative sum
  for (dlong n=0;n<A->Ncols-A->NlocalCols;n++)
    ToffdRowOffsets[n+1] += ToffdRowOffsets[n];


  // The next step to compute D^{-1}*A*T is to multiply each entry A(i,j) by the
  // row T(j,:), store the all the results, sort them by row+col, and compress
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
  for(dlong i=0; i<A->Nrows; i++) {
    /*Start with entries for T*/
    rowStarts[i+1]+=T->diag.rowStarts[i+1]-T->diag.rowStarts[i] +
                    T->offd.rowStarts[i+1]-T->offd.rowStarts[i];

    /*Then add entries from A*T*/
    dlong Jstart = A->diag.rowStarts[i];
    dlong Jend   = A->diag.rowStarts[i+1];
    for(dlong jj=Jstart; jj<Jend; jj++){
      const dlong col = A->diag.cols[jj];
      rowStarts[i+1]+=T->diag.rowStarts[col+1]-T->diag.rowStarts[col] +
                      T->offd.rowStarts[col+1]-T->offd.rowStarts[col];
    }
    //non-local entries
    Jstart = A->offd.rowStarts[i];
    Jend   = A->offd.rowStarts[i+1];
    for (dlong jj=Jstart;jj<Jend;jj++) {
      const dlong col = A->offd.cols[jj]-A->NlocalCols;
      rowStarts[i+1]+= ToffdRowOffsets[col+1] - ToffdRowOffsets[col];
    }
  }

  /*Cumulative sum*/
  for(dlong i=1; i<A->Nrows+1; i++) {
    rowStarts[i] += rowStarts[i-1];
  }

  dlong NNZ = rowStarts[A->Nrows];

  nonZero_t *Ptmp = new nonZero_t[NNZ];

  //count total number of nonzeros we find
  dlong nnz =0;

  // Fill the intermediate form of P
  // #pragma omp parallel for
  for (dlong i=0;i<A->Nrows;i++) {
    const dlong cStart = rowStarts[i];
    dlong& c = rowCounts[i];

    /*Start with P=T entries*/

    //local T entries
    dlong start = T->diag.rowStarts[i];
    dlong end   = T->diag.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      Ptmp[cStart+c].row = i + T->rowOffsetL;
      Ptmp[cStart+c].col = T->diag.cols[j] + T->colOffsetL; //global id
      Ptmp[cStart+c].val = T->diag.vals[j];
      c++;
    }
    //non-local T entries
    start = T->offd.rowStarts[i];
    end   = T->offd.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      Ptmp[cStart+c].row = i + T->rowOffsetL;
      Ptmp[cStart+c].col = T->colMap[T->offd.cols[j]];
      Ptmp[cStart+c].val = T->offd.vals[j];
      c++;
    }

    /*Then P -= omega*invD*A*T*/

    //local A entries
    start = A->diag.rowStarts[i];
    end   = A->diag.rowStarts[i+1];

    const dfloat invDi = 1.0/A->diagA[i];

    for (dlong j=start;j<end;j++) {
      const dlong col = A->diag.cols[j];
      const dfloat Aval = -omega*invDi*A->diag.vals[j];

      //local T entries
      dlong Tstart = T->diag.rowStarts[col];
      dlong Tend   = T->diag.rowStarts[col+1];
      for (dlong jj=Tstart;jj<Tend;jj++) {
        Ptmp[cStart+c].row = i + A->rowOffsetL;
        Ptmp[cStart+c].col = T->diag.cols[jj] + T->colOffsetL; //global id
        Ptmp[cStart+c].val = Aval*T->diag.vals[jj];
        c++;
      }
      //non-local T entries
      Tstart = T->offd.rowStarts[col];
      Tend   = T->offd.rowStarts[col+1];
      for (dlong jj=Tstart;jj<Tend;jj++) {
        Ptmp[cStart+c].row = i + A->rowOffsetL;
        Ptmp[cStart+c].col = T->colMap[T->offd.cols[jj]]; //global id
        Ptmp[cStart+c].val = Aval*T->offd.vals[jj];
        c++;
      }
    }
    //non-local A entries
    start = A->offd.rowStarts[i];
    end   = A->offd.rowStarts[i+1];
    for (dlong j=start;j<end;j++) {
      const dlong col = A->offd.cols[j]-A->NlocalCols;
      const dfloat Aval = -omega*invDi*A->offd.vals[j];

      // entries from recived rows of T
      dlong Tstart = ToffdRowOffsets[col];
      dlong Tend   = ToffdRowOffsets[col+1];
      for (dlong jj=Tstart;jj<Tend;jj++) {
        Ptmp[cStart+c].row = i + A->rowOffsetL;
        Ptmp[cStart+c].col = ToffdRows[jj].col; //global id
        Ptmp[cStart+c].val = Aval*ToffdRows[jj].val;
        c++;
      }
    }

    //sort entries in this row by col id
    std::sort(Ptmp+cStart, Ptmp+cStart+c,
              [](const nonZero_t& a, const nonZero_t& b) {
                return a.col < b.col;
              });

    /*Count how many actual nonzeros will be in this row*/
    dlong nnzRow=0;
    if (c>0) nnzRow++;
    for (dlong j=1;j<c;j++) {
      if ((Ptmp[cStart+j].col!=Ptmp[cStart+j-1].col)) nnzRow++;
    }

    nnz+=nnzRow; //Add to total
  }
  delete[] ToffdRowOffsets;
  delete[] ToffdRows;

  delete[] rowStarts;
  delete[] rowCounts;

  // cooP.nnz = nnz;
  nonZero_t *entries = new nonZero_t[nnz];

  //compress nonzeros
  nnz = 0;
  if (NNZ) entries[nnz++] = Ptmp[0];
  for (dlong i=1;i<NNZ;i++) {
    if ((Ptmp[i].row!=Ptmp[i-1].row)||
        (Ptmp[i].col!=Ptmp[i-1].col)) {
      entries[nnz++] = Ptmp[i];
    } else {
      entries[nnz-1].val += Ptmp[i].val;
    }
  }
  //clean up
  delete[] Ptmp;

  //build P from coo matrix
  parCSR* P = new parCSR(A->Nrows, T->NlocalCols,
                         nnz, entries,
                         A->platform, A->comm);

  delete[] entries;

  return P;
}

} //namespace paradogs

} //namespace libp
