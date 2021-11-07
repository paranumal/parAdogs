/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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
#include <random>

namespace paradogs {

extern std::mt19937 RNG;

//------------------------------------------------------------------------
//
//  parCSR matrix
//
//------------------------------------------------------------------------

void parCSR::SpMV(const dfloat alpha, dfloat *x,
                  const dfloat beta, dfloat *y) {

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  // #pragma omp parallel for
  for(dlong i=0; i<Nrows; i++){ //local
    dfloat result = 0.0;
    for(dlong jj=diag.rowStarts[i]; jj<diag.rowStarts[i+1]; jj++)
      result += diag.vals[jj]*x[diag.cols[jj]];

    if (beta!=0.0)
      y[i] = alpha*result + beta*y[i];
    else
      y[i] = alpha*result;
  }

  // halo->Exchange(x, 1, ogs::Dfloat);

  // // #pragma omp parallel for
  // for(dlong i=0; i<offd.nzRows; i++){ //local
  //   const dlong row = offd.rows[i];
  //   dfloat result = 0.0;
  //   for(dlong jj=offd.mRowStarts[i]; jj<offd.mRowStarts[i+1]; jj++)
  //     result += offd.vals[jj]*x[offd.cols[jj]];

  //   y[row] += alpha*result;
  // }
}

void parCSR::SpMV(const dfloat alpha, dfloat *x,
                  const dfloat beta, const dfloat *y, dfloat *z) {

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  // #pragma omp parallel for
  for(dlong i=0; i<Nrows; i++){ //local
    dfloat result = 0.0;
    for(dlong jj=diag.rowStarts[i]; jj<diag.rowStarts[i+1]; jj++)
      result += diag.vals[jj]*x[diag.cols[jj]];

    z[i] = alpha*result + beta*y[i];
  }

  // halo->Exchange(x, 1, ogs::Dfloat);

  // for(dlong i=0; i<offd.nzRows; i++){ //local
  //   const dlong row = offd.rows[i];
  //   dfloat result = 0.0;
  //   for(dlong jj=offd.mRowStarts[i]; jj<offd.mRowStarts[i+1]; jj++)
  //     result += offd.vals[jj]*x[offd.cols[jj]];

  //   z[row] += alpha*result;
  // }
}

//------------------------------------------------------------------------
//
//  parCSR matrix setup
//
//------------------------------------------------------------------------

//build a parCSR matrix from a distributed COO matrix
parCSR::parCSR(parCOO& A)       // number of nonzeros on this rank
  // platform(A.platform),
  // comm(A.comm) 
  {

  // int rank;
  // int size;
  // MPI_Comm_rank(comm, &rank);
  // MPI_Comm_size(comm, &size);

  // //copy global partition
  // globalRowStarts = (hlong *) calloc(size+1,sizeof(hlong));
  // globalColStarts = (hlong *) calloc(size+1,sizeof(hlong));
  // memcpy(globalRowStarts, A.globalRowStarts, (size+1)*sizeof(hlong));
  // memcpy(globalColStarts, A.globalColStarts, (size+1)*sizeof(hlong));

  // const hlong globalRowOffset = globalRowStarts[rank];
  // const hlong globalColOffset = globalColStarts[rank];

  // Nrows = (dlong)(globalRowStarts[rank+1]-globalRowStarts[rank]);
  // Ncols = (dlong)(globalColStarts[rank+1]-globalColStarts[rank]);

  // diag.rowStarts = (dlong *) calloc(Nrows+1, sizeof(dlong));
  // offd.rowStarts = (dlong *) calloc(Nrows+1, sizeof(dlong));

  // //count the entries in each row
  // for (dlong n=0;n<A.nnz;n++) {
  //   const dlong row = (dlong) (A.entries[n].row - globalRowOffset);
  //   if (   (A.entries[n].col < globalColOffset)
  //       || (A.entries[n].col > globalColOffset+Ncols-1))
  //     offd.rowStarts[row+1]++;
  //   else
  //     diag.rowStarts[row+1]++;
  // }

  // offd.nzRows=0;

  // // count how many rows are shared
  // for(dlong i=0; i<Nrows; i++)
  //   if (offd.rowStarts[i+1]>0) offd.nzRows++;

  // offd.rows       = (dlong *) calloc(offd.nzRows, sizeof(dlong));
  // offd.mRowStarts = (dlong *) calloc(offd.nzRows+1, sizeof(dlong));

  // // cumulative sum
  // dlong cnt=0;
  // for(dlong i=0; i<Nrows; i++) {
  //   if (offd.rowStarts[i+1]>0) {
  //     offd.rows[cnt] = i; //record row id
  //     offd.mRowStarts[cnt+1] = offd.mRowStarts[cnt] + offd.rowStarts[i+1];
  //     cnt++;
  //   }
  //   diag.rowStarts[i+1] += diag.rowStarts[i];
  //   offd.rowStarts[i+1] += offd.rowStarts[i];
  // }
  // diag.nnz = diag.rowStarts[Nrows];
  // offd.nnz = offd.rowStarts[Nrows];

  // // Halo setup
  // cnt=0;
  // hlong *colIds = (hlong *) malloc(offd.nnz*sizeof(hlong));
  // for (dlong n=0;n<A.nnz;n++) {
  //   if ( (A.entries[n].col < globalColOffset)
  //     || (A.entries[n].col > globalColOffset+Ncols-1))
  //     colIds[cnt++] = A.entries[n].col;
  // }
  // haloSetup(colIds); //setup halo, and transform colIds to a local indexing

  // //fill the CSR matrices
  // diag.cols = (dlong *)  calloc(diag.nnz, sizeof(dlong));
  // offd.cols = (dlong *)  calloc(offd.nnz, sizeof(dlong));
  // diag.vals = (pfloat *) calloc(diag.nnz, sizeof(pfloat));
  // offd.vals = (pfloat *) calloc(offd.nnz, sizeof(pfloat));
  // dlong diagCnt = 0;
  // dlong offdCnt = 0;
  // for (dlong n=0;n<A.nnz;n++) {
  //   if ( (A.entries[n].col < globalColOffset)
  //     || (A.entries[n].col > globalColOffset+NlocalCols-1)) {
  //     offd.cols[offdCnt] = colIds[offdCnt];
  //     offd.vals[offdCnt] = A.entries[n].val;
  //     offdCnt++;
  //   } else {
  //     diag.cols[diagCnt] = (dlong) (A.entries[n].col - globalColOffset);
  //     diag.vals[diagCnt] = A.entries[n].val;
  //     diagCnt++;
  //   }
  // }
  // free(colIds);
}

//------------------------------------------------------------------------
//
//  parCSR halo setup
//
//------------------------------------------------------------------------

typedef struct {

  dlong localId;
  hlong globalId;

  dlong newId;

} parallelId_t;


void parCSR::haloSetup(hlong *colIds) {

  // int rank;
  // MPI_Comm_rank(comm, &rank);

  // const hlong globalOffset = globalColStarts[rank];

  // //collect the unique nonlocal column ids
  // parallelId_t*  parIds = (parallelId_t*) malloc(offd.nnz*sizeof(parallelId_t));

  // for (dlong n=0;n<offd.nnz;n++) {
  //   parIds[n].localId  = n;
  //   parIds[n].globalId = colIds[n];
  // }

  // //sort by global index
  // std::sort(parIds, parIds+offd.nnz,
  //           [](const parallelId_t& a, const parallelId_t& b) {
  //             if(a.globalId < b.globalId) return true;
  //             if(a.globalId > b.globalId) return false;

  //             return (a.localId < b.localId);
  //           });

  // //count unique nonlocal column ids
  // dlong Noffdcols = 0; //number of unique columns
  // if(offd.nnz) parIds[0].newId = Noffdcols;
  // for (dlong n=1;n<offd.nnz;n++) {
  //   if (parIds[n].globalId != parIds[n-1].globalId)
  //     Noffdcols++;

  //   parIds[n].newId = Noffdcols;
  // }
  // if(offd.nnz) Noffdcols++;

  // //record the global ids of the unique columns
  // hlong *offdcols = (hlong *) malloc(Noffdcols*sizeof(hlong));
  // Noffdcols = 0;
  // if(offd.nnz) offdcols[Noffdcols++] = parIds[0].globalId;
  // for (dlong n=1;n<offd.nnz;n++)
  //   if (parIds[n].globalId != parIds[n-1].globalId)
  //     offdcols[Noffdcols++] = parIds[n].globalId;

  // //sort back to local order
  // std::sort(parIds, parIds+offd.nnz,
  //           [](const parallelId_t& a, const parallelId_t& b) {
  //             if(a.localId < b.localId) return true;
  //             if(a.localId > b.localId) return false;

  //             return (a.globalId < b.globalId);
  //           });

  // // be careful to make sure Ncols is set at this point
  // NlocalCols = Ncols;
  // Ncols += Noffdcols;

  // //make an array of all the column ids required on this rank (local first)
  // colMap = (hlong*) malloc(Ncols*sizeof(hlong));
  // for (dlong n=0; n<NlocalCols; n++)      colMap[n] = n+globalOffset+1; //local rows
  // for (dlong n=NlocalCols; n<Ncols; n++)  colMap[n] = -(offdcols[n-NlocalCols]+1);    //nonlocal rows

  // //make a halo exchange to share column entries and an ogs for gsops accross columns
  // bool verbose = false;
  // halo = new ogs::halo_t(platform);
  // halo->Setup(Ncols, colMap, comm, ogs::Auto, verbose);

  // //shift back to 0-indexed
  // for (dlong n=0; n<Ncols; n++) colMap[n]=abs(colMap[n])-1;

  // //update column numbering
  // for (dlong n=0;n<offd.nnz;n++)
  //   colIds[n] = NlocalCols + parIds[n].newId;

  // free(parIds);
}

parCSR::~parCSR() {
  if (diag.rowStarts) {delete[] diag.rowStarts; diag.rowStarts=nullptr;}
  if (diag.cols)      {delete[] diag.cols; diag.cols=nullptr;}
  if (diag.vals)      {delete[] diag.vals; diag.vals=nullptr;}

  if (offd.rowStarts)  {delete[] offd.rowStarts; offd.rowStarts=nullptr;}
  if (offd.mRowStarts) {delete[] offd.mRowStarts; offd.mRowStarts=nullptr;}
  if (offd.rows) {delete[] offd.rows; offd.rows=nullptr;}
  if (offd.cols) {delete[] offd.cols; offd.cols=nullptr;}
  if (offd.vals) {delete[] offd.vals; offd.vals=nullptr;}

  if (diagA)   {delete[] diagA; diagA=nullptr;}
  if (diagInv) {delete[] diagInv; diagInv=nullptr;}

  if (globalRowStarts) {delete[] globalRowStarts; globalRowStarts=nullptr;}
  if (globalColStarts) {delete[] globalColStarts; globalColStarts=nullptr;}
  if (colMap) {delete[] colMap; colMap=nullptr;}

  if (halo)   {halo->Free(); halo = nullptr;}
}

//------------------------------------------------------------------------
//
//  parCSR Estimate max Eigenvalue of diagA^{-1}*A
//
//------------------------------------------------------------------------

dfloat parCSR::rhoDinvA(dfloat null[]){

  // int size;
  // MPI_Comm_size(comm, &size);

  int k = 10;

  // hlong Ntotal = globalRowStarts[size];
  hlong Ntotal = Nrows;
  if(k > Ntotal) k = (int) Ntotal;

  // do an arnoldi

  // allocate memory for Hessenberg matrix
  double *H = new double[k*k];
  for(int n=0; n<k*k; n++) H[n] = 0.0;

  // allocate memory for basis
  dfloat *V = new dfloat[(k+1)*Nrows];
  dfloat *Vx = new dfloat[Ncols];

  /*Create rng*/
  std::uniform_real_distribution<dfloat> distrib(-0.5, 0.5);

  // generate a random vector for initial basis vector
  for(dlong n=0; n<Nrows; n++) Vx[n] = distrib(RNG);

  /*Project out null vector*/
  dfloat nulldot =0.0;
  for(dlong n=0; n<Nrows; n++) nulldot += null[n]*Vx[n];
  for(dlong n=0; n<Nrows; n++) Vx[n] -= nulldot*null[n];

  // dfloat norm_vo = vectorNorm(Nrows,Vx, comm);
  dfloat norm_vo=0.0, gnorm_vo=0.0;
  for(dlong n=0; n<Nrows; n++) norm_vo += Vx[n]*Vx[n];
  // MPI_Allreduce(&norm_vo, &gnorm_vo, 1, MPI_DFLOAT, MPI_SUM, comm);
  gnorm_vo = norm_vo;
  norm_vo = sqrt(gnorm_vo);

  // vectorScale(Nrows, 1.0/norm_vo, Vx);
  for(dlong n=0; n<Nrows; n++) Vx[n] *= (1.0/norm_vo);

  //V[0] = Vx
  std::copy(Vx, Vx+Nrows, V);

  for(int j=0; j<k; j++){
    dfloat *Vj   = V+j*Nrows;
    dfloat *Vjp1 = V+(j+1)*Nrows;

    //Vx = V[j]
    std::copy(Vj, Vj+Nrows, Vx);

    // v[j+1] = invD*(A*v[j])
    SpMV(1.0, Vx, 0., Vjp1);
    // vectorDotStar(Nrows, diagInv, V[j+1]);
    for(dlong n=0; n<Nrows; n++) Vjp1[n] *= diagInv[n];

    // modified Gram-Schmidth
    for(int i=0; i<=j; i++){
      dfloat *Vi = V+i*Nrows;
      // H(i,j) = v[i]'*A*v[j]
      // dfloat hij = vectorInnerProd(Nrows, V[i], V[j+1],comm);
      dfloat local_hij=0.0, hij=0.0;
      for(dlong n=0; n<Nrows; n++) local_hij += Vi[n]*Vjp1[n];
      // MPI_Allreduce(&local_hij, &hij, 1, MPI_DFLOAT, MPI_SUM, comm);
      hij = local_hij;

      // v[j+1] = v[j+1] - hij*v[i]
      // vectorAdd(Nrows,-hij, V[i], 1.0, V[j+1]);
      for(dlong n=0; n<Nrows; n++) Vjp1[n] += -hij*Vi[n];

      H[i + j*k] = (double) hij;
    }

    if(j+1 < k){

      // dfloat norm_vj = vectorNorm(Nrows,V[j+1],comm);
      dfloat norm_vj=0.0, gnorm_vj=0.0;
      for(dlong n=0; n<Nrows; n++) norm_vj += Vjp1[n]*Vjp1[n];
      // MPI_Allreduce(&norm_vj, &gnorm_vj, 1, MPI_DFLOAT, MPI_SUM, comm);
      gnorm_vj = norm_vj;
      norm_vj = sqrt(gnorm_vj);

      H[j+1+ j*k] = (double) norm_vj;

      // vectorScale(Nrows, 1./H[j+1 + j*k], V[j+1]);
      for(dlong n=0; n<Nrows; n++) Vjp1[n] *= (1./H[j+1 + j*k]);
    }
  }

  double *WR = new double[k];
  double *WI = new double[k];

  matrixEigenValues(k, H, WR, WI);

  double RHO = 0.;

  for(int i=0; i<k; i++){
    double RHO_i  = sqrt(WR[i]*WR[i] + WI[i]*WI[i]);

    if(RHO < RHO_i) {
      RHO = RHO_i;
    }
  }

  // free memory
  delete[] H;
  delete[] WR;
  delete[] WI;

  delete[] Vx;
  delete[] V;

  // printf("weight = %g \n", RHO);

  return RHO;
}


} //namespace paradogs
