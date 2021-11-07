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

/****************************************/
/* Construct Multilevel Laplacian       */
/****************************************/
void graph_t::MultigridSetup() {

  /*Target size for coarsest graph*/
  const int coarseSize = 100;

  /*Coarsening tolerance. If a coarse graph isn't at least
    coarseTol times smaller than a fine graph, we consider the coarsening
    to be stalling*/
  const float coarseTol = 0.8;

  do{
    /*Get coarsest level*/
    mgLevel_t& Lf = L[Nlevels-1];

    /*If the graph is small enough, we're done*/
    if (Lf.A.Nrows <= coarseSize) {
      coarseSolver.Setup(Lf.A, Lf.null);
      break;
    }

    if (Nlevels>=PARADOGS_MAX_LEVELS)
      LIBP_ABORT("Paradogs: Max levels exceeded in coarse graph creation. Increase PARADOGS_MAX_LEVELS.")

    Lf.SetupSmoother();

    /*Construct next level via coarsening*/
    mgLevel_t& Lc = L[Nlevels];    
    Lc.CoarsenLevel(Lf);
    Nlevels++;
    
    /*Check for stalls*/
    if (Lc.A.Nrows > coarseTol*Lf.A.Nrows) {
      stringstream ss;
      ss << "Paradogs: Graph coarsening stalling. Coarse graph has " << Lc.A.Nrows << " nodes.";
      LIBP_WARNING(ss.str())
      coarseSolver.Setup(Lc.A, Lc.null);
      break;
    }
  } while(true);
}

void mgLevel_t::SetupSmoother() {
  A.diagInv = new dfloat[A.Nrows];

  for (dlong n=0;n<A.Nrows;++n) {
    A.diagInv[n] = 1.0/A.diagA[n];
  }

  // estimate rho(invD * A)
  dfloat rho = A.rhoDinvA(null);

  /*Smoothing params*/
  lambda1 = rho;
  lambda0 = rho/10.;

  /*Scratch space*/
  scratch = new dfloat[2*A.Ncols];
}

/*Coarsen a graph using an aggregation*/
void mgLevel_t::CoarsenLevel(mgLevel_t& Lf) {

  /*Create a FineToCoarse mapping*/
  const dlong Nf = Lf.A.Nrows;

  /*Create a vertex matching*/
  Lf.FineToCoarse = new dlong[Nf];
  Lf.A.Aggregate(Lf.Nc, Lf.FineToCoarse);
  dlong* FtoC = Lf.FineToCoarse;

  /* Tentative prolongation operator*/
  Lf.P = new dfloat[Nf];

  /* Each entry is the null vector entry*/
  for (dlong v=0;v<Nf;++v) Lf.P[v] = Lf.null[v];

  /*Create coarse null*/
  null = new dfloat[Lf.Nc];
  for (dlong v=0;v<Lf.Nc;++v) null[v] = 0.0;

  /*Sum columns of P*/
  //add local nonzeros
  for (dlong v=0;v<Nf;++v)
    null[FtoC[v]] += Lf.P[v] * Lf.P[v];

  //add nonlocal nonzeros
  // for(dlong i=0; i<P->offd.nnz; i++)
  //   null[P->offd.cols[i]] += P->offd.vals[i] * P->offd.vals[i];

  //add the halo values to their origins
  // P->halo->Combine(null, 1, ogs::Dfloat);

  for (dlong n=0;n<Lf.Nc;++n)
    null[n] = sqrt(null[n]);

  //share the results
  // P->halo->Exchange(null, 1, ogs::Dfloat);

  for(dlong v=0; v<Nf; v++)
    Lf.P[v] /= null[FtoC[v]];
  // for(dlong v=0; v<Nf; v++)
  //   P[v] /= cgraph.null[FtoC[v]];

  A.GalerkinProduct(Lf.A, Lf.Nc, FtoC, Lf.P);

  /*Space for Fiedler*/
  Fiedler = new dfloat[A.Nrows];

  /*Multigrid buffers*/
  RHS = new dfloat[A.Nrows];
  X = new dfloat[A.Nrows];
  RES = new dfloat[A.Nrows];
}

void coarseSolver_t::Setup(parCSR &A, dfloat null[]) {

  // comm = A->comm;
  // MPI_Comm_rank(comm, &rank);
  // MPI_Comm_size(comm, &size);

  //copy the global coarse partition as ints
  // coarseOffsets = (int* ) calloc(size+1,sizeof(int));
  // for (int r=0;r<size+1;r++) coarseOffsets[r] = (int) A->globalRowStarts[r];

  // coarseTotal   = coarseOffsets[size];
  // coarseOffset  = coarseOffsets[rank];

  N = (int) A.Nrows;
  Nrows = A.Nrows;
  Ncols = A.Ncols;

  coarseTotal=N;

  // coarseCounts = (int*) calloc(size,sizeof(int));

  // int sendNNZ = (int) (A->diag.nnz+A->offd.nnz);

  // if((rank==0)&&(settings.compareSetting("VERBOSE","TRUE")))
  //   {printf("Setting up coarse solver...");fflush(stdout);}

  // parCOO::nonZero_t *sendNonZeros = (parCOO::nonZero_t *) calloc(sendNNZ, sizeof(parCOO::nonZero_t));

  //populate matrix
  // int cnt = 0;
  // for (int n=0;n<N;n++) {
  //   const int start = (int) A->diag.rowStarts[n];
  //   const int end   = (int) A->diag.rowStarts[n+1];
  //   for (int m=start;m<end;m++) {
  //     sendNonZeros[cnt].row = n + coarseOffset;
  //     sendNonZeros[cnt].col = A->diag.cols[m] + coarseOffset;
  //     sendNonZeros[cnt].val = A->diag.vals[m];
  //     cnt++;
  //   }
  // }

  // for (int n=0;n<A->offd.nzRows;n++) {
  //   const int row   = (int) A->offd.rows[n];
  //   const int start = (int) A->offd.mRowStarts[n];
  //   const int end   = (int) A->offd.mRowStarts[n+1];
  //   for (int m=start;m<end;m++) {
  //     sendNonZeros[cnt].row = row + coarseOffset;
  //     sendNonZeros[cnt].col = A->colMap[A->offd.cols[m]];
  //     sendNonZeros[cnt].val = A->offd.vals[m];
  //     cnt++;
  //   }
  // }

  // //get the nonzero counts from all ranks
  // int *recvNNZ    = (int*) calloc(size,sizeof(int));
  // int *NNZoffsets = (int*) calloc(size+1,sizeof(int));
  // MPI_Allgather(&sendNNZ, 1, MPI_INT,
  //                recvNNZ, 1, MPI_INT, comm);

  // int totalNNZ = 0;
  // for (int r=0;r<size;r++) {
  //   totalNNZ += recvNNZ[r];
  //   NNZoffsets[r+1] = NNZoffsets[r] + recvNNZ[r];
  // }

  // parCOO::nonZero_t *recvNonZeros = (parCOO::nonZero_t *) calloc(totalNNZ, sizeof(parCOO::nonZero_t));

  // MPI_Allgatherv(sendNonZeros, sendNNZ,             MPI_NONZERO_T,
  //                recvNonZeros, recvNNZ, NNZoffsets, MPI_NONZERO_T, comm);

  // //gather null vector
  // dfloat *nullTotal = (dfloat*) calloc(coarseTotal,sizeof(dfloat));

  // for (int r=0;r<size;r++)
  //   coarseCounts[r] = coarseOffsets[r+1]-coarseOffsets[r];

  // MPI_Allgatherv(nullVector,            N,                MPI_DFLOAT,
  //                 nullTotal, coarseCounts, coarseOffsets, MPI_DFLOAT,
  //                comm);

  // //clean up
  // MPI_Barrier(comm);
  // free(sendNonZeros);
  // free(NNZoffsets);
  // free(recvNNZ);

  //assemble the full matrix
  dfloat *coarseA = new dfloat[coarseTotal*coarseTotal];
  // for (int i=0;i<totalNNZ;i++) {
  //   int n = recvNonZeros[i].row;
  //   int m = recvNonZeros[i].col;
  //   coarseA[n*coarseTotal+m] = recvNonZeros[i].val;
  // }

  // if (nullSpace) { //A is dense due to nullspace augmentation
  //   for (int n=0;n<coarseTotal;n++) {
  //     for (int m=0;m<coarseTotal;m++) {
  //       coarseA[n*coarseTotal+m] += nullSpacePenalty*nullTotal[n]*nullTotal[m];
  //     }
  //   }
  // }

  // free(recvNonZeros);
  // free(nullTotal);

  for (int n=0;n<coarseTotal;n++) {
    for (int m=0;m<coarseTotal;m++) {
      coarseA[n*coarseTotal+m] = 0.0;
    }
  }
  for (int n=0;n<A.Nrows;n++) {
    const dlong start = A.diag.rowStarts[n];
    const dlong end   = A.diag.rowStarts[n+1];
    for (int j=start;j<end;++j) {
      coarseA[n*N + A.diag.cols[j]] = A.diag.vals[j];
    }
  }

  /*Augment null space*/
  for (int n=0;n<coarseTotal;n++) {
    for (int m=0;m<coarseTotal;m++) {
      coarseA[n*coarseTotal+m] += null[n]*null[m];
    }
  }

  matrixInverse(coarseTotal, coarseA);

  //determine size of offd piece
  // offdTotal = coarseTotal - N;

  // //shift offsets for MPI_AllGatherv of offd pieces
  // for (int r=rank+1;r<=size;r++)
  //   coarseOffsets[r]-= N;

  // //dont copy the local piece in MPI_AllGatherv
  // coarseCounts[rank]=0;

  // //counts for all-to-all
  // sendCounts = (int* ) calloc(size,sizeof(int));
  // sendOffsets = (int* ) calloc(size,sizeof(int));
  // for (int r=0;r<size;r++) {
  //   sendCounts[r] = N;
  //   sendOffsets[r] = 0;
  // }
  // sendCounts[rank] = 0;

  // //diag piece of invA
  // diagInvAT = (dfloat *) calloc(N*N,sizeof(dfloat));
  // for (int n=0;n<N;n++) {
  //   for (int m=0;m<N;m++) {
  //     diagInvAT[n+m*N] = coarseA[(n+coarseOffset)*coarseTotal+(m+coarseOffset)];
  //   }
  // }

  // //offd piece of invA
  // offdInvAT = (dfloat *) calloc(N*offdTotal,sizeof(dfloat));
  // for (int n=0;n<N;n++) {
  //   for (int m=0;m<coarseOffset;m++) {
  //     offdInvAT[n+m*N] = coarseA[(n+coarseOffset)*coarseTotal+m];
  //   }
  //   for (int m=coarseOffset+N;m<coarseTotal;m++) {
  //     offdInvAT[n+(m-N)*N] = coarseA[(n+coarseOffset)*coarseTotal+m];
  //   }
  // }

  // o_diagInvAT = platform.malloc(N*N*sizeof(dfloat), diagInvAT);
  // o_offdInvAT = platform.malloc(N*offdTotal*sizeof(dfloat), offdInvAT);

  // diagRhs = (dfloat*) calloc(N,sizeof(dfloat));
  // offdRhs = (dfloat*) calloc(offdTotal,sizeof(dfloat));

  // o_offdRhs = platform.malloc(offdTotal*sizeof(dfloat));

  // if((rank==0)&&(settings.compareSetting("VERBOSE","TRUE"))) printf("done.\n");


  //diag piece of invA
  diagInvA = new dfloat[N*N];
  for (int n=0;n<N;n++) {
    for (int m=0;m<N;m++) {
      diagInvA[n*N+m] = coarseA[n*N+m];
    }
  }

  delete[] coarseA;
}

coarseSolver_t::~coarseSolver_t() {
  if(coarseOffsets) {delete[] coarseOffsets; coarseOffsets=nullptr;}
  if(coarseCounts)  {delete[] coarseCounts; coarseCounts=nullptr;}
  if(sendOffsets)   {delete[] sendOffsets; sendOffsets=nullptr;}
  if(sendCounts)    {delete[] sendCounts; sendCounts=nullptr;}
  if(diagInvA)      {delete[] diagInvA; diagInvA=nullptr; }
  if(offdInvA)      {delete[] offdInvA; offdInvA=nullptr;}
  if(diagRhs)       {delete[] diagRhs; diagRhs=nullptr; }
  if(offdRhs)       {delete[] offdRhs; offdRhs=nullptr;}
}

}


