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

/*Create a vertex matching using distance-2 aggregation*/
void parCSR::Aggregate(dlong& Nc,
                       dlong FineToCoarse[]) {

  /*Stength threashold*/
  const dfloat theta=0.08;

  /*Create rng*/
  std::uniform_real_distribution<> distrib(-0.25, 0.25);

  parCSR strong(Nrows, Ncols);
  strong.diag.rowStarts = new dlong[Nrows+1];

  #pragma omp parallel for
  for(dlong i=0; i<Nrows+1; i++) {
    strong.diag.rowStarts[i]=0;
  }

  #pragma omp parallel for
  for(dlong i=0; i<Nrows; i++){
    int strong_per_row = 0;

    const dfloat Aii = diagA[i];

    //local entries
    dlong Jstart = diag.rowStarts[i];
    dlong Jend   = diag.rowStarts[i+1];
    for(dlong jj= Jstart; jj<Jend; jj++){
      const dlong col = diag.cols[jj];
      if (col==i) continue;

      const dfloat Ajj = diagA[col];

      if(-diag.vals[jj] > theta*(sqrt(Aii*Ajj)))
        strong_per_row++;
    }
    //non-local entries
    // Jstart = A->offd.rowStarts[i];
    // Jend   = A->offd.rowStarts[i+1];
    // for(dlong jj= Jstart; jj<Jend; jj++){
    //   const dlong col = A->offd.cols[jj];
    //   const dfloat Ajj = fabs(diagA[col]);

    //   if(fabs(A->offd.vals[jj]) > theta*(sqrt(Aii*Ajj)))
    //     strong_per_row++;
    // }

    strong.diag.rowStarts[i+1] = strong_per_row;
  }

  // cumulative sum
  for(dlong i=1; i<Nrows+1 ; i++) {
    strong.diag.rowStarts[i] += strong.diag.rowStarts[i-1];
  }
  strong.diag.nnz = strong.diag.rowStarts[Nrows];
  strong.diag.cols = new dlong[strong.diag.nnz];
  strong.diag.vals = new dfloat[strong.diag.nnz];

  // fill in the columns for strong connections
  // #pragma omp parallel for
  for(dlong i=0; i<Nrows; i++){
    const dfloat Aii = diagA[i];

    dlong counter = strong.diag.rowStarts[i];

    //local entries
    dlong Jstart = diag.rowStarts[i];
    dlong Jend   = diag.rowStarts[i+1];
    for(dlong jj= Jstart; jj<Jend; jj++){
      const dlong col = diag.cols[jj];
      if (col==i) continue;

      const dfloat Ajj = diagA[col];

      if(-diag.vals[jj] > theta*(sqrt(Aii*Ajj))) {
        strong.diag.cols[counter] = col;
        strong.diag.vals[counter++] = -diag.vals[jj] + distrib(paradogs::RNG);
      }
    }
    //non-local entries
    // Jstart = A->offd.rowStarts[i];
    // Jend   = A->offd.rowStarts[i+1];
    // for(dlong jj= Jstart; jj<Jend; jj++){
    //   const dlong col = A->offd.cols[jj];

    //   const dfloat Ajj = fabs(diagA[col]);

    //   if(fabs(A->offd.vals[jj]) > theta*(sqrt(Aii*Ajj)))
    //     strong.cols[counter++] = col;
    // }
  }

  int  *state = new int[Nrows];
  float* rand = new float[Nrows];
  int   *Ts = new int[Nrows];
  float *Tr = new float[Nrows];
  dlong *Tn = new dlong[Nrows];

  /*Initialize state array*/
  /*  0 - Undecided */
  /* -1 - Not MIS */
  /*  1 - MIS */
  #pragma omp parallel for
  for (dlong n=0;n<Nrows;++n) state[n] = 0;

  /*Use vertex degree with random noise to break ties*/
  // #pragma omp parallel for
  for (dlong n=0;n<Nrows;++n) {
    rand[n] = strong.diag.rowStarts[n+1]
              - strong.diag.rowStarts[n]
              + distrib(paradogs::RNG);
  }

  do {
    // first neighbours
    #pragma omp parallel for
    for(dlong n=0; n<Nrows; n++){
      int    smax = state[n];

      if (smax==1) continue;

      float  rmax = rand[n];
      dlong  nmax = n;

      for(dlong j=strong.diag.rowStarts[n];j<strong.diag.rowStarts[n+1];j++){
        const dlong k  = strong.diag.cols[j];
        const int   sk = state[k];
        const float rk = rand[k];
        if ((sk>smax)              || /*If neighbor is MIS node*/
           ((sk==smax)&&(rk>rmax)) || /*Else if it has a bigger weight*/
           ((sk==smax)&&(rk==rmax)&&(k>nmax))) { /*Rare, but just in case, break tie with index number*/
          smax = sk;
          rmax = rk;
          nmax = k;
        }
      }
      Ts[n] = smax;
      Tr[n] = rmax;
      Tn[n] = nmax;
    }

    // second neighbours
    #pragma omp parallel for
    for(dlong n=0; n<Nrows; n++){
      if (state[n]!=0) continue;

      int   smax = Ts[n];
      float rmax = Tr[n];
      dlong nmax = Tn[n];

      for(dlong j=strong.diag.rowStarts[n];j<strong.diag.rowStarts[n+1];j++){
        const dlong k = strong.diag.cols[j];
        const int   sk = Ts[k];
        const float rk = Tr[k];
        const dlong nk = Tn[k];
        if ((sk>smax)              || /*If neighbor is MIS node*/
           ((sk==smax)&&(rk>rmax)) || /*Else if it has a bigger weight*/
           ((sk==smax)&&(rk==rmax)&&(nk>nmax))) { /*Rare, but just in case, break tie with index number*/
          smax = sk;
          rmax = rk;
          nmax = nk;
        }
      }

      // if I am the strongest among all the 1 and 2 ring neighbours
      // I am an MIS node
      if(nmax == n) state[n] = 1;

      // if there is an MIS node within distance 2, I am removed
      if(smax>0) state[n] = -1;
    }

    // if number of undecided nodes = 0, algorithm terminates
    dlong cnt = 0;
    for (dlong n=0;n<Nrows;n++) if (state[n]==0) cnt++;

    if (cnt==0) break;

  } while(true);

  delete[] Ts;
  delete[] Tr;
  delete[] Tn;
  delete[] rand;

  /*Initialize Matching array*/
  Nc=0;
  for(dlong i=0; i<Nrows; i++) {
    if(state[i] == 1) {
      FineToCoarse[i] = Nc++;
    } else {
      FineToCoarse[i] = -1;
    }
  }

  // first neighbours
  #pragma omp parallel for
  for(dlong n=0; n<Nrows; n++){
    if (FineToCoarse[n]==-1) {
      for(dlong j=strong.diag.rowStarts[n];j<strong.diag.rowStarts[n+1];j++){
        const dlong k  = strong.diag.cols[j];
        const int   sk = FineToCoarse[k];

        /*If this node is an MIS node, join the aggregate*/
        if (state[k]==1) {
          FineToCoarse[n] = sk;
          break;
        }
      }
    }
  }

  // second neighbours
  #pragma omp parallel for
  for(dlong n=0; n<Nrows; n++){
    if (FineToCoarse[n]==-1) { //If we're still undecided
      int   smax = -1;
      float rmax = -1.0;
      dlong kmax = -1;

      for(dlong j=strong.diag.rowStarts[n];j<strong.diag.rowStarts[n+1];j++){
        const dlong k = strong.diag.cols[j];
        const int   sk = FineToCoarse[k];
        if (sk!=-1) {
          // const float rk = rand[k];
          const float rk = strong.diag.vals[j];
          if( (rk>rmax)            || /*If edge is strongest*/
             ((rk==rmax)&&(k>kmax))) { /*Rare, but just in case, break tie with index number*/
            smax = sk;
            rmax = rk;
            kmax = k;
          }
        }
      }
      FineToCoarse[n] = smax;
    }
  }

  delete[] state;
}

}
