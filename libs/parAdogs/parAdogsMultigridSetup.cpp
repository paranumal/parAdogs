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
/* Construct Multigrid Hierarchy        */
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


  // estimate rho(invD * A)
  A.rho = A.rhoDinvA(null);

  /*Smoothing params*/
  lambda1 = A.rho;
  lambda0 = A.rho/10.;

  /*Scratch space*/
  scratch = new dfloat[2*A.Ncols];
}

/*Coarsen a graph using an aggregation*/
void mgLevel_t::CoarsenLevel(mgLevel_t& Lf) {

  /*Create a FineToCoarse mapping*/
  const dlong Nf = Lf.A.Nrows;

  /*Create a vertex matching*/
  dlong Nc=0;
  dlong *FineToCoarse = new dlong[Nf];
  Lf.A.Aggregate(Nc, FineToCoarse);

  /*Create coarse nullvector*/
  null = new dfloat[Nc];

  /* Tentative prolongation operator*/
  parCSR T = TentativeProlongator(Nf, Nc, FineToCoarse,
                                  Lf.null, null);
  delete[] FineToCoarse;

  /* Smoothed prologontion */
  Lf.P = SmoothProlongator(Lf.A, T);
  T.Free();

  /* R = P^T*/
  Lf.R = Transpose(Lf.P);

  /*Galerkin product*/
  parCSR AP = SpMM(Lf.A, Lf.P);
  A = SpMM(Lf.R, AP);
  AP.Free();
  // A.GalerkinProduct(Lf.A, Lf.P);

  /*fill diagonal*/
  A.diagA = new dfloat[A.Nrows];
  A.diagInv = new dfloat[A.Nrows];

  #pragma omp parallel for
  for (dlong i=0;i<A.Nrows;i++) {
    const dlong start = A.diag.rowStarts[i];
    const dlong end   = A.diag.rowStarts[i+1];

    for (dlong j=start;j<end;j++) {
      //record the diagonal
      if (A.diag.cols[j]==i) {
        A.diagA[i] = A.diag.vals[j];
        A.diagInv[i] = 1.0/A.diagA[i];
      }
    }
  }

  /*Space for Fiedler*/
  Fiedler = new dfloat[A.Nrows];

  /*Multigrid buffers*/
  RHS = new dfloat[A.Nrows];
  X = new dfloat[A.Nrows];
  RES = new dfloat[A.Nrows];
}

/*Free coarse levels of hierarchy*/
void graph_t::MultigridDestroy() {
  coarseSolver.Free();
  for (int n=Nlevels-1;n>0;--n) L[n].Free();
  Nlevels=1;
}

}


