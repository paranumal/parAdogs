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
#include "parAdogs/parAdogsPartition.hpp"

namespace paradogs {

/****************************************/
/* Multigrid vcycle                     */
/****************************************/
void graph_t::vcycle(const int l,
                     dfloat r[],
                     dfloat x[]) {

  //check for base level
  if(l==Nlevels-1) {
    coarseSolver.Solve(r, x);
    return;
  }

  mgLevel_t& Lf = L[l];
  dfloat *res = Lf.RES;

  mgLevel_t& Lc = L[l+1];
  dfloat *rC = Lc.RHS;
  dfloat *xC = Lc.X;

  //Pre smooth and then compute res = rhs-Ax
  Lf.Smooth(r, x, true);
  Lf.Residual(r, x, res);

  // rhsC = P^T res
  Lf.Coarsen(res, rC);

  // Recursive call
  vcycle(l+1, rC, xC);

  // x = x + P xC
  Lf.Prolongate(xC, x);

  // Post smooth
  Lf.Smooth(r, x, false);
}

}


