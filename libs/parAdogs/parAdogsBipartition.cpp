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
#include <algorithm>

#ifdef GLIBCXX_PARALLEL_SORT
#include <parallel/algorithm>
using __gnu_parallel::sort;
#else
using std::sort;
#endif

namespace paradogs {

struct keyVal_t {
  dfloat key;
  dlong val;

  keyVal_t() {};
  keyVal_t(const dfloat _key, const dlong _val):
    key{_key}, val{_val} {}

  friend bool operator<(const keyVal_t& a, const keyVal_t& b) {
    return a.key < b.key;
  }
};

/****************************************/
/* Serial Multilevel Spectral Bisection */
/****************************************/
void Bipartition(graph_t& graph,
                 const dfloat targetFraction[2],
                 int partition[]) {

  /*Compute Fiedler vector */
  dfloat *Fiedler = graph.FiedlerVector();

  /*Use Fiedler vector to bipartion graph*/
  keyVal_t *F = new keyVal_t[graph.Nverts];
  for (dlong n=0;n<graph.Nverts;++n) {
    F[n] = keyVal_t(Fiedler[n], n);
  }

  sort(F, F+graph.Nverts);

  dlong Nverts0 = std::ceil(targetFraction[0]*graph.Nverts);

  for (dlong n=0;n<Nverts0;++n) {
    partition[F[n].val] = 0;
  }
  for (dlong n=Nverts0;n<graph.Nverts;++n) {
    partition[F[n].val] = 1;
  }

  delete[] F;
}

}

