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

#ifndef PARADOGS_GRAPH_HPP
#define PARADOGS_GRAPH_HPP 1

#include "parAdogs.hpp"
#include "parAdogs/parAdogsDefines.h"
#include "parAdogs/parAdogsMatrix.hpp"
#include "parAdogs/parAdogsMultigrid.hpp"

namespace paradogs {

#define MAX_NVERTS 8
#define MAX_NFACES 6
#define MAX_NFACEVERTS 4

class graph_t {
private:
  MPI_Comm gcomm=MPI_COMM_NULL;
  MPI_Comm comm=MPI_COMM_NULL;

  int rank, size;
  dlong Nverts=0;
  hlong NVertsGlobal=0;
  hlong VoffsetL=0, VoffsetU=0;

  int grank, gsize;
  hlong gNVertsGlobal=0;
  hlong gVoffsetL=0, gVoffsetU=0;
  
  /*Mesh data*/
  dlong Nelements=0;
  int dim=0;
  int Nfaces=0;
  int NelementVerts=0;
  int NfaceVerts=0;
  struct element_t {
    dfloat EX[MAX_NVERTS]; //x coordinates of verts
    dfloat EY[MAX_NVERTS]; //y coordinates of verts
    dfloat EZ[MAX_NVERTS]; //z coordinates of verts
    hlong V[MAX_NVERTS];   //Global Vertex Ids of verts

    hlong E[MAX_NFACES];   //Global element ids of neighbors
    int F[MAX_NFACES];     //Face ids of neighbors
  } *elements=nullptr;

  int faceVerts[MAX_NFACES*MAX_NFACEVERTS];

  /*Multilevel Laplacian (for spectral partitioning)*/
  int Nlevels=0;
  mgLevel_t L[PARADOGS_MAX_LEVELS];
  coarseSolver_t coarseSolver;

public:
  /*Build a graph from mesh connectivity info*/
  graph_t(const dlong _Nelements,
          const int _dim,
          const int _Nverts,
          const int _Nfaces,
          const int _NfaceVerts,
          const int* _faceVerts,
          const hlong EToV[],
          const dfloat EX[],
          const dfloat EY[],
          const dfloat EZ[],
          MPI_Comm _comm);

  ~graph_t();

  void InertialPartition();

  void SpectralPartition();

  void Connect();

  void CuthillMckee();

  void Report();

  void ExtractMesh(dlong &Nelements_,
                   hlong*  &EToV,
                   hlong*  &EToE,
                   int*    &EToF,
                   dfloat* &EX,
                   dfloat* &EY,
                   dfloat* &EZ);

private:
  void InertialBipartition(const dfloat targetFraction[2]);
  void SpectralBipartition(const dfloat targetFraction[2]);


  /*Divide graph into two pieces according to a bisection*/
  void Split(const int partition[]);

  void CreateLaplacian();

  /*Compute Fiedler vector of graph */
  dfloat* FiedlerVector();

  /*Improve a Fiedler vector*/
  void Refine(const int level);

  /* Solve A_{l}*x = b*/
  int Solve(const int level,
            const dfloat TOL,
            dfloat r[],
            dfloat x[],
            dfloat scratch[]);

  /*Create multilevel heirarchy*/
  void MultigridSetup();

  void MultigridVcycle(const int l,
                      dfloat r[],
                      dfloat x[]);

  /*Clear multilevel heirarchy*/
  void MultigridDestroy();
};

}

#endif

