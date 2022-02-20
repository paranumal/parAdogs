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
#include <random>

namespace paradogs {

extern std::mt19937 RNG;

void MeshPartition(platform_t &platform,
                   settings_t &settings,
                   dlong &Nelements,
                   const  int dim,
                   const  int Nverts,
                   const  int Nfaces,
                   const  int NfaceVertices,
                   const  int* faceVertices,
                   hlong*  &EToV,
                   hlong*  &EToE,
                   int*    &EToF,
                   dfloat* &EX,
                   dfloat* &EY,
                   dfloat* &EZ,
                   MPI_Comm comm) {

  /* Create RNG*/
  int rank;
  int size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  RNG = std::mt19937(rank);

  /* Create graph from mesh info*/
  graph_t graph(platform,
                Nelements,
                dim,
                Nverts,
                Nfaces,
                NfaceVertices,
                faceVertices,
                EToV,
                EX,
                EY,
                EZ,
                comm);

  double timeStart = MPI_Wtime();

  if (settings.compareSetting("PARADOGS PARTITIONING", "INERTIAL")) {
    /*Inertial partitioning*/
    graph.InertialPartition();
  } else if (settings.compareSetting("PARADOGS PARTITIONING", "SPECTRAL")) {
    /*Connect element faces before partitioning*/
    if (size>1) graph.Connect();

    /*Spectral partitioning*/
    graph.SpectralPartition();
  }

  /*Connect element faces after partitioning*/
  graph.Connect();

  /*Reorder rank-local element list for better locality*/
  graph.CuthillMckee();

  double timeEnd = MPI_Wtime();
  double elaplsed = timeEnd-timeStart;
  MPI_Allreduce(MPI_IN_PLACE, &elaplsed, 1, MPI_DOUBLE, MPI_MAX, comm);

  /*Print some stats about the partitioning*/
  graph.Report();

  if (rank==0) {
    printf("   Partitioning time:  %5.2f seconds                                                          |\n",
           elaplsed);
    printf("-----------------------------------------------------------------------------------------------\n");
  }

  /*Get the new mesh data*/
  graph.ExtractMesh(Nelements,
                    EToV,
                    EToE,
                    EToF,
                    EX,
                    EY,
                    EZ);
}

}

