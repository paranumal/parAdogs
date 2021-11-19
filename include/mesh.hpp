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

#ifndef MESH_HPP
#define MESH_HPP 1

#include "core.hpp"
#include "settings.hpp"
#include "platform.hpp"

#define TRIANGLES 3
#define QUADRILATERALS 4
#define TETRAHEDRA 6
#define HEXAHEDRA 12

class meshSettings_t: public settings_t {
public:
  meshSettings_t(MPI_Comm& _comm);
  void report();
};

class mesh_t {
public:
  platform_t& platform;
  meshSettings_t& settings;

  occa::properties props;

  MPI_Comm comm;
  int rank, size;

  int dim;
  int Nverts, Nfaces, NfaceVertices;
  int *faceVertices; // list of mesh vertices on each face

  // indices of vertex nodes
  int *vertexNodes;

  int elementType;

  hlong Nnodes=0; //global number of element vertices
  dfloat *EX=nullptr; // coordinates of vertices for each element
  dfloat *EY=nullptr;
  dfloat *EZ=nullptr;

  dlong Nelements=0;       //local element count
  hlong NelementsGlobal=0; //global element count
  hlong *EToV=nullptr; // element-to-vertex connectivity
  hlong *EToE=nullptr; // element-to-element connectivity
  int   *EToF=nullptr; // element-to-(local)face connectivity
  int   *EToB=nullptr; // element-to-boundary condition type

  hlong *elementInfo=nullptr; //type of element

  // boundary faces
  hlong NboundaryFaces=0; // number of boundary faces
  hlong *boundaryInfo=nullptr; // list of boundary faces (type, vertex-1, vertex-2, vertex-3)

  int    plotNverts=0;    // number of vertices for each plot element
  int    plotNelements=0; // number of "plot elements" per element
  int    *plotEToV=nullptr;       // triangulation of plot nodes

  mesh_t() = delete;
  mesh_t(platform_t& _platform, meshSettings_t& _settings,
         MPI_Comm _comm);

  virtual ~mesh_t();

  // generic mesh setup
  static mesh_t& Setup(platform_t& _platform, meshSettings_t& _settings,
                       MPI_Comm _comm);

  // box mesh
  virtual void SetupBox() = 0;

  // mesh reader
  virtual void ParallelReader(const char *fileName) = 0;

  // repartition elements
  void Partition();

  void Plot(const dfloat* q);
};

#endif

