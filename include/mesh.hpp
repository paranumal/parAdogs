/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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
#include "ogs.hpp"

namespace libp {

class meshSettings_t: public settings_t {
public:
  meshSettings_t(comm_t _comm);
  void report();
};

class mesh_t {
public:
  platform_t platform;
  settings_t settings;
  properties_t props;

  comm_t comm;
  int rank, size;

  /*************************/
  /* Element Data          */
  /*************************/
  int dim;
  int Nverts, Nfaces, NfaceVertices;
  int elementType;

  // indices of vertex nodes
  memory<int> vertexNodes;

  hlong Nnodes=0; //global number of element vertices
  memory<dfloat> EX; // coordinates of vertices for each element
  memory<dfloat> EY;
  memory<dfloat> EZ;

  dlong Nelements=0;       //local element count
  hlong NelementsGlobal=0; //global element count
  memory<hlong> EToV; // element-to-vertex connectivity
  memory<hlong> EToE; // element-to-element connectivity
  memory<int>   EToF; // element-to-(local)face connectivity
  memory<int>   EToP; // element-to-partition/process connectivity
  memory<int>   EToB; // element-to-boundary condition type

  memory<hlong> elementInfo; //type of element

  memory<dlong> VmapM;  // list of vertices on each face
  memory<dlong> VmapP;  // list of vertices that are paired with face vertices

  memory<int> faceVertices; // list of mesh vertices on each face

  // boundary faces
  hlong NboundaryFaces=0; // number of boundary faces
  memory<hlong> boundaryInfo; // list of boundary faces (type, vertex-1, vertex-2, vertex-3)

  int plotNverts=0;            // number of vertices for each plot element
  int plotNelements=0;         // number of "plot elements" per element
  memory<int> plotEToV;  // triangulation of plot nodes

  mesh_t()=default;
  mesh_t(platform_t& _platform,
         meshSettings_t& _settings,
         comm_t _comm) {
    Setup(_platform, _settings, _comm);
  }

  ~mesh_t() = default;

  void Setup(platform_t& _platform,
             meshSettings_t& _settings,
             comm_t _comm);

private:
  /*Element types*/
  static constexpr int TRIANGLES     =3;
  static constexpr int QUADRILATERALS=4;
  static constexpr int TETRAHEDRA    =6;
  static constexpr int HEXAHEDRA     =12;

  /*Set the type of mesh*/
  void SetElementType(const int eType);

  // box mesh
  void SetupBox() {
    switch (elementType) {
      case TRIANGLES:
        SetupBoxTri2D();
        break;
      case QUADRILATERALS:
        SetupBoxQuad2D();
        break;
      case TETRAHEDRA:
        SetupBoxTet3D();
        break;
      case HEXAHEDRA:
        SetupBoxHex3D();
        break;
    }
  }
  void SetupBoxTri2D();
  void SetupBoxQuad2D();
  void SetupBoxTet3D();
  void SetupBoxHex3D();

  // mesh reader
  void ParallelReader(const std::string fileName) {
    switch (elementType) {
      case TRIANGLES:
        ParallelReaderTri2D(fileName);
        break;
      case QUADRILATERALS:
        ParallelReaderQuad2D(fileName);
        break;
      case TETRAHEDRA:
        ParallelReaderTet3D(fileName);
        break;
      case HEXAHEDRA:
        ParallelReaderHex3D(fileName);
        break;
    }
  }
  void ParallelReaderTri2D(const std::string fileName);
  void ParallelReaderQuad2D(const std::string fileName);
  void ParallelReaderTet3D(const std::string fileName);
  void ParallelReaderHex3D(const std::string fileName);

  // repartition elements
  void Partition();

  void Plot(const memory<dfloat> q);
};

} //namespace libp

#endif

