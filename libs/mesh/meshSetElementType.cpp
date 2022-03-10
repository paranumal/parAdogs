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

#include "mesh.hpp"

namespace libp {

void mesh_t::SetElementType(const int eType) {

  if (eType==TRIANGLES) {
    elementType = TRIANGLES;

    dim = 2;           // mesh dimension
    Nverts = 3;        // number of vertices per element
    Nfaces = 3;        // number of faces per element
    NfaceVertices = 2; // number of vertices per face

    // vertices on each face
    int _faceVertices[4][2] = {{0,1},{1,2},{2,0}};

    faceVertices.malloc(NfaceVertices*Nfaces);
    faceVertices.copyFrom(_faceVertices[0]);

    // Plotting
    plotNelements = 1;
    plotNverts = 3;
    plotEToV.malloc(plotNelements*plotNverts);
    plotEToV[0] = 0;
    plotEToV[1] = 1;
    plotEToV[2] = 2;

  } else if (eType==QUADRILATERALS) {
    elementType = QUADRILATERALS;

    dim = 2;           // mesh dimension
    Nverts = 4;        // number of vertices per element
    Nfaces = 4;        // number of faces per element
    NfaceVertices = 2; // number of vertices per face

    // vertices on each face
    int _faceVertices[4][2] = {{0,1},{1,2},{2,3},{3,0}};

    faceVertices.malloc(NfaceVertices*Nfaces);
    faceVertices.copyFrom(_faceVertices[0]);

    // Plotting
    plotNelements = 2;
    plotNverts = 3;
    plotEToV.malloc(plotNelements*plotNverts);

    plotEToV[0*plotNverts+0] = 0;
    plotEToV[0*plotNverts+1] = 1;
    plotEToV[0*plotNverts+2] = 2;

    plotEToV[1*plotNverts+0] = 0;
    plotEToV[1*plotNverts+1] = 3;
    plotEToV[1*plotNverts+2] = 2;

  } else if (eType==TETRAHEDRA) {
    elementType = TETRAHEDRA;

    dim = 3;           // mesh dimension
    Nverts = 4;        // number of vertices per element
    Nfaces = 4;        // number of faces per element
    NfaceVertices = 3; // number of vertices per face

    // vertices on each face
    int _faceVertices[4][3] = {{0,1,2},{0,1,3},{1,2,3},{2,0,3}};

    faceVertices.malloc(NfaceVertices*Nfaces);
    faceVertices.copyFrom(_faceVertices[0]);

    // Plotting
    plotNelements = 1;
    plotNverts = 4;
    plotEToV.malloc(plotNelements*plotNverts);

    plotEToV[0] = 0;
    plotEToV[1] = 1;
    plotEToV[2] = 2;
    plotEToV[3] = 3;

  } else if (eType==HEXAHEDRA) {
    elementType = HEXAHEDRA;

    dim = 3;           // mesh dimension
    Nverts = 8;        // number of vertices per element
    Nfaces = 6;        // number of faces per element
    NfaceVertices = 4; // number of vertices per face

    // vertices on each face
    int _faceVertices[6][4] =
      {{0,1,2,3},{0,4,5,1},{1,5,6,2},{2,6,7,3},{0,3,7,4},{4,7,6,5}};

    faceVertices.malloc(NfaceVertices*Nfaces);
    faceVertices.copyFrom(_faceVertices[0]);

    // Plotting
    plotNelements = 6;
    plotNverts = 4;
    plotEToV.malloc(plotNelements*plotNverts);

    //tet 1 (0,3,2,7)
    plotEToV[0*plotNverts+0] = 0;
    plotEToV[0*plotNverts+1] = 2;
    plotEToV[0*plotNverts+2] = 3;
    plotEToV[0*plotNverts+3] = 6;
    //tet 2 (0,1,3,7)
    plotEToV[1*plotNverts+0] = 0;
    plotEToV[1*plotNverts+1] = 1;
    plotEToV[1*plotNverts+2] = 2;
    plotEToV[1*plotNverts+3] = 6;
    //tet 3 (0,2,6,7)
    plotEToV[2*plotNverts+0] = 0;
    plotEToV[2*plotNverts+1] = 3;
    plotEToV[2*plotNverts+2] = 7;
    plotEToV[2*plotNverts+3] = 6;
    //tet 4 (0,6,4,7)
    plotEToV[3*plotNverts+0] = 0;
    plotEToV[3*plotNverts+1] = 7;
    plotEToV[3*plotNverts+2] = 4;
    plotEToV[3*plotNverts+3] = 6;
    //tet 5 (0,5,1,7)
    plotEToV[4*plotNverts+0] = 0;
    plotEToV[4*plotNverts+1] = 5;
    plotEToV[4*plotNverts+2] = 1;
    plotEToV[4*plotNverts+3] = 6;
    //tet 6 (0,4,5,7)
    plotEToV[5*plotNverts+0] = 0;
    plotEToV[5*plotNverts+1] = 4;
    plotEToV[5*plotNverts+2] = 5;
    plotEToV[5*plotNverts+3] = 6;
  } else {
    LIBP_FORCE_ABORT("Unknown element type: " << eType);
  }
}

} //namespace libp
