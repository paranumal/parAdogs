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

#include "mesh.hpp"
#include "mesh/mesh2D.hpp"
#include "mesh/mesh3D.hpp"

//makeing a mesh object requires it to be bound to a device and communicator
mesh_t::mesh_t(platform_t& _platform, meshSettings_t& _settings, MPI_Comm _comm):
  platform(_platform), settings(_settings), comm(_comm) {
  props = platform.props;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
}

mesh2D::mesh2D(platform_t& _platform, meshSettings_t& _settings, MPI_Comm _comm):
  mesh_t(_platform, _settings, _comm) {}

mesh3D::mesh3D(platform_t& _platform, meshSettings_t& _settings, MPI_Comm _comm):
  mesh_t(_platform, _settings, _comm) {}

meshTri2D::meshTri2D(platform_t& _platform, meshSettings_t& _settings, MPI_Comm _comm):
  mesh2D(_platform, _settings, _comm) {

  plotNelements = 1;
  plotNverts = 3;
  plotEToV = (int*) malloc(plotNelements*plotNverts*sizeof(int));
  plotEToV[0] = 0;
  plotEToV[1] = 1;
  plotEToV[2] = 2;
}

meshQuad2D::meshQuad2D(platform_t& _platform, meshSettings_t& _settings, MPI_Comm _comm):
  mesh2D(_platform, _settings, _comm) {

  plotNelements = 2;
  plotNverts = 3;
  plotEToV = (int*) malloc(plotNelements*plotNverts*sizeof(int));

  plotEToV[0*plotNverts+0] = 0;
  plotEToV[0*plotNverts+1] = 1;
  plotEToV[0*plotNverts+2] = 2;

  plotEToV[1*plotNverts+0] = 0;
  plotEToV[1*plotNverts+1] = 3;
  plotEToV[1*plotNverts+2] = 2;
}

meshTri3D::meshTri3D(platform_t& _platform, meshSettings_t& _settings, MPI_Comm _comm):
  mesh3D(_platform, _settings, _comm) {}

meshQuad3D::meshQuad3D(platform_t& _platform, meshSettings_t& _settings, MPI_Comm _comm):
  mesh3D(_platform, _settings, _comm) {}

meshTet3D::meshTet3D(platform_t& _platform, meshSettings_t& _settings, MPI_Comm _comm):
  mesh3D(_platform, _settings, _comm) {

  plotNelements = 1;
  plotNverts = 4;
  plotEToV = (int*) malloc(plotNelements*plotNverts*sizeof(int));

  plotEToV[0] = 0;
  plotEToV[1] = 1;
  plotEToV[2] = 2;
  plotEToV[3] = 3;
}

meshHex3D::meshHex3D(platform_t& _platform, meshSettings_t& _settings, MPI_Comm _comm):
  mesh3D(_platform, _settings, _comm) {

  plotNelements = 6;
  plotNverts = 4;
  plotEToV = (int*) malloc(plotNelements*plotNverts*sizeof(int));

  //Tensor product
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
}

mesh_t::~mesh_t() {
  if (plotEToV) {free(plotEToV); plotEToV=nullptr;}
  if (faceVertices) {free(faceVertices); faceVertices=nullptr;}
  if (boundaryInfo) {free(boundaryInfo); boundaryInfo=nullptr;}
  if (elementInfo) {free(elementInfo); elementInfo=nullptr;}
  if (EToV) {free(EToV); EToV=nullptr;}
  if (EToE) {free(EToE); EToE=nullptr;}
  if (EToF) {free(EToF); EToF=nullptr;}
  if (EX) {free(EX); EX=nullptr;}
  if (EY) {free(EY); EY=nullptr;}
  if (EZ) {free(EZ); EZ=nullptr;}
}
