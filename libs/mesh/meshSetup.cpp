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

namespace libp {

void mesh_t::Setup(platform_t& _platform, meshSettings_t& _settings,
                   MPI_Comm _comm){

  platform = _platform;
  settings = _settings;
  props = platform.props();

  MPI_Comm_dup(_comm, &comm);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  std::string fileName;
  settings.getSetting("MESH FILE", fileName);
  settings.getSetting("ELEMENT TYPE", elementType);
  settings.getSetting("MESH DIMENSION", dim);

  SetElementType(elementType);

  if (settings.compareSetting("MESH FILE","BOX")) {
    //build a box mesh
    SetupBox();
  } else {
    // read chunk of elements from file
    ParallelReader(fileName);
  }

  // partition mesh among processes
  Partition();
}

} //namespace libp
