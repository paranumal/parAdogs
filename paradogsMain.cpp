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
#include "paradogsMain.hpp"

int main(int argc, char **argv){

  // start up MPI
  comm_t::Init(argc, argv);

  {
    comm_t comm(comm_t::world().Dup());

    LIBP_ABORT(std::string("Usage: ./paradogsMain setupfile"),
               argc!=2);

    //create default settings
    platformSettings_t platformSettings(comm);
    meshSettings_t meshSettings(comm);
    paradogsSettings_t paradogsSettings(comm);

    //load settings from file
    paradogsSettings.parseFromFile(platformSettings, meshSettings,
                                   argv[1]);

    // set up platform
    platform_t platform(platformSettings);

    platformSettings.report();
    meshSettings.report();
    paradogsSettings.report();

    // set up mesh
    mesh_t mesh(platform, meshSettings, comm);
  }

  // close down MPI
  comm_t::Finalize();
  return LIBP_SUCCESS;
}
