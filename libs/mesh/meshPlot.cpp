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

//plot field at vertices
void mesh_t::Plot(const memory<dfloat> q){

  FILE *fp;

  char fileName[BUFSIZ];
  sprintf(fileName, "mesh_%04d.vtu", rank);

  fp = fopen(fileName, "w");

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");
  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
          Nelements*Nverts,
          Nelements*plotNelements);

  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  if (dim==2) {
    for(dlong e=0;e<Nelements;++e){
      for(int n=0;n<Nverts;++n){
        fprintf(fp, "       ");
        fprintf(fp, "%g %g 0.0\n", EX[n+e*Nverts],
                                   EY[n+e*Nverts]);
      }
    }
  } else if (dim==3) {
    for(dlong e=0;e<Nelements;++e){
      for(int n=0;n<Nverts;++n){
        fprintf(fp, "       ");
        fprintf(fp, "%g %g %g\n", EX[n+e*Nverts],
                                  EY[n+e*Nverts],
                                  EZ[n+e*Nverts]);
      }
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");

  // write out q
  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"q\" Format=\"ascii\">\n");
  for(dlong e=0;e<Nelements;++e){
    for(int n=0;n<Nverts;++n){
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", q[n+e*Nverts]);
    }
  }
  fprintf(fp, "       </DataArray>\n");
  fprintf(fp, "     </PointData>\n");

  fprintf(fp, "    <Cells>\n");
  fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");

  for(dlong e=0;e<Nelements;++e){
    for(int n=0;n<plotNelements;++n){
      fprintf(fp, "       ");
      for(int m=0;m<plotNverts;++m){
        fprintf(fp, "%d ", e*Nverts + plotEToV[n*plotNverts+m]);
      }
      fprintf(fp, "\n");
    }
  }
  fprintf(fp, "        </DataArray>\n");

  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  dlong cnt = 0;
  for(dlong e=0;e<Nelements;++e){
    for(int n=0;n<plotNelements;++n){
      cnt += plotNverts;
      fprintf(fp, "       ");
      fprintf(fp, "%d\n", cnt);
    }
  }
  fprintf(fp, "       </DataArray>\n");

  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(dlong e=0;e<Nelements;++e){
    for(int n=0;n<plotNelements;++n){
      if(dim==2)
        fprintf(fp, "5\n");
      else
        fprintf(fp, "10\n");
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

} //namespace libp
