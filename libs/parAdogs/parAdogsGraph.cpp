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

namespace libp {

namespace paradogs {

/*Build a graph from mesh connectivity info*/
graph_t::graph_t(platform_t &_platform,
                 const dlong _Nelements,
                 const int _dim,
                 const int _Nverts,
                 const int _Nfaces,
                 const int _NfaceVerts,
                 const libp::memory<int>& faceVertices,
                 const libp::memory<hlong>& EToV,
                 const libp::memory<dfloat>& EX,
                 const libp::memory<dfloat>& EY,
                 const libp::memory<dfloat>& EZ,
                 MPI_Comm _comm):
  platform(_platform),
  Nverts(_Nelements),
  Nelements(_Nelements),
  dim(_dim),
  Nfaces(_Nfaces),
  NelementVerts(_Nverts),
  NfaceVerts(_NfaceVerts) {

  MPI_Comm_dup(_comm, &gcomm);
  MPI_Comm_rank(gcomm, &grank);
  MPI_Comm_size(gcomm, &gsize);

  MPI_Comm_dup(_comm, &comm);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  for (int n=0;n<Nfaces*NfaceVerts;++n)
    faceVerts[n] = faceVertices[n];

  /*Global number of elements*/
  hlong localNverts=static_cast<hlong>(Nverts);
  MPI_Allreduce(&localNverts, &NVertsGlobal, 1, MPI_HLONG, MPI_SUM, comm);

  /*Get global element count offsets*/
  MPI_Scan(&localNverts, &VoffsetU, 1, MPI_HLONG, MPI_SUM, comm);
  VoffsetL = VoffsetU-Nverts;

  gNVertsGlobal = NVertsGlobal;
  gVoffsetL = VoffsetL;
  gVoffsetU = VoffsetU;

  /*Create array of packed element data*/
  elements = new element_t[Nelements];

  if (dim==2) {
    for (dlong e=0;e<Nelements;++e) {
      for (int v=0;v<NelementVerts;++v) {
        elements[e].EX[v] = EX[v+e*NelementVerts];
        elements[e].EY[v] = EY[v+e*NelementVerts];

        elements[e].V[v] = EToV[v+e*NelementVerts];
      }
      for (int f=0;f<Nfaces;++f) {
        elements[e].E[f] = -1;
        elements[e].F[f] = -1;
      }
    }
  } else {
    for (dlong e=0;e<Nelements;++e) {
      for (int v=0;v<NelementVerts;++v) {
        elements[e].EX[v] = EX[v+e*NelementVerts];
        elements[e].EY[v] = EY[v+e*NelementVerts];
        elements[e].EZ[v] = EZ[v+e*NelementVerts];

        elements[e].V[v] = EToV[v+e*NelementVerts];
      }
      for (int f=0;f<Nfaces;++f) {
        elements[e].E[f] = -1;
        elements[e].F[f] = -1;
      }
    }
  }
}

/*Globally divide graph into two pieces according to a bipartition*/
void graph_t::Split(const int partition[]) {

  /*Count how much of each partition we have locally*/
  dlong Nverts0=0;
  dlong Nverts1=0;
  for (dlong n=0;n<Nverts;++n) {
    if (partition[n]==0) Nverts0++;
    else                 Nverts1++;
  }

  hlong localNverts0 = static_cast<hlong>(Nverts0);
  hlong localNverts1 = static_cast<hlong>(Nverts1);
  hlong globalNverts0=0;
  hlong globalNverts1=0;
  MPI_Allreduce(&localNverts0, &globalNverts0, 1, MPI_HLONG, MPI_SUM, comm);
  MPI_Allreduce(&localNverts1, &globalNverts1, 1, MPI_HLONG, MPI_SUM, comm);

  /*Get offsets of partitions on each rank*/
  hlong *starts0 = new hlong[size+1];
  hlong *starts1 = new hlong[size+1];
  starts0[0]=0;
  starts1[0]=0;
  MPI_Allgather(&localNverts0, 1, MPI_HLONG, starts0+1, 1,  MPI_HLONG, comm);
  MPI_Allgather(&localNverts1, 1, MPI_HLONG, starts1+1, 1,  MPI_HLONG, comm);

  for(int r=0;r<size;++r) {
    starts0[r+1] += starts0[r];
    starts1[r+1] += starts1[r];
  }

  /*Determine number of ranks to hold left and right partitions*/
  const int size0 = (size+1)/2;
  const int size1 = size-size0;

  const hlong chunk0 = globalNverts0/size0;
  const hlong chunk1 = globalNverts1/size1;

  const int remainder0 = static_cast<int>(globalNverts0 - chunk0*size0);
  const int remainder1 = static_cast<int>(globalNverts1 - chunk1*size1);

  // Make the MPI_ELEMENT_T data type
  MPI_Datatype MPI_ELEMENT_T;
  MPI_Datatype dtype[6] = {MPI_DFLOAT, MPI_DFLOAT, MPI_DFLOAT,
                           MPI_HLONG, MPI_HLONG, MPI_INT};
  int blength[6] = {MAX_NVERTS, MAX_NVERTS, MAX_NVERTS,
                    MAX_NVERTS, MAX_NFACES, MAX_NFACES};
  MPI_Aint addr[6], displ[6];
  MPI_Get_address ( &(elements[0]      ), addr+0);
  MPI_Get_address ( &(elements[0].EY[0]), addr+1);
  MPI_Get_address ( &(elements[0].EZ[0]), addr+2);
  MPI_Get_address ( &(elements[0].V[0] ), addr+3);
  MPI_Get_address ( &(elements[0].E[0] ), addr+4);
  MPI_Get_address ( &(elements[0].F[0] ), addr+5);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  displ[3] = addr[3] - addr[0];
  displ[4] = addr[4] - addr[0];
  displ[5] = addr[5] - addr[0];
  MPI_Type_create_struct (6, blength, displ, dtype, &MPI_ELEMENT_T);
  MPI_Type_commit (&MPI_ELEMENT_T);

  int *Nsend0 = new int[size];
  int *Nsend1 = new int[size];
  int *Nrecv0 = new int[size];
  int *Nrecv1 = new int[size];
  int *sendOffsets0 = new int[size];
  int *sendOffsets1 = new int[size];
  int *recvOffsets0 = new int[size];
  int *recvOffsets1 = new int[size];

  hlong *newIds = new hlong[Nverts+Nhalo];

  for (int r=0;r<size;++r) {
    Nsend0[r]=0;
    Nsend1[r]=0;
  }

  /*Determine new ids and send counts*/
  dlong cnt0=0;
  dlong cnt1=0;
  for(dlong e=0;e<Nverts;++e){
    if (partition[e]==0) {
      // new global element index
      const hlong ep = starts0[rank]+cnt0++;
      newIds[e] = ep;

      // 0, chunk+1, 2*(chunk+1) ..., remainder*(chunk+1), remainder*(chunk+1) + chunk
      int r;
      if(ep<remainder0*(chunk0+1))
        r = ep/(chunk0+1);
      else
        r = remainder0 + ((ep-remainder0*(chunk0+1))/chunk0);

      ++Nsend0[r];
    } else {
      // new global element index
      const hlong ep = starts1[rank]+cnt1++;
      newIds[e] = ep;

      // 0, chunk+1, 2*(chunk+1) ..., remainder*(chunk+1), remainder*(chunk+1) + chunk
      int r;
      if(ep<remainder1*(chunk1+1))
        r = ep/(chunk1+1);
      else
        r = remainder1 + ((ep-remainder1*(chunk1+1))/chunk1);

      ++Nsend1[r+size0];
    }
  }

  delete[] starts0;
  delete[] starts1;

  if (L[0].A) {
    /*If we have connected the elements, share the newIds*/
    L[0].A->halo.Exchange(newIds, 1, ogs::Hlong);

    /*Then update the connectivity*/
    dlong cnt=0;
    for(dlong e=0;e<Nverts;++e){
      const int part = partition[e];
      for (int f=0;f<Nfaces;++f) {
        const hlong gE = elements[e].E[f];
        if (gE!=-1) {
          dlong eN;
          if (gE>=VoffsetL && gE<VoffsetU) { /*local neighbor*/
            eN = static_cast<dlong>(gE-VoffsetL);
          } else { /*halo neighbor*/
            eN = colIds[cnt++]; /*Get the local id in the halo (we make this when building the Laplacian)*/
          }

          const int partN = partition[eN];
          if (partN==part) { /*If both elements are in the same partition*/
            elements[e].E[f] = newIds[eN]; /*Re index*/
          } else {
            elements[e].E[f] = -1;/*else break connections across the partitions*/
          }
        }
      }
    }
  }
  delete[] newIds;

  // find send offsets
  sendOffsets0[0]=0;
  sendOffsets1[0]=0;
  for(int r=1;r<size;++r) {
    sendOffsets0[r] = sendOffsets0[r-1] + Nsend0[r-1];
    sendOffsets1[r] = sendOffsets1[r-1] + Nsend1[r-1];
  }
  int NsendTotal0=0;
  int NsendTotal1=0;
  for(int r=0;r<size;++r) {
    NsendTotal0 += Nsend0[r];
    NsendTotal1 += Nsend1[r];
  }

  // exchange counts
  MPI_Alltoall(Nsend0, 1, MPI_INT, Nrecv0, 1, MPI_INT, comm);
  MPI_Alltoall(Nsend1, 1, MPI_INT, Nrecv1, 1, MPI_INT, comm);

  // find recv offsets
  recvOffsets0[0]=0;
  recvOffsets1[0]=0;
  for(int r=1;r<size;++r) {
    recvOffsets0[r] = recvOffsets0[r-1] + Nrecv0[r-1];
    recvOffsets1[r] = recvOffsets1[r-1] + Nrecv1[r-1];
  }

  // count incoming clusters
  dlong newNverts = 0;

  if (rank<size0) {
    for(int r=0;r<size;++r) {
      newNverts += Nrecv0[r];
    }
  } else {
    for(int r=0;r<size;++r) {
      newNverts += Nrecv1[r];
    }
  }

  /*make send buffers*/
  element_t *sendElements0 = new element_t[NsendTotal0];
  element_t *sendElements1 = new element_t[NsendTotal1];

  cnt0=0;
  cnt1=0;
  for(dlong e=0;e<Nverts;++e){
    if (partition[e]==0) {
      sendElements0[cnt0++] = elements[e];
    } else {
      sendElements1[cnt1++] = elements[e];
    }
  }

  /*free old element list*/
  delete[] elements;

  /*make new list*/
  Nverts = newNverts;
  Nelements = newNverts;
  elements = new element_t[Nverts];

  // exchange elements
  if (rank<size0) {
    MPI_Alltoallv(sendElements0, Nsend0, sendOffsets0, MPI_ELEMENT_T,
                  elements, Nrecv0, recvOffsets0, MPI_ELEMENT_T, comm);
    MPI_Alltoallv(sendElements1, Nsend1, sendOffsets1, MPI_ELEMENT_T,
                  NULL, Nrecv1, recvOffsets1, MPI_ELEMENT_T, comm);
  } else {
    MPI_Alltoallv(sendElements0, Nsend0, sendOffsets0, MPI_ELEMENT_T,
                  NULL, Nrecv0, recvOffsets0, MPI_ELEMENT_T, comm);
    MPI_Alltoallv(sendElements1, Nsend1, sendOffsets1, MPI_ELEMENT_T,
                  elements, Nrecv1, recvOffsets1, MPI_ELEMENT_T, comm);
  }

  MPI_Barrier(comm);
  MPI_Type_free(&MPI_ELEMENT_T);

  delete[] sendElements0;
  delete[] sendElements1;

  delete[] Nsend0;
  delete[] Nsend1;
  delete[] Nrecv0;
  delete[] Nrecv1;
  delete[] sendOffsets0;
  delete[] sendOffsets1;
  delete[] recvOffsets0;
  delete[] recvOffsets1;

  MPI_Comm newComm;
  MPI_Comm_split(comm, rank<size0, rank, &newComm);
  MPI_Comm_free(&comm);
  comm = newComm;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  /*Global number of elements*/
  hlong localNverts=static_cast<hlong>(Nverts);
  MPI_Allreduce(&localNverts, &NVertsGlobal, 1, MPI_HLONG, MPI_SUM, comm);

  /*Get global element count offsets*/
  MPI_Scan(&localNverts, &VoffsetU, 1, MPI_HLONG, MPI_SUM, comm);
  VoffsetL = VoffsetU-Nverts;
}

void graph_t::Report() {

  /* Min,Avg,Max Element counts*/
  hlong globalNverts = static_cast<hlong>(Nverts);
  MPI_Allreduce(MPI_IN_PLACE, &globalNverts, 1, MPI_HLONG, MPI_SUM, gcomm);
  dfloat avgNverts = static_cast<dfloat>(globalNverts)/gsize;

  dlong minNverts=0;
  dlong maxNverts=0;
  MPI_Allreduce(&Nverts, &minNverts, 1, MPI_DLONG, MPI_MIN, gcomm);
  MPI_Allreduce(&Nverts, &maxNverts, 1, MPI_DLONG, MPI_MAX, gcomm);


  dlong cut=0.0;
  for (dlong n=0;n<Nverts;++n) {
    for (int f=0;f<Nfaces;++f) {
      const hlong eN = elements[n].E[f];
      if (eN!=-1) {
        if ((eN<gVoffsetL) || (eN>=gVoffsetU) ) {
          cut++;
        }
      }
    }
  }

  hlong gCut = static_cast<hlong>(cut);
  MPI_Allreduce(MPI_IN_PLACE, &gCut, 1, MPI_HLONG, MPI_SUM, gcomm);
  hlong avgCut = gCut/gsize;

  dlong minCut=0;
  dlong maxCut=0;
  MPI_Allreduce(&cut, &minCut, 1, MPI_DLONG, MPI_MIN, gcomm);
  MPI_Allreduce(&cut, &maxCut, 1, MPI_DLONG, MPI_MAX, gcomm);

  if(grank==0) {
    printf("--------------------------------------ParAdogs Report------------------------------------------\n");
    printf("-----------------------------------------------------------------------------------------------\n");
    printf("   Nranks   |    Elements   |   Per Rank Elements   |   Halo Faces   |   Per Rank Halo Faces  |\n");
    printf("            |               |       (min,avg,max)   |                |         (min,avg,max)  |\n");
    printf("-----------------------------------------------------------------------------------------------\n");
    printf(      "%9d   | %11lld   |       %13lld   | %12lld   |         %13lld  |\n",
            gsize,
            static_cast<long long int>(globalNverts),
            static_cast<long long int>(minNverts),
            static_cast<long long int>(gCut),
            static_cast<long long int>(minCut));
    printf("            |               |       %13lld   |                |         %13lld  |\n",
            static_cast<long long int>(avgNverts),
            static_cast<long long int>(avgCut));
    printf("            |               |       %13lld   |                |         %13lld  |\n",
            static_cast<long long int>(maxNverts),
            static_cast<long long int>(maxCut));
    printf("-----------------------------------------------------------------------------------------------\n");
  }
}

void graph_t::ExtractMesh(dlong &Nelements_,
                          libp::memory<hlong>& EToV,
                          libp::memory<hlong>& EToE,
                          libp::memory<int>& EToF,
                          libp::memory<dfloat>& EX,
                          libp::memory<dfloat>& EY,
                          libp::memory<dfloat>& EZ) {

  /*Destroy any exiting mesh data and create new data from current graph*/
  Nelements_ = Nelements;

  EToV.malloc(Nelements*NelementVerts);
  EToE.malloc(Nelements*NelementVerts);
  EToF.malloc(Nelements*NelementVerts);

  EX.malloc(Nelements*NelementVerts);
  EY.malloc(Nelements*NelementVerts);
  if (dim==3)
    EZ.malloc(Nelements*NelementVerts);

  if (dim==2) {
    for (dlong e=0;e<Nelements;++e) {
      for (int v=0;v<NelementVerts;++v) {
        EToV[v+e*NelementVerts] = elements[e].V[v];
        EX[v+e*NelementVerts] = elements[e].EX[v];
        EY[v+e*NelementVerts] = elements[e].EY[v];
      }
      for (int f=0;f<Nfaces;++f) {
        EToE[f+e*Nfaces] = elements[e].E[f];
        EToF[f+e*Nfaces] = elements[e].F[f];
      }
    }
  } else {
    for (dlong e=0;e<Nelements;++e) {
      for (int v=0;v<NelementVerts;++v) {
        EToV[v+e*NelementVerts] = elements[e].V[v];
        EX[v+e*NelementVerts] = elements[e].EX[v];
        EY[v+e*NelementVerts] = elements[e].EY[v];
        EZ[v+e*NelementVerts] = elements[e].EZ[v];
      }
      for (int f=0;f<Nfaces;++f) {
        EToE[f+e*Nfaces] = elements[e].E[f];
        EToF[f+e*Nfaces] = elements[e].F[f];
      }
    }
  }
}

graph_t::~graph_t() {
  if (gcomm!=MPI_COMM_NULL) MPI_Comm_free(&gcomm);
  if (comm!=MPI_COMM_NULL) MPI_Comm_free(&comm);

  if (elements) {delete[] elements; elements=nullptr;}
  if (colIds) {delete[] colIds; colIds=nullptr;}
}

} //namespace paradogs

} //namespace libp
