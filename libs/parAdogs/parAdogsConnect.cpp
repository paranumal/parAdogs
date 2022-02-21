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

typedef struct {
  hlong v[graph_t::MAX_NFACEVERTS]; // vertices on face
  hlong element, elementN;
  int face, faceN;    // face info
  int rank;

}parallelFace_t;


void graph_t::Connect(){

  /*Global number of elements*/
  hlong localNverts=static_cast<hlong>(Nverts);
  MPI_Allreduce(&localNverts, &gNVertsGlobal, 1, MPI_HLONG, MPI_SUM, gcomm);

  /*Get global element count offsets*/
  MPI_Scan(&localNverts, &gVoffsetU, 1, MPI_HLONG, MPI_SUM, gcomm);
  gVoffsetL = gVoffsetU-Nverts;

  /* build list of faces */
  libp::memory<parallelFace_t> faces(Nelements*Nfaces);

  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      const dlong id = f+Nfaces*e;

      for(int n=0;n<NfaceVerts;++n){
        dlong vid = faceVerts[f*NfaceVerts+n];
        faces[id].v[n] = elements[e].V[vid];
      }

      std::sort(faces[id].v, faces[id].v+NfaceVerts,
                std::less<hlong>());

      faces[id].element = e + gVoffsetL;
      faces[id].face = f;

      faces[id].elementN= -1;
      faces[id].faceN = -1;
    }
  }

  /* sort faces by their vertex number pairs */
  std::sort(faces.ptr(), faces.ptr()+Nelements*Nfaces,
            [&](const parallelFace_t& a, const parallelFace_t& b) {
              return std::lexicographical_compare(a.v, a.v+NfaceVerts,
                                                  b.v, b.v+NfaceVerts);
            });

  /* scan through sorted face lists looking for adjacent
     faces that have the same vertex ids */
  for(dlong n=0;n<Nelements*Nfaces-1;++n){
    if(std::equal(faces[n].v, faces[n].v+NfaceVerts,
                  faces[n+1].v)){
      // match
      faces[n].elementN = faces[n+1].element;
      faces[n].faceN = faces[n+1].face;

      faces[n+1].elementN = faces[n].element;
      faces[n+1].faceN = faces[n].face;
      ++n;
    }
  }

  /* resort faces back to the original element/face ordering */
  std::sort(faces.ptr(), faces.ptr()+Nelements*Nfaces,
            [](const parallelFace_t& a, const parallelFace_t& b) {
              if(a.element < b.element) return true;
              if(a.element > b.element) return false;

              return (a.face < b.face);
            });

  /* extract the element to element and element to face connectivity */

  // count # of elements to send to each rank based on
  // minimum {vertex id % gsize}
  libp::memory<int> Nsend(gsize, 0);
  libp::memory<int> Nrecv(gsize);
  libp::memory<int> sendOffsets(gsize);
  libp::memory<int> recvOffsets(gsize);

  int allNsend=0;
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      const dlong id = f+Nfaces*e;
      if (faces[id].elementN>-1) { /*matched face*/
        elements[e].E[f] = faces[id].elementN; //global id
        elements[e].F[f] = faces[id].faceN;
      } else { /*unmatched*/
        elements[e].E[f] = -1; //global id
        elements[e].F[f] = -1; /*mark face*/

        // find rank of destination for sorting based on min(face vertices)%gsize
        int destRank = static_cast<int>(faces[id].v[0]%gsize);

        // increment send gsize for
        ++Nsend[destRank];
        ++allNsend;
      }
    }
  }

  // find send offsets
  sendOffsets[0]=0;
  for(int rr=1;rr<gsize;++rr)
    sendOffsets[rr] = sendOffsets[rr-1] + Nsend[rr-1];

  // reset counters
  for(int rr=0;rr<gsize;++rr)
    Nsend[rr] = 0;

  // buffer for outgoing data
  libp::memory<parallelFace_t> sendFaces(allNsend);

  // Make the MPI_PARALLELFACE_T data type
  MPI_Datatype MPI_PARALLELFACE_T;
  MPI_Datatype dtype[6] = {MPI_HLONG, MPI_HLONG, MPI_HLONG, MPI_INT,
                            MPI_INT, MPI_INT};
  int blength[6] = {MAX_NFACEVERTS, 1, 1, 1, 1, 1};
  MPI_Aint addr[6], displ[6];
  MPI_Get_address ( &(sendFaces[0]              ), addr+0);
  MPI_Get_address ( &(sendFaces[0].element      ), addr+1);
  MPI_Get_address ( &(sendFaces[0].elementN     ), addr+2);
  MPI_Get_address ( &(sendFaces[0].face         ), addr+3);
  MPI_Get_address ( &(sendFaces[0].faceN        ), addr+4);
  MPI_Get_address ( &(sendFaces[0].rank         ), addr+5);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  displ[3] = addr[3] - addr[0];
  displ[4] = addr[4] - addr[0];
  displ[5] = addr[5] - addr[0];
  MPI_Type_create_struct (6, blength, displ, dtype, &MPI_PARALLELFACE_T);
  MPI_Type_commit (&MPI_PARALLELFACE_T);

  // pack face data
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      const dlong id = f+Nfaces*e;
      if (faces[id].elementN==-1) { /*unmatched face*/

        // find rank of destination for sorting based on min(face vertices)%gsize
        int destRank = static_cast<int>(faces[id].v[0]%gsize);

        // populate face to send out staged in segment of sendFaces array
        const int sid = sendOffsets[destRank]+Nsend[destRank];
        sendFaces[sid] = faces[id];
        sendFaces[sid].rank = grank;
        ++Nsend[destRank];
      }
    }
  }

  faces.free();

  // exchange byte counts
  MPI_Alltoall(Nsend.ptr(), 1, MPI_INT,
               Nrecv.ptr(), 1, MPI_INT,
               gcomm);

  // count incoming faces
  int allNrecv = 0;
  for(int rr=0;rr<gsize;++rr)
    allNrecv += Nrecv[rr];

  // find offsets for recv data
  recvOffsets[0]=0;
  for(int rr=1;rr<gsize;++rr)
    recvOffsets[rr] = recvOffsets[rr-1] + Nrecv[rr-1]; // byte offsets

  // buffer for incoming face data
  libp::memory<parallelFace_t> recvFaces(allNrecv);

  // exchange parallel faces
  MPI_Alltoallv(sendFaces.ptr(), Nsend.ptr(), sendOffsets.ptr(), MPI_PARALLELFACE_T,
                recvFaces.ptr(), Nrecv.ptr(), recvOffsets.ptr(), MPI_PARALLELFACE_T,
                gcomm);

  // local sort allNrecv received faces
  std::sort(recvFaces.ptr(), recvFaces.ptr()+allNrecv,
            [&](const parallelFace_t& a, const parallelFace_t& b) {
              return std::lexicographical_compare(a.v, a.v+NfaceVerts,
                                                  b.v, b.v+NfaceVerts);
            });

  // find matches
  for(int n=0;n<allNrecv-1;++n){
    // since vertices are ordered we just look for pairs
    if(std::equal(recvFaces[n].v, recvFaces[n].v+NfaceVerts,
                  recvFaces[n+1].v)){
      recvFaces[n].elementN = recvFaces[n+1].element;
      recvFaces[n].faceN = recvFaces[n+1].face;

      recvFaces[n+1].elementN = recvFaces[n].element;
      recvFaces[n+1].faceN = recvFaces[n].face;
      ++n;
    }
  }

  // sort back to original ordering
  std::sort(recvFaces.ptr(), recvFaces.ptr()+allNrecv,
            [](const parallelFace_t& a, const parallelFace_t& b) {
              if(a.rank < b.rank) return true;
              if(a.rank > b.rank) return false;

              if(a.element < b.element) return true;
              if(a.element > b.element) return false;

              return (a.face < b.face);
            });

  // send faces back from whence they came
  MPI_Alltoallv(recvFaces.ptr(), Nrecv.ptr(), recvOffsets.ptr(), MPI_PARALLELFACE_T,
                sendFaces.ptr(), Nsend.ptr(), sendOffsets.ptr(), MPI_PARALLELFACE_T,
                gcomm);

  // extract connectivity info
  for(int cnt=0;cnt<allNsend;++cnt){
    dlong e = static_cast<dlong>(sendFaces[cnt].element-gVoffsetL);
    hlong eN = sendFaces[cnt].elementN;
    int f = sendFaces[cnt].face;
    int fN = sendFaces[cnt].faceN;

    if(eN>=0 && fN>=0){ /*match found*/
      elements[e].E[f] = eN;
      elements[e].F[f] = fN;
    }
  }

  MPI_Barrier(gcomm);
  MPI_Type_free(&MPI_PARALLELFACE_T);
}

} //namespace paradogs

} //namespace libp
