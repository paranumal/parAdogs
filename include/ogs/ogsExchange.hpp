/*

The MIT License (MIT)

Copyright (c) 2017-2021 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#ifndef OGS_EXCHANGE_HPP
#define OGS_EXCHANGE_HPP

#include "ogs.hpp"
#include "ogs/ogsOperator.hpp"

namespace libp {

namespace ogs {

//virtual base class to perform MPI exchange of gatherScatter
class ogsExchange_t {
public:
  platform_t platform;
  MPI_Comm comm;
  int rank, size;

  dlong Nhalo, NhaloP;

  char* haloBuf;
  occa::memory o_haloBuf, h_haloBuf;

  static occa::stream dataStream;
  static occa::kernel extractKernel[4];

#ifdef GPU_AWARE_MPI
  bool gpu_aware=true;
#else
  bool gpu_aware=false;
#endif

  ogsExchange_t(platform_t &_platform, MPI_Comm _comm):
    platform(_platform),
    comm(_comm) {
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (!dataStream.isInitialized())
      dataStream = platform.device.createStream();
  }
  virtual ~ogsExchange_t() {}

  virtual void Start(const int k,
                     const Type type,
                     const Op op,
                     const Transpose trans,
                     const bool host=false)=0;
  virtual void Finish(const int k,
                      const Type type,
                      const Op op,
                      const Transpose trans,
                      const bool host=false)=0;

  virtual void AllocBuffer(size_t Nbytes)=0;

  friend void InitializeKernels(platform_t& platform, const Type type, const Op op);
};

//MPI communcation via single MPI_Alltoallv call
class ogsAllToAll_t: public ogsExchange_t {
private:

  dlong NsendN=0, NsendT=0;
  libp::memory<dlong> sendIdsN, sendIdsT;
  occa::memory o_sendIdsN, o_sendIdsT;

  ogsOperator_t postmpi;

  char* sendBuf;
  occa::memory o_sendBuf;
  occa::memory h_sendBuf;

  libp::memory<int> mpiSendCountsN;
  libp::memory<int> mpiSendCountsT;
  libp::memory<int> mpiRecvCountsN;
  libp::memory<int> mpiRecvCountsT;
  libp::memory<int> mpiSendOffsetsN;
  libp::memory<int> mpiSendOffsetsT;
  libp::memory<int> mpiRecvOffsetsN;
  libp::memory<int> mpiRecvOffsetsT;

  libp::memory<int> sendCounts;
  libp::memory<int> recvCounts;
  libp::memory<int> sendOffsets;
  libp::memory<int> recvOffsets;

public:
  ogsAllToAll_t(dlong Nshared,
               libp::memory<parallelNode_t> &sharedNodes,
               ogsOperator_t &gatherHalo,
               MPI_Comm _comm,
               platform_t &_platform);

  virtual void Start(const int k,
                     const Type type,
                     const Op op,
                     const Transpose trans,
                     const bool host=false);
  virtual void Finish(const int k,
                      const Type type,
                      const Op op,
                      const Transpose trans,
                      const bool host=false);

  virtual void AllocBuffer(size_t Nbytes);

};

//MPI communcation via pairwise send/recvs
class ogsPairwise_t: public ogsExchange_t {
private:

  dlong NsendN=0, NsendT=0;
  libp::memory<dlong> sendIdsN, sendIdsT;
  occa::memory o_sendIdsN, o_sendIdsT;

  ogsOperator_t postmpi;

  char* sendBuf;
  occa::memory o_sendBuf;
  occa::memory h_sendBuf;

  int NranksSendN=0, NranksRecvN=0;
  int NranksSendT=0, NranksRecvT=0;
  libp::memory<int> sendRanksN;
  libp::memory<int> sendRanksT;
  libp::memory<int> recvRanksN;
  libp::memory<int> recvRanksT;
  libp::memory<int> sendCountsN;
  libp::memory<int> sendCountsT;
  libp::memory<int> recvCountsN;
  libp::memory<int> recvCountsT;
  libp::memory<int> sendOffsetsN;
  libp::memory<int> sendOffsetsT;
  libp::memory<int> recvOffsetsN;
  libp::memory<int> recvOffsetsT;
  libp::memory<MPI_Request> requests;
  libp::memory<MPI_Status> statuses;

public:
  ogsPairwise_t(dlong Nshared,
               libp::memory<parallelNode_t> &sharedNodes,
               ogsOperator_t &gatherHalo,
               MPI_Comm _comm,
               platform_t &_platform);

  virtual void Start(const int k,
                     const Type type,
                     const Op op,
                     const Transpose trans,
                     const bool host=false);
  virtual void Finish(const int k,
                      const Type type,
                      const Op op,
                      const Transpose trans,
                      const bool host=false);

  virtual void AllocBuffer(size_t Nbytes);
};

//MPI communcation via Crystal Router
class ogsCrystalRouter_t: public ogsExchange_t {
private:

  struct crLevel {
    int Nmsg;
    int partner;

    int Nsend, Nrecv0, Nrecv1;
    dlong recvOffset;

    libp::memory<dlong> sendIds;
    occa::memory o_sendIds;

    ogsOperator_t gather;
  };

  int buf_id=0;
  occa::memory o_buf[2];
  occa::memory h_buf[2];
  char* buf[2];

  MPI_Request request[3];
  MPI_Status status[3];

  int Nlevels=0;
  libp::memory<crLevel> levelsN;
  libp::memory<crLevel> levelsT;

  int NsendMax=0, NrecvMax=0;
  char* sendBuf;
  char* recvBuf;
  occa::memory o_sendBuf;
  occa::memory h_sendBuf;
  occa::memory o_recvBuf;

public:
  ogsCrystalRouter_t(dlong Nshared,
                   libp::memory<parallelNode_t> &sharedNodes,
                   ogsOperator_t &gatherHalo,
                   MPI_Comm _comm,
                   platform_t &_platform);

  virtual void Start(const int k,
                     const Type type,
                     const Op op,
                     const Transpose trans,
                     const bool host=false);
  virtual void Finish(const int k,
                      const Type type,
                      const Op op,
                      const Transpose trans,
                      const bool host=false);

  virtual void AllocBuffer(size_t Nbytes);
};

} //namespace ogs

} //namespace libp

#endif
