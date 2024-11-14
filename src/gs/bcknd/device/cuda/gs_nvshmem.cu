/*
 Copyright (c) 2024, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

   * Redistributions in binary form must reproduce the above
     copyright notice, this list of conditions and the following
     disclaimer in the documentation and/or other materials provided
     with the distribution.

   * Neither the name of the authors nor the names of its
     contributors may be used to endorse or promote products derived
     from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
*/

#include <device/device_config.h>
#include <device/cuda/check.h>
#include <mpi.h>
#include <nvshmem.h>
#include <nvshmemx.h>
#include <comm/comm.h>
#include "gs_kernels.h"
#include "gs_nvshmem_kernels.h"

extern "C" {

  void cudamalloc_nvshmem(void** ptr, size_t size)
  {
    *ptr = nvshmem_malloc(size);
    CUDA_CHECK(cudaGetLastError());
    cudaMemset(*ptr, 0, size);
    CUDA_CHECK(cudaGetLastError());
  }    
  
  void cudafree_nvshmem(void** ptr, size_t size)
  {
    nvshmem_free(*ptr);
    CUDA_CHECK(cudaGetLastError());
  }
  
  void cuda_gs_pack_and_push(void *u_d, void *sbuf_d, void *sdof_d,
                             int soffset, int n, cudaStream_t stream,
                             int srank,  void *rbuf_d, int roffset, int* remote_offset,
			     int rrank, uint64_t counter, void* notifyDone, void* notifyReady,
			     int iter)
  {
    
    if(remote_offset[iter-1] == -1)
    {
      MPI_Sendrecv(&roffset, 1, MPI_INT,
                   rrank, 0,
                   &(remote_offset[iter-1]), 1, MPI_INT,
                   srank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    const int nthrds = 1024;
    const int nblcks = (n+nthrds-1)/nthrds;
      
    // TO DO investigate merging following 2 kernels (and also unpack).  
    gs_pack_kernel<real>
      <<<nblcks, nthrds, 0, stream>>>((real *) u_d, (real *) sbuf_d + soffset,
                                      (int *) sdof_d + soffset, n);
    
    pushShmemKernel<real>
      <<<nblcks,nthrds,0,stream>>>((real *) rbuf_d + remote_offset[iter-1],
                                   (real *) sbuf_d + soffset,
                                   (int *) sdof_d + soffset,
                                   srank, rrank, n, counter,
                                   (uint64_t*) notifyDone,
                                   (uint64_t*) notifyReady);
    CUDA_CHECK(cudaGetLastError());
  }

  void cuda_gs_pack_and_push_wait(cudaStream_t stream, int counter, void* notifyDone)
  {
    pushShmemKernelWait<<<1,1,0,stream>>>(counter,(uint64_t*) notifyDone);
    CUDA_CHECK(cudaGetLastError());
  }
}