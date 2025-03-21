/*
 Copyright (c) 2021-2025, The Neko Authors
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

#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "math_kernel.h"

extern "C" {

#include <math/bcknd/device/device_mpi_reduce.h>
#include <math/bcknd/device/device_mpi_op.h>

#ifdef HAVE_RCCL
#include <math/bcknd/device/device_nccl_reduce.h>
#include <math/bcknd/device/device_nccl_op.h>
#endif
  
  /** Fortran wrapper for copy
   * Copy a vector \f$ a = b \f$
   */
  void hip_copy(void *a, void *b, int *n) {
    HIP_CHECK(hipMemcpyAsync(a, b, (*n) * sizeof(real),
                             hipMemcpyDeviceToDevice,
                             (hipStream_t) glb_cmd_queue));
  }

  /** Fortran wrapper for masked copy
   * Copy a vector \f$ a(mask) = b(mask) \f$
   */
  void hip_masked_copy(void *a, void *b, void *mask, int *n, int *m) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(masked_copy_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b, (int *) mask, *n, *m);

    HIP_CHECK(hipGetLastError());

  }
  
 
  /** Fortran wrapper for masked reduced copy
   * Copy a vector \f$ a = b(mask) \f$
   */
  void hip_masked_red_copy(void *a, void *b, void *mask, int *n, int *m) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(masked_red_copy_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b, (int *) mask, *n, *m);

    HIP_CHECK(hipGetLastError());

  } 

  /** Fortran wrapper for masked atomic reduction
   * update a vector \f$ a += b(mask) \f$ where mask is not unique
   */
  void hip_masked_atomic_reduction(void *a, void *b, void *mask, int *n, int *m) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(masked_atomic_reduction_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b, (int *) mask, *n, *m);

    HIP_CHECK(hipGetLastError());

  } 

  /** Fortran wrapper for cfill_mask
   * Fill a scalar to vector \f$ a_i = s, for i \in mask \f$
   */
  void hip_cfill_mask(void* a, real* c, int* size, void* mask, int* mask_size) {

      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks(((*mask_size) + 1024 - 1) / 1024, 1, 1);
      
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cfill_mask_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real*)a, *c, *size, (int*)mask, *mask_size);

    HIP_CHECK(hipGetLastError());
  }

  /** Fortran wrapper for rzero
   * Zero a real vector
   */
  void hip_rzero(void *a, int *n) {
    HIP_CHECK(hipMemsetAsync(a, 0, (*n) * sizeof(real),
                             (hipStream_t) glb_cmd_queue));
  }

  /** Fortran wrapper for cmult
   * Multiplication by constant c \f$ a = c \cdot a \f$
   */
  void hip_cmult(void *a, real *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cmult_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, *c, *n);
    HIP_CHECK(hipGetLastError());
    
  }

  /** Fortran wrapper for cmult
   * Multiplication by constant c \f$ a = c \cdot b \f$
   */
  void hip_cmult2(void *a, void *b, real *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cmult2_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a,(real *) b, *c, *n);
    HIP_CHECK(hipGetLastError());
    
  }
  /** Fortran wrapper for cadd
   * Add a scalar to vector \f$ a = \sum a_i + s \f$
   */
  void hip_cadd(void *a, real *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cadd_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, *c, *n);
    HIP_CHECK(hipGetLastError()); 
  }

  /**
   * Fortran wrapper for cadd2
   * Add a scalar to a vector \f$ a = b + s \f$
   */
  void hip_cadd2(void *a, void *b, real *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(cadd2_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b, *c, *n);
    HIP_CHECK(hipGetLastError());
  }

  /** Fortran wrapper for cfill
   * Multiplication by constant c \f$ a = c \cdot a \f$
   */
  void hip_cfill(void *a, real *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    if (*n > 0) { 
      hipLaunchKernelGGL(HIP_KERNEL_NAME(cfill_kernel<real>),
                         nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                         (real *) a, *c, *n);
      HIP_CHECK(hipGetLastError());
    }
  }

  /**
   * Fortran wrapper for add2
   * Vector addition \f$ a = a + b \f$
   */
  void hip_add2(void *a, void *b, int *n) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(add2_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b, *n);
    HIP_CHECK(hipGetLastError());
  }
  
  /**
   * Fortran wrapper for add4
   * Vector addition \f$ a = b + c + d\f$
   */
  void hip_add4(void *a, void *b, void *c, void *d, int *n) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(add4_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b, (real *) c, (real *) d, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for add2s1
   * Vector addition with scalar multiplication \f$ a = c_1 a + b \f$
   * (multiplication on first argument) 
   */
  void hip_add2s1(void *a, void *b, real *c1, int *n) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(add2s1_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b,
                       *c1, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for add2s2
   * Vector addition with scalar multiplication \f$ a = a + c_1 b \f$
   * (multiplication on second argument) 
   */
  void hip_add2s2(void *a, void *b, real *c1, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(add2s2_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b,
                       *c1, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for addsqr2s2
   * Vector addition with scalar multiplication \f$ a = a + c_1 (b * b) \f$
   * (multiplication on second argument) 
   */
  void hip_addsqr2s2(void *a, void *b, real *c1, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(addsqr2s2_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b,
                       *c1, *n);
    HIP_CHECK(hipGetLastError());
  }
  
  /**
   * Fortran wrapper for add3s2
   * Vector addition with scalar multiplication \f$ a = c_1 b + c_2 c \f$
   * (multiplication on second argument) 
   */
  void hip_add3s2(void *a, void *b, void *c, real *c1, real *c2, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(add3s2_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b, (real *) c,
                       *c1, *c2, *n);
    HIP_CHECK(hipGetLastError());
  }
  
  /**
   * Fortran wrapper for invcol1
   * Invert a vector \f$ a = 1 / a \f$
   */
  void hip_invcol1(void *a, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(invcol1_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for invcol2
   * Vector division \f$ a = a / b \f$
   */
  void hip_invcol2(void *a, void *b, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(invcol2_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b, *n);
    HIP_CHECK(hipGetLastError());
  }
  
  /**
   * Fortran wrapper for col2
   * Vector multiplication with 2 vectors \f$ a = a \cdot b \f$
   */
  void hip_col2(void *a, void *b, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(col2_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b, *n);
    HIP_CHECK(hipGetLastError());
  }
  
  /**
   * Fortran wrapper for col3
   * Vector multiplication with 3 vectors \f$ a = b \cdot c \f$
   */
  void hip_col3(void *a, void *b, void *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(col3_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b, (real *) c, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for subcol3
   * Vector multiplication with 3 vectors \f$ a = b \cdot c \f$
   */
  void hip_subcol3(void *a, void *b, void *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(subcol3_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b, (real *) c, *n);
    HIP_CHECK(hipGetLastError());
  }
  
  /**
   * Fortran wrapper for sub2
   * Vector subtraction \f$ a = a - b \f$
   */
  void hip_sub2(void *a, void *b, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(sub2_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for sub3
   * Vector subtraction \f$ a = b - c \f$
   */
  void hip_sub3(void *a, void *b, void *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(sub3_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b, (real *) c, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for addcol3
   * \f$ a = a + b * c \f$
   */
  void hip_addcol3(void *a, void *b, void *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(addcol3_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b, (real *) c, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for addcol4
   * \f$ a = a + b * c * d \f$
   */
  void hip_addcol4(void *a, void *b, void *c, void *d, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(addcol4_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b, (real *) c, (real *) d, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for vdot3
   * \f$ dot = u \cdot v \f$
   */
  void hip_vdot3(void *dot, void *u1, void *u2, void *u3,
                 void *v1, void *v2, void *v3, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(vdot3_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) dot, (real *) u1, (real *) u2, (real *) u3,
                       (real *) v1, (real *) v2, (real *) v3, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper for vdot3
   * \f$ u = v \times w \f$
   */
  void hip_vcross(void *u1, void *u2, void *u3,
                 void *v1, void *v2, void *v3,
                 void *w1, void *w2, void *w3, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(vcross_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) u1, (real *) u2, (real *) u3,
                       (real *) v1, (real *) v2, (real *) v3, 
                       (real *) w1, (real *) w2, (real *) w3, *n);
    HIP_CHECK(hipGetLastError());
  }


  /*
   * Reduction buffer
   */
  int red_s = 0;
  real * bufred = NULL;
  real * bufred_d = NULL;
  
  /**
   * Fortran wrapper vlsc3
   * Compute multiplication sum \f$ dot = u \cdot v \cdot w \f$
   */
  real hip_vlsc3(void *u, void *v, void *w, int *n) {
        
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;
    
    if ( nb > red_s){
      red_s = nb;
      if (bufred != NULL) {
        HIP_CHECK(hipHostFree(bufred));
        HIP_CHECK(hipFree(bufred_d));        
      }
      HIP_CHECK(hipHostMalloc(&bufred,nb*sizeof(real),hipHostMallocDefault));
      HIP_CHECK(hipMalloc(&bufred_d, nb*sizeof(real)));
    }
     
    hipLaunchKernelGGL(HIP_KERNEL_NAME(glsc3_kernel<real>),
                       nblcks, nthrds, 0, stream,
                       (real *) u, (real *) v,
                       (real *) w, bufred_d, *n);
    HIP_CHECK(hipGetLastError());
    hipLaunchKernelGGL(HIP_KERNEL_NAME(reduce_kernel<real>),
                       1, 1024, 0, stream, bufred_d, nb);
    HIP_CHECK(hipGetLastError());

    HIP_CHECK(hipMemcpyAsync(bufred, bufred_d, sizeof(real),
                             hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipStreamSynchronize(stream));

    return bufred[0];
  }
  

  /**
   * Fortran wrapper glsc3
   * Weighted inner product \f$ a^T b c \f$
   */
  real hip_glsc3(void *a, void *b, void *c, int *n) {
        
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;
    
    if ( nb > red_s){
      red_s = nb;
      if (bufred != NULL) {
        HIP_CHECK(hipHostFree(bufred));
        HIP_CHECK(hipFree(bufred_d));        
      }
      HIP_CHECK(hipHostMalloc(&bufred,nb*sizeof(real),hipHostMallocDefault));
      HIP_CHECK(hipMalloc(&bufred_d, nb*sizeof(real)));
    }
     
    hipLaunchKernelGGL(HIP_KERNEL_NAME(glsc3_kernel<real>),
                       nblcks, nthrds, 0, stream,
                       (real *) a, (real *) b,
                       (real *) c, bufred_d, *n);
    HIP_CHECK(hipGetLastError());
    hipLaunchKernelGGL(HIP_KERNEL_NAME(reduce_kernel<real>),
                       1, 1024, 0, stream, bufred_d, nb);
    HIP_CHECK(hipGetLastError());

#ifdef HAVE_RCCL
    device_nccl_allreduce(bufred_d, bufred_d, 1, sizeof(real),
                          DEVICE_NCCL_SUM, stream);
    HIP_CHECK(hipMemcpyAsync(bufred, bufred_d, sizeof(real),
                             hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipStreamSynchronize(stream));
#elif HAVE_DEVICE_MPI
    HIP_CHECK(hipStreamSynchronize(stream));
    device_mpi_allreduce(bufred_d, bufred, 1, sizeof(real), DEVICE_MPI_SUM);
#else
    HIP_CHECK(hipMemcpyAsync(bufred, bufred_d, sizeof(real),
                             hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipStreamSynchronize(stream));
#endif
    return bufred[0];
  }

  /**
   * Fortran wrapper for doing an reduction to an array
   * Weighted inner product \f$ w^T v(n,1:j) c \f$
   */
  void hip_glsc3_many(real *h, void * w, void *v,void *mult, int *j, int *n){ 
    int pow2 = 1;
    while(pow2 < (*j)){
      pow2 = 2*pow2;
    }
    const int nt = 1024/pow2;   
    const dim3 nthrds(pow2, nt, 1);
    const dim3 nblcks(((*n)+nt - 1)/nt, 1, 1);
    const dim3 nthrds_red(1024,1,1);
    const dim3 nblcks_red( (*j),1,1);
    const int nb = ((*n) + nt - 1)/nt;
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;
    
    if((*j)*nb>red_s){
      red_s = (*j)*nb;
      if (bufred != NULL) {
        HIP_CHECK(hipHostFree(bufred));
        HIP_CHECK(hipFree(bufred_d));
      }      
      HIP_CHECK(hipHostMalloc(&bufred,(*j)*nb*sizeof(real),hipHostMallocDefault));
      HIP_CHECK(hipMalloc(&bufred_d, (*j)*nb*sizeof(real)));
    }
    hipLaunchKernelGGL(HIP_KERNEL_NAME(glsc3_many_kernel<real>),
                       nblcks, nthrds, 0, stream,
                       (const real *) w, (const real **) v,
                       (const real *)mult, bufred_d, *j, *n);
    HIP_CHECK(hipGetLastError());
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(glsc3_reduce_kernel<real>),
                       nblcks_red, nthrds_red, 0, stream,
                       bufred_d, nb, *j);
    HIP_CHECK(hipGetLastError());

#ifdef HAVE_RCCL
    device_nccl_allreduce(bufred_d, bufred_d, (*j), sizeof(real),
                          DEVICE_NCCL_SUM, stream);
    HIP_CHECK(hipMemcpyAsync(h, bufred_d, (*j)* sizeof(real),
                             hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipStreamSynchronize(stream));
#elif HAVE_DEVICE_MPI
    HIP_CHECK(hipStreamSynchronize(stream));
    device_mpi_allreduce(bufred_d, h, (*j), sizeof(real), DEVICE_MPI_SUM);
#else
    HIP_CHECK(hipMemcpyAsync(h, bufred_d, (*j)* sizeof(real),
                             hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipStreamSynchronize(stream));
#endif    
  }

  /**
   * Fortran wrapper for add2s2
   * Vector addition with scalar multiplication 
   * \f$ x = x + c_1 p1 + c_2p2 + ... + c_jpj \f$
   * (multiplication on second argument) 
   */
  void hip_add2s2_many(void *x, void **p, void *alpha, int *j, int *n) {
        
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(add2s2_many_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) x, (const real **) p, (real *) alpha, *j, *n);
    HIP_CHECK(hipGetLastError());

  }

  /**
   * Fortran wrapper for add3
   * Vector addition \f$ a = b + c \f$
   */
  void hip_add3(void *a, void *b, void *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(add3_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, (real *) b, (real *) c, *n);
    HIP_CHECK(hipGetLastError());
  }

  /**
   * Fortran wrapper glsc2
   * Weighted inner product \f$ a^T b \f$
   */
  real hip_glsc2(void *a, void *b, int *n) {
        
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;
        
    if ( nb > red_s){
      red_s = nb;
      if (bufred != NULL) {
        HIP_CHECK(hipHostFree(bufred));
        HIP_CHECK(hipFree(bufred_d));        
      }
      HIP_CHECK(hipHostMalloc(&bufred,nb*sizeof(real),hipHostMallocDefault));
      HIP_CHECK(hipMalloc(&bufred_d, nb*sizeof(real)));
    }
     
    hipLaunchKernelGGL(HIP_KERNEL_NAME(glsc2_kernel<real>),
                       nblcks, nthrds, 0, stream,
                       (real *) a, (real *) b, bufred_d, *n);
    HIP_CHECK(hipGetLastError());
    hipLaunchKernelGGL(HIP_KERNEL_NAME(reduce_kernel<real>),
                       1, 1024, 0, stream, bufred_d, nb);
    HIP_CHECK(hipGetLastError());

#ifdef HAVE_RCCL
    device_nccl_allreduce(bufred_d, bufred_d, 1, sizeof(real),
                          DEVICE_NCCL_SUM, stream);
    HIP_CHECK(hipMemcpyAsync(bufred, bufred_d,sizeof(real),
                             hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipStreamSynchronize(stream));
#elif HAVE_DEVICE_MPI
    HIP_CHECK(hipStreamSynchronize(stream));
    device_mpi_allreduce(bufred_d, bufred, 1, sizeof(real), DEVICE_MPI_SUM);
#else
    HIP_CHECK(hipMemcpyAsync(bufred, bufred_d, sizeof(real),
                             hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipStreamSynchronize(stream));
#endif
    return bufred[0];
  }

  /** 
   * Fortran wrapper glsum
   * Sum a vector of length n
   */
  real hip_glsum(void *a, int *n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;
    
    if ( nb > red_s){
      red_s = nb;
      if (bufred != NULL) {
        HIP_CHECK(hipHostFree(bufred));
        HIP_CHECK(hipFree(bufred_d));        
      }
      HIP_CHECK(hipHostMalloc(&bufred,nb*sizeof(real),hipHostMallocDefault));
      HIP_CHECK(hipMalloc(&bufred_d, nb*sizeof(real)));
    }
    if( *n > 0) { 
      hipLaunchKernelGGL(HIP_KERNEL_NAME(glsum_kernel<real>),
                       nblcks, nthrds, 0, stream,
                       (real *) a, bufred_d, *n);
      HIP_CHECK(hipGetLastError());
      hipLaunchKernelGGL(HIP_KERNEL_NAME(reduce_kernel<real>),
                       1, 1024, 0, stream, bufred_d, nb);
      HIP_CHECK(hipGetLastError());
    }

#ifdef HAVE_RCCL
    device_nccl_allreduce(bufred_d, bufred_d, 1, sizeof(real),
                          DEVICE_NCCL_SUM, stream);
    HIP_CHECK(hipMemcpyAsync(bufred, bufred_d,sizeof(real),
                             hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipStreamSynchronize(stream));
#elif HAVE_DEVICE_MPI
    HIP_CHECK(hipStreamSynchronize(stream));
    device_mpi_allreduce(bufred_d, bufred, 1, sizeof(real), DEVICE_MPI_SUM);
#else
    HIP_CHECK(hipMemcpyAsync(bufred, bufred_d,sizeof(real),
                             hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipStreamSynchronize(stream));
#endif    
    return bufred[0];
  }

  /** Fortran wrapper for absval
   * \f$ a = abs(a) \f$
   */
  void hip_absval(void *a, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(absval_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a, *n);
    HIP_CHECK(hipGetLastError());
    
  }
}
