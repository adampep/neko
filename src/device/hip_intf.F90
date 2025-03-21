! Copyright (c) 2021-2022, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Fortran HIP interface
module hip_intf
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_HIP

  !> Enum @a hipError_t
  enum, bind(c)
    enumerator :: hipSuccess = 0
    enumerator :: hipErrorInvalidContext = 1
    enumerator :: hipErrorInvalidKernelFile = 2
    enumerator :: hipErrorMemoryAllocation = 3
    enumerator :: hipErrorInitializationError = 4
    enumerator :: hipErrorLaunchFailure = 5
    enumerator :: hipErrorLaunchOutOfResources = 6
    enumerator :: hipErrorInvalidDevice = 7
    enumerator :: hipErrorInvalidValue = 8
    enumerator :: hipErrorInvalidDevicePointer = 9
    enumerator :: hipErrorInvalidMemcpyDirection = 10
    enumerator :: hipErrorUnknown = 11
    enumerator :: hipErrorInvalidResourceHandle = 12
    enumerator :: hipErrorNotReady = 13
    enumerator :: hipErrorNoDevice = 14
    enumerator :: hipErrorPeerAccessAlreadyEnabled = 15
    enumerator :: hipErrorPeerAccessNotEnabled = 16
    enumerator :: hipErrorRuntimeMemory = 17
    enumerator :: hipErrorRuntimeOther = 18
    enumerator :: hipErrorHostMemoryAlreadyRegistered = 19
    enumerator :: hipErrorHostMemoryNotRegistered = 20
    enumerator :: hipErrorMapBufferObjectFailed = 21
    enumerator :: hipErrorTbd = 22
  end enum

  !> Enum @a hipMemcpyKind
  enum, bind(c)
    enumerator :: hipMemcpyHostToHost = 0
    enumerator :: hipMemcpyHostToDevice = 1
    enumerator :: hipMemcpyDeviceToHost = 2
    enumerator :: hipMemcpyDevicetoDevice = 3
    enumerator :: hipMemcpyDefault = 4
  end enum

  interface
     integer(c_int) function hipMalloc(ptr_d, s) &
          bind(c, name = 'hipMalloc')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: ptr_d
       integer(c_size_t), value :: s
     end function hipMalloc

     integer(c_int) function hipFree(ptr_d) &
          bind(c, name = 'hipFree')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr_d
     end function hipFree

     integer(c_int) function hipMemcpy(ptr_dst, ptr_src, s, dir) &
          bind(c, name = 'hipMemcpy')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr_dst, ptr_src
       integer(c_size_t), value :: s
       integer(c_int), value :: dir
     end function hipMemcpy

     integer(c_int) function hipMemcpyAsync(ptr_dst, ptr_src, s, dir, stream) &
          bind(c, name = 'hipMemcpyAsync')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr_dst, ptr_src, stream
       integer(c_size_t), value :: s
       integer(c_int), value :: dir
     end function hipMemcpyAsync

     integer(c_int) function hipDeviceSynchronize() &
          bind(c, name = 'hipDeviceSynchronize')
       use, intrinsic :: iso_c_binding
       implicit none
     end function hipDeviceSynchronize

     integer(c_int) function hipDeviceGetName(name, len, device) &
          bind(c, name = 'hipDeviceGetName')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: name
       integer(c_int), value :: len
       integer(c_int), value :: device
     end function hipDeviceGetName

     integer(c_int) function hipGetDeviceCount(amount) &
          bind(c, name = 'hipGetDeviceCount')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: amount
     end function hipGetDeviceCount

     integer(c_int) function hipStreamCreate(stream) &
          bind(c, name = 'hipStreamCreate')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: stream
     end function hipStreamCreate

     integer(c_int) function hipStreamCreateWithFlags(stream, flags) &
          bind(c, name = 'hipStreamCreateWithFlags')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: stream
       integer(c_int), value :: flags
     end function hipStreamCreateWithFlags

     integer(c_int) function hipStreamCreateWithPriority(stream, flags, prio) &
          bind(c, name = 'hipStreamCreateWithPriority')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: stream
       integer(c_int), value :: flags, prio
     end function hipStreamCreateWithPriority

     integer(c_int) function hipStreamDestroy(steam) &
          bind(c, name = 'hipStreamDestroy')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: steam
     end function hipStreamDestroy

     integer(c_int) function hipStreamSynchronize(stream) &
          bind(c, name = 'hipStreamSynchronize')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: stream
     end function hipStreamSynchronize

     integer(c_int) function hipStreamWaitEvent(stream, event, flags) &
          bind(c, name = 'hipStreamWaitEvent')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: stream, event
       integer(c_int), value :: flags
     end function hipStreamWaitEvent

     integer(c_int) function hipDeviceGetStreamPriorityRange(low_prio, &
                                                             high_prio) &
          bind(c, name = 'hipDeviceGetStreamPriorityRange')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: low_prio, high_prio
     end function hipDeviceGetStreamPriorityRange

     integer(c_int) function hipEventCreate(event) &
          bind(c, name = 'hipEventCreate')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: event
     end function hipEventCreate

     integer(c_int) function hipEventDestroy(event) &
          bind(c, name = 'hipEventDestroy')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: event
     end function hipEventDestroy

     integer(c_int) function hipEventCreateWithFlags(event, flags) &
          bind(c, name = 'hipEventCreateWithFlags')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: event
       integer(c_int), value :: flags
     end function hipEventCreateWithFlags

     integer(c_int) function hipEventRecord(event, stream) &
          bind(c, name = 'hipEventRecord')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: event, stream
     end function hipEventRecord

     integer(c_int) function hipEventSynchronize(event) &
          bind(c, name = 'hipEventSynchronize')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: event
     end function hipEventSynchronize
  end interface

contains

  subroutine hip_init(glb_cmd_queue, aux_cmd_queue, &
       STRM_HIGH_PRIO, STRM_LOW_PRIO)
    type(c_ptr), intent(inout) :: glb_cmd_queue
    type(c_ptr), intent(inout) :: aux_cmd_queue
    integer, intent(inout) :: STRM_HIGH_PRIO
    integer, intent(inout) :: STRM_LOW_PRIO

    if (hipDeviceGetStreamPriorityRange(STRM_LOW_PRIO, STRM_HIGH_PRIO) &
        .ne. hipSuccess) then
       call neko_error('Error retrieving stream priority range')
    end if

    if (hipStreamCreateWithPriority(glb_cmd_queue, 1, STRM_HIGH_PRIO) &
        .ne. hipSuccess) then
       call neko_error('Error creating main stream')
    end if

    if (hipStreamCreateWithPriority(aux_cmd_queue, 1, STRM_LOW_PRIO) &
        .ne. hipSuccess) then
       call neko_error('Error creating main stream')
    end if
  end subroutine hip_init

  subroutine hip_finalize(glb_cmd_queue, aux_cmd_queue)
    type(c_ptr), intent(inout) :: glb_cmd_queue
    type(c_ptr), intent(inout) :: aux_cmd_queue
    
    if (hipStreamDestroy(glb_cmd_queue) .ne. hipSuccess) then
       call neko_error('Error destroying main stream')
    end if

    if (hipStreamDestroy(aux_cmd_queue) .ne. hipSuccess) then
       call neko_error('Error destroying aux stream')
    end if
  end subroutine hip_finalize

  subroutine hip_device_name(name)
    character(len=*), intent(inout) :: name
    character(kind=c_char, len=1024), target :: c_name
    integer :: end_pos

    if (hipDeviceGetName(c_loc(c_name), 1024, 0) .ne. hipSuccess) then
       call neko_error('Failed to query device')
    end if

    end_pos = scan(c_name, C_NULL_CHAR)
    if (end_pos .ge. 2) then
       name(1:end_pos-1) = c_name(1:end_pos-1)
    end if

  end subroutine hip_device_name

  !> Return the number of available HIP devices
  integer function hip_device_count()
    type(c_ptr) :: hip_count_ptr
    integer :: amount

    if (hipGetDeviceCount(amount) .ne. hipSuccess) then
       call neko_error('Failed to query device count')
    end if

    hip_device_count = amount
  end function hip_device_count

#endif

end module hip_intf
