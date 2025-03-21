AC_PREREQ([2.69])
AC_INIT([neko],[0.9.99])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AM_MAINTAINER_MODE
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST
AC_LANG(Fortran)
AC_PROG_INSTALL
AC_PROG_SED
AC_PROG_MAKE_SET
LT_INIT([disable-shared])
LT_LANG([Fortran])

AC_ARG_ENABLE(contrib,
              AS_HELP_STRING([--enable-contrib],[Compile various tools]), 
              [enable_contrib=${enableval}], [enable_contrib=yes])

AC_ARG_ENABLE(real,
	      AS_HELP_STRING([--enable-real=Xp  Specify working precision of REAL types:
                    ssp -- REAL(kind=4) xp -- sp 
                    sp -- REAL(kind=4) xp -- dp
                    dp -- REAL(kind=8) xp -- dp (default)
                    qp -- REAL(kind=16) xp -- qp]),
	        [enable_real=${enableval}], [enable_real=dp])

AC_ARG_ENABLE(blk_size,
              AS_HELP_STRING([--enable-blk_size=number of values Specify size of fast memory]),
              [blk_size=${enableval}], [blk_size="10000"])

AC_ARG_ENABLE(device-mpi,
	      AS_HELP_STRING([--enable-device-mpi],[Enable device aware MPI]),
	      [enable_device_mpi=${enableval}], [enable_device_mpi=no])

AC_ARG_ENABLE(openmp, AC_HELP_STRING([--enable-openmp], [Enable OpenMP]),
              [enable_openmp=${enableval}], [enable_openmp=no])

# Test for a sane fortran environment (^-^)
AC_LANG(Fortran)
AC_PROG_FC(,90)
AX_COARRAY
AX_DTYPE_IO
AX_REAL128

# Test for a working MPI compiler
AX_MPI([have_mpi=yes], [have_mpi=no])
if test "x${have_mpi}" != xno; then
   FC="$MPIFC"
   LIBS="$MPILIBS $LIBS"
   AC_LANG(C)
   AX_MPI([have_mpi=yes],[have_mpi=no])
   if  test "x${have_mpi}" != xno; then
       CC="$MPICC"
   else
       AC_MSG_ERROR([Can't find a suitable C MPI compiler])
   fi
   AC_LANG(C++)
   AX_MPI([have_mpi=yes],[have_mpi=no])
   if  test "x${have_mpi}" != xno; then
       CXX="$MPICXX"
   fi
   AC_LANG(Fortran)
else
   AC_MSG_ERROR([Can't find a suitable MPI compiler])   
fi
AX_MPIF08
AX_MPI_DTYPE

# Test for OpenMP (if requested)
if test "x${enable_openmp}" != xno; then
   AX_OPENMP([have_openmp=yes], [have_openmp=no])
   if test "x${have_openmp}" != xno; then
      FCFLAGS="$FCFLAGS $OPENMP_FCFLAGS"
      AC_LANG(C)
      AX_OPENMP([have_openmp=yes], [have_openmp=no])
      if test "x${have_openmp}" != xno; then
         CFLAGS="$CFLAGS $OPENMP_CFLAGS"
      fi
      AC_LANG(Fortran)
   fi
fi

AC_FC_PP_SRCEXT([F90])
AX_COMPILER_VENDOR

        
              
# Checks for platform headers
AC_LANG_PUSH([C])
AC_CHECK_HEADERS([cpuid.h])
AC_LANG_POP([C])  

# Set necessary flags to fix issues with various compilers
if test "x$ax_cv_fc_compiler_vendor" = xportland; then
   FCFLAGS="$FCFLAGS  -Mbackslash"
fi

if test "x$ax_cv_fc_compiler_vendor" = xnvhpc; then
   FCFLAGS="$FCFLAGS  -Mbackslash"
fi

# Store build information (date, host, FC)
AC_SUBST(NEKO_BUILD_INFO,
        ["(build: $(date +%Y-%m-%d) on $host using $ax_cv_fc_compiler_vendor)"], [])

# Checks for Doxygen
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -n "$DOXYGEN"; then
  AC_CHECK_PROGS([DOXYGEN_DEPS], [dot])
fi
AM_CONDITIONAL([ENABLE_DOXYGEN], [test -n "$DOXYGEN_DEPS"])

# Checks for makedepf90
AC_CHECK_PROGS([MAKEDEPF90], [makedepf90])
AM_CONDITIONAL([ENABLE_MAKEDEPF90], [test -n "$MAKEDEPF90"])

# Checks for flinter
AC_CHECK_PROGS([FLINT], [flint])
AM_CONDITIONAL([ENABLE_FLINT], [test -n "$FLINT"])

# Checks for various platforms
AX_CRAY
AX_HPE_CRAY
AX_SX

# Checks for pFUnit
AX_PFUNIT

# Check for blas and lapack
AS_IF([test "x${is_cray}" = xyes || test "x${is_hpe_cray}" = xyes],
       AX_CRAY_LIBSCI
       if test "x${have_cray_libsci}" = xyes; then
         have_blas_lapack="yes"
       else
         have_blas_lapack="no"
       fi, have_blas_lapack="no")

if test "x${have_blas_lapack}" != xyes; then
  AX_LAPACK
  if test "x${ax_lapack_ok}" != xno; then
     LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS"
     have_blas_lapack=yes
  else
     have_blas_lapack=no
  fi
fi

if test "x${have_blas_lapack}" != xyes; then
   AC_MSG_ERROR([Can't find a suitable BLAS/LAPACK library])
fi

if test "x${have_darshan}" != xyes; then
    AX_DARSHAN
    if test "x${ax_darshan_ok}" != xno; then
        LIBS="$DARSHAN_LIBS $LIBS"
        have_darshan=yes
    else
        have_darshan=no
    fi
fi

found_domain_decomp="no"

# Checks for ParMetis
if test "x${is_cray}" = xyes; then
   AX_CRAY_PARMETIS
   if test "x${have_cray_parmetis}" = xyes; then
      found_domain_decomp="yes"
   fi
fi

if test "x${found_domain_decomp}" != xyes; then
   AX_PARMETIS
   if test "${PARMETIS_LIBS}"; then
       LIBS="$PARMETIS_LIBS $LIBS"
       found_domain_decomp="yes"
   fi
fi

# Checks for ADIOS2
AX_ADIOS2

# Checks for gslib
AX_GSLIB

# Select Working precision of REAL
if test "x$enable_real" = xssp; then
   AC_SUBST(NEKO_REAL_TYPE, "sp")
   AC_SUBST(NEKO_C_REAL_TYPE, "c_float")
   AC_SUBST(NEKO_MPI_REAL_TYPE, "MPI_REAL")
   AC_SUBST(NEKO_DEV_REAL_TYPE, "float")
   AC_SUBST(NEKO_XREAL_TYPE, "sp")
   AC_SUBST(NEKO_C_XREAL_TYPE, "c_float")
   AC_SUBST(NEKO_MPI_XREAL_TYPE, "MPI_REAL")
   AC_SUBST(NEKO_DEV_XREAL_TYPE, "float")
elif test "x$enable_real" = xsp; then
   AC_SUBST(NEKO_REAL_TYPE, "sp")
   AC_SUBST(NEKO_C_REAL_TYPE, "c_float")
   AC_SUBST(NEKO_MPI_REAL_TYPE, "MPI_REAL")
   AC_SUBST(NEKO_DEV_REAL_TYPE, "float")
   AC_SUBST(NEKO_XREAL_TYPE, "dp")
   AC_SUBST(NEKO_C_XREAL_TYPE, "c_double")
   AC_SUBST(NEKO_MPI_XREAL_TYPE, "MPI_DOUBLE_PRECISION")
   AC_SUBST(NEKO_DEV_XREAL_TYPE, "double") 
elif test "x$enable_real" = xdp; then
   AC_SUBST(NEKO_REAL_TYPE, "dp")
   AC_SUBST(NEKO_C_REAL_TYPE, "c_double")
   AC_SUBST(NEKO_MPI_REAL_TYPE, "MPI_DOUBLE_PRECISION")
   AC_SUBST(NEKO_DEV_REAL_TYPE, "double")
   AC_SUBST(NEKO_XREAL_TYPE, "dp")
   AC_SUBST(NEKO_C_XREAL_TYPE, "c_double")
   AC_SUBST(NEKO_MPI_XREAL_TYPE, "MPI_DOUBLE_PRECISION")
   AC_SUBST(NEKO_DEV_XREAL_TYPE, "double")   				
elif test "x$enable_real" = xqp; then
   AC_SUBST(NEKO_REAL_TYPE, "qp")
   AC_SUBST(NEKO_C_REAL_TYPE, "c_long_double")
   AC_SUBST(NEKO_MPI_REAL_TYPE, "MPI_REAL16")
   AC_SUBST(NEKO_DEV_REAL_TYPE, "long double")
   AC_SUBST(NEKO_XREAL_TYPE, "qp")
   AC_SUBST(NEKO_C_XREAL_TYPE, "c_long_double")
   AC_SUBST(NEKO_MPI_XREAL_TYPE, "MPI_REAL16")
   AC_SUBST(NEKO_DEV_XREAL_TYPE, "long double")
else
   AC_MSG_ERROR([Invalid working precision])
fi

if test "x$have_real128" = xyes; then
   AC_SUBST(NEKO_REAL128_TYPE, "REAL128")
else
   AC_SUBST(NEKO_REAL128_TYPE, "selected_real_kind(1)")
   if test "x$enable_real" = xqp; then
      AC_MSG_ERROR([REAL128 not supported])
   fi
fi                          

AC_SUBST(blk_size, ${blk_size})

# Checks for libxsmm backend
AX_LIBXSMM

# Checks for HIP backend
if test "x${is_hpe_cray}" = xyes; then
   AX_CRAY_ROCM
else
   AX_HIP
fi

# Checks for CUDA backend
if test "x${is_cray}" = xyes; then
  AX_CRAY_CUDATOOLKIT
else
  AX_CUDA
fi

# Checks for OpenCL backend
AX_OPENCL

# Check if device aware MPI is requested
if (test "x${have_cuda}" = xyes ||
    test "x${have_hip}" = xyes  ||
    test "x${have_opencl}" = xyes); then
   if test "x${enable_device_mpi}" = xyes; then

      # Check that necessary accelerator modules are loaded
      if (test "x${is_cray}" = xyes || test "x${is_hpe_cray}" = xyes); then    
      	 AX_CRAY_ACCEL
      fi

      AC_SUBST(device_mpi, .true.)
      AC_DEFINE(HAVE_DEVICE_MPI, 1, [Define if device MPI is requested.])
      
      if test "x${have_cuda}" = xyes; then
         CUDA_CFLAGS="$CUDA_CFLAGS -DHAVE_DEVICE_MPI=1"

         # Make sure we have support for atomicAdd
      	 AS_IF([test "$CUDA_ARCH"],[],[CUDA_ARCH="-arch sm_60"])         
      fi      

      if  test "x${have_hip}" = xyes; then
         HIP_HIPCC_FLAGS="$HIP_HIPCC_FLAGS -DHAVE_DEVICE_MPI=1"
      fi
   else
      AC_SUBST(device_mpi, .false.)
   fi
else
   AC_SUBST(device_mpi, .false.)
fi

# Checks for NCCL backend
AX_NCCL

# Checks for RCCL backend
AX_RCCL

# Checks for NVSHMEM backend
AX_NVSHMEM

# Checks for NVTX profiling
AX_NVTX

# Checks for ROCTX profiling
AX_ROCTX

# Checks for HDF5
if (test "x${is_cray}" = xyes || test "x${is_hpe_cray}" = xyes); then
   AX_CRAY_HDF5_PARALLEL
else
   AX_HDF5
fi


# JSON-Fortran
PKG_CHECK_MODULES([JSON_Fortran],[json-fortran >= 7.1.0])
LIBS="$JSON_Fortran_LIBS $LIBS"
FCFLAGS="$JSON_Fortran_CFLAGS $FCFLAGS"

AM_CONDITIONAL([ENABLE_PARMETIS], [test "x$found_domain_decomp" = xyes])
AM_CONDITIONAL([ENABLE_ADIOS2], [test "x${have_adios2}" = xyes])
AM_CONDITIONAL([ENABLE_CUDA], [test "x${have_cuda}" = xyes])
AM_CONDITIONAL([ENABLE_HIP], [test "x${have_hip}" = xyes])
AM_CONDITIONAL([ENABLE_OPENCL], [test "x${have_opencl}" = xyes])
AM_CONDITIONAL([ENABLE_CONTRIB], [test "x${enable_contrib}" = xyes])
AM_CONDITIONAL([ENABLE_NVSHMEM], [test "x${have_nvshmem}" = xyes])

# Set device dependent flags
AC_SUBST(CUDA_ARCH)
AC_SUBST(CUDA_CFLAGS)
AC_SUBST(HIP_HIPCC_FLAGS)

# Set device backend if requested (and found)
if (test "x${have_cuda}" = xyes ||
    test "x${have_hip}" = xyes  ||
    test "x${have_opencl}" = xyes); then
    AC_SUBST(device_bcknd, "1")
else
    AC_SUBST(device_bcknd, "0")
fi

AC_CONFIG_FILES([Makefile\
		 src/Makefile\
		 tests/unit/Makefile\
		 src/comm/comm.F90\
	 	 src/config/neko_config.f90\
		 src/config/num_types.f90\
		 src/device/device_config.h\
		 src/scripts/cltostring.sh\
		 makeneko\
		 neko.pc])

# Config tests/unit
AC_CONFIG_FILES([tests/unit/stack/Makefile\
		 tests/unit/tuple/Makefile\
 		 tests/unit/point/Makefile\
		 tests/unit/htable/Makefile\
		 tests/unit/uset/Makefile\
		 tests/unit/quad/Makefile\
		 tests/unit/hex/Makefile\
		 tests/unit/tet/Makefile\
		 tests/unit/tri/Makefile\
		 tests/unit/math/Makefile\
		 tests/unit/mesh/Makefile\
		 tests/unit/space/Makefile\
		 tests/unit/dofmap/Makefile\
		 tests/unit/field/Makefile\
		 tests/unit/mean_field/Makefile\
		 tests/unit/mean_sqr_field/Makefile\
		 tests/unit/gather_scatter/Makefile\
		 tests/unit/krylov/Makefile\
		 tests/unit/bc/Makefile\
		 tests/unit/signal/Makefile\
		 tests/unit/jobctrl/Makefile\
		 tests/unit/device/Makefile\
		 tests/unit/device_math/Makefile\
		 tests/unit/mathops/Makefile\
		 tests/unit/device_mathops/Makefile\
		 tests/unit/octree/Makefile\
                 tests/unit/vector/Makefile\
                 tests/unit/matrix/Makefile\
                 tests/unit/scratch_registry/Makefile\
                 tests/unit/fast3d/Makefile\
                 tests/unit/time_scheme/Makefile\
                 tests/unit/time_scheme_controller/Makefile\
                 tests/unit/runge_kutta_scheme/Makefile\
                 tests/unit/time_based_controller/Makefile\
                 tests/unit/point_interpolation/Makefile])

if test "x${enable_contrib}" = xyes; then
   AC_CONFIG_FILES([contrib/Makefile\
		    contrib/rea2nbin/Makefile\
		    contrib/genmeshbox/Makefile\
		    contrib/average_fields_in_time/Makefile\
		    contrib/calc_lift_from_field/Makefile\
		    contrib/postprocess_fluid_stats/Makefile\
		    contrib/average_field_in_space/Makefile\
		    contrib/map_to_equidistant_1d/Makefile\
		    contrib/mesh_checker/Makefile\
		    contrib/prepart/Makefile])
fi
# Doxygen
AC_CONFIG_FILES([doc/Doxyfile doc/Makefile])

AC_OUTPUT
