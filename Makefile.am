pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = neko.pc

SUBDIRS = src

if ENABLE_CONTRIB
SUBDIRS += contrib
endif

if ENABLE_PFUNIT
SUBDIRS += tests
endif

if ENABLE_DOXYGEN
SUBDIRS += doc
endif


bin_SCRIPTS = makeneko
EXTRA_DIST = \
	makeneko.in\
	regen.sh\
	flinter_rc.yml\
	bench/ax/setup.f90\
	bench/ax/Makefile\
	bench/ax/axkernel.f90\
	bench/ax/driver.f90\
	bench/gs/Makefile\
	bench/gs/driver.f90\
	bench/tgv32/tgv.case\
	bench/tgv32/README\
	bench/tgv32/extract.sh\
	bench/tgv32/tgv.f90\
	examples/turb_channel/turb_channel.case\
	examples/turb_channel/README.md\
	examples/turb_channel/turb_channel.f90\
	examples/turb_channel/box.nmsh\
	examples/turbpipe/pipe.case\
	examples/turbpipe/pipe.f90\
	examples/turbpipe/pipe.nmsh\
	examples/turbpipe/README.md\
	examples/turbpipe/pipe_source.case\
	examples/hemi/hemi.nmsh\
	examples/hemi/README.md\
	examples/hemi/hemi.case\
	examples/lid/lid.f90\
	examples/lid/lid2d.nmsh\
	examples/lid/lid2d.f90\
	examples/lid/lid.nmsh\
	examples/lid/README.md\
	examples/lid/lid2d.box\
	examples/lid/lid2d.case\
	examples/lid/lid.case\
	examples/lid/lid.box\
	examples/poisson/ax_poisson.f90\
	examples/poisson/setup.f90\
	examples/poisson/Makefile\
	examples/poisson/README.md\
	examples/poisson/data/4096.nmsh\
	examples/poisson/data/16384.nmsh\
	examples/poisson/data/256.nmsh\
	examples/poisson/data/2048.nmsh\
	examples/poisson/data/128.nmsh\
	examples/poisson/data/32768.nmsh\
	examples/poisson/data/8192.nmsh\
	examples/poisson/data/512.nmsh\
	examples/poisson/data/1024.nmsh\
	examples/poisson/driver.f90\
	examples/2d_cylinder/ext_cyl.nmsh\
	examples/2d_cylinder/2d_cylinder.case\
	examples/2d_cylinder/README.md\
	examples/cyl_boundary_layer/cyl_bl_rot.case\
	examples/cyl_boundary_layer/cyl_bl.f90\
	examples/cyl_boundary_layer/README.md\
	examples/cyl_boundary_layer/cyl_bl_user.case\
	examples/cyl_boundary_layer/cyl_bl_basic.case\
	examples/cyl_boundary_layer/cyl.nmsh\
	examples/rayleigh-benard/rayleigh.f90\
	examples/rayleigh-benard/README.md\
	examples/rayleigh-benard/box.nmsh\
	examples/rayleigh-benard/rayleigh.case\
	examples/tgv/262144.nmsh\
	examples/tgv/tgv.case\
	examples/tgv/README.md\
	examples/tgv/32768.nmsh\
	examples/tgv/tgv.f90\
	examples/tgv/512.nmsh\
	examples/rayleigh-benard-cylinder/probes.csv\
	examples/rayleigh-benard-cylinder/rayleigh.f90\
	examples/rayleigh-benard-cylinder/README.md\
	examples/rayleigh-benard-cylinder/cylinder.nmsh\
	examples/rayleigh-benard-cylinder/rayleigh.case\
	examples/rayleigh-benard-cylinder/postprocess_probes.py\
	reframe/checks.py\
	reframe/settings.py\
	reframe/README.md\
	reframe/src/get-tgv-dns.sh\
	reframe/src/minihemi.case.template\
	reframe/src/rayleigh.case.template\
	reframe/src/rayleigh.f90\
	reframe/src/rocm_select_gpu_device\
	reframe/src/tgv.case.template\
	reframe/src/tgv.f90\
	patches/nvhpc_bge.patch

if ENABLE_FLINT
lint:
	$(FLINT) score -d 0 src/ -r flinter_rc.yml

check-local: lint

endif
