TOPDIR	= .
include ./Makefile.system

BLASDIRS = interface driver/level2 driver/level3 driver/others

ifneq ($(DYNAMIC_ARCH), 1)
BLASDIRS += kernel 
endif

ifdef UTEST_CHECK
SANITY_CHECK = 1
endif

ifdef SANITY_CHECK
BLASDIRS += reference
endif

ifndef PREFIX
PREFIX = /opt/OpenBLAS
endif

SUBDIRS	= $(BLASDIRS)
ifneq ($(NO_LAPACK), 1)
SUBDIRS	+= lapack
endif

SUBDIRS_ALL = $(SUBDIRS) test ctest utest exports benchmark ../laswp ../bench

.PHONY : all libs netlib test ctest shared install
.NOTPARALLEL : all libs prof lapack-test install

all :: libs netlib tests shared
	@echo
	@echo " OpenBLAS build complete."
	@echo
	@echo "  OS               ... $(OSNAME)             "
	@echo "  Architecture     ... $(ARCH)               "
ifndef BINARY64
	@echo "  BINARY           ... 32bit                 "
else
	@echo "  BINARY           ... 64bit                 "
endif
ifdef INTERFACE64
	@echo "  Use 64 bits int    (equivalent to \"-i8\" in Fortran)      "
endif
	@echo "  C compiler       ... $(C_COMPILER)  (command line : $(CC))"
	@echo "  Fortran compiler ... $(F_COMPILER)  (command line : $(FC))"
ifneq ($(OSNAME), AIX)
	@echo -n "  Library Name     ... $(LIBNAME)"
else
	@echo "  Library Name     ... $(LIBNAME)"
endif

ifndef SMP
	@echo " (Single threaded)  "
else
	@echo " (Multi threaded; Max num-threads is $(NUM_THREADS))"
endif

ifeq ($(USE_OPENMP), 1)
	@echo
	@echo " Use OpenMP in the multithreading. Becasue of ignoring OPENBLAS_NUM_THREADS and GOTO_NUM_THREADS flags, "
	@echo " you should use OMP_NUM_THREADS environment variable to control the number of threads."
	@echo
endif

ifeq ($(OSNAME), Darwin)
	@echo "WARNING: If you plan to use the dynamic library $(LIBDYNNAME), you must run:"
	@echo
	@echo "\"make PREFIX=/your_installation_path/ install\"."
	@echo
	@echo "(or set PREFIX in Makefile.rule and run make install."
	@echo "If you want to move the .dylib to a new location later, make sure you change"
	@echo "the internal name of the dylib with:"
	@echo
	@echo "install_name_tool -id /new/absolute/path/to/$(LIBDYNNAME) $(LIBDYNNAME)"
endif
	@echo
	@echo "To install the library, you can run \"make PREFIX=/path/to/your/installation install\"."
	@echo

shared :
ifndef NO_SHARED
ifeq ($(OSNAME), Linux)
	@$(MAKE) -C exports so
	@-ln -fs $(LIBSONAME) $(LIBPREFIX).so
	@-ln -fs $(LIBSONAME) $(LIBPREFIX).so.$(MAJOR_VERSION)
endif
ifeq ($(OSNAME), FreeBSD)
	@$(MAKE) -C exports so
	@-ln -fs $(LIBSONAME) $(LIBPREFIX).so
endif
ifeq ($(OSNAME), NetBSD)
	@$(MAKE) -C exports so
	@-ln -fs $(LIBSONAME) $(LIBPREFIX).so
endif
ifeq ($(OSNAME), Darwin)
	@$(MAKE) -C exports dyn
	@-ln -fs $(LIBDYNNAME) $(LIBPREFIX).dylib
endif
ifeq ($(OSNAME), WINNT)
	@$(MAKE) -C exports dll
endif
ifeq ($(OSNAME), CYGWIN_NT)
	@$(MAKE) -C exports dll
endif
endif

tests :
ifndef NOFORTRAN
ifndef TARGET
ifndef CROSS
	touch $(LIBNAME)
ifndef NO_FBLAS
	$(MAKE) -C test all
ifdef UTEST_CHECK
	$(MAKE) -C utest all
endif
endif
ifndef NO_CBLAS
	$(MAKE) -C ctest all
endif
endif
endif
endif

libs :
ifeq ($(CORE), UNKOWN)
	$(error OpenBLAS: Detecting CPU failed. Please set TARGET explicitly, e.g. make TARGET=your_cpu_target. Please read README for the detail.)
endif
ifeq ($(NOFORTRAN), 1)
	$(error OpenBLAS: Detecting fortran compiler failed. Please install fortran compiler, e.g. gfortran, ifort, openf90.)
endif
	@-ln -fs $(LIBNAME) $(LIBPREFIX).$(LIBSUFFIX)
	@for d in $(SUBDIRS) ; \
	do if test -d $$d; then \
	  $(MAKE) -C $$d $(@F) || exit 1 ; \
	fi; \
	done
#Save the config files for installation
	@cp Makefile.conf Makefile.conf_last
	@cp config.h config_last.h
ifdef QUAD_PRECISION
	@echo "#define QUAD_PRECISION">> config_last.h
endif
ifeq ($(EXPRECISION), 1)
	@echo "#define EXPRECISION">> config_last.h
endif
## 
ifeq ($(DYNAMIC_ARCH), 1)
	@$(MAKE) -C kernel commonlibs || exit 1
	@for d in $(DYNAMIC_CORE) ; \
	do  $(MAKE) GOTOBLAS_MAKEFILE= -C kernel TARGET_CORE=$$d kernel || exit 1 ;\
	done
	@echo DYNAMIC_ARCH=1 >> Makefile.conf_last
endif
ifdef USE_THREAD
	@echo USE_THREAD=$(USE_THREAD) >>  Makefile.conf_last
endif
	@touch lib.grd

prof : prof_blas prof_lapack

prof_blas :
	ln -fs $(LIBNAME_P) $(LIBPREFIX)_p.$(LIBSUFFIX)
	for d in $(SUBDIRS) ; \
	do if test -d $$d; then \
	  $(MAKE) -C $$d prof || exit 1 ; \
	fi; \
	done
ifeq ($(DYNAMIC_ARCH), 1)
	  $(MAKE) -C kernel commonprof || exit 1
endif

blas :
	ln -fs $(LIBNAME) $(LIBPREFIX).$(LIBSUFFIX)
	for d in $(BLASDIRS) ; \
	do if test -d $$d; then \
	  $(MAKE) -C $$d libs || exit 1 ; \
	fi; \
	done

hpl : 
	ln -fs $(LIBNAME) $(LIBPREFIX).$(LIBSUFFIX)
	for d in $(BLASDIRS) ../laswp exports ; \
	do if test -d $$d; then \
	  $(MAKE) -C $$d $(@F) || exit 1 ; \
	fi; \
	done
ifeq ($(DYNAMIC_ARCH), 1)
	  $(MAKE) -C kernel commonlibs || exit 1
	for d in $(DYNAMIC_CORE) ; \
	do  $(MAKE) GOTOBLAS_MAKEFILE= -C kernel TARGET_CORE=$$d kernel || exit 1 ;\
	done
endif

hpl_p :
	ln -fs $(LIBNAME_P) $(LIBPREFIX)_p.$(LIBSUFFIX)
	for d in $(SUBDIRS) ../laswp exports ; \
	do if test -d $$d; then \
	  $(MAKE) -C $$d $(@F) || exit 1 ; \
	fi; \
	done

ifeq ($(NO_LAPACK), 1)
netlib : 

else
netlib : lapack_prebuild
ifndef NOFORTRAN
	@$(MAKE) -C $(NETLIB_LAPACK_DIR) lapacklib
endif
ifndef NO_LAPACKE
	@$(MAKE) -C $(NETLIB_LAPACK_DIR) lapackelib
endif
endif

prof_lapack : lapack_prebuild
	@$(MAKE) -C $(NETLIB_LAPACK_DIR) lapack_prof

lapack_prebuild :
ifndef NOFORTRAN
	-@echo "FORTRAN     = $(FC)" > $(NETLIB_LAPACK_DIR)/make.inc
	-@echo "OPTS        = $(FFLAGS)" >> $(NETLIB_LAPACK_DIR)/make.inc
	-@echo "POPTS       = $(FPFLAGS)" >> $(NETLIB_LAPACK_DIR)/make.inc
	-@echo "NOOPT       = $(FFLAGS) -O0" >> $(NETLIB_LAPACK_DIR)/make.inc
	-@echo "PNOOPT      = $(FPFLAGS) -O0" >> $(NETLIB_LAPACK_DIR)/make.inc
	-@echo "LOADOPTS    = $(FFLAGS) $(EXTRALIB)" >> $(NETLIB_LAPACK_DIR)/make.inc
	-@echo "CC          = $(CC)" >> $(NETLIB_LAPACK_DIR)/make.inc
	-@echo "override CFLAGS      = $(LAPACK_CFLAGS)" >> $(NETLIB_LAPACK_DIR)/make.inc
	-@echo "ARCH        = $(AR)" >> $(NETLIB_LAPACK_DIR)/make.inc
	-@echo "ARCHFLAGS   = -ru" >> $(NETLIB_LAPACK_DIR)/make.inc
	-@echo "RANLIB      = $(RANLIB)" >> $(NETLIB_LAPACK_DIR)/make.inc
	-@echo "LAPACKLIB   = ../$(LIBNAME)" >> $(NETLIB_LAPACK_DIR)/make.inc
	-@echo "LAPACKELIB  = ../$(LIBNAME)" >> $(NETLIB_LAPACK_DIR)/make.inc
	-@echo "LAPACKLIB_P = ../$(LIBNAME_P)" >> $(NETLIB_LAPACK_DIR)/make.inc
	-@echo "SUFFIX      = $(SUFFIX)" >> $(NETLIB_LAPACK_DIR)/make.inc
	-@echo "PSUFFIX     = $(PSUFFIX)" >> $(NETLIB_LAPACK_DIR)/make.inc
	-@echo "CEXTRALIB   = $(EXTRALIB)" >> $(NETLIB_LAPACK_DIR)/make.inc
	-@cat  make.inc >> $(NETLIB_LAPACK_DIR)/make.inc
endif

lapack-3.4.2 : lapack-3.4.2.tgz
ifndef NOFORTRAN
ifndef NO_LAPACK
	@if test `$(MD5SUM) $< | $(AWK) '{print $$1}'` = 61bf1a8a4469d4bdb7604f5897179478; then \
		echo $(TAR) zxf $< ;\
		$(TAR) zxf $< && (cd $(NETLIB_LAPACK_DIR); $(PATCH) -p1 < ../patch.for_lapack-3.4.2) ;\
		rm -f $(NETLIB_LAPACK_DIR)/lapacke/make.inc ;\
	else \
		rm -rf $(NETLIB_LAPACK_DIR) ;\
		echo "	Cannot download lapack-3.4.2.tgz or the MD5 check sum is wrong (Please use orignal)."; \
		exit 1; \
	fi
endif
endif

LAPACK_URL=http://www.netlib.org/lapack/lapack-3.4.2.tgz

lapack-3.4.2.tgz :
ifndef NOFORTRAN
#http://stackoverflow.com/questions/7656425/makefile-ifeq-logical-or
ifeq ($(OSNAME), $(filter $(OSNAME),Darwin NetBSD))
	curl -O $(LAPACK_URL);
else
ifeq ($(OSNAME), FreeBSD)
	fetch $(LAPACK_URL);
else
	wget -O $@ $(LAPACK_URL);
endif
endif
endif

large.tgz : 
ifndef NOFORTRAN
	if [ ! -a $< ]; then
	-wget http://www.netlib.org/lapack/timing/large.tgz;
	fi
endif

timing.tgz :
ifndef NOFORTRAN
	if [ ! -a $< ]; then
	-wget http://www.netlib.org/lapack/timing/timing.tgz;
	fi
endif

lapack-timing : large.tgz timing.tgz
ifndef NOFORTRAN
	(cd $(NETLIB_LAPACK_DIR); $(TAR) zxf ../timing.tgz TIMING)
	(cd $(NETLIB_LAPACK_DIR)/TIMING; $(TAR) zxf ../../large.tgz )
	make -C $(NETLIB_LAPACK_DIR) tmglib
	make -C $(NETLIB_LAPACK_DIR)/TIMING
endif


lapack-test :
	$(MAKE) -C $(NETLIB_LAPACK_DIR) tmglib
	$(MAKE) -C $(NETLIB_LAPACK_DIR)/TESTING xeigtstc xeigtstd xeigtsts xeigtstz xlintstc xlintstd xlintstds xlintsts xlintstz xlintstzc
	@rm	-f $(NETLIB_LAPACK_DIR)/TESTING/*.out
	$(MAKE) -j 1 -C $(NETLIB_LAPACK_DIR)/TESTING
	$(GREP) failed $(NETLIB_LAPACK_DIR)/TESTING/*.out

dummy :

install :
	$(MAKE) -f Makefile.install install

clean ::
	@for d in $(SUBDIRS_ALL) ; \
	do if test -d $$d; then \
	  $(MAKE) -C $$d $(@F) || exit 1 ; \
	fi; \
	done
#ifdef DYNAMIC_ARCH
	@$(MAKE) -C kernel clean
#endif
	@$(MAKE) -C reference clean
	@rm -f *.$(LIBSUFFIX) *.so *~ *.exe getarch getarch_2nd *.dll *.lib *.$(SUFFIX) *.dwf $(LIBPREFIX).$(LIBSUFFIX) $(LIBPREFIX)_p.$(LIBSUFFIX) $(LIBPREFIX).so.$(MAJOR_VERSION) *.lnk myconfig.h
ifeq ($(OSNAME), Darwin)
	@rm -rf getarch.dSYM getarch_2nd.dSYM
endif
	@rm -f Makefile.conf config.h cblas_noconst.h Makefile_kernel.conf config_kernel.h st* *.dylib
	@touch $(NETLIB_LAPACK_DIR)/make.inc
	@$(MAKE) -C $(NETLIB_LAPACK_DIR) clean
	@rm -f $(NETLIB_LAPACK_DIR)/make.inc $(NETLIB_LAPACK_DIR)/lapacke/include/lapacke_mangling.h
	@rm -f *.grd Makefile.conf_last config_last.h
	@echo Done.
