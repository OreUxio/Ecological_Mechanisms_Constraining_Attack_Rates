# $Id: Makefile 2476 2016-10-23 11:37:49Z axel $
# you may want to change this part as you work:

default: NewWorld Inspect Miniweb CompetiWeb IntervalTest TAGS

PROFILE = no # yes or no
DEBUG = no  # yes or no
DEP_ON_MAKEFILE = yes # yes or no, default depends on machine
PARALLEL = yes # yes or no (this option is defunct)
STATIC_BINARIES = no

#FEEDBACK= -fprofile-generate
#FEEDBACK= -fprofile-use
#FEEDBACK= --coverage
#FEEDBACK=-fprofile-arcs -ftest-coverage

#SX8 = yes

#FEEDBACK= -pg

# hardware/system specific setup ##################################
ifdef SX
OS = SUPER-UX
CPU = SX
else
OS = $(shell uname -s)
CPU = $(shell uname -m)
endif
HOSTNAME = $(shell hostname)

ARCH = $(OS).$(CPU)

CFLAGS+=-D_XOPEN_SOURCE=600

ifeq ($(OS), Linux)
ifneq ($(findstring XXdell,$(HOSTNAME)),)
	CXXFLAGS += -DUSE_OPENMP
	CXX_COMPILER := icc -g -openmp -O2 $(FEEDBACK)
	CC_COMPILE := icc -g -openmp -O2 $(FEEDBACK)
else
# ifneq ($(findstring login00,$(HOSTNAME)),)
# 	CXXFLAGS += -DUSE_OPENMP
# 	CXX_COMPILER := icc -g -openmp -O2 $(FEEDBACK)
# 	CC_COMPILE := icc -g -openmp -O2 $(FEEDBACK)
# else
	CXX_COMPILER := g++ -g -O2  $(FEEDBACK)
	CC_COMPILE := gcc -g -O2  $(FEEDBACK)
# endif
endif
	GSL = yes
	PREFIX := $(HOME)
	ifdef SYSTEM_FOUND
	ERROR := $(error "confusion identifying system or hardware")
	else
	SYSTEM_FOUND = true
	endif
	ifndef DEP_ON_MAKEFILE
	DEP_ON_MAKEFILE = yes
	endif
	LIBS += gfortran pthread
ifneq ($(findstring login00,$(HOSTNAME)),) 
	LIB_DIRS += /usr/lib64
endif
ifeq ($(findstring dell,$(HOSTNAME)),dell)
	INC_DIRS += /home/arossberg/include
	LIB_DIRS += /home/arossberg/lib
	INC_DIRS += $(wildcard /cm/shared/apps/[^a]*/include)
	LIB_DIRS += $(wildcard /cm/shared/apps/[^a]*/lib)
	LIB_DIRS += /cm/shared/apps/gotoblas/core/64/1.26
endif
	CXX_OPTIM_FLAGS += -O2 # -Wdisabled-optimization
ifeq ($(findstring login00,$(HOSTNAME)),) 
	CXX_OPTIM_FLAGS += -msse -msse2
	CXX_OPTIM_FLAGS += -maccumulate-outgoing-args 
	CXX_OPTIM_FLAGS += -fmodulo-sched-allow-regmoves
	CXX_OPTIM_FLAGS += -fno-signed-zeros
	CXX_OPTIM_FLAGS += -fmodulo-sched -fmodulo-sched-allow-regmoves
	CXX_OPTIM_FLAGS += -fgcse-las -fgcse-after-reload
	CXX_OPTIM_FLAGS += -fstrict-aliasing
# for Zoom profiler:
	CXXFLAGS += -fno-omit-frame-pointer 
	CXX_OPTIM_FLAGS += -ftree-vectorize #-ftree-vectorizer-verbose=5
endif
##	CXX_OPTIM_FLAGS += -fipa-matrix-reorg
#	CXX_OPTIM_FLAGS += -ftree-loop-im
#	CXX_OPTIM_FLAGS += -ftree-loop-ivcanon
# 	this make code MUCH SLOWER, but may help with larger webs:
#	CXX_OPTIM_FLAGS += -funroll-loops -fvariable-expansion-in-unroller -fprefetch-loop-arrays -fpredictive-commoning
#	CXX_OPTIM_FLAGS += -freorder-blocks-and-partition
	CXXFLAGS += -Wno-unused-result
ifeq ($(HOSTNAME),almond)
	CXX_OPTIM_FLAGS += -msse4.1 -msse4.2
else
	CXX_OPTIM_FLAGS +=  -momit-leaf-frame-pointer
endif
endif
ifeq ($(OS), Darwin)
	CXX_COMPILER := /sw/bin/g++-4 -g -Wno-deprecated -Xlinker -bind_at_load $(FEEDBACK)
	CC_COMPILE := /sw/bin/gcc-4 -g $(FEEDBACK)
	GSL = yes
	PREFIX := /sw/
	ifdef SYSTEM_FOUND
	ERROR := $(error "confusion identifying system or hardware")
	else
	SYSTEM_FOUND = true
	endif
	ifndef DEP_ON_MAKEFILE
	DEP_ON_MAKEFILE = yes
	endif
	LDFLAGS += -lgfortran
	LIB_DIRS += /sw/lib
	INC_DIRS += /sw/include
	INC_DIRS += /sw/lib/gcc4.3/include
	CXXFLAGS += -D NO_CONST_WEBNAME_FILTER
	CXX_OPTIM_FLAGS += -O2 -Wdisabled-optimization
	CXX_OPTIM_FLAGS += -ffast-math -mfpmath=sse -march=native 
	CXX_OPTIM_FLAGS += -maccumulate-outgoing-args -momit-leaf-frame-pointer
	CXX_OPTIM_FLAGS += -fmodulo-sched -fmodulo-sched-allow-regmoves
	CXX_OPTIM_FLAGS += -fgcse-las -fgcse-after-reload
	CXX_OPTIM_FLAGS += -fsee
##	CXX_OPTIM_FLAGS += -fipa-matrix-reorg
#	CXX_OPTIM_FLAGS += -ftree-loop-im
#	CXX_OPTIM_FLAGS += -ftree-loop-ivcanon
# 	this make code MUCH SLOWER, but may help with larger webs:
#	CXX_OPTIM_FLAGS += -funroll-loops -fvariable-expansion-in-unroller -fprefetch-loop-arrays -fpredictive-commoning
#	CXX_OPTIM_FLAGS += -freorder-blocks-and-partition
	CXX_OPTIM_FLAGS += -fstrict-aliasing
	CXX_OPTIM_FLAGS += -fno-signed-zeros
endif
ifeq ($(OS), SUPER-UX)
	CXX_COMPILER := . /SX/etc/$(SX).sh;sxc++ -$(SX) -K exceptions -T local -O fullmsg -pi fullmsg -pvctl fullmsg
	CC_COMPILE :=  . /SX/etc/$(SX).sh;sxcc -$(SX)
	CXXFLAGS +=  -R diaglist -R objlist
	CCFLAGS +=  -R diaglist -R objlist
	GSL = no
	LDFLAGS += -f90lib #for lapack
	PREFIX := $(HOME)/axel/$(SX)
	ifdef SYSTEM_FOUND
	ERROR := $(error "confusion identifying system or hardware")
	else
	SYSTEM_FOUND = true
	endif
	ifndef DEP_ON_MAKEFILE
	DEP_ON_MAKEFILE = no
	endif
endif
ifeq ($(OS), IRIX64)
	CXX_COMPILER := CC -no_prelink
	CC_COMPILE := cc -no_prelink
	GSL = no
	PREFIX := $(HOME)/axel
	CXXFLAGS += -woff 3333 # no warning on LEDA memory management
	CXXFLAGS += -woff 1460 # redeclared "inline" after being called
	CXXFLAGS += -woff 1355 # The extra semicolon (";") is ignored.
	CXXFLAGS += -woff 1195 # floating-point operation out of range.
	CCFLAGS += -woff 1185 # enumerated type mixed with other type
	ifdef SYSTEM_FOUND
	ERROR := $(error "confusion identifying system or hardware")
	else
	SYSTEM_FOUND = true
	endif
	ifndef DEP_ON_MAKEFILE
	DEP_ON_MAKEFILE = no
	endif
endif

ifeq ($(OS),Darwin)
BOOST_LIB_SUFFIX = -mt
else 
ifneq ($(findstring dell,$(HOSTNAME)),)
BOOST_LIB_SUFFIX = 
else
ifneq ($(findstring Linux,$(OS)),)
BOOST_LIB_SUFFIX = 
else
BOOST_LIB_SUFFIX = -gcc43-mt
endif
endif
endif

# ... insert other specific setups here
ifdef SYSTEM_FOUND
else
	ERROR := $(error "system or hardware could not be identified")
endif

ifeq ($(findstring yes,$(DEP_ON_MAKEFILE)),yes)
	MAKEFILE_DEP=Makefile
endif

################# end of hardware/system specific setup ###########
###################################################################
# executables: ####################################################

sources =  $(wildcard *.cc *.c)

# we use some heuristics to guess the executable names:
executable_sources := $(shell grep -l '^int main' $(sources))
executable_c_sources := $(filter %.c, $(executable_sources))
executable_cc_sources := $(filter %.cc, $(executable_sources))
executables = $(shell echo $(executable_sources) | sed '/s/\.\(cc|c\)$//')

################# end of executables ##############################
###################################################################
# libs and includes: ##############################################

## sundials is now a system library, no special paths required (and no
## profiling possible!)
# ifeq ($(findstring yes,$(PROFILE)),yes)
# SUNDIALS_DIR = $(PREFIX)/sundials-profile
# else
# SUNDIALS_DIR = $(PREFIX)/sundials
# endif

#LIB_DIRS += /usr/lib/atlas /usr/local/atlas/lib
#INC_DIRS += /usr/local/atlas/include

## sundials is now a system library, no special paths required (and no
## profiling possible!)
#LIB_DIRS += $(SUNDIALS_DIR)/lib
#INC_DIRS += $(SUNDIALS_DIR)/include

LIB_DIRS += $(PREFIX)/lib
INC_DIRS += $(PREFIX)/include
INC_DIRS += /home/af06/tmp/numeric_bindings

ifeq ($(OS), IRIX64) 
LIB_DIRS += /usr/freeware/lib
INC_DIRS += /usr/freeware/include
endif

#deepgeen module include and lib folders 
#guessing the $HOSTNAME of deepgreen contains "deep". 
ifneq ($(findstring deep,$(HOSTNAME)),)
# You may try to put this sequence into ~/.profile instead (worked for me):
## # Add ~/lib to LD_LIBRARY_PATH
## export LD_LIBRARY_PATH="${HOME}/lib:${LD_LIBRARY_PATH}"
## # For each /lib in LD_LIBRARY_PATH put the corresponding /include in CPATH
## # CPATH is the compiler's general include path
## export CPATH=$(echo "${LD_LIBRARY_PATH}:"|sed 's+/lib:+/include:+g'|sed 's+:$++' )
## export LIBRARY_PATH="${LIBRARY_PATH}:${LD_LIBRARY_PATH}"
# --- end of code for ~/.profile
LIB_DIRS += /apps/boost/1.54.0/lib 
INC_DIRS += /apps/boost/1.54.0/include
LIB_DIRS += /apps/clhep/2.1.3.1/lib 
INC_DIRS += /apps/clhep/2.1.3.1/include 
LIB_DIRS += /apps/sundials/2.5.0/lib 
INC_DIRS += /apps/sundials/2.5.0/include 
LIB_DIRS += /apps/nlopt/2.3/lib 
INC_DIRS += /apps/nlopt/2.3/include 
LIB_DIRS += /apps/xerces/3.1.1/lib 
INC_DIRS += /apps/xerces/3.1.1/include 
LIB_DIRS += /apps/fftw/1.54.0/lib 
INC_DIRS += /apps/fftw/3.3.3/include 
LIB_DIRS += /scratch/apps/boost/1.54.0/lib 
INC_DIRS += /scratch/apps/boost/1.54.0/include 
LIB_DIRS += /scratch/apps/clhep/2.1.3.1/lib 
INC_DIRS += /scratch/apps/clhep/2.1.3.1/include 
LIB_DIRS += /scratch/apps/sundials/2.5.0/lib 
INC_DIRS += /scratch/apps/sundials/2.5.0/include 
LIB_DIRS += /scratch/apps/nlopt/2.3/lib 
INC_DIRS += /scratch/apps/nlopt/2.3/include 
LIB_DIRS += /scratch/apps/xerces/3.1.1/lib 
INC_DIRS += /scratch/apps/xerces/3.1.1/include 
LIB_DIRS += /scratch/apps/fftw/1.54.0/lib 
INC_DIRS += /scratch/apps/fftw/3.3.3/include 
LIB_DIRS += /usr/lib64/atlas
endif

# INC_DIRS += /usr/include/xercesc-2.7.0
# add also library directory /usr/lib64/xerces-c-2.7.? when it exists:
# LIB_DIRS += $(subst /usr/lib64/xerces-c-2.7.?,,$(wildcard /usr/lib64/xerces-c-2.7.?))

packed_simulation_asm.o : packed_simulation_asm.s $(MAKEFILE_DEP)
	gcc -c -masm=intel packed_simulation_asm.s

# we are not so picky to choose just the right libs:
MERPSpectrum SizeSpectrum : LIBS += fftw3_threads fftw3 pthread plotter Xaw Xmu Xt SM ICE Xext X11
Miniweb : LIBS += boost_iostreams$(BOOST_LIB_SUFFIX)
LIBS += sundials_cvode sundials_kinsol #rt sundials_shared
                        # sundial ODE solver
LIBS += xerces-c	# xerces-c XML parser 

LIBS += armadillo

ifeq (1,0)
ifeq ($(findstring g++,$(CXX_COMPILER)),g++)
ifeq ($(findstring dell,$(HOSTNAME)),dell)
CLHEP_SUFFIX = -2.0.4.6
endif
ifeq ($(findstring almond,$(HOSTNAME)),almond)
CLHEP_SUFFIX = -2.1.0.1
endif
ifneq ($(findstring adrian-VirtualBox,$(HOSTNAME)),)
CLHEP_SUFFIX = -2.1.2.3
endif
# ifeq ($(findstring pine,$(HOSTNAME)),pine)
# LIBS += CLHEP-Matrix$(CLHEP_SUFFIX) # *some* matrices (should better use boost or gsl)
# LIBS += CLHEP-Vector$(CLHEP_SUFFIX) # 
# LIBS += CLHEP-Random$(CLHEP_SUFFIX) # for random (should better use boost or gsl)
# else
LIBS += CLHEP
# endif
else
ifeq ($(findstring sxc++,$(CXX_COMPILER)),sxc++)
ifeq ($(SX),sx5)
LIBS += CLHEP-sx5c++	# for random, evaluator, and some matrices
else
LIBS += CLHEP-sxc++	# for random, evaluator, and some matrices
endif
else
ifeq ($(findstring dell,$(HOSTNAME)),dell)
CLHEP_SUFFIX = -2.0.4.6
# LIBS += CLHEP-Matrix$(CLHEP_SUFFIX) # *some* matrices 
LIBS += CLHEP-Vector$(CLHEP_SUFFIX) # 
LIBS += CLHEP-Random$(CLHEP_SUFFIX) # for random
else
LIBS += CLHEP-cc	# for random, evaluator, and some matrices
endif
endif
endif
endif #eq (1,0)

ifeq ($(findstring yes,$(GSL)),yes)
LIBS += gsl # for some special function (gslcblas may be obsolete)
endif

ifeq ($(OS), IRIX64)
LIBS += scs 
else
ifeq ($(findstring dell,$(HOSTNAME)),dell)
LIBS += lapack cblas goto atlas
else
LIBS += openblas lapack gslcblas blas
#LIBS += lapack gslcblas blas
endif
#ifeq ($(findstring yes,$(STATIC_BINARIES)),yes)
#LIBS += g2c
#endif
endif

# ODE.o : HOMPACK90.o hompack90_calls.o
# LDFLAGS += HOMPACK90.o hompack90_calls.o

NewWorld Inspect: LIBS += nlopt

ifeq ($(findstring login00,$(HOSTNAME)),) 
LIBS += bz2
else
NewWorld Inspect : LDFLAGS += /usr/lib64/libbz2.so.1 
endif

LIBS += boost_serialization$(BOOST_LIB_SUFFIX) boost_thread$(BOOST_LIB_SUFFIX) boost_system$(BOOST_LIB_SUFFIX)

#LIBS += m 		# floating-point exception handling might require this

################# end of libs and includes: #######################
###################################################################
# now, so some calculations: ######################################


ifeq ($(findstring CC,$(CXX_COMPILER)),CC)

ifeq ($(findstring yes,$(DEBUG)),yes)
	DEBUG_OPTS += -g2 -DDEBUGGING -fno-inline 
	LDFLAGS += -W l,-v
else
	CXXFLAGS += -O3 -apo -IPA -apolist -mp -mplist
	CCFLAGS += -O3 -apo -IPA
	LDFLAGS += -O3 -apo -IPA -apolist -mp -mplist
endif

endif # CC compiler

ifeq ($(findstring sxc++,$(CXX_COMPILER)),sxc++)

ifeq ($(findstring yes,$(PROFILE)),yes)
#	DEBUG_OPTS += -p
	DEBUG_OPTS += -ftrace
endif

ifeq ($(findstring yes,$(DEBUG)),yes)
	DEBUG_OPTS += -DDEBUGGING #-C debug # debug compiles to slow
else
	CXXFLAGS += -C aopt -pi,copy_arg  
#	CXXFLAGS += -Orestrict=all #dangerous!
#	CXXFLAGS += -alias,pointer,type,variable #dangerous!
#	CXXFLAGS += -pvctl,nodep #dangerous!
#	CXXFLAGS += -pi,nocopy_arg #dangerous!
	CXXFLAGS += -pi,func_size=500000,times=1000,max_depth=1000,max_size=500000 #very agressive inlining
#	CXXFLAGS += -O extendreorder #increases compilation time
endif

endif # sxc++ compiler

ifeq ($(findstring g++,$(CXX_COMPILER)),g++)

ifeq ($(findstring yes,$(PROFILE)),yes)
	DEBUG_OPTS += -pg -g -O2 -fprofile-generate
endif

ifeq ($(findstring yes,$(DEBUG)),yes)
	DEBUG_OPTS += -O0 -g -DDEBUGGING -fno-inline -finline-limit=0
else
	CXXFLAGS += $(CXX_OPTIM_FLAGS)
endif

endif # g++ compiler

CC := $(CC_COMPILE)
CXX := $(CXX_COMPILER)

LDLIBS  += $(patsubst %, -l%, $(LIBS))

LDFLAGS +=  $(patsubst %, -L%, $(LIB_DIRS))

CPPFLAGS += $(patsubst %, -I%, $(INC_DIRS))
ifeq ($(OS),Linux)
 CPPFLAGS+= -Wno-deprecated #for hash_map
endif

ifeq ($(findstring yes,$(PARALLEL)),yes)
	CPPFLAGS += -D PARALLEL
endif

ifeq ($(findstring CC,$(CXX_COMPILER)),g++)
  CPPFLAGS+= -D HAS_TYPEOF
endif
ifeq ($(OS),IRIX64)
  CCFLAGS+= -LANG:std 
  CPPFLAGS+= -LANG:std 
  LDFLAGS+= -LANG:std 
  CPPFLAGS += -D ON_SX5FSV
endif
ifeq ($(findstring yes,$(GSL)),yes)
  CPPFLAGS += -D WITH_GSL
endif

#CPPFLAGS += -D __BUILTIN_BOOL__ # avoid redefining type bool in LEDA
#CPPFLAGS += -D AVOID_CONFLICTS_WITH_NEW # also for LEDA

# these should be last to overwrite others:
CXXFLAGS += $(DEBUG_OPTS)
CCFLAGS += $(DEBUG_OPTS)
LDFLAGS += $(DEBUG_OPTS)

SHELL = bash

################# end of calculations: ############################
###################################################################
# now, do the actual work: ########################################

.PHONY: run

ifdef SX
ifeq ($(SX),sx8)
ifeq ($(findstring yes,$(PARALLEL)),yes)
QUEUENAME:=Ap
else
QUEUENAME:=As
endif
else
ifeq ($(SX),sx5)
ifeq ($(findstring yes,$(PARALLEL)),yes)
ifeq ($(findstring yes,$(PROFILE)),yes)
QUEUENAME:=Sp
else
QUEUENAME:=Mp
endif
else
ifeq ($(findstring yes,$(PROFILE)),yes)
QUEUENAME:=Ss
else
QUEUENAME:=Ms
endif
endif
endif
endif

bat.installed :  $(wildcard *.bat)
	rsync -u $(wildcard *.bat) /hwork0/d6/ishii/
	touch bat.installed

installed : NewWorld NewWorld.cfg
	mkdir -p /hwork0/d6/ishii/$(QUEUENAME)
	rsync -u NewWorld NewWorld.cfg /hwork0/d6/ishii/$(QUEUENAME)
ifeq ($(findstring yes,$(PROFILE)),yes)
#	-rm /hwork0/d6/ishii/$(QUEUENAME)/web*.xml
endif
#	cd /hwork0/d6/ishii/$(QUEUENAME);. /SX/etc/$(SX).sh;sxstrip NewWorld
	touch installed

run : installed bat.installed
	@echo runnning in queue $(QUEUENAME)
ifeq ($(QUEUENAME),Sp)
	cd /hwork0/d6/ishii;\
	qsub -q $(QUEUENAME) /hwork0/d6/ishii/bftrace.bat
#	qsub -q As /hwork0/d6/ishii/bprof	
else
	cd /hwork0/d6/ishii;\
	qsub -q $(QUEUENAME) /hwork0/d6/ishii/bc.bat
endif
	sleep 2
	qstatr | grep $(QUEUENAME)
endif



TAGS: $(wildcard *.cc *.c *.h) $(MAKEFILE_DEP)
ifeq ($(OS), Linux)
# etags -m ??
	@etags --language="c++" `ls -t *.cc *.h` \
	      --language="c" `ls -t *.c` ||\
	etags --language="c++" `ls -t *.cc *.h` \
	      --language="c" `ls -t *.c`
else
	@etags --c++ `ls -t *.cc *.h`  `ls -t *.c`
endif

# extra compiler options
CXXFLAGS += $(X_CXXFLAGS)

ifdef MAKEFILE_DEP
%.o : %.c Makefile
	$(CC_COMPILE) -c $(CPPFLAGS) $(CCFLAGS) $< -o $@ 
%.o : %.cc Makefile
	$(CXX_COMPILER) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@ 
%.o : %.f Makefile
	gfortran -g -c $< -o $@ 
else
%.o : %.c
	$(CC_COMPILE) -c $(CPPFLAGS) $(CCFLAGS) $< -o $@ 
%.o : %.cc
	$(CXX_COMPILER) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@ 
%.o : %.f
	gfortran -g -c $< -o $@ 
endif

ifeq ($(OS), IRIX64)
CC = $(CXX_COMPILER)
else 
ifeq ($(OS), SUPER-UX)
CC = $(CXX_COMPILER)
else
# for linking, this is sufficient:
$(patsubst %.cc, %, $(executable_cc_sources)) : CC = $(CXX_COMPILER)
#but use CC_COMPILE for normal c-files
$(patsubst %.c, %.o, $(wildcard *.c)) : CC = $(CC_COMPILER)
LINK.cc = $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(TARGET_ARCH)
STATIC_LINK.cc = $(CXX) -static $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(TARGET_ARCH)

#  commands to execute (built-in):
%: %.o nvector_blas.o
	$(CC) compile_time.cc -c
	$(LINK.cc) $^ compile_time.o $(LOADLIBES) $(LDLIBS) $(STATICS) -o $@
ifeq ($(findstring dell,$(HOSTNAME)),dell)
ifeq ($(findstring rossberg,$(LD_LIBRARY_PATH)),)
	$(warning Before executing, do 'export LD_LIBRARY_PATH="$$LD_LIBRARY_PATH:/home/arossberg/lib:/cm/shared/apps/clhep/lib:/cm/shared/apps/gsl/lib:/cm/shared/apps/sundials/lib:/cm/shared/apps/xerces/lib"' or so)
endif
endif
ifeq ($(findstring yes,$(STATIC_BINARIES)),yes)
	$(warning $(LDLIBS))
	$(warning $(subst sundials_ ,sundials_f , $(LDLIBS) ))
	$(STATIC_LINK.cc) $^ compile_time.o $(LOADLIBES) $(STATICS) $(subst sundials_,sundials_f, $(LDLIBS) ) -o $@.static
endif
endif
endif

NewWorld: NewWorld.o nvector_blas.o
	$(CC) compile_time.cc -c
	$(LINK.cc) $^ compile_time.o $(LOADLIBES) $(LDLIBS) $(STATICS) -o $@

Inspect: Inspect.o nvector_blas.o
	$(CC) compile_time.cc -c
	$(LINK.cc) $^ compile_time.o $(LOADLIBES) $(LDLIBS) $(STATICS) -o $@




ifndef X_CXXFLAGS
ifeq ($(findstring sxc++,$(CXX_COMPILER)),sxc++)
#we need some more complier option to get perfect vectorization
hotspots.o : hotspots.cc 
	@echo doing hotspots.o with extra options
	$(MAKE) X_CXXFLAGS="-pvctl,nodep -O extendreorder -dir inline -pi fullmsg" hotspots.o
NewWeb.o : NewWeb.cc 
	@echo doing NewWeb.o with extra options
	$(MAKE) X_CXXFLAGS="-pvctl,nodep -dir inline -pi fullmsg" NewWeb.o 
NetworkAnalysis.o : NetworkAnalysis.cc 
	@echo doing NetworkAnalysis.o with debug options to avoid compiler bug
	$(MAKE) DEBUG=yes X_CXXFLAGS="-C debug" NetworkAnalysis.o
endif
endif

# we need to put cfgList.o in first position in the list of executables
# so we do it by hand for now:
ifeq ($(OS), IRIX64)
NewWorld : NewWorld.o ODE.o error.o scheduler.o parsecfg.o evaluate.o XMLStore.o remember.o matrix_transformers.o NewWeb.o NetworkAnalysis.o CMatrix.o NewSpecies.o sequence.o NetworkHelpers.o indices.o random_pick_field.o standardize.o cfgPermanent.o random.o Statistics.o linpack_eigen.o
	$(CXX_COMPILER) $(LDFLAGS) cfgList.o random.o NewWorld.o ODE.o error.o scheduler.o parsecfg.o evaluate.o XMLStore.o remember.o matrix_transformers.o NewWeb.o NetworkAnalysis.o CMatrix.o NewSpecies.o sequence.o NetworkHelpers.o indices.o random_pick_field.o standardize.o cfgPermanent.o Statistics.o linpack_eigen.o mp_sug_numthreads.o  $(LDLIBS) -o NewWorld
endif

################# end of actual work: ############################
##################################################################
# automatic dependencies: ########################################

include $(shell echo $(sources) ' '|sed 's,\.cc ,.d ,g'|sed 's,\.c ,.d ,g')



%.o: %.d

%.d1: %.cc $(MAKEFILE_DEP)
# Automatically make and load dependencies
# Automatically collect required object files for executables
	@$(CXX) -M $(CPPFLAGS) -DGET_INDIRECT_DEPENDENCIES $< |\
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' |\
	sed 's,\([: ]\)/usr/[^ ]*,\1,g'|\
	grep -v '^ *\\$$' \
	> $@

.PRECIOUS : %.d1

%.d: %.d1 $(MAKEFILE_DEP)
	@cp $*.d1 $*.d
	@echo $* : $(filter $(wildcard *.h),$(shell cat $*.d1))' ' |\
	sed 's,\.h ,.o ,g' \
	>> $*.d


%.d: %.c $(MAKEFILE_DEP)
# Automatically make and load dependencies
# Automatically collect required object files for executables
	@$(CC) -M $(CPPFLAGS) -DGET_INDIRECT_DEPENDENCIES $< |\
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' |\
	sed 's,\([: ]\)/usr/[^ ]*,\1,g'|\
	grep -v '^ *\\$$' \
	> $@

################# end of automatic dependencies: #################
##################################################################
# others: ########################################################

test : test.cc Makefile
	echo $$SX $$SX_BASE_PSUITE
	$(CXX_COMPILER) test.cc -o test

.PHONY: doxygen

doxygen: 
	doxygen

doxygen-publish: doxygen
	ncftpput -zf ~/.rossberg.net foodwebs/NewWeb-doc/ `find ../doc/html -newer ../doc/.last_upload`
	touch ../doc/.last_upload

### Local Variables: 
### End: 
