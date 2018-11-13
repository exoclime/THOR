# for release build and SM=35
# $ make release -j8 SM=35
# for debug build and default SM
# $ make debug -j8
# default build builds release and SM=30
# $ make -j8

# to test localy while forcing g++-5:
#  make release -j8 ccbin='-ccbin g++-5'
#
# to show commands echoed
# $ make VERBOSE=1

# SM number can also be set in Makefile.conf

# shell colors
RED := "\033[1;31m"
YELLOW := "\033[1;33m"
GREEN := "\033[1;32m"
BLUE := "\033[1;34m"
MAGENTA := "\033[1;35m"
CYAN := "\033[1;36m"
END := "\033[0m"



#######################################################################
# get local config
ifeq ($(wildcard Makefile.conf),Makefile.conf)
$(info Using local configuration from Makefile.conf )
include Makefile.conf
else
$(info No local configuration )
endif


#######################################################################
# Set up variables for files and compilation



# Builds THOR executable
SM ?= 35 # Streaming Multiprocessor version

COMP ?= nvcc




# objects
obj_cuda   := esp.o grid.o esp_initial.o planet.o thor_driver.o profx_driver.o esp_output.o debug_helpers.o profx_conservation.o reduction_add.o phy_modules_device.o

obj_cpp := storage.o binary_test.o config_file.o cmdargs.o directories.o log_writer.o iteration_timer.o
obj := $(obj_cpp) $(obj_cuda)

#objects for tests
obj_tests_cmdargs := cmdargs_test.o cmdargs.o
obj_tests_config := config_test.o config_file.o
obj_tests_storage := storage_test.o storage.o
obj_tests_directories := directories_test.o directories.o
obj_tests_gen_init := gen_init.o storage.o grid.o planet.o
obj_tests_reduction_add := reduction_add_test.o reduction_add.o




#######################################################################
# flags


CUDA_PATH := /usr/lib/cuda/
CUDA_LIBS := /usr/lib/x86-64-linux-gnu/

ifeq ($(COMP), nvcc)
	# define specific compiler for nvcc. if if fails on newer installations, get it to use g++-5
	CC = nvcc
	ccbin :=
	# ccbin := -ccbin g++-5
	CDB = none
	arch := -arch sm_$(SM)
	dependencies_flags = --generate-dependencies

	# define common flags
	cpp_flags := $(ccbin)  --compiler-options  -Wall -std=c++11 -DDEVICE_SM=$(SM)
	cuda_flags := $(ccbin) --compiler-options  -Wall -std=c++11 -DDEVICE_SM=$(SM)

	cpp_dep_flags := $(ccbin) -std=c++11
	cuda_dep_flags := $(ccbin) -std=c++11
	link_flags = $(ccbin)
else
	# need to compile with clang for compilation database
	CC := $(COMP)
	CDB = -MJ
	arch := --cuda-gpu-arch=sm_$(SM)
	dependencies_flags := -MM

	# define common flags
	cpp_flags := -Wall -std=c++11 -DDEVICE_SM=$(SM)
	cuda_flags := -Wall -std=c++11 -DDEVICE_SM=$(SM) --cuda-path=$(CUDA_PATH)

	cpp_dep_flags := -std=c++11
	cuda_dep_flags := -std=c++11 --cuda-path=$(CUDA_PATH)
	link_flags = --cuda-path=$(CUDA_PATH) -L$(CUDA_LIBS) -lcudart_static -ldl -lrt -pthread
endif





# define debug flags
debug_flags := -g -G

# define release flags
release_flags := -O3

# define profiling flags
profiling_flags := -pg -lineinfo --default-stream per-thread

#######################################################################
# define where to find sources
source_dirs := src src/grid src/initial src/thor src/profx src/output src/devel src/input src/files src/utils src/test

vpath %.cu $(source_dirs)
vpath %.cpp $(source_dirs)

# Path where the hdf5 lib was installed,
h5libs := $(shell h5c++ -show -shlib | awk -v ORS=" " '{ for ( n=1; n<=NF; n++ ) if ($$n ~ "^-lh") print $$n  }')
h5libdir := $(shell h5c++ -show -shlib | awk -v ORS=" " '{ for ( n=1; n<=NF; n++ ) if ($$n ~ "^-L") print $$n  }')
h5include := $(shell h5c++ -show -shlib | awk -v ORS=" " '{ for ( n=1; n<=NF; n++ ) if ($$n ~ "^-I") print $$n  }')
$(info h5libs="$(h5libs)")
$(info h5libdir="$(h5libdir)")
$(info h5include="$(h5include)")

includehdf = $(h5include)
includedir = src/headers

#######################################################################
# directory names
OBJDIR = obj
BINDIR = bin
RESDIR = results
TESTDIR = tests

# modules search directory.
# set if not set from command line or config file
# will run the Makefile in that directory, passing it this makefile's
# variables
MODULES_SRC ?= src/physics/managers/empty/

.PHONY: all clean

#######################################################################
# compute targets

$(info goals: $(MAKECMDGOALS))

# set some build values depending on target
# default target value, falls back to release below
OUTPUTDIR := UNDEF
MODE := UNDEF

# profiling target
ifeq "$(findstring prof, $(MAKECMDGOALS))" "prof"
	cuda_flags += $(profiling_flags) -DBUILD_LEVEL="\"profiling\""
	cpp_flags += $(profiling_flags) -DBUILD_LEVEL="\"profiling\""
	MODE := prof
	OUTPUTDIR := prof
endif

# debug target
ifeq "$(findstring debug, $(MAKECMDGOALS))" "debug"
	cuda_flags += $(debug_flags) -DBUILD_LEVEL="\"debug\""
	cpp_flags += $(debug_flags) -DBUILD_LEVEL="\"debug\""
	MODE := debug
	OUTPUTDIR := debug
endif

# release target
ifeq "$(findstring release, $(MAKECMDGOALS))" "release"
	cuda_flags += $(release_flags) -DBUILD_LEVEL="\"release\""
	cpp_flags += $(release_flags) -DBUILD_LEVEL="\"release\""
	MODE := release
	OUTPUTDIR := release
endif

# test target
ifeq "$(findstring tests, $(MAKECMDGOALS))" "tests"
	cuda_flags += $(debug_flags) -DBUILD_LEVEL="\"test\""
	cpp_flags += $(debug_flags) -DBUILD_LEVEL="\"test\""
	MODE := tests
	OUTPUTDIR := debug
endif

# by default, build release target
ifeq "$(MODE)" "UNDEF"
	cuda_flags += $(release_flags) -DBUILD_LEVEL="\"release\""
	cpp_flags += $(release_flags) -DBUILD_LEVEL="\"release\""
	MODE := release
	OUTPUTDIR := release
endif

debug: symlink
release: symlink
prof: symlink
cdb:


#######################################################################
# main binary
all: symlink

$(info Compiler: $(CC))
$(info Compile mode: $(MODE))
$(info Output objects to: $(OBJDIR)/$(OUTPUTDIR))
$(info Output tests to: $(BINDIR)/$(TESTDIR))
$(info cuda_flags: $(cuda_flags) )
$(info cpp_flags: $(cpp_flags) )
$(info link_flags: $(link_flags) )
$(info arch: $(arch) )
$(info MODULES_SRC: $(MODULES_SRC))

ifndef VERBOSE
.SILENT:
endif


#######################################################################
# create object directory if missing
$(BINDIR) $(RESDIR) $(OBJDIR):
	mkdir $@

$(BINDIR)/${OUTPUTDIR}: $(BINDIR)
	mkdir -p $(BINDIR)/$(OUTPUTDIR)

$(OBJDIR)/${OUTPUTDIR}: $(OBJDIR)
	mkdir -p $(OBJDIR)/$(OUTPUTDIR)

$(BINDIR)/${TESTDIR}: $(BINDIR)
	mkdir -p $(BINDIR)/$(TESTDIR)

#######################################################################
# make dependency targets
# this generates obj/(debug|release)/*.d files containing dependencies targets.
# it uses the compiler to find all the header files a CUDA/C++ file depends on
# those files are then included at the end
# make is run twice, once to create the dependency targets and once to run the compilation

# for CUDA files
$(OBJDIR)/${OUTPUTDIR}/%.d: %.cu | $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
	@echo $(BLUE)computing dependencies $@ $(END)
	set -e; rm -f $@; \
	$(CC) $(dependencies_flags) $(arch) $(cuda_dep_flags) $(h5include) -I$(includedir)  $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$


# for C++ files
$(OBJDIR)/${OUTPUTDIR}/%.d: %.cpp | $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
	@echo $(BLUE)computing dependencies $@ $(END)
	set -e; rm -f $@; \
	$(CC) $(dependencies_flags) $(arch) $(cpp_dep_flags) $(h5include) -I$(includedir) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

#######################################################################
# build objects
# CUDA files
$(OBJDIR)/${OUTPUTDIR}/%.o: %.cu $(OBJDIR)/$(OUTPUTDIR)/%.d| $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
	@echo $(YELLOW)creating $@ $(END)
	if test $$CDB = "-MJ" ; then \
		$(CC) -c $(arch)  $(cuda_flags) $(h5include) -I$(includedir) $(CDB) $@.json -o $@ $<; \
	else \
		$(CC) -c $(arch)  $(cuda_flags) $(h5include) -I$(includedir) -o $@ $<; \
	fi

# C++ files
$(OBJDIR)/${OUTPUTDIR}/%.o: %.cpp $(OBJDIR)/$(OUTPUTDIR)/%.d| $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
	@echo $(YELLOW)creating $@ $(END)
	if test $$CDB = "-MJ" ; then \
		$(CC) -c $(arch) $(cpp_flags) $(h5include) -I$(includedir) $(CDB) $@.json -o $@ $<; \
	else \
		$(CC) -c $(arch) $(cpp_flags) $(h5include) -I$(includedir) -o $@ $<; \
	fi


# link *.o objects
$(BINDIR)/${OUTPUTDIR}/esp: $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj)) $(MODULES_SRC)/libphy_modules.a | $(BINDIR) $(RESDIR) $(BINDIR)/$(OUTPUTDIR)  $(OBJDIR)
	@echo $(YELLOW)creating $@ $(END)
	$(CC) -o $(BINDIR)/$(OUTPUTDIR)/esp $(arch) $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj)) -L$(MODULES_SRC) -lphy_modules $(h5libdir) $(h5libs) $(link_flags)
	if test $$CDB = "-MJ" ; then \
	rm -f compile_commands.json; \
	sed -e '1s/^/[\n/' -e '$$s/,$$/\n]/' $(OBJDIR)/$(OUTPUTDIR)/*.o.json $(MODULES_SRC)/obj/*.o.json > compile_commands.json; \
	sed -i.bak s/-xcuda/-xc++/g compile_commands.json; \
	rm -f compile_commands.json.bak; \
	fi

# phony so that it will always be run
.PHONY: symlink
symlink: $(BINDIR)/$(OUTPUTDIR)/esp
	@echo $(BLUE)make link from $(BINDIR)/$(OUTPUTDIR)/esp to $(BINDIR)/esp  $(END)
	rm -f $(BINDIR)/esp
	ln -s $(BINDIR)/$(OUTPUTDIR)/esp -r -t bin

#######################################################################
# include physics modules

export
# always call submakefile for modules
.PHONY: $(MODULES_SRC)/libphy_modules.a
$(MODULES_SRC)/libphy_modules.a:
	@echo $(MAGENTA)Creating physics module from subdir $(MODULES_SRC)$(END)
	$(MAKE) -f Makefile -C $(MODULES_SRC)

#######################################################################
# Build Tests

#define build_test =
#	$(CC) $(arch) $(flags) $(debug_flags) -o $(BINDIR)/$(TESTDIR)/cmdargs_test $(addprefix #$(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_cmdargs))
#endef

tests: ${BINDIR}/${TESTDIR}/cmdargs_test ${BINDIR}/${TESTDIR}/config_test ${BINDIR}/${TESTDIR}/storage_test ${BINDIR}/${TESTDIR}/directories_test ${BINDIR}/${TESTDIR}/gen_init  ${BINDIR}/${TESTDIR}/reduction_add_test

$(BINDIR)/$(TESTDIR)/cmdargs_test:  $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_cmdargs)) | $(BINDIR)/${OUTPUTDIR} $(BINDIR)/$(TESTDIR) $(BINDIR) $(RESDIR)
	@echo $(YELLOW)creating $@ $(END)
	$(CC) $(CDB) $(arch) $(link_flags) -o $(BINDIR)/$(TESTDIR)/cmdargs_test $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_cmdargs))

$(BINDIR)/$(TESTDIR)/config_test:  $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_config)) | $(BINDIR)/${OUTPUTDIR} $(BINDIR)/$(TESTDIR) $(BINDIR) $(RESDIR)
	@echo $(YELLOW)creating $@ $(END)
	$(CC) $(CDB) $(arch) $(link_flags) -o $(BINDIR)/$(TESTDIR)/config_test $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_config))

$(BINDIR)/$(TESTDIR)/directories_test:  $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_directories)) | $(BINDIR)/${OUTPUTDIR} $(BINDIR)/$(TESTDIR) $(BINDIR) $(RESDIR)
	@echo $(YELLOW)creating $@ $(END)
	$(CC) $(CDB) $(arch) $(link_flags) -o $(BINDIR)/$(TESTDIR)/directories_test $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_directories))

$(BINDIR)/$(TESTDIR)/storage_test:  $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_storage)) | $(BINDIR)/${OUTPUTDIR} $(BINDIR)/$(TESTDIR) $(BINDIR) $(RESDIR)
	@echo $(YELLOW)creating $@ $(END)
	$(CC) $(CDB) $(arch) $(link_flags) -o $(BINDIR)/$(TESTDIR)/storage_test $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_storage))  $(h5libdir) $(h5libs)

$(BINDIR)/$(TESTDIR)/gen_init:  $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_gen_init)) | $(BINDIR)/${OUTPUTDIR} $(BINDIR)/$(TESTDIR) $(BINDIR) $(RESDIR)
	@echo $(YELLOW)creating $@ $(END)
	$(CC) $(CDB) $(arch) $(link_flags) -o $(BINDIR)/$(TESTDIR)/gen_init $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_gen_init))  $(h5libdir) $(h5libs)

$(BINDIR)/$(TESTDIR)/reduction_add_test:  $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_reduction_add)) | $(BINDIR)/${OUTPUTDIR} $(BINDIR)/$(TESTDIR) $(BINDIR) $(RESDIR)
	@echo $(YELLOW)creating $@ $(END)
	$(CC) $(CDB) $(arch) $(link_flags) -o $(BINDIR)/$(TESTDIR)/reduction_add_test $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_reduction_add))  $(h5libdir) $(h5libs)

#######################################################################
# Cleanup
.phony: clean,ar
clean:
	@echo $(CYAN)clean up binaries $(END)
	-$(RM) $(BINDIR)/debug/esp
	-$(RM) $(BINDIR)/release/esp
	-$(RM) $(BINDIR)/prof/esp
	@echo $(CYAN)clean up objects files and dependencies $(END)
	-$(RM) $(addprefix $(OBJDIR)/debug/,$(obj)) $(obj:%.o=$(OBJDIR)/debug/%.d) $(obj:%.o=$(OBJDIR)/debug/%.o.json)
	-$(RM) $(addprefix $(OBJDIR)/release/,$(obj)) $(obj:%.o=$(OBJDIR)/release/%.d) $(obj:%.o=$(OBJDIR)/release/%.o.json)
	-$(RM) $(addprefix $(OBJDIR)/prof/,$(obj)) $(obj:%.o=$(OBJDIR)/prof/%.d) $(obj:%.o=$(OBJDIR)/prof/%.o.json)
	@echo $(CYAN)clean up tests binaries $(END)
	-$(RM) $(BINDIR)/tests/cmdargs_test
	-$(RM) $(BINDIR)/tests/storage_test
	-$(RM) $(BINDIR)/tests/config_test
	-$(RM) $(BINDIR)/tests/directories_test
	-$(RM) $(BINDIR)/tests/gen_init
	-$(RM) $(BINDIR)/tests/test_reduction_add
	@echo $(CYAN)clean up test object files $(END)
	-$(RM) $(addprefix $(OBJDIR)/debug/,$(obj_tests_storage))
	-$(RM) $(addprefix $(OBJDIR)/debug/,$(obj_tests_config))
	-$(RM) $(addprefix $(OBJDIR)/debug/,$(obj_tests_cmdargs))
	-$(RM) $(addprefix $(OBJDIR)/debug/,$(obj_tests_directories))
	-$(RM) $(addprefix $(OBJDIR)/debug/,$(obj_tests_gen_init))
	@echo $(CYAN)clean up test dependencies $(END)
	-$(RM) $(obj_tests_cmdargs:%.o=$(OBJDIR)/debug/%.d)
	-$(RM) $(obj_tests_storage:%.o=$(OBJDIR)/debug/%.d)
	-$(RM) $(obj_tests_config:%.o=$(OBJDIR)/debug/%.d)
	-$(RM) $(obj_tests_directories:%.o=$(OBJDIR)/debug/%.d)
	-$(RM) $(obj_tests_gen_init:%.o=$(OBJDIR)/debug/%.d)
	-$(RM) $(obj_tests_reduction_add:%.o=$(OBJDIR)/debug/%.d)
	@echo $(CYAN)clean up symlink $(END)
	-$(RM) $(BINDIR)/esp
	@echo $(CYAN)clean up modules directory $(END)
	$(MAKE) -C $(MODULES_SRC) clean
	@echo $(CYAN)clean up directories $(END)
	-$(RM) -d $(BINDIR)/debug $(BINDIR)/release $(BINDIR)/prof
	-$(RM) -d $(BINDIR)/tests
	-$(RM) -d $(OBJDIR)/debug
	-$(RM) -d $(OBJDIR)/release
	-$(RM) -d $(OBJDIR)/prof
	-$(RM) -d $(BINDIR)
	-$(RM) -d $(OBJDIR)


#######################################################################
# dependencies includes
ifneq ($(MAKECMDGOALS),clean)
ifeq "${MODE}" "tests"
include $(obj_tests_config:%.o=$(OBJDIR)/$(OUTPUTDIR)/%.d)
include $(obj_tests_cmdargs:%.o=$(OBJDIR)/$(OUTPUTDIR)/%.d)
include $(obj_tests_storage:%.o=$(OBJDIR)/$(OUTPUTDIR)/%.d)
include $(obj_tests_directories:%.o=$(OBJDIR)/$(OUTPUTDIR)/%.d)
include $(obj_tests_gen_init:%.o=$(OBJDIR)/$(OUTPUTDIR)/%.d)
include $(obj_tests_reduce_add:%.o=$(OBJDIR)/$(OUTPUTDIR)/%.d)
else
include $(obj:%.o=$(OBJDIR)/$(OUTPUTDIR)/%.d)
endif
endif
