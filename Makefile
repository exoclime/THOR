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
RED := \033[1;31m
YELLOW := \033[1;33m
GREEN := \033[1;32m
BLUE := \033[1;34m
MAGENTA := \033[1;35m
CYAN := \033[1;36m
END := \033[0m



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
obj_cuda   := esp.o grid.o esp_initial.o simulation_setup.o thor_driver.o profx_driver.o esp_output.o debug_helpers.o profx_globdiag.o reduction_add.o phy_modules_device.o ultrahot_thermo.o profx_sponge.o cuda_device_memory.o insolation.o

obj_cpp := storage.o binary_test.o config_file.o cmdargs.o directories.o log_writer.o iteration_timer.o
obj := $(obj_cpp) $(obj_cuda)

#objects for tests
obj_tests_cmdargs := cmdargs_test.o cmdargs.o
obj_tests_config := config_test.o config_file.o
obj_tests_storage := storage_test.o storage.o log_writer.o directories.o
obj_tests_directories := directories_test.o directories.o log_writer.o
obj_tests_gen_init := gen_init.o storage.o grid.o simulation_setup.o directories.o log_writer.o
obj_tests_reduction_add := reduction_add_test.o reduction_add.o  directories.o log_writer.o




#######################################################################
# flags


CUDA_PATH := /usr/lib/cuda/
CUDA_LIBS := /usr/lib/x86-64-linux-gnu/

ifeq ($(COMP), nvcc)
	# define specific compiler for nvcc. if if fails on newer installations, get it to use g++-5
	CC = nvcc
	CC_comp_flag = -dc
	ccbin :=
	# ccbin := -ccbin g++-5
	CDB :=
	HAS_CBD := none
	arch := -arch sm_$(SM)
	dependencies_flags = --compiler-options -MMD,-MP,"-MT $(OBJDIR)/$(OUTPUTDIR)/$(notdir $(basename $@)).o","-MF $(OBJDIR)/$(OUTPUTDIR)/$(DEPDIR)/$(notdir $(basename $@)).d"

	cuda_dependencies_flags = --generate-nonsystem-dependencies -MF $(OBJDIR)/$(OUTPUTDIR)/$(DEPDIR)/$(notdir $(basename $@)).d -MT "$(OBJDIR)/$(OUTPUTDIR)/$(notdir $(basename $@)).o $(OBJDIR)/$(OUTPUTDIR)/$(DEPDIR)/$(notdir $(basename $@)).d" 


	# define common flags
	cpp_flags := $(ccbin)  --compiler-options  -Wall -std=c++14 -DDEVICE_SM=$(SM)
	cuda_flags := $(ccbin) --compiler-options  -Wall -std=c++14 -DDEVICE_SM=$(SM)

	cpp_dep_flags := $(ccbin) -std=c++14
	cuda_dep_flags := $(ccbin) -std=c++14
	link_flags = $(ccbin)
else
	# Compiolatiopn with clang to create compilation database for external tools
	# need to compile with clang for compilation database
	CC := $(COMP)
	CC_comp_flag = -c
	CDB = -MJ $@.json
	HAS_CBD := -MJ
	arch := --cuda-gpu-arch=sm_$(SM)
	dependencies_flags := -MP -MMD -MF $(OBJDIR)/$(OUTPUTDIR)/$.d
	cuda_dependencies_flags := -MP -MMD -MF $(OBJDIR)/$(OUTPUTDIR)/$.d

	# define common flags
	cpp_flags := -Wall -std=c++14 -DDEVICE_SM=$(SM)
	cuda_flags := -Wall -std=c++14 -DDEVICE_SM=$(SM) --cuda-path=$(CUDA_PATH)

	cpp_dep_flags := -std=c++14
	cuda_dep_flags := -std=c++14 --cuda-path=$(CUDA_PATH)
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
source_dirs := src src/ESP src/utils src/test

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
DEPDIR = deps
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
# test for alfrodull
ifeq ($(wildcard Alfrodull),)
	ALFRODULL_FLAGS =
	ALFRODULL_LINK_FLAGS =
else 

	ALFRODULL_FLAGS = -DHAS_ALFRODULL -IAlfrodull/thor_module/inc/
    ALFRODULL_LINK_FLAGS =  -LAlfrodull -lalfrodull
	ALFRODULL_TARGET = Alfrodull/libalfrodull.a
export
# always call submakefile for alf
.PHONY: Alfrodull/libalfrodull.a
Alfrodull/libalfrodull.a:
	@echo -e '$(MAGENTA)Creating Alfrodull from subdir Alfrodull$(END)'
	$(MAKE) -f Makefile -C Alfrodull

endif 




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
$(info CC compile flag: $(CC_comp_flag))

ifndef VERBOSE
.SILENT:
endif

GITREV_FILE = $(OBJDIR)/$(OUTPUTDIR)/git-rev.h


.PHONY: $(GITREV_FILE)
$(GITREV_FILE): | $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
	@echo -e '$(GREEN)Creating version file $(END)'
	@echo -n "#define GIT_BRANCH_RAW " > $(GITREV_FILE)
	@echo \"$(shell git rev-parse --abbrev-ref HEAD)\" >> $(GITREV_FILE)
	@echo -n "#define GIT_HASH_RAW " >> $(GITREV_FILE)
	@echo \"$(shell git describe --always --dirty --abbrev=40 --match="NoTagWithThisName")\" >> $(GITREV_FILE)


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

$(OBJDIR)/${OUTPUTDIR}/$(DEPDIR): $(OBJDIR)/${OUTPUTDIR} $(OBJDIR)
	mkdir -p $(OBJDIR)/$(OUTPUTDIR)/$(DEPDIR)

#######################################################################
# build objects
# CUDA files
$(OBJDIR)/${OUTPUTDIR}/esp.o: esp.cu $(OBJDIR)/${OUTPUTDIR}/${DEPDIR}/esp.d $(GITREV_FILE) | $(OBJDIR)/${OUTPUTDIR}/$(DEPDIR) $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
	@echo -e '$(BLUE)creating dependencies for $@ $(END)'
	$(CC) $(cuda_dependencies_flags) $(CC_comp_flag) $(arch)  $(cuda_flags) $(h5include) -I$(includedir)  -I$(OBJDIR)/$(OUTPUTDIR) $(CDB) $(ALFRODULL_FLAGS) -o $@ $<
	@echo -e '$(YELLOW)creating object file for $@ $(END)'
	$(CC) $(CC_comp_flag) $(arch)  $(cuda_flags) $(h5include) -I$(includedir) -I$(OBJDIR)/$(OUTPUTDIR) $(CDB) $(ALFRODULL_FLAGS) -o $@ $<


$(OBJDIR)/${OUTPUTDIR}/%.o: %.cu $(OBJDIR)/${OUTPUTDIR}/${DEPDIR}/%.d  | $(OBJDIR)/${OUTPUTDIR}/$(DEPDIR) $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR) 
	@echo -e '$(BLUE)creating dependencies for $@ $(END)'
	$(CC) $(cuda_dependencies_flags) $(CC_comp_flag) $(arch)  $(cuda_flags) $(h5include) -I$(includedir) $(CDB) $(ALFRODULL_FLAGS) -o $@ $<
	@echo -e '$(YELLOW)creating object file for $@ $(END)'
	$(CC) $(CC_comp_flag) $(arch)  $(cuda_flags) $(h5include) -I$(includedir) $(CDB) $(ALFRODULL_FLAGS) -o $@ $<

# C++ files
$(OBJDIR)/${OUTPUTDIR}/%.o: %.cpp $(OBJDIR)/${OUTPUTDIR}/${DEPDIR}/%.d  | $(OBJDIR)/${OUTPUTDIR}/$(DEPDIR) $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
	@echo -e '$(YELLOW)creating dependencies and object file for $@  $(END)'
	$(CC) $(dependencies_flags) $(CC_comp_flag) $(arch) $(cpp_flags) $(h5include) -I$(includedir) $(CDB) $(ALFRODULL_FLAGS) -o $@ $<

# Target for dependencies, so that they are always defined
$(OBJDIR)/${OUTPUTDIR}/${DEPDIR}/%.d: ;

# link *.o objects
.PHONY: 
$(BINDIR)/${OUTPUTDIR}/esp: $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj)) $(MODULES_SRC)/libphy_modules.a $(ALFRODULL_TARGET) | $(BINDIR) $(RESDIR) $(BINDIR)/$(OUTPUTDIR)  $(OBJDIR)
	@echo -e '$(YELLOW)creating $@ $(END)'
	$(CC) -o $(BINDIR)/$(OUTPUTDIR)/esp $(arch) $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj)) -L$(MODULES_SRC) -lphy_modules $(h5libdir) $(h5libs) $(link_flags) $(ALFRODULL_LINK_FLAGS)
	if [ "$(HAS_CDB)" = "-MJ" ] ; then \
		rm -f compile_commands.json; \
		sed -e '1s/^/[\n/' -e '$$s/,$$/\n]/' $(OBJDIR)/$(OUTPUTDIR)/*.o.json $(MODULES_SRC)/obj/*.o.json > compile_commands.json; \
		sed -i.bak s/-xcuda/-xc++/g compile_commands.json; \
		rm -f compile_commands.json.bak; \
	fi


# phony so that it will always be run
.PHONY: symlink
symlink: $(BINDIR)/$(OUTPUTDIR)/esp   | $(OBJDIR) $(GITREV_FILE)
	@echo -e '$(BLUE)make link from $(BINDIR)/$(OUTPUTDIR)/esp to $(BINDIR)/esp  $(END)'
	rm -f $(BINDIR)/esp
	ln -s $(BINDIR)/$(OUTPUTDIR)/esp -r -t bin

# dependencies
DEPFILES := $($(obj):%.o=$(OBJDIR)/$(OUTPUTDIR)/$(DEPDIR)/%.d)
$(DEPFILES):

#######################################################################
# include physics modules

export
# always call submakefile for modules
.PHONY: $(MODULES_SRC)/libphy_modules.a
$(MODULES_SRC)/libphy_modules.a:
	@echo -e '$(MAGENTA)Creating physics module from subdir $(MODULES_SRC)$(END)'
	$(MAKE) -f Makefile -C $(MODULES_SRC)

#######################################################################
# Build Tests

#define build_test =
#	$(CC) $(arch) $(flags) $(debug_flags) -o $(BINDIR)/$(TESTDIR)/cmdargs_test $(addprefix #$(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_cmdargs))
#endef

tests: ${BINDIR}/${TESTDIR}/cmdargs_test ${BINDIR}/${TESTDIR}/config_test ${BINDIR}/${TESTDIR}/storage_test ${BINDIR}/${TESTDIR}/directories_test ${BINDIR}/${TESTDIR}/gen_init  ${BINDIR}/${TESTDIR}/reduction_add_test

$(BINDIR)/$(TESTDIR)/cmdargs_test:  $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_cmdargs)) | $(BINDIR)/${OUTPUTDIR} $(BINDIR)/$(TESTDIR) $(BINDIR) $(RESDIR)
	@echo -e '$(YELLOW)creating $@ $(END)'
	$(CC) $(arch) $(link_flags) -o $(BINDIR)/$(TESTDIR)/cmdargs_test $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_cmdargs))

$(BINDIR)/$(TESTDIR)/config_test:  $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_config)) | $(BINDIR)/${OUTPUTDIR} $(BINDIR)/$(TESTDIR) $(BINDIR) $(RESDIR)
	@echo -e '$(YELLOW)creating $@ $(END)'
	$(CC)  $(arch) $(link_flags) -o $(BINDIR)/$(TESTDIR)/config_test $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_config))

$(BINDIR)/$(TESTDIR)/directories_test:  $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_directories)) | $(BINDIR)/${OUTPUTDIR} $(BINDIR)/$(TESTDIR) $(BINDIR) $(RESDIR)
	@echo -e '$(YELLOW)creating $@ $(END)'
	$(CC) $(arch) $(link_flags) -o $(BINDIR)/$(TESTDIR)/directories_test $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_directories))

$(BINDIR)/$(TESTDIR)/storage_test:  $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_storage)) | $(BINDIR)/${OUTPUTDIR} $(BINDIR)/$(TESTDIR) $(BINDIR) $(RESDIR)
	@echo -e '$(YELLOW)creating $@ $(END)'
	$(CC) $(arch) $(link_flags) -o $(BINDIR)/$(TESTDIR)/storage_test $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_storage))  $(h5libdir) $(h5libs)

$(BINDIR)/$(TESTDIR)/gen_init:  $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_gen_init)) | $(BINDIR)/${OUTPUTDIR} $(BINDIR)/$(TESTDIR) $(BINDIR) $(RESDIR)
	@echo -e '$(YELLOW)creating $@ $(END)'
	$(CC) $(arch) $(link_flags) -o $(BINDIR)/$(TESTDIR)/gen_init $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_gen_init))  $(h5libdir) $(h5libs)

$(BINDIR)/$(TESTDIR)/reduction_add_test:  $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_reduction_add)) | $(BINDIR)/${OUTPUTDIR} $(BINDIR)/$(TESTDIR) $(BINDIR) $(RESDIR)
	@echo -e '$(YELLOW)creating $@ $(END)'
	$(CC) $(arch) $(link_flags) -o $(BINDIR)/$(TESTDIR)/reduction_add_test $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_reduction_add))  $(h5libdir) $(h5libs)

#######################################################################
# Cleanup
.phony: clean,ar
clean:
	@echo -e '$(CYAN)clean up binaries $(END)'
	-$(RM) $(BINDIR)/debug/esp
	-$(RM) $(BINDIR)/release/esp
	-$(RM) $(BINDIR)/prof/esp
	echo -e '$(CYAN)clean up objects files and dependencies $(END)'
	-$(RM) $(OBJDIR)/debug/*.o $(OBJDIR)/debug/$(DEPDIR)/*.d $(OBJDIR)/debug/$(DEPDIR)/*.d.* $(OBJDIR)/debug/*.o.json
	-$(RM) $(OBJDIR)/release/*.o $(OBJDIR)/release/$(DEPDIR)/*.d $(OBJDIR)/release/$(DEPDIR)/*.d.* $(OBJDIR)/release/*.o.json
	-$(RM) $(OBJDIR)/prof/*.o $(OBJDIR)/prof/$(DEPDIR)/*.d $(OBJDIR)/prof/$(DEPDIR)/*.d.* $(OBJDIR)/prof/*.o.json
	-$(RM) $(GITREV_FILE)
	@echo -e '$(CYAN)clean up tests binaries $(END)'
	-$(RM) $(BINDIR)/tests/cmdargs_test
	-$(RM) $(BINDIR)/tests/storage_test
	-$(RM) $(BINDIR)/tests/config_test
	-$(RM) $(BINDIR)/tests/directories_test
	-$(RM) $(BINDIR)/tests/gen_init
	-$(RM) $(BINDIR)/tests/test_reduction_add
	@echo -e '$(CYAN)clean up symlink $(END)'
	-$(RM) $(BINDIR)/esp
	@echo -e '$(CYAN)clean up modules directory $(END)'
	$(MAKE) -C $(MODULES_SRC) clean
	@echo -e '$(CYAN)clean up directories $(END)'
	-$(RM) -d $(BINDIR)/debug $(BINDIR)/release $(BINDIR)/prof
	-$(RM) -d $(OBJDIR)/debug/$(DEPDIR)
	-$(RM) -d $(OBJDIR)/release/$(DEPDIR)
	-$(RM) -d $(OBJDIR)/prof/$(DEPDIR)
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
include $(obj_tests_config:%.o=$(OBJDIR)/$(OUTPUTDIR)/$(DEPDIR)/%.d)
include $(obj_tests_cmdargs:%.o=$(OBJDIR)/$(OUTPUTDIR)/$(DEPDIR)/%.d)
include $(obj_tests_storage:%.o=$(OBJDIR)/$(OUTPUTDIR)/$(DEPDIR)/%.d)
include $(obj_tests_directories:%.o=$(OBJDIR)/$(OUTPUTDIR)/$(DEPDIR)/%.d)
include $(obj_tests_gen_init:%.o=$(OBJDIR)/$(OUTPUTDIR)/$(DEPDIR)/%.d)
include $(obj_tests_reduce_add:%.o=$(OBJDIR)/$(OUTPUTDIR)/$(DEPDIR)/%.d)
else
include $(obj:%.o=$(OBJDIR)/$(OUTPUTDIR)/$(DEPDIR)/%.d)
endif
endif
