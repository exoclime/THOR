# for release build and SM=35
# $ make release -j8 SM=35
# for debug build and default SM
# $ make debug -j8
# default build builds release and SM=30
# $ make -j8

# to show commands echoed
# $ make VERBOSE=1

# shell colors
RED := "\033[1;31m"
YELLOW := "\033[1;33m"
GREEN := "\033[1;32m"
BLUE := "\033[1;34m"
MAGENTA := "\033[1;35m"
CYAN := "\033[1;36m"
END := "\033[0m"

# Builds THOR executable
CC = nvcc

SM:=35 # Streaming Multiprocessor version
arch := -arch sm_$(SM)

# objects
obj_cuda   := esp.o grid.o esp_initial.o planet.o thor_driver.o profx_driver.o esp_output.o
obj_cpp := storage.o binary_test.o config_file.o cmdargs.o directories.o
obj := $(obj_cpp) $(obj_cuda)

#objects for tests
obj_tests_cmdargs := cmdargs_test.o cmdargs.o
obj_tests_config := config_test.o config_file.o
obj_tests_storage := storage_test.o storage.o
obj_tests_directories := directories_test.o directories.o



#######################################################################
# flags
# define specific compiler. if if fails on newer installations, get it to use g++-5
ccbin :=
# ccbin := -ccbin g++-5

# define common flags
flags := $(ccbin) --compiler-options -Wall -std=c++11
dep_flags := $(ccbin) -std=c++11
link_flags := $(ccbin)
# define debug flags
debug_flags := -g -G

# define release flags
release_flags := -O3

# define profiling flags
profiling_flags := -pg

#######################################################################
# define where to find sources
source_dirs := src src/grid src/initial src/thor src/profx src/output src/devel src/input src/files src/test

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
	flags += $(profiling_flags) -DBUILD_LEVEL="\"profiling\""
	MODE := prof
	OUTPUTDIR := prof
endif

# debug target
ifeq "$(findstring debug, $(MAKECMDGOALS))" "debug"
	flags += $(debug_flags) -DBUILD_LEVEL="\"debug\""
	MODE := debug
	OUTPUTDIR := debug
endif

# release target
ifeq "$(findstring release, $(MAKECMDGOALS))" "release"
	flags += $(release_flags) -DBUILD_LEVEL="\"release\""
	MODE := release
	OUTPUTDIR := release
endif

# test target
ifeq "$(findstring tests, $(MAKECMDGOALS))" "tests"
	flags += $(debug_flags) -DBUILD_LEVEL="\"test\""
	MODE := tests
	OUTPUTDIR := debug
endif

# by default, build release target
ifeq "$(MODE)" "UNDEF"
	flags += $(release_flags) -DBUILD_LEVEL="\"release\""
	MODE := release
	OUTPUTDIR := release
endif

debug: symlink
release: symlink
prof: symlink

#######################################################################
# main binary
all: symlink

$(info Compile mode: $(MODE))
$(info Output objects to: $(OBJDIR)/$(OUTPUTDIR))
$(info Output tests to: $(BINDIR)/$(TESTDIR))
$(info flags: $(flags) )
$(info arch: $(arch) )


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
	@set -e; rm -f $@; \
	$(CC) $(arch) $(dep_flags) $(h5include) $(h5libdir) -I$(includedir) --generate-dependencies $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

# for C++ files
$(OBJDIR)/${OUTPUTDIR}/%.d: %.cpp | $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
	@echo $(BLUE)computing dependencies $@ $(END)
	@set -e; rm -f $@; \
	$(CC) $(arch) $(dep_flags) $(h5include) $(h5libdir) -I$(includedir)  --generate-dependencies  $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

#######################################################################
# build objects
# CUDA files
$(OBJDIR)/${OUTPUTDIR}/%.o: %.cu $(OBJDIR)/$(OUTPUTDIR)/%.d| $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
	@echo $(YELLOW)creating $@ $(END)
	$(CC) $(arch)  $(flags) $(h5include) $(h5libdir) -I$(includedir) -dc -o $@ $<

# C++ files
$(OBJDIR)/${OUTPUTDIR}/%.o: %.cpp $(OBJDIR)/$(OUTPUTDIR)/%.d| $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
	@echo $(YELLOW)creating $@ $(END)
	$(CC) $(arch)  $(flags) $(h5include) $(h5libdir) -I$(includedir) -dc -o $@ $<


# link *.o objects
$(BINDIR)/${OUTPUTDIR}/esp: $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj)) | $(BINDIR) $(RESDIR) $(BINDIR)/$(OUTPUTDIR)  $(OBJDIR)
	@echo $(YELLOW)creating $@ $(END)
	$(CC) $(arch) $(link_flags) -o $(BINDIR)/$(OUTPUTDIR)/esp $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj))  $(h5libdir) $(h5libs)

# phony so that it will always be run
.PHONY: symlink
symlink: $(BINDIR)/$(OUTPUTDIR)/esp
	@echo $(BLUE)make link from $(BINDIR)/$(OUTPUTDIR)/esp to $(BINDIR)/esp  $(END)
	rm -f $(BINDIR)/esp
	ln -s $(BINDIR)/$(OUTPUTDIR)/esp -r -t bin


#######################################################################
# Build Tests

#define build_test =
#	$(CC) $(arch) $(flags) $(debug_flags) -o $(BINDIR)/$(TESTDIR)/cmdargs_test $(addprefix #$(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_cmdargs))
#endef

tests: ${BINDIR}/${TESTDIR}/cmdargs_test ${BINDIR}/${TESTDIR}/config_test ${BINDIR}/${TESTDIR}/storage_test ${BINDIR}/${TESTDIR}/directories_test

$(BINDIR)/$(TESTDIR)/cmdargs_test:  $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_cmdargs)) | $(BINDIR)/${OUTPUTDIR} $(BINDIR)/$(TESTDIR) $(BINDIR) $(RESDIR)
	@echo $(YELLOW)creating $@ $(END)
	$(CC) $(arch) $(flags) $(debug_flags) -o $(BINDIR)/$(TESTDIR)/cmdargs_test $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_cmdargs))

$(BINDIR)/$(TESTDIR)/config_test:  $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_config)) | $(BINDIR)/${OUTPUTDIR} $(BINDIR)/$(TESTDIR) $(BINDIR) $(RESDIR)
	@echo $(YELLOW)creating $@ $(END)
	$(CC) $(arch) $(flags) $(debug_flags) -o $(BINDIR)/$(TESTDIR)/config_test $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_config))

$(BINDIR)/$(TESTDIR)/directories_test:  $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_directories)) | $(BINDIR)/${OUTPUTDIR} $(BINDIR)/$(TESTDIR) $(BINDIR) $(RESDIR)
	@echo $(YELLOW)creating $@ $(END)
	$(CC) $(arch) $(flags) $(debug_flags) -o $(BINDIR)/$(TESTDIR)/directories_test $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_directories))

$(BINDIR)/$(TESTDIR)/storage_test:  $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_storage)) | $(BINDIR)/${OUTPUTDIR} $(BINDIR)/$(TESTDIR) $(BINDIR) $(RESDIR)
	@echo $(YELLOW)creating $@ $(END)
	$(CC) $(arch) $(flags) $(debug_flags) -o $(BINDIR)/$(TESTDIR)/storage_test $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests_storage))  $(h5libdir) $(h5libs)


#######################################################################
# Cleanup
.phony: clean,ar
clean:
	@echo $(CYAN)clean up binaries $(END)
	-$(RM) $(BINDIR)/debug/esp
	-$(RM) $(BINDIR)/release/esp
	-$(RM) $(BINDIR)/prof/esp
	@echo $(CYAN)clean up objects files and dependencies $(END)
	-$(RM) $(addprefix $(OBJDIR)/debug/,$(obj)) $(obj:%.o=$(OBJDIR)/debug/%.d)
	-$(RM) $(addprefix $(OBJDIR)/release/,$(obj)) $(obj:%.o=$(OBJDIR)/release/%.d)
	-$(RM) $(addprefix $(OBJDIR)/prof/,$(obj)) $(obj:%.o=$(OBJDIR)/prof/%.d)
	@echo $(CYAN)clean up tests binaries $(END)
	-$(RM) $(BINDIR)/tests/cmdargs_test
	-$(RM) $(BINDIR)/tests/storage_test
	-$(RM) $(BINDIR)/tests/config_test
	-$(RM) $(BINDIR)/tests/directories_test
	@echo $(CYAN)clean up test object files $(END)
	-$(RM) $(addprefix $(OBJDIR)/debug/,$(obj_tests_storage))
	-$(RM) $(addprefix $(OBJDIR)/debug/,$(obj_tests_config))
	-$(RM) $(addprefix $(OBJDIR)/debug/,$(obj_tests_cmdargs))
	-$(RM) $(addprefix $(OBJDIR)/debug/,$(obj_tests_directories))
	@echo $(CYAN)clean up test dependencies $(END)
	-$(RM) $(obj_tests_cmdargs:%.o=$(OBJDIR)/debug/%.d)
	-$(RM) $(obj_tests_storage:%.o=$(OBJDIR)/debug/%.d)
	-$(RM) $(obj_tests_config:%.o=$(OBJDIR)/debug/%.d)
	-$(RM) $(obj_tests_directories:%.o=$(OBJDIR)/debug/%.d)
	@echo $(CYAN)clean up symlink $(END)
	-$(RM) $(BINDIR)/esp
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
else
include $(obj:%.o=$(OBJDIR)/$(OUTPUTDIR)/%.d)
endif
endif
