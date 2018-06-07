# for release build and SM=35
# $ make release -j8 SM=35
# for debug build and default SM
# $ make debug -j8
# default build builds release and SM=30
# $ make -j8

# Builds THOR executable
CC = nvcc

SM:=30 # Streaming Multiprocessor version
arch := -arch sm_$(SM)

path_hd  := $(shell pwd)/src/headers
path_src := $(shell pwd)/src
obj_cuda   := esp.o grid.o esp_initial.o planet.o thor_driver.o profx_driver.o esp_output.o 
obj_cpp := storage.o binary_test.o config_file.o cmdargs.o
obj_tests := cmdargs_test.o cmdargs.o
obj := $(obj_cpp) $(obj_cuda)
headers := $(path_hd)/define.h $(path_hd)/grid.h $(path_hd)/planet.h $(path_hd)/esp.h \
           $(path_hd)/dyn/thor_fastmodes.h $(path_hd)/dyn/thor_adv_cor.h $(path_hd)/dyn/thor_auxiliary.h \
           $(path_hd)/dyn/thor_vertical_int.h $(path_hd)/dyn/thor_slowmodes.h $(path_hd)/dyn/thor_diff.h \
           $(path_hd)/dyn/thor_div.h $(path_hd)/phy/profx_auxiliary.h $(path_hd)/phy/profx_held_suarez.h \
           $(path_hd)/storage.h $(path_hd)/binary_test.h $(path_hd)/debug.h $(path_hd)/config_file.h \
           $(path_hd)/config_file.h

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

# define where to find sources
source_dirs := src src/grid src/initial src/thor src/profx src/output src/devel src/input src/test

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
includedir = ${path_hd}

# directory names 
OBJDIR = obj
BINDIR = bin
RESDIR = results
TESTDIR = tests

.PHONY: all clean

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

debug: all
release: all
prof: all

#######################################################################
# main binary
all: $(BINDIR)/$(OUTPUTDIR)/esp $(BINDIR)/esp

$(info Compile mode: $(MODE))
$(info Output objects to: $(OBJDIR)/$(OUTPUTDIR))
$(info Output tests to: $(BINDIR)/$(TESTDIR))
$(info flags: $(flags))
$(info arch: $(arch))


# create object directory if missing
$(BINDIR) $(RESDIR) $(OBJDIR):
	mkdir $@

$(BINDIR)/${OUTPUTDIR}: $(BINDIR)
	mkdir -p $(BINDIR)/$(OUTPUTDIR)

$(OBJDIR)/${OUTPUTDIR}: $(OBJDIR)
	mkdir -p $(OBJDIR)/$(OUTPUTDIR)

$(BINDIR)/${TESTDIR}: $(BINDIR)
	mkdir -p $(BINDIR)/$(TESTDIR)


# make dependency lists
# for CUDA files
$(OBJDIR)/${OUTPUTDIR}/%.d: %.cu | $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
	@set -e; rm -f $@; \
	$(CC) $(arch) $(dep_flags) $(h5include) $(h5libdir) -I$(includedir) --generate-dependencies $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

# for C++ files
$(OBJDIR)/${OUTPUTDIR}/%.d: %.cpp | $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
	@set -e; rm -f $@; \
	$(CC) $(arch) $(dep_flags) $(h5include) $(h5libdir) -I$(includedir)  --generate-dependencies  $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

# build objects
# CUDA files
$(OBJDIR)/${OUTPUTDIR}/%.o: %.cu $(OBJDIR)/$(OUTPUTDIR)/%.d| $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
	@echo creating $@
	$(CC) $(arch)  $(flags) $(h5include) $(h5libdir) -I$(includedir) -dc -o $@ $<
	@echo done $@

# C++ files
$(OBJDIR)/${OUTPUTDIR}/%.o: %.cpp $(OBJDIR)/$(OUTPUTDIR)/%.d| $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
	@echo creating $@
	$(CC) $(arch)  $(flags) $(h5include) $(h5libdir) -I$(includedir) -dc -o $@ $< 
	@echo done $@


# link *.o objects
$(BINDIR)/${OUTPUTDIR}/esp: $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj)) | $(BINDIR) $(RESDIR) $(BINDIR)/$(OUTPUTDIR)  $(OBJDIR)
	@echo creating $@
	$(CC) $(arch) $(link_flags) -o $(BINDIR)/$(OUTPUTDIR)/esp $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj))  $(h5libdir) $(h5libs)
	rm -f bin/esp
	ln -s bin/$(OUTPUTDIR)/esp -r -t bin


#######################################################################
# Build Tests

tests: ${BINDIR}/${TESTDIR}/cmdargs_test

$(BINDIR)/$(TESTDIR)/cmdargs_test:  $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests)) | $(BINDIR)/${OUTPUTDIR} $(BINDIR)/$(TESTDIR) $(BINDIR) $(RESDIR) 
	@echo creating $@
	$(CC) $(arch) $(flags) $(debug_flags) -o $(BINDIR)/$(TESTDIR)/cmdargs_test $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj_tests))  $(h5libdir) $(h5libs)

#######################################################################
# Cleanup 
.phony: clean,ar
clean:
	-rm -f bin/debug/esp $(addprefix $(OBJDIR)/debug/,$(obj)) $(obj:%.o=$(OBJDIR)/debug/%.d)
	-rm -f bin/release/esp $(addprefix $(OBJDIR)/release/,$(obj)) $(obj:%.o=$(OBJDIR)/release/%.d)
	-rm -f bin/prof/esp $(addprefix $(OBJDIR)/prof/,$(obj)) $(obj:%.o=$(OBJDIR)/prof/%.d)
	-rm -f bin/tests/cmdargs_test $(addprefix $(OBJDIR)/debug/,$(obj_tests)) $(obj_tests:%.o=$(OBJDIR)/debug/%.d)
	-rm -f bin/esp

#######################################################################
# dependencies includes
ifeq "${MODE}" "tests"
include $(obj_tests:%.o=$(OBJDIR)/$(OUTPUTDIR)/%.d)
else
include $(obj:%.o=$(OBJDIR)/$(OUTPUTDIR)/%.d)
endif
