# Builds THOR executable
CC = nvcc

sm:=61 # Streaming Multiprocessor version
arch := -arch sm_$(sm)

path_hd  := $(shell pwd)/src/headers
path_src := $(shell pwd)/src
obj_cuda   := esp.o grid.o esp_initial.o planet.o thor_driver.o profx_driver.o esp_output.o 
obj_cpp := storage.o binary_test.o
obj := $(obj_cpp) $(obj_cuda)
headers := $(path_hd)/define.h $(path_hd)/grid.h $(path_hd)/planet.h $(path_hd)/esp.h \
           $(path_hd)/dyn/thor_fastmodes.h $(path_hd)/dyn/thor_adv_cor.h $(path_hd)/dyn/thor_auxiliary.h \
           $(path_hd)/dyn/thor_vertical_int.h $(path_hd)/dyn/thor_slowmodes.h $(path_hd)/dyn/thor_diff.h \
           $(path_hd)/dyn/thor_div.h $(path_hd)/phy/profx_auxiliary.h $(path_hd)/phy/profx_held_suarez.h \
           $(path_hd)/storage.h $(path_hd)/binary_test.h $(path_hd)/debug.h

# define specific compiler. if if fails on newer installations, get it to use g++-5
ccbin := 
# ccbin := -ccbin g++-5

# define common flags
flags := $(ccbin) --compiler-options -Wall -std=c++11
dep_flags := $(ccbin) -std=c++11
# define debug flags
debug_flags := -g -G

# define release flags
release_flags := -O3

# define profiling flags
profiling_flags := -pg

# define where to find sources
source_dirs := src src/grid src/initial src/thor src/profx src/output src/devel

vpath %.cu $(source_dirs)
vpath %.cpp $(source_dirs)

h5libs := $(shell h5c++ -show -shlib | awk -v ORS=" " '{ for ( n=1; n<=NF; n++ ) if ($$n ~ "^-lh") print $$n  }')
h5libdir := $(shell h5c++ -show -shlib | awk -v ORS=" " '{ for ( n=1; n<=NF; n++ ) if ($$n ~ "^-L") print $$n  }')
h5include := $(shell h5c++ -show -shlib | awk -v ORS=" " '{ for ( n=1; n<=NF; n++ ) if ($$n ~ "^-I") print $$n  }')
$(info h5libs="$(h5libs)")
$(info h5libdir="$(h5libdir)")
$(info h5include="$(h5include)")

includehdf = $(h5include)
includedir = ${path_hd}
# Path where the hdf5 lib was installed,
# should contain `include` and `lib`
# hdf_dir := <hdf5 directory>
#
# h5lib = -L$(hdf_dir)/lib -lhdf5
# includehdf = -I$(hdf_dir)/include

OBJDIR = obj
BINDIR = bin
RESDIR = results
.PHONY: all clean

$(info goals: $(MAKECMDGOALS))

# default target
OUTPUTDIR := default
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
debug: all
release: all
prof: all


all: $(BINDIR)/$(OUTPUTDIR)/esp

$(info Compile mode: $(MODE))
$(info Output objects to: $(OBJDIR)/$(OUTPUTDIR))
$(info flags: $(flags))

# create object directory if missing
$(BINDIR) $(RESDIR) $(OBJDIR):
	mkdir $@

$(BINDIR)/${OUTPUTDIR}: $(BINDIR)
	mkdir -p $(BINDIR)/$(OUTPUTDIR)

$(OBJDIR)/${OUTPUTDIR}: $(OBJDIR)
	mkdir -p $(OBJDIR)/$(OUTPUTDIR)

# make dependency lists
$(OBJDIR)/${OUTPUTDIR}/%.d: %.cu | $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
	@set -e; rm -f $@; \
	$(CC) $(arch) $(dep_flags) $(h5include) $(h5libdir) -I$(includedir) --generate-dependencies $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

$(OBJDIR)/${OUTPUTDIR}/%.d: %.cpp | $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
	@set -e; rm -f $@; \
	$(CC) $(arch) $(dep_flags) $(h5include) $(h5libdir) -I$(includedir)  --generate-dependencies  $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

# build *.cu CUDA objects
$(OBJDIR)/${OUTPUTDIR}/%.o: %.cu $(OBJDIR)/$(OUTPUTDIR)/%.d| $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
	@echo creating $@
	$(CC) $(arch)  $(flags) $(h5include) $(h5libdir) -I$(includedir) -dc -o $@ $<
	@echo done $@

$(OBJDIR)/${OUTPUTDIR}/%.o: %.cpp $(OBJDIR)/$(OUTPUTDIR)/%.d| $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
	@echo creating $@
	$(CC) $(arch)  $(flags) $(h5include) $(h5libdir) -I$(includedir) -dc -o $@ $< 
	@echo done $@

# link *.o objects
$(BINDIR)/${OUTPUTDIR}/esp: $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj)) | $(BINDIR) $(RESDIR) $(BINDIR)/$(OUTPUTDIR)  $(OBJDIR)
	@echo creating $@
	$(CC) $(arch) $(flags) -o $(BINDIR)/$(OUTPUTDIR)/esp $(addprefix $(OBJDIR)/$(OUTPUTDIR)/,$(obj))  $(h5libdir) $(h5libs)
	rm -f bin/esp
	ln -s bin/$(OUTPUTDIR)/esp -r -t bin

.phony: clean,ar
clean:
	-rm -f bin/default/esp $(addprefix $(OBJDIR)/default/,$(obj)) $(obj:%.o=$(OBJDIR)/default/%.d)
	-rm -f bin/debug/esp $(addprefix $(OBJDIR)/debug/,$(obj)) $(obj:%.o=$(OBJDIR)/debug/%.d)
	-rm -f bin/release/esp $(addprefix $(OBJDIR)/release/,$(obj)) $(obj:%.o=$(OBJDIR)/release/%.d)
	-rm -f bin/prof/esp $(addprefix $(OBJDIR)/prof/,$(obj)) $(obj:%.o=$(OBJDIR)/prof/%.d)
	-rm -f bin/esp

include $(obj:%.o=$(OBJDIR)/$(OUTPUTDIR)/%.d)

