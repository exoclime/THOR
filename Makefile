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
# ccbin := 
# ccbin := -ccbin g++-5
flag := $(ccbin) -g -G -std=c++11

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

all: $(BINDIR)/esp

# create object directory if missing
$(OBJDIR) $(BINDIR) $(RESDIR):
	mkdir $@

# make dependency lists
$(OBJDIR)/%.d: %.cu |$(OBJDIR)
	@set -e; rm -f $@; \
	$(CC) $(arch) $(flags) $(h5include) $(h5libdir) -I$(includedir) --generate-dependencies $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

$(OBJDIR)/%.d: %.cpp |$(OBJDIR)
	@set -e; rm -f $@; \
	$(CC) $(arch) $(flag) $(h5include) $(h5libdir) -I$(includedir)  --generate-dependencies  $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

# build *.cu CUDA objects
$(OBJDIR)/%.o: %.cu $(OBJDIR)/%.d|$(OBJDIR)
	@echo creating $@
	$(CC) $(arch)  $(flag) $(h5include) $(h5libdir) -I$(includedir) -dc -o $@ $<
	@echo done $@

$(OBJDIR)/%.o: %.cpp $(OBJDIR)/%.d|$(OBJDIR)
	@echo creating $@
	$(CC) $(arch)  $(flag) $(h5include) $(h5libdir) -I$(includedir) -dc -o $@ $< 
	@echo done $@

# link *.o objects
$(BINDIR)/esp: $(addprefix $(OBJDIR)/,$(obj)) | $(BINDIR) $(RESDIR)
	@echo creating $@
	$(CC) $(arch) $(flag) -o $(BINDIR)/esp $(addprefix $(OBJDIR)/,$(obj))  $(h5libdir) $(h5libs)


.phony: clean,ar
clean:
	-rm bin/esp $(addprefix $(OBJDIR)/,$(obj)) $(obj:%.o=$(OBJDIR)/%.d)

include $(obj:%.o=$(OBJDIR)/%.d)

