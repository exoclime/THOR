# makefile for RT physics modules, called from main makefile
# must create libphy_modules.a in its root directory

$(info -------------------------------------------------------------- )
$(info Modules Alfrodull Makefile )

# set some variables if not set
includedir ?= unset
h5include ?= unset
cpp_flags ?= unset
cuda_flags ?= unset
arch ?= unset
CC_comp_flag ?= unset
MODE ?= unset

$(info Some variables inherited from parent makefiles)
$(info includes: $(includedir))
$(info h5includes: $(h5include))
$(info cpp_flags: $(cpp_flags))
$(info cuda_flags: $(cuda_flags))
$(info CC compile flag: $(CC_comp_flag))
$(info dependeny flags: $(dependency_flags))
$(info arch: $(arch))
$(info MODE: $(MODE))

######################################################################
# Directories
THOR_ROOT = ../

# Includes
LOCAL_INCLUDE = src/inc 

# shared modules
SHARED_MODULES_INCLUDE = $(THOR_ROOT)src/physics/modules/inc/

# thor root include if we want to use code from there
THOR_INCLUDE = $(THOR_ROOT)src/headers

# source dirs
LOCAL_SOURCES = src
SHARED_MODULES_DIR = $(THOR_ROOT)src/physics/modules/src/
# object directory
BUILDDIR = obj
OUTPUTDIR = $(MODE)

.PHONY: all clean


######################################################################
$(info Sub Makefile variables)
$(info THOR root from submakefile: $(THOR_ROOT))

######################################################################
# all: libphy_modules.a libalfrodull.a
all: libalfrodull.a

# path to local module code
vpath %.cu src src/kernels src/opacities 

vpath %.cpp $(LOCAL_SOURCES) 
vpath %.h $(LOCAL_INCLUDE) 
# path to thor headers
vpath %.h $(THOR_INCLUDE)
# path to phy_modules
vpath %.h $(SHARED_MODULES_INCLUDE)
vpath %.cu $(SHARED_MODULES_DIR)
vpath %.cpp $(SHARED_MODULES_DIR)



# objects
obj_cuda := alfrodull_engine.o alfrodullib.o planck_table.o opacities.o cloud_opacities.o calculate_physics.o integrate_flux.o interpolate_values.o math_helpers.o two_streams_radiative_transfer.o 
obj_cpp := gauss_legendre_weights.o 
obj :=  $(obj_cuda) $(obj_cpp)


ifndef VERBOSE
.SILENT:
endif

dependencies_flags = --generate-dependencies
cuda_dependencies_flags = --generate-dependencies

#######################################################################
# create directory

$(BUILDDIR):
	mkdir $@

$(BUILDDIR)/${OUTPUTDIR}: $(BUILDDIR)
	mkdir -p $(BUILDDIR)/$(OUTPUTDIR)

$(BUILDDIR)/${OUTPUTDIR}/$(DEPDIR): $(BUILDDIR)/${OUTPUTDIR} $(BUILDDIR)
	mkdir -p $(BUILDDIR)/$(OUTPUTDIR)/$(DEPDIR)

$(info Build dir: $(BUILDDIR))
$(info Output dir: $(BUILDDIR)/${OUTPUTDIR})
$(info Dependencies dir: $(BUILDDIR)/${OUTPUTDIR}/$(DEPDIR))
#######################################################################
# build objects

INCLUDE_DIRS = -I$(SHARED_MODULES_INCLUDE) -I$(THOR_INCLUDE) -I$(LOCAL_INCLUDE)


#######################################################################
# build objects
# CUDA files

$(BUILDDIR)/${OUTPUTDIR}/%.o: %.cu $(BUILDDIR)/${OUTPUTDIR}/${DEPDIR}/%.d  | $(BUILDDIR)/${OUTPUTDIR}/$(DEPDIR) $(BUILDDIR)/$(OUTPUTDIR) $(BUILDDIR) 
	@echo -e '$(BLUE)creating dependencies for $@ $(END)'
	set -e; rm -f $(BUILDDIR)/${OUTPUTDIR}/${DEPDIR}/$*.d; \
	$(CC) $(cuda_dependencies_flags) $(arch)  $(cuda_flags) $(h5include) $(INCLUDE_DIRS)   -I$(includedir) $(CDB) $< > $(BUILDDIR)/${OUTPUTDIR}/${DEPDIR}/$*.d.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $(BUILDDIR)/${OUTPUTDIR}/${DEPDIR}/$*.d.$$$$ > $(BUILDDIR)/${OUTPUTDIR}/${DEPDIR}/$*.d; \
	rm -f $(BUILDDIR)/${OUTPUTDIR}/${DEPDIR}/$*.d.$$$$
	@echo -e '$(YELLOW)creating object file for $@ $(END)'
	$(CC) $(CC_comp_flag) $(arch)  $(cuda_flags) $(h5include) -I$(includedir) $(INCLUDE_DIRS) $(CDB) -o $@ $<

# C++ files
$(BUILDDIR)/${OUTPUTDIR}/%.o: %.cpp $(BUILDDIR)/${OUTPUTDIR}/${DEPDIR}/%.d  | $(BUILDDIR)/${OUTPUTDIR}/$(DEPDIR) $(BUILDDIR)/$(OUTPUTDIR) $(BUILDDIR)
	@echo -e '$(BLUE)creating dependencies for $@ $(END)'
	set -e; rm -f $(BUILDDIR)/${OUTPUTDIR}/${DEPDIR}/$*.d; \
	$(CC) $(dependencies_flags) $(arch) $(cpp_flags) $(h5include) $(h5include) $(INCLUDE_DIRS) -I$(includedir) $(CDB) $< > $(BUILDDIR)/${OUTPUTDIR}/${DEPDIR}/$*.d.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $(BUILDDIR)/${OUTPUTDIR}/${DEPDIR}/$*.d.$$$$ > $(BUILDDIR)/${OUTPUTDIR}/${DEPDIR}/$*.d; \
	rm -f $(BUILDDIR)/${OUTPUTDIR}/${DEPDIR}/$*.d.$$$$
	@echo -e '$(YELLOW)creating object file for $@  $(END)'
	$(CC) $(CC_comp_flag) $(arch) $(cpp_flags) $(h5include) $(INCLUDE_DIRS) -I$(includedir) $(CDB) -o $@ $<

libalfrodull.a: $(addprefix $(BUILDDIR)/$(OUTPUTDIR)/,$(obj)) | $(BUILDDIR)/$(OUTPUTDIR) $(BUILDDIR)
	@echo -e '$(YELLOW)creating $@ $(END)'
	@echo -e '$(GREEN)Linking Class into static lib $(END)'
	ar rcs $@ $(addprefix $(BUILDDIR)/$(OUTPUTDIR)/,$(obj))

DEPFILES := $(obj:%.o=$(BUILDDIR)/$(OUTPUTDIR)/$(DEPDIR)/%.d)
$(DEPFILES):

#######################################################################
# Cleanup
.phony: clean,ar
clean:
	@echo -e '$(CYAN)clean up library $(END)'
	-$(RM) libalfrodull.a
	@echo -e '$(CYAN)clean up modules objects $(END)'
	-$(RM) $(BUILDDIR)/debug/*.o
	-$(RM) $(BUILDDIR)/debug/*.o.json
	-$(RM) $(BUILDDIR)/release/*.o
	-$(RM) $(BUILDDIR)/release/*.o.json
	-$(RM) $(BUILDDIR)/prof/*.o
	-$(RM) $(BUILDDIR)/prof/*.o.json
	@echo -e '$(CYAN)clean up dependencies $(END)'
	-$(RM) $(BUILDDIR)/debug/$(DEPDIR)/*.d
	-$(RM) $(BUILDDIR)/debug/$(DEPDIR)/*.d.*
	-$(RM) $(BUILDDIR)/release/$(DEPDIR)/*.d
	-$(RM) $(BUILDDIR)/release/$(DEPDIR)/*.d.*
	-$(RM) $(BUILDDIR)/prof/$(DEPDIR)/*.d
	-$(RM) $(BUILDDIR)/prof/$(DEPDIR)/*.d.*
	-$(RM) -d $(BUILDDIR)/debug/$(DEPDIR)/
	-$(RM) -d $(BUILDDIR)/release/$(DEPDIR)/
	-$(RM) -d $(BUILDDIR)/prof/$(DEPDIR)/
	-$(RM) -d $(BUILDDIR)/debug/
	-$(RM) -d $(BUILDDIR)/release/
	-$(RM) -d $(BUILDDIR)/prof/
	@echo -e '$(CYAN)remove modules object dir $(END)'
	-$(RM) -d $(BUILDDIR)
$(info -------------------------------------------------------------- )
