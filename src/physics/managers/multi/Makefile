# makefile for RT physics modules, called from main makefile
# must create libphy_modules.a in its root directory

$(info Modules RT Makefile )

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
$(info arch: $(arch))
$(info MODE: $(MODE))

######################################################################
# Directories
THOR_ROOT = ../../../../

# Includes
LOCAL_INCLUDE = inc

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


ifndef VERBOSE
.SILENT:
endif

######################################################################
$(info Sub Makefile variables)
$(info THOR root from submakefile: $(THOR_ROOT))

######################################################################
all: libphy_modules.a

# path to local module code
vpath %.cu $(LOCAL_SOURCES)
vpath %.cpp $(LOCAL_SOURCES)
vpath %.h $(LOCAL_INCLUDE)
# path to thor headers
vpath %.h $(THOR_INCLUDE)
# path to phy_modules
vpath %.h $(SHARED_MODULES_INCLUDE)
vpath %.cu $(SHARED_MODULES_DIR)
vpath %.cpp $(SHARED_MODULES_DIR)


#######################################################################
# create directory

$(BUILDDIR):
	mkdir $@

$(BUILDDIR)/${OUTPUTDIR}: $(BUILDDIR)
		mkdir -p $(BUILDDIR)/$(OUTPUTDIR)

#######################################################################
# build objects

INCLUDE_DIRS = -I$(SHARED_MODULES_INCLUDE) -I$(THOR_INCLUDE) -I$(LOCAL_INCLUDE)

$(BUILDDIR)/${OUTPUTDIR}/chemistry.o: chemistry.cu chemistry.h phy_modules.h phy_module_base.h thor_chemistry.h | $(BUILDDIR)
	@echo -e '$(YELLOW)creating $@ $(END)'
	if test $$CDB = "-MJ" ; then \
		$(CC) $(CC_comp_flag) $(arch) $(cuda_flags) $(h5include) $(h5libdir) $(INCLUDE_DIRS) $(CDB) $@.json -o $@ $<; \
	else \
		$(CC) $(CC_comp_flag) $(arch) $(cuda_flags) $(h5include) $(h5libdir) $(INCLUDE_DIRS) -o $@ $<; \
	fi

$(BUILDDIR)/${OUTPUTDIR}/radiative_transfer.o: radiative_transfer.cu radiative_transfer.h profx_RT.h phy_modules.h phy_module_base.h | $(BUILDDIR)/${OUTPUTDIR}
	@echo -e '$(YELLOW)creating $@ $(END)'
	if test $$CDB = "-MJ" ; then \
		$(CC) $(CC_comp_flag) $(arch)  $(cuda_flags) $(h5include) $(h5libdir) $(INCLUDE_DIRS) $(CDB) $@.json -o $@ $<; \
	else \
		$(CC) $(CC_comp_flag) $(arch)  $(cuda_flags) $(h5include) $(h5libdir) $(INCLUDE_DIRS) -o $@ $<; \
	fi

$(BUILDDIR)/${OUTPUTDIR}/boundary_layer.o: boundary_layer.cu boundary_layer.h phy_modules.h phy_module_base.h | $(BUILDDIR)/${OUTPUTDIR}
	@echo -e '$(YELLOW)creating $@ $(END)'
	if test $$CDB = "-MJ" ; then \
		$(CC) $(CC_comp_flag) $(arch)  $(cuda_flags) $(h5include) $(h5libdir) $(INCLUDE_DIRS) $(CDB) $@.json -o $@ $<; \
	else \
		$(CC) $(CC_comp_flag) $(arch)  $(cuda_flags) $(h5include) $(h5libdir) $(INCLUDE_DIRS) -o $@ $<; \
	fi

$(BUILDDIR)/${OUTPUTDIR}/phy_modules.o: phy_modules.cu phy_modules.h phy_module_base.h | $(BUILDDIR)/${OUTPUTDIR}
	@echo -e '$(YELLOW)creating $@ $(END)'
	if test $$CDB = "-MJ" ; then \
		$(CC) $(CC_comp_flag) $(arch)  $(cuda_flags) $(h5include) $(INCLUDE_DIRS) $(CDB) $@.json -o $@ $<; \
	else \
		$(CC) $(CC_comp_flag) $(arch)  $(cuda_flags) $(h5include) $(INCLUDE_DIRS) -o $@ $<; \
	fi


libphy_modules.a: $(BUILDDIR)/${OUTPUTDIR}/chemistry.o $(BUILDDIR)/${OUTPUTDIR}/radiative_transfer.o $(BUILDDIR)/${OUTPUTDIR}/boundary_layer.o $(BUILDDIR)/${OUTPUTDIR}/phy_modules.o | $(BUILDDIR)
	@echo -e '$(YELLOW)creating $@ $(END)'
	@echo -e '$(GREEN)Linking Modules into static lib $(END)'
	ar rcs $@ $(BUILDDIR)/${OUTPUTDIR}/phy_modules.o $(BUILDDIR)/${OUTPUTDIR}/radiative_transfer.o $(BUILDDIR)/${OUTPUTDIR}/boundary_layer.o $(BUILDDIR)/${OUTPUTDIR}/chemistry.o

#######################################################################
# Cleanup
.phony: clean,ar
clean:
	@echo -e '$(CYAN)clean up library $(END)'
	-$(RM) libphy_modules.a
	@echo -e '$(CYAN)clean up modules objects $(END)'
	-$(RM) $(BUILDDIR)/debug/*.o
	-$(RM) $(BUILDDIR)/debug/*.o.json
	-$(RM) $(BUILDDIR)/release/*.o
	-$(RM) $(BUILDDIR)/release/*.o.json
	-$(RM) $(BUILDDIR)/prof/*.o
	-$(RM) $(BUILDDIR)/prof/*.o.json
	-$(RM) -d $(BUILDDIR)/debug/
	-$(RM) -d $(BUILDDIR)/release/
	-$(RM) -d $(BUILDDIR)/prof/
	@echo -e '$(CYAN)remove modules object dir $(END)'
	-$(RM) -d $(BUILDDIR)
