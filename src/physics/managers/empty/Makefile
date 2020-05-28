# makefile for phy modules tempalte, called from main makefile
# must create libphy_modules.a in its root directory

$(info Modules Template Makefile )

# set some variables if not set 
includedir ?= unset
h5include ?= unset
cpp_flags ?= unset
cuda_flags ?= unset
arch ?= unset

$(info Some variables inherited from parent makefiles)
$(info includes: $(includedir))
$(info h5includes: $(h5include))
$(info cpp_flags: $(cpp_flags))
$(info cuda_flags: $(cuda_flags))
$(info arch: $(arch))
$(info CC compile flag: $(CC_comp_flag))
$(info Sub Makefile variables)

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


ifndef VERBOSE
.SILENT:
endif

#######################################################################
# create directory
$(BUILDDIR):
	mkdir $@
#######################################################################
# build objects

INCLUDE_DIRS = -I$(SHARED_MODULES_INCLUDE) -I$(THOR_INCLUDE) -I$(LOCAL_INCLUDE)

$(BUILDDIR)/phy_modules.o: phy_modules.cu phy_modules.h phy_module_base.h | $(BUILDDIR)
	@echo -e '$(YELLOW)creating $@ $(END)'
	# test if we build with compilation database (clang) or not (nvcc)	
	if test $$CDB = "-MJ" ; then \
		$(CC) $(CC_comp_flag) $(arch)  $(cuda_flags) $(h5include) $(INCLUDE_DIRS) $(CDB) $@.json -o $@ $<; \
	else \
		$(CC) $(CC_comp_flag) $(arch)  $(cuda_flags) $(h5include) $(INCLUDE_DIRS) -o $@ $<; \
	fi


libphy_modules.a: $(BUILDDIR)/phy_modules.o | $(BUILDDIR)
	@echo -e '$(YELLOW)creating $@ $(END)'
	@echo -e '$(GREEN)Linking Modules into static lib $(END)'
	ar rcs $@ $(BUILDDIR)/phy_modules.o 

#######################################################################
# Cleanup
.phony: clean,ar
clean:
	@echo -e '$(CYAN)clean up library $(END)'
	-$(RM) libphy_modules.a
	@echo -e '$(CYAN)clean up modules objects $(END)'
	-$(RM) $(BUILDDIR)/phy_modules.o
	-$(RM) $(BUILDDIR)/phy_modules.o.json
	@echo -e '$(CYAN)remove modules object dir $(END)'
	-$(RM) -d $(BUILDDIR)
