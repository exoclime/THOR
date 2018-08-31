# makefile for phy modules, called from main makefile
# must create phy_modules.o in MODULES_DIR from THOR_ROOT

$(info Modules RT Makefile )

# set some variables if not set 

$(info Some variables inherited from parent makefiles)
$(info includes $(includedir))
$(info h5includes $(h5include))
$(info flags $(flags))
$(info arch $(arch))

$(info Sub Makefile variables)

######################################################################
# Directories
THOR_ROOT = ../../

# Includes
LOCAL_INCLUDE = inc

# shared modules
SHARED_MODULES_INCLUDE = $(THOR_ROOT)src/phy_modules/inc/

# thor root include if we want to use code from there
THOR_INCLUDE = $(THOR_ROOT)src/headers

# source dirs
LOCAL_SOURCES = src
SHARED_MODULES_DIR = $(THOR_ROOT)src/phy_modules/src/

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


#######################################################################
# create directory
$(BUILDDIR):
	mkdir $@
#######################################################################
# build objects

INCLUDE_DIRS = -I$(SHARED_MODULES_INCLUDE) -I$(THOR_INCLUDE) -I$(LOCAL_INCLUDE)

$(BUILDDIR)/radiative_transfer.o: radiative_transfer.cu radiative_transfer.h phy_modules.h phy_module_base.h | $(BUILDDIR)
	@echo $(YELLOW)creating $@ $(END)
	$(CC) $(arch)  $(flags) $(h5include) $(h5libdir) $(INCLUDE_DIRS) -dc -o $@ $<


$(BUILDDIR)/phy_modules.o: phy_modules.cu phy_modules.h phy_module_base.h | $(BUILDDIR)
	@echo $(YELLOW)creating $@ $(END)
	$(CC) $(arch)  $(flags) $(h5include) $(h5libdir) $(INCLUDE_DIRS) -dc -o $@ $<


libphy_modules.a: $(BUILDDIR)/radiative_transfer.o $(BUILDDIR)/phy_modules.o | $(BUILDDIR)
	@echo $(YELLOW)creating $@ $(END)
	@echo $(GREEN)Linking Modules into static lib $(END)
	ar rcs $@ $(BUILDDIR)/phy_modules.o $(BUILDDIR)/radiative_transfer.o 

#######################################################################
# Cleanup
.phony: clean,ar
clean:
	@echo $(CYAN)clean up library $(END)
	-$(RM) libphy_modules.a
	@echo $(CYAN)clean up objects $(END)
	-$(RM) $(BUILDDIR)/phy_modules_main.o
	-$(RM) $(BUILDDIR)/radiative_transfer.o 
	@echo $(CYAN)remove object dir $(END)
	-$(RM) -d $(BUILDDIR)
