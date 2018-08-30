# makefile for phy modules, called from main makefile
# must create phy_modules.o in MODULES_DIR from THOR_ROOT

$(info Modules RT Makefile )

$(info Some variables inherited from parent makefiles)
$(info includes $(includedir))
$(info h5includes $(h5include))
$(info flags $(flags))
$(info arch $(arch))
$(info MODULES_DIR: $(MODULES_DIR))

$(info Sub Makefile variables)

THOR_ROOT = ../../
MODULES_INCLUDE = $(THOR_ROOT)src/phy_modules/inc/


BUILDDIR = $(THOR_ROOT)$(MODULES_DIR)

link_flags = --compiler-options -Wall,-std=c++11,-flinker-output=rel

$(info Linker Flags: $(link_flags))

all: $(BUILDDIR)/libphy_modules.a

# path to local module code
vpath %.cu src
vpath %.cpp src
vpath %.h inc
# path to thor headers
vpath %.h $(THOR_ROOT)src/headers/
# path to phy_modules
vpath %.h $(THOR_ROOT)src/phy_modules/inc/
vpath %.cu $(THOR_ROOT)src/phy_modules/src/
vpath %.cpp $(THOR_ROOT)src/phy_modules/src/
#######################################################################
# build objects



$(BUILDDIR)/radiative_transfer.o: radiative_transfer.cu radiative_transfer.h phy_modules.h phy_module_base.h 
	@echo $(YELLOW)creating $@ $(END)
	$(CC) $(arch)  $(flags) $(h5include) $(h5libdir)  -I$(THOR_ROOT)$(includedir) -I$(MODULES_INCLUDE) -Iinc -dc -o $@ $<



$(BUILDDIR)/phy_modules_main.o: phy_modules.cu phy_modules.h phy_module_base.h 
	@echo $(YELLOW)Building Module $(END)
	$(CC) $(arch)  $(flags) $(h5include) $(h5libdir) -I$(THOR_ROOT)$(includedir) -I$(MODULES_INCLUDE) -dc -o $(BUILDDIR)/phy_modules_main.o src/phy_modules.cu


$(BUILDDIR)/libphy_modules.a: $(BUILDDIR)/phy_modules_main.o $(BUILDDIR)/radiative_transfer.o 
	@echo $(YELLOW)Linking Modules into static lib $(END)
	ar rcs $(BUILDDIR)/libphy_modules.a $(BUILDDIR)/phy_modules_main.o $(BUILDDIR)/radiative_transfer.o 

#	$(CC) $(arch) $(link_flags) $(h5include) $(h5libdir) -I$(THOR_ROOT)$(includedir) -I$(MODULES_INCLUDE) -dc -o $(BUILDDIR)/phy_modules.o $(BUILDDIR)/phy_modules_main.o $(BUILDDIR)/radiative_transfer.o 
