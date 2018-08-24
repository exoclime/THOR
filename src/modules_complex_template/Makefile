# makefile for phy modules, called from main makefile

$(info Modules Makefile )
$(info includes $(includedir))
$(info h5includes $(h5include))
$(info flags $(flags))
$(info arch $(arch))


THOR_ROOT = ../../
MODULES_INCLUDE = inc

link_flags = --compiler-options -Wall,-std=c++11,-flinker-output=rel

all: $(THOR_ROOT)$(MODULES_DIR)/phy_modules.o

vpath %.cu src
vpath %.cpp src
vpath %.h inc
vpath %.h $(THOR_ROOT)src/headers/
#######################################################################
# build objects



$(THOR_ROOT)$(MODULES_DIR)/radiative_transfer.o: radiative_transfer.cu radiative_transfer.h phy_modules.h phy_module_base.h 
	@echo $(YELLOW)creating $@ $(END)
	$(CC) $(arch)  $(flags) $(h5include) $(h5libdir)  -I$(THOR_ROOT)$(includedir) -I$(MODULES_INCLUDE) -dc -o $@ $<



$(THOR_ROOT)$(MODULES_DIR)/phy_modules_main.o: phy_modules.cu phy_modules.h phy_module_base.h 
	@echo "Building Module"
	$(CC) $(arch)  $(flags) $(h5include) $(h5libdir) -I$(THOR_ROOT)$(includedir) -I$(MODULES_INCLUDE) -dc -o $(THOR_ROOT)$(MODULES_DIR)/phy_modules_main.o src/phy_modules.cu


$(THOR_ROOT)$(MODULES_DIR)/phy_modules.o: $(THOR_ROOT)$(MODULES_DIR)/phy_modules_main.o $(THOR_ROOT)$(MODULES_DIR)/radiative_transfer.o 
	@echo "Linking Modules"
	$(CC) $(arch) $(link_flags) $(h5include) $(h5libdir) -I$(THOR_ROOT)$(includedir) -I$(MODULES_INCLUDE) -dc -o $(THOR_ROOT)$(MODULES_DIR)/phy_modules.o $(THOR_ROOT)$(MODULES_DIR)/phy_modules_main.o $(THOR_ROOT)$(MODULES_DIR)/radiative_transfer.o 
