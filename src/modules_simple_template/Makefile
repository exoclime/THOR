# makefile for phy modules, called from main makefile

$(info Modules Makefile )
$(info includes $(includedir))
$(info h5includes $(h5include))
$(info flags $(flags))
$(info arch $(arch))

all: phy_modules.o

THOR_ROOT = ../../
MODULES_INCLUDE = inc

vpath %.cu src

#######################################################################
# build objects
# CUDA files
# $(OBJDIR)/${OUTPUTDIR}/%.o: %.cu $(OBJDIR)/$(OUTPUTDIR)/%.d| $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
# 	@echo $(YELLOW)creating $@ $(END)
# 	$(CC) $(arch)  $(flags) $(h5include) $(h5libdir) -I$(includedir) -dc -o $@ $<

# # C++ files
# $(OBJDIR)/${OUTPUTDIR}/%.o: %.cpp $(OBJDIR)/$(OUTPUTDIR)/%.d| $(OBJDIR)/$(OUTPUTDIR) $(OBJDIR)
# 	@echo $(YELLOW)creating $@ $(END)
# 	$(CC) $(arch)  $(flags) $(h5include) $(h5libdir) -I$(includedir) -dc -o $@ $<



phy_modules.o: phy_modules.cu | $(THOR_ROOT)$(MODULES_DIR)
	@echo "Building Modules"
	$(CC) $(arch)  $(flags) $(h5include) $(h5libdir) -I$(THOR_ROOT)$(includedir) -I$(MODULES_INCLUDE) -dc -o $(THOR_ROOT)$(MODULES_DIR)/$@ $<
