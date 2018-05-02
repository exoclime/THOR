# Builds THOR executable

sm:=50 # Streaming Multiprocessor version
arch := -arch sm_$(sm)

path_hd  := $(shell pwd)/src/headers
path_src := $(shell pwd)/src
obj   := esp.o grid.o esp_initial.o planet.o thor_driver.o profx_driver.o esp_output.o
headers := $(path_hd)/define.h $(path_hd)/grid.h $(path_hd)/planet.h $(path_hd)/esp.h \
           $(path_hd)/dyn/thor_fastmodes.h $(path_hd)/dyn/thor_adv_cor.h $(path_hd)/dyn/thor_auxiliary.h \
           $(path_hd)/dyn/thor_vertical_int.h $(path_hd)/dyn/thor_slowmodes.h $(path_hd)/dyn/thor_diff.h \
           $(path_hd)/dyn/thor_div.h $(path_hd)/phy/profx_auxiliary.h $(path_hd)/phy/profx_held_suarez.h
flag := -g -G

vvv := $(shell grep '^prefix' `which h5cc`)
$(eval $(vvv))

vvv1 := $(shell grep '^exec_prefix' `which h5cc`)
$(eval $(vvv1))

vvv2 := $(shell grep '^libdir' `which h5cc`)
$(eval $(vvv2))
h5lib = -L$(libdir) -lhdf5

vvv3 := $(shell grep '^includedir' `which h5cc`)
$(eval $(vvv3))

includehdf = -I$(includedir)
# Path where the hdf5 lib was installed,
# should contain `include` and `lib`
# hdf_dir := <hdf5 directory>
#
# h5lib = -L$(hdf_dir)/lib -lhdf5
# includehdf = -I$(hdf_dir)/include

esp: $(obj) $(headers)
	nvcc $(arch) $(includehdf) $(h5lib) -o bin/esp $(obj) $(flag)
	mv *.o obj/

grid.o: $(path_src)/grid/grid.cu
	nvcc $(arch) $(includehdf) $(h5lib) -dc $(path_src)/grid/grid.cu $(flag)

esp_initial.o: $(path_src)/initial/esp_initial.cu
	nvcc $(arch) $(includehdf) $(h5lib) -dc $(path_src)/initial/esp_initial.cu $(flag)

planet.o: $(path_src)/initial/planet.cu
	nvcc $(arch) $(includehdf) $(h5lib) -dc $(path_src)/initial/planet.cu $(flag)

thor_driver.o: $(path_src)/thor/thor_driver.cu
	nvcc $(arch) $(includehdf) $(h5lib) -dc $(path_src)/thor/thor_driver.cu $(flag)

profx_driver.o: $(path_src)/profx/profx_driver.cu
	nvcc $(arch) $(includehdf) $(h5lib) -dc $(path_src)/profx/profx_driver.cu $(flag)

esp.o: $(path_src)/esp.cu
	nvcc $(arch) $(includehdf) $(h5lib) -dc $(path_src)/esp.cu $(flag)

esp_output.o: $(path_src)/output/esp_output.cu
	nvcc $(arch) $(includehdf) $(h5lib) -dc $(path_src)/output/esp_output.cu $(flag)

.phony: clean, ar
clean:
	rm bin/esp obj/esp.o obj/grid.o obj/esp_initial.o obj/thor_driver.o \
	                     obj/profx_driver.o obj/esp_output.o obj/planet.o
