# specify where the installed copy of ROSE is located. 
# Essentially the --prefix path used with configure

## Your Plugin source files
Plugin=ProsePlugin
Plugin_SOURCE=$(Plugin).cpp

## Input testcode for your plugin
TESTCODE=test1.cpp

# Standard C++ compiler stuff (see rose-config --help)
comma   := ,
CXX      = $(shell $(ROSE_EXE_PATH)/rose-config cxx)
CPPFLAGS = $(shell $(ROSE_EXE_PATH)/rose-config cppflags) -I.
CXXFLAGS = $(shell $(ROSE_EXE_PATH)/rose-config cxxflags)
LIBDIRS  = $(shell $(ROSE_EXE_PATH)/rose-config libdirs)
LDFLAGS  = $(shell $(ROSE_EXE_PATH)/rose-config ldflags) -lboost_graph -L. \
		$(addprefix -Wl$(comma)-rpath -Wl$(comma), $(subst :, , $(LIBDIRS)))

#-------------------------------------------------------------
# Makefile Targets
#-------------------------------------------------------------

all: $(Plugin).so

# compile the plugin and generate a shared library
# -g is recommended to be used by default to enable debugging your code
$(Plugin).so: $(Plugin_SOURCE)
	$(CXX) -g $(Plugin_SOURCE) -Wno-deprecated-declarations -fpic -shared $(CPPFLAGS) $(LDFLAGS) -o $@

# test the plugin
check: $(Plugin).so
	$(ROSE_EXE_PATH)/rose-compiler -c -rose:plugin_lib $(Plugin).so -rose:plugin_action print-names -rose:plugin_arg_print-names op1 -I. -I$(ROSE_INSTALL)/include $(TESTCODE) 

clean:
	rm -rf $(Plugin).so
