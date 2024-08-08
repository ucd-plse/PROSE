ifneq ($(GPTL), )

## entire
# FFLAGS += -g -finstrument-functions
# FFLAGS := $(filter-out -O3, $(FFLAGS))
# FFLAGS += -O1
# CPPFLAGS += -g -finstrument-functions
# CPPFLAGS := $(filter-out -O3, $(CPPFLAGS))
# CPPFLAGS += -O1
# LDFLAGS += -L /glade/u/home/jdvanover/gptl-8.1.1/install-intel-pmpi/lib -lgptl -lgptlf -rdynamic
# mpas.o			: private override FFLAGS += -I /glade/u/home/jdvanover/gptl-8.1.1/install-intel-pmpi/include -DGPTL
# mpas.o			: private override FFLAGS := $(filter-out -finstrument-functions, $(FFLAGS))

## atm_time_integration
mpas_atm_time_integration.o 		: private override FFLAGS += -I /glade/u/home/jdvanover/gptl-8.1.1/install-intel/include -DGPTL
mpas_atm_core.o						: private override FFLAGS += -I /glade/u/home/jdvanover/gptl-8.1.1/install-intel/include -DGPTL
mpas_timer.o						: private override FFLAGS += -I /glade/u/home/jdvanover/gptl-8.1.1/install-intel/include -DGPTL
mpas_dmpar.o						: private override FFLAGS += -I /glade/u/home/jdvanover/gptl-8.1.1/install-intel/include -DGPTL
mpas_atm_iau.o						: private override FFLAGS += -I /glade/u/home/jdvanover/gptl-8.1.1/install-intel/include -DGPTL
mpas_log.o							: private override FFLAGS += -I /glade/u/home/jdvanover/gptl-8.1.1/install-intel/include -DGPTL
mpas_pool_routines.o				: private override FFLAGS += -I /glade/u/home/jdvanover/gptl-8.1.1/install-intel/include -DGPTL
mpas_atmphys_driver_microphysics.o	: private override FFLAGS += -I /glade/u/home/jdvanover/gptl-8.1.1/install-intel/include -DGPTL
mpas_atmphys_todynamics.o			: private override FFLAGS += -I /glade/u/home/jdvanover/gptl-8.1.1/install-intel/include -DGPTL
mpas_atm_boundaries.o				: private override FFLAGS += -I /glade/u/home/jdvanover/gptl-8.1.1/install-intel/include -DGPTL
mpas_timekeeping.o					: private override FFLAGS += -I /glade/u/home/jdvanover/gptl-8.1.1/install-intel/include -DGPTL
mpas_vector_reconstruction.o		: private override FFLAGS += -I /glade/u/home/jdvanover/gptl-8.1.1/install-intel/include -DGPTL
mpas.o	: private override FFLAGS += -I /glade/u/home/jdvanover/gptl-8.1.1/install-intel/include -DGPTL

endif