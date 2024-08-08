FC = gfortran
FFLAGS = -cpp -fdefault-real-8

default: main

main: main.f90 target_module

	$(FC) $(FFLAGS) -o main main.f90 target_module.o

target_module: target_module.f90
	$(FC) $(FFLAGS) -c target_module.f90

clean:
	rm -rf *.mod *.o main

reset: clean
	rm  -rf *.rmod prose_* *postprocessed*
	git checkout *90