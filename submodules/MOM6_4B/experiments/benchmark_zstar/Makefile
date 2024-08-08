setup:
	cp -r ../../prose_preprocessed_src/* ../../src/MOM6/src/

clean:
	rm -f *.nc *.nc.* MOM_parameter_doc.* mom6.standalone.o* ocean_solo.res exitcode available_diags.000000 ocean.stats time_stamp.out logfile.000000.out CPU_stats

rose-clean: clean
	# remove rose generated files
	rm -rf prose_logs prose_workspace *postprocessed* *.rmod rose_*

reset: rose-clean
	rm ../../build/* -rf
	cd ../../src/MOM6 && git clean -i . && git checkout --recurse-submodules . && cd -
	cd ../../src/FMS && git clean -i . && git checkout --recurse-submodules . && cd -
	cd ../../src/MOM6/pkg/CVMix-src && git clean -i . && cd -