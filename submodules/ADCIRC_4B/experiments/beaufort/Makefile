clean:
	rm -f fort.63 fort.64 maxele.63 maxvel.63

rose-clean: clean
	# remove rose generated files
	rm -rf prose_logs prose_workspace *postprocessed* *.rmod rose_*

reset: rose-clean
	rm ../../work_rose/odir4/* -rf
	rm ../../work_rose/*.f -f
	rm ../../work_rose/*.F -f
	rm ../../work_rose/*.rmod -f
	rm ../../work/odir4/* -rf
	git checkout ../../src
	git checkout ../../thirdparty
	git checkout ../../wind