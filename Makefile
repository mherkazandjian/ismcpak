.PHONY: doc clean

doc:
	cd doc_source && make html && cd ..
clean:
	rm -rvf *.pyc *~ .*~ *.log
	cd doc_source && make clean && cd ..
	cd galaxies && make clean
	cd ismrad && make clean
	cd paper2 && make clean
	cd data-for-collaborators && make clean

