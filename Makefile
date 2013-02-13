.PHONY: doc clean

doc:
	cd doc_source && make html && cd ..
clean:
	rm -rvf *.pyc
	cd doc_source && make clean && cd ..

