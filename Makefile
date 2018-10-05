AMUSE_DIR=~/amuse

.PHONY: clean

nothing:

pdr_install:
	cp -fvr oneSided ${AMUSE_DIR}/src/amuse/community/pdr
	cd ${AMUSE_DIR}/src/amuse/community/pdr && make
pdr_clean:
	@rm -fvr ${AMUSE_DIR}/src/amuse/community/pdr

clean:pdr_clean
	rm -rvf *.pyc *~ .*~ *.log 

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

