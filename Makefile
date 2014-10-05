AMUSE_DIR=~/amuse

.PHONY: clean

nothing:

clean:
	rm -rvf *.pyc *~ .*~ *.log

pdr_install:
	cp -fvr oneSided ${AMUSE_DIR}/src/amuse/community/pdr
	cd ${AMUSE_DIR}/src/amuse/community/pdr && make
pdr_clean:
	@rm -fvr ${AMUSE_DIR}/src/amuse/community/pdr
