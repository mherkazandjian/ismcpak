AMUSE_DIR=~/amuse

.PHONY: clean

nothing:

pdr_install:
	cp -fvr oneSided ${AMUSE_DIR}/src/amuse/community/pdr
	cd ${AMUSE_DIR}/src/amuse/community/pdr && make
pdr_clean:
	@rm -fvr ${AMUSE_DIR}/src/amuse/community/pdr

singularity:
	@singularity build container.sif Singularity

jupyterlab:
	@singularity exec --scratch /run/user container.sif jupyter-lab

clean:pdr_clean
	rm -rvf *.pyc *~ .*~ *.log 
