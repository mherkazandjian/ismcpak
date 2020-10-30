[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/2343)

Quickstart
==========
To run a quick example, the following container can be used:

    $ git clone -b alpha-master https://github.com/mherkazandjian/ismcpak.git ~/ismcpak
    $ cd ~/ismcpak/tests
    $ singularity exec library://mher/default/ismcpak:master mpiexec python run_singleMesh.py 

This is a package which implements some utilities useful for modelling and
analyzing simulation output of PDRs.

Build the container locally
---------------------------

The following command build the singularity container on a local machine. The
only prerequisite is to have singularity installed and to have sudo access. 

    $ sudo make singularity

Prerequisites
=============
amuse  - mpich
PyQt4
ipython


Installing the PDR code
=======================
The PDR code should be copied into:
      amuse/src/amuse/community/pdr

Compiling the PDR code
======================
The PDR code can be compiled using:
    ~> cd amuse/src/amuse/community/pdr
    ~> make all
The generates the libpdr.a library

Setting up the working environment
==================================
The path to ismcpak should be added to the PYTHONPATH environment variable. For
bash, the following line should be added:

    export PYTHONPATH=/PATH/TO/ismcpak:$PYTHONPATH
   
to tcsh :

    setenv PYTHONPATH /PATH/TO/ismcpak:$PYTHONPATH

Running the code
================
  Basic test - single model
  -------------------------
  The PDR code can only be run through the AMUSE ( http://amusecode.org ).
  Depending on the mpi environment installed with AMUSE, it might be 
  necessary to launch the mpd deamon before executing either:
  
    ~> mpirun -np 1 python run_singleMesh.py
    
  or via
  
    ~> python run_singleMesh.py

Running a Grid of models
========================

# setup the working environment variables
1) source setdev

# install the pdr code into amuse (make sure the correct 
# path of amuse is set in setenv
2) make pdr_install

# after these two steps, the tests 
   run_singleMesh.py
   chemical_network_pdr_code.py
should run without errors


3) to run a grid, use the following under ismcpak:
   ~>  ipython --pylab=qt tests/run_oneSidedGrid.py

4) after the model data is written to 
      tests/oneSidedGrid/meshes
   we need to construct the database files .db using constructReadArchive.py
   ~> ipython --pylab=qt constructReadArchive.py

   after the database is constructed we must have the file 
         meshes.db  meshes.db.info
   in the output directory and a message 
         archive integrity test passed
   must be displayed

5) after creating the database, a reference file must be generated which 
   stores information about the parameters which have been used in 
   generating the data. A template of this file is located under
        runs/tests/templateDir/used_params.py
   where the parameters used by run_oneSidedGrid.py should be filled in
   by hand. Once the values are changed :
       ~> python used_parms.py
   generates the pickle file

6) set the desired display parameters in analyzeArchive.py and invoke :
     ~> ipython --pylab=qt analyzeArchive.py

7) to generate the radex databases, the bottom part of analyzeArchive.py should be enabled to
   allow radex databases to be computed and written do disk. Set the desired values of 
   Av to compute and the species whose emission will be computed and re-run: 
     ~> ipython --pylab=qt analyzeArchive.py
   As a check, the data in 
        tests/oneSidedGrid/radexDbs
   should have directories with the Avs we have set and each directory should 
   have files for each species we have specified.

8) after producing the radex database files, we can convert that data to ascii data using :
     ~> ipython ismcpak2Ascii.py   
   
Disclaimer
==========
THIS SOFTWARE IS PROVIDED UNDER THE GPL LICENSE BY THE COPYRIGHT HOLDERS AND 
CONTRIBUTORS “AS IS” AND DOES NOT EXPRESS OR PROVIDE IMPLIED WARRANTIES, 
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND F
ITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT
, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY 
WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH 
DAMAGE. 

See LICENSE.txt for more information about the GPL license.

Please cite the following papers if any part of this package is used in your 
research. 

   Kazandjian, M. V., Meijerink, R., Pelupessy, I., Israel, F. P., Spaans, M.,
   arXiv:1403.7000

   Kazandjian, M. V., Meijerink, R., Pelupessy, I., Israel, F. P., Spaans, M.,
   2012, A&A, 542, A65, 26

   Meijerink, R., Spaans, M., & Israel, F. P. 2007, A&A, 461, 793

   Meijerink, R. & Spaans, M. 2005, A&A, 436, 397


 
