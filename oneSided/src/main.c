#include "prototypes.h"  /* function headers */
#include "vars.h"        /* global variables */

int main(int argc, char **argv)
{
  mesh msh;
  char outFile[256];

  if(readPfile(argv[1],&pv)); else errMsg(1); /* reading the parameter file                              */
  init(argc, argv);                           /* initilize main run parameters and read stuff from files */
  initConstantsAndAbundances();               /* initialize some arrays and constants                    */
  setMaxAllowedDensities();                   /* set the maximum allowed densities                       */
  initMesh( &msh );                           /* initialize the mesh, the first slab for now             */

  find_mesh_equilibrium(&msh);

  sprintf(outFile, "%smesh.dat",pv.outputDir);
  writeDataOutputBinary(&msh, chemNet.nSpec, outFile);
  if( verbose & VERB1_MASK) fprintf(outputFd, "%sDumped the mesh to file %s\n", INDENT1, outFile);
  
  finalize(&msh);

  fprintf(stdout, " total number of errors = %d\n", nErrors);
  
  return 0;
}
