#include "src/prototypes.h"
#include "src/vars.h"
#include "interface.h"
#include <math.h>
#include <mpi.h>

#define MAX_NUMBER_MESHES 1000000

#define EXISTS      0x80
#define INITIALIZED 0x40
#define IS_MY_MESH  0x20
#define SUCCESS     0x10

mesh **meshes;
int nMeshes;
int highest_id;
unsigned short int *flags;
int *errFlags;
char outFile[256];
char command[256];
int foo;

initialConditions ics0;

int nProcs=-1, rank=-1;

double *rho_all, *T_init_all, *G0_FUV_all, *Lmech_all;
int *ids_all;

/* MPI recv buffers */
int    *recvBuffi;
double *recvBuffd;

MPI_Status mpStt;


/*----------------------------------------------------------------------------------------------------*/
/*   initialize : mpi environment, chemical network, read self sheilding and rovib tables             */
/*                all meshe parameter and flag buffers                                                */
/*----------------------------------------------------------------------------------------------------*/
int initialize()
{
   char outFname[256];
   int i;

   /* intializing MPI stuff */
   MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   printf("AMUSE: my rank is %d out of %d \n", rank, nProcs-1);

   sprintf(outFname, "%smesh.out-id-%d", pv.outputDir, rank);
   set_diagnostics_parameters( outFname );

   /* reading and initializing IO related stuff */
   chemical_network_data();   /* read and alloc all info about the chemical network (alloc and setup chemNet) */
   read_selfshlding_CO();     /* read and alloc stuff needed for sheilding                                    */
   read_vibrational_data();   /* read data of vibrational modes                                               */
   read_rotational_data();    /* read data of rotational modes                                                */

   /* making a copy of the guess abundances read in chemical_network_data() */
   MAKE_1D_ARRAY(ics0.abun   , chemNet.nBaseSpec+1, double); /* Elemental relative abundances */
   MAKE_1D_ARRAY(ics0.maxAbun, chemNet.nSpec+1    , double); /* max allowed abundances        */
   for(i = 0; i < chemNet.nBaseSpec+1; i++) ics0.abun[i] = ics.abun[i];
   for(i = 0; i < chemNet.nSpec+    1; i++) ics0.maxAbun[i] = ics.maxAbun[i];

   meshes   = (mesh **)malloc(MAX_NUMBER_MESHES*sizeof(mesh *));
   flags    = (unsigned short int *)malloc(MAX_NUMBER_MESHES*sizeof(unsigned short int));
   errFlags = (int *)malloc(MAX_NUMBER_MESHES*sizeof(int));

   recvBuffi = (int *  )malloc(MAX_NUMBER_MESHES * sizeof(int));
   recvBuffd = (double *)malloc(MAX_NUMBER_MESHES * sizeof(double));

   rho_all    = malloc(MAX_NUMBER_MESHES * sizeof(double));
   T_init_all = malloc(MAX_NUMBER_MESHES * sizeof(double));
   G0_FUV_all = malloc(MAX_NUMBER_MESHES * sizeof(double));
   Lmech_all  = malloc(MAX_NUMBER_MESHES * sizeof(double));
   ids_all    = malloc(MAX_NUMBER_MESHES * sizeof(int));

   nMeshes=0;
   for(i = 0; i < MAX_NUMBER_MESHES; i++) {
      flags[i] = 0x00;
      errFlags[i] = 0;
      meshes[i] = (mesh *)malloc(sizeof(mesh));
   }

   highest_id = 0;

   if(rank == 0)
      dumpVersionFile();    /* copying the version file to the output Directory */
   fflush(stdout);
   MPI_Barrier(MPI_COMM_WORLD);

   /* initializing the guess database in the G0, rho parameter space*/
   init_database();

   if(rank == 0) printf("AMUSE: all initialized\n");
   MPI_Barrier(MPI_COMM_WORLD);

   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
/*      adds all the meshes to the stack and allocates local buffers for meshes local to procceses    */
/*----------------------------------------------------------------------------------------------------*/
int add_mesh(int *id, double rho, double T_init, double G0_FUV, double Lmech)
{
   *id=highest_id++;
   //printf("AMUSE: proccess %d adding a new mesh to the que, id = %d\n", rank, *id);

   rho_all[*id]    = rho;
   T_init_all[*id] = T_init;
   G0_FUV_all[*id] = G0_FUV;
   Lmech_all[*id]  = Lmech;

   nMeshes++;
   flags[*id]=EXISTS;
   ids_all[*id]=*id;

   /*------------------------------------------------------------------------------*/
   /* setting up and allocating meshes on each proccess according to the mesh ID   */
   /* each proccess gets (allocates) a certain number of meshes. Here there values */
   /* for the first slab are also initialized                                      */
   /*------------------------------------------------------------------------------*/
   if( isMyMesh(*id) ) {
      flags[*id] = flags[*id] | IS_MY_MESH;
   }

   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int evolve_mesh( int i )
{
   FILE *fd;
   int j, k, idx;

   double hBefore=heapUse;


   double *AvTmp; // temporary array holding the ending positions in Av of all the slabs
   AvTmp = (double *)malloc(sizeof(double)*5000);

   meshes[i]->nx = descretizeSlab1( AvTmp );
   meshes[i]->nx += 1; // taking the zeroth slab into account

   setupAllocMesh( meshes[i], meshes[i]->nx, chemNet.nSpec);
   flags[i] = flags[i] | INITIALIZED;

   printf("AMUSE: process %d evolving mesh %d : log(G0) = %.2f , log(dens) = %.2f, log(Lmech) = %e\n",
         rank, i, G0_FUV_all[i], rho_all[i], Lmech_all[i]);

   pv.dens0      = pow(10e0, rho_all[i]);
   pv.G0         = pow(10e0, G0_FUV_all[i]);
   pv.gamma_mech = Lmech_all[i];
   pv.T_init     = T_init_all[i];

   nNewton=nBisect=0;
   tmr1.dt=tmr2.dt=tmr3.dt=tmr4.dt=tmr1.tot=tmr2.tot=tmr3.tot=tmr4.tot=0.0;
   nErrors=0;
   errCode=0;

   TIC(tmr1);

   T0     = 12.2 * pow(pv.G0, 0.2);
   tau100 = 2.7e2 * pv.G0 / pow(T0, 5e0) * pv.metalicity;

   meshes[i]->AvL[0] = 0.0; /* setting the position/size of the first slab to 0               */
   meshes[i]->xs[0]  = 0.0;
   meshes[i]->xc[0]  = 0.0;
   meshes[i]->xe[0]  = 0.0;
   meshes[i]->dx[0]  = 0.0;

   meshes[i]->NH2All_L[0] = 0.0; /* setting the initial colon densities due to the slabs on the  */
   meshes[i]->NH2_L[0]    = 0.0; /* left to zero                                                 */
   meshes[i]->NCO_L[0]    = 0.0;
   meshes[i]->NC_L[0]     = 0.0;
   meshes[i]->NH2O_L[0]   = 0.0;

   /* set the self sheilding factor of the first slabs on the right and left */
   meshes[i]->betaH2_L[0]= meshes[i]->betaCO_L[0] = meshes[i]->beta13CO_L[0] = 1.0;
   meshes[i]->betaH2_R[0]= meshes[i]->betaCO_R[0] = meshes[i]->beta13CO_R[0] = 1.0;

   for(j = 0; j < chemNet.nSpec; j++)
      for( k = 0; k < meshes[i]->nx; k++)
         meshes[i]->abun[j][k] = 0.0;

   /* setting the precomputed locations of Av from the temporary one */
   for( j = 1; j < meshes[i]->nx; j++)
      meshes[i]->AvL[j] = AvTmp[j-1];
   free(AvTmp);

   setInitialGuessAbundances( ); /* ;;; see if this is neccessarry.. and remove the stuff in it that are redundant */
   /*     after usuing the databse                                                   */

   errFlags[i] = find_mesh_equilibrium(meshes[i]);

   //sprintf(command, "mv %sflux-00.out %sflux-00.out-id-%02d", pv.outputDir, pv.outputDir, i);
   //system(command);

   if( errFlags[i]==0) {
      sprintf(outFile, "%smesh.dat-id-%06d", pv.outputDir, i);
      printf("AMUSE: dumping binary file to %s\n", outFile);
      writeDataOutputBinary(meshes[i], chemNet.nSpec, outFile );
      flags[i] = flags[i] | SUCCESS;
   } else {
      //  AMUSE_MESH_FAILED:
      printf("AMUSE ERROR: mesh %d on process %d failed - error Code %d\n", i, rank, errFlags[i]);
      printf("AMUSE ERROR: mesh parameters : id = %d rho0 = %f, G0 = %f, gamma_mech = %e\n", i, rho_all[i], G0_FUV_all[i], Lmech_all[i]);
   }

   finalize(meshes[i]);
   freeMesh(meshes[i], chemNet.nSpec );

   printf("########################################################################################\n");

   //printf("===========> heap use after  = %f MB\n", (double)heapUse/((double)(1024.0*1024.0)) );
   //printf("===========> heap leak       = %f MB\n", (double)(heapUse-hBefore)/((double)(1024.0*1024.0)) );

   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int calc_equilibrium() // modify this to check for load balancing
{
   int i;
   timer tWall;            /* tmr1: wall time from start of simulation */

   /* printing the share of meshes of each proccess */
   /*-----------------------------------------------*/
   fflush(stdout);
   MPI_Barrier(MPI_COMM_WORLD);
   if(rank!=0)        MPI_Recv(&foo, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &mpStt);
   printf("AMUSE: process %3d will evolve meshes : ", rank);
   for( i=0; i<MAX_NUMBER_MESHES; i++) {
      if( flags[i] & IS_MY_MESH ) {
         printf("%d ", ids_all[i]);
      }
   }
   printf("\n");
   fflush(stdout);
   if(rank!=nProcs-1) MPI_Send(&foo, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);

   fflush(stdout);
   MPI_Barrier(MPI_COMM_WORLD);
   /*-----------------------------------------------*/

   /* evolving the meshes */
   /*---------------------*/
   tWall.dt=0.0;

   TIC(tWall);

   for( i=0; i<MAX_NUMBER_MESHES; i++) {

      if( flags[i] &  IS_MY_MESH  ) {
         evolve_mesh(i);
      }
   }
   TOC(tWall);
   /*---------------------*/

   printf("AMUSE: Wall time evolving all meshes : %.2f s\n", tWall.dt);
   return 0;
}
/*----------------------------- get the error flags --------------------------------------------------*/
/* when called from amuse, errFlagsRet is already allocated, so we do not have to allocate it before  */
/* it is used                                                                                         */
/*----------------------------------------------------------------------------------------------------*/
int get_errorFlags(const int *mesh_ids, const int n, int *errFlagsRet)
{
   int i, idx;

   for(i=0; i<n; i++) {

      idx = mesh_ids[i];

      if( rank !=0 &&  (flags[idx] & IS_MY_MESH) )
         MPI_Send( &errFlags[idx], 1, MPI_INT, 0, idx, MPI_COMM_WORLD);

      if( rank==0 && (flags[idx] & EXISTS) ) {
         if ( flags[idx] & IS_MY_MESH )
            errFlagsRet[idx] = errFlags[idx];
         else
            MPI_Recv( &errFlagsRet[idx], 1, MPI_INT, MPI_ANY_SOURCE, idx, MPI_COMM_WORLD, &mpStt);
      }
   }

   return 0;
}
/*------------------------------get statisitcs on temperate-------------------------------------------*/
int get_temperature(const int *mesh_ids, const int n, double *temperatureRet)
{
   double tempSend=-2, stat;
   int i, idx;

   stat=0.0;

   for(i=0; i<n; i++) {

      idx = mesh_ids[i];

      if( rank !=0 &&  (flags[idx] & IS_MY_MESH) ) {

         if( errFlags[idx]==0 ) /* setting the temperature to be sent */
            tempSend = meshes[idx]->gasT[0];
         else
            tempSend = -1.0;

         MPI_Send( &tempSend, 1, MPI_DOUBLE, 0, idx, MPI_COMM_WORLD);
      }

      if( rank==0 && (flags[idx] & EXISTS) ) {

         if ( flags[idx] & IS_MY_MESH ) {
            if( errFlags[idx]==0 )
               temperatureRet[idx] = meshes[idx]->gasT[0];
            else
               temperatureRet[idx] = -1.0;
         }
         else
            MPI_Recv( &temperatureRet[idx], 1, MPI_DOUBLE, MPI_ANY_SOURCE, idx, MPI_COMM_WORLD, &mpStt);
      }
   }

   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int setInitialGuessAbundances( )
{
   int i;

   /* copyting original guesses */
   for(i=0; i<chemNet.nBaseSpec+1; i++) ics.abun[i]=ics0.abun[i];
   for(i=0; i<chemNet.nSpec+    1; i++) ics.maxAbun[i]=ics0.maxAbun[i];

   /* setting up the intial guesses for the abundances of base species */
   initConstantsAndAbundances( );
   /* setting the maximum allowed absolute abundances */
   setMaxAllowedDensities();

   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int remove_mesh( int id )
{
   freeMesh( meshes[id], chemNet.nSpec );
   nMeshes--;
   flags[id]=0;

   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int isMyMesh( int id )
{
   if( id % nProcs == rank )
      return 1;
   else
      return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int get_errorFlags2(const int *ids, const int n, int *errFlagsRet)
{
   int i;
   for(i=0;i<n;i++) {
      printf("%d\n", ids[i]);
      errFlagsRet[i]=-i;
   }

   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int set_outputDir( char *dirName )
{
   sprintf(pv.outputDir, "%s", dirName);
   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int set_species_fName( char *fName )
{
   sprintf(pv.species_fName, "%s", fName);
   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int set_underUbundant_fName( char *fName )
{
   sprintf(pv.underUbundant_fName, "%s", fName);
   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int set_rate99_fName( char *fName )
{
   sprintf(pv.rate99_fName, "%s", fName);
   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int set_selfSheilding_CO_fName( char *fName )
{
   sprintf(pv.selfSheilding_CO_fName, "%s", fName);
   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int set_rotationalCooling_baseName( char *baseName )
{
   sprintf(pv.rotationalCooling_baseName, "%s", baseName );
   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int set_vibrationalCooling_baseName( char *baseName )
{
   sprintf(pv.vibrationalCooling_baseName, "%s", baseName );
   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int set_zeta( double zeta )
{
   pv.zeta = zeta;
   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int set_S_depletion( double S_depletion )
{
   pv.S_depletion = S_depletion;
   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int set_TTol( double TTol )
{
   pv.TTol = TTol;
   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int set_CTol( double CTol )
{
   pv.CTol = CTol;
   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int set_metalicity( double metalicity )
{
   pv.metalicity = metalicity;
   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int set_AvMax( double AvMax )
{
   pv.AvMax = AvMax;
   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int set_slabSizeCrit( double slabSizeCrit )
{
   pv.slabSizeCrit = slabSizeCrit;
   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int set_max_deltaAv( double max_deltaAv )
{
   pv.max_deltaAv = max_deltaAv;
   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int set_min_deltaAv( double min_deltaAv )
{
   pv.min_deltaAv = min_deltaAv;
   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int set_maxSlabs( int maxSlabs )
{
   pv.maxSlabs = maxSlabs;
   return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int set_database_fName( char *fName )
{
   sprintf(pv.database_fName, "%s", fName);
   return 0;
}
