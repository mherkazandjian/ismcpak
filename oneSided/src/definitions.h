#ifndef __DEFINITIONS_H
#define __DEFINITIONS_H

/******************** Service macros ********************************/
#define max(x,y)     ( (x)>(y) ? (x) : (y) )
#define min(x,y)     ( (x)<(y) ? (x) : (y) )
#define dabs(x)      ( (x) < 0 ? (-x) : (x) )
#define sign(x)      ( (x) < 0 ? -1.0 : 1.0 )
#define SQ(x)        ( (x)*(x) )

/******************** Resource Strings *******************************/
#define ALLDONE    "All done"
#define FAILOPEN   "I cannot open this file"
#define PROCESSING "Working on file"
#define WRITERROR  "Write error"
/******************** Memory macros **********************************/
#define MAKE_1D_ARRAY(name, elements, type)	\
		{name=(type *)calloc(elements, sizeof(type));\
		if (name==(type *)NULL) errMsg(15);\
		heapUse+=(elements)*sizeof(type);}
#define FREE_1D_ARRAY(name, elements, type)	\
		{free(name);\
		heapUse-=(elements)*sizeof(type);}

/* creates a 2D array A of m rows and n colomns  */
/* the entries of the array can be accessed with */
/* A[i][j], where i is the  row index and j is   */
/* the colon index. Like MATLAB notation         */
/* example :                                     */
/* we create an array of 5 rows and 10 colons,   */
/*                                               */
/*    MAKE_2D_ARRAY(A, 5, 10, double)            */
/*                                               */
/* to set all the elements of the 3rd row to -1  */
/*                                               */
/*   for(i=0;i<10;i++)                           */
/*       A[2][i] = -1                            */
/* the memory allocation scheme is done such     */
/* colons are allocated as contigueous sections  */
/* in memeory. So to write the k^th colon, it    */
/* can be done thisway :                         */
/* fwrite(&(A[k][0]), sizeof(double), n, fd)     */

#define	MAKE_2D_ARRAY(rowPtr, m, n, type) \
		{register type *dataPtr;\
		int i;\
		dataPtr=(type *)calloc((m)*(n),sizeof(type));	\
		if (dataPtr==(type *)NULL) errMsg(15);\
		rowPtr=(type **)calloc(m, sizeof (type *));\
		if (rowPtr==(type **)NULL) errMsg(15);\
		for (i=0; i<m; i++) {\
		  rowPtr[i]=dataPtr;\
		  dataPtr+=(n);}\
		heapUse+=(m)*(n)*sizeof(type);}
#define	FREE_2D_ARRAY(rowPtr, m, n, type) \
		{free(*rowPtr);\
		free(rowPtr);\
		heapUse-=(m)*(n)*sizeof(type);}

#define	MAKE_3D_ARRAY(gridPtr, grid, rows, columns, type) \
		{register type *data_Ptr; register type **rowPtr;\
		int row_count;\
		data_Ptr=(type *)calloc((grid)*(rows)*(columns), sizeof(type));\
		if (data_Ptr==(type *)NULL) errMsg(15);\
		rowPtr=(type **)calloc((grid)*(rows), sizeof (type *));\
		if (rowPtr==(type **)NULL) errMsg(15);\
		gridPtr=(type ***)calloc(grid, sizeof (type **));\
		if (gridPtr==(type ***)NULL) errMsg(15);\
		for (row_count=0; row_count<(grid)*(rows); row_count++) {\
		  rowPtr[row_count]=data_Ptr;\
		  data_Ptr+=(columns);}\
		for (row_count=0; row_count<grid; row_count++) {\
		  gridPtr[row_count]=rowPtr;\
		  rowPtr+=(rows);}\
		heapUse+=(grid)*(rows)*(columns)*sizeof(type);}
#define	FREE_3D_ARRAY(gridPtr, grid, rows, columns, type) \
		{free(**gridPtr);\
		free(*gridPtr);\
		free(gridPtr);\
		heapUse-=(grid)*(rows)*(columns)*sizeof(type);}
#define	MAKE_4D_ARRAY(blockPtr, blocks, grid, rows, columns, type) \
		{register type *data_Ptr; register type **rowPtr; register type ***gridPtr;\
		int row_count;\
		data_Ptr=(type *)calloc((blocks)*(grid)*(rows)*(columns), sizeof(type));\
		if (data_Ptr==(type *)NULL) errMsg(15);\
		rowPtr=(type **)calloc((blocks)*(grid)*(rows), sizeof (type *));\
		if (rowPtr==(type **)NULL) errMsg(15);\
		gridPtr=(type ***)calloc((blocks)*(grid), sizeof (type **));\
		if (gridPtr==(type ***)NULL) errMsg(15);\
                blockPtr=(type ****)calloc(blocks, sizeof (type ***));\
		if (blockPtr==(type ****)NULL) errMsg(15);\
		for (row_count=0; row_count<(blocks)*(grid)*(rows); row_count++) {\
		  rowPtr[row_count]=data_Ptr;\
		  data_Ptr+=(columns);}\
		for (row_count=0; row_count<(blocks)*(grid); row_count++) {\
		  gridPtr[row_count]=rowPtr;\
		  rowPtr+=(rows);}\
                for (row_count=0; row_count< blocks; row_count++) {\
		  blockPtr[row_count]=gridPtr;\
		  gridPtr+=(grid);}\
		heapUse+=(blocks)*(grid)*(rows)*(columns)*sizeof(type);}
#define	FREE_4D_ARRAY(blockPtr, blocks, grid, rows, columns, type) \
		{free(***blockPtr);\
                free(**blockPtr);\
		free(*blockPtr);\
		free(blockPtr);\
		heapUse-=(blocks)*(grid)*(rows)*(columns)*sizeof(type);}

#define	READ_SPACE(nChars, file)\
                  fgets(__SPACE_BUFFER__, nChars, file);

/******************** Control constants *****************************/
enum errorcodes
  {NO_ERROR, DIALOGUE_ERROR, FILE_OPEN_ERROR, OUTPUT_ERROR, 
   EOF_ERROR, FPRINTF_ERROR, WRONG_TYPE_ERROR, MEMORY_ERROR,
   SINGULAR_MATRIX};

enum chemical_reactions
  {REACT, PHOTION, IONMOL, ELREC, NEUTRAL, CHARGEX, RADASS, ACTIBAR, CRPHOTON, REACTION};

enum finestruct
  {CII, CI, OI, SPECIES};

enum   coollines
  {H2, H2_low_T, H2O, H2O_ortho, H2O_para, CO, CO_low_T, COOLSPEC};

enum self_shld
  {_12CO, _13CO, _C_18O, _13C_18O, MOL};

/* for non void functions */
#define CHECK_ERR if(nErrors>0){\
                       fprintf(stderr, "caught error in file \"%s\" function \"%s\" line %d... going up one level\n",\
                                         __FILE__, __FUNCTION__,__LINE__);\
		       return -errCode;					\
                       }

/* for void functions */
#define CHECK_ERRV if(nErrors>0){\
                       fprintf(stderr, "caught error in file \"%s\" function \"%s\" line %d... going up one level\n",\
                           __FILE__, __FUNCTION__,__LINE__);\
		       return;\
                    }


#endif
