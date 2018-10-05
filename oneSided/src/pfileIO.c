#include "vars.h"
#include "prototypes.h"

#define MAX_LINE_LENGTH  1024
#define MAX_TOKEN_LENGTH 256

const char DELIMITERS1[] = "\n\0";         /* new line delimiters */
const char DELIMITERS2[] = " \f\n\r\t\v";  /* isspace characters  */
const char COMMENT       = '#';            /* comment symbol      */

/*----------- headers ----------------*/
char * copyFileToBuffer(const char *fName, int *chrBuffLen);
int    getParamAndValue(char *line, char *pname, char *pvalue);
void   setParamValueToStruct(const char *pname, const char *pvalue,  parameters *parms, short int *readFlag);
int    dumpParametersToFile(const parameters *parms);
/*------------------------------------*/

int readPfile(char *fname, parameters *parms)
{
    char *buffer, *line, pname[MAX_TOKEN_LENGTH] , pvalue[MAX_TOKEN_LENGTH];
    int  buffLen, current, shift, status;
    short int *readFlag, i;
    
    /*allocate and init array which holds tag found flags*/
    readFlag=malloc(sizeof(short int)*NTAGS);
    for(i=0;i<NTAGS;i++) readFlag[i]=0;
    
    /*copying the parameter file into a local buffer*/
    buffer = copyFileToBuffer(fname, &buffLen);    
    
    /*parsing the parameter file*/
    shift = current = 0;
    
    while( (current += shift) < buffLen-1 ) {
	
	/* get a line */
	buffer  += shift;
	line = strtok(buffer, DELIMITERS1); 
	
	if(line == NULL) { 
	    shift = 1;                /* empty line, move to next one ..    */
	    continue;
	} else 
	    shift = strlen(line)+1;
	
	/* parse each line */
	status = getParamAndValue(line, pname, pvalue);  
	if( status == 0 ) continue;
	
	/* assigninig the read parameter value to the correcponding value in the parm struct */
	setParamValueToStruct(pname, pvalue, parms, readFlag);
    }
    
    /*checking for missing parameters read values*/
    for(i=0;i<NTAGS;i++) {
	if(readFlag[i]==0) {
	    fprintf(stderr, "value of \"%s\" not set : ",tagNames[i]);
	    errMsg(6);
	}
    }
    
    dumpParametersToFile(parms);   /*dumping the read paramters into a file*/
    
    return TRUE;
}
/*--------------------------------------------------------------------------------------*/
/* copy the file (fName) into a memory buffer (chrBuffer) and return the number of      */
/* characters in the file                                                               */
char * copyFileToBuffer(const char *fName, int *chrBuffLen)
{
    FILE *fp;
    int nc1=0, nc2=0;
    char c, *chrBuffer;

    /* counting the number of charachters and allocating the char buffer*/
    if((fp=fopen(fName,"r"))) {

	while((c=getc(fp) )!= EOF) 
	    nc1++;
	
	chrBuffer=(char *)malloc(sizeof(char)*(nc1+1));
	fflush(fp);
	fclose(fp);
    } else {
	fprintf(stderr, "%s:%d : file %s not found\n", __FILE__, __LINE__, fName);
	return NULL;
    }

    /* copying the file to the char buffer */
    fp=fopen(fName,"r");
    while((c=getc(fp) )!= EOF) {
	chrBuffer[nc2]=c;
	nc2++;
    }
    fflush(fp);
    fclose(fp);
    
    chrBuffer[nc2+1]='\0'; /* appending with the null char */
    
    if(nc2!=nc1) 
	return FALSE;
    else {
	*chrBuffLen = nc1;
	return chrBuffer;
    }
}
/*------------------------------------------------------------------------------------------*/
void setParamValueToStruct(const char *pname, const char *pvalue, parameters *parms, short int *readFlag)
{
    int i, found=0;
	
    for(i=0;i<NTAGS;i++) {
	
	if(strcmp(tagNames[i], pname)==0) {
	    
	    found++;
	    switch(i) {
               /*add your parameters here*/
		case 0  : parms->dens0=atof(pvalue);                          break;
		case 1  : parms->G0=atof(pvalue);                             break;
		case 2  : parms->gamma_mech=atof(pvalue);                     break;
		case 3  : parms->T_init=atof(pvalue);                         break;
		case 4  : strcpy(parms->outputDir, pvalue);                   break;
		case 5  : strcpy(parms->species_fName, pvalue);               break;
		case 6  : strcpy(parms->underUbundant_fName, pvalue);         break;
		case 7  : strcpy(parms->rate99_fName, pvalue);                break;
		case 8  : strcpy(parms->selfSheilding_CO_fName, pvalue);      break;
		case 9  : strcpy(parms->rotationalCooling_baseName, pvalue);  break;
		case 10 : strcpy(parms->vibrationalCooling_baseName, pvalue); break;
		case 11 : parms->zeta=atof(pvalue);                           break;
		case 12 : parms->S_depletion=atof(pvalue);                    break;
	        case 13 : parms->TTol=atof(pvalue);                           break;
		case 14 : parms->CTol=atof(pvalue);                           break;
		case 15 : parms->metalicity=atof(pvalue);                     break;
		case 16 : parms->AvMax=atof(pvalue);                          break;
		case 17 : parms->slabSizeCrit=atof(pvalue);                   break;
		case 18 : parms->max_deltaAv=atof(pvalue);                    break;
		case 19 : parms->min_deltaAv=atof(pvalue);                    break;
	        case 20 : parms->maxSlabs=atoi(pvalue);                       break;		  
	        case 21 : strcpy(parms->database_fName, pvalue);              break;		  
	    }
	    readFlag[i]++;
	}
    }

    if( found==0 ) {
	fprintf(stderr, "%s:%d: parameter |%s| not found\n", __FILE__, __LINE__, pname);
	exit(-1);
    }
}
/*------------------------------------------------------------------------------------------*/
int dumpParametersToFile(const parameters *parms)
{
    int i=0;
    char *outFname;
    FILE *fp;

    /* if parms->outputDir is not define, just name the file parmaeter */
    /* and dump it in the current working dir                          */
    outFname=malloc(strlen(parms->outputDir)+strlen("paramters\0")); 
    sprintf(outFname,"%s%s%c",parms->outputDir,"parameters",'\0');

    if((fp=fopen(outFname,"w"))) {
	
	fprintf(fp,"%-20s %e\n", tagNames[i++], parms->dens0);
	fprintf(fp,"%-20s %e\n", tagNames[i++], parms->G0);
	fprintf(fp,"%-20s %e\n", tagNames[i++], parms->gamma_mech);
	fprintf(fp,"%-20s %e\n", tagNames[i++], parms->T_init);
	fprintf(fp,"%-20s %s\n", tagNames[i++], parms->outputDir);
	fprintf(fp,"%-20s %s\n", tagNames[i++], parms->species_fName);
	fprintf(fp,"%-20s %s\n", tagNames[i++], parms->underUbundant_fName);
	fprintf(fp,"%-20s %s\n", tagNames[i++], parms->rate99_fName);
	fprintf(fp,"%-20s %s\n", tagNames[i++], parms->selfSheilding_CO_fName);
	fprintf(fp,"%-20s %s\n", tagNames[i++], parms->rotationalCooling_baseName);
	fprintf(fp,"%-20s %s\n", tagNames[i++], parms->vibrationalCooling_baseName);
	fprintf(fp,"%-20s %e\n", tagNames[i++], parms->zeta);
	fprintf(fp,"%-20s %e\n", tagNames[i++], parms->S_depletion);
	fprintf(fp,"%-20s %e\n", tagNames[i++], parms->TTol);
	fprintf(fp,"%-20s %e\n", tagNames[i++], parms->CTol);
	fprintf(fp,"%-20s %e\n", tagNames[i++], parms->metalicity);
	fprintf(fp,"%-20s %e\n", tagNames[i++], parms->AvMax);
	fprintf(fp,"%-20s %e\n", tagNames[i++], parms->slabSizeCrit);
	fprintf(fp,"%-20s %e\n", tagNames[i++], parms->max_deltaAv);
	fprintf(fp,"%-20s %e\n", tagNames[i++], parms->min_deltaAv);
	fprintf(fp,"%-20s %d\n", tagNames[i++], parms->maxSlabs);
	fprintf(fp,"%-20s %s\n", tagNames[i++], parms->database_fName);

	fflush(fp);
	fclose(fp);
	free(outFname);
	return TRUE;
    } else {
	fprintf(stderr, "%s:%d : file %s not found\n", __FILE__, __LINE__, outFname);
	return FALSE;
    }
}
/*------------------------------------------------------------------------------------------*/
/* fetch the parameter/value pair, retuns 1 if found and zero if no pair is found           */
int getParamAndValue(char *line, char *pname, char *pvalue) 
{
    char *word, lineBckp[MAX_LINE_LENGTH];

    /* make a copy of the original line   */
    strcpy(lineBckp, line);           
    
    /* extract first token from the line (either a params or a commented line) */
    word = strtok(line, DELIMITERS2); 
    if( word == NULL ) {
	return 0;                          /* blank line, move to next line */
    } else {
	if(word[0]==COMMENT) {
	    return 0;                      /*commented  line, move to the next */
	} else {
	    if(strchr(word, COMMENT)!=NULL) {  /* parameter should not contain comment symbol */
		fprintf(stderr, "%s:%d: parameter name cannot contain comment sybmol :\n    %s\n",
			__FILE__, __LINE__, lineBckp);
		exit(-1);
	    }
	}
    }
    strcpy(pname, word);

    /*parameter velue expected*/
    word = strtok(NULL, DELIMITERS2);
    if( word == NULL ) {
	fprintf(stderr, "%s:%d: missing parameter value in parameter file line :\n    %s\n",
		__FILE__, __LINE__, lineBckp);
	exit(-1);
    } else {
	if( strchr(word, COMMENT) != NULL ) {
	    fprintf(stderr, "%s:%d: comment symbol found where parameter value expected:\n    %s\n",
		    __FILE__, __LINE__, lineBckp);
	    exit(-1);
	}
    }
    strcpy( pvalue, word);
    
    /* getting the commented part */
    word = strtok(NULL, DELIMITERS2);  
    if(word!=NULL) {
	if(word[0]!=COMMENT) {
	  fprintf(stderr, "%s:%d: comment part doesnt start with comment symbol:\n    %s\n",
		  __FILE__, __LINE__, lineBckp);
	  exit(-1);
	}
    }
	
    return 1;
}
