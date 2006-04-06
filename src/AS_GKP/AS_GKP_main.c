
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received (LICENSE.txt) a copy of the GNU General Public 
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/
static char CM_ID[] = "$Id: AS_GKP_main.c,v 1.7 2006-04-06 18:53:14 brianwalenz Exp $";

/*************************************************
* Module:  AS_GKP_main.c
* Description:
*    Gatekeeper main
* 
*    Reference: GateKeeper.rtf
*
*    Command Line Interface:
*        gatekeeper [-G] [-O][-af] -e <errorThreshhold> <gateKeeperStorePath> <inputFile>.<ext>
*
*        -a append to the gateKeeperStore (it it doesn't exist, create it )
*        -f force the creation of a new store 
*        -e <errorThreshhold> Set the error threshold to the value (default = 1)
*        -O ignore messages that are intended only for Grande (shredded BACs/Bactigs)
*        -G ignore messages that are intended only for overlay (FULLBacs and BACtigs)
* 
*        gatekeeperStorePath   Used to find/create a directory that will house 4 files:
*           gkp.phash  Symbol table for UID mapping and to keep track of RPT_IDs
*                  Each UIDand RPT_ID has a reference count.  PHash also keeps track
*                  of IIDs that have been assigned for each class of input.
*                  See UTL_PHash.[ch] for a more complete description.
*           gkp.frg A GateKeeperFragmentStore.  An array of GateKeeperFragmentRecords (see AS_PER_GkpStore.h)
*                  One record per fragment.  Stores IID->UID map, and links to gkplstore for
*                  links.
*           gkp.lnk A GateKeeperLinkStore.  An array of GateKeeperLinkRecords (see AS_PER_GkpStore.h)
*                  One record per LKG record, and one per set of Join messages.  This is
*                  maintained as an array of linked records.
*           gkp.btg A GateKeeperBacTigStore.  An array of GateKeeperBactigRecords (see AS_PER_GkpStore.h)
*                  One record per BTG subrecord.
*           gkp.loc A GateKeeperLocaleStore.  An array of GateKeeperLocaleRecords (see AS_PER_GkpStore.h)
*                  One record per BAC record.
*           gkp.s_loc A GateKeeperLocaleStore.  An array of GateKeeperLocaleRecords (see AS_PER_GkpStore.h)
*                  One record for each redefinition of a BAC.
*           gkp.bat A GateKeeperBatchStore.  An array of GateKeeperBatchRecords (see AS_PER_GkpStore.h)
*                  One record per BAT record.
*           gkp.dst A GateKeeperDistanceStore.  An array of GateKeeperDistanceRecords (see AS_PER_GkpStore.h)
*                  One record per DST record.
*           gkp.s_dst A GateKeeperDistanceStore.  An array of GateKeeperDistanceRecords (see AS_PER_GkpStore.h)
*                  One record for each redefinition of a DST.
*           gkp.scn A GateKeeperScreenStore.  An array of GateKeeperScreenRecords (see AS_PER_GkpStore.h)
*                  One record per SCN record.
*           gkp.rpt A GateKeeperRepeatStore.  An array of GateKeeperRepeatRecords (see AS_PER_GkpStore.h)
*                  One record per RPT record.
*           gkp.seq A GateKeeperSequenceStore.  An array of GateKeeperSequenceRecords (see AS_PER_GkpStore.h)
*                  One record per SEQ record.
*           gkp.pla A GateKeeperPlateStore.  An array of GateKeeperPlateRecords (see AS_PER_GkpStore.h)
*                  One record per PLA record.
*           gkp.wel A GateKeeperWellStore.  An array of GateKeeperWellRecords (see AS_PER_GkpStore.h)
*                  One record per WEL record.  These are NOT indexed, just collected.  They unltimatedly
*                  will be associated with all of the Celera reads.

*       gatekeeper processes the input file, reporting errors and warnings.  Any message that causes
*       an error warning is output to stderr, along with the associated warning/error messages.
*       Messages that are ignored due to not being appropriate are sent to <inputFile>.ign
*       If more than <errorThreshhold> errors are detected (default 10), no output is produced.
*       If less than <errorThreshhold> errors are detected, the output data set will be found
*       in <inputFile>.inp, <inputFile.ign> and the gateKeeperStore will be updated.
*
*       See GateKeeper.rtf or AS_GKP_include.h for descriptions of the checks performed.
*
*    Programmer:  S. Kravitz
*       Written:  Jan 1999
*       Revised:  Feb 2000
* 
*************************************************/

//#define DEBUG 1
//#define DEBUGIO 1

/******************************************************************************/
/* AS_GKP_main
 *    Main for gatekeeper.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_PHash.h"
#include "AS_UTL_version.h"
#include "AS_MSG_pmesg.h"
#include "AS_GKP_include.h"
#include "AS_UTL_param_proc.h"

int  nerrs = 0;   // Number of errors in current run
int maxerrs = -1; // Number of errors allowed before we punt
static void usage(void);

MesgReader Reader;
MesgWriter Writer, ErrorWriter;
FILE *Infp, *Outfp, *Msgfp, *Ignfp, *Ibcfp, *Errfp;
GateKeeperStore GkpStore;
GateKeeperStore GkpStore_input;
char *GlobalParamText  = NULL;
/********************************************/
SequenceBucketArrayT *LinkerDetector_READ;   // Used to collect stats for 1-mers - 8-mers on 5-prime ends of sequence
SequenceBucketArrayT *Linker3pDetector_READ;   // Used to collect stats for 1-mers - 8-mers on 3-prime ends of sequence
SequenceBucketArrayT *LinkerDetector_EBAC;   // Used to collect stats for 1-mers - 8-mers on 5-prime ends of sequence
SequenceBucketArrayT *LinkerDetector_LBAC;   // Used to collect stats for 1-mers - 8-mers on 5-prime ends of sequence

SequenceBucketArrayT *SanityDetector_READ;   // Used to collect stats for 1-mers - 8-mers on clr_rng_end -50bp of sequence
SequenceBucketArrayT *SanityDetector_EBAC;   // Used to collect stats for 1-mers - 8-mers on clr_rng_end -50bp of sequence
SequenceBucketArrayT *SanityDetector_LBAC;   // Used to collect stats for 1-mers - 8-mers on clr_rng_end -50bp of sequence

SequenceBucketT *SequenceProbabilities; // Used to collect stats on all sequence within clear ranges


SequenceLengthHistogramT *Linker5pHistogram = NULL;
SequenceLengthHistogramT *Linker3pHistogram = NULL;
SequenceLengthHistogramT *LinkerSanityHistogram= NULL;


/********************************************/

  char  Output_File_Name_Prefix [FILENAME_MAX];
  char  File_Name_Prefix [FILENAME_MAX];
  char  Bac_File_Name [FILENAME_MAX];
  char  Input_Store_Name [FILENAME_MAX];
  char  Output_Store_Name [FILENAME_MAX];
  char  Ignore_File_Name [FILENAME_MAX];
  char  Error_File_Name [FILENAME_MAX];
  char  Input_File_Name [FILENAME_MAX];
  char  Output_File_Name [FILENAME_MAX];
  char  Parameter_File_Name [FILENAME_MAX];

#define NMER_SIGMA_THRESH (60.0)


int  main(int argc, char * argv [])

{
  int status = 0;
  char cmd[FILENAME_MAX * 4];
  char tmpFilePath[FILENAME_MAX];
  int strict = FALSE;
  int  illegal, create, append, force, input_exists, output_exists, check_qvs, verbose, check_nmers,
    batchNumFileNames, experimental, compatibility, inputStoreSpecified, outputStoreSpecified;
  int assembler = AS_ASSEMBLER_GRANDE;
  char *projectName = NULL;
  char *param_projectName = NULL;
  char *param_maxerrs = NULL;
  double threshhold = NMER_SIGMA_THRESH;
  char *parameterPath;
  OutputType output;
  mode_t   mode = S_IRWXU | S_IRWXG | S_IROTH;
  CDS_CID_t currentBatchID = NULLINDEX;
  VA_TYPE(int32) *addedBactigs = NULL;
  VA_TYPE(PtrT) *nmersToCheck = CreateVA_PtrT(100);
  int matchBAC_DstsInLkgMsgs = 1;

  experimental = 0;
  projectName = NULL;
  parameterPath = NULL;
  verbose = 0;
  illegal = 0;
  create = 1;
  append = 0;
  force = 0;
  input_exists = 0;
  output_exists = 0;
  check_qvs = 1;
  check_nmers = 1;
  strict = FALSE;
  compatibility = FALSE;
  batchNumFileNames = TRUE;
  output = AS_BINARY_OUTPUT;
  inputStoreSpecified = outputStoreSpecified = FALSE;





  /**************** Process Command Line Arguments *********************/
  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    optind = 1;
    while (!errflg && ((ch = getopt(argc, argv, "bfe:i:o:p:t:n:achCNGsOPQTXvd:l")) != EOF)){
      switch(ch) {
      case 'l':
        matchBAC_DstsInLkgMsgs = 0;
        break;
      case 'a':
	append = 1;
	create = 0;
	force  = 0;
	break;
      case 'b':
	batchNumFileNames = TRUE;
	break;
      case 'c':
	create = 1;
	break;
      case 'd':
	fprintf(stderr,"* option -d %s  length = " F_SIZE_T "\n", optarg, strlen(optarg));
	if(strlen(optarg) != 8 ||
	   (strspn(optarg,"actg") != 8)){
	  fprintf(stderr,"* -d option must specify 8-mer from alphabet [actg]*\n");
	  usage();
	}
	fprintf(stderr,"* Appending %s to nmersToCheck " F_SIZE_T "\n",
		optarg, GetNumVA_PtrT(nmersToCheck));

	{
	  char *copy = strdup(optarg);
	  
	  AppendVA_PtrT(nmersToCheck, (const void *)&copy);
	}
	break;
      case 'e':
	fprintf(stderr,"* -e option optarg: %s\n", optarg);
	maxerrs = atoi(optarg);
	fprintf(stderr,"* maxerrs set to %d\n", maxerrs );
	break;
      case 'f':
	force  = 1;
	append = 0;
	nerrs = 0;
	break;

      case 'h':
	usage();
	fprintf(stderr,"The following failures are detected:\n");
	printAllGKPErrors(stderr);
	exit(1);
	break;

      case 'i':
	create = 0;
	fprintf(stderr,"* -i %s\n", optarg);
	strcpy(Input_Store_Name, optarg);
	inputStoreSpecified = TRUE;
	break;
      case 'o':
	fprintf(stderr,"* -o %s\n", optarg);
	strcpy(Output_Store_Name, optarg);
	outputStoreSpecified = TRUE;
	break;
      case 'p':
	{
	  char *allParams;
	fprintf(stderr,"* -p %s\n", optarg);
	strcpy(Parameter_File_Name, optarg);

	/* Load parameters file */
	loadParams(Parameter_File_Name);
	getAllParams("gatekeeper", &allParams);
	if(GlobalParamText){
	  free(GlobalParamText);
	}
	GlobalParamText = (char *)malloc(strlen(allParams)+100);
	sprintf(GlobalParamText,"Parameters from parameter file are:\n%s\n",
		allParams);
	fprintf(stderr,"* Loaded parameters %s\n", allParams);
	
	/* Load all command line parameters  from parameters file */
	param_projectName = getParam("gatekeeper.project_name");
	if(param_projectName){
	  fprintf(stderr,"* Read Project Name: %s from parameter file\n",param_projectName);
	}
	param_maxerrs = getParam("gatekeeper.maxErrors");
	if(param_maxerrs){
	  fprintf(stderr,"* Read maxErrors: %s from parameter file\n",param_maxerrs);
	}

	}
	break;
      case 'n':
	fprintf(stderr,"* Project name specified <%s> with -n argument.\n", optarg);
	projectName = optarg;
	break;
      case 's':
	strict = TRUE;
	break;
      case 't':
	threshhold = atof(optarg);
	fprintf(stderr,"** n-mer screening threshhold set to %f sigma\n", threshhold);
	break;
      case 'v':
	verbose = 1;
	break;
      case 'C':
	compatibility = TRUE;
	break;
      case 'G':
	assembler = AS_ASSEMBLER_GRANDE;
	break;
      case 'N':
	fprintf(stderr,"** n-mer screening disabled\n");
	check_nmers = 0;
	break;
      case 'O':
	assembler = AS_ASSEMBLER_OVERLAY;
	break;
      case 'P':
	output = AS_PROTO_OUTPUT;
	break;
      case 'Q':
	check_qvs = 0;
	break;
      case 'T':
        assembler = AS_ASSEMBLER_OBT;
        break;
      case 'X':
        experimental = 1;
	break;
      case '?':
	fprintf(stderr,"Unrecognized option -%c",optopt);
      default :
	errflg++;
	illegal = 1;
      }
  }

    switch (assembler) {
      case AS_ASSEMBLER_GRANDE:  fprintf(stderr, "* Gatekeeper for Assembler Grande\n"); break;
      case AS_ASSEMBLER_OVERLAY: fprintf(stderr, "* Gatekeeper for Overlay Assembler\n"); break;
      case AS_ASSEMBLER_OBT:     fprintf(stderr, "* Gatekeeper for Assembler Grande with Overlap Based Trimming\n"); break;
    }

    if(force == 1 && append  == 1)
      {
	fprintf (stderr,
		 "* Illegal combination of command line flags"
		 "-- they are mutually exclusive\n");
	illegal = 1;
      }
     
    if((illegal == 1) || (argc - optind != (compatibility?2:1) ))
      {
	fprintf(stderr,"* argc = %d optind = %d compatibility = %d %s\n",
		argc, optind, compatibility, argv[optind]);
	usage();
      }

    if(compatibility){
      batchNumFileNames = FALSE;
      strcpy(Output_Store_Name,argv[optind++]);
    }
    if(optind < argc) 
      {
	char *suffix;

	suffix = strrchr(argv[optind],(int)'.');
       
	fprintf(stderr,"Input file is %s suffix is %s\n",argv[optind], suffix);
	strcpy(Input_File_Name, argv[optind]);
	Infp = File_Open (Input_File_Name, "r", TRUE);     // frg file
	Reader = (MesgReader)InputFileType_AS(Infp);

	if(suffix)
	  *suffix = '\0';

	strcpy(File_Name_Prefix,argv[optind]);

	Writer = (MesgWriter)OutputFileType_AS(output);
	ErrorWriter = (MesgWriter)OutputFileType_AS(AS_PROTO_OUTPUT);


	optind++;
      }
    /* End of command line parsing */
  }
   

  if(projectName == NULL &&
     param_projectName){
    fprintf(stderr,"* Using project name from param file: %s\n", param_projectName);
    projectName = param_projectName;
  }
  if(maxerrs == -1){
    if(!param_maxerrs){
      maxerrs = 1;
      fprintf(stderr,"* Using default maxerrs: %d\n", maxerrs);
    }else{
      maxerrs = atol(param_maxerrs);
      fprintf(stderr,"* Using maxerrs from param file: %d\n", maxerrs);
    }
  }

  

  if(projectName == NULL &&
     batchNumFileNames == TRUE){
    fprintf(stderr,"* Project Name must be specified on command line or through parameter file..exiting\n");
    exit(1);
  }
  if(!experimental &&
     compatibility){
    fprintf(stderr,"* Compatiblity mode is only valid with -X option..exiting\n");
    exit(1);
  }
     

    if(compatibility){
      /**************** Open or Create Files *********************/
      InitGateKeeperStore(&GkpStore, Output_Store_Name);
      output_exists = TestOpenGateKeeperStore(&GkpStore);
      strcpy(Output_File_Name_Prefix,File_Name_Prefix);

      if(force && output_exists){
	fprintf(stderr,"* Gatekeeper Store %s exists ...nuking\n", Output_Store_Name);
	RemoveGateKeeperStoreFiles(&GkpStore);
	create = 1;
	append = 0;
      }

      if(append){
	if( output_exists == 0){
	  fprintf(stderr,"* Directory %s DOES NOT exist ...creating before append\n", Output_Store_Name);
	  create = 1;
	  append = 0;
	}else if(output_exists == -1){
	  fprintf(stderr,"* Directory %s DOES exist ...but not all files are present...exiting\n", 
		  Output_Store_Name);
	  exit(1);
	}else {
	  fprintf(stderr,"* Directory %s DOES  exist ...\n", Output_Store_Name);
	  append = 1;
	  create = 0;
	}
      }

      /* We need to make backup copies of the old stores, and
      restore them if we hit fatal errors */

      if(append){
	fprintf(stderr,"* Appending to gateKeeperStore %s\n", Output_Store_Name);
	sprintf(tmpFilePath,"%s___tmp", Output_Store_Name);

	/* Remove temporary files if any */
	sprintf(cmd,"rm -rf %s", tmpFilePath);
        if(system(cmd) != 0) assert(0);
	mkdir(tmpFilePath, mode);
	CopyGateKeeperStoreFiles(&GkpStore, tmpFilePath);
     
	InitGateKeeperStore(&GkpStore, tmpFilePath);
	OpenGateKeeperStore(&GkpStore);
      }else  if(create){
	char buffer[FILENAME_MAX];

	fprintf(stderr,"* output_exists = %d\n", output_exists);
	if(output_exists == 0){
	  if(mkdir(Output_Store_Name, mode)){
	    sprintf(buffer,"gateKeeper: Failure to create directory %s", Output_Store_Name);
	    perror(buffer);
	    exit(1);
	  }
	}
	fprintf (stderr, "Creating NEW Store! Output_Store_Name = %s\n", Output_Store_Name);
	CreateGateKeeperStore(&GkpStore);
      }else{
	fprintf(stderr,"** Serious error...bye\n");
	exit(1);
      }
      currentBatchID = getNumGateKeeperBatchs(GkpStore.batStore);

    }else{
      /* First check valid store modes */
      if((inputStoreSpecified && create) ||
	 (outputStoreSpecified && append) ||
	 (!outputStoreSpecified && !append) ||
	 (!inputStoreSpecified && !create) ||
	 (!inputStoreSpecified && !outputStoreSpecified ) ||
	 (create && append) ||
	 (inputStoreSpecified && outputStoreSpecified && !strcmp(Input_Store_Name, Output_Store_Name))){
	fprintf(stderr,"* Illegal I/O combination create %d append %d input %d output %d\n",
		create, append, inputStoreSpecified, outputStoreSpecified);
	exit(1);
      }
      /* Open or create input Store */
      if(create){
	InitGateKeeperStore(&GkpStore, Output_Store_Name);
	output_exists = TestOpenGateKeeperStore(&GkpStore);

	if(output_exists){
	  fprintf(stderr,"* Gatekeeper store %s exists...\n", Output_Store_Name);
	  if(!force){
	    fprintf(stderr,"* Output store exists...exiting\n");
	    exit(1);
	  }else{
	    fprintf(stderr,"* Nuking\n");
	    RemoveGateKeeperStoreFiles(&GkpStore);
	  }
	}
	if(output_exists == 0){
	  fprintf(stderr,"* Creating directory %s in cwd of %s\n", Output_Store_Name, getcwd(NULL, 256));
	  mkdir(Output_Store_Name, mode);
	}
	CreateGateKeeperStore(&GkpStore);
      }else if(append){
	InitGateKeeperStore(&GkpStore_input, Input_Store_Name);
	strcpy(Output_Store_Name, Input_Store_Name);
	fprintf(stderr,"* Appending to gateKeeperStore %s\n", Input_Store_Name);
	sprintf(tmpFilePath,"%s___tmp", Input_Store_Name);

	/* Remove temporary files if any */
	sprintf(cmd,"rm -rf %s", tmpFilePath);
        if(system(cmd) != 0) assert(0);
	fprintf(stderr,"* Creating directory %s in cwd of %s\n", tmpFilePath, getcwd(NULL, 256));
	mkdir(tmpFilePath, mode);
	CopyGateKeeperStoreFiles(&GkpStore_input, tmpFilePath);
	InitGateKeeperStore(&GkpStore, tmpFilePath);
	OpenGateKeeperStore(&GkpStore);

      }else{ // normal case, reading from input writing to output
	InitGateKeeperStore(&GkpStore, Output_Store_Name);
	output_exists = TestOpenGateKeeperStore(&GkpStore);

	if(output_exists){
	  fprintf(stderr,"* Gatekeeper store %s exists...\n", Output_Store_Name);
	  if(!force){
	    fprintf(stderr,"* Output store exists...exiting\n");
	    exit(1);
	  }else{
	    fprintf(stderr,"* Nuking\n");
	    RemoveGateKeeperStoreFiles(&GkpStore);
	  }
	}else{

	  mkdir(Output_Store_Name, mode);

	}
	fprintf(stderr,"* TestOpen of input: %s\n", Input_Store_Name);
	InitGateKeeperStore(&GkpStore_input, Input_Store_Name);
	input_exists = TestOpenGateKeeperStore(&GkpStore_input);

	if(input_exists != 1){
	  fprintf(stderr,"* Can't open input store...exiting\n");
	  exit(1);
	}
	CopyGateKeeperStoreFiles(&GkpStore_input,Output_Store_Name);
	OpenGateKeeperStore(&GkpStore);

      }
      currentBatchID = getNumGateKeeperBatchs(GkpStore.batStore);
      fprintf(stderr,"* Current batch is " F_CID "\n", currentBatchID);
      
      sprintf(Output_File_Name_Prefix,"%s_%05" F_CIDP,projectName,currentBatchID + 1);
	

    }   
  /**************** Open Files ***********************/
  
  {


  sprintf(Output_File_Name,"%s.inp",Output_File_Name_Prefix);

  fprintf(stderr,"Output file is %s \n",Output_File_Name);
  Outfp = fopen (Output_File_Name, "w");     // inp file

  sprintf(Ignore_File_Name,"%s.ign",Output_File_Name_Prefix);
  fprintf(stderr,"Ignore file is %s \n",Ignore_File_Name);
  Ignfp = fopen (Ignore_File_Name, "w");     // ign file

  sprintf(Error_File_Name,"%s.err",Output_File_Name_Prefix);
  fprintf(stderr,"Error file is %s \n",Error_File_Name);
  Errfp = fopen (Error_File_Name, "w");     // ign file

  }
  Msgfp = Errfp;

  if(assembler == AS_ASSEMBLER_OVERLAY){
    sprintf(Bac_File_Name,"%s.ibc",Output_File_Name_Prefix);
    fprintf(stderr,"ibc file is %s \n",Bac_File_Name);
    Ibcfp = fopen(Bac_File_Name, "w");     // .ibc
  }else{
    Ibcfp = NULL;
  }
    addedBactigs = CreateVA_int32(100);

  if(Outfp && Ignfp && Errfp){
    status = GATEKEEPER_SUCCESS;
  }else{
    fprintf(stderr,"* Failed to open output files! Exiting...\n");
    status = GATEKEEPER_FAILURE;
  }	 

  if(status == GATEKEEPER_SUCCESS){
    char buffer[2048];

   /**************** Set up Linker Detectors *************
    * Each class of read-like data has its own detector  *
    * so that a multitude of data in one class won't     *
    * mask the problems in a different class             *
    * These are used in checkFrag                        *
    ******************************************************/
    LinkerDetector_READ = CreateSequenceBucketArray(8);
    Linker3pDetector_READ = CreateSequenceBucketArray(8);
    LinkerDetector_LBAC = CreateSequenceBucketArray(8);
    LinkerDetector_EBAC = CreateSequenceBucketArray(8);
    SanityDetector_READ = CreateSequenceBucketArray(8);
    SanityDetector_LBAC = CreateSequenceBucketArray(8);
    SanityDetector_EBAC = CreateSequenceBucketArray(8);

    sprintf(buffer,"%s_%s", Output_File_Name_Prefix, "5p");
    Linker5pHistogram = CreateSequenceLengthHistogram(8,buffer);
    sprintf(buffer,"%s_%s", Output_File_Name_Prefix, "3p");
    Linker3pHistogram = CreateSequenceLengthHistogram(8,buffer);
    sprintf(buffer,"%s_%s", Output_File_Name_Prefix, "Sanity");
    LinkerSanityHistogram = CreateSequenceLengthHistogram(8,buffer);


    {
      size_t i;

      for(i = 0; i < GetNumVA_PtrT(nmersToCheck); i++){
	char *nmer = *(char **)GetVA_PtrT(nmersToCheck, i);
	fprintf(stderr,"* ActivatingSequenceLengthHistogram for nmer %s\n", nmer);
	ActivateSequenceLengthHistogram(Linker5pHistogram,nmer);
	ActivateSequenceLengthHistogram(Linker3pHistogram,nmer);
	ActivateSequenceLengthHistogram(LinkerSanityHistogram,nmer);
      }
    }


    /**************** CreateTable  *********************/

    InitQualityToFractionError();



    /* This is used to collect the global probabilities of each
      base in the input set.  Used as input to check the sanity
      of the data collected in the bucket arrays */
    SequenceProbabilities = CreateSequenceBucket(1);

    /**************** Process Input  *********************/

    status = ReadFile(check_qvs,
		      check_nmers,
		      addedBactigs,
		      currentBatchID,
		      assembler,
		      strict,
		      argc, argv,
		      verbose,
                      matchBAC_DstsInLkgMsgs);


    if(check_nmers){
      /* now check that the n-mer probablities we've collected are OK */
      CheckNmerProbabilities(Msgfp, threshhold);
    }
    /* If everything looks OK, check that the input in this batch is self contained
      with respect to the overlay assembler.  OA expects that a Bactig and its associated
      fragment are in the same batch */
    if(status == GATEKEEPER_SUCCESS &&
       assembler == AS_ASSEMBLER_OVERLAY){
      int bactigStatus = CheckOverlayInputComplete(addedBactigs, verbose);
     
      if(bactigStatus != GATEKEEPER_SUCCESS){
	status = GATEKEEPER_FAILURE;
      }
       
    }

  }

   /**************** Close files   *********************/
   if(status == GATEKEEPER_SUCCESS){ //  OK to update persistent data and generate output
     /* Remove temporary files if any */
     fprintf(stderr,"#  Successful run with %d errors < %d maxerrs ..removing temp backup files\n",
	     nerrs , maxerrs);
     if(append){
       CloseGateKeeperStore(&GkpStore);
       sprintf(cmd,"mv  %s/* %s", tmpFilePath, Output_Store_Name);
       fprintf(stderr,"* %s\n", cmd);
       if(system(cmd) != 0) assert(0);
       sprintf(cmd,"rm -rf %s ", tmpFilePath);
       fprintf(stderr,"* Removing temp output store: %s\n",cmd);
       if(system(cmd) != 0) assert(0);
     }else{
       CloseGateKeeperStore(&GkpStore);
     }
   }else{
     fprintf(stderr, "# ReadFile() failed - see %s for details\n",
	     Error_File_Name);
     fprintf(stderr,"# Too Many Errors -- removing output and exiting "
	     "output_exists = %d\n", output_exists);

     if (unlink(Output_File_Name) < 0) {
       fprintf(stderr, "%s:%d - unlink(%s) failed: %s\n",
	       __FILE__, __LINE__, Output_File_Name,
	       strerror(errno));
       exit(1);
     }

     if (unlink(Ignore_File_Name) < 0) {
       fprintf(stderr, "%s:%d - unlink(%s) failed: %s\n",
	       __FILE__, __LINE__, Ignore_File_Name,
	       strerror(errno));
       exit(1);
     }

     if(append){
       sprintf(cmd,"rm -rf %s ", tmpFilePath);
       fprintf(stderr,"* Removing temp output store: %s\n",cmd);
       if(system(cmd) != 0) assert(0);
     }else{
       CloseGateKeeperStore(&GkpStore);
       RemoveGateKeeperStoreFiles(&GkpStore);
       if ((output_exists == 0)) {
         fprintf(stderr,"* Removing output store: %s ", Output_Store_Name);
         if (rmdir(Output_Store_Name) < 0) {
           fprintf(stderr, "%s:%d - rmdir(%s) failed: %s\n",
              __FILE__, __LINE__, Output_Store_Name, strerror(errno));
           exit(1);
         }
       }
     }
   }
   fclose (Errfp);
   fclose (Ignfp);
   fclose (Infp);
   fclose (Outfp);

   exit(status != GATEKEEPER_SUCCESS);
}


/******************************************************************************/



int incrementErrors(int num, FILE *msgFile){
  //  fprintf(msgFile,"* incrementErrors by %d -- (%d,%d)\n",
  //  num, nerrs, maxerrs);
  nerrs += num;
  if(nerrs >= maxerrs){
    fprintf(msgFile, "GateKeeper: max allowed errors reached %d > %d...bye\n",
	    nerrs, maxerrs);
    return(GATEKEEPER_FAILURE);
  }
  //fprintf(stderr,"* incrementErrors returning SUCCESS\n");
  return(GATEKEEPER_SUCCESS);
}

/********************************************************************************/
/********************************************************************************/
/* function ReadFile:
   Read the input stream, and invoke the appropriate check routines.

 */
int ReadFile(int check_qvs,
	     int check_nmers,
	     VA_TYPE(int32) *addedBactigs,
	     CDS_CID_t currentBatchID,
	     int32 assembler, // grande or overlay
	     int32 strict, 
	     int argc, char **argv,
	     int32 verbose,
             int matchBAC_DstsInLkgMsgs){


  MessageType imesgtype;
  GenericMesg   *pmesg;
  time_t currentTime = time(0);
  int messageCount = 0;
  int shreddedFragmentsIgnored = 0; // For overlay assembler
  int bacFragmentsIgnored = 0;      // For grande

  /* Read maxFrags fragments from the stream, adding their accession numbers, and
     those of the DST records to their respective hash tables */

  while(  EOF != Reader(Infp, &pmesg)) {

    messageCount++;
    imesgtype = pmesg->t;

    if(messageCount == 1 && imesgtype != MESG_BAT){
      fprintf(Msgfp,"# First message must be {BAT\n");
      ErrorWriter(Msgfp,pmesg);
      return GATEKEEPER_FAILURE;
    }

    switch(imesgtype){

    case MESG_BAT:
      {
	InternalBatchMesg iba_mesg;
	BatchMesg *bat_mesg = (BatchMesg *)pmesg->m;
	int gkp_result;

	

	if(messageCount != 1){
	  fprintf(Msgfp,"\n\n# Line %d of input (message %d)\n", GetProtoLineNum_AS(), messageCount);
	  printGKPError(Msgfp, GKPError_FirstMessageBAT);
	  ErrorWriter(Msgfp,pmesg);
	  return GATEKEEPER_FAILURE;
	}

	gkp_result = Check_BatchMesg(bat_mesg, &iba_mesg, currentTime, verbose);

	switch(gkp_result){

	case GATEKEEPER_SUCCESS:
	  pmesg->t = MESG_IBA;
	  pmesg->m = &iba_mesg;
	  Writer((assembler == AS_ASSEMBLER_OVERLAY?Ignfp:Outfp),pmesg);      
	  currentBatchID = iba_mesg.iaccession;
	  fprintf(stderr,"Gatekeeper reading batch " F_IID "\n",
		  iba_mesg.iaccession);
	  break;
	case GATEKEEPER_WARNING:
	case GATEKEEPER_FAILURE:
	  fprintf(Msgfp,"# Invalid BAT message at Line %d of input...exiting\n", GetProtoLineNum_AS());
	  ErrorWriter(Msgfp,pmesg);      
	  return GATEKEEPER_FAILURE;
	  break;
	default:
	  assert(0);
	}

      }
      break;
    case MESG_RPT:
      {
	RepeatItemMesg *rpt_mesg = (RepeatItemMesg *)pmesg->m;
	InternalRepeatItemMesg irp_mesg;
	int gkp_result;

	gkp_result = Check_RepeatItemMesg(rpt_mesg, &irp_mesg, currentBatchID,  verbose);

	switch(gkp_result){
	case GATEKEEPER_WARNING:
	  fprintf(Msgfp,"# Line %d of input\n", GetProtoLineNum_AS());
	  ErrorWriter(Msgfp,pmesg);      
	case GATEKEEPER_SUCCESS:
	  pmesg->t = MESG_IRP;
          pmesg->m = &irp_mesg;
	  Writer(Outfp, pmesg);
	  break;
	case GATEKEEPER_FAILURE:
	  fprintf(Msgfp,"# Line %d of input\n", GetProtoLineNum_AS());
	  ErrorWriter(Msgfp,pmesg);      
	  if(incrementErrors(1, Msgfp) == GATEKEEPER_FAILURE){
	    return GATEKEEPER_FAILURE;
	  }
	  break;
	default:
	  assert(0);
	}
      }
      break;
    case MESG_SCN:
      {
	ScreenItemMesg *scn_mesg = (ScreenItemMesg *)pmesg->m;
	InternalScreenItemMesg isn_mesg;

#ifdef DEBUGIO
	fprintf(Msgfp,"Read SCN message with eaccession " F_UID "\n", 
		scn_mesg->eaccession);
#endif

	if(GATEKEEPER_SUCCESS == 
	   Check_ScreenItemMesg(scn_mesg, &isn_mesg, currentBatchID,  verbose)){
	   pmesg->m = &isn_mesg;
	   pmesg->t = MESG_ISN;
	  Writer(Outfp,pmesg);
	}else{
	  fprintf(Msgfp,"# Line %d of input\n", GetProtoLineNum_AS());
	   ErrorWriter(Msgfp,pmesg);      
	  if(incrementErrors(1, Msgfp) == GATEKEEPER_FAILURE){
	    return GATEKEEPER_FAILURE;
	  }
	}
      }
      break;


    case MESG_DST:
      {
	DistanceMesg  *dst_mesg = (DistanceMesg *)pmesg->m;
	InternalDistMesg idt_mesg;

#ifdef DEBUGIO
	fprintf(Msgfp,"Read DST message with accession " F_UID " (%d)\n", 
		dst_mesg->eaccession,
		dst_mesg->action);
#endif

	if(GATEKEEPER_SUCCESS == 
	   Check_DistanceMesg(dst_mesg, &idt_mesg, currentBatchID,  verbose)){
	   pmesg->m = &idt_mesg;
	   pmesg->t = MESG_IDT;
	  Writer(Outfp,pmesg);      
	}else{
	  fprintf(Msgfp,"# Line %d of input\n", GetProtoLineNum_AS());
	   ErrorWriter(Msgfp,pmesg);      
	  if(incrementErrors(1, Msgfp) == GATEKEEPER_FAILURE){
	    return GATEKEEPER_FAILURE;
	  }
	}
      }
      break;

    case MESG_FRG:
      {
	/* Put the record where it belongs in the array.
	   This array is indexed by the overlaps. */
	FragMesg   *frg_mesg = (FragMesg *)pmesg->m;
	InternalFragMesg ifg_mesg;
	
#ifdef DEBUGIO
	fprintf(Msgfp,"Read FRG message with accession " F_UID "\n",
		frg_mesg->eaccession);
#endif
	if(GATEKEEPER_SUCCESS == 
	   Check_FragMesg(frg_mesg, &ifg_mesg, check_nmers, check_qvs, currentBatchID, currentTime, assembler,  verbose)){

	  pmesg->m = &ifg_mesg;
	  pmesg->t = MESG_IFG;
	  if((assembler == AS_ASSEMBLER_GRANDE) || (assembler == AS_ASSEMBLER_OBT)){
	     
	    if(ifg_mesg.type == AS_BACTIG || ifg_mesg.type == AS_FULLBAC){
	      bacFragmentsIgnored++;
	      Writer(Ignfp,pmesg);      
	    }else{
	      Writer(Outfp,pmesg);      
	    }
	  }else {
	    if(ifg_mesg.type == AS_UBAC || ifg_mesg.type == AS_FBAC || ifg_mesg.type == AS_LBAC){
	      shreddedFragmentsIgnored++;
	      Writer(Ignfp,pmesg);      
	    }else{
	      Writer(Outfp,pmesg);      
	    }
	  }
	}else{
	  fprintf(Msgfp,"# Line %d of input\n", GetProtoLineNum_AS());
	  ErrorWriter(Msgfp,pmesg);      
	  if(incrementErrors(1, Msgfp) == GATEKEEPER_FAILURE){
	    return GATEKEEPER_FAILURE;
	  }
	}

      }
      break;

    case MESG_LKG:
      {
	/* Put the record where it belongs in the array.
	   This array is indexed by the overlaps. */
	LinkMesg   *lnk_mesg = (LinkMesg *)pmesg->m;
	InternalLinkMesg ilk_mesg;
	int gkp_result;


#ifdef DEBUGIO
	fprintf(Msgfp,"Read LNK message with accessions "
		F_UID " , " F_UID "\n", 
		lnk_mesg->frag1, lnk_mesg->frag2);
#endif
	gkp_result = 
	  Check_LinkMesg(lnk_mesg, &ilk_mesg,  currentBatchID, currentTime,  verbose, matchBAC_DstsInLkgMsgs);

	switch(gkp_result){
	case GATEKEEPER_WARNING:
	  fprintf(Msgfp,"# Line %d of input\n", GetProtoLineNum_AS());
	  ErrorWriter(Msgfp,pmesg);      
	case GATEKEEPER_SUCCESS:
#ifdef OLD
	  /**************************** We DON'T output ILK messages anymore *********************************/
	  pmesg->m = &ilk_mesg;
	  pmesg->t = MESG_ILK;
	  Writer((assembler == AS_ASSEMBLER_OVERLAY?Ignfp:Outfp),pmesg);      
#endif
	  break;
	case GATEKEEPER_FAILURE:
	  fprintf(Msgfp,"# Line %d of input\n", GetProtoLineNum_AS());
	  ErrorWriter(Msgfp,pmesg);      
	  if(incrementErrors(1, Msgfp) == GATEKEEPER_FAILURE){
	    return GATEKEEPER_FAILURE;
	  }
	  break;
	}
      }
      break;




    case MESG_ADT:
      {
	AuditLine auditLine;
	AuditMesg *adt_mesg;
        char *params = (char *) malloc(200);
	char * emptyString = "";
        sprintf(params,
         "\nQV_MAX_ERROR: %5.4f\nQV_WINDOW_WIDTH: %3d\nQV_WINDOW_ERROR: %5.4f\n",
         GATEKEEPER_MAX_ERROR_RATE,
         GATEKEEPER_QV_WINDOW_WIDTH,
         GATEKEEPER_QV_WINDOW_THRESH);
	adt_mesg = (AuditMesg *)(pmesg->m);
	pmesg->t = MESG_ADT;

	VersionStampADTWithCommentAndVersion(
	    adt_mesg, 
            argc, 
            argv, 
	    // Next param assumed const -- Jason, 6/01.
            (GlobalParamText?GlobalParamText:emptyString),
            "(blank)");
        AppendAuditLine_AS(
            adt_mesg, 
	    &auditLine, 
	    currentTime, 
	    argv[0], "",
	    params);

	Writer((assembler == AS_ASSEMBLER_OVERLAY?Ignfp:Outfp),pmesg);
        /*
	fprintf(stderr,"# GateKeeper $Revision: 1.7 $\n");
	ErrorWriter(Msgfp,pmesg);
        */
        free(params);
	

      }
      break;
    case MESG_BAC:
      {
	BacMesg *bac_mesg = (BacMesg *)(pmesg->m);
	InternalBacMesg ibc_mesg;

#ifdef DEBUGIO
	fprintf(Msgfp,"Read BAC message with accession " F_UID "\n", 
		bac_mesg->ebac_id);
#endif

	if(GATEKEEPER_SUCCESS == 
	   Check_BacMesg(bac_mesg, &ibc_mesg, addedBactigs,  currentBatchID, currentTime, assembler, strict, verbose)){
	   pmesg->m = &ibc_mesg;
	   pmesg->t = MESG_IBC;
	   if(assembler == AS_ASSEMBLER_OVERLAY){
	     Writer((bac_mesg->type == AS_EBAC?Ignfp:Ibcfp),pmesg);      
	   }else{
	     Writer(Outfp,pmesg);      
	   }
	}else{
	  fprintf(Msgfp,"# Line %d of input\n", GetProtoLineNum_AS());
	   ErrorWriter(Msgfp,pmesg);      
	  if(incrementErrors(1, Msgfp) == GATEKEEPER_FAILURE){
	    return GATEKEEPER_FAILURE;
	  }
	}
      }
      break;

    case MESG_LIB:
      {
	LibDonorMesg *lib_mesg = (LibDonorMesg *) pmesg->m;

#ifdef DEBUGIO
	fprintf(Msgfp,"Read LIB message with accession " F_UID "\n", 
		lib_mesg->elibrary);
#endif

	if(GATEKEEPER_SUCCESS == 
	   Check_LibDonorMesg(lib_mesg, currentBatchID, currentTime, assembler, strict, verbose)){
          // nothing
	}else{
	  fprintf(Msgfp,"# Line %d of input\n", GetProtoLineNum_AS());
	   ErrorWriter(Msgfp,pmesg);      
	  if(incrementErrors(1, Msgfp) == GATEKEEPER_FAILURE){
	    return GATEKEEPER_FAILURE;
	  }
	}
      }
      break;
    case MESG_PLA:
      {
	PlateMesg *pla_mesg = (PlateMesg *) pmesg->m;

#ifdef DEBUGIO
	fprintf(Msgfp,"Read PLA message with accession " F_UID "\n",
                pla_mesg->eaccession );
#endif

	if(GATEKEEPER_SUCCESS == 
	   Check_PlateMesg(pla_mesg, currentBatchID, currentTime, assembler, strict, verbose)){
          // nothing
	}else{
	  fprintf(Msgfp,"# Line %d of input\n", GetProtoLineNum_AS());
	   ErrorWriter(Msgfp,pmesg);      
	  if(incrementErrors(1, Msgfp) == GATEKEEPER_FAILURE){
	    return GATEKEEPER_FAILURE;
	  }
	}
      }
      break;
    case MESG_LKP:
      {
        LinkPlateMesg * lkp_mesg = (LinkPlateMesg *) pmesg->m;
        
#ifdef DEBUGIO
	fprintf(Msgfp,"Read LKP message with forward accession " F_UID ", reverse accession " F_UID "\n",
                lkp_mesg->eplate_for, lkp_mesg->eplate_rev );
#endif

	if(GATEKEEPER_SUCCESS == 
	   Check_LinkPlateMesg(lkp_mesg, currentBatchID, currentTime, assembler, strict, verbose)){
          // nothing
	}else{
	  fprintf(Msgfp,"# Line %d of input\n", GetProtoLineNum_AS());
	   ErrorWriter(Msgfp,pmesg);      
	  if(incrementErrors(1, Msgfp) == GATEKEEPER_FAILURE){
	    return GATEKEEPER_FAILURE;
	  }
	}
      }
      break;
    default:
      fprintf(Msgfp,"# ERROR: Read Message with type %s line %d...skipping\n", MessageTypeName[imesgtype],
	      GetProtoLineNum_AS());
	  ErrorWriter(Msgfp,pmesg);      
	  if(incrementErrors(1, Msgfp) == GATEKEEPER_FAILURE){
	    return GATEKEEPER_FAILURE;
	  }
      break;
    }
  }

    if(assembler == AS_ASSEMBLER_OVERLAY)
      fprintf(stderr,"* Ignored %d shredded BAC Fragments\n", shreddedFragmentsIgnored);
    else
      fprintf(stderr,"* Ignored %d BACTIG/BAC Fragments\n", bacFragmentsIgnored);

  return(GATEKEEPER_SUCCESS);
}




/***************************************************************************************/
void printGKPError(FILE *fout, GKPErrorType type){

  switch(type){
  default:
#if 1
    fprintf(stderr,"#### printGKPError: error type %d\n", (int)type);
#endif
    break;

  case GKPError_Invalid:
  fprintf(stderr,"# printGKPError: Invalid error type %d\n", (int)type);
  break;

 case GKPError_FirstMessageBAT:
    fprintf(fout,"# GKP Error %d: First message MUST be BAT\n",(int)type);
    break;

    break;
  case GKPError_BadUniqueBAT:
    fprintf(fout,"# GKP Error %d: UID of batch definition was previously seen\n",(int)type);
    break;
  case GKPError_BadUniqueFRG:
    fprintf(fout,"# GKP Error %d: UID of fragment definition was previously seen\n",(int)type);
    break;
  case GKPError_BadUniqueDST:
    fprintf(fout,"# GKP Error %d: UID of distance definition was previously seen\n",(int)type);
    break;
  case GKPError_BadUniqueBAC:
    fprintf(fout,"# GKP Error %d: UID of BAC definition was previously seen\n",(int)type);
    break;
  case GKPError_BadUniqueBTG:
    fprintf(fout,"# GKP Error %d: UID of Bactig definition was previously seen\n",(int)type);
    break;
  case GKPError_BadUniqueSEQ:
    fprintf(fout,"# GKP Error %d: UID of Sequence definition was previously seen\n",(int)type);
    break;
  case GKPError_BadUniqueRPT:
    fprintf(fout,"# GKP Error %d: UID of Repeat definition was previously seen\n",(int)type);
    break;
  case GKPError_BadUniqueSCN:
    fprintf(fout,"# GKP Error %d: UID of Screen definition was previously seen\n",(int)type);
    break;
  case GKPError_BadUniqueWEL:
    fprintf(fout,"# GKP Error %d: UID of Well definition was previously seen\n",(int)type);
    break;
  case GKPError_BadUniquePLA:
    fprintf(fout,"# GKP Error %d: UID of Plate definition was previously seen\n",(int)type);
    break;

    
  case GKPError_MissingFRG:
    fprintf(fout,"# GKP Error %d: Fragment not previously defined\n",(int)type);
    break;

  case GKPError_MissingDST:
    fprintf(fout,"# GKP Error %d: Distance not previously defined\n",(int)type);
    break;
  case GKPError_MissingBAC:
    fprintf(fout,"# GKP Error %d: BAC not previously defined\n",(int)type);
    break;
  case GKPError_MissingBTG:
    fprintf(fout,"# GKP Error %d: Bactig not previously defined\n",(int)type);
    break;
  case GKPError_MissingSEQ:
    fprintf(fout,"# GKP Error %d: Sequence not previously defined\n",(int)type);
    break;
  case GKPError_MissingPLA:
    fprintf(fout,"# GKP Error %d: Plate not previously defined\n",(int)type);
    break;
  case GKPError_MissingRPT:
    fprintf(fout,"# GKP Error %d: Repeat not previously defined\n",(int)type);
    break;

  case GKPError_DeleteFRG:
    fprintf(fout,"# GKP Error %d: Can't delete Fragment\n",(int)type);
    break;
  case GKPError_DeleteDST:
    fprintf(fout,"# GKP Error %d: Can't delete Distance\n",(int)type);
    break;
  case GKPError_DeleteBAC:
    fprintf(fout,"# GKP Error %d: Can't delete BAC\n",(int)type);
    break;
  case GKPError_DeleteLNK:
    fprintf(fout,"# GKP Error %d: Can't delete LNK\n",(int)type);
    break;


  case GKPError_Time:
    fprintf(fout,"# GKP Error %d: Invalid creation time\n",(int)type);
    break;

  case GKPError_Action:
    fprintf(fout,"# GKP Error %d: Invalid action\n",(int)type);
    break;

  case GKPError_Scalar:
    fprintf(fout,"# GKP Error %d: Invalid scalar\n",(int)type);
    break;

  case GKPError_FRGSequence:
    fprintf(fout,"# GKP Error %d: Invalid fragment sequence characters\n",(int)type);
    break;

  case GKPError_FRGQuality:
    fprintf(fout,"# GKP Error %d: Invalid fragment quality characters\n",(int)type);
    break;
  case GKPError_FRGLength:
    fprintf(fout,"# GKP Error %d: Invalid fragment length\n",(int)type);
    break;
  case GKPError_FRGClrRange:
    fprintf(fout,"# GKP Error %d: Invalid fragment clear range must be 0<=clr1<clr2<=length\n",(int)type);
    break;

  case GKPError_FRGLocalPos:
    fprintf(fout,"# GKP Error %d: Invalid fragment locale pos\n",(int)type);
    break;

  case GKPError_FRGAccession:
    fprintf(fout,"# GKP Error %d: Invalid accession for BAC or Bactig fragment\n",(int)type);
    break;

  case GKPError_FRGBacSeqMismatch:
    fprintf(fout,"# GKP Error %d: Fragment BacID/SeqID mismatch\n",(int)type);
    break;

  case GKPError_LNKFragtypeMismatch:
    fprintf(fout,"# GKP Error %d: Link fragment type mismatch\n",(int)type);
    break;

  case GKPError_LNKLocaleMismatch:
    fprintf(fout,"# GKP Error %d: Link fragment locale mismatch\n",(int)type);
    break;

  case GKPError_LNKLocaleDistanceMismatch:
    fprintf(fout,"# GKP Error %d: Link fragment distance mismatch with associated locale\n",(int)type);
    break;

  case GKPError_LNKOneLink:
    fprintf(fout,"# GKP Error %d: Violation of unique mate/bacend link per fragment\n",(int)type);
    break;

  case GKPError_BACNumBactigs:
    fprintf(fout,"# GKP Error %d: Unfinished BAC must have >=1 Bactigs\n",(int)type);
    break;
  case GKPError_BACNumBactigsShouldBeZero:
    fprintf(fout,"# GKP Error %d: BAC of this type should have 0 Bactigs\n",(int)type);
    break;

  case GKPError_BACRedefinition:
    fprintf(fout,"# GKP Error %d: Error in BAC redefinition\n",(int)type);
    break;

  case GKPError_BACWrongTypeForOverlay:
    fprintf(fout,"# GKP Error %d: Overlay Assembler only accepts BACs of type UNFINISHED\n",(int)type);
    break;

  case GKPError_FRGWrongTypeForOverlay:
    fprintf(fout,"# GKP Error %d: Overlay Assembler only accepts BAC FRGs of type BACTIG\n",(int)type);
    break;

  case GKPError_DSTValues:
    fprintf(fout,"# GKP Error %d: DST mean,stddev must be >0 and mean must be >= 3 * stddev\n",(int)type);
    break;
  
  case GKPError_SCNVariation:
    fprintf(fout,"# GKP Error %d: SCN variation out of bounds\n",(int)type);
    break;
  case GKPError_SCNminLength:
    fprintf(fout,"# GKP Error %d: SCN min length out of bounds\n",(int)type);
    break;

  case GKPError_RPTLength:
    fprintf(fout,"# GKP Error %d: RPT length out of bounds\n",(int)type);
    break;

  case GKPError_IncompleteOAInput:
    fprintf(fout,"# GKP Error %d: Overlay Assembler Requires that a Bactig and its\n",(int)type);
    fprintf(fout,"#               associated fragment must appear in the same batch\n");
    break;
    }
}


void printAllGKPErrors(FILE *fout){
  int i;
  for(i =1; i <= MAX_GKPERROR; i++){
    printGKPError(stderr, (GKPErrorType)i);
 } 
}


/***********************************************************************************/
/* CheckOverlayInputComplete                                                       */
/***********************************************************************************/
int CheckOverlayInputComplete(VA_TYPE(int32) *addedBactigs, int verbose){
  size_t i;
  int status = GATEKEEPER_SUCCESS;

  for(i = 0; i < GetNumVA_int32(addedBactigs); i++){
    GateKeeperBactigRecord gkpbactig;
    int32 bactigIid = *GetVA_int32(addedBactigs, i);

    getGateKeeperBactigStore(GkpStore.btgStore,bactigIid, &gkpbactig);

    if(gkpbactig.hasSequence == FALSE){
      status = GATEKEEPER_FAILURE;
      printGKPError(Msgfp,GKPError_IncompleteOAInput);
      fprintf(Msgfp,"# Bactig (" F_UID "," F_S32 ") does not have an associated FRG\n",
	      gkpbactig.UID, bactigIid);
    }
  }
  return status;
}


int32   CheckNmerProbabilities(FILE *fout, double threshhold){
  float32 *probs;

  ComputeBucketActualRates(SequenceProbabilities);
  probs = SequenceProbabilities->arate;

  fprintf(stderr,"* Global probabilities  P(a) = %6g P(c) = %6g P(g) = %6g P(t) = %6g\n",
	  probs[0], probs[1], probs[2], probs[3]);

  CheckSequenceBucketArraySanity(LinkerDetector_READ, SanityDetector_READ, probs, threshhold, fout,"#linker5p READ");
  CheckSequenceBucketArraySanity(Linker3pDetector_READ, SanityDetector_READ, probs, threshhold, fout,"#linker3p READ");
  //  CheckSequenceBucketArraySanity(LinkerDetector_LBAC, SanityDetector_LBAC, probs, threshhold, fout, "#linker LBAC");
  //  CheckSequenceBucketArraySanity(LinkerDetector_EBAC, SanityDetector_EBAC, probs, threshhold, fout,"#linker EBAC");

  return GATEKEEPER_SUCCESS;


}

 static void usage(void){
	fprintf (stderr, "USAGE:  gatekeeper [-aiefnopsCGPNQX] <Output_Store_Name> <Input>.<ext>\n"
		 "Opens <Input>.<ext> to read .frg input\n"
		 "Creates GateKeeperFragmentStore in <Output_Store_Name>\n"
		 "Writes output to <InputFileName>.inp\n"
		 "  -a  append to Store\n"
		 "  -b  batchNumFileNames on\n"
		 "  -c  create new Store\n"
		 "  -d <8-mer>  Do clear range length histograms for fragments with the nmer at their 3p,5p,Sanity\n"
		 "  -e <errorThreshhold>  set error threshhold\n"
		 "  -f  force new Store\n"
		 "  -h  print usage and error messages\n"
		 "  -i  <input store>\n"
		 "  -n  <projectname>\n"
		 "  -o  <output store>\n"
		 "  -s  strict enforcement of -G -O rules\n"
		 "  -t <float>  threshhold for frequent n-mer check in units of std deviations\n"
		 "  -G  gatekeeper for assembler Grande (default)\n"
		 "  -O  gatekeeper for overlay assembler\n"
		 "  -T  gatekeeper for assembler Grande with Overlap Based Trimming\n"
		 "  -C  compatiblity mode\n"
		 "  -N  don't check n-mer frequencies\n"
		 "  -Q  don't check quality-based data quality\n"
		 "  -X  enable experimental switches\n"
                 "  -l  don't check consistency of EBAC link distance UIDs\n\n"

		 );
	exit (EXIT_FAILURE);
 }
