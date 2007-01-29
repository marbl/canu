
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
static char CM_ID[] = "$Id: PopulateFragStore.c,v 1.10 2007-01-29 20:41:20 brianwalenz Exp $";

/*************************************************
* Module:  PopulateFragStore.c
* Description:
*    Populate Frag Store with Batch of Fragments.  This is basically equivalent to "overlap -n",
*    but with the added benefit of adhering to the new Data Flow spec.
* 
*    Reference: DataFlow document
*
*    Command Line Interface:
*        PopulateFragStore [-fcav] [-i <InputStorePath>] [-o <OutputStorePath>] <inputFile>.<ext>
*
*        -a append to the fragment store (it it doesn't exist, create it )
*        -A same as -a, but no backup is made
*        -f force clobber an existing store if it exists
*        -c create
*        -i <InputStorePath>
*        -o <OutputStorePath>
*
*    Programmer:  S. Kravitz
*       Written:  Mar 2000
* 
*************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <ctype.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_distStore.h"
#include "AS_PER_fragStore.h"
#include "AS_UTL_version.h"


/**** Static Global Variables ****/

static Batch_ID  Batch_Msg_UID = 0;
static IntBatch_ID  Batch_Msg_IID = 0;


/**** Utility Functions ****/

static void usage(void);
static void  stripWhiteSpace (char *target, char *source, int maxlen);

static FILE *  File_Open (const char * Filename, const char * Mode);
static int  ReadFrags(int maxFrags, 
		      FragStoreHandle oldFragStore, DistStore oldDistStore, 
		      FragStoreHandle newFragStore, DistStore newDistStore, 
		      FILE *In_Stream, 
		      FILE *Out_Stream, 
		      int argc, char **argv,
		      int * fragsRead, int *distsRead);

FragStoreHandle  Frag_Store_Open
    (char * storename, char * mode);
DistStore  Dist_Store_Open
    (char * storename, char * mode);



/* Buffers for file names */
char  Input_Store_Name [FILENAME_MAX];
char  Output_Store_Name [FILENAME_MAX];
char  Input_File_Name [FILENAME_MAX];
static char  Output_File_Name [FILENAME_MAX] = "";
char  Stats_File_Name [FILENAME_MAX];
char  File_Name_Prefix [FILENAME_MAX];
char  Parameter_File_Name [FILENAME_MAX];


int  main(int argc, char * argv [])

{
  char cmd[FILENAME_MAX * 4] = "";
  char tmpFilePath[FILENAME_MAX] = "";
  FILE *In_fp = NULL, *Out_fp = NULL;
  FragStoreHandle tmpFragStore, targetFragStore;
  DistStore tmpDistStore, targetDistStore;

  int verbose = 0;
  int create = 0;
  int append = 0;
  int force = 0;
  int backup = 1;
  int input_exists = 0;
  int output_exists = 0;
  int inputStoreSpecified = FALSE;
  int outputStoreSpecified = FALSE;
  int64 firstFrag, lastFrag;
  int64 firstDist, lastDist;


  /**************** Process Command Line Arguments *********************/
  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    optind = 1;
    while (!errflg && ((ch = getopt(argc, argv, "aAchfi:o:Pp:V:vX")) != EOF)){
      switch(ch) {
      case 'a':
	append = 1;
	create = 0;
	force  = 0;
	backup = 1;
	break;
      case 'A':
	append = 1;
	backup = 0;
	create = 0;
	force  = 0;
	break;
      case 'c':
	backup = 0;
	create = 1;
	break;
      case 'f':
	backup = 0;
	force  = 1;
	append = 0;
	break;
      case 'h':
	usage();
	exit(0);
	break;
      case 'i':
	fprintf(stderr,"* -i %s\n", optarg);
	strcpy(Input_Store_Name, optarg);
	inputStoreSpecified = TRUE;
	//errflg = inputStoreSpecified;
	break;
      case 'o':
	fprintf(stderr,"* -o %s\n", optarg);
	strcpy(Output_Store_Name, optarg);
	outputStoreSpecified = TRUE;
	//errflg = inputStoreSpecified;
	break;
      case 'P':
	fprintf(stderr, "-P is no longer needed.\n");
	break;
      case 'p':
	fprintf(stderr,"* -p %s\n", optarg);
	strcpy(Parameter_File_Name, optarg);
        fprintf(stderr,"** The parameter file is currently ignored!!!\n");
	break;
      case 'v':
	verbose = 1;
	break;
      case  'V' :
        assert (strlen (optarg) < FILENAME_MAX - 1);
        strcpy (Output_File_Name, optarg);
        break;
      case 'X':
        // Activate the expert options.
	break;
      case '?':
      default :
	fprintf(stderr,"Unrecognized option -%c",optopt);
	errflg++;
	exit(1);
      }
    }
    /* End of command line parsing */
       
    /*** Check -afcio options for meaningful combinations ***/
    if(force == 1 && append  == 1) {
      fprintf (stderr,
               "* Illegal combination of command line flags"
               "-- they are mutually exclusive\n");
      errflg = 1;
    }
    
    /* If we've seen an error, or there isn't
       one and only one argument left...exit. */
    if(errflg || (argc - optind != 1 )) {
      fprintf(stderr,
              "If we've seen %d errors, or there isn't "
              "one and only one argument (%d) left...exit\n",
              errflg,argc - optind);
      usage();
      exit(1);
    }

  /**** Open input and output files ****/
  {
    char *suffix;

    suffix = strrchr(argv[optind],(int)'.');
       
    fprintf(stderr,"Input file is %s suffix is %s\n",argv[optind], suffix);
    strcpy(Input_File_Name, argv[optind]);

    In_fp = File_Open (Input_File_Name, "r");     // frg file
    AssertPtr(In_fp);

    if(suffix)
      *suffix = '\0';

    strcpy(File_Name_Prefix,argv[optind]);

    if  (Output_File_Name [0] == '\0')
        sprintf (Output_File_Name, "%s.ovl", File_Name_Prefix);
    fprintf(stderr,"* Output_File_Name %s\n",Output_File_Name);

    Out_fp = File_Open (Output_File_Name, "w");     // ovl file
    AssertPtr(Out_fp);

    optind++;
  }
}

{
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
    usage();
    exit(1);
  }
  /* Open or create input Store */
  if(create){
    output_exists = (existsFragStore(Output_Store_Name) != 0);

    if(output_exists){
      fprintf(stderr,"* Frag store %s exists...\n", Output_Store_Name);
      if(!force){
	fprintf(stderr,"* Output store exists. You must use the -f option to clobber it...exiting\n");
	exit(1);
      }else{
	fprintf(stderr,"* Nuking\n");
	removeFragStore (Output_Store_Name);
	sprintf (cmd,"rm -f %s/db.dst", Output_Store_Name);
	system (cmd);
	append = 0;
      }
    }
    if(! output_exists){ /* The directory doesn't exist */
      fprintf(stderr,"* Creating directory %s in cwd of %s\n", Output_Store_Name, getcwd(NULL, 256));
      mkdir(Output_Store_Name, S_IRWXU | S_IRWXG | S_IROTH);
    }
    fprintf(stderr,"* Creating fragment store %s\n", Output_Store_Name);
    targetFragStore = createFragStore(Output_Store_Name,"",1);
    {
      char buffer[2048];
      sprintf(buffer,"%s/db.dst", Output_Store_Name);
      targetDistStore = createDistStore(buffer,1);
    }
    firstFrag = 1;
    firstDist = 1;
  }else if(append){
    strcpy(Output_Store_Name, Input_Store_Name);
    fprintf(stderr,"* Appending to fragment store %s\n", Input_Store_Name);

    if(backup){
      sprintf(tmpFilePath,"%s___tmp", Input_Store_Name);

      /* Remove temporary files if any */
      sprintf(cmd,"rm -rf %s", tmpFilePath);
      if(system(cmd) != 0) assert(0);
      fprintf(stderr,"* Creating directory %s in cwd of %s\n", tmpFilePath, getcwd(NULL, 256));
      mkdir(tmpFilePath, S_IRWXU | S_IRWXG | S_IROTH);

      /* Copy the frag and dst stores into the tmp directory and open them */
      copyFragStore(Input_Store_Name, tmpFilePath,FALSE);
      targetFragStore = Frag_Store_Open (tmpFilePath, "r+");
      {
	char buffer[2048];
	sprintf(buffer,"%s/db.dst", Output_Store_Name);
	sprintf(cmd,"cp %s/db.dst %s",Input_Store_Name,tmpFilePath);
	if(system(cmd) != 0) assert(0);

	sprintf(buffer,"%s/db.dst", tmpFilePath);
        fprintf(stderr,"* Opening dist store %s\n", tmpFilePath);
	targetDistStore = Dist_Store_Open(buffer,"r+");
      }
    }else{ // no backup

      /* Open the input store for r/w */
      targetFragStore = Frag_Store_Open (Input_Store_Name, "r+");
      {
	char buffer[2048];
	sprintf(buffer,"%s/db.dst", Input_Store_Name);
	targetDistStore = Dist_Store_Open(buffer,"r+");
      }

    }

    /* Remember the state of the stores before we start */
        {
          StoreStat stats;

          statsFragStore(targetFragStore, &stats);
          firstFrag = (stats.lastElem >= stats.firstElem?stats.lastElem + 1:1);
          statsStore(targetDistStore, &stats);
          firstDist = (stats.lastElem >= stats.firstElem?stats.lastElem + 1:1);
#if 0
	firstDist = getLastElemStore(targetDistStore) + 1;
	firstFrag = getLastElemFragStore(targetFragStore) + 1;
#endif
        }

  }else{ // normal case, reading from input writing to output
	output_exists = (testOpenFragStore(Output_Store_Name,"r") == 0);

	if(output_exists){
	  fprintf(stderr,"* Frag store %s exists...\n", Output_Store_Name);
	  if(!force){
	    fprintf(stderr,"* Output store exists...exiting\n");
	    exit(1);
	  }else{
	    fprintf(stderr,"* Nuking\n");
	    removeFragStore(Output_Store_Name);
	  }
	}else{

	  mkdir(Output_Store_Name, S_IRWXU | S_IRWXG | S_IROTH);

	}
	fprintf(stderr,"* TestOpen of input: %s\n", Input_Store_Name);
	input_exists = testOpenFragStore(Input_Store_Name,"r");

	if(input_exists != 0){
	  fprintf(stderr,"* Can't open input store...exiting\n");
	  exit(1);
	}
	copyFragStore(Input_Store_Name,Output_Store_Name, FALSE);
	sprintf(cmd,"cp %s/db.dst %s", Input_Store_Name, Output_Store_Name);
	fprintf(stderr,"* %s\n", cmd);
	if(system(cmd) != 0) assert(0);

	targetFragStore = Frag_Store_Open(Output_Store_Name,"r+");
	{
	  char buffer[2048];
	  sprintf(buffer,"%s/db.dst", Output_Store_Name);
	  targetDistStore = Dist_Store_Open(buffer,"r+");
	}
        {
          StoreStat stats;

          statsFragStore(targetFragStore, &stats);
          firstFrag = (stats.lastElem >= stats.firstElem?stats.lastElem + 1:1);
          statsStore(targetDistStore, &stats);
          firstDist = (stats.lastElem >= stats.firstElem?stats.lastElem + 1:1);
#if 0
	firstDist = getLastElemStore(targetDistStore) + 1;
	firstFrag = getLastElemFragStore(targetFragStore) + 1;
#endif
        }
      }

  }

  fprintf(stderr,
          "* Expecting fragment " F_S64 " and distance " F_S64 "\n",
          firstFrag, firstDist);

  tmpFragStore = createFragStore(NULL,"tmp",firstFrag);
  tmpDistStore = createDistStore(NULL,firstDist);


  {
    int distRead, fragRead;

#define FRAGS_PER_BATCH 50000

  while(ReadFrags(FRAGS_PER_BATCH, 
		  targetFragStore, targetDistStore,
		  tmpFragStore, tmpDistStore,
		  In_fp, Out_fp, 
		  argc,argv, 
		  &fragRead, &distRead)){
    if(fragRead){
      fprintf(stderr,"* Fragsread %d\n", fragRead);
                concatFragStore (targetFragStore, tmpFragStore);
                fprintf(stderr,"* Resetting fragstore to " F_S64 "\n",
                        getLastElemFragStore(targetFragStore) + 1);
		resetFragStore (tmpFragStore, getLastElemFragStore(targetFragStore) + 1);
    }

    if(distRead){
      fprintf(stderr,"* Distsread %d\n", distRead);
                concatStore (targetDistStore, tmpDistStore);
                fprintf(stderr,"* Resetting distStore to " F_S64 "\n",
                        getLastElemStore(targetDistStore) + 1);
		resetDistStore (tmpDistStore, getLastElemStore (targetDistStore) + 1);
    }
  }

   
}
  lastFrag = getLastElemFragStore(targetFragStore);
  lastDist = getLastElemStore(targetDistStore);

   if  (firstFrag >= 0)
     {
       char buffer[2048];
       sprintf(buffer,"%s.range",File_Name_Prefix);
       {
         FILE  * fp = File_Open (buffer, "w");

         fprintf (fp, "batch:" F_UID "\n", Batch_Msg_UID);
         fprintf (fp, "lo:" F_S64 "  hi:" F_S64 "\n", firstFrag, lastFrag);
         fclose (fp);
       }
     }

  closeFragStore(targetFragStore);
  closeDistStore(targetDistStore);

  fclose(In_fp);
  fclose(Out_fp);

  if(append && backup){
    removeFragStore (Output_Store_Name);
    sprintf (cmd,"rm -f %s/db.dst", Output_Store_Name);
    fprintf(stderr,"* %s\n", cmd);
    if(system(cmd) != 0) assert(0);
    copyFragStore(tmpFilePath,Output_Store_Name, TRUE);
    sprintf(cmd,"cp %s/db.dst %s", tmpFilePath, Output_Store_Name);
    fprintf(stderr,"* %s\n", cmd);
    if(system(cmd) != 0) assert(0);
    sprintf(cmd,"rm -rf %s", tmpFilePath);
    if(system(cmd) != 0) assert(0);
  }

 return  0;
}



/************************************************************************************
 *  ReadFrags:
 *  Read up to maxFrags, but maybe a little less.
 *
 ************************************************************************************/

static int  ReadFrags(int maxFrags, 
		      FragStoreHandle oldFragStore, DistStore oldDistStore, 
		      FragStoreHandle newFragStore, DistStore newDistStore, 
		      FILE *In_Stream,    FILE *Out_Stream, 
		      int argc, char **argv,
		      int * fragsRead, int *distsRead)

  {
  ReadStructp  myRead = new_ReadStruct();
  int numFragsRead = 0;
  int numDistsRead = 0;
  StoreStat stats;
  MessageType imesgtype;
  GenericMesg   *pmesg; 
  //OverlapMesg   *ovl_mesg;
  int64  total_len = 0;
  int64 firstFrag, nextFrag;
  int64 firstDist, nextDist;

  AssertPtr(In_Stream);

  statsFragStore(oldFragStore, &stats);
  nextFrag = firstFrag = (stats.lastElem >= stats.firstElem?stats.lastElem + 1:stats.firstElem);
  statsFragStore(newFragStore, &stats);

  if(firstFrag != stats.firstElem) {
    fprintf(stderr,"firstFrag = " F_S64 " stats.firstElem=" F_S64 "\n",
            firstFrag, stats.firstElem);
  }
  assert(firstFrag == stats.firstElem);

  statsStore(oldDistStore, &stats);
  fprintf(stderr,"* oldDistStore first:" F_S64 " last:" F_S64 "\n",
          stats.firstElem, stats.lastElem);
  nextDist = firstDist = (stats.lastElem >= stats.firstElem?stats.lastElem + 1:stats.firstElem);
  statsStore(newDistStore, &stats);
  fprintf(stderr,"* newDistStore first:" F_S64 " last:" F_S64 "\n",
          stats.firstElem, stats.lastElem);
  if(firstDist != stats.firstElem) {
    fprintf(stderr,"firstDist = " F_S64 " stats.firstElem=" F_S64 "\n",
            firstDist, stats.firstElem);
  }
  assert(firstDist == stats.firstElem );

  AssertPtr(In_Stream);

  /* Read maxFrags fragments from the stream */

  while  ((numFragsRead < maxFrags)    && 
	  (EOF != ReadProtoMesg_AS (In_Stream, & pmesg)))
    {
    imesgtype = pmesg->t;
    switch(imesgtype){
    case MESG_IDT:
      {
	DistRecord distRecord;
	InternalDistMesg  *idt_mesg = (InternalDistMesg *) pmesg->m;
#ifdef DEBUG
	fprintf(stderr,"Read DST message %c ( " F_UID ", " F_IID " ) nextDist = " F_S64 "\n", 
		idt_mesg->action,
		idt_mesg->eaccession,
		idt_mesg->iaccession,
		nextDist);
#endif
	distRecord.UID = idt_mesg->eaccession;
	distRecord.IID = idt_mesg->iaccession;
	switch(idt_mesg->action){
	case AS_ADD:
	  {
	    if(distRecord.IID != nextDist){
	      fprintf(stderr,"*** Fatal Error -- distance record should have IID " F_S64 " not " F_IID "\n",
		      nextDist, distRecord.IID);
	      exit(1);
	    }

	    nextDist++;
	    numDistsRead++;
	    distRecord.mean = idt_mesg->mean;
	    distRecord.stddev = idt_mesg->stddev;

            /*
	    fprintf(stderr,"* Distance Message (" F_UID "," F_IID ")\n", 
                    distRecord.UID, distRecord.IID);
            */

	    appendDistStore(newDistStore, &distRecord);
	  }
	  break;
	case AS_REDEFINE:
	  {
	    distRecord.mean = idt_mesg->mean;
	    distRecord.stddev = idt_mesg->stddev;

            /*
	    fprintf(stderr,"* Distance Message (" F_UID "," F_IID ")\n", 
                    distRecord.UID, distRecord.IID);
            */

	  if  (distRecord . IID > nextDist){
	       fprintf (stderr, "*** Fatal Error -- distance record " F_IID " doesn't exist\n",
		        idt_mesg -> iaccession);
	       exit (EXIT_FAILURE);
	  } else if  (idt_mesg -> iaccession < firstDist){
	    setDistStore(oldDistStore, distRecord.IID,&distRecord);
	  }else{
	    setDistStore(newDistStore, distRecord.IID,&distRecord);
	  }
	  }
	  break;
	case AS_DELETE:
	  if  (distRecord . IID > nextDist){
	       fprintf (stderr, "*** Fatal Error -- distance record " F_IID " doesn't exist\n",
		        idt_mesg -> iaccession);
	       exit (EXIT_FAILURE);
	  } else if  (idt_mesg -> iaccession < firstDist){
	    deleteDistStore (oldDistStore, idt_mesg -> iaccession);
	  }else{
	    deleteDistStore (newDistStore, idt_mesg -> iaccession);
	  }
          break;

        default:
          assert(0);
        }

	WriteProtoMesg_AS (Out_Stream,pmesg);
      }
      break;

    case MESG_IFG:
      {

	InternalFragMesg *ifg_mesg = (InternalFragMesg *) pmesg->m;

        OFGMesg ofg_mesg;
        Transfer_IFG_to_OFG_AS (ifg_mesg, &ofg_mesg);
        pmesg->m = &ofg_mesg;
        pmesg->t = MESG_OFG;

#ifdef DEBUG
        fprintf(stderr,"Read IFG message %c ( " F_UID ", " F_IID " ) \n",
                ifg_mesg->action,
                ifg_mesg->eaccession,
                ifg_mesg->iaccession);
#endif
          switch(ifg_mesg->action){
          case AS_ADD:
            numFragsRead++;
#ifdef DEBUG
            fprintf(stderr,
                    "Read message %c (" F_UID ", " F_IID ") " F_S64 "\n", 
                    ifg_mesg->action, ifg_mesg->eaccession,
                    ifg_mesg->iaccession, nextFrag);
#endif
            if(ifg_mesg->iaccession != nextFrag){
              fprintf(stderr,"*** Fatal Error -- fragment record should have IID " F_S64 " not " F_IID "\n",
                      nextFrag, ifg_mesg->iaccession);
              exit(1);
            }

            nextFrag++;

            {
              char buffer1[AS_READ_MAX_LEN * 2];
              char buffer2[AS_READ_MAX_LEN * 2];

              /* Add it to the fragment Store, as well */
              setAccID_ReadStruct(myRead, ofg_mesg.eaccession);
              setReadIndex_ReadStruct(myRead, ofg_mesg.iaccession);
              setReadType_ReadStruct(myRead, ofg_mesg.type);
              stripWhiteSpace(buffer1, ifg_mesg->sequence, AS_READ_MAX_LEN * 2);
              stripWhiteSpace(buffer2, ifg_mesg->quality, AS_READ_MAX_LEN * 2);
              setSequence_ReadStruct(myRead, buffer1, buffer2);
              setSource_ReadStruct(myRead, ifg_mesg->source);
              setEntryTime_ReadStruct(myRead, ofg_mesg.entry_time);
              setClearRegion_ReadStruct(myRead, ofg_mesg.clear_rng.bgn, ofg_mesg.clear_rng.end, READSTRUCT_ORIGINAL);
              setLocalePos_ReadStruct(myRead,ofg_mesg.locale_pos.bgn, ofg_mesg.locale_pos.end);
              setLocID_ReadStruct(myRead,ofg_mesg.ilocale);
              total_len += ofg_mesg.clear_rng.end - ofg_mesg.clear_rng.bgn;;

              appendFragStore(newFragStore, myRead);
            }
          break;

        case AS_DELETE:
          if  (ifg_mesg -> iaccession < firstFrag)
              deleteFragStore (oldFragStore, ifg_mesg -> iaccession);
            else
              deleteFragStore (newFragStore, ifg_mesg -> iaccession);
	  break;

	default:
          fprintf (stderr,
                   "ERROR:  Bad action = %d  on frag iid = " F_IID "\n",
                   (int) ifg_mesg -> action, ifg_mesg -> iaccession);
	  assert(0);
	  break;
	  }

	  WriteProtoMesg_AS (Out_Stream, pmesg);
	}
	break;

      case MESG_ADT:
	{
	  AuditMesg *adt_mesg = (AuditMesg *) pmesg->m;
	  pmesg->t = MESG_ADT;

	  VersionStampADT(adt_mesg, argc, argv);

	  WriteProtoMesg_AS (Out_Stream, pmesg);
	}
	break;

    case MESG_IBA:
       {
        InternalBatchMesg  * iba_mesg = (InternalBatchMesg  *) pmesg -> m;
        Batch_Msg_UID = iba_mesg -> eaccession;
        Batch_Msg_IID = iba_mesg -> iaccession;

        // fall through
       }

    case MESG_IBC:

#ifdef DEBUG
	fprintf(stderr,"Read IBC/IBA\n");
#endif

	WriteProtoMesg_AS (Out_Stream, pmesg);
	break;


      default:
	fprintf(stderr,"* Oops: Read Message with type imesgtype = %d\n",
                imesgtype);
	WriteProtoMesg_AS(stderr,pmesg);      
	
	exit(1);
	}
      }


 *fragsRead = numFragsRead;
 *distsRead = numDistsRead;
  delete_ReadStruct(myRead);
  return(numFragsRead + numDistsRead);
}


/******************************************************************************/

/* Function stripWhiteSpace:
   Input:  source   string of maximum length maxlen
           maxlen   maximum length of source
   Output: target
   
   Description:
     Copy non-space characters from source to target.
*/

void  stripWhiteSpace
    (char *target, char *source, int maxlen)

  {
  int i = 0;
  *target = '\0';
  while(i < maxlen){
    if(!isspace(*source)){
      *target++ = *source;
      i++;
    }
    if(*source == '\0')
      break;
    source++;
  }

}



#define  FILE_OPEN_FAILURE   134

FragStoreHandle  Frag_Store_Open
    (char * storename, char * mode)

/* Open  storename  in  mode  and return a pointer to its control
*  block.  If fail, print a message and exit. */

  {
   FragStoreHandle  fp;
   
   fp = openFragStore (storename, mode);
   if  (fp == NULLSTOREHANDLE)
       {
        fprintf (stderr, "ERROR %d:  Could not open fragstore  %s \n",
                 errno, storename);
        perror (strerror (errno));
       }

   return  fp;
  }




FILE *  File_Open (const char * Filename, const char * Mode)

/* Open  Filename  in  Mode  and return a pointer to its control
*  block.  If fail, print a message and exit. */

  {
   FILE  *  fp;

   fp = fopen (Filename, Mode);
   if  (fp == NULL)
       {
        fprintf (stderr, "ERROR %d:  Could not open file  %s \n",
                 errno, Filename);
        perror (strerror (errno));
       }

   return  fp;
  }



int File_Exists (const char * Filename){

  /* Test for filename existence */

   FILE  *  fp;

   fp = fopen (Filename, "r+");
   if  (fp)
       {
	 fclose(fp);
	 return 1;
       }

   return  0;

}


void usage(void){
  fprintf (stderr,
           "usage: Populator [-aAcfX] [-i <inputStore>] [-o <outputStore>]\n"
           "                 [-V <ovlfilename>] <inputFile>.<ext>\n"
	   " Populates a fragment store with the fragments and distance messages in the input file. \n"
	  " Outputs two files:\n"
	  "\t<inputFile>.ovl    All proto messages, with IFGs and SFGs converted to OFGs\n"
	  "\t<inputFile>.range  Range of fragment IIDs processed in the <inputFile>"
	  "Options:\n"
	  "  -i    specifies an input fragStore\n"
	  "  -o    specifies an output fragStore\n"
	  "  -a    append with backup (requires -i)\n"
	  "  -A    append with no backup (requires -i)\n"
	  "  -c    create (requires -o)\n"
	  "  -f    force.  If output exists, nuke it.\n"
          "  -V    specifies name of output (.ovl) file\n"
          "  -X    activate expert options.\n");
}

DistStore  Dist_Store_Open
    (char * storename, char * mode)

/* Open  storename  in  mode  and return a pointer to its control
*  block.  If fail, print a message and exit. */

  {
   DistStore  fp;

   fp = openDistStore (storename, mode);
   if  (fp == NULLSTOREHANDLE)
       {
        fprintf (stderr, "ERROR %d:  Could not open diststore  %s \n",
                 errno, storename);
        perror (strerror (errno));
       }

   return  fp;
  }

