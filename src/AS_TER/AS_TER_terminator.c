
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
/*********************************************************************		Module: AS_TER_terminator.c

 Description: Assembly terminator module. It is the backend of the
assembly pipeline and replaces internal accession numbers by external
accession numbers.

 Assumptions: Input meets specifications in the ProtoIO documents.

**********************************************************************/

/* RCS info
 * $Id: AS_TER_terminator.c,v 1.2 2004-09-23 20:25:29 mcschatz Exp $
 * $Revision: 1.2 $
 */


static char CM_ID[] = "$Id: AS_TER_terminator.c,v 1.2 2004-09-23 20:25:29 mcschatz Exp $";


/*************************************************************************/
/* System include files */
/*************************************************************************/

#include  <stdlib.h>
#include  <stdio.h>
#include  <unistd.h>
#include  <assert.h>
#include  <sys/types.h>
#include  <string.h>

/*************************************************************************/
/* Local include files */
/*************************************************************************/

#include "AS_global.h"
#include "AS_TER_terminator.h"
#include "AS_TER_utils.h" 
#include "AS_TER_terminator_funcs.h"
#include "AS_TER_alloc.h"


int main (int argc, char *argv[]) {
  char *inputFileName   = NULL;
  char *outputFileName  = NULL;
  char *mapFileName     = NULL;
  char *fragStoreName   = NULL;
  char *bactigStoreName = NULL;
  char *gkpStoreName    = NULL;
  char *euidServerNames = NULL;
  /*Variable defintions for parsing the commandline */

  /* The default behaviour of the Terminator. Should
     be changed to binary */
  OutputType  output = AS_BINARY_OUTPUT;
  int32 illegal = 0;
  int32 realUIDs = FALSE;
  int32 help     = FALSE;
  uint32 random  = FALSE;
  uint32 quiet   = FALSE;
  uint64 uidStart = 0;

  { /* Parse the argument list using "man 3 getopt". */
    int ch,errflg=FALSE;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv, "Pub:f:g:i:o:m:rs:E:hQ")) != EOF))
      switch(ch) 
	{
	case 'P':
	  output = AS_PROTO_OUTPUT;
	  break;
	case 'u':	  
	  realUIDs = TRUE;
	  break;
	case 'r':	  
	  random = TRUE;
	  break;
	case 'Q':	  
	  quiet = TRUE;
	  break;
	case 's':	  
	  uidStart = strtoul(optarg,(char**) NULL,10);
	  break;
	case 'f':	  
	  fragStoreName = strdup(optarg);
	  assert(fragStoreName != NULL);
	  break;
	case 'b':	  
	  bactigStoreName = strdup(optarg);
	  assert(bactigStoreName != NULL);
	  break;
	case 'g':	  
	  gkpStoreName = strdup(optarg);
	  assert(gkpStoreName != NULL);	
	  break;
	case 'i':	  
	  inputFileName = strdup(optarg);
	  assert(inputFileName != NULL);
#if 0
	  fprintf(stderr,"**** -i option is no longer supported!  Please specify input file(s) at end of command line\n");
	  fprintf(stderr,"**** to read from stdin, specify no files at end of command line\n");
	  exit(1);
#endif
	  break;
	case 'o':	  
	  outputFileName = strdup(optarg);
	  assert(outputFileName != NULL);
	  break;
	case 'm':	  
	  mapFileName = strdup(optarg);
	  assert(mapFileName != NULL);
	  break;
#ifdef USE_SOAP_UID
	case 'E':
	  realUIDs = TRUE;
	  SYS_UIDset_euid_server(optarg);
	  break;
#endif
	case 'h':	  
	  help = TRUE;
	  break;
	default:
	  fprintf(stderr,"Unrecognized option -%c\n",optopt);
	  errflg++;
	  illegal = 1;
	  break;
	}
  }
  
  if((illegal == 1)  || help || (realUIDs != TRUE && uidStart == 0 && random)){
    fprintf (stderr,
	     "USAGE:  terminator [-P -u -b<bactig_store_directory> -g <gatekeeper_store_directory] <input files>*\n"
	     "specify zero or more input files at end of command line.  Zero files means read from stdin\n"
	     "-f <frag_store_directory> [-i <input_file>] -o <output_file> -m <map_file> \n"
	     "-s <uid_start> sets dummy UIDs start has to be > 0 (only required for -r)\n"
	     "-r puts terminator in random access mode\n"
	     "-P forces ASCII output\n"
#ifdef USE_SOAP_UID
	     "-E <server:port> changes default EUID server from tools.tigr.org:8190 (implies -u)\n"
#endif
	     "-u forces real UIDs\n"
	     "-Q runs also for simulator\n");
    exit (AS_TER_EXIT_FAILURE);
  }
  
  {
    char *stdinInput = "-";
    char **inputList;
    int32 numInputs;
    if(inputFileName != NULL){
      inputList = &inputFileName;
      numInputs = 1;
      if(argc - optind > 0){
	int i;
	fprintf(stderr,"* -i option specified, ignoring argments: ");
	for(i = optind; i < argc; i++){
	  fprintf(stderr,"%s ", argv[i]);
	}
	fprintf(stderr,"\n");
      }
    }else if(argc - optind < 1 ){
// handle the special case where NO files are specified on the command line. 
// This means read from stdin.  Since output_snapshot already understands a file with name '-' to
// mean read from stdin, we create an input file list with a single file named '-'.
      inputList = &stdinInput;
      numInputs = 1;
    }else{
// This is the normal case, where there are one or more files on the command line
//
      inputList = argv + optind;
      numInputs = argc - optind;
    }
    
  /* Start the main routines */  
  check_environment();
  output_snapshot(fragStoreName,bactigStoreName,gkpStoreName,inputList,numInputs,outputFileName,mapFileName,output,300,realUIDs,quiet,random,uidStart,argc,argv);
  
  }



  return (int)AS_TER_EXIT_SUCCESS;
} 
