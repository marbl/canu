
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
/* RCS info
 * $Id: AS_SIM_tobatches.c,v 1.1.1.1 2004-04-14 13:53:41 catmandew Exp $
 */


/********************************************************************** 

1. Overview

The Celera shotgun genome simulator needs to produce output the is
suitable for incremental assembly.  This routine reads a monolithic
fragment read message file and produces a family of files that
represent daily batches of fragment reads.

2. Memory Usage
Minimal.  

3. Interface

4. Design

This executable takes a fragment read message file in the Celera
Assembler's Prototype Assembler Message Routine format and produces a
sequence of smaller files with the same content.

This filter counts the number of messages it has read. When a
threshold number of message is encountered, then the output starts to a
a new file in the sequence.

5. Limitations

There must be a three letter file name extenstion.  This file filter
does not discriminate about message type when chopping the file. It
just counts the number of message that it has read.

6. Status

This executable is incorporated into the celsim AWK script.

**********************************************************************/

/********************************************************************** 

Module:

tobatch [-P] [-b <batch size>] filename.ext

Description: 

This executable takes a proto-msg file and chops it into batches
of about <batch size> in number of records. 
A batch may be smaller if there is not enough data.

-P : Specifies ASCII output. The default is binary.
-b <batch size>: Specifies the desired batch size.  The default batch
size is 200,000 records.

filename: This is the file name to be read.

Assumptions: 

The input file must use the protomesg package.

**************************************************************************/

/*************************************************************************/
/* System include files */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <unistd.h> /* man 3 getopt */
/*************************************************************************/
/* Local include files */
#include "AS_global.h"
#include "AS_MSG_pmesg.h"

/*************************************************************************/
/* Conditional compilation */
/*************************************************************************/
/* Constant definitions; Macro definitions; type definitions */
#define DEFAULT_BATCH_SIZE 200000
/* This is the expected daily batch size expected in April 1999. */
/*************************************************************************/
/* Static Globals */

/*************************************************************************/
/* Function prototypes for internal static functions */
/*************************************************************************/
/* Utility functions */
void outputBAT(FILE *fout, int binary_output_mode, char *name,
               char *comment, CDS_UID_t eid)
{
  
  GenericMesg outMesg;
  BatchMesg batchMesg;

  fprintf(stderr,"* OutputBAT %s %s " F_UID "\n", name, comment, eid);
  batchMesg.name = name;
  batchMesg.comment = comment;
  batchMesg.created = CREATION_TIME;
  batchMesg.eaccession = eid;

  outMesg.m = &batchMesg;
  outMesg.t = MESG_BAT;

  if(binary_output_mode) {
    WriteBinaryMesg_AS(fout,&outMesg);
  } else {
    WriteProtoMesg_AS(fout, &outMesg);
  }


}

/*************************************************************************/
/*************************************************************************/

int main(int argc, char *argv[])
{
  char *File_Prefix = NULL,*File_Suffix = NULL,File_Name[80],File_Name_Tmp[80];
  char strtmp[1024];
  int ifile,imsg,notdone;
  int batch_size = DEFAULT_BATCH_SIZE;
  char strlabel[100];
  int binary_output_mode=TRUE;
  /* int32 seconds from sometime in 1970 */
  time_t tp1,tp2;

  MesgReader   reader;
  MessageType    imesgtype;
  GenericMesg   *pmesg;

  FILE *fin,*fout;
  
   { /* Parse the argument list using "man 3 getopt". */ 
     int ch,errflg=0,illegal = 0;
     optarg = NULL;
     while (!errflg && ((ch = getopt(argc, argv, "Pb:")) != EOF))
       switch(ch) {
       case 'P':
	 binary_output_mode=FALSE;
	 break;
       case 'b':
	 batch_size = atoi(optarg);
	 assert(batch_size>0);
	 break;
       case '?':
	 fprintf(stderr,"Unrecognized option -%c",optopt);
       default :
	 errflg++;
       }

     if(optind < argc) {
       strcpy(File_Name,argv[optind++]);
       fprintf(stderr,"File_Name = %s\n",File_Name);
       strcpy(File_Name_Tmp,File_Name);
       File_Prefix = File_Name_Tmp;
       File_Suffix = strrchr(File_Name_Tmp,'.');
       *File_Suffix = '\0'; File_Suffix++;
       fprintf(stderr,"File_Prefix = %s\n",File_Prefix);
       fprintf(stderr,"File_Suffix = %s\n",File_Suffix);
     }
     if(illegal == 1)
       {
	 fprintf (stderr, 
		  "USAGE: %s [-P] [-b <batch size>] <dataset_name>\n", 
		  argv [0]);
	 exit (EXIT_FAILURE);
       }
     /* End of command line parsing */
   }
   
  fprintf(stderr,__FILE__ " "  __DATE__ " " __TIME__ "\n");
  fprintf(stderr,"$Id: AS_SIM_tobatches.c,v 1.1.1.1 2004-04-14 13:53:41 catmandew Exp $\n");
  fprintf(stderr,"Batch size = %d\n",batch_size);

  fin = fopen(File_Name,"r");
  assert(fin != NULL);
  reader = InputFileType_AS(fin);

  notdone = TRUE;
  fprintf(stderr,"zoot\n");
  for(ifile=1;(ifile<10000)&&notdone;ifile++){
    int okay_to_chop_here =TRUE;
    strcpy(strtmp,File_Prefix);
    sprintf(strlabel,"_%05d.",ifile); /* This width of 4 constrains us 
					to 10000 files per dataset batch.*/
    strcat(strtmp,strlabel);
    strcat(strtmp,File_Suffix);
    fout = fopen(strtmp,"w");
    if(ifile > 1)
    {
      fprintf(stderr,"* Printing outputBAT\n");
      outputBAT(fout, binary_output_mode, strlabel, "(created by tobatches)",
#if 0
		(1L<<63) - ifile
#else
		(((CDS_UID_t)1)<<63) - ifile
#endif
		);
    }
    time(&tp1); fprintf(stderr,"Begin writing %s\n",strtmp);
    fflush(NULL);
    
    for(imsg=0; (imsg<batch_size); imsg++) {
      if((imsg >= batch_size ) && (okay_to_chop_here)) break;

      notdone = (EOF != reader(fin, &pmesg));
      //fprintf(stderr,"ifile=%d,imsg=%d\n",ifile,imsg);
      /* if(!notdone)then(okay_to_chop_here) */
      assert(notdone || okay_to_chop_here);
      if(!notdone) break;

      imesgtype = pmesg->t;
      if(imesgtype != EOF) {
	if(binary_output_mode) {
	  WriteBinaryMesg_AS(fout,pmesg);
	} else {
	  WriteProtoMesg_AS(fout,pmesg);
	  }
	fflush(fout);
      }
    }
    fclose(fout);
    time(&tp2); fprintf(stderr,"Finished %s, walltime=%10d sec\n",strtmp,
		       (int)(tp2-tp1));
    if(!notdone) break;
  }
  fclose(fin);
  return 0;
}
