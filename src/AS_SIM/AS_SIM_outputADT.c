
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
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> /* man 3 getopt */
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"

#define NSTRLEN 1024 // Maximum string length quantum.

static int trim_eol(char * const tool,const size_t nstrlen)
{
  int ii;
  if(NULL == tool) { return FALSE;} 
  for(ii=0;ii<nstrlen;ii++) {
    char ch = tool[ii];
    if(
       (ch == '\0') ||
       (ch == '\n') ||
       (ch == '\r')
       ) {
      tool[ii] = '\0';
      return TRUE;
    }
  }
  return FALSE;
}

void outputADT(char *tool, char *version, char *commentString)
{
  
  GenericMesg outMesg;
  AuditMesg auditMesg;
  AuditLine auditLine;

  fprintf(stderr,"tool = <%s>  version = <%s>\n", tool, version);
  auditLine.complete = CREATION_TIME;
  auditLine.name = tool;
  auditLine.version = version;
  auditLine.comment = commentString;
  auditLine.next = NULL;

  auditMesg.list = &auditLine;
  outMesg.m = &auditMesg;
  outMesg.t = MESG_ADT;

  WriteProtoMesg_AS(stdout,&outMesg);
}

static void outputBAT
(char * name,
 char * comment,
 CDS_UID_t  eid)
{
  
  GenericMesg outMesg;
  BatchMesg batchMesg;

  batchMesg.name = name;
  batchMesg.comment = comment;
  batchMesg.created = CREATION_TIME;
  batchMesg.eaccession = eid;

  outMesg.m = &batchMesg;
  outMesg.t = MESG_BAT;

  WriteProtoMesg_AS(stdout,&outMesg);
}


/* >>>> MAIN / TOP <<<< */

int main(int argc, char *argv[])
{
  // UID of first dros read in first input file...
  char eid_string[NSTRLEN] = "";
  char tool[NSTRLEN] = "";
  char genomeLength[NSTRLEN] = "";
  char commentString[256 * NSTRLEN] = "";
  int c;
  char *p = commentString;
  CDS_UID_t eid = FIRST_UID;
  char * ret = NULL;
  int iret = 0;

  /**************** Process Command Line Arguments *********************/
  { /* Parse the argument list using "man 3 getopt". */ 
    int ch;
    optarg = NULL;
    while (
	   ((ch = getopt(argc, argv, 
			 "b:c:l:"
                         )) != EOF)) {
      switch(ch) {
        /* The required command line options: */
      case 'b':
	strcpy(eid_string,optarg);
	break;
      case 'l':
	strcpy(genomeLength,optarg);
	break;
      case 'c':
	strcpy(commentString,optarg);
	break;
      case '?':
      default :
	fprintf(stderr,"Unrecognized option -%c\n",optopt);
        assert(0);
      }
    }
  }

  if(eid_string[0] != '\0')
    eid = STR_TO_UID(eid_string, NULL, 10);
  ret = fgets(tool,NSTRLEN,stdin);
  assert(NULL != ret);
  iret = trim_eol(tool,NSTRLEN);
  
  ret = fgets(genomeLength,NSTRLEN,stdin);
  assert(NULL != ret);
  iret = trim_eol(genomeLength,NSTRLEN);

  // ret = fgets(commentString,NSTRLEN,stdin);
  // assert(NULL != ret);
  // iret = trim_eol(tool,NSTRLEN);
  while((c = getchar()) != EOF){
    *p++ = c;
  }
  *p = '\0';

  outputBAT("celsim output","(No comment)", eid);
  outputADT(tool, genomeLength, commentString);

  exit(0);
}

