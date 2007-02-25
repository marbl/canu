
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
static const char CM_ID[] = "$Id: AS_UTL_version.c,v 1.11 2007-02-25 08:13:38 brianwalenz Exp $";

#include "AS_UTL_version.h"

#define SAFE_FOPEN( adt_tmp_file, adt_tmp_name, the_mode) \
  assert(NULL == adt_tmp_file); adt_tmp_file = fopen(adt_tmp_name, the_mode); assert(NULL != adt_tmp_file);

#define  SAFE_FCLOSE(adt_tmp_file) \
  assert(NULL != adt_tmp_file); fclose(adt_tmp_file); adt_tmp_file = NULL;

VA_DEF(char);

int VersionStamp(int argc, char *argv[]) {
  time_t t;
  int i;
  int rc;
  
  char *command_name = NULL;

  command_name = (char *)safe_malloc(sizeof(char) * (strlen(argv[0])+15));

  fprintf(stderr,"\n+++++++++++++++++++++++++ VERSION INFO +++++++++++++++++++++++++\n");
  sprintf(command_name,"ident `which %s`",argv[0]);
  fprintf(stderr,"Version <%s>\n", command_name);
  rc=system(command_name);
  safe_free(command_name);
  
  fprintf(stderr,"\nComplete call: %s ",argv[0]);
  for (i=1;i<argc;i++) {
     fprintf(stderr,"%s ",argv[i]);
  }
  fprintf(stderr,"\n");
  t = time(0);
  fprintf(stderr,"Started: %s",ctime(&t));
  fprintf(stderr,"Working directory: ");
  system("pwd");
  fprintf(stderr,"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  return rc;
} 


int VersionStampADT(AuditMesg *adt_mesg, int argc, char *argv[]) {

  return VersionStampADTWithCommentAndVersion(adt_mesg, argc, argv,"","(blank)");
}

int VersionStampADTWithCommentAndVersion(AuditMesg *adt_mesg, int argc, char *argv[], char *comment, char *version) {

  time_t t0 = time(0);
  int i;
  int rc;
  char *command_name = NULL;
  char adt_tmp_name[] = "tmp_XXXXXX";
  int fd;
  FILE *adt_tmp_file = NULL;

  fd = mkstemp(adt_tmp_name);
  unlink(adt_tmp_name);

  command_name = (char *)safe_malloc(sizeof(char) * (strlen(argv[0])+20+FILENAME_MAX));
  sprintf(command_name,"ident `which %s` > %s",argv[0],adt_tmp_name);
  //  fprintf(stderr,"* Calling: %s\n", command_name);
  rc=system(command_name);

  adt_tmp_file = fdopen( fd, "a");
  assert(NULL != adt_tmp_file);

  fprintf(adt_tmp_file,"\nComplete call: %s ",argv[0]);
  for (i=1;i<argc;i++) {
    fprintf(adt_tmp_file,"%s ",argv[i]);
  }
  fprintf(adt_tmp_file,"\n");
  fprintf(adt_tmp_file,"Started: %s",ctime(&t0));
  fprintf(adt_tmp_file,"Working directory: ");
  SAFE_FCLOSE(adt_tmp_file);
  
  sprintf(command_name,"pwd >> %s",adt_tmp_name);
  //  fprintf(stderr,"* Calling: %s\n", command_name);
  rc = system(command_name);
  safe_free(command_name);
  // Now, read back into char *
  {
    int c;
    char ch;
    int ci = 0;
    //AuditLine *auditLine = (AuditLine *)malloc(sizeof(AuditLine));
    static AuditLine auditLine_memory;
    AuditLine *auditLine = &auditLine_memory;
    char *input = NULL;
    char *adt_char = NULL;
    char *startOfInterestingPart = NULL;
    size_t len;
    VA_TYPE(char) * adt_ident = CreateVA_char(5000);
    char * NoVersionWarning =
      "WARNING! No source code version keywords detected with ident!\n\n"
      "Please consider recompiling with different flags so that this\n"
      "useful information is not stripped from the program binary.\n\n";

    SAFE_FOPEN(adt_tmp_file, adt_tmp_name, "r");
    while( c = fgetc(adt_tmp_file), ch = c, c != EOF) {
      AppendVA_char(adt_ident,&ch);
      ci++;
    }
    // fprintf(stderr,"* Read in %d characters ... %d\n",    ci, GetNumchars(adt_ident));

    unlink(adt_tmp_name); 
    input = GetVA_char(adt_ident,0);
    //    fprintf(stderr,"* input = %s\n",input);
    // We're not interested in the first part of the ident output

    len = strlen(comment);
    startOfInterestingPart = strstr(input,"     $Id: ");
    // NOTE: gcc versions 3.4.0 to 3.4.2 (?) strip out ident strings
    //       and this ought to be checked anyway
    startOfInterestingPart = (startOfInterestingPart == NULL ?
                              NoVersionWarning : startOfInterestingPart);
    len += strlen(startOfInterestingPart);
    
    // fprintf(stderr,"* len = %d\n", len);
    adt_char = safe_malloc(sizeof(char) * (len+3));
    sprintf(adt_char,"\n%s\n", comment);
    strcat(adt_char,startOfInterestingPart);
    
    /*        fprintf(stderr,"* version:%s argc:%d argv:%s adt_char = %s\n", version, argc,argv[0],adt_char);
	      fflush(NULL); */
    AppendAuditLine_AS(adt_mesg, auditLine, t0, argv[0], version, adt_char);
    Delete_VA(adt_ident);
    // SAFE_FREE(adt_char); // Proto-IO needs this memory leak!!
    // free(auditLine);
    SAFE_FCLOSE(adt_tmp_file);
  }
  return rc;
} 
