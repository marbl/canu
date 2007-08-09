
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
static char CM_ID[]= "$Id: AS_MSG_pmesg.c,v 1.41 2007-08-09 16:55:34 brianwalenz Exp $";

#include "AS_MSG_pmesg_internal.h"

AS_MSG_global_t *AS_MSG_globals = NULL;



// Make sure there is a block of size bytes left in memory buffer
// starting at an index that is a multiple of boundary.  Return the
// *index* into the array, so that the realloc does not blow
// structures in the process of being built.  All pointers in such
// structures are saved as integers and converted to pointers after
// all allocation has taken place.
//
void
MakeSpace(const int size){
  size_t newsize=1;
  char *newbufr;

  if(AS_MSG_globals->MemMax> size)
    return;

  if(AS_MSG_globals->MemMax < 0)
    newsize = 2 * 2048 * 2048; // This may be excessive, but this code is BRITTLE!
  else
    newsize = 2 * size;

  newsize = ROUNDUP(newsize,8);
  newbufr = (char *)safe_realloc(AS_MSG_globals->MemBuffer,newsize);
  AS_MSG_globals->MemBuffer = newbufr;
  AS_MSG_globals->MemMax    = newsize;
}


long
MoreSpace(const int size, const int boundary) { 
  AS_MSG_globals->MemTop = ROUNDUP(AS_MSG_globals->MemTop,boundary);

  MakeSpace(AS_MSG_globals->MemTop + size);
  { int alloc;

    alloc   = AS_MSG_globals->MemTop;
    AS_MSG_globals->MemTop += size;
    return (alloc);
  }
}


char *
ReadLine(FILE *fin) {

  AS_MSG_globals->CurLine[MAX_LINE_LEN-2] = '\n';
  AS_MSG_globals->CurLine[MAX_LINE_LEN-1] = 0;

  errno = 0;

  AS_MSG_globals->LineNum++;

#ifdef FGETS_IS_BROKEN
  int p=0;

  for (p=0; p<MAX_LINE_LEN; p++)
    AS_MSG_globals->CurLine[p] = 0;
  AS_MSG_globals->CurLine[MAX_LINE_LEN-2] = '\n';
  AS_MSG_globals->CurLine[MAX_LINE_LEN-1] = 0;

  for (p=0; p<MAX_LINE_LEN-1; p++) {
    AS_MSG_globals->CurLine[p]   = fgetc(fin);
    AS_MSG_globals->CurLine[p+1] = 0;
    if (AS_MSG_globals->CurLine[p] == 0)
      AS_MSG_globals->CurLine[p] = '\n';
    if (AS_MSG_globals->CurLine[p] == '\n') {
      AS_MSG_globals->CurLine[p+1] = 0;
      break;
    }
  }

  //  Whatever.  We can print the line, or we can iterate it.  If we
  //  don't, we die.  OK, we can't.  Gotta print.
  //
  //AS_MSG_globals->CurLine[MAX_LINE_LEN-2] = '\n';
  //AS_MSG_globals->CurLine[MAX_LINE_LEN-1] = 0;
  //
  //for (p=0; AS_MSG_globals->CurLine[p]; p++)
  //  ;
  //fprintf(stdout, "%*s", AS_MSG_globals->CurLine);

#else
  if (fgets(AS_MSG_globals->CurLine, MAX_LINE_LEN-1, fin) == NULL) {
    fprintf(stderr,"ERROR: AS_MSG_pmesg.c::ReadLine()-- Premature end of input at line %d (%s)\n", AS_MSG_globals->LineNum, AS_MSG_globals->Mcode);
    fprintf(stderr,"       '%s'\n", AS_MSG_globals->CurLine);
    exit(1);
  }
#endif

  if (errno) {
    fprintf(stderr,"ERROR: AS_MSG_pmesg.c::ReadLine()-- Read error at line %d: '%s'\n", AS_MSG_globals->LineNum, strerror(errno));
    fprintf(stderr,"       '%s'\n", AS_MSG_globals->CurLine);
    exit(1);
  }

  if (AS_MSG_globals->CurLine[MAX_LINE_LEN-2] != '\n') {
    fprintf(stderr,"ERROR: Input line %d is too long (%s)\n", AS_MSG_globals->LineNum, AS_MSG_globals->Mcode);
    fprintf(stderr,"       '%s'\n", AS_MSG_globals->CurLine);
    exit(1);
  }

  //fprintf(stderr, "READLINE -- %d %d %s", AS_MSG_globals->CurLine[0], AS_MSG_globals->CurLine[1], AS_MSG_globals->CurLine);

  return(AS_MSG_globals->CurLine);
}


char *
GetLine(FILE *fin, int skipComment) {
  /* Get next input line (there must be one). */

  do {
    ReadLine(fin);
  } while (skipComment && AS_MSG_globals->CurLine[0] == '#');
  return (AS_MSG_globals->CurLine);
}

/* Found an enum out-of-range error. */

void
MtypeError(const char * const name) {
  fprintf(stderr,"ERROR: Illegal %s type value \"%c\" (%s) at line %d\n",
          name,AS_MSG_globals->CurLine[4],AS_MSG_globals->Mcode, AS_MSG_globals->LineNum);
  exit (1);
}

/* Found a bad 3-code field name. */

void
MtagError(const char * const tag) {
  fprintf(stderr,"ERROR: Illegal tag \"%s\" (expected \"%s\") (%s) at line %d\n",
          AS_MSG_globals->CurLine, tag, AS_MSG_globals->Mcode, AS_MSG_globals->LineNum);
  exit (1);
}

/* Field content area did not have correct syntax. */

void
MfieldError(const char * const mesg) {
  int len;

  len = strlen(AS_MSG_globals->CurLine)-1;
  if (AS_MSG_globals->CurLine[len] == '\n')
    AS_MSG_globals->CurLine[len] = 0;
  fprintf(stderr,"ERROR: %s \"%s\" (%s) at line %d\n",
          mesg,AS_MSG_globals->CurLine,AS_MSG_globals->Mcode,AS_MSG_globals->LineNum);
  exit (1);
}

/* General error message exit. */

void
MgenError(const char * const mesg) {
  fprintf(stderr,"ERROR: %s (%s) at line %d\n",
          mesg,AS_MSG_globals->Mcode,AS_MSG_globals->LineNum);
  exit (1);
}

/* Get a text field item: syntax "tag\n(%s\n)*.\n".
   Tricky part is length is not known ahead of time.  */

long
GetText(const char * const tag, FILE *fin, const int delnewlines) {
  long text, idx; 
  int len;
 
  if (strncmp(GetLine(fin, TRUE),tag,4) != 0)
    MtagError(tag);
  text = AS_MSG_globals->MemTop;
  while (1)
    { 
      char *line = GetLine(fin, FALSE);
      if ((line[0] == '.') && (line[1] == '\n'))
        break;
      len = strlen(AS_MSG_globals->CurLine);
      if (delnewlines && AS_MSG_globals->CurLine[len-1] == '\n') len -= 1;
      idx = MoreSpace(len,1);
      strncpy(AS_MSG_globals->MemBuffer+idx,AS_MSG_globals->CurLine,len);
    }
  idx = MoreSpace(1,1);
  AS_MSG_globals->MemBuffer[idx] = 0;
#if 0
  while ((text != AS_MSG_globals->MemBuffer+idx) && isspace(AS_MSG_globals->MemBuffer[idx-1])) {
    idx--;
    AS_MSG_globals->MemBuffer[idx] = 0;
  }
#endif
  return (text);
}

/* Get a string field item: syntax "tag%s\n".
   Tricky part is length is not known ahead of time.  */

long
GetString(const char * const tag, FILE *fin) {
  char *str;
  int   eos, len, text, idx;

  errno = 0;

  text = AS_MSG_globals->MemTop;
  do
    {
      ReadLine(fin);
    }
  while (AS_MSG_globals->CurLine[0] == '#');
  if (strncmp(AS_MSG_globals->CurLine,tag,4) != 0)
    MtagError(tag);
  str = AS_MSG_globals->CurLine + 4;

  while (1)
    { len = strlen(str);
      eos = (str[len-1] == '\n');
      if (eos) len -= 1;
      idx = MoreSpace(len,1);
      strncpy(AS_MSG_globals->MemBuffer+idx,str,len);
      if (eos) break;
      str = ReadLine(fin);
    }
  idx = MoreSpace(1,1);
  AS_MSG_globals->MemBuffer[idx] = 0;
  return (text);
}

/* Output text field item with 3-code field-name "tag". */

void
PutText(FILE *fout, const char * const tag, 
        char * text, const int format) {
  // Note that the data of "text" is modified!!!
  int i, len;
  fprintf(fout,"%s\n",tag);
  if ((text != NULL) && (text[0] != 0)) {
    len = strlen(text);
    if (format) { 
      for (i = 0; i < len; i += 70)
        { 
          fprintf(fout,"%.*s\n",70,text);
          text += 70;
        }
      fprintf(fout,".\n");
    } else{ 
      if (text[len-1] == '\n')     /* Strip trailing new line if prez. */
        text[len-1] = 0;
      fprintf(fout,"%s\n.\n", text);
    }
  } else{
    fprintf(fout,".\n");
  }
}






void
AppendAuditLine_AS(AuditMesg *adt, AuditLine *adl,
                   time_t t, char *name, char *version, char *comment) {
  AuditLine *list;

  // Skip over the entire list, stopping at the last element
  if(adt->list) {
    for (list = adt->list; list->next != NULL; list = list->next);

    list->next = adl;
  }else{
    adt->list = adl;
  }

  adl->complete = t;
  adl->name     = name;
  adl->version  = version;
  adl->comment  = comment;
  adl->next     = NULL;
}




static
void
AS_MSG_globalsInitialize(void) {
  if (AS_MSG_globals == NULL) {
    AS_MSG_globals = (AS_MSG_global_t *)safe_calloc(1, sizeof(AS_MSG_global_t));
    AS_MSG_globals->MemMax        = -1;
    AS_MSG_globals->MemTop        = -1;
    AS_MSG_setFormatVersion(1);
  }
}

void
AS_MSG_setFormatVersion(int format) {
  AS_MSG_globalsInitialize();
  switch (format) {
    case 1:
      AS_MSG_setFormatVersion1();
      break;
    case 2:
      AS_MSG_setFormatVersion2();
      break;
    default:
      fprintf(stderr, "AS_MSG_SetFormatVersion()--  Unknown format %d/\n", format);
      exit(1);
      break;
  }
}


int
GetMessageType(char *string){
  int t;
  AS_MSG_globalsInitialize();
  for(t = 1; t < NUM_OF_REC_TYPES;t++){
    if ((AS_MSG_globals->CallTable[t].header) &&
        (strncmp(string, AS_MSG_globals->CallTable[t].header + 1, 3) == 0))
      return t;
  }
  return 0;
}


const char *
GetMessageName(int type){
  AS_MSG_globalsInitialize();
  if (type >= NUM_OF_REC_TYPES || type < 1)
    return(NULL);
  if (AS_MSG_globals->CallTable[type].header == NULL)
    return(NULL);
  return(AS_MSG_globals->CallTable[type].header + 1);
}


int
GetProtoLineNum_AS(void) {
  AS_MSG_globalsInitialize();
  return (AS_MSG_globals->LineNum);
}



int
ReadProtoMesg_AS(FILE *fin, GenericMesg **pmesg) {
  int t;

  AS_MSG_globalsInitialize();

  *pmesg = &AS_MSG_globals->ReadMesg;

  errno = 0;

  AS_MSG_globals->CurLine[MAX_LINE_LEN-2] = '\n';
  AS_MSG_globals->MemTop    = 0;
  do {
    //  Can't use ReadLine() here, because we want to return EOF if we
    //  read an empty line.
    AS_MSG_globals->LineNum++;
    if (fgets(AS_MSG_globals->CurLine,MAX_LINE_LEN,fin) == NULL)
      return (EOF);
  } while (AS_MSG_globals->CurLine[0] == '#');

  if (errno)
    fprintf(stderr, "ERROR: Read Failure looking for message type: %s\n",
            strerror(errno)), exit(1);

  for(t = 1; t <= NUM_OF_REC_TYPES; t++)
    if ((AS_MSG_globals->CallTable[t].header) &&
        (strncmp(AS_MSG_globals->CurLine, AS_MSG_globals->CallTable[t].header, 4) == 0))
      break;
  if (t > NUM_OF_REC_TYPES) {
    int len = strlen(AS_MSG_globals->CurLine)-1;
    if (AS_MSG_globals->CurLine[len] == '\n')
      AS_MSG_globals->CurLine[len] = 0;
    fprintf(stderr,"ERROR: Unrecognized message type in \"%s\" at line %d\n",
            AS_MSG_globals->CurLine,AS_MSG_globals->LineNum);
    exit (1);
  }
  AS_MSG_globals->Mcode = AS_MSG_globals->CallTable[t].header+1;

  (*pmesg)->t = (MessageType) t;
  (*pmesg)->m = AS_MSG_globals->CallTable[t].reader(fin);
  (*pmesg)->s = AS_MSG_globals->MemTop;

  if (errno) {
    fprintf(stderr, "ERROR: Read Failure reading message %s: %s\n",
            AS_MSG_globals->CallTable[t].header+1, strerror(errno));
    exit(1);
  }

  return (0);
}

int
WriteProtoMesg_AS(FILE *fout, GenericMesg *pmesg) {

  AS_MSG_globalsInitialize();

  errno = 0;
  AS_MSG_globals->CallTable[pmesg->t].writer(fout,pmesg->m);
  if (errno) {
    fprintf(stderr, "ERROR: Write Failure: %s\n", strerror(errno));
    exit(1);
  }
  return (0);
}
