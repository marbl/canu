
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
static char CM_ID[]= "$Id: AS_MSG_pmesg.c,v 1.45 2008-06-27 06:29:17 brianwalenz Exp $";

#include "AS_MSG_pmesg_internal.h"

AS_MSG_global_t *AS_MSG_globals = NULL;



char *
GetMemory(size_t nbytes) {
  if (AS_MSG_globals->msgLen % 8)
    AS_MSG_globals->msgLen += 8 - (AS_MSG_globals->msgLen % 8);
  char *ret = AS_MSG_globals->msgBuffer + AS_MSG_globals->msgLen;
  AS_MSG_globals->msgLen += nbytes;
  return(ret);
}



char *
ReadLine(FILE *fin, int skipComment) {

  do {
    AS_MSG_globals->curLine[MAX_LINE_LEN-2] = '\n';
    AS_MSG_globals->curLine[MAX_LINE_LEN-1] = 0;

    errno = 0;

    AS_MSG_globals->curLineNum++;

#ifdef FGETS_IS_BROKEN
    int p=0;

    for (p=0; p<MAX_LINE_LEN; p++)
      AS_MSG_globals->curLine[p] = 0;
    AS_MSG_globals->curLine[MAX_LINE_LEN-2] = '\n';
    AS_MSG_globals->curLine[MAX_LINE_LEN-1] = 0;

    for (p=0; p<MAX_LINE_LEN-1; p++) {
      AS_MSG_globals->curLine[p]   = fgetc(fin);
      AS_MSG_globals->curLine[p+1] = 0;
      if (AS_MSG_globals->curLine[p] == 0)
        AS_MSG_globals->curLine[p] = '\n';
      if (AS_MSG_globals->curLine[p] == '\n') {
        AS_MSG_globals->curLine[p+1] = 0;
        break;
      }
    }
#else
    if (fgets(AS_MSG_globals->curLine, MAX_LINE_LEN-1, fin) == NULL) {
      fprintf(stderr,"ERROR: AS_MSG_pmesg.c::ReadLine()-- Premature end of input at line %d (%s)\n", AS_MSG_globals->curLineNum, AS_MSG_globals->msgCode);
      fprintf(stderr,"       '%s'\n", AS_MSG_globals->curLine);
      exit(1);
    }
#endif

    if (errno) {
      fprintf(stderr,"ERROR: AS_MSG_pmesg.c::ReadLine()-- Read error at line %d: '%s'\n", AS_MSG_globals->curLineNum, strerror(errno));
      fprintf(stderr,"       '%s'\n", AS_MSG_globals->curLine);
      exit(1);
    }

    if (AS_MSG_globals->curLine[MAX_LINE_LEN-2] != '\n') {
      fprintf(stderr,"ERROR: Input line %d is too long (%s)\n", AS_MSG_globals->curLineNum, AS_MSG_globals->msgCode);
      fprintf(stderr,"       '%s'\n", AS_MSG_globals->curLine);
      exit(1);
    }
  } while (skipComment && AS_MSG_globals->curLine[0] == '#');

  //fprintf(stderr, "ReadLine()-- '%s", AS_MSG_globals->curLine);

  return(AS_MSG_globals->curLine);
}





// Found an enum out-of-range error.
void
MtypeError(const char * const name) {
  fprintf(stderr,"ERROR: Illegal %s type value '%c' (%s) at line %d\n",
          name,AS_MSG_globals->curLine[4],AS_MSG_globals->msgCode, AS_MSG_globals->curLineNum);
  exit (1);
}

// Found a bad 3-code field name.
void
MtagError(const char * const tag) {
  fprintf(stderr,"ERROR: Illegal tag '%s' (expected '%s') (%s) at line %d\n",
          AS_MSG_globals->curLine, tag, AS_MSG_globals->msgCode, AS_MSG_globals->curLineNum);
  exit (1);
}

// Field content area did not have correct syntax.
void
MfieldError(const char * const mesg) {
  int len;

  len = strlen(AS_MSG_globals->curLine)-1;
  if (AS_MSG_globals->curLine[len] == '\n')
    AS_MSG_globals->curLine[len] = 0;
  fprintf(stderr,"ERROR: %s '%s' (%s) at line %d\n",
          mesg,AS_MSG_globals->curLine,AS_MSG_globals->msgCode,AS_MSG_globals->curLineNum);
  exit (1);
}

// General error message exit.
void
MgenError(const char * const mesg) {
  fprintf(stderr,"ERROR: %s (%s) at line %d\n",
          mesg,AS_MSG_globals->msgCode,AS_MSG_globals->curLineNum);
  exit (1);
}







// Get a text field item: syntax "tag\n(%s\n)*.\n".
char *
GetText(const char * const tag, FILE *fin, const int delnewlines) {
  char *ret = AS_MSG_globals->msgBuffer + AS_MSG_globals->msgLen;

  ReadLine(fin, TRUE);

  if (strncmp(AS_MSG_globals->curLine,tag,4) != 0)
    MtagError(tag);

  while (1) {
    ReadLine(fin, FALSE);

    if ((AS_MSG_globals->curLine[0] == '.') &&
        (AS_MSG_globals->curLine[1] == '\n'))
      break;

    int len = strlen(AS_MSG_globals->curLine);

    if (delnewlines && AS_MSG_globals->curLine[len-1] == '\n') {
      len -= 1;
      AS_MSG_globals->curLine[len] = 0;
    }

    memcpy(AS_MSG_globals->msgBuffer + AS_MSG_globals->msgLen, AS_MSG_globals->curLine, len);
    AS_MSG_globals->msgLen += len;
  }

  AS_MSG_globals->msgBuffer[AS_MSG_globals->msgLen++] = 0;

  return(ret);
}


// Get a string field item: syntax "tag%s\n".
char *
GetString(const char * const tag, FILE *fin) {
  char *ret = AS_MSG_globals->msgBuffer + AS_MSG_globals->msgLen;
  char *str = ReadLine(fin, TRUE);
  int   len = 0;

  if ((tag[0] != str[0]) ||
      (tag[1] != str[1]) ||
      (tag[2] != str[2]) ||
      (tag[3] != str[3]))
    MtagError(tag);

  str += 4;
  len  = strlen(str);

  while (isspace(str[len-1])) {
    len -= 1;
    str[len] = 0;
  }

  memcpy(AS_MSG_globals->msgBuffer + AS_MSG_globals->msgLen, str, len);
  AS_MSG_globals->msgLen += len;

  AS_MSG_globals->msgBuffer[AS_MSG_globals->msgLen++] = 0;

  return(ret);
}


char
GetType(char *format, char *name, FILE *fin) {
  char value[2];
  ReadLine(fin, TRUE);
  if (sscanf(AS_MSG_globals->curLine, format, value) != 1)
    MtypeError(name);
  return(value[0]);
}



void
GetEOM(FILE *fin) {
  if (ReadLine(fin,TRUE)[0] != '}')
    MgenError("Expecting end of message");
}


AS_UID
GetUID(char *tag, FILE *fin) {
  return(AS_UID_load(GetString(tag, fin)));
}


AS_UID
GetUIDIID(char *tag, AS_IID *iid, FILE *fin) {
  char   *uidstr = GetString(tag, fin) + 1;
  char   *iidstr = uidstr;

  //  Skip over the initial '(' (that's the +1 above), then move
  //  iidstr to the ',' separator and make it a null terminating byte,
  //  then move iidstr to the beginning of the iid.

  while (*iidstr != ',')
    iidstr++;
  *iidstr = 0;
  iidstr++;

  *iid = strtoul(iidstr, NULL, 10);

  return(AS_UID_load(uidstr));
}





// Output text field item with 3-code field-name "tag".
//
// Note that the data of "text" is modified!!!
void
PutText(FILE *fout,
        const char * const tag,
        char * text,
        const int format) {

  int i, len;

  fprintf(fout,"%s\n",tag);

  if ((text != NULL) && (text[0] != 0)) {
    len = strlen(text);

    if (format) {
      for (i = 0; i < len; i += 70) {
        fprintf(fout,"%.*s\n",70,text);
        text += 70;
      }

    } else{
      // Strip trailing new line if prez.
      if (text[len-1] == '\n')
        text[len-1] = 0;
      fprintf(fout,"%s\n", text);
    }
  }
  fprintf(fout,".\n");
}






static
void
AS_MSG_globalsInitialize(void) {
  if (AS_MSG_globals == NULL) {
    AS_MSG_globals = (AS_MSG_global_t *)safe_calloc(1, sizeof(AS_MSG_global_t));

    AS_MSG_globals->msgMax    = MAX_MESG_LEN;
    AS_MSG_globals->msgLen    = 0;
    AS_MSG_globals->msgBuffer = (char *)safe_malloc(sizeof(char) * AS_MSG_globals->msgMax);

    AS_MSG_globals->curLine   = (char *)safe_malloc(sizeof(char) * MAX_LINE_LEN);

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
  return (AS_MSG_globals->curLineNum);
}



int
ReadProtoMesg_AS(FILE *fin, GenericMesg **pmesg) {
  int t;

  AS_MSG_globalsInitialize();

  *pmesg = &AS_MSG_globals->readMesg;

  errno = 0;

  AS_MSG_globals->msgLen = 0;

  //  Can't use ReadLine() here, because we want to return EOF if we
  //  read an empty line.
  //
  do {
    AS_MSG_globals->curLineNum++;
    if (fgets(AS_MSG_globals->curLine,MAX_LINE_LEN,fin) == NULL)
      return (EOF);

    //  Brute force skip ADT messages.
    //
    //  Some legacy pipelines still produce these, but we don't care.  Format is:
    //
    //  {ADT
    //  {ADL
    //  who:name
    //  ctm:1206563351
    //  vsn:1.00
    //  com:
    //  comments
    //  .
    //  }
    //  {ADL
    //  ...
    //  }
    //  {ADL
    //  ...
    //  }
    //  .
    //  }
    //
    if ((AS_MSG_globals->curLine[0] == '{') &&
        (AS_MSG_globals->curLine[1] == 'A') &&
        (AS_MSG_globals->curLine[2] == 'D') &&
        (AS_MSG_globals->curLine[3] == 'T')) {
      while (strncmp(ReadLine(fin,TRUE),"{ADL",4) == 0) {
        while (strncmp(ReadLine(fin,TRUE),".",1) != 0)
          //  Do nothing.
          ;
        GetEOM(fin);
      }
      GetEOM(fin);

      //  Force us to read another message
      AS_MSG_globals->curLine[0] = '#';
    }
  } while (AS_MSG_globals->curLine[0] == '#');

  if (errno)
    fprintf(stderr, "ERROR: Read Failure looking for message type: %s\n",
            strerror(errno)), exit(1);

  for (t=1; t<=NUM_OF_REC_TYPES; t++)
    if ((AS_MSG_globals->CallTable[t].header) &&
        (AS_MSG_globals->curLine[0] == AS_MSG_globals->CallTable[t].header[0]) &&
        (AS_MSG_globals->curLine[1] == AS_MSG_globals->CallTable[t].header[1]) &&
        (AS_MSG_globals->curLine[2] == AS_MSG_globals->CallTable[t].header[2]) &&
        (AS_MSG_globals->curLine[3] == AS_MSG_globals->CallTable[t].header[3]))
      break;

  if (t > NUM_OF_REC_TYPES) {
    int len = strlen(AS_MSG_globals->curLine)-1;
    if (AS_MSG_globals->curLine[len] == '\n')
      AS_MSG_globals->curLine[len] = 0;
    fprintf(stderr,"ERROR: Unrecognized message type in '%s' at line %d\n",
            AS_MSG_globals->curLine, AS_MSG_globals->curLineNum);
    exit (1);
  }

  AS_MSG_globals->msgCode = AS_MSG_globals->CallTable[t].header + 1;

  (*pmesg)->t = (MessageType)t;
  (*pmesg)->m = AS_MSG_globals->CallTable[t].reader(fin);

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

  //  This is a quirk with printf; redirecting output to /dev/null
  //  generates ENOTTY.

  errno = 0;
  AS_MSG_globals->CallTable[pmesg->t].writer(fout,pmesg->m);
  if ((errno) && (errno != ENOTTY)) {
    fprintf(stderr, "ERROR: Write Failure: %s\n", strerror(errno));
    exit(1);
  }
  return (0);
}
