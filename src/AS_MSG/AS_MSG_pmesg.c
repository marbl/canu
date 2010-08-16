
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
static char *rcsid= "$Id: AS_MSG_pmesg.c,v 1.52 2010-08-16 07:22:28 brianwalenz Exp $";

#include "AS_MSG_pmesg_internal.h"

AS_MSG_global_t *AS_MSG_globals = NULL;

char *
GetMemory(size_t nbytes) {

  //  Pad to a multiple of 8.
  if (nbytes % 8)
    nbytes += 8 - (nbytes % 8);

  assert((nbytes % 8) == 0);

  return((char *)GetHeapItems_AS(AS_MSG_globals->msgHeap, nbytes));
}


char *
ReadLine(FILE *fin, int skipComment) {

  //  Do until we get a complete non-comment line.
  do {
    errno = 0;

    AS_MSG_globals->curLineNum++;

    uint64   cloffset = 0;

    //  Do until we get a complete line.
    do {
      AS_MSG_globals->curLine[AS_MSG_globals->curLineMax - 2] = '\n';
      AS_MSG_globals->curLine[AS_MSG_globals->curLineMax - 1] = 0;

      //  Read as much of the line as possible into the current line buffer.
      //
      if (fgets(AS_MSG_globals->curLine + cloffset, AS_MSG_globals->curLineMax - cloffset - 1, fin) == NULL) {
        fprintf(stderr,"ERROR: AS_MSG_pmesg.c::ReadLine()-- Premature end of input at line " F_U64 " (%s)\n", AS_MSG_globals->curLineNum, AS_MSG_globals->msgCode);
        fprintf(stderr,"       '%s'\n", AS_MSG_globals->curLine);
        exit(1);
      }

      //  Detect other errors, exit.
      //
      if (errno) {
        fprintf(stderr,"ERROR: AS_MSG_pmesg.c::ReadLine()-- Read error at line " F_U64 ": '%s'\n", AS_MSG_globals->curLineNum, strerror(errno));
        fprintf(stderr,"       '%s'\n", AS_MSG_globals->curLine);
        exit(1);
      }

      //  If the line didn't fit completely, make the buffer larger.
      //
      if (AS_MSG_globals->curLine[AS_MSG_globals->curLineMax - 2] != '\n') {
#if 0
        fprintf(stderr, "WARNING: Input line "F_U64" is long (%s), resizing.\n", AS_MSG_globals->curLineNum, AS_MSG_globals->msgCode);
        fprintf(stderr, "         '%s'\n", AS_MSG_globals->curLine);
        fprintf(stderr, "         length = %d\n", strlen(AS_MSG_globals->curLine));
        fprintf(stderr, "         cloffset = %d\n", cloffset);
        fprintf(stderr, "         max = %d\n", AS_MSG_globals->curLineMax);
#endif
        cloffset = AS_MSG_globals->curLineMax - 2;

        AS_MSG_globals->curLineMax *= 2;
        AS_MSG_globals->curLine     = (char *)safe_realloc(AS_MSG_globals->curLine, sizeof(char) * AS_MSG_globals->curLineMax);

        AS_MSG_globals->curLine[AS_MSG_globals->curLineMax - 2] = 0;  //  Reset the end-of-line mark so we read more of the line
        AS_MSG_globals->curLine[AS_MSG_globals->curLineMax - 1] = 0;
      }
    } while (AS_MSG_globals->curLine[AS_MSG_globals->curLineMax - 2] != '\n');
  } while (skipComment && AS_MSG_globals->curLine[0] == '#');

  return(AS_MSG_globals->curLine);
}





// Field content area did not have correct syntax.
void
MfieldError(const char * const mesg) {
  int len;

  len = strlen(AS_MSG_globals->curLine)-1;
  if (AS_MSG_globals->curLine[len] == '\n')
    AS_MSG_globals->curLine[len] = 0;
  fprintf(stderr,"ERROR: %s '%s' (%s) at line " F_U64 "\n",
          mesg,AS_MSG_globals->curLine,AS_MSG_globals->msgCode,AS_MSG_globals->curLineNum);
  exit (1);
}





// Get a text field item: syntax "tag\n(%s\n)*.\n".
char *
GetText(const char * const tag, FILE *fin, const int delnewlines) {
  char *ret = AS_MSG_globals->msgBuffer + AS_MSG_globals->msgLen;

  ReadLine(fin, TRUE);

  if (strncmp(AS_MSG_globals->curLine,tag,4) != 0) {
    fprintf(stderr,"ERROR: Illegal tag '%s' (expected '%s') (%s) at line " F_U64 "\n",
            AS_MSG_globals->curLine, tag, AS_MSG_globals->msgCode, AS_MSG_globals->curLineNum);
    exit (1);
  }

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

  assert(AS_MSG_globals->msgLen <= AS_MSG_globals->msgMax);

  //if (AS_MSG_globals->curLineNum > 258000000)
  //  fprintf(stderr, "TEXT (len:"F_U64"): %s\n", AS_MSG_globals->msgLen, ret);

  return(ret);
}


// Get a string field item: syntax "tag%s\n".
char *
GetString(const char * const tag, FILE *fin) {
  char *ret = NULL;
  char *str = ReadLine(fin, TRUE);
  int   len = 0;

  if ((tag[0] != str[0]) ||
      (tag[1] != str[1]) ||
      (tag[2] != str[2]) ||
      (tag[3] != str[3])) {
    fprintf(stderr,"ERROR: Illegal tag '%s' (expected '%s') (%s) at line " F_U64 "\n",
            AS_MSG_globals->curLine, tag, AS_MSG_globals->msgCode, AS_MSG_globals->curLineNum);
    exit (1);
  }

  str += 4;
  len  = strlen(str);

  while (isspace(str[len-1])) {
    len -= 1;
    str[len] = 0;
  }

  ret = GetMemory(len + 1);
  memcpy(ret, str, len + 1);

  return(ret);
}


char
GetType(char *format, char *name, FILE *fin) {
  char value[2];
  ReadLine(fin, TRUE);
  if (sscanf(AS_MSG_globals->curLine, format, value) != 1) {
    fprintf(stderr,"ERROR: Illegal %s type value '%c' (%s) at line " F_U64 " \n",
            name,AS_MSG_globals->curLine[4],AS_MSG_globals->msgCode, AS_MSG_globals->curLineNum);
    exit(1);
  }
  return(value[0]);
}



void
GetEOM(FILE *fin) {
  if (ReadLine(fin,TRUE)[0] != '}')
    fprintf(stderr, "ERROR: ??? expecting end of message '}' at line "F_U64", got '%s' instead.\n",
            AS_MSG_globals->curLineNum, AS_MSG_globals->curLine), exit(1);
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

    AS_MSG_globals->msgMax     = 256;
    AS_MSG_globals->msgMax    *= 1024;
    AS_MSG_globals->msgMax    *= 1024;
    //AS_MSG_globals->msgMax    *= 1024;  //  8GB needed for GOSIII

    AS_MSG_globals->msgLen     = 0;
    AS_MSG_globals->msgBuffer  = (char *)safe_malloc(sizeof(char) * AS_MSG_globals->msgMax);

    AS_MSG_globals->msgHeap    = AllocateHeap_AS(1, 128 * 1024 * 1024);

    AS_MSG_globals->curLineMax = 1024;
    AS_MSG_globals->curLine    = (char *)safe_malloc(sizeof(char) * AS_MSG_globals->curLineMax);

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


uint64
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

  //  Our memory is now round-robin.  This is VERY important for VAR messages.  Terminator (before
  //  tigStore) would read an IMV, and attempt to write a VAR based on that memory.  Resetting the
  //  msgLen pointer to zero caused the IMV to get overwritten.
  //
  //  AFTER tigStore, we can switch back to resetting msgLen on reads.
  //
  AS_MSG_globals->msgLen = 0;

  ClearHeap_AS(AS_MSG_globals->msgHeap);

  //  Can't use ReadLine() here, because we want to return EOF if we
  //  read an empty line.
  //
  do {
    AS_MSG_globals->curLineNum++;
    if (fgets(AS_MSG_globals->curLine, AS_MSG_globals->curLineMax, fin) == NULL)
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
    fprintf(stderr,"ERROR: Unrecognized message type in '%s' at line " F_U64 "\n",
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

void
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

  //  After we have written the message, reset the buffer.
  AS_MSG_globals->msgLen = 0;

  ClearHeap_AS(AS_MSG_globals->msgHeap);
}
