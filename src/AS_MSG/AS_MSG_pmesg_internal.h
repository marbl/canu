
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

#ifndef AS_MSG_PMESG_INTERNAL_H
#define AS_MSG_PMESG_INTERNAL_H

//  FreeBSD 6.1 fgets() sporadically replaces \n with \0, which
//  horribly breaks this reader.  Defined this to replace
//  fgets() with fgetc().
//#define FGETS_IS_BROKEN

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <ctype.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"


#define ROUNDUP(n,u) ((((n)-1)/(u) + 1)*(u))  /* Round n up to nearest multiple of u */


typedef struct {
  const char          *header;
  void              *(*reader)(FILE *);
  void               (*writer)(FILE *, void *);
  size_t               size;
} AS_MSG_callrecord;


typedef struct {

  //  Current line number of the reader
  int         LineNum;

  //  3-code of current read/write routine
  const char *Mcode;

  //  Memory allocation buffer for messages, and the current ceiling/top.
  char       *MemBuffer;
  int         MemMax;
  int         MemTop;

  //  Where we read messages into
  GenericMesg ReadMesg;

  //  The current line
#define MAX_LINE_LEN (128 * 1024)
  char        CurLine[MAX_LINE_LEN];

  //  The current calling table
  AS_MSG_callrecord CallTable[NUM_OF_REC_TYPES+1];
} AS_MSG_global_t;

extern AS_MSG_global_t  *AS_MSG_globals;






void    MakeSpace(const int size);
long    MoreSpace(const int size, const int boundary);
char   *ReadLine(FILE *fin);
char   *GetLine(FILE *fin, int skipComment);
void    MtypeError(const char * const name);
void    MtagError(const char * const tag);
void    MfieldError(const char * const mesg);
void    MgenError(const char * const mesg);
long    GetText(const char * const tag, FILE *fin, const int delnewlines);
long    GetString(const char * const tag, FILE *fin);
void    PutText(FILE *fout, const char * const tag, char * text, const int format);

//  These set the call function table in AS_MSG_global_t to the
//  appropriate functions for that format.
//
void    AS_MSG_setFormatVersion1(void);
void    AS_MSG_setFormatVersion2(void);



#define GET_TYPE(lvalue,format,name)                 \
{ char value[2];                                     \
  if (sscanf(GetLine(fin, TRUE),format,value) != 1)  \
    MtypeError(name);                                \
  lvalue =  *value;                                  \
}

/** DO NOT USE FOR SETTING ENUM TYPES FROM CHARACTER FLAGS!!!! **/
#define GET_FIELD(lvalue,format,emesg)                        \
{ if (sscanf(GetLine(fin, TRUE),format,&(lvalue)) != 1)       \
    MfieldError(emesg);                                       \
}

#define GET_PAIR(lvalue1,lvalue2,format,emesg)                        \
{ if (sscanf(GetLine(fin,TRUE),format,&(lvalue1),&(lvalue2)) != 2)    \
    MfieldError(emesg);                                               \
}

// lvalue3 is _ALWAYS_ type char coupled with a %1 format
#define GET_TRIPLE(lvalue1,lvalue2,lvalue3,format,emesg)                   \
{                                                                          \
  char value[2];                                                           \
  if(sscanf(GetLine(fin,TRUE),format,&(lvalue1),&(lvalue2), value) != 3)   \
    MfieldError(emesg);                                                    \
  lvalue3 = *value;                                                        \
}

#define GET_EOM                              \
{ if (GetLine(fin,TRUE)[0] != '}')           \
    MgenError("Expecting end of message");   \
}

#endif
