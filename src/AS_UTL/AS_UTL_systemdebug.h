
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
#ifndef AS_UTL_SYSTEM_DEBUG
#define AS_UTL_SYSTEM_DEBUG

#include "string.h"
#include "stdio.h"
size_t utl_strspn(const char *cs, const char *ct) {return strspn(cs,ct); }
size_t utl_strcspn(const char *cs, const char *ct) {return strcspn(cs,ct); }
size_t utl_strlen(const char *cs) {return strlen(cs);}

char 	*utl_strtok(char *s, const char *ct) { return strtok(s,ct); }
char	*utl_strchr(const char *ct, char c) { return strchr(ct,c); }
char	*utl_strpbrk(const char *cs, const char *ct) { return strpbrk(cs,ct); }
char	*utl_strstr(const char *cs, const char *ct) { return strstr(cs,ct); }
char	*utl_strrchr(const char *ct, char c) { return strrchr(ct,c); }

FILE *utl_stderr(void) { return stderr; }
FILE *utl_stdin(void) { return stdin; }
FILE *utl_stdout(void) { return stdout; }

void utl_showstring(FILE *out,const char *cs, int width) {
   int len=strlen(cs);
   int s=0;
   const char *p;
   while( s<len ) {
    p=cs+s;
    fprintf(out,"%.*s\n",width,p);
    s+=width;
   }
}
      

#endif
