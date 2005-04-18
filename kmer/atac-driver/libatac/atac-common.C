// This file is part of A2Amapper.
// Copyright (c) 2005 J. Craig Venter Institute
// Author: Brian Walenz
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "bio++.H"


void
readHeader(char *inLine, FILE *in, char *file1, char *file2, FILE *out) {

  fgets(inLine, 1024, in);
  while (!feof(in) && (inLine[0] != 'M')) {

    //  Copy the line to the output, if there is an output.
    //
    if (out)
      fputs(inLine, out);

    if (file1) {
      if (strncmp(inLine, "/assemblyFile1", 14) == 0) {
        //  Skip any whitespace, the =, and more whitespace.  Copy.  Nuke whitespace at the end.
        char *tmp = inLine + 14;
        while (isspace(*tmp))  tmp++;
        while (*tmp == '=') tmp++;
        while (isspace(*tmp))  tmp++;
        strcpy(file1, tmp);
        chomp(file1);
      }
    }
    if (file2) {
      if (strncmp(inLine, "/assemblyFile2", 14) == 0) {
        char *tmp = inLine + 14;
        while (isspace(*tmp))  tmp++;
        while (*tmp == '=') tmp++;
        while (isspace(*tmp))  tmp++;
        strcpy(file2, tmp);
        chomp(file2);
      }
    }

    fgets(inLine, 1024, in);
  }

  if (((file1) && (file1[0] == 0)) || ((file2) && (file2[0] == 0))) {
    fprintf(stderr, "I didn't find /assemblyFile1= or /assemblyFile2=\n");
    exit(1);
  }
}





bool
decodeMatch(splitToWords &W,
            u32bit &iid1, u32bit &pos1, u32bit &len1, u32bit &ori1,
            u32bit &iid2, u32bit &pos2, u32bit &len2, u32bit &ori2) {

  if ((W[0][0] == 'M') && (W[1][0] == 'u')) {
    char *tmp = W[4];
    while (*tmp && (*tmp != ':'))
      tmp++;
    if (*tmp)
      iid1 = strtou32bit(tmp+1, 0L);

    tmp = W[8];
    while (*tmp && (*tmp != ':'))
      tmp++;
    if (*tmp)
      iid2 = strtou32bit(tmp+1, 0L);

    pos1 = strtou32bit(W[5], 0L);
    len1 = strtou32bit(W[6], 0L);
    ori1 = (W[7][0] == '-') ? 0 : 1;

    pos2 = strtou32bit(W[9], 0L);
    len2 = strtou32bit(W[10], 0L);
    ori2 = (W[11][0] == '-') ? 0 : 1;

    return(true);
  }

  return(false);
}

