// This file is part of A2Amapper.
// Copyright (c) 2004 Applera Corporation
// Copyright (c) 2005 The J. Craig Venter Institute
// Author: Clark Mobarry
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

#include "heavychains.H"

#define BUFFERSIZE 1024

int main(int argc, char *argv[]) {
  int    beVerbose    = 0;
  char  *assemblyId1  = 0L;
  char  *assemblyId2  = 0L;
  double minScore     = 100.0;   // Default minimum of bp filled in a good run.
  int    maxJump      = 100000;  // Default maximum intra-run jump allowed in a good run.
  char  *inFileName   = 0L;
  char  *outFileName  = 0L;

  int  arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-v") == 0) {
      beVerbose++;
    } else if (strcmp(argv[arg], "-1") == 0) {
      assemblyId1 = argv[++arg];
    } else if (strcmp(argv[arg], "-2") == 0) {
      assemblyId2 = argv[++arg];
    } else if (strcmp(argv[arg], "-s") == 0) {
      minScore = atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-j") == 0) {
      maxJump = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-i") == 0) {
      inFileName = argv[++arg];
    } else if (strcmp(argv[arg], "-o") == 0) {
      outFileName = argv[++arg];
    } else {
      fprintf(stderr,"%s : unknown flag '-%s'\n", argv[0], *argv);
    }
  }

  FILE *inpF = fopen(inFileName, "r");
  FILE *outF = fopen(outFileName, "w");

  fprintf(outF,"! format atac 1.0\n");

  int    old_stra1  = -1;
  int    old_stra2  = -1; // True strand ordinals are non-negative.
  char   linebuffer[BUFFERSIZE] = {0};
  long   matchid    = 0;

  StrandPair *sp = new StrandPair(beVerbose, assemblyId1, assemblyId2, maxJump, minScore);
  TheStats   *ts = new TheStats(beVerbose, assemblyId1, assemblyId2, maxJump, minScore);

  bool  endOfInput = false;

  while (!endOfInput) {
    endOfInput = true;

    int   new_stra1   = -1;
    int   new_stra2   = -1;
    int   xln         = 0;
    int   yln         = 0;
    int   tmp_xlo     = 0;
    int   tmp_ylo     = 0;
    int   tmp_filled  = 0;  //  This is never changed!
    char  tmp_ori     = 0;

    if (fgets(linebuffer, BUFFERSIZE, inpF)) {
      endOfInput = false;

      if(linebuffer[0] == 'M') {
        char  classCode;
        char  subtype;
        char  selfId[100];
        char  parentId[100];
        char  new_ass1[100];
        char  new_ass2[100];
        int   xfl;
        int   yfl;

        if (12 != sscanf(linebuffer,
                         "%c %c %s %s %s %d %d %d %s %d %d %d\n",
                         &classCode,
                         &subtype,
                         selfId,
                         parentId,
                         new_ass1,
                         &tmp_xlo,
                         &xln,
                         &xfl,
                         new_ass2,
                         &tmp_ylo,
                         &yln,
                         &yfl)) {
          fprintf(stderr, "WARNING: short read on '%s'\n", linebuffer);
        }

#if 0
        printf("classCode=%c\n", classCode);
        printf("subtype  =%c\n", subtype);
        printf("selfId   =%s\n", selfId);
        printf("parentId =%s\n", parentId);
        printf("new_ass1 =%s\n", new_ass1);
        printf("xfl      =%d\n", xfl);
        printf("new_ass2 =%s\n", new_ass2);
        printf("yfl      =%d\n", yfl);
#endif

        if ((xfl != 1 && xfl != -1) ||
            (yfl != 1 && yfl != -1)) {
          fprintf(stderr, "ERROR: orientation wrong.\n%s\n", linebuffer);
          exit(1);
        }

        tmp_ori = (xfl == yfl ? 'f' : 'r');

        //  Parse the IID out of the ID
        //
        for (char *p = new_ass1; *p; p++)
          if (*p == ':')
            new_stra1 = atoi(p+1);

        for (char *p = new_ass2; *p; p++)
          if (*p == ':')
            new_stra2 = atoi(p+1);

      } else if ((linebuffer[0] == '#') || (linebuffer[0] == '!') || (linebuffer[0] == '/')) {
        fprintf(stderr,"%s",linebuffer);
      } else {
        fprintf(stderr, "UNRECOGNIZED: %s", linebuffer);
      }
    }

    if ((new_stra1 != old_stra1) ||
        (new_stra2 != old_stra2) || endOfInput) {
      if (sp->size() > 0) {
        sp->process();
        matchid = sp->print(outF, matchid);

        ts->add(sp);
      }
      delete sp;
      sp = new StrandPair(beVerbose, assemblyId1, assemblyId2, maxJump, minScore);
    }

    // Add the hit to the sp if we just read a point
    //
    if (linebuffer[0] == 'M') {
      sp->addHit(tmp_ori,
                 new_stra1, tmp_xlo, xln,
		 new_stra2, tmp_ylo, yln,
		 tmp_filled);

      old_stra1 = new_stra1;
      old_stra2 = new_stra2;
    }
  }

  ts->add(sp);
  ts->show(outF);

  delete sp;
  delete ts;

  fclose(inpF);
  fclose(outF);
}
