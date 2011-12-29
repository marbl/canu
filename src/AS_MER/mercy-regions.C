
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2006-2007, J. Craig Venter Institute
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

const char *mainid = "$Id: mercy-regions.C,v 1.6 2011-12-29 09:26:03 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "AS_global.h"
#include "AS_UTL_splitToWords.H"
#include "AS_UTL_intervalList.H"

//  This reads the assembly frgctg, varctg and merQC badmers, computes
//  the number and location of bad-mer, bad-var regions, and their
//  depth, in contig space.
//
//  File paths are hardcoded.
//  This code ONLY works on 64-bit hardware, but it's easy to fix.

using namespace std;
#include <map>


void
readDepth(char *depthname, map<uint64,intervalDepth*> &lowCoverage) {
  char                         line[1024] = {0};
  map<uint64,intervalList*>    ILs;

  fprintf(stderr, "Reading depth from '%s'\n", depthname);

  errno = 0;
  FILE *F = fopen(depthname, "r");
  if (errno)
    fprintf(stderr, "failed to open '%s': %s\n", depthname, strerror(errno)), exit(1);

  uint32 i=0;

  fgets(line, 1024, F);
  while (!feof(F)) {
    splitToWords   W(line);

    uint64  uid = strtoul(W[1], 0L, 10);
    uint32  beg = strtoul(W[2], 0L, 10);
    uint32  end = strtoul(W[3], 0L, 10);

    if (beg > end)
      fprintf(stderr, "ERROR: l="F_U32" h="F_U32"\n", beg, end);

    if (ILs[uid] == 0L)
      ILs[uid] = new intervalList();
    ILs[uid]->add(beg, end - beg);

    i++;

    fgets(line, 1024, F);
  }

  fclose(F);
  fprintf(stderr, " "F_U32" lines.\n", i);

  map<uint64,intervalList*>::iterator    it = ILs.begin();
  map<uint64,intervalList*>::iterator    ed = ILs.end();

  while (it != ed) {
    lowCoverage[it->first] = new intervalDepth(*it->second);
    delete it->second;
    it->second = 0L;
    it++;
  }
}


void
readVariation(char *depthname, map<uint64,intervalList*> &variation) {
  char                         line[1024 * 1024] = {0};

  fprintf(stderr, "Reading variation from '%s'\n", depthname);

  errno = 0;
  FILE *F = fopen(depthname, "r");
  if (errno)
    fprintf(stderr, "failed to open '%s': %s\n", depthname, strerror(errno)), exit(1);

  uint32 i=0;

  fgets(line, 1024 * 1024, F);
  while (!feof(F)) {
    splitToWords   W(line);

    uint64  uid = strtoul(W[1], 0L, 10);
    uint32  beg = strtoul(W[2], 0L, 10);
    uint32  end = strtoul(W[3], 0L, 10);

    if (variation[uid] == 0L)
      variation[uid] = new intervalList();
    variation[uid]->add(beg, end - beg);

    i++;

    fgets(line, 1024 * 1024, F);
  }

  fclose(F);
  fprintf(stderr, " "F_U32" lines.\n", i);
}


void
readBadMers(char *depthname, map<uint64,intervalList*> &badMers) {
  char                         line[1024] = {0};

  fprintf(stderr, "Reading badMers from '%s'\n", depthname);

  errno = 0;
  FILE *F = fopen(depthname, "r");
  if (errno)
    fprintf(stderr, "failed to open '%s': %s\n", depthname, strerror(errno)), exit(1);

  uint32 i=0;

  fgets(line, 1024, F);
  while (!feof(F)) {
    splitToWords   W(line);

    //  Change every non-digit to a space in the first word.
    for (uint32 z=strlen(W[0])-1; z--; )
      if (!isdigit(W[0][z]))
        W[0][z] = ' ';

    uint64  uid = strtoul(W[0], 0L, 10);
    uint32  beg = strtoul(W[3], 0L, 10);
    uint32  end = strtoul(W[4], 0L, 10);

    if (badMers[uid] == 0L)
      badMers[uid] = new intervalList();
    badMers[uid]->add(beg, end - beg);

    i++;

    fgets(line, 1024, F);
  }

  fclose(F);
  fprintf(stderr, " "F_U32" lines.\n", i);
}



int
main(int argc, char **argv) {
  map<uint64,intervalList*>    badMers;
  map<uint64,intervalList*>    variation;
  map<uint64,intervalDepth*>   lowCoverage;

  bool  showDepthIntersect    = false;
  bool  showVariantIntersect  = false;
  bool  showVarDepthIntersect = false;

  argc = AS_configure(argc, argv);

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-D") == 0) {

    } else if (strcmp(argv[arg], "-pd") == 0) {
      showDepthIntersect = true;
    } else if (strcmp(argv[arg], "-pv") == 0) {
      showVariantIntersect = true;
    } else if (strcmp(argv[arg], "-pvd") == 0) {
      showVarDepthIntersect = true;
    } else {
      fprintf(stderr, "usage: %s [-D debugfile] [-pd] [-pv] [-pvd]\n", argv[0]);
      fprintf(stderr, " -pd    print bad mers regions isect depth\n");
      fprintf(stderr, " -pv    print bad mers regions isect variants\n");
      fprintf(stderr, " -pvd   print bad mers regions isect both variants and depth\n");
      exit(1);
    }
    arg++;
  }

#if 1
  //  HuRef6, in the assembly directory.
  //
  readDepth    ("/project/huref6/assembly/h6/9-terminator/h6.posmap.frgctg", lowCoverage);
  readVariation("/project/huref6/assembly/h6/9-terminator/h6.posmap.varctg", variation);
  readBadMers  ("/project/huref6/assembly/h6-mer-validation/h6-ms22-allfrags-normalcontigs.badmers.0.singlecontig.zerofrag.badmers", badMers);
#endif

#if 0
  //  HuRef6, ws=25, in the assembly directory.
  //
  readDepth    ("/project/huref6/assembly/h6/9-terminator-ws25/h6.posmap.frgctg", lowCoverage);
  readVariation("/project/huref6/assembly/h6/9-terminator-ws25/h6.posmap.varctg", variation);
  readBadMers  ("/project/huref6/assembly/h6-mer-validation/h6-version4-ws25/h6-ms22-allfrags-normalcontigs.badmers.0.singlecontig.zerofrag.badmers", badMers);
#endif

#if 0
  //  Our scratch huref
  //
  readDepth    ("/project/huref6/redo_consensus-gennady/mer-validation/h6tmp.posmap.frgctg", lowCoverage);
  readVariation("/project/huref6/redo_consensus-gennady/mer-validation/h6tmp.posmap.varctg", variation);
  readBadMers  ("/project/huref6/redo_consensus-gennady/mer-validation/h6tmp-ms22-allfrags-allcontigs.badmers.0.singlecontig.zerofrag.badmers", badMers);
#endif

  uint32   badBegDepth[1024] = {0};
  uint32   badEndDepth[1024] = {0};

  uint32   badDepth[32][32];
  for (uint32 i=0; i<32; i++)
    for (uint32 j=0; j<32; j++)
      badDepth[i][j] = 0;

  map<uint64,intervalList*>::iterator    it = badMers.begin();
  map<uint64,intervalList*>::iterator    ed = badMers.end();
  while (it != ed) {
    uint64         uid        = it->first;

    intervalList  *Iv = variation[uid];
    intervalList  *Ib = badMers[uid];
    intervalList  *Ii = 0L;
    intervalDepth *Id = lowCoverage[uid];

    if (Iv)
      Iv->merge();
    if (Ib)
      Ib->merge();

    if (Iv && Ib) {
      Ii = new intervalList();
      Ii->intersect(*Iv, *Ib);
    }


    if (Ii) {
      uint32 ii = 0;
      uint32 id = 0;

      while ((ii < Ii->numberOfIntervals()) &&
             (id < Id->numberOfIntervals())) {

        //  We want to count the number of times a badmer region
        //  begins/ends in some depth.

        //fprintf(stderr, "testing beg        "F_U32" "F_U32" -- "F_U32" "F_U32"\n",
        //        Ii->lo(ii), Ii->hi(ii), Id->lo(id), Id->hi(id));

        uint32  beg = 0;
        uint32  end = 0;

        //  Low points are not allowed to be equal to high points, skip to the next
        while ((id < Id->numberOfIntervals()) &&
               (Id->hi(id) <= Ii->lo(ii))) {
          id++;
          //fprintf(stderr, "testing beg (m)     "F_U32" "F_U32" -- "F_U32" "F_U32"\n",
          //        Ii->lo(ii), Ii->hi(ii), Id->lo(id), Id->hi(id));
        }
        if (id < Id->numberOfIntervals()) {
          uint32 lo = Id->lo(id);
          uint32 hi = Id->hi(id);

          //  Low points are not allowed to be equal to high points.
          if ((lo <= Ii->lo(ii)) && (Ii->lo(ii) < hi)) {
            beg = Id->de(id);
          } else {
            fprintf(stderr, "failed to find begin "F_U64" "F_U64" -- "F_U64" "F_U64" "F_U32"\n",
                    Ii->lo(ii), Ii->hi(ii), Id->lo(id), Id->hi(id), Id->de(id));
            if (id > 0)
              fprintf(stderr, "                     "F_U64" "F_U64" -- "F_U64" "F_U64" "F_U32"\n",
                      Ii->lo(ii), Ii->hi(ii), Id->lo(id-1), Id->hi(id-1), Id->de(id-1));
            //exit(1);
          }
        }

        //fprintf(stderr, "testing end        "F_U64" "F_U64" -- "F_U64" "F_U64"\n",
        //        Ii->lo(ii), Ii->hi(ii), Id->lo(id), Id->hi(id));

        //  High points can be equal.
        while ((id < Id->numberOfIntervals()) &&
               (Id->hi(id) < Ii->hi(ii))) {
          id++;
          //fprintf(stderr, "testing end (m)    "F_U64" "F_U64" -- "F_U64" "F_U64"\n",
          //        Ii->lo(ii), Ii->hi(ii), Id->lo(id), Id->hi(id));
        }
        if (id < Id->numberOfIntervals()) {
          uint32 lo = Id->lo(id);
          uint32 hi = Id->hi(id);

          //  High points aren't allowed to be equal to lo, but can be equal to hi.
          if ((lo < Ii->hi(ii)) && (Ii->hi(ii) <= hi)) {
            end = Id->de(id);
          } else {
            fprintf(stderr, "failed to find end "F_U64" "F_U64" -- "F_U64" "F_U64" "F_U32"\n",
                    Ii->lo(ii), Ii->hi(ii), Id->lo(id), Id->hi(id), Id->de(id));
            if (id > 0)
              fprintf(stderr, "                     "F_U64" "F_U64" -- "F_U64" "F_U64" "F_U32"\n",
                      Ii->lo(ii), Ii->hi(ii), Id->lo(id-1), Id->hi(id-1), Id->de(id-1));
            //exit(1);
          }
        }

        badBegDepth[beg]++;
        badEndDepth[end]++;

        fprintf(stdout, F_U64"\t"F_U64"\t"F_U64"\tdepth="F_U32","F_U32"\n",
                uid, Ii->lo(ii), Ii->hi(ii), beg, end);

        if ((beg < 32) && (end < 32))
          badDepth[beg][end]++;

        ii++;
      }
    }

    it++;
  }

  uint32 bb = 0;
  uint32 be = 0;
  for (uint32 x=0; x<32; x++) {
    fprintf(stdout, F_U32"\t"F_U32"\t"F_U32"\n", x, badBegDepth[x], badEndDepth[x]);
    bb += badBegDepth[x];
    be += badEndDepth[x];
  }
  fprintf(stdout, "total\t"F_U32"\t"F_U32"\n", bb, be);

  for (uint32 i=0; i<30; i++) {
    for (uint32 j=0; j<30; j++)
      fprintf(stdout, "%5u", badDepth[i][j]);
    fprintf(stdout, "\n");
  }

  return(0);
}
