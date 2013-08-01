
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * Copyright (C) 2007, J. Craig Venter Institute.
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

const char *mainid = "$Id$";

#include "analyzePosMap.H"




//  Describes surrogates per library (globally sounds more useful).
//
//  numSurrogates
//  numInstances        (histogram; x surrogates had 50 instances)
//  numResolved         (histogram; fraction reads resolved per surrogate or per instance)
//  length              (histogram)
//  length-vs-instances
//
//  For salmon, we'd like this broken down by class of unitig.  Need some mechanism of coloring the unitigs.

void
analyzeSurrogates(char *outPrefix) {
}


void
analyzeDegenerates(char *outPrefix) {
}



#define ANALYZE_NOTHING              0x0000
#define ANALYZE_READS_IN_GAPS        0x0001  //  Report on unplaced reads that could close a gap
#define ANALYZE_LIBRARY_FATE         0x0002  //  Report on the fraction of each library that happy, stretched, misorient, etc mated
#define ANALYZE_SURROGATES           0x0004
#define ANALYZE_DEGENERATES          0x0008

int
main(int argc, char **argv) {
  char    *mapPrefix   = NULL;
  char    *outPrefix   = NULL;
  char    *gkpName     = NULL;
  uint32   analyze     = ANALYZE_NOTHING;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-p") == 0) {
      mapPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      outPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-A") == 0) {
      arg++;

      if        (strcmp(argv[arg], "readsingaps") == 0) {
        analyze |= ANALYZE_READS_IN_GAPS;

      } else if (strcmp(argv[arg], "libraryfate")    == 0) {
        analyze |= ANALYZE_LIBRARY_FATE;

      } else if (strcmp(argv[arg], "surrogates") == 0) {
        analyze |= ANALYZE_SURROGATES;

      } else {
        fprintf(stderr, "WARNING: unknown analysis '%s'\n", argv[arg]);
      }

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if (mapPrefix == NULL)
    err++;
  if (outPrefix == NULL)
    err++;
  if (analyze == ANALYZE_NOTHING)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s [...]\n", argv[0]);
    fprintf(stderr, "  -p posmap-prefix     prefix of posmap files (e.g., posmap-prefix.posmap.frgctg)\n");
    fprintf(stderr, "  -o output-prefix     prefix of output files\n");
    fprintf(stderr, "  -g gkpStore          path to gkpStore (not used)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -A analysis          select an analysis (multiple -A allowed)\n");
    fprintf(stderr, "                         readsingaps - probability that a gap can be filled with a read\n");
    fprintf(stderr, "                         libraryfate - per library detail of where each read ended up\n");
    fprintf(stderr, "                                     - and the status of each mate\n");
    fprintf(stderr, "                         surrogates  - \n");
    fprintf(stderr, "\n");

    if  (mapPrefix == NULL)
      fprintf(stderr, "ERROR: No posmap prefix (-p) supplied.\n");

    if  (outPrefix == NULL)
      fprintf(stderr, "ERROR: No output prefix (-p) supplied.\n");

    if (analyze == ANALYZE_NOTHING)
      fprintf(stderr, "ERROR: No analysis (-A) selected.\n");

    exit(1);
  }

  
  loadPosMap(mapPrefix, gkpName);




  //  From some single scaffold, count the number of mates that are external to it

#if 0
  for (uint32 si=0; si < scfNames.size(); si++) {
    uint32  internal = 0;
    uint32  external = 0;

    for (uint32 fi=0; fi < frgNames.size(); fi++) {
      if (frgDat[fi].con != si)
        continue;

      uint32  mi = frgMate[fi];

      if (mi == 0)
        continue;
      
      if (frgDat[fi].con == frgDat[mi].con)
        internal++;
      else
        external++;
    }

    fprintf(stderr, "scaffold "F_U32" %s %c internal "F_U32" external "F_U32"\n",
            si,
            scfNam[si].c_str(),
            scfDat[si].typ,
            internal,
            external);
  }
#endif


  if (analyze & ANALYZE_READS_IN_GAPS)
    analyzeGapFillProbability(outPrefix);

  if (analyze & ANALYZE_LIBRARY_FATE)
    analyzeLibraryFate(outPrefix);

  if (analyze & ANALYZE_SURROGATES)
    analyzeSurrogates(outPrefix);

  if (analyze & ANALYZE_DEGENERATES)
    analyzeDegenerates(outPrefix);

  exit(0);
}
