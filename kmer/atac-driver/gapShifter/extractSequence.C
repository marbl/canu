// This file is part of A2Amapper.
// Copyright (c) 2006 J. Craig Venter Institute
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

#include "atac.H"
#include "bio++.H"
#include "seqCache.H"

//  Reads a set of matches and outputs sequence that was mapped.  Filters matches, etc.

void
extractA(seqCache *A, seqCache *B,
         FILE *Aoutput, FILE *Boutput,
         uint32 Aiid,
         uint32 Biid,
         atacMatchList &ML) {

}


void
usage(char *name) {
  fprintf(stderr, "usage: %s [-OP output.fasta] [-t trfile] -m matches\n", name);
  fprintf(stderr, "   OP\n");
  fprintf(stderr, "   -a        extract all unmapped sequence in A\n");
  fprintf(stderr, "   -b        extract all unmapped sequence in B\n");
  fprintf(stderr, "   -ar       extract within run unmapped sequence in A\n");
  fprintf(stderr, "   -br       extract within run unmapped sequence in B\n");
  fprintf(stderr, "             BOTH -ar and -br need to be specified!\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   -t        mask out tandem repeats listed in trfile\n");
}

FILE *
openOutputFile(char *name) {
  errno = 0;
  FILE *R = fopen(name, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", name, strerror(errno)), exit(1);
  return(R);
}

int
main(int argc, char *argv[]) {
  char         *matchesFile = 0L;
  FILE         *Aoutput = 0L;
  FILE         *Boutput = 0L;
  uint32        Aiid = ~uint32ZERO;
  uint32        Biid = ~uint32ZERO;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-m") == 0) {
      matchesFile = argv[++arg];
    } else if (strcmp(argv[arg], "-a") == 0) {
      Aoutput = openOutputFile(argv[++arg]);
    } else if (strcmp(argv[arg], "-b") == 0) {
      Boutput = openOutputFile(argv[++arg]);

    } else if (strcmp(argv[arg], "-1") == 0) {
      Aiid = strtouint32(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-2") == 0) {
      Biid = strtouint32(argv[++arg], 0L);
    } else {
      usage(argv[0]);
      exit(1);
    }
    arg++;
  }

  if (matchesFile == 0L)
    usage(argv[0]), exit(1);

  atacFile       AF(matchesFile);
  atacMatchList &ML = *AF.matches();

  seqCache  *A = new seqCache(AF.assemblyFileA(), 0, true);
  seqCache  *B = new seqCache(AF.assemblyFileB(), 0, true);

  A->loadAllSequences();
  B->loadAllSequences();

  for (uint32 x=0; x<ML.numMatches(); x++) {
    atacMatch *m = ML[x];

    if (((Aiid == ~uint32ZERO) || (Aiid == m->iid1)) &&
        ((Biid == ~uint32ZERO) || (Biid == m->iid2))) {

      if (Aoutput) {
        seqInCore *S = A->getSequenceInCore(m->iid1);

        fprintf(Aoutput, "%s extracted from iid "uint32FMT" pos "uint32FMT" "uint32FMT" match %s(%s)\n",
                S->header(), S->getIID(),
                m->pos1, m->pos1 + m->len1,
                m->matchuid, m->parentuid);
        fwrite(S->sequence() + m->pos1, sizeof(char), m->len1, Aoutput);
        fprintf(Aoutput, "\n");
      }

      if (Boutput) {
        seqInCore *S = B->getSequenceInCore(m->iid2);

        fprintf(Boutput, "%s extracted from iid "uint32FMT" pos "uint32FMT" "uint32FMT" match %s(%s)\n",
                S->header(), S->getIID(),
                m->pos2, m->pos2 + m->len2,
                m->matchuid, m->parentuid);
        fwrite(S->sequence() + m->pos2, sizeof(char), m->len2, Boutput);
        fprintf(Boutput, "\n");
      }
    }
  }

  if (Aoutput)  fclose(Aoutput);
  if (Boutput)  fclose(Boutput);

  return(0);
}
