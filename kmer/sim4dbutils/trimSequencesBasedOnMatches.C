#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include "bri++.H"
#include "sim4reader.h"

int
main(int argc, char ** argv) {
  FastA        *seqs = 0L;
  FastABuffer   seqsbuffer;
  FILE         *pfile = 0L;
  sim4polish   *p = 0L;

  if (argc == 1) {
    fprintf(stderr, "usage: %s -sequence s.fasta -polishes p.polished\n", argv[0]);
    exit(1);
  }

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-sequence", 2) == 0) {
      seqs = new FastA(argv[++arg], true);
    } else if (strncmp(argv[arg], "-polishes", 2) == 0) {
      errno = 0;
      pfile = fopen(argv[++arg], "r");
      if (errno) {
        fprintf(stderr, "Can't open '%s': %s\n", argv[arg], strerror(errno));
        exit(1);
      }
    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[arg]);
    }

    arg++;
  }

  if (seqs == 0L) {
    fprintf(stderr, "error: you need to specify '-sequence s.fasta'\n");
    exit(1);
  }

  if (pfile == 0L) {
    fprintf(stderr, "error: you need to specify '-polishes p.polished'\n");
    exit(1);
  }

  uint32   numseqs = seqs->numberOfSequences();
  uint32  *lrange = new uint32 [numseqs];
  uint32  *hrange = new uint32 [numseqs];

  for (uint32 i=0; i<numseqs; i++) {
    lrange[i] = ~uint32ZERO;
    hrange[i] =  uint32ZERO;
  }

  uint32  numRead = 0;

  while ((p = readPolish(pfile)) != 0L) {
    if (lrange[p->estID] > p->exons[0].estFrom-1)
      lrange[p->estID] = p->exons[0].estFrom-1;

    if (hrange[p->estID] < p->exons[p->numExons-1].estTo)
      hrange[p->estID] = p->exons[p->numExons-1].estTo;

    numRead++;
    if ((numRead & 0xff) == 0) {
      fprintf(stderr, "Reading matches: %u\r", numRead);
      fflush(stderr);
    }

    destroyPolish(p);
  }

  fprintf(stderr, "\n");

  uint32  seqcopylen = 128 * 1024;
  char   *seqcopy    = new char [seqcopylen + 1];
  char   *defcopy    = new char [128 * 1024];

  seqs->first(seqsbuffer);

  for (uint32 i=0; i<numseqs; i++, seqs->next(seqsbuffer)) {
    
    //  If there is no polish for the sequence, just write the whole
    //  thing out.  This is a hack, so that svi will run.
    //
    if (lrange[i] >= hrange[i]) {
      lrange[i] = 0;
      hrange[i] = seqsbuffer.sequenceLength();
    }

    if (lrange[i] < hrange[i]) {
      //seqs->seek(seqsbuffer, i);

      if (seqsbuffer.sequenceLength() > seqcopylen) {
        delete [] seqcopy;

        seqcopylen = seqsbuffer.sequenceLength() + 128 * 1024;
        seqcopy    = new char [seqcopylen + 1];
      }

      for (uint32 j=0, k=lrange[i]; k<hrange[i]; j++, k++)
        seqcopy[j] = seqsbuffer.sequence()[k];

      seqcopy[hrange[i] - lrange[i]] = 0;

      //  Mangle the defline
      //
      uint32 j = 0;
      for (j=0; !isspace(seqsbuffer.header()[j]) && j<seqsbuffer.headerLength(); j++)
        defcopy[j] = seqsbuffer.header()[j];

      defcopy[j] = 0;

#if 0
      fprintf(stdout, "%u] Trim from 0:%u to %u:%u\n",
              i,
              seqsbuffer.sequenceLength(),
              lrange[i],
              hrange[i]);
#endif

      fprintf(stdout, "%s trimmed to %u:%u\n%s\n", defcopy, lrange[i], hrange[i], seqcopy);
    }

    if ((i & 0x1ff) == 0) {
      fprintf(stderr, "Writing trimmed sequences: %u\r", i);
      fflush(stderr);
    }
  }

  fprintf(stderr, "\n");

  return(0);
}
