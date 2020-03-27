
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "runtime.H"
#include "sqStore.H"
#include "ovStore.H"

#include <vector>

using namespace std;


int
main(int argc, char **argv) {
  char                  *seqStoreName = NULL;
  sqStore               *seqStore = NULL;

  ovOverlapDisplayType   dt = ovOverlapAsCoords;
  vector<char *>         files;


  int32     arg = 1;
  int32     err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-coords") == 0) {
      dt = ovOverlapAsCoords;

    } else if (strcmp(argv[arg], "-hangs") == 0) {
      dt = ovOverlapAsHangs;

    } else if (strcmp(argv[arg], "-unaligned") == 0) {
      dt = ovOverlapAsUnaligned;

    } else if (fileExists(argv[arg])) {
      files.push_back(argv[arg]);

    } else {
      fprintf(stderr, "ERROR:  invalid arg '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if ((seqStoreName == NULL) && (dt == ovOverlapAsCoords))
    err++;

  if ((err) || (files.size() == 0)) {
    fprintf(stderr, "usage: %s [options] file.ovb[.gz]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S             seqStore (needed for -coords, the default)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -coords        output coordiantes on reads\n");
    fprintf(stderr, "  -hangs         output hangs on reads\n");
    fprintf(stderr, "  -unaligned     output unaligned regions on each read\n");
    fprintf(stderr, "\n");

    if ((seqStoreName == NULL) && (dt == ovOverlapAsCoords))
      fprintf(stderr, "ERROR:  -coords mode requires a seqStore (-S)\n");

    if (files.size() == 0)
      fprintf(stderr, "ERROR:  no overlap files supplied\n");

    exit(1);
  }

  if (seqStoreName)
    seqStore = new sqStore(seqStoreName);

  char  *ovStr = new char [1024];

  for (uint32 ff=0; ff<files.size(); ff++) {
    ovFile      *of = new ovFile(seqStore, files[ff], ovFileFull);
    ovOverlap   ov;

    while (of->readOverlap(&ov))
      fputs(ov.toString(ovStr, dt, true), stdout);

    delete of;
  }

  delete [] ovStr;

  delete seqStore;

  exit(0);
}
