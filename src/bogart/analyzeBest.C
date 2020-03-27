
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
#include "AS_PER_seqStore.H"
#include "strings.H"


//  Read the 'best.edges' and 'best.contains' outputs from BOG, compute
//  how many singletons, spurs and contains are present per library.
//
//  Assumes it is run from 4-unitigger.

int
main(int argc, char **argv) {
  char  *seqName = 0L;
  char  *bEdge   = "best.edges";
  char  *bCont   = "best.contains";
  char  *bSing   = "best.singletons";

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-e") == 0) {
      bEdge = argv[++arg];

    } else if (strcmp(argv[arg], "-c") == 0) {
      bCont = argv[++arg];

    } else if (strcmp(argv[arg], "-s") == 0) {
      bSing = argv[++arg];

    } else {
      err++;
    }

    arg++;
  }
  if ((err) || (seqName == 0L)) {
    fprintf(stderr, "usage: %s -g seqName [-e best.edges] [-c best.contains]\n", argv[0]);
    exit(1);
  }

  fprintf(stderr, "Opening best edges and contains.\n");

  FILE *be = AS_UTL_openInputFile(bEdge);
  FILE *bc = AS_UTL_openInputFile(bCont);
  FILE *bs = AS_UTL_openInputFile(bSing);

  fprintf(stderr, "Loading read to library mapping.\n");

  sqStore    *seq = new sqStore(seqName, false, false);
  gkStream   *str = new gkStream(seq, 0, 0, GKFRAGMENT_INF);
  sqRead      fr;

  uint32   numFrg = seq->sqStore_getNumFragments();
  uint32   numLib = seq->sqStore_getNumLibraries();

  uint32  *frgToLib  = new uint32 [numFrg + 1];  memset(frgToLib, 0, sizeof(uint32) * (numFrg + 1));

  uint64  *readPerLib = new uint64 [numLib + 1];  memset(readPerLib, 0, sizeof(uint64) * (numLib + 1));
  uint64  *deldPerLib = new uint64 [numLib + 1];  memset(deldPerLib, 0, sizeof(uint64) * (numLib + 1));

  while (str->next(&fr)) {
    frgToLib[fr.sqRead_getReadIID()] = fr.sqRead_getLibraryIID();

    if (fr.sqRead_getIsDeleted())
      deldPerLib[fr.sqRead_getLibraryIID()]++;
    else
      readPerLib[fr.sqRead_getLibraryIID()]++;
  }

  delete str;

  uint64      *cntdPerLib = new uint64 [numLib + 1];  memset(cntdPerLib, 0, sizeof(uint64) * (numLib + 1));
  uint64      *cntrPerLib = new uint64 [numLib + 1];  memset(cntrPerLib, 0, sizeof(uint64) * (numLib + 1));
  uint32      *cntr       = new uint32 [numFrg + 1];  memset(cntr,       0, sizeof(uint32) * (numFrg + 1));

  uint64      *singPerLib = new uint64 [numLib + 1];  memset(singPerLib, 0, sizeof(uint64) * (numLib + 1));
  uint64      *spu5PerLib = new uint64 [numLib + 1];  memset(spu5PerLib, 0, sizeof(uint64) * (numLib + 1));
  uint64      *spu3PerLib = new uint64 [numLib + 1];  memset(spu3PerLib, 0, sizeof(uint64) * (numLib + 1));
  uint64      *dovePerLib = new uint64 [numLib + 1];  memset(dovePerLib, 0, sizeof(uint64) * (numLib + 1));

  int32          lineMax = 1048576;
  char          *line    = new char [lineMax];
  splitToWords   W;


  fprintf(stderr, "Processing best.singletons.\n");

  while (!feof(bs)) {
    fgets(line, lineMax, bs);  chomp(line);
    W.split(line);

    if (line[0] == '#')
      continue;

    singPerLib[frgToLib[W(0)]]++;
  }


  fprintf(stderr, "Processing best.contains.\n");

  while (!feof(bc)) {
    fgets(line, lineMax, bc);  chomp(line);
    W.split(line);

    if (line[0] == '#')
      continue;

    cntdPerLib[frgToLib[W(0)]]++;

    cntr[W(3)]++;
  }

  for (uint32 i=1; i<=numFrg; i++) {
    if (cntr[i] > 0)
      cntrPerLib[frgToLib[i]]++;
  }


  fprintf(stderr, "Processing best.edges.\n");

  while (!feof(be)) {
    fgets(line, lineMax, be);  chomp(line);
    W.split(line);

    if (line[0] == '#')
      continue;

    if      ((W(2) == 0) && (W(4) == 0))
      singPerLib[frgToLib[W(0)]]++;
    else if (W(2) == 0)
      spu5PerLib[frgToLib[W(0)]]++;
    else if (W(4) == 0)
      spu3PerLib[frgToLib[W(0)]]++;
    else
      dovePerLib[frgToLib[W(0)]]++;
  }

  fprintf(stderr, "libIID   libUID                             #frg     #del             cnt'd             cnt'r              sing             spur5             spur3              dove\n");

  for (uint32 i=0; i<numLib+1; i++) {
    double tot = readPerLib[i] + deldPerLib[i];

    fprintf(stderr, "%-8u %-30s %8" F_U64P " %8" F_U64P " (%4.1f%%)  %8" F_U64P " (%4.1f%%)  %8" F_U64P " (%4.1f%%)  %8" F_U64P " (%4.1f%%)  %8" F_U64P " (%4.1f%%)  %8" F_U64P " (%4.1f%%)  %8" F_U64P " (%4.1f%%)\n",
            i,
            seq->sqStore_getLibrary(i)->libraryName,
            readPerLib[i],
            deldPerLib[i], (deldPerLib[i] == 0) ? 0.0 : 100.0 * deldPerLib[i] / tot,
            cntdPerLib[i], (cntdPerLib[i] == 0) ? 0.0 : 100.0 * cntdPerLib[i] / tot,
            cntrPerLib[i], (cntrPerLib[i] == 0) ? 0.0 : 100.0 * cntrPerLib[i] / tot,
            singPerLib[i], (singPerLib[i] == 0) ? 0.0 : 100.0 * singPerLib[i] / tot,
            spu5PerLib[i], (spu5PerLib[i] == 0) ? 0.0 : 100.0 * spu5PerLib[i] / tot,
            spu3PerLib[i], (spu3PerLib[i] == 0) ? 0.0 : 100.0 * spu3PerLib[i] / tot,
            dovePerLib[i], (dovePerLib[i] == 0) ? 0.0 : 100.0 * dovePerLib[i] / tot);
  }

  delete seq;

  AS_UTL_closeFile(be, bEdge);
  AS_UTL_closeFile(bc, bCont);
  AS_UTL_closeFile(bs, bSing);

  delete [] frgToLib;

  delete [] readPerLib;
  delete [] deldPerLib;

  delete [] cntdPerLib;
  delete [] cntrPerLib;
  delete [] cntr;

  delete [] singPerLib;
  delete [] spu5PerLib;
  delete [] spu3PerLib;
  delete [] dovePerLib;

  delete [] line;

  return(0);
}

