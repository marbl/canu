
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_GKP/fastqSample.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-FEB-22 to 2013-AUG-01
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-AUG-06 to 2015-FEB-24
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include <vector>
#include <algorithm>

using namespace std;

#define MAXLEN 1024*1024

class aRead {
public:
  aRead() {
    memset(a, 0, sizeof(char) * MAXLEN);
    memset(b, 0, sizeof(char) * MAXLEN);
    memset(c, 0, sizeof(char) * MAXLEN);
    memset(d, 0, sizeof(char) * MAXLEN);
  };

  bool  read(FILE *F) {

    a[0] = 0;
    b[0] = 0;
    c[0] = 0;
    d[0] = 0;

    if (F == NULL)
      return(false);

    fgets(a, MAXLEN, F);
    fgets(b, MAXLEN, F);
    fgets(c, MAXLEN, F);
    fgets(d, MAXLEN, F);

    if (feof(F) == true)
      return(false);

    if ((a[0] != '@') || (c[0] != '+')) {
      fprintf(stderr, "ERROR:  Not FastQ.  Read lines:\n");
      fprintf(stderr, "  %s", a);
      fprintf(stderr, "  %s", b);
      fprintf(stderr, "  %s", c);
      fprintf(stderr, "  %s", d);
      exit(1);
    }

    return(true);
  };
  void  write(FILE *F) {
    if (F == NULL)
      return;
    fputs(a, F);
    fputs(b, F);
    fputs(c, F);
    fputs(d, F);
  };

  uint32  length(void) {
    return(strlen(b) - 1);  //  Newline still exists
  };

private:
  char  a[MAXLEN];
  char  b[MAXLEN];
  char  c[MAXLEN];
  char  d[MAXLEN];
};



class anInput {
public:
  anInput() {
    id   = 0;
    len  = 0;
  };
  anInput(uint64 id_, uint32 len1_, uint32 len2_) {
    id   = id_;
    len  = len1_ + len2_;
  };

  uint64  id;
  uint32  len;
};


inline
bool
anInputByLongest(const anInput &a, const anInput &b) {
  return(a.len > b.len);
};


int
main(int argc, char **argv) {
  aRead    *Ar = new aRead, *Br = new aRead;
  FILE     *Ai = NULL,      *Bi = NULL;
  FILE     *Ao = NULL,      *Bo = NULL;

  vector<anInput>   ids;
  vector<bool>      sav;

  char     *INPNAME     = NULL;
  char     *OUTNAME     = NULL;
  bool      AUTONAME    = false;

  uint64    NUMINPUT    = 0;        //  Number of pairs in the input
  uint64    READLENGTH  = 0;        //  For mated reads, 2x read size
  bool      isMated     = true;

  uint32    MINLEN      = 0;
  uint32    LONGEST     = false;

  uint64    GENOMESIZE  = 0;        //  Size of the genome in bp
  double    COVERAGE    = 0;        //  Desired coverage in output

  uint64    NUMOUTPUT   = 0;        //  Number of pairs to output

  double    FRACTION    = 0.0;      //  Desired fraction of the input

  uint64    BASES       = 0;        //  Desired amount of sequence

  char      path1[FILENAME_MAX];
  char      path2[FILENAME_MAX];

  srand48(time(NULL));

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-I") == 0) {
      INPNAME = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      OUTNAME = argv[++arg];

    } else if (strcmp(argv[arg], "-A") == 0) {
      AUTONAME = true;

    } else if (strcmp(argv[arg], "-T") == 0) {
      NUMINPUT = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-L") == 0) {
      READLENGTH = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-U") == 0) {
      isMated = false;

    } else if (strcmp(argv[arg], "-m") == 0) {
      MINLEN = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-max") == 0) {
      LONGEST = true;


    } else if (strcmp(argv[arg], "-g") == 0) {
      GENOMESIZE = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-c") == 0) {
      COVERAGE = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-p") == 0) {
      NUMOUTPUT = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-f") == 0) {
      FRACTION = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-b") == 0) {
      BASES = atoi(argv[++arg]);

    } else {
      err++;
    }

    arg++;
  }
  if (INPNAME == NULL)
    err++;
  if (OUTNAME == NULL)
    OUTNAME = INPNAME;
  if ((GENOMESIZE == 0) && (COVERAGE != 0))
    err++;
  if ((COVERAGE == 0) && (NUMOUTPUT == 0) && (FRACTION == 0.0) && (BASES == 0))
    err++;
  if (err) {
    fprintf(stderr, "\n");
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    fprintf(stderr, "  Input Specification\n");
    fprintf(stderr, "    -I NAME  input name (prefix) of the reads\n");
    fprintf(stderr, "    -T T     total number of mate pairs in the input (if not supplied, will be counted)\n");
    fprintf(stderr, "    -L L     length of a single read (if not supplied, will be determined)\n");
    fprintf(stderr, "    -U       reads are unmated, expected in *.u.fastq\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Output Specification\n");
    fprintf(stderr, "    -O NAME  output name (prefix) of the reads (default is same as -I)\n");
    fprintf(stderr, "    -A       automatically include coverage or number of reads in the output name\n");
    fprintf(stderr, "    -m L     ignore reads shorter than L bases\n");
    fprintf(stderr, "    -max     don't sample randomly, pick the longest reads\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Method 1: specify desired output coverage:\n");
    fprintf(stderr, "    -g G     genome size\n");
    fprintf(stderr, "    -c C     desired coverage in the output reads\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Method 2: specify desired number of output pairs\n");
    fprintf(stderr, "    -p N     for mated reads, output 2N reads, or N pairs of reads\n");
    fprintf(stderr, "             for unmated reads, output N reads\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Method 3: specify a desired fraction of the input:\n");
    fprintf(stderr, "    -f F     output F * T pairs of reads (T as above in -t option)\n");
    fprintf(stderr, "             0.0 < F <= 1.0\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Method 4: specify a desired total length\n");
    fprintf(stderr, "    -b B     output reads/pairs until B bases is exceeded\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Samples reads from paired Illumina reads NAME.1.fastq and NAME.2.fastq and outputs:\n");
    fprintf(stderr, "    NAME.Cx.1.fastq and N.Cx.2.fastq (for coverage based sampling)\n");
    fprintf(stderr, "    NAME.n=N.1.fastq and N.n=N.2.fastq (for coverage based sampling)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "If -T is not supplied, the number of reads will be counted for you.\n");
    fprintf(stderr, "\n");
    if (INPNAME == NULL)
      fprintf(stderr, "ERROR: no name supplied with -I.\n");
    if ((GENOMESIZE == 0) && (COVERAGE != 0))
      fprintf(stderr, "ERROR: no genome size supplied with -g (when using -c)\n");
    if ((COVERAGE == 0) && (NUMOUTPUT == 0) && (FRACTION == 0.0) && (BASES == 0))
      fprintf(stderr, "ERROR: no method supplied with -c, -p, -f or -b\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  //
  //  We know not enough about the reads, and are forced to scan the entire
  //  inputs.
  //

  uint64   totBasesInInput = 0;
  uint64   totPairsInInput = 0;

  if ((NUMINPUT == 0) || (READLENGTH == 0)) {
    uint64  Ac = 0;
    uint64  Bc = 0;

    fprintf(stderr, "Counting the number of reads in the input.\n");

    sprintf(path1, "%s.%c.fastq", INPNAME, (isMated == true) ? '1' : 'u');
    sprintf(path2, "%s.%c.fastq", INPNAME, (isMated == true) ? '2' : 'u');

    errno = 0;
    Ai = fopen(path1, "r");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", path1, strerror(errno)), exit(1);

    if (isMated == true) {
      errno = 0;
      Bi = fopen(path2, "r");
      if (errno)
        fprintf(stderr, "Failed to open '%s': %s\n", path2, strerror(errno)), exit(1);
    }

    bool  moreA = Ar->read(Ai);
    bool  moreB = Br->read(Bi);

    while (moreA || moreB) {
      uint32  lA = Ar->length();
      uint32  lB = Br->length();

      if (lA > 0)  Ac++;
      if (lB > 0)  Bc++;

      ids.push_back(anInput(totPairsInInput, lA, lB));
      sav.push_back(false);

      totPairsInInput += 1;
      totBasesInInput += lA + lB;

      moreA = Ar->read(Ai);
      moreB = Br->read(Bi);
    }

    if (Ai)  fclose(Ai);
    if (Bi)  fclose(Bi);

    fprintf(stderr, "Found "F_U64" bases and "F_U64" reads in '%s'\n",
            totBasesInInput, totPairsInInput, path1);

    if (Ac != Bc) {
      fprintf(stderr, "ERROR:  Number of reads in the .1 and .2 files must be the same.\n");
      exit(1);
    }
  }

  //  Otherwise, both NUMINPUT and READLENGTH are defined, so we can fill out ids
  //  with defaults.
  else {
    if (isMated == false)
      for (uint64 ii=0; ii<NUMINPUT; ii++) {
        ids.push_back(anInput(ii, READLENGTH, 0));
        sav.push_back(false);
      }
    else
      for (uint64 ii=0; ii<NUMINPUT; ii++) {
        ids.push_back(anInput(ii, READLENGTH, READLENGTH));
        sav.push_back(false);
      }

    totBasesInInput = NUMINPUT * READLENGTH;
    totPairsInInput = NUMINPUT;
  }

  //
  //  Decide how much data to output, either as a read or base limit.
  //

  uint64   nBasesToOutput = 0;
  uint64   nPairsToOutput = 0;


  if (GENOMESIZE > 0) {
    nBasesToOutput = (uint64)(COVERAGE * GENOMESIZE);
  }

  if (NUMOUTPUT > 0) {
    nPairsToOutput = NUMOUTPUT;
  }

  if (FRACTION > 0) {
    nPairsToOutput = (uint64)(FRACTION * totPairsInInput);
  }

  if (BASES > 0) {
    nBasesToOutput = BASES;
  }

  if (totBasesInInput < nBasesToOutput)
    fprintf(stderr, "ERROR: not enough reads, "F_U64" bp in input, "F_U64" needed for desired .....\n",
            totBasesInInput, nBasesToOutput),
      exit(1);

  if (totPairsInInput < nPairsToOutput)
    fprintf(stderr, "ERROR: not enough reads, "F_U64" %s in input, "F_U64" needed for desired ......\n",
            totPairsInInput, (isMated) ? "pairs" : "reads", nPairsToOutput),
      exit(1);

  //fprintf(stderr, "OUTPUT: %lu bases\n", nBasesToOutput);
  //fprintf(stderr, "OUTPUT: %lu pairs\n", nPairsToOutput);


  //
  //  Randomize the ID list, or sort by length.
  //

  if (LONGEST) {
    fprintf(stderr, "Sorting by length\n");
    sort(ids.begin(), ids.end(), anInputByLongest);

  } else {
    fprintf(stderr, "Shuffling sequences\n");
    for (uint64 i=0; i<totPairsInInput; i++) {
      uint64   p = (uint64)lrand48() % totPairsInInput;
      anInput  a = ids[p];

      ids[p] = ids[i];
      ids[i] = a;
    }
  }


  //
  //  Decide what to output.
  //

  if (nPairsToOutput > 0) {
    uint64  nPairs = 0;

    for (uint64 i=0; i<totPairsInInput; i++) {
      if (nPairs < nPairsToOutput)
        sav[ids[i].id] = true;
      else
        break;

      nPairs++;
    }
  }

  if (nBasesToOutput > 0) {
    uint64  nBases = 0;

    for (uint64 i=0; i<totPairsInInput; i++) {
      if (nBases < nBasesToOutput)
        sav[ids[i].id] = true;
      else
        break;

      nBases += ids[i].len;
    }
  }


  //
  //  Count.
  //

  nPairsToOutput = 0;
  nBasesToOutput = 0;

  for (uint64 i=0; i<totPairsInInput; i++) {
    if (sav[ids[i].id] == true) {
      nPairsToOutput += 1;
      nBasesToOutput += ids[i].len;
    }
  }


  //
  //  Do the output
  //

  sprintf(path1, "%s.%c.fastq", INPNAME, (isMated == true) ? '1' : 'u');
  sprintf(path2, "%s.%c.fastq", INPNAME, (isMated == true) ? '2' : 'u');

  errno = 0;
  Ai = fopen(path1, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", path1, strerror(errno)), exit(1);

  if (isMated == true) {
    errno = 0;
    Bi = fopen(path2, "r");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", path2, strerror(errno)), exit(1);
  }


  if (AUTONAME == false) {
    sprintf(path1, "%s.%c.fastq", OUTNAME, (isMated == true) ? '1' : 'u');
    sprintf(path2, "%s.%c.fastq", OUTNAME, (isMated == true) ? '2' : 'u');

  } else if (GENOMESIZE > 0) {
    sprintf(path1, "%s.x=%07.3f.n=%09"F_U64P".%c.fastq", OUTNAME, (double)nBasesToOutput / GENOMESIZE, nPairsToOutput, (isMated == true) ? '1' : 'u');
    sprintf(path2, "%s.x=%07.3f.n=%09"F_U64P".%c.fastq", OUTNAME, (double)nBasesToOutput / GENOMESIZE, nPairsToOutput, (isMated == true) ? '2' : 'u');

  } else {
    sprintf(path1, "%s.x=UNKNOWN.n=%09"F_U64P".%c.fastq", OUTNAME, nPairsToOutput, (isMated == true) ? '1' : 'u');
    sprintf(path2, "%s.x=UNKNOWN.n=%09"F_U64P".%c.fastq", OUTNAME, nPairsToOutput, (isMated == true) ? '2' : 'u');
  }

  errno = 0;
  Ao = fopen(path1, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", path1, strerror(errno)), exit(1);

  if (isMated == true) {
    errno = 0;
    Bo = fopen(path2, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", path2, strerror(errno)), exit(1);
  }

  uint64 i=0;
  uint64 s=0;


  if (isMated == true) {

    if (nPairsToOutput > 0)
      fprintf(stderr, "Extracting "F_U64" mate pairs into %s and %s\n",
              nPairsToOutput, path1, path2);
    else
      fprintf(stderr, "Extracting "F_U64" bases of mate pairs into %s and %s\n",
              nBasesToOutput, path1, path2);

    for (; Ar->read(Ai) && Br->read(Bi); i++) {
      if ((i < totPairsInInput) && (sav[i])) {
        Ar->write(Ao);
        Br->write(Bo);
        s++;
      }
    }

    fclose(Ai);
    fclose(Bi);

    fclose(Ao);
    fclose(Bo);

  } else {

    if (nPairsToOutput > 0)
      fprintf(stderr, "Extracting "F_U64" reads into %s\n",
              nPairsToOutput, path1);
    else
      fprintf(stderr, "Extracting "F_U64" bases of reads into %s\n",
              nBasesToOutput, path1);

    for (; Ar->read(Ai); i++) {
      if ((i < totPairsInInput) && (sav[i])) {
        Ar->write(Ao);
        s++;
      }
    }

    fclose(Ai);
    fclose(Ao);
  }

  delete Ar;
  delete Br;

  if (i > totPairsInInput) {
    fprintf(stderr, "WARNING:  There are "F_U64" %s in the input; you claimed there are "F_U64" (-t option) %s.\n",
            i,               (isMated) ? "mates" : "reads",
            totPairsInInput, (isMated) ? "mates" : "reads");
    fprintf(stderr, "WARNING:  Result is not a random sample of the input file.\n");
  }

  if (i < totPairsInInput) {
    fprintf(stderr, "WARNING:  There are "F_U64" %s in the input; you claimed there are "F_U64" (-t option) %s.\n",
            i,               (isMated) ? "mates" : "reads",
            totPairsInInput, (isMated) ? "mates" : "reads");
    fprintf(stderr, "WARNING:  Result is only %f X coverage.\n", (double)s * READLENGTH / GENOMESIZE);
  }

  return(0);
}
