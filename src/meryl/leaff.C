
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
 *    kmer/leaff/leaff.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2003-JAN-02 to 2003-OCT-14
 *      are Copyright 2003 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-FEB-20 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-MAR-06 to 2014-APR-11
 *      are Copyright 2005-2009,2011-2012,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Liliana Florea on 2011-NOV-16
 *      are Copyright 2011 Liliana Florea, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-AUG-22 to 2015-JAN-13
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-DEC-07
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "seqCache.H"
#include "seqStore.H"

#include "md5.H"
#include "mt19937ar.H"
#include "dnaAlphabets.H"

//  Analysis functions
//
void dumpBlocks(char *filename);
void stats(char *filename, uint64 refLen);
void partitionBySize(char *prefix, uint64 partitionSize, char *filename);
void partitionByBucket(char *prefix, uint64 partitionSize, char *filename);
void partitionBySegment(char *prefix, uint64 numSegments, char *filename);
void simseq(char *,char *, uint32, uint32, uint32, uint32, double);
void computeGCcontent(char *name);
void findDuplicates(char *filename);
void mapDuplicates(char *filea, char *fileb);

void processFile(char  *filename);
void processArray(int argc, char **argv);

bool                   doReverse         = false;
bool                   doComplement      = false;
bool                   withDefLine       = true;
char                  *specialDefLine    = 0L;
uint32                 withLineBreaks    = 0;

bool                   toUppercase       = false;
char                   translate[256]    = {0};

seqCache              *fasta             = 0L;

uint32                 begPos            =  (uint32)0;
uint32                 endPos            = ~(uint32)0;

uint32                 endExtract        = ~(uint32)0;

mtRandom               MT;


static
void
failIfNoSource(void) {
  if (fasta == 0L)
    fprintf(stderr, "No source file specified.\n"), exit(1);
}

static
void
failIfNotRandomAccess(void) {
  if (fasta->randomAccessSupported() == false)
    fprintf(stderr, "Algorithm required random access; soruce file not supported.\n"), exit(1);
}


static
void
helpStandard(char *program) {
    fprintf(stderr, "usage: %s [-f fasta-file] [options]\n", program);
    fprintf(stderr, "\n");
    fprintf(stderr, "SOURCE FILES\n");
    fprintf(stderr, "   -f file:     use sequence in 'file' (-F is also allowed for historical reasons)\n");
    fprintf(stderr, "   -A file:     read actions from 'file'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "SOURCE FILE EXAMINATION\n");
    fprintf(stderr, "   -d:          print the number of sequences in the fasta\n");
    fprintf(stderr, "   -i name:     print an index, labelling the source 'name'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OUTPUT OPTIONS\n");
    fprintf(stderr, "   -6 <#>:      insert a newline every 60 letters\n");
    fprintf(stderr, "                  (if the next arg is a number, newlines are inserted every\n");
    fprintf(stderr, "                  n letters, e.g., -6 80.  Disable line breaks with -6 0,\n");
    fprintf(stderr, "                  or just don't use -6!)\n");
    fprintf(stderr, "   -e beg end:  Print only the bases from position 'beg' to position 'end'\n");
    fprintf(stderr, "                  (space based, relative to the FORWARD sequence!)  If\n");
    fprintf(stderr, "                  beg == end, then the entire sequence is printed.  It is an\n");
    fprintf(stderr, "                  error to specify beg > end, or beg > len, or end > len.\n");
    fprintf(stderr, "   -ends n      Print n bases from each end of the sequence.  One input\n");
    fprintf(stderr, "                  sequence generates two output sequences, with '_5' or '_3'\n");
    fprintf(stderr, "                  appended to the ID.  If 2n >= length of the sequence, the\n");
    fprintf(stderr, "                  sequence itself is printed, no ends are extracted (they\n");
    fprintf(stderr, "                  overlap).\n");
    fprintf(stderr, "   -C:          complement the sequences\n");
    fprintf(stderr, "   -H:          DON'T print the defline\n");
    fprintf(stderr, "   -h:          Use the next word as the defline (\"-H -H\" will reset to the\n");
    fprintf(stderr, "                  original defline\n");
    fprintf(stderr, "   -R:          reverse the sequences\n");
    fprintf(stderr, "   -u:          uppercase all bases\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "SEQUENCE SELECTION\n");
    fprintf(stderr, "   -G n s l:    print n randomly generated sequences, 0 < s <= length <= l\n");
    fprintf(stderr, "   -L s l:      print all sequences such that s <= length < l\n");
    fprintf(stderr, "   -N l h:      print all sequences such that l <= %% N composition < h\n");
    fprintf(stderr, "                  (NOTE 0.0 <= l < h < 100.0)\n");
    fprintf(stderr, "                  (NOTE that you cannot print sequences with 100%% N\n");
    fprintf(stderr, "                   This is a useful bug).\n");
    fprintf(stderr, "   -q file:     print sequences from the seqid list in 'file'\n");
    fprintf(stderr, "   -r num:      print 'num' randomly picked sequences\n");
    fprintf(stderr, "   -s seqid:    print the single sequence 'seqid'\n");
    fprintf(stderr, "   -S f l:      print all the sequences from ID 'f' to 'l' (inclusive)\n");
    fprintf(stderr, "   -W:          print all sequences (do the whole file)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "LONGER HELP\n");
    fprintf(stderr, "   -help analysis\n");
    fprintf(stderr, "   -help examples\n");
}


static
void
helpAnalysis(char *program) {
  fprintf(stderr, "usage: %s [-f <fasta-file>] [options]\n", program);
  fprintf(stderr, "\n");
  fprintf(stderr, "   --findduplicates a.fasta\n");
  fprintf(stderr, "                Reports sequences that are present more than once.  Output\n");
  fprintf(stderr, "                is a list of pairs of deflines, separated by a newline.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   --mapduplicates a.fasta b.fasta\n");
  fprintf(stderr, "                Builds a map of IIDs from a.fasta and b.fasta that have\n");
  fprintf(stderr, "                identical sequences.  Format is \"IIDa <-> IIDb\"\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   --md5 a.fasta:\n");
  fprintf(stderr, "                Don't print the sequence, but print the md5 checksum\n");
  fprintf(stderr, "                (of the entire sequence) followed by the entire defline.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   --partition     prefix [ n[gmk]bp | n ] a.fasta\n");
  fprintf(stderr, "   --partitionmap         [ n[gmk]bp | n ] a.fasta\n");
  fprintf(stderr, "                Partition the sequences into roughly equal size pieces of\n");
  fprintf(stderr, "                size nbp, nkbp, nmbp or ngbp; or into n roughly equal sized\n");
  fprintf(stderr, "                parititions.  Sequences larger that the partition size are\n");
  fprintf(stderr, "                in a partition by themself.  --partitionmap writes a\n");
  fprintf(stderr, "                description of the partition to stdout; --partiton creates\n");
  fprintf(stderr, "                a fasta file 'prefix-###.fasta' for each partition.\n");
  fprintf(stderr, "                Example: -F some.fasta --partition parts 130mbp\n");
  fprintf(stderr, "                         -F some.fasta --partition parts 16\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   --segment prefix n a.fasta\n");
  fprintf(stderr, "                Splits the sequences into n files, prefix-###.fasta.\n");
  fprintf(stderr, "                Sequences are not reordered.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   --gccontent a.fasta\n");
  fprintf(stderr, "                Reports the GC content over a sliding window of\n");
  fprintf(stderr, "                3, 5, 11, 51, 101, 201, 501, 1001, 2001 bp.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   --testindex a.fasta\n");
  fprintf(stderr, "                Test the index of 'file'.  If index is up-to-date, leaff\n");
  fprintf(stderr, "                exits successfully, else, leaff exits with code 1.  If an\n");
  fprintf(stderr, "                index file is supplied, that one is tested, otherwise, the\n");
  fprintf(stderr, "                default index file name is used.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   --dumpblocks a.fasta\n");
  fprintf(stderr, "                Generates a list of the blocks of N and non-N.  Output\n");
  fprintf(stderr, "                format is 'base seq# beg end len'.  'N 84 483 485 2' means\n");
  fprintf(stderr, "                that a block of 2 N's starts at space-based position 483\n");
  fprintf(stderr, "                in sequence ordinal 84.  A '.' is the end of sequence\n");
  fprintf(stderr, "                marker.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   --errors L N C P a.fasta\n");
  fprintf(stderr, "                For every sequence in the input file, generate new\n");
  fprintf(stderr, "                sequences including simulated sequencing errors.\n");
  fprintf(stderr, "                L -- length of the new sequence.  If zero, the length\n");
  fprintf(stderr, "                     of the original sequence will be used.\n");
  fprintf(stderr, "                N -- number of subsequences to generate.  If L=0, all\n");
  fprintf(stderr, "                     subsequences will be the same, and you should use\n");
  fprintf(stderr, "                     C instead.\n");
  fprintf(stderr, "                C -- number of copies to generate.  Each of the N\n");
  fprintf(stderr, "                     subsequences will have C copies, each with different\n");
  fprintf(stderr, "                     errors.\n");
  fprintf(stderr, "                P -- probability of an error.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "                HINT: to simulate ESTs from genes, use L=500, N=10, C=10\n");
  fprintf(stderr, "                         -- make C=10 sequencer runs of N=10 EST sequences\n");
  fprintf(stderr, "                            of length 500bp each.\n");
  fprintf(stderr, "                      to simulate mRNA from genes, use L=0, N=10, C=10\n");
  fprintf(stderr, "                      to simulate reads from genomes, use L=800, N=10, C=1\n");
  fprintf(stderr, "                         -- of course, N= should be increased to give the\n");
  fprintf(stderr, "                            appropriate depth of coverage\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   --stats a.fasta [refLen]\n");
  fprintf(stderr, "                Reports size statistics; number, N50, sum, largest.\n");
  fprintf(stderr, "                If 'refLen' is supplied, N50 is based on this size.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   --seqstore out.seqStore\n");
  fprintf(stderr, "                Converts the input file (-f) to a seqStore file.\n");
}


static
void
helpExamples(char *program) {
  fprintf(stderr, "usage: %s [-f <fasta-file>] [options]\n", program);
  fprintf(stderr, "\n");
  fprintf(stderr, "Options are ORDER DEPENDENT.  Sequences are printed whenever an ACTION occurs\n");
  fprintf(stderr, "on the command line.  SEQUENCE OPTIONS are not reset when a sequence is printed.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "SEQUENCES are numbered starting at ZERO, not one.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   Print the first 10 bases of the fourth sequence in file 'genes':\n");
  fprintf(stderr, "       -f genes -e 0 10 -s 3\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   Print the first 10 bases of the fourth and fifth sequences:\n");
  fprintf(stderr, "       -f genes -e 0 10 -s 3 -s 4\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   Print the fourth and fifth sequences reverse complemented, and the sixth\n");
  fprintf(stderr, "   sequence forward.  The second set of -R -C toggle off reverse-complement:\n");
  fprintf(stderr, "       -f genes -R -C -s 3 -s 4 -R -C -s 5\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   Convert file 'genes' to a seqStore 'genes.seqStore'.  The seqStore\n");
  fprintf(stderr, "   provides better performance with the kmer tools.\n");
  fprintf(stderr, "       -f genes --seqstore genes.seqStore\n");
}


static
void
printSequence(char        *def,
              char        *seq,
              uint32       beg,
              uint32       end) {

  if (beg >= end)
    return;

  if ((endExtract != ~uint32ZERO) &&
      (endExtract + endExtract < end - beg)) {
    char    d[1024];
    uint32  l = strlen(seq);

    sprintf(d, "%s_5", def);
    printSequence(d, seq, 0, endExtract);

    sprintf(d, "%s_3", def);
    printSequence(d, seq, l-endExtract, l);

    return;
  }

  if (specialDefLine)
    def = specialDefLine;

  if (withDefLine == false)
    def = 0L;

  uint32    limit = end - beg;
  char     *n = new char [end - beg + 1];
  char     *m;

  if        ((doReverse == false) && (doComplement == false)) {
    m    = n;
    seq += beg;
    while (limit--)
      *(m++) = translate[*(seq++)];

  } else if ((doReverse == true) && (doComplement == false)) {
    m    = n + limit - 1;
    seq += beg;
    while (limit--)
      *(m--) = translate[*(seq++)];

  } else if ((doReverse == false) && (doComplement == true)) {
    m    = n;
    seq += beg;
    while (limit--)
      *(m++) = alphabet.complementSymbol(translate[*(seq++)]);

  } else if ((doReverse == true) && (doComplement == true)) {
    m    = n + limit - 1;
    seq += beg;
    while (limit--)
      *(m--) = alphabet.complementSymbol(translate[*(seq++)]);
  }

  n[end-beg] = 0;

  if (def)
    fprintf(stdout, ">%s\n", def);

  if (withLineBreaks) {
    char      *t = n;
    char      *a = new char [withLineBreaks+1];

    while (*t) {
      uint32 i=0;
      while ((*t) && (i < withLineBreaks))
        a[i++] = *(t++);
      a[i++] = '\n';
      a[i]   = 0;
      fprintf(stdout, "%s", a);
    }

    delete [] a;
  } else {
    fprintf(stdout, "%s\n", n);
  }

  delete [] n;
}


static
void
printSequence(seqInCore *sic) {
  printSequence(sic->header(), sic->sequence(), (begPos!=(uint32)0) ? begPos:0, (endPos!=~uint32(0)) ? endPos:sic->sequenceLength());
}


static
void
printSequence(uint32 sid) {
  seqInCore *sic   = fasta->getSequenceInCore(sid);
  if (sic == 0L)
    fprintf(stderr, "WARNING: Didn't find sequence with iid '"F_U32"'\n", sid);
  else
    printSequence(sic);
  delete sic;
}


static
void
printSequence(char *sid) {
  seqInCore  *sic = fasta->getSequenceInCore(sid);
  if (sic == 0L)
    fprintf(stderr, "WARNING: Didn't find sequence with name/iid '%s'\n", sid);
  else
    printSequence(sic);
  delete sic;
}


static
void
printIDsFromFile(char *name) {
  uint32      idLen = 0;
  uint32      idMax = 63;
  char       *id    = new char [idMax+1];

  readBuffer  B(name);
  char        x = B.read();

  //  For optimal performance, we should sort the list of ID's given
  //  by their IID, but the user might have a good reason for wanting
  //  them unsorted.

  while (B.eof() == false) {
    while (alphabet.isWhitespace(x) && (B.eof() == false))
      x = B.read();

    if (B.eof() == false) {
      idLen = 0;

      while (!alphabet.isWhitespace(x) && (B.eof() == false)) {
        id[idLen++] = x;
        x = B.read();

        if (idLen >= idMax) {
          idMax *= 2;
          char *newid = new char [idMax+1];
          memcpy(newid, id, sizeof(char) * idLen);
          delete [] id;
          id = newid;
        }
      }

      id[idLen] = 0;

      seqInCore  *S = fasta->getSequenceInCore(id);

      if (S == 0L)
        fprintf(stderr, "WARNING: Didn't find sequence with name/iid '%s'\n", id);
      else
        printSequence(S);
    }
  }

  delete [] id;
}


void
processArray(int argc, char **argv) {

  int arg = 1;
  while (arg < argc) {

    if       ((strcmp(argv[arg], "-f") == 0) ||
              (strcmp(argv[arg], "-F") == 0)) {
      delete fasta;
      fasta = new seqCache(argv[++arg]);

    } else if (strcmp(argv[arg], "-i") == 0) {

      failIfNoSource();

      ++arg;
      if ((argv[arg] == 0L) || (argv[arg][0] == '-'))
        fprintf(stderr, "ERROR: next arg to -i should be 'name', I got '%s'\n",
                (argv[arg] == 0L) ? "(nullpointer)" : argv[arg]), exit(1);

      for (uint32 s=0; s<fasta->getNumberOfSequences(); s++)
        fprintf(stdout, "G\tseq\t%s:"F_U32"\t"F_U32"\t%s\n",
                argv[arg], s, fasta->getSequenceLength(s), ">unimplemented");

    } else if (strcmp(argv[arg], "-d") == 0) {
      failIfNoSource();
      printf(F_U32"\n", fasta->getNumberOfSequences());

    } else if (strcmp(argv[arg], "-L") == 0) {
      uint32 small = strtouint32(argv[++arg]);
      uint32 large = strtouint32(argv[++arg]);

      failIfNoSource();

      for (uint32 s=0; s<fasta->getNumberOfSequences(); s++)
        if ((small <= fasta->getSequenceLength(s)) && (fasta->getSequenceLength(s) < large))
          printSequence(s);

    } else if (strcmp(argv[arg], "-N") == 0) {
      double small = atof(argv[++arg]);
      double large = atof(argv[++arg]);

      failIfNoSource();

      for (uint32 s=0; s<fasta->getNumberOfSequences(); s++) {
        seqInCore *S   = fasta->getSequenceInCore(s);
        uint32     Ns  = 0;
        uint32     len = S->sequenceLength();
        char      *seq = S->sequence();

        for (uint32 i=begPos; i<len && i<endPos; i++)
          if ((seq[i] == 'n') || (seq[i] == 'N'))
            Ns++;

        double Np = 100.0 * Ns / len;

        if ((small <= Np) && (Np < large))
          printSequence(S);

        delete S;
      }

    } else if (strcmp(argv[arg], "-W") == 0) {
      failIfNoSource();

      for (uint32 s=0; s<fasta->getNumberOfSequences(); s++)
        printSequence(s);

    } else if (strcmp(argv[arg], "-G") == 0) {
      uint32 n = strtouint32(argv[++arg]);
      uint32 s = strtouint32(argv[++arg]);
      uint32 l = strtouint32(argv[++arg]);

      char      bases[4] = {'A', 'C', 'G', 'T'};
      char     *def      = new char [1024];
      char     *seq      = new char [l + 1];

      if (s == 0)
        s = 1;
      if (s > l)
        fprintf(stderr, "leaff: usage: -G num-seqs min-length max-length\n"), exit(1);

      for (uint32 i=0; i<n; i++) {
        uint32 j = s + ((l-s == 0) ? 0 : (MT.mtRandom32() % (l-s)));
        uint32 p = 0;

        while (p < j)
          seq[p++] = bases[MT.mtRandom32() & 0x3];
        seq[p] = 0;

        sprintf(def, "random%06"F_U32P, i);

        printSequence(def, seq, 0, j);
      }

      delete [] seq;
      delete [] def;

    } else if (strcmp(argv[arg], "-s") == 0) {
      failIfNoSource();
      failIfNotRandomAccess();  //  Easy to fix, just read the first N sequences
      printSequence(argv[++arg]);

    } else if (strcmp(argv[arg], "-S") == 0) {
      failIfNoSource();
      failIfNotRandomAccess();  //  Easy to fix, just read the first N sequences

      uint32 lowID  = fasta->getSequenceIID(argv[++arg]);
      uint32 highID = fasta->getSequenceIID(argv[++arg]);

      if (lowID > highID) {
        uint32 t = lowID;
        lowID    = highID;
        highID   = t;
      }

      for (uint32 s=lowID; (s <= highID) && (s <= fasta->getNumberOfSequences()); s++)
        printSequence(s);

    } else if (strcmp(argv[arg], "-r") == 0) {
      uint32 num = strtouint32(argv[++arg]);

      failIfNoSource();
      failIfNotRandomAccess();  //  Impossible to fix, or load whole thing into memory

      if (num >= fasta->getNumberOfSequences())
        num = fasta->getNumberOfSequences();

      uint32  *seqs = new uint32 [fasta->getNumberOfSequences()];

      for (uint32 i=0; i<fasta->getNumberOfSequences(); i++)
        seqs[i] = i;

      for (uint32 i=0; i<fasta->getNumberOfSequences(); i++) {
        uint32 j = MT.mtRandom32() % (fasta->getNumberOfSequences() - i) + i;
        uint32 t = seqs[j];
        seqs[j] = seqs[i];
        seqs[i] = t;
      }

      for (uint32 i=0; i<num; i++)
        printSequence(seqs[i]);

      delete [] seqs;

    } else if (strcmp(argv[arg], "-q") == 0) {
      failIfNoSource();
      failIfNotRandomAccess();  //  Impossible to fix, or load whole thing into memory
      printIDsFromFile(argv[++arg]);

    } else if (strcmp(argv[arg], "-6") == 0) {
      withLineBreaks = 60;
      if ((argv[arg+1] != 0L) && (argv[arg+1][0] != '-'))
        withLineBreaks = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-w") == 0) {
      toUppercase = !toUppercase;
      for (int z=0; z<256; z++)
        translate[z] = (toUppercase) ? alphabet.toUpper(z) : (char)z;

    } else if (strcmp(argv[arg], "-R") == 0) {
      doReverse = !doReverse;

    } else if (strcmp(argv[arg], "-C") == 0) {
      doComplement = !doComplement;

    } else if (strcmp(argv[arg], "-H") == 0) {
      withDefLine    = !withDefLine;
      specialDefLine = 0L;

    } else if (strcmp(argv[arg], "-h") == 0) {
      withDefLine    = true;
      specialDefLine = argv[++arg];

    } else if (strcmp(argv[arg], "-e") == 0) {
      begPos = strtouint32(argv[++arg]);
      endPos = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-ends") == 0) {
      endExtract = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-A") == 0) {
      processFile(argv[++arg]);



    } else if (strcmp(argv[arg], "--findduplicates") == 0) {
      findDuplicates(argv[++arg]);
      exit(0);

    } else if (strcmp(argv[arg], "--mapduplicates") == 0) {
      mapDuplicates(argv[arg+1], argv[arg+2]);
      exit(0);

    } else if (strcmp(argv[arg], "--md5") == 0) {
      md5_s     md5;
      char      sum[33];

      fasta = new seqCache(argv[++arg]);

      for (uint32 s=0; s<fasta->getNumberOfSequences(); s++) {
        seqInCore *S = fasta->getSequenceInCore(s);
        fprintf(stdout, "%s %s\n",
                md5_toascii(md5_string(&md5, S->sequence(), S->sequenceLength()), sum),
                S->header());
        delete S;
      }
      delete fasta;
      exit(0);

    } else if ((strcmp(argv[arg], "--partition") == 0) ||
               (strcmp(argv[arg], "--partitionmap") == 0)) {

      char *prefix = 0L;
      if (strcmp(argv[arg], "--partition") == 0)
        prefix = argv[++arg];

      //  does the next arg end with gbp, mbp, kbp or bp?  If so,
      //  partition by length, else partition into buckets.
      //
      int     al = strlen(argv[arg+1]);
      uint64  ps = strtouint64(argv[arg+1]);

      char a3 = (al<3) ? '0' : alphabet.toLower(argv[arg+1][al-3]);
      char a2 = (al<2) ? '0' : alphabet.toLower(argv[arg+1][al-2]);
      char a1 = (al<1) ? '0' : alphabet.toLower(argv[arg+1][al-1]);

      //  partition!

      if (!isdigit(a1) || !isdigit(a2) || !isdigit(a3)) {
        if        ((a3 == 'g') && (a2 == 'b') && (a1 == 'p')) {
          ps *= 1000000000;
        } else if ((a3 == 'm') && (a2 == 'b') && (a1 == 'p')) {
          ps *= 1000000;
        } else if ((a3 == 'k') && (a2 == 'b') && (a1 == 'p')) {
          ps *= 1000;
        } else if (isdigit(a3) && (a2 == 'b') && (a1 == 'p')) {
          ps *= 1;
        } else {
          fprintf(stderr, "Unknown partition size option '%s'\n", argv[arg+1]), exit(1);
        }

        if (ps == 0)
          fprintf(stderr, "Unknown or zero partition size '%s'\n", argv[arg+1]), exit(1);
        partitionBySize(prefix, ps, argv[arg+2]);
      } else {
        if (ps == 0)
          fprintf(stderr, "Unknown or zero partition size '%s'\n", argv[arg+1]), exit(1);
        partitionByBucket(prefix, ps, argv[arg+2]);
      }
      exit(0);

    } else if (strcmp(argv[arg], "--segment") == 0) {
      partitionBySegment(argv[arg+1], strtouint32(argv[arg+2]), argv[arg+3]);
      exit(0);

    } else if (strcmp(argv[arg], "--gccontent") == 0) {
      computeGCcontent(argv[++arg]);
      exit(0);

    } else if (strcmp(argv[arg], "--dumpblocks") == 0) {
      dumpBlocks(argv[++arg]);
      exit(0);

    } else if (strcmp(argv[arg], "--stats") == 0) {
      stats(argv[arg+1], (argv[arg+2] != 0L) ? strtouint64(argv[arg+2]) : 0);
      exit(0);

    } else if (strcmp(argv[arg], "--errors") == 0) {
      uint32  L = strtouint32(argv[++arg]);     //  Desired length
      uint32  l = 0;                            //  min of desired length, length of sequence
      uint32  N = strtouint32(argv[++arg]);     //  number of copies per sequence
      uint32  C = strtouint32(argv[++arg]);     //  number of mutations per copy
      double  P = atof(argv[++arg]);            //  probability of mutation
      uint32  i = 0;

      fasta = new seqCache(argv[++arg]);

      seqInCore *S = fasta->getSequenceInCore(i++);
      while (S) {
        char   *seq = S->sequence();
        char   *hdr = S->header();
        uint32  len = S->sequenceLength();

        l = len;
        if ((L > 0) && (L < len))
          l = L;

        simseq(seq, hdr, len, N, l, C, P);

        delete S;
        S = fasta->getSequenceInCore(i++);
      }
      delete fasta;
      exit(0);

    } else if (strcmp(argv[arg], "--seqstore") == 0) {
      constructSeqStore(argv[++arg], fasta);
      exit(0);

    } else if (strcmp(argv[arg], "-help") == 0) {
      if      ((argv[arg+1]) && (strcmp(argv[arg+1], "analysis") == 0))
        helpAnalysis(argv[0]);
      else if ((argv[arg+1]) && (strcmp(argv[arg+1], "examples") == 0))
        helpExamples(argv[0]);
      else
        helpStandard(argv[0]);
      exit(0);

    } else {
      helpStandard(argv[0]);
      fprintf(stderr, "Unknown option '%s'\n", argv[arg]);
      exit(1);
    }

    arg++;
  }

  delete fasta;
  fasta = 0L;
}


void
processFile(char  *filename) {
  FILE   *F    = NULL;

  if (strcmp(filename, "-") == 0) {
    F = stdin;
  } else {
    errno = 0;
    F = fopen(filename, "r");
    if (errno)
      fprintf(stderr, "Couldn't open '%s': %s\n", filename, strerror(errno)), exit(1);
  }

  uint64  max  = 16 * 1024 * 1024;
  uint64  pos  = 0;
  size_t  len  = 0;

  char   *data = new char [max];

  //  Suck the file into 'data'

  while (!feof(F)) {
    errno = 0;
    len = fread(data+pos, 1, max - pos, F);
    if (errno)
      fprintf(stderr, "Couldn't read "F_U64" bytes from '%s': %s\n",
              (uint64)(max-pos), filename, strerror(errno)), exit(1);

    pos += len;

    if (pos >= max) {
      max += 16 * 1024 * 1024;
      char *tmpd = new char [max];
      memcpy(tmpd, data, pos);
      delete [] data;
      data = tmpd;
    }
  }

  if (strcmp(filename, "-") != 0)
    fclose(F);

  len = pos;

  //  (over)count the number of words; we start at two, since the
  //  first arg is the name of the program, and if there is only one
  //  word and no whitespace in the file, the below loop fails to
  //  count the second word.

  int     argc = 2;
  char  **argv = 0L;

  for (uint32 i=0; i<len; i++) {
    if (isspace(data[i])) {
      argc++;
      data[i] = 0;
    }
  }

  //  Allocate space for word pointers, then set them.  First arg in
  //  argv[] is the name of the program -- we'll substitute the name
  //  of the file instead.

  argv = new char * [argc];

  argv[0] = filename;
  argc    = 1;

  //  Three steps: Skip leading whitespace; save the arg if it's a
  //  real arg (and not the end of the file; then skip the word.

  for (uint32 pos=0; pos<len; pos++) {
    while ((data[pos] == 0) && (pos < len))
      pos++;
    if (pos < len)
      argv[argc++] = data+pos;
    while ((data[pos] != 0) && (pos < len))
      pos++;
  }

  processArray(argc, argv);

  delete [] argv;
  delete [] data;
}


int
main(int argc, char **argv) {

  if (argc < 2) {
    helpStandard(argv[0]);
    exit(1);
  }

  for (int z=0; z<256; z++)
    translate[z] = (char)z;

  processArray(argc, argv);

  delete fasta;
}




