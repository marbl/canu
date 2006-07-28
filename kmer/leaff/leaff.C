#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>

//  Linux needs to include time.h; others can use sys/time.h.
//  Tru64 is ok with time.h
//
#include <time.h>

#include "bio++.H"

void          simseq(char *,char *,int,int,int,int,double);

#include "buildinfo-leaff.h"
#include "buildinfo-libbio.h"
#include "buildinfo-libutil.h"


const char *usage =
"usage: %s [--buildinfo] [-f|-F|-I <fasta-file>] [options]\n"
"\n"
"ENTERTAINMENT OPTIONS\n"
"       -V:          whenever a sequence is read, print the defline to stderr\n"
"\n"
"SOURCE FILE\n"
"       -f file:     use 'file' as an UN-INDEXED source file\n"
"       -Ft[c] file: use 'file' as an INDEXED source file (the index is built if\n"
"                    it doesn't exist), where 't' is the type of index to build:\n"
"                      i:  internal id's only (the default if t is not specified)\n"
"                      n:  names only (the first word on defline)\n"
"                      d:  full deflines\n"
"                    The optional character c will build an index with checksums,\n"
"                    and --checksum will verify the checksums.\n"
"\n"
"OUTPUT OPTIONS\n"
"       -6:          insert a newline every 60 letters\n"
"                      (if the next arg is a number, newlines are inserted every\n"
"                       n letters, e.g., -6 80.  Disable line breaks with -6 0,\n"
"                       or just don't use -6!)\n"
"       -u:          uppercase all bases\n"
"       -R:          print the sequence reversed\n"
"       -C:          print the sequence complemented\n"
"       -H:          DON'T print the defline\n"
"       -h:          Use the next word as the defline\n"
"                      (use \"-H -H\" to resume using the original defline)\n"
"       -e beg end:  Print only the bases from position 'beg' to position 'end'\n"
"                    (space based, relative to the FORWARD sequence!)  If\n"
"                    beg == end, then the entire sequence is printed.  It is an\n"
"                    error to specify beg > end, or beg > len, or end > len.\n"
"       -md5:        Don't print the sequence, but print the md5 checksum\n"
"                    (of the entire sequence) followed by the entire defline.\n"
"                    (this should really be an ANALYSIS OPTION)\n"
"\n"
"SEQUENCE SELECTION (no index needed)\n"
"       -L s l:      print all sequences such that s <= length < l\n"
"       -N l h:      print all sequences such that l <= % N composition < h\n"
"                      (NOTE 0.0 <= l < h < 100.0)\n"
"                      (NOTE that you cannot print sequences with 100% N\n"
"                       This is a useful bug).\n"
"       -G n s l:    print n randomly generated sequences, 0 < s <= length <= l\n"
"       -W:          print all sequences (do the whole file)\n"
"\n"
"SEQUENCE SELECTION (index needed)\n"
"       -s seqid:    print the single sequence 'seqid'\n"
"       -Ii:         'seqid' is an internal id (an integer), this is the default\n"
"                      for -Fi\n"
"       -Ie:         'seqid' is an external id (a word), this is the default for\n"
"                      -Fn or -Fd, and is an error for -Fi\n"
"       -d:          print the number of sequences in the fasta\n"
"       -i name:     print an ATA compliant index, labelling the source 'name'\n"
"       -ii:         print the index in an almost-human readable format\n"
"       -S f l:      print all the sequences from IID 'f' to 'l' (inclusive)\n"
"       -r num:      print 'num' randomly picked sequences\n"
"       -q file:     print sequences from the seqid list in 'file'\n"
"\n"
"ANALYSIS OPTIONS\n"
"\n"
"       --findduplicates a.fasta\n"
"                    Reports sequences that are present more than once.  Output\n"
"                    is a list of pairs of deflines, separated by a newline.\n"
"\n"
"       --mapduplicates a.fasta b.fasta\n"
"                    Builds a map of IIDs from a.fasta and b.fasta that have\n"
"                    identical sequences.  Format is \"IIDa <-> IIDb\"\n"
"\n"
"       --partition prefix [ n[gmk]bp | n ]\n"
"       --partitionmap [ n[gmk]bp | n ]\n"
"                    Partition the sequences into roughly equal size pieces of\n"
"                    size nbp, nkbp, nmbp or ngbp; or into n roughly equal sized\n"
"                    parititions.  Sequences larger that the partition size are\n"
"                    in a partition by themself.  --partitionmap writes a\n"
"                    description of the partition to stdout; --partiton creates\n"
"                    a fasta file 'prefix-###.fasta' for each partition.\n"
"                    Example: -F some.fasta --partition parts 130mbp\n"
"                             -F some.fasta --partition parts 16\n"
"\n"
"       --checksum a.fasta\n"
"                    One of three actions:\n"
"                    1) If no fastaidx file exists, this is equivalent to\n"
"                       \"-Fc a.fasta\".\n"
"                    2) If a fastaidx file exists, but doesn't have checksums,\n"
"                       checksums are added to that index file.\n"
"                    3) If a fastaidx file exists and checksums are present,\n"
"                       the checksums are verified.\n"
"\n"
"       --testindex a.fasta\n"
"                    Test the index of 'file'.  If index is up-to-date, leaff\n"
"                    exits successfully, else, leaff exits with code 1.  If an\n"
"                    index file is supplied, that one is tested, otherwise, the\n"
"                    default index file name is used.\n"
"\n"
"       --dumpblocks\n"
"                    Generates a list of the blocks of N and non-N.  Output\n"
"                    format is 'base seq# beg end len'.  'N 84 483 485 2' means\n"
"                    that a block of 2 N's starts at space-based position 483\n"
"                    in sequence ordinal 84.  A '.' is the end of sequence\n"
"                    marker.\n"
"\n"
"       --errors L N C P\n"
"                    For every sequence in the input file, generate new\n"
"                    sequences including simulated sequencing errors.\n"
"                    L -- length of the new sequence.  If zero, the length\n"
"                         of the original sequence will be used.\n"
"                    N -- number of subsequences to generate.  If L=0, all\n"
"                         subsequences will be the same, and you should use\n"
"                         C instead.\n"
"                    C -- number of copies to generate.  Each of the N\n"
"                         subsequences will have C copies, each with different\n"
"                         errors.\n"
"                    P -- probability of an error.\n"
"\n"
"                    HINT: to simulate ESTs from genes, use L=500, N=10, C=10\n"
"                          to simulate mRNA from genes, use L=0, N=10, C=10\n"
"\n"
"EXPERT OPTIONS\n"
"       -A:  Read actions from 'file'\n"
"\n"
"EXAMPLES\n"
"       Options are ORDER DEPENDENT.  Sequences are printed whenever an\n"
"       ACTION occurs on the command line.  SEQUENCE OPTIONS are not reset\n"
"       when a sequence is printed.\n"
"\n"
"       -F file -e 0 10 -s 3       Print the first 10 bases of the fourth\n"
"                                  sequence in file.\n"
"       -F file -e 0 10 -s 3 -s 4  Print the first 10 bases of the fourth\n"
"                                  and fifth sequences in file.\n"
"       -F file -R -C -s 3 -s 4    Print the fourth and fifth sequences\n"
"                                  reverse complemented.\n"
"\n";


//  This will do our character translation (e.g., to upper)
//
char translate[256];

void  processFile(char  *filename);
void  processArray(int argc, char **argv);

void
printSequence(FastASequenceInCore *b,
              u32bit               beg,
              u32bit               end,
              u32bit               withLineBreaks=0,
              bool                 reverse=false,
              bool                 complement=false) {
  u32bit           style  = 0;
  if (reverse)     style += 1;
  if (complement)  style += 2;

  u32bit l = b->sequenceLength();

  if (beg == end) {
    beg = 0;
    end = l;
  }

#if 0
  if ((beg > l) || (end > l) || (beg > end)) {
    fprintf(stderr, "WARNING:  Printing "u32bitFMT" to "u32bitFMT" of sequence "u32bitFMT" (len="u32bitFMT") is out of bounds -- NO SEQUENCE PRINTED!\n",
            beg, end, b->getIID(), l);
    return;
  }
#else
  if (beg > l)
    beg = 0;
  if (end > l)
    end = l;
#endif

  u32bit    limit = end - beg;
  char     *n = new char [end - beg + 1];
  char     *m;
  char     *s = b->sequence();

  switch (style) {
    case 0:
      //  forward
      m  = n;
      s += beg;

      while (limit--)
        *(m++) = translate[(int)*(s++)];
      break;
    case 1:
      //  reverse
      m  = n + limit - 1;
      s += beg;

      while (limit--)
        *(m--) = translate[(int)*(s++)];
      break;
    case 2:
      //  complement
      m  = n;
      s += beg;

      while (limit--)
        *(m++) = complementSymbol[(int)translate[(int)*(s++)]];
      break;
    case 3:
      //  reverse complement
      m  = n + limit - 1;
      s += beg;

      while (limit--)
        *(m--) = complementSymbol[(int)translate[(int)*(s++)]];
      break;
  }

  n[end-beg] = 0;

  if (withLineBreaks) {
    char      *t = n;
    char      *b = new char [withLineBreaks+1];

    while (*t) {
      u32bit i=0;
      while ((*t) && (i < withLineBreaks))
        b[i++] = *(t++);
      b[i++] = '\n';
      b[i]   = 0;
      fprintf(stdout, "%s", b);
    }

    delete [] b;
  } else {
    fprintf(stdout, "%s\n", n);
  }

  delete [] n;
}


int
comp(const void *a, const void *b) {
  const u32bit A = *((const u32bit *)a);
  const u32bit B = *((const u32bit *)b);
  if (A < B)
    return(-1);
  if (A > B)
    return(1);
  return(0);
}


//  Global stuff for processArray
//
//  XXX:  ideally, this should be a structure, so
//  that the user can select between "new state" and
//  "old state"
//
bool                   printOnLoad       = false;
bool                   reverse           = false;
bool                   complement        = false;
bool                   withDefLine       = true;
char                  *specialDefLine    = 0L;
u32bit                 withLineBreaks    = 0;
bool                   toUppercase       = false;
char                  *sourceFile        = 0L;
FastAWrapper          *fasta             = 0L;
FastACache            *cache             = 0L;
char                   seqIDtype         = 'i';
u32bit                 begPos            = ~(u32bit)0;
u32bit                 endPos            = ~(u32bit)0;
bool                   printMD5          = false;
FastASequenceInCore   *lastSeq           = 0L;
mt_s                  *mtctx             = 0L;

void
failIfNoSource(void) {
  if (fasta == 0L) {
    fprintf(stderr, "No source file specified.\n");
    exit(1);
  }
}

void
failIfNotRandomAccess(void) {
  if (fasta->isRandomAccess() == false) {
    fprintf(stderr, "No index exists for %s\n", sourceFile);
    exit(1);
  }
}



//  Just so we can support the printOnLoad flag, this
//  routine will call loadSequence on a FastAWrapper.
//
FastASequenceInCore *loadSequence(FastAWrapper *F) {
  FastASequenceInCore *s = F->getSequence();

  if (printOnLoad)
    fprintf(stderr, "%s\n", s->header());

  return(s);
}




FastAWrapper*
openNewFile(char *name, char *arg) {

  if (fasta)
    delete fasta;

  fasta = new FastAWrapper(name);

  seqIDtype = 'i';

  if (arg[1] == 'F') {
    u32bit  indextype = FASTA_INDEX_ONLY;
    u32bit  md5type   = 0;

    for (int ap=2; arg[ap]; ap++) {
      if        (arg[ap] == 'i') {
        seqIDtype = 'i';
        indextype = FASTA_INDEX_ONLY;
      } else if (arg[ap] == 'n') {
        seqIDtype = 'e';
        indextype = FASTA_INDEX_PLUS_IDS;
      } else if (arg[ap] == 'd') {
        seqIDtype = 'e';
        indextype = FASTA_INDEX_PLUS_DEFLINES;
      } else if (arg[ap] == 'c') {
        md5type   = FASTA_INDEX_MD5;
      } else {
        fprintf(stderr, "Unknown option for '-F': '%c'\n", arg[ap]);
      }
    }

    fasta->openIndex(indextype | md5type);
  }

  delete lastSeq;
  lastSeq = 0L;

  return(fasta);
}


void
printIID(u32bit iid, FastASequenceInCore *s=0L) {
  bool  mySeq = false;

  if (s == 0L) {
    mySeq = true;
    fasta->find(iid);
    s = loadSequence(fasta);
  }

  if (printMD5) {
    md5_s     md5;
    char      sum[33];

    md5_toascii(md5_string(&md5, s->sequence(), s->sequenceLength()), sum);

    fprintf(stdout, "%s %s\n", sum, s->header());
  } else {
    if (withDefLine)
      if (specialDefLine)
        fprintf(stdout, ">%s\n", specialDefLine);
      else
        fprintf(stdout, ">%s\n", s->header()+1);

    printSequence(s, begPos, endPos, withLineBreaks, reverse, complement);
  }

  if (mySeq)
    delete s;
}



void
printSequenceBetweenSize(u32bit small, u32bit large) {

  //  XXX: Would probably be faster if we used the index

  while (!fasta->eof()) {
    FastASequenceInCore *s = loadSequence(fasta);

    if ((small <= s->sequenceLength()) && (s->sequenceLength() < large))
      printIID(0, s);

    delete s;
  }
}

void
printSequenceBetweenNComposition(double small, double large) {

  while (!fasta->eof()) {
    FastASequenceInCore *s = loadSequence(fasta);

    u32bit   Ns  = 0;
    u32bit   len = s->sequenceLength();
    char    *seq = s->sequence();
    for (u32bit i=0; i<len; i++)
      if ((seq[i] == 'n') || (seq[i] == 'N'))
        Ns++;

    double Np = 100.0 * Ns / len;

    if ((small <= Np) && (Np < large))
      printIID(0, s);

    delete s;
  }
}

void
printRandomlyGeneratedSequence(u32bit n, u32bit s, u32bit l) {
  char      bases[4] = {'A', 'C', 'G', 'T'};
  char     *seq      = new char [l + 1];

  if (s > l) {
    fprintf(stderr, "leaff: usage: -G num-seqs min-length max-length\n");
    exit(1);
  }

  if (s == 0)
    s = 1;

  for (u32bit i=0; i<n; i++) {
    u32bit j = s + (mtRandom32(mtctx) % (l-s));

    seq[j] = 0;
    while (j)
      seq[--j] = bases[mtRandom32(mtctx) & 0x3];            

    if (withDefLine)
      if (specialDefLine)
        fprintf(stdout, ">%s\n", specialDefLine);
      else
        fprintf(stdout, ">"u32bitFMT"\n", i);

    fprintf(stdout, "%s\n", seq);
  }

  delete [] seq;
}

void
findSequenceAndPrint(char *id) {

  bool found = false;

  if (seqIDtype == 'i')
    found = fasta->find(strtou32bit(id, 0L));
  else
    found = fasta->find(id);

  if (found) {
    if ((lastSeq == 0L) ||
        (lastSeq->getIID() != fasta->currentIID())) {
      delete lastSeq;
      lastSeq = loadSequence(fasta);
    }

    printIID(fasta->currentIID(), lastSeq);
  } else {
    fprintf(stderr, "WARNING: Didn't find %s id '%s'\n",
            seqIDtype == 'i' ? "internal" : "external",
            id);
  }
}

void
printRangeOfSequences(char *argl, char *argh) {
  u32bit lowID  = 0;
  u32bit highID = 0;
  bool   fail   = false;

  if (seqIDtype == 'i') {
    lowID  = strtou32bit(argl, 0L);
    if (lowID >= fasta->getNumberOfSequences()) {
      fprintf(stderr, "ERROR: Internal id of "u32bitFMT" for starting sequence is too large; only "u32bitFMT" sequences.\n",
              lowID, fasta->getNumberOfSequences());
      fail = true;
    }

    highID = strtou32bit(argh, 0L);
    if (highID >= fasta->getNumberOfSequences()) {
      fprintf(stderr, "ERROR: Internal id of "u32bitFMT" for ending sequence is too large; only "u32bitFMT" sequences.\n",
              highID, fasta->getNumberOfSequences());
      fail = true;
    }
  } else {
    if (fasta->find(argl)) {
      lowID = fasta->currentIID();
    } else {
      fprintf(stderr, "ERROR: Can't find external id '%s' of starting sequence.\n", argl);
      fail = true;
    }
    
    if (fasta->find(argh)) {
      highID = fasta->currentIID();
    } else {
      fprintf(stderr, "ERROR: Can't find external id '%s' of ending sequence.\n", argh);
      fail = true;
    }
  }

  if (fail)
    exit(1);

  //  Silently swap the ranges if needed.
  //
  if (lowID > highID) {
    u32bit t = lowID;
    lowID    = highID;
    highID   = t;
  }

  for (u32bit seqID=lowID; seqID <= highID; seqID++)
    printIID(seqID);
}


void
findAndPrintRandomSequences(u32bit num) {

  if (num >= fasta->getNumberOfSequences()) {
    fprintf(stderr, "WARNING: file has "u32bitFMT" sequences, and you asked for "u32bitFMT".\n",
            fasta->getNumberOfSequences(), num);
    printSequenceBetweenSize(u32bitZERO, ~u32bitZERO);
    return;
  }

  u32bit  *seqs = new u32bit [fasta->getNumberOfSequences()];

  for (u32bit i=0; i<fasta->getNumberOfSequences(); i++)
    seqs[i] = i;

  for (u32bit i=0; i<fasta->getNumberOfSequences(); i++) {
    u32bit j = (unsigned int)(mtRandom32(mtctx) % fasta->getNumberOfSequences());
    u32bit t = seqs[j];
    seqs[j] = seqs[i];
    seqs[i] = t;
  }

  qsort(seqs, num, sizeof(unsigned int), comp);

  for (u32bit i=0; i<num; i++)
    printIID(seqs[i]);

  delete [] seqs;
}


void
printIDsFromFile(char *name) {
  u32bit      idLen = 0;
  u32bit      idMax = 63;
  char       *id    = new char [idMax+1];

  readBuffer  B(name);

  //  For optimal performance, we should sort the list of ID's given
  //  by their IID, but the user might have a good reason for wanting
  //  them unsorted.

  while (B.eof() == false) {
    while (isspace(B.get()) && (B.eof() == false))
      B.next();

    if (B.eof() == false) {
      idLen = 0;

      while (!isspace(B.get()) && (B.eof() == false)) {
        id[idLen++] = B.get();

        if (idLen >= idMax) {
          idMax *= 2;
          char *newid = new char [idMax+1];
          memcpy(newid, id, sizeof(char) * idLen);
          delete [] id;
          id = newid;
        }

        B.next();
      }

      id[idLen] = 0;

      findSequenceAndPrint(id);
    }
  }

  delete [] id;
}


md5_s *
computeMD5ForEachSequence(FastAWrapper *F) {
  u32bit   numSeqs = fasta->getNumberOfSequences();
  md5_s   *result  = new md5_s [numSeqs];

  F->find((u32bit)0);

  for (u32bit idx=0; idx < numSeqs; idx++) {
    FastASequenceInCore *s1 = loadSequence(F);
    md5_string(result+idx, s1->sequence(), s1->sequenceLength());
    result[idx].i = s1->getIID();
    delete s1;
  }

  return(result);
}

void
findDuplicates(char *filename) {
  FastASequenceInCore  *s1 = 0L;
  FastASequenceInCore  *s2 = 0L;
  FastAWrapper         *A = new FastAWrapper(filename);

  A->openIndex(FASTA_INDEX_ONLY);

  u32bit numSeqs = A->getNumberOfSequences();

  fprintf(stderr, "Computing MD5's for each sequence in '%s'.\n", filename);
  md5_s *result = computeMD5ForEachSequence(A);

  fprintf(stderr, "Sorting MD5's.\n");
  qsort(result, numSeqs, sizeof(md5_s), md5_compare);

  fprintf(stderr, "Verifying identity, and output\n");
  for (u32bit idx=1; idx<numSeqs; idx++) {
    if (md5_compare(result+idx-1, result+idx) == 0) {
      if (result[idx-1].i == result[idx].i) {
        fprintf(stderr, "Internal error: found two copies of the same sequence iid ("u32bitFMT")!\n", result[idx].i);
        exit(1);
      }

      A->find(result[idx-1].i);
      s1 = loadSequence(A);

      A->find(result[idx].i);
      s2 = loadSequence(A);

      if (strcmp(s1->sequence(), s2->sequence()) == 0) {
        fprintf(stdout, u32bitFMT":%s\n"u32bitFMT":%s\n\n",
                result[idx-1].i, s1->header(),
                result[idx  ].i, s2->header());
      } else {
        fprintf(stderr, "COLLISION DETECTED BETWEEN IID "u32bitFMT" AND "u32bitFMT"!\nPLEASE REPORT THIS TO bri@walenz.org!\n",
                result[idx-1].i, result[idx].i);
      }

      delete s1;
      delete s2;
    }
  }

  delete [] result;
  delete    A;
}



void
mapDuplicates_Print(char *filea, FastASequenceInCore *sa,
                    char *fileb, FastASequenceInCore *sb) {
  if (strcmp(sa->sequence(), sb->sequence()) == 0) {
    fprintf(stdout, u32bitFMT" <-> "u32bitFMT"\n", sa->getIID(), sb->getIID());
  } else {
    fprintf(stderr, "COLLISION DETECTED BETWEEN %s:"u32bitFMT" AND %s:"u32bitFMT"!\nPLEASE REPORT THIS TO bri@walenz.org!\n",
            filea, sa->getIID(), fileb, sb->getIID());
  }
}



void
mapDuplicates(char *filea, char *fileb) {
  FastAWrapper   *A = new FastAWrapper(filea);
  A->openIndex(FASTA_INDEX_ONLY);
  fprintf(stderr, "Computing MD5's for each sequence in '%s'.\n", filea);
  md5_s *resultA = computeMD5ForEachSequence(A);

  FastAWrapper   *B = new FastAWrapper(fileb);
  B->openIndex(FASTA_INDEX_ONLY);
  fprintf(stderr, "Computing MD5's for each sequence in '%s'.\n", fileb);
  md5_s *resultB = computeMD5ForEachSequence(B);

  u32bit  numSeqsA = A->getNumberOfSequences();
  u32bit  numSeqsB = B->getNumberOfSequences();
  u32bit  idxA = 0;
  u32bit  idxB = 0;

  fprintf(stderr, "Sorting MD5's.\n");
  qsort(resultA, numSeqsA, sizeof(md5_s), md5_compare);
  qsort(resultB, numSeqsB, sizeof(md5_s), md5_compare);

  fprintf(stderr, "Finding duplicates.\n");
  while ((idxA<numSeqsA) && (idxB<numSeqsB)) {
    int res = md5_compare(resultA+idxA, resultB+idxB);

    if (res == 0) {
      A->find(resultA[idxA].i);
      FastASequenceInCore *sa = loadSequence(A);

      B->find(resultB[idxB].i);
      FastASequenceInCore *sb = loadSequence(B);

      mapDuplicates_Print(filea, sa, fileb, sb);

      //  While the B sequence matches the current A sequence, output a match
      //
      u32bit idxBb = idxB+1;
      int resb = md5_compare(resultA+idxA, resultB+idxBb);
      while (resb == 0) {
        B->find(resultB[idxBb].i);
        FastASequenceInCore *sbb = loadSequence(B);

        mapDuplicates_Print(filea, sa, fileb, sbb);

        delete sbb;

        idxBb++;
        resb = md5_compare(resultA+idxA, resultB+idxBb);
      }

      //  And likewise for A
      //
      u32bit idxAa = idxA+1;
      int resa = md5_compare(resultA+idxAa, resultB+idxB);
      while (resa == 0) {
        A->find(resultA[idxAa].i);
        FastASequenceInCore *saa = loadSequence(A);

        mapDuplicates_Print(filea, saa, fileb, sb);

        delete saa;

        idxAa++;
        resa = md5_compare(resultA+idxAa, resultB+idxB);
      }

      delete sa;
      delete sb;

      idxA++;
      idxB++;
    } else {
      if (res < 0)
        idxA++;
      else
        idxB++;
    }
  }

  delete A;
  delete B;
}











void
checksum(char *name) {
  char            stra[33];
  char            strb[33];
  FastAWrapper   *A = new FastAWrapper(name);
  A->openIndex(FASTA_INDEX_ONLY | FASTA_INDEX_MD5);

  md5_s *result   = computeMD5ForEachSequence(A);
  u32bit  numSeqs = A->getNumberOfSequences();

  for (u32bit idx=0; idx < numSeqs; idx++) {
    if (md5_compare(result+idx, A->getMD5(idx)) != 0) {
      fprintf(stderr, "Checksum mismatch for iid "u32bitFMT": computed %s saved %s\n",
              idx,
              md5_toascii(result+idx, stra),
              md5_toascii(A->getMD5(idx), strb));
    }
  }

  delete A;
}









struct partition_s {
  u32bit  length;
  u32bit  index;
  u32bit  used;
};

int
partition_s_compare(const void *A, const void *B) {
  partition_s *a = (partition_s *)A;
  partition_s *b = (partition_s *)B;
  if (a->length < b->length)
    return(1);
  if (a->length > b->length)
    return(-1);
  return(0);
}

partition_s *loadPartition(void) {
  u32bit        n = fasta->getNumberOfSequences();
  partition_s  *p = new partition_s [n];

  for (u32bit i=0; i<n; i++) {
    p[i].length = fasta->sequenceLength(i);
    p[i].index  = i;
    p[i].used   = 0;
  }

  qsort(p, n, sizeof(partition_s), partition_s_compare);

  return(p);
}

void
outputPartition(char *prefix, partition_s *p, u32bit openP, u32bit n) {
  char  filename[1024];

  //  Check that everything has been partitioned
  //
  for (u32bit i=0; i<n; i++)
    if (p[i].used == 0)
      fprintf(stderr, "ERROR: Failed to partition "u32bitFMT"\n", i);

  if (prefix) {

    //  This rewrites the source fasta file into partitioned fasta files
    //
    for (u32bit o=1; o<=openP; o++) {
      sprintf(filename, "%s-"u32bitFMTW(03)".fasta", prefix, o);

      errno = 0;
      FILE *file = fopen(filename, "w");
      if (errno)
        fprintf(stderr, "Couldn't open '%s' for write: %s\n", filename, strerror(errno));

      for (u32bit i=0; i<n; i++)
        if (p[i].used == o) {
          fasta->find(p[i].index);
          FastASequenceInCore *S = fasta->getSequence();
          fprintf(file, "%s\n", S->header());
          fwrite(S->sequence(), sizeof(char), S->sequenceLength(), file);
          fprintf(file, "\n");

          if (S->sequenceLength() != p[i].length) {
            fprintf(stderr, "Huh?  '%s' "u32bitFMT" != "u32bitFMT"\n", S->header(), S->sequenceLength(), p[i].length);
          }

          delete S;
        }

      fclose(file);
    }

  } else {

    //  This dumps the partition information to stdout.
    //
    fprintf(stdout, u32bitFMT"\n", openP);
    for (u32bit o=1; o<=openP; o++) {
      u32bit  sizeP = 0;
      for (u32bit i=0; i<n; i++)
        if (p[i].used == o)
          sizeP += p[i].length;
      fprintf(stdout, u32bitFMT"]("u32bitFMT")", o, sizeP);
      for (u32bit i=0; i<n; i++)
        if (p[i].used == o)
          fprintf(stdout, " "u32bitFMT"("u32bitFMT")", p[i].index, p[i].length);
      fprintf(stdout, "\n");
    }

  }
}


void
partitionBySize(char *prefix, u64bit partitionSize) {
  u32bit        n = fasta->getNumberOfSequences();
  partition_s  *p = loadPartition();

  u32bit  openP = 1;  //  Currently open partition
  u32bit  sizeP = 0;  //  Size of open partition
  u32bit  seqsP = n;  //  Number of sequences to partition

  //  For any sequences larger than partitionSize, create
  //  partitions containing just one sequence
  //
  for (u32bit i=0; i<n; i++) {
    if (p[i].length > partitionSize) {
      p[i].used = openP++;
      seqsP--;
    }
  }

  //  For the remaining, iterate through the list,
  //  greedily placing the longest sequence that fits
  //  into the open partition
  //
  while (seqsP > 0) {
    for (u32bit i=0; i<n; i++) {
      if ((p[i].used == 0) &&
          (p[i].length + sizeP < partitionSize)) {
        p[i].used = openP;
        sizeP += p[i].length;
        seqsP--;
      }
    }

    openP++;
    sizeP = 0;
  }

  outputPartition(prefix, p, openP-1, n);
  delete [] p;
}


void
partitionByBucket(char *prefix, u64bit partitionSize) {
  u32bit        n = fasta->getNumberOfSequences();
  partition_s  *p = loadPartition();

  if (partitionSize > n)
    partitionSize = n;

  //  The size, in bases, of each partition
  //
  u32bit       *s = new u32bit [partitionSize];
  for (u32bit i=0; i<partitionSize; i++)
    s[i] = 0;

  //  For each sequence
  //
  for (u32bit nextS=0; nextS<n; nextS++) {

    //  find the smallest partition
    //
    u32bit openP = 0;
    for (u32bit i=0; i<partitionSize; i++)
      if (s[i] < s[openP])
        openP = i;

    //  add the next largest sequence to the open partition
    //
    s[openP] += p[nextS].length;
    p[nextS].used = openP+1;
  }

  outputPartition(prefix, p, (u32bit)partitionSize, n);
  delete [] p;
}




void
dumpBlocks(void) {
  FastASequenceOnDisk  *S     = 0L;
  u32bit                seqno = 0;

  failIfNoSource();
  failIfNotRandomAccess();

  bool                  V[256];
  for (u32bit i=0; i<256; i++)
    V[i] = false;
  V[(int)'n'] = true;
  V[(int)'N'] = true;

  S = fasta->getSequenceOnDisk();
  while (S) {
    char    seq    = S->get();
    u32bit  len    = S->sequenceLength();
    bool    nnn    = V[(int)seq];
    char    begseq = seq;
    u32bit  begpos = 0;

    u32bit pos = 1;
    for ( ; pos<len; pos++) {
      S->next();
      seq = S->get();
      if (nnn != V[(int)seq]) {
        fprintf(stdout, "%c "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
                begseq, seqno, begpos, pos, pos - begpos);
        nnn = V[(int)seq];
        begpos = pos;
        begseq = seq;
      }
    }

    fprintf(stdout, "%c "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
            begseq, seqno, begpos, pos, pos - begpos);
    fprintf(stdout, ". "u32bitFMT" "u32bitFMT" "u32bitFMT"\n", seqno, pos, u32bitZERO);

    delete S;
    S = fasta->getSequenceOnDisk();

    seqno++;
  }
}





void
processArray(int argc, char **argv) {

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "--buildinfo", 3) == 0) {
      buildinfo_leaff(stderr);
      buildinfo_libbio(stderr);
      buildinfo_libutil(stderr);
      exit(1);
    } else if (strncmp(argv[arg], "--findduplicates", 3) == 0) {
      arg++;
      findDuplicates(argv[arg]);
    } else if (strncmp(argv[arg], "--mapduplicates", 3) == 0) {
      arg += 2;
      mapDuplicates(argv[arg-1], argv[arg]);
    } else if ((strncmp(argv[arg], "--partition", 3) == 0) ||
               (strncmp(argv[arg], "--partitionmap", 3) == 0)) {
      failIfNoSource();
      failIfNotRandomAccess();

      char *prefix = 0L;
      if (strcmp(argv[arg], "--partition") == 0) {
        arg++;
        prefix = argv[arg];
      }

      arg++;
      //  does the next arg end with gbp, mbp, kbp or bp?  If so,
      //  partition by length, else partition into buckets.
      //
      int     al = (int)strlen(argv[arg]);
      u64bit  ps = (u32bit)strtou32bit(argv[arg], 0L);

      char a3 = (al<3) ? '0' : (char)toLower[argv[arg][al-3]];
      char a2 = (al<2) ? '0' : (char)toLower[argv[arg][al-2]];
      char a1 = (al<1) ? '0' : (char)toLower[argv[arg][al-1]];

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
          fprintf(stderr, "Unknown partition option '%s'\n", argv[arg]), exit(1);
        }

        if (ps == 0)
          fprintf(stderr, "Unknown or zero partition size '%s'\n", argv[arg]), exit(1);
        partitionBySize(prefix, ps);
      } else {
        if (ps == 0)
          fprintf(stderr, "Unknown or zero partition size '%s'\n", argv[arg]), exit(1);
        partitionByBucket(prefix, ps);
      }
    } else if (strncmp(argv[arg], "--checksum", 3) == 0) {
      arg++;
      checksum(argv[arg]);
    } else if (strncmp(argv[arg], "--testindex", 3) == 0) {
      fasta = new FastAWrapper(argv[arg+1]);
      if (fasta->isIndexValid())
        exit(0);
      exit(1);
    } else if (strncmp(argv[arg], "--dumpblocks", 3) == 0) {
      dumpBlocks();
    } else if (strncmp(argv[arg], "--errors", 3) == 0) {
      int    L = strtou32bit(argv[++arg], 0L);  //  Desired length
      int    l = 0;                             //  min of desired length, length of sequence
      int    N = strtou32bit(argv[++arg], 0L);  //  number of copies per sequence
      int    C = strtou32bit(argv[++arg], 0L);  //  number of mutations per copy
      double P = atof(argv[++arg]);             //  probability of mutation
      
      FastASequenceInCore *S = fasta->getSequence();
      while (S) {
        char   *seq = S->sequence();
        char   *hdr = S->header();
        int     len = S->sequenceLength();

        l = len;
        if ((L > 0) && (L < len))
          l = L;

        simseq(seq, hdr, len, N, l, C, P);

        delete S;
        S = fasta->getSequence();
      }
    } else {
      switch(argv[arg][1]) {
        case 'V':
          printOnLoad = !printOnLoad;
          break;
        case 'f':
        case 'F':
          sourceFile = argv[++arg];
          fasta = openNewFile(sourceFile, argv[arg-1]);
          break;
        case 'I':
          switch (argv[arg][2]) {
            case 'i':
              seqIDtype = 'i';
              break;
            case 'e':
              seqIDtype = 'e';
              break;
            default:
              fprintf(stderr, "WARNING: unknown id type '%c'\n", argv[arg][2]);
              break;
          }
          break;
        case 'i':
          switch (argv[arg][2]) {
            case 0:
              failIfNoSource();
              failIfNotRandomAccess();
              fasta->printATADescription(stdout, argv[++arg]);
              break;
            case 'i':
              failIfNoSource();
              failIfNotRandomAccess();
              fasta->printTextDescription(stdout);
              break;
            default:
              fprintf(stderr, "WARNING: unknown option '%s'\n", argv[arg]);
              break;
          }
          break;
        case 'd':
          failIfNoSource();
          failIfNotRandomAccess();
          printf(u32bitFMT"\n", fasta->getNumberOfSequences());
          break;
        case 'L':
          failIfNoSource();
          printSequenceBetweenSize(strtou32bit(argv[arg+1], 0L),
                                   strtou32bit(argv[arg+2], 0L));
          arg += 2;
          break;
        case 'N':
          failIfNoSource();
          printSequenceBetweenNComposition(atof(argv[arg+1]),
                                           atof(argv[arg+2]));
          arg += 2;
          break;
        case 'W':
          failIfNoSource();
          printSequenceBetweenSize(u32bitZERO, ~u32bitZERO);
          break;
        case 'G':
          printRandomlyGeneratedSequence(strtou32bit(argv[arg+1], 0L),
                                         strtou32bit(argv[arg+2], 0L),
                                         strtou32bit(argv[arg+3], 0L)+1);
          arg += 3;
          break;
        case 's':
          failIfNoSource();
          failIfNotRandomAccess();
          fasta->optimizeRandomAccess();
          findSequenceAndPrint(argv[++arg]);
          break;
        case 'S':
          failIfNoSource();
          failIfNotRandomAccess();
          printRangeOfSequences(argv[arg+1], argv[arg+2]);
          arg += 2;
          break;
        case 'r':
          failIfNoSource();
          failIfNotRandomAccess();
          fasta->optimizeRandomAccess();
          findAndPrintRandomSequences(strtou32bit(argv[++arg], 0L));
          break;
        case 'q':
          failIfNoSource();
          failIfNotRandomAccess();
          printIDsFromFile(argv[++arg]);
          break;
        case '6':
          //  If there is a next argument, and it doesn't begin with a
          //  '-', assume it's the value for line breaking, otherwise,
          //  use the default 60.
          //
          withLineBreaks = 60;
          if ((argv[arg+1] != 0L) && (argv[arg+1][0] != '-')) {
            arg++;
            withLineBreaks = strtou32bit(argv[arg], 0L);
          }
          break;
        case 'u':
          toUppercase = !toUppercase;

          if (toUppercase)
            for (int z=0; z<256; z++)
              translate[z] = (char)toUpper[z];
          else
            for (int z=0; z<256; z++)
              translate[z] = (char)z;
          break;
        case 'R':
          reverse = !reverse;
          break;
        case 'C':
          complement = !complement;
          break;
        case 'H':
          withDefLine    = !withDefLine;
          specialDefLine = 0L;
          break;
        case 'h':
          withDefLine    = true;
          specialDefLine = argv[++arg];
          break;
        case 'e':
          begPos = strtou32bit(argv[arg+1], 0L);
          endPos = strtou32bit(argv[arg+2], 0L);
          arg += 2;
          break;
        case 'm':
          printMD5 = !printMD5;
          break;
        case 'A':
          processFile(argv[++arg]);
          break;
      }
    }
    arg++;
  }

  if (fasta) {
    delete fasta;
    fasta = 0L;
  }
}



//  Reads commands from a file, calls processArray()
//
void
processFile(char  *filename) {
  int     argc = 0;
  char   *data = 0L;
  char  **argv = 0L;
  size_t  len  = 0;


  //  If we read from a file, preallocate, otherwise, make it grow
  //
  if (strcmp(filename, "-") == 0) {
    size_t  pos = 0;
    size_t  max = 1024 * 1024;

    data = new char [max];

    while (!feof(stdin)) {
      errno = 0;
      len = fread(data+pos, 1, max - pos, stdin);
      if (errno) {
#ifdef TRUE64BIT
        fprintf(stderr, "error: Couldn't read %lu bytes from '%s'\n%s\n",
                max-pos, filename, strerror(errno));
#else
        fprintf(stderr, "error: Couldn't read %d bytes from '%s'\n%s\n",
                max-pos, filename, strerror(errno));
#endif
        exit(1);
      }
      pos += len;
      
      if (pos >= max) {
        max += (size_t)floor(0.2 * max);
        char *tmpd = new char [max];
        memcpy(tmpd, data, pos);
        delete [] data;
        data = tmpd;
      }
    }

    //  save the length of the thing we read in
    len = pos;
  } else {
    len = sizeOfFile(filename);

    data = new char [len+1];

    errno = 0;
    FILE *F = fopen(filename, "r");
    if (errno) {
      fprintf(stderr, "error: Couldn't open '%s'\n%s\n", filename, strerror(errno));
      exit(1);
    }
    fread(data, 1, len, F);
    if (errno) {
#ifdef TRUE64BIT
      fprintf(stderr, "error: Couldn't read %lu bytes from '%s'\n%s\n",
              len, filename, strerror(errno));
#else
      fprintf(stderr, "error: Couldn't read %d bytes from '%s'\n%s\n",
              len, filename, strerror(errno));
#endif
      exit(1);
    }
    fclose(F);
  }


  //  (over)count the number of words
  //
  for (u32bit i=0; i<len; i++) {
    if (isspace(data[i])) {
      argc++;
      data[i] = 0;
    }
  }


  //  Allocate space for word pointers
  //
  argv = new char * [argc+1];


  //  Set word pointers -- first arg is the name of the "program"
  //
  argv[0] = filename;
  argc = 1;

  for (u32bit pos=0; pos<len; pos++) {

    //  Skip leading whitespace
    while (data[pos] == 0 && pos < len)
      pos++;

    //  save the arg if it's a real arg
    if (pos < len)
      argv[argc++] = data+pos;

    //  skip the word
    while (data[pos] != 0 && pos < len)
      pos++;
  }
 
  processArray(argc, argv);

  delete [] argv;
  delete [] data;
}




int
main(int argc, char **argv) {

  if (argc < 2) {
    fprintf(stderr, usage, argv[0]);
    exit(1);
  }

  for (int z=0; z<256; z++)
    translate[z] = (char)z;

  mtctx = mtInit(getpid() * time(NULL));

  processArray(argc, argv);

  if (fasta)
    delete fasta;
}
