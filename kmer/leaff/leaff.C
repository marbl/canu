#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>

//  Linux needs to include time.h; others can use sys/time.h.
//  Tru64 is ok with time.h
//
#include <time.h>

#include "libbri.H"
#include "fasta.H"

#include "buildinfo-leaff.h"
#include "buildinfo-libfasta.h"
#include "buildinfo-libbri.h"

#ifdef TRUE64BIT
const char *idxMsg   = "%8u] %5u %10u\n";
const char *descMsg  = "%u\n";
#else
const char *idxMsg   = "%8lu] %5lu %10lu\n";
const char *descMsg  = "%lu\n";
#endif


const char *usage =
"usage: %s [--buildinfo] [-f|-F|-I <fasta-file>] [options]\n"
"\n"
"ENTERTAINMENT OPTIONS\n"
"       -V:         whenever a sequence is read, print the defline to stderr\n"
"\n"
"SOURCE FILE\n"
"       -f file:    use 'file' as an UN-INDEXED source file\n"
"       -Ft file:   use 'file' as an INDEXED source file (the index is built if\n"
"                   it doesn't exist), where 't' is the type of index to build:\n"
"                     i:  internal id's only (the default if t is not specified)\n"
"                     n:  names only (the first word on defline)\n"
"                     d:  full deflines\n"
"       -Ii:        internal 'seqid' (an integer), this is the default for -Fi\n"
"       -Ie:        external 'seqid' (a word), this is the default for -Fn or\n"
"                   -Fd, and is an error for -Fi\n"
"\n"
"ACTIONS (no index needed)\n"
"       -L s l:     print all sequences such that:  s <= length < l\n"
"       -G n s l:   print n random sequences of length 0 < s <= length <= l)\n"
"       -W:         print all sequences (do the whole file)\n"
"\n"
"ACTIONS (index needed)\n"
"       -d:             print the number of sequences in the fasta\n"
"       -i name:        print an ATA compliant index, labelling the source 'name'\n"
"       -ii:            print the index in an almost-human readable format\n"
"       -s seqid:       print the single sequence 'seqid'\n"
"       -S f l:         print all the sequences from IID 'f' to 'l' (inclusive)\n"
"       -r num:         print 'num' randomly picked sequences\n"
"       -q file:        print sequences from the seqid list in 'file'\n"
"\n"
"SEQUENCE OPTIONS\n"
"       -6:          insert a newline every 60 letters\n"
"                      (if the next arg is a number, newlines are inserted every\n"
"                       n letters -- e.g., -6 80 -- disable line breaks with -6 0)\n"
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
"                    NOTE: The coordinates are always relative to the unreversed\n"
"                          sequence, even if -R is specified\n"
"       -md5:        Don't print the sequence, but print the md5 checksum\n"
"                    (of the entire sequence) followed by the entire defline.\n"
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
"       --partition n\n"
"                    Partition the sequences into n roughly equal size pieces\n"
"\n"
"       --partition n[gmk]bp\n"
"                    Partition the sequences into roughly equal size pieces of\n"
"                    size nbp, nkbp, nmbp or ngbp.  Sequences larger that the\n"
"                    partition size are in a partition by themself.\n"
"                    Example: --partition 130mbp\n"
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
"       -F file -R -C -s 3 -s 4    Print thethe fourth and fifth sequences\n"
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

  if ((beg > l) || (end > l) || (beg > end)) {
#ifdef TRUE64BIT
    fprintf(stderr, "WARNING:  Printing %u to %u of sequence %u (len=%u) is out of bounds -- NO SEQUENCE PRINTED!\n",
            beg, end, b->getIID(), l);
#else
    fprintf(stderr, "WARNING:  Printing %lu to %lu of sequence %lu (len=%lu) is out of bounds -- NO SEQUENCE PRINTED!\n",
            beg, end, b->getIID(), l);
#endif
  }

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
        *(m++) = translate[*(s++)];
      break;
    case 1:
      //  reverse
      m  = n + limit - 1;
      s += beg;

      while (limit--)
        *(m--) = translate[*(s++)];
      break;
    case 2:
      //  complement
      m  = n;
      s += beg;

      while (limit--)
        *(m++) = complementSymbol[translate[*(s++)]];
      break;
    case 3:
      //  reverse complement
      m  = n + limit - 1;
      s += beg;

      while (limit--)
        *(m--) = complementSymbol[translate[*(s++)]];
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
FastAWrapper          *f                 = 0L;
char                   seqIDtype         = 'i';
u32bit                 begPos            = ~(u32bit)0;
u32bit                 endPos            = ~(u32bit)0;
bool                   printMD5          = false;
FastASequenceInCore   *lastSeq           = 0L;


void
failIfNoSource(void) {
  if (f == 0L) {
    fprintf(stderr, "No source file specified.\n");
    exit(1);
  }
}

void
failIfNotRandomAccess(void) {
  if (f->isRandomAccess() == false) {
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

  if (f)
    delete f;

  f = new FastAWrapper(name);

  if (arg[1] == 'F') {
    switch (arg[2]) {
      case 'n':
        seqIDtype = 'e';
        f->openIndex(FASTA_INDEX_PLUS_IDS);
        break;
      case 'd':
        seqIDtype = 'e';
        f->openIndex(FASTA_INDEX_PLUS_DEFLINES);
        break;
      default:
        seqIDtype = 'i';
        f->openIndex(FASTA_INDEX_ONLY);
        break;
    }
  }

  delete lastSeq;
  lastSeq = 0L;

  return(f);
}


void
printIID(u32bit iid, FastASequenceInCore *s=0L) {
  bool  mySeq = false;

  if (s == 0L) {
    mySeq = true;
    f->find(iid);
    s = loadSequence(f);
  }

  if (printMD5) {
    md5_s     md5;
    char      sum[33];

    md5_toascii(md5_string(&md5,
                           s->sequence(), s->sequenceLength()),
                sum);

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

  //  XXX: could be faster (maybe) with an index

  while (!f->eof()) {
    FastASequenceInCore *s = loadSequence(f);

    if ((small <= s->sequenceLength()) &&
        (s->sequenceLength() < large))
      printIID(0, s);

    delete s;
  }
}

void
printRandomlyGeneratedSequence(u32bit n, u32bit s, u32bit l) {
  char      bases[4] = {'A', 'C', 'G', 'T'};
  char     *seq      = new char [l + 1];

  if (s > l) {
    u32bit t = s;
    s = l;
    l = t;
  }

  if (s == 0)
    s = 1;

  for (u32bit i=0; i<n; i++) {
    u32bit j = s + (random() % (l-s));
    seq[j] = 0;

    while (j)
      seq[--j] = bases[random() & 0x3];            

    if (withDefLine)
      if (specialDefLine)
        fprintf(stdout, ">%s\n", specialDefLine);
      else
        fprintf(stdout, ">%lu\n", i);

    fprintf(stdout, "%s\n", seq);
  }

  delete [] seq;
}

void
findSequenceAndPrint(char *id) {

  bool found = false;

  if (seqIDtype == 'i')
    found = f->find((u32bit)atoi(id));
  else
    found = f->find(id);

  if (found) {
    if ((lastSeq == 0L) ||
        (lastSeq->getIID() != f->currentIID())) {
      delete lastSeq;
      lastSeq = loadSequence(f);
    }

    printIID(f->currentIID(), lastSeq);
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
    lowID  = atoi(argl);
    if (lowID >= f->getNumberOfSequences()) {
      fprintf(stderr, "ERROR: Internal id of %lu for starting sequence is too large; only %lu sequences.\n",
              lowID, f->getNumberOfSequences());
      fail = true;
    }

    highID = atoi(argh);
    if (highID >= f->getNumberOfSequences()) {
      fprintf(stderr, "ERROR: Internal id of %lu for ending sequence is too large; only %lu sequences.\n",
              highID, f->getNumberOfSequences());
      fail = true;
    }
  } else {
    if (f->find(argl)) {
      lowID = f->currentIID();
    } else {
      fprintf(stderr, "ERROR: Can't find external id '%s' of starting sequence.\n", argl);
      fail = true;
    }
    
    if (f->find(argh)) {
      highID = f->currentIID();
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

  if (num >= f->getNumberOfSequences()) {
    fprintf(stderr, "WARNING: file has %lu sequences, and you asked for %lu.\n",
            f->getNumberOfSequences(), num);
    printSequenceBetweenSize(u32bitZERO, ~u32bitZERO);
    return;
  }

  u32bit  *seqs = new u32bit [f->getNumberOfSequences()];

  for (u32bit i=0; i<f->getNumberOfSequences(); i++)
    seqs[i] = i;

  for (u32bit i=0; i<f->getNumberOfSequences(); i++) {
    u32bit j = (unsigned int)(random() % f->getNumberOfSequences());
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
  u32bit   numSeqs = F->getNumberOfSequences();
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
#ifdef TRUE64BIT
        fprintf(stderr, "Internal error: found two copies of the same sequence iid (%u)!\n", result[idx].i);
#else
        fprintf(stderr, "Internal error: found two copies of the same sequence iid (%lu)!\n", result[idx].i);
#endif
        exit(1);
      }

      A->find(result[idx-1].i);
      s1 = loadSequence(A);

      A->find(result[idx].i);
      s2 = loadSequence(A);

      if (strcmp(s1->sequence(), s2->sequence()) == 0) {
        fprintf(stdout, "%s\n%s\n\n", s1->header(), s2->header());
      } else {
#ifdef TRUE64BIT
        fprintf(stderr, "COLLISION DETECTED BETWEEN IID %u AND %u!\nPLEASE REPORT THIS TO brian.walenz@celera.com!\n",
                result[idx-1].i, result[idx].i);
#else
        fprintf(stderr, "COLLISION DETECTED BETWEEN IID %lu AND %lu!\nPLEASE REPORT THIS TO brian.walenz@celera.com!\n",
                result[idx-1].i, result[idx].i);
#endif
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
#ifdef TRUE64BIT
    fprintf(stdout, "%u <-> %u\n", sa->getIID(), sb->getIID());
#else
    fprintf(stdout, "%lu <-> %lu\n", sa->getIID(), sb->getIID());
#endif
  } else {
#ifdef TRUE64BIT
    fprintf(stderr, "COLLISION DETECTED BETWEEN %s:%u AND %s:%u!\nPLEASE REPORT THIS TO brian.walenz@celera.com!\n",
            filea, sa->getIID(), fileb, sb->getIID());
#else
    fprintf(stderr, "COLLISION DETECTED BETWEEN %s:%lu AND %s:%lu!\nPLEASE REPORT THIS TO brian.walenz@celera.com!\n",
            filea, sa->getIID(), fileb, sb->getIID());
#endif
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
  u32bit        n = f->getNumberOfSequences();
  partition_s  *p = new partition_s [n];

  for (u32bit i=0; i<n; i++) {
    p[i].length = f->sequenceLength(i);
    p[i].index  = i;
    p[i].used   = 0;
  }

  qsort(p, n, sizeof(partition_s), partition_s_compare);

  return(p);
}

void
outputPartition(partition_s *p, u32bit openP, u32bit n) {

  //  Check that everything has been partitioned
  //
  for (u32bit i=0; i<n; i++)
    if (p[i].used == 0)
      fprintf(stderr, "ERROR: Failed to partition %u\n", i);

  fprintf(stdout, "%u\n", openP);
  for (u32bit o=1; o<=openP; o++) {
    u32bit  sizeP = 0;
    for (u32bit i=0; i<n; i++)
      if (p[i].used == o)
        sizeP += p[i].length;
    fprintf(stdout, "%u](%u)", o, sizeP);
    for (u32bit i=0; i<n; i++)
      if (p[i].used == o)
        fprintf(stdout, " %u(%u)", p[i].index, p[i].length);
    fprintf(stdout, "\n");
  }
}


void
partitionBySize(u64bit partitionSize) {
  u32bit        n = f->getNumberOfSequences();
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

  outputPartition(p, openP-1, n);
  delete [] p;
}


void
partitionByBucket(u64bit partitionSize) {
  u32bit        n = f->getNumberOfSequences();
  partition_s  *p = loadPartition();

  u32bit       *s = new u32bit [partitionSize];

  for (u32bit i=0; i<partitionSize; i++)
    s[i] = 0;

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

  outputPartition(p, (u32bit)partitionSize, n);
  delete [] p;
}



void
processArray(int argc, char **argv) {

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "--buildinfo", 3) == 0) {
      buildinfo_leaff(stderr);
      buildinfo_libfasta(stderr);
      buildinfo_libbri(stderr);
      exit(1);
    } else if (strncmp(argv[arg], "--findduplicates", 3) == 0) {
      arg++;
      findDuplicates(argv[arg]);
    } else if (strncmp(argv[arg], "--mapduplicates", 3) == 0) {
      arg += 2;
      mapDuplicates(argv[arg-1], argv[arg]);
    } else if (strncmp(argv[arg], "--partition", 3) == 0) {
      failIfNoSource();
      failIfNotRandomAccess();

      arg++;
      //  does the next arg end with gbp, mbp, kbp or bp?  If so,
      //  partition by length, else partition into buckets.
      //
      int     al = (int)strlen(argv[arg]);
      u64bit  ps = (u32bit)atoi(argv[arg]);

      char a3 = (al<3) ? '0' : (char)tolower(argv[arg][al-3]);
      char a2 = (al<2) ? '0' : (char)tolower(argv[arg][al-2]);
      char a1 = (al<1) ? '0' : (char)tolower(argv[arg][al-1]);

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
          fprintf(stderr, "Unknown partition option '%s'\n", argv[arg]);
          exit(1);
        }

        partitionBySize(ps);
      } else {
        partitionByBucket(ps);
      }
    } else {
    switch(argv[arg][1]) {
      case 'V':
        printOnLoad = !printOnLoad;
        break;
      case 'f':
      case 'F':
        sourceFile = argv[++arg];
        f = openNewFile(sourceFile, argv[arg-1]);
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
            f->printATADescription(stdout, argv[++arg]);
            break;
          case 'i':
            failIfNoSource();
            failIfNotRandomAccess();
            f->printTextDescription(stdout);
            break;
          default:
            fprintf(stderr, "WARNING: unknown option '%s'\n", argv[arg]);
            break;
        }
        break;
      case 'd':
        failIfNoSource();
        failIfNotRandomAccess();
        printf(descMsg, f->getNumberOfSequences());
        break;
      case 'L':
        failIfNoSource();
        printSequenceBetweenSize(atoi(argv[arg+1]),
                                 atoi(argv[arg+2]));
        arg += 2;
        break;
      case 'W':
        failIfNoSource();
        printSequenceBetweenSize(u32bitZERO, ~u32bitZERO);
        break;
      case 'G':
        printRandomlyGeneratedSequence(atoi(argv[arg+1]), atoi(argv[arg+2]), atoi(argv[arg+3])+1);
        arg += 3;
        break;
      case 's':
        failIfNoSource();
        failIfNotRandomAccess();
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
        findAndPrintRandomSequences(atoi(argv[++arg]));
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
          withLineBreaks = atoi(argv[arg]);
        }
        break;
      case 'u':
        toUppercase = !toUppercase;

        if (toUppercase)
          for (int z=0; z<256; z++)
            translate[z] = (char)toupper(z);
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
        begPos = atoi(argv[arg+1]);
        endPos = atoi(argv[arg+2]);
        arg += 2;
        break;
      case 'm':
        printMD5 = !printMD5;
        break;
      case 'A':
        processFile(argv[++arg]);
        break;
    }
    arg++;
  }
  }

  if (f) {
    delete f;
    f = 0L;
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
        fprintf(stderr, "error: Couldn't read %llu bytes from '%s'\n%s\n",
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
      fprintf(stderr, "error: Couldn't read %llu bytes from '%s'\n%s\n",
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

  srandom(getpid() * time(NULL));

  processArray(argc, argv);

  if (f)
    delete f;
}
