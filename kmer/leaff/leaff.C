#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

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
"       -6:          insert a newline every 60 bases\n"
"       -u:          uppercase all bases\n"
"       -R:          print the sequence reversed\n"
"       -C:          print the sequence complemented\n"
"       -H:          DON'T print the defline\n"
"       -h:          Use the next word as the defline\n"
"                      (use \"-H -H\" to resume using the original defline)\n"
"       -e beg end:  Print only the bases from position 'beg' to position 'end'\n"
"                    (space based!)  If beg == end, then the entire sequence is\n"
"                    printed.  It is an error to specify beg > end, or beg > len,\n"
"                    or end > len.\n"
"       -md5:        Don't print the sequence, but print the md5 checksum\n"
"                    (of the entire sequence) followed by the entire defline.\n"
"\n"
"           XXXXXXXXXXX: how does -e handle reverse complement??\n"
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
              bool                 withLineBreaks=false,
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
    char       b[62];
    int        i = 0;

    while (*t) {
      for (i=0; (*t) && (i < 60); )
        b[i++] = *(t++);
      b[i++] = '\n';
      b[i]   = 0;
      fprintf(stdout, "%s", b);
    }
  } else {
    fprintf(stdout, "%s", n);
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
bool                   reverse           = false;
bool                   complement        = false;
bool                   withDefLine       = true;
char                  *specialDefLine    = 0L;
bool                   withLineBreaks    = false;
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
    s = f->getSequence();
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
    fprintf(stdout, "\n");
  }

  if (mySeq)
    delete s;
}



void
printSequenceBetweenSize(u32bit small, u32bit large) {

  //  XXX: could be faster (maybe) with an index

  while (!f->eof()) {
    FastASequenceInCore *s = f->getSequence();

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
      lastSeq = f->getSequence();
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
    FastASequenceInCore *s1 = F->getSequence();
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
      s1 = A->getSequence();

      A->find(result[idx].i);
      s2 = A->getSequence();

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
      FastASequenceInCore *sa = A->getSequence();

      B->find(resultB[idxB].i);
      FastASequenceInCore *sb = B->getSequence();

      mapDuplicates_Print(filea, sa, fileb, sb);

      //  While the B sequence matches the current A sequence, output a match
      //
      u32bit idxBb = idxB+1;
      int resb = md5_compare(resultA+idxA, resultB+idxBb);
      while (resb == 0) {
        B->find(resultB[idxBb].i);
        FastASequenceInCore *sbb = B->getSequence();

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
        FastASequenceInCore *saa = A->getSequence();

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
    } else {
    switch(argv[arg][1]) {
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
        withLineBreaks = !withLineBreaks;
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
