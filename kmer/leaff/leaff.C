#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#define WITH_MD5

#ifdef WITH_MD5
#include "../external/md5lib/global.h"
#include "../external/md5lib/md5.h"
#endif

//  Linux needs to include time.h; others can use sys/time.h.
//  Tru64 is ok with time.h
//
#include <time.h>


//  only for complementSymbol
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
"       -Ft file:   use 'file' as an INDEXED source file (the index is built if it doesn't exist)\n"
"                   where 't' is the type of index to build:\n"
"                     i:  internal id's only (the default if t is not specified)\n"
"                     n:  names only (the first word on defline)\n"
"                     d:  full deflines\n"
"       -Ii:        internal 'seqid' (an integer), the default for -Fi\n"
"       -Ie:        external 'seqid' (a word), the default for -Fn or -Fd, error for -Fi\n"
"\n"
"ACTIONS (no index needed)\n"
"       -L s l:     print all sequences such that:  s <= length < l\n"
"       -G n s l:   print n random sequences of length between s and l (0 < s <= l)\n"
"       -W:         print all sequences (do the whole file)\n"
"\n"
"ACTIONS (index needed)\n"
"       -d:             print the number of sequences in the fasta\n"
"       -i name:        print an ATA compliant index, labelling the source 'name'\n"
"       -ii:            print the index in an almost-human readable format\n"
"       -s seqid:       print the single sequence 'seqid'\n"
"       -S first last:  print all the sequences from 'first' to 'last' (inclusive)\n"
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
"       -e beg end:  Print only the bases from position 'beg' to position 'end' (space\n"
"                    based!)  If beg == end, then the entire sequence is printed.  It\n"
"                    is an error to specify beg > end, or beg > len, or end > len.\n"
#ifdef WITH_MD5
"       -md5:        Don't print the sequence, but print the md5 checksum\n"
"                    (of the entire sequence) followed by the entire defline.\n"
#endif
"\n"
"           XXXXXXXXXXX: how does -e handle reverse complement??\n"
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

#if 0
  if (beg > l)  beg = 0;
  if (end > l)  end = l;

  if (beg > end) {
    beg = 0;
    end = l;
  }
#else
  if ((beg > l) || (end > l) || (beg > end)) {
    fprintf(stderr, "WARNING:  Printing %u to %u of sequence %u (len=%u) is out of bounds -- NO SEQUENCE PRINTED!\n",
            beg, end, b->getIID(), l);
  }
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

#ifdef WITH_MD5
  if (printMD5) {
    MD5_CTX         ctx;
    unsigned char   dig[16];
    char           *prt     = "0123456789abcdef";
    char            sum[33];

    MD5Init(&ctx);
    MD5Update(&ctx, (unsigned char *)s->sequence(), s->sequenceLength());
    MD5Final(dig, &ctx);

    for (int i=0; i<16; i++) {
      sum[2*i  ] = prt[(dig[i] >> 4) & 0x0f];
      sum[2*i+1] = prt[(dig[i])      & 0x0f];
    }
    sum[32] = 0;

    fprintf(stdout, "%s %s\n", sum, s->header());
  } else {
#endif
    if (withDefLine)
      if (specialDefLine)
        fprintf(stdout, ">%s\n", specialDefLine);
      else
        fprintf(stdout, ">%s\n", s->header()+1);

    printSequence(s, begPos, endPos, withLineBreaks, reverse, complement);
    fprintf(stdout, "\n");
#ifdef WITH_MD5
  }
#endif

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
              lowID, f->getNumberOfSequences());
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



void
processArray(int argc, char **argv) {

  int arg = 1;
  while (arg < argc) {
    switch(argv[arg][1]) {
      case '-':  //  Ick!  Must be --buildinfo
        buildinfo_leaff(stderr);
        buildinfo_libfasta(stderr);
        buildinfo_libbri(stderr);
        exit(1);
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
#ifdef WITH_MD5
      case 'm':
        printMD5 = !printMD5;
        break;
#endif
      case 'A':
        processFile(argv[++arg]);
        break;
    }
    arg++;
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
        fprintf(stderr, "error: Couldn't read %llu bytes from '%s'\n%s\n",
                max-pos, filename, strerror(errno));
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
      fprintf(stderr, "error: Couldn't read %llu bytes from '%s'\n%s\n", len, filename, strerror(errno));
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
