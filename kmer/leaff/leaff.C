#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

//  Linux needs to include time.h; others can use sys/time.h.
//  Tru64 is ok with time.h
#include <time.h>
#include "libbri.H"

#include "buildinfo-leaff.h"
#include "buildinfo-libbri.h"

#ifdef TRUE64BIT
const char *idxMsg   = "%8u] %5u %10u\n";
const char *descMsg  = "%u\n";
#else
const char *idxMsg   = "%8lu] %5lu %10lu\n";
const char *descMsg  = "%lu\n";
#endif


//  Tue Jun  4 11:11:11 EDT 2002
//  Made -i use the index if it exists.  Run time from 215 seconds to
//  13 seconds (dbEST).
//
//  Tue Sep  3 13:26:57 EDT 2002
//  Modified -S to limit highID to the number of sequenecs (it used to just exit)
//  and to warn when lowID >= highID



const char *usage =
"usage: %s [--buildinfo] [-f|-F|-I <fasta-file>] [options]\n"
"\n"
"SOURCE FILE\n"
"       -f:  use 'fasta-file' as the source file\n"
"       -F:  use 'fasta-file' as the source file, building the index if needed\n"
"       -I:  only build the index\n"
"       -v:  be verbose\n"
"\n"
"ACTIONS (no index needed)\n"
"       -i:  print the index in an almost-human readable format\n"
"       -d:  print the number of sequences in the fasta\n"
"       -L:  print all sequences such that:  smallest <= length < largest\n"
"       -G:  print n random sequences of length between l and h\n"
"       -W:  print all sequences (do the whole file)\n"
"\n"
"ACTIONS (index needed)\n"
"       -s:  print the single sequence 'seqid'; 0 <= seqid < max\n"
"       -S:  print all the sequences such that:  first <= seqid < last\n"
"       -r:  print n randomly picked sequences; 0 <= n < max\n"
"       -q:  print sequences from the list in 'file'; 0 <= n < max\n"
"\n"
"SEQUENCE OPTIONS\n"
"       -6:  reformat the input fasta with newlines every 60 bases\n"
"       -u:  uppercase all bases\n"
"       -R:  print the sequence reversed\n"
"       -C:  print the sequence complemented\n"
"       -H:  DON'T print the defline\n"
"       -h:  Use the next word as the defline\n"
"            (use \"-H -H\" to resume using the original defline)\n"
"       -e:  Print only the bases from position 'begin' to position 'end'\n"
"          counting starts at zero!\n"
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
printSequence(FastABuffer &b,
              unsigned int beg, unsigned int end,
              bool withLineBreaks=false,
              bool reverse=false, bool complement=false) {
  unsigned int    style = 0;

  if (reverse)     style += 1;
  if (complement)  style += 2;

  u32bit l = b.sequenceLength();

  if (beg > l)  beg = 0;
  if (end > l)  end = l;

  if (beg > end) {
    beg = 0;
    end = l;
  }

  unsigned int    limit = end - beg;
  unsigned char  *n = new unsigned char [end - beg + 1];
  unsigned char  *m;
  unsigned char  *s = b.sequence();

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
    unsigned char *t = n;
    unsigned char  b[62];
    int            i = 0;

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

  delete n;
}


int
comp(const void *a, const void *b) {
  const unsigned int A = *((const unsigned int *)a);
  const unsigned int B = *((const unsigned int *)b);
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
bool           beVerbose         = false;
bool           reverse           = false;
bool           complement        = false;
bool           withDefLine       = true;
char          *specialDefLine    = 0L;
bool           withLineBreaks    = false;
bool           toUppercase       = false;
char          *sourceFile        = 0L;
FastA         *f                 = 0L;
FastABuffer    b;
u32bit         bid;
unsigned int   begPos           = ~(unsigned int)0;
unsigned int   endPos           = ~(unsigned int)0;



void
processArray(int argc, char **argv) {

  int arg = 1;
  while (arg < argc) {
    switch(argv[arg][1]) {
      case '-':  //  Ick!  Must be --buildinfo
        buildinfo_leaff(stderr);
        buildinfo_libbri(stderr);
        exit(1);
        break;
      case 'v':
        beVerbose = true;
        break;
      case 'f':
        if (f)
          delete f;
        sourceFile = argv[++arg];
        f = new FastA(sourceFile, false, beVerbose);
        break;
      case 'F':
      case 'I':
        if (f)
          delete f;
        sourceFile = argv[++arg];
        f = new FastA(sourceFile, true, beVerbose);
        break;
      case 'i':
        if (f == 0L) {
          fprintf(stderr, "No source file specified.\n");
          exit(1);
        } else {
          //  See if the index exists; if so, use a fast method
          //
          if (f->numberOfSequences() > 0) {
            for (u32bit i=0; i < f->numberOfSequences(); i++)
              fprintf(stdout, idxMsg, i, f->headerLength(i), f->sequenceLength(i));
          } else {
            u32bit i;
            for (i=0, f->first(b); f->eof() == false; f->next(b), i++)
              fprintf(stdout, idxMsg, i, b.headerLength(), b.sequenceLength());
          }
        }
        break;
      case 'd':
        if (f == 0L) {
          fprintf(stderr, "No source file specified.\n");
          exit(1);
        } else {
          printf(descMsg, f->numberOfSequences());
        }
        break;
      case 'L':
        if (f == 0L) {
          fprintf(stderr, "No source file specified.\n");
          exit(1);
        } else {
          u32bit small = atoi(argv[++arg]);
          u32bit large = atoi(argv[++arg]);

          for (f->first(b); f->eof() == false; f->next(b))
            if ((small <= b.sequenceLength()) && (b.sequenceLength() < large)) {
              if (specialDefLine)
                fprintf(stdout, ">%s\n", specialDefLine);
              else
                if (withDefLine)
                  fprintf(stdout, ">%s\n", b.header()+1);
              printSequence(b, begPos, endPos, withLineBreaks, reverse, complement);
              fprintf(stdout, "\n");
            }
        }
        break;
      case 'W':
        if (f == 0L) {
          fprintf(stderr, "No source file specified.\n");
          exit(1);
        } else {
          for (f->first(b); f->eof() == false; f->next(b)) {
              if (specialDefLine)
                fprintf(stdout, ">%s\n", specialDefLine);
              else
                if (withDefLine)
                  fprintf(stdout, ">%s\n", b.header()+1);
            printSequence(b, 0, b.sequenceLength(), withLineBreaks, reverse, complement);
            fprintf(stdout, "\n");
          }
        }
        break;
      case 'G':
        //  Give the variables below their own scope.
        {
          unsigned int    n, s, l;
          char            bases[4] = {'A', 'C', 'G', 'T'};

          n = atoi(argv[++arg]);
          s = atoi(argv[++arg]);
          l = atoi(argv[++arg]) + 1;

          char           *seq      = new char [l + 1];

          srandom(getpid() * time(NULL));

          for (unsigned int i=0; i<n; i++) {
            unsigned int j = s + (random() % (l-s));
            seq[j] = 0;
            while (j)
              seq[--j] = bases[random() & 0x3];            
            if (specialDefLine)
              fprintf(stdout, ">%s\n", specialDefLine);
            else
              if (withDefLine)
                fprintf(stdout, ">%s\n", b.header()+1);
            fprintf(stdout, "%s\n", seq);
          }
        }
        break;
      case 's':
        if (f == 0L) {
          fprintf(stderr, "No source file specified.\n");
          exit(1);
        } else {
          u32bit seqID = atoi(argv[++arg]);

          if (seqID < f->numberOfSequences()) {
            f->seek(b, seqID);
            if (specialDefLine)
              fprintf(stdout, ">%s\n", specialDefLine);
            else
              if (withDefLine)
                fprintf(stdout, ">%s\n", b.header()+1);
            printSequence(b, begPos, endPos, withLineBreaks, reverse, complement);
            fprintf(stdout, "\n");
          }
        }
        break;
      case 'S':
        if (f == 0L) {
          fprintf(stderr, "No source file specified.\n");
          exit(1);
        } else {
          u32bit lowID  = atoi(argv[++arg]);
          u32bit highID = atoi(argv[++arg]);

          if (highID > f->numberOfSequences())
            highID = f->numberOfSequences();

          if (lowID >= highID) {
            fprintf(stderr, "Make sure that lowID < highID!  You gave lowID = %d and highID = %d\n", lowID, highID);
          } else {
            if ((lowID < f->numberOfSequences()) && (highID <= f->numberOfSequences())) {
              unsigned int i;
              for (i=lowID, f->seek(b, lowID); (f->eof() == false) && (i<highID); f->next(b), i++) {
                if (specialDefLine)
                  fprintf(stdout, ">%s\n", specialDefLine);
                else
                  if (withDefLine)
                    fprintf(stdout, ">%s\n", b.header()+1);
                printSequence(b, begPos, endPos, withLineBreaks, reverse, complement);
                fprintf(stdout, "\n");
              }
            }
          }
        }
        break;
      case 'r':
        if (f == 0L) {
          fprintf(stderr, "No source file specified.\n");
          exit(1);
        } else {
          unsigned int numberSeqs = atoi(argv[++arg]);

          if (numberSeqs < f->numberOfSequences()) {
            unsigned int  *useThisSequence = new unsigned int [f->numberOfSequences()];
            unsigned int   i, j, t;

	    srandom(getpid() * time(NULL));

            for (i=0; i<f->numberOfSequences(); i++)
              useThisSequence[i] = i;

            for (i=0; i<f->numberOfSequences(); i++) {
              j = (unsigned int)(random() % f->numberOfSequences());
              t = useThisSequence[j];
              useThisSequence[j] = useThisSequence[i];
              useThisSequence[i] = t;
            }

            qsort(useThisSequence, numberSeqs, sizeof(unsigned int), comp);

            for (i=0; i<numberSeqs; i++) {
              f->seek(b, useThisSequence[i]);
              if (specialDefLine)
                fprintf(stdout, ">%s\n", specialDefLine);
              else
                if (withDefLine)
                  fprintf(stdout, ">%s\n", b.header()+1);
              printSequence(b, begPos, endPos, withLineBreaks, reverse, complement);
              fprintf(stdout, "\n");
            }

            delete useThisSequence;
          }
        }
        break;
      case 'q':
        if (f == 0L) {
          fprintf(stderr, "No source file specified!\n");
          exit(1);
        } else {
          FILE *In = fopen(argv[++arg], "r");
          if (In) {
            while (!feof(In)) {
              int i;

              if (fscanf(In, " %d ", &i) == 1) {
                f->seek(b, i);
                if (specialDefLine)
                  fprintf(stdout, ">%s\n", specialDefLine);
                else
                  if (withDefLine)
                    fprintf(stdout, ">%s\n", b.header()+1);
                printSequence(b, begPos, endPos, withLineBreaks, reverse, complement);
                fprintf(stdout, "\n");
              }
            }

            fclose(In);
          } else {
            fprintf(stderr, "Couldn't open '%s'\n", argv[arg]);
          }
        }
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
        withDefLine = !withDefLine;
        if (withDefLine)
          specialDefLine = 0L;
        break;
      case 'h':
        specialDefLine = argv[++arg];
        break;
      case 'e':
        begPos = atoi(argv[arg+1]);
        endPos = atoi(argv[arg+2]);
        arg += 2;
        break;
      case 'A':
        processFile(argv[++arg]);
        break;
    }
    arg++;
  }

  if (f)
    delete f;
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
        fprintf(stderr, "error: Couldn't read %lu bytes from '%s'\n%s\n",
                max-pos, filename, strerror(errno));
        exit(1);
      }
      pos += len;
      
      if (pos >= max) {
        max += 0.2 * max;
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
      fprintf(stderr, "error: Couldn't read %lu bytes from '%s'\n%s\n", len, filename, strerror(errno));
      exit(1);
    }
    fclose(F);
  }


  //  (over)count the number of words
  //
  for (int i=0; i<len; i++) {
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

  for (int pos=0; pos<len; pos++) {

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

  processArray(argc, argv);
}
