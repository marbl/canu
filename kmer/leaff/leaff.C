#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>

//  Linux needs to include time.h; others can use sys/time.h.
//  Tru64 is ok with time.h
//
#include <time.h>

#include <algorithm>

#include "bio++.H"
#include "seqCache.H"
#include "seqStore.H"  //  constructSeqStore()


void          simseq(char *,char *,int,int,int,int,double);


const char *usage =
"usage: %s [-f <fasta-file>] [options]\n"
"\n"
"SOURCE FILE\n"
"       -f file:     use sequence in 'file' (-F is also allowed for historical reasons)\n"
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
"       -ends n      Print n bases from each end of the sequence.  One input\n"
"                    sequence generates two output sequences, with '_5' or '_3'\n"
"                    appended to the ID.  If 2n >= length of the sequence, the\n"
"                    sequence itself is printed, no ends are extracted (they\n"
"                    overlap).\n"
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
"       --segment prefix n\n"
"                    Splits the sequences into n files, prefix-###.fasta.\n"
"                    Sequences are not reordered; the first n sequences are in\n"
"                    the first file, the next n in the second file, etc.\n"
"\n"
"       --gccontent a.fasta\n"
"                    Reports the GC content over a sliding window of\n"
"                    3, 5, 11, 51, 101, 201, 501, 1001, 2001 bp.\n"
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
"                             -- make C=10 sequencer runs of N=10 EST sequences\n"
"                                of length 500bp each.\n"
"                          to simulate mRNA from genes, use L=0, N=10, C=10\n"
"                          to simulate reads from genomes, use L=800, N=10, C=1\n"
"                             -- of course, N= should be increased to give the\n"
"                                appropriate depth of coverage\n"
"\n"
"       --stats\n"
"                    Reports size statistics; number, N50, sum, largest.\n"
"\n"
"\n"
"EXPERT OPTIONS\n"
"       -A:  Read actions from 'file'\n"
"\n"
"\n"
"CONVERSION\n"
"       --seqstore out.seqStore\n"
"                    Converts the input file (-f) to a seqStore file.\n"
"\n"
"\n"
"EXAMPLES\n"
"       Options are ORDER DEPENDENT.  Sequences are printed whenever an\n"
"       ACTION occurs on the command line.  SEQUENCE OPTIONS are not reset\n"
"       when a sequence is printed.\n"
"\n"
"\n"
"       Print the first 10 bases of the fourth sequence in file:\n"
"           -f file -e 0 10 -s 3\n"
"\n"
"       Print the first 10 bases of the fourth and fifth sequences in file:\n"
"           -f file -e 0 10 -s 3 -s 4\n"
"\n"
"       Print the fourth and fifth sequences reverse complemented:\n"
"           -f file -R -C -s 3 -s 4\n"
"\n"
"       Convert 'file' to a seqStore 'out.seqStore':\n"
"           -f file --seqstore out.seqStore\n"
"\n";


//  This will do our character translation (e.g., to upper)
//
char translate[256];

void  processFile(char  *filename);
void  processArray(int argc, char **argv);

void
printSequence(char        *d,
              seqInCore   *b,
              u32bit       beg,
              u32bit       end,
              u32bit       withLineBreaks=0,
              bool         reverse=false,
              bool         complement=false) {
  u32bit           style  = 0;
  if (reverse)     style += 1;
  if (complement)  style += 2;

  u32bit l = b->sequenceLength();

  if (beg == end)
    return;
  if (beg > l)
    return;
  if (end > l)
    end = l;

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

  if (d)
    fprintf(stdout, ">%s\n", d);

  if (withLineBreaks) {
    char      *t = n;
    char      *a = new char [withLineBreaks+1];

    while (*t) {
      u32bit i=0;
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


int
u32bit_compare(const void *a, const void *b) {
  const u32bit A = *((const u32bit *)a);
  const u32bit B = *((const u32bit *)b);
  if (A < B)
    return(-1);
  if (A > B)
    return(1);
  return(0);
}


bool                   reverse           = false;
bool                   complement        = false;
bool                   withDefLine       = true;
char                  *specialDefLine    = 0L;
u32bit                 withLineBreaks    = 0;
bool                   toUppercase       = false;
char                  *sourceFile        = 0L;
seqCache              *fasta             = 0L;
seqCache              *cache             = 0L;
u32bit                 begPos            =  (u32bit)0;
u32bit                 endPos            = ~(u32bit)0;
u32bit                 endExtract        = ~(u32bit)0;
bool                   printMD5          = false;
mt_s                  *mtctx             = 0L;

void
failIfNoSource(void) {
  if (fasta == 0L) {
    fprintf(stderr, "No source file specified.\n");
    exit(1);
  }
}



void
printIID(u32bit iid, seqInCore *s=0L) {
  bool  mySeq = false;

  if (s == 0L) {
    mySeq = true;
    s = fasta->getSequenceInCore(iid);
  }

  if (printMD5) {
    md5_s     md5;
    char      sum[33];

    md5_toascii(md5_string(&md5, s->sequence(), s->sequenceLength()), sum);

    fprintf(stdout, "%s %s\n", sum, s->header());
  } else if (endExtract + endExtract < s->sequenceLength()) {
    char d[1024];

    sprintf(d, "%s_5", s->header());
    printSequence(d, s, 0, endExtract, withLineBreaks, reverse, complement);

    sprintf(d, "%s_3", s->header());
    printSequence(d, s, s->sequenceLength()-endExtract, s->sequenceLength(), withLineBreaks, reverse, complement);
  } else {
    char *d = 0L;

    if (withDefLine)
      d = (specialDefLine) ? specialDefLine : s->header();

    printSequence(d, s, begPos, endPos, withLineBreaks, reverse, complement);
  }

  if (mySeq)
    delete s;
}


void
printDescription(char *label) {
  for (u32bit s=0; s<fasta->getNumberOfSequences(); s++) {
    fprintf(stdout, "G\tseq\t%s:"u32bitFMT"\t"u32bitFMT"\t%s\n",
            label, s, fasta->getSequenceLength(s), ">unimplemented");
  }
}



void
printSequenceBetweenSize(u32bit small, u32bit large) {

  for (u32bit s=0; s<fasta->getNumberOfSequences(); s++) {
    if ((small <= fasta->getSequenceLength(s)) && (fasta->getSequenceLength(s) < large))
      printIID(s, 0L);
  }
}

void
printSequenceBetweenNComposition(double small, double large) {

  for (u32bit s=0; s<fasta->getNumberOfSequences(); s++) {
    seqInCore *S = fasta->getSequenceInCore(s);

    u32bit   Ns  = 0;
    u32bit   len = S->sequenceLength();
    char    *seq = S->sequence();
    for (u32bit i=0; i<len; i++)
      if ((seq[i] == 'n') || (seq[i] == 'N'))
        Ns++;

    double Np = 100.0 * Ns / len;

    if ((small <= Np) && (Np < large))
      printIID(0, S);

    delete S;
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
    u32bit p = 0;

    seq[j] = '\n';

    while (p < j)
      seq[p++] = bases[mtRandom32(mtctx) & 0x3];            

    if (withDefLine)
      if (specialDefLine)
        fprintf(stdout, ">%s\n", specialDefLine);
      else
        fprintf(stdout, ">"u32bitFMT"\n", i);

    fwrite(seq, sizeof(char), j+1, stdout);
  }

  delete [] seq;
}

void
findSequenceAndPrint(char *id) {
  seqInCore  *S = fasta->getSequenceInCore(id);

  if (S == 0L)
    fprintf(stderr, "WARNING: Didn't find sequence with name/iid '%s'\n", id);
  else
    printIID(0, S);
}

void
printRangeOfSequences(char *argl, char *argh) {
  u32bit lowID  = fasta->getSequenceIID(argl);
  u32bit highID = fasta->getSequenceIID(argh);

  if (lowID > highID) {
    u32bit t = lowID;
    lowID    = highID;
    highID   = t;
  }

  for (u32bit seqID=lowID; seqID <= highID; seqID++)
    printIID(seqID, 0L);
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

  qsort(seqs, num, sizeof(u32bit), u32bit_compare);

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
  char        x = B.read();

  //  For optimal performance, we should sort the list of ID's given
  //  by their IID, but the user might have a good reason for wanting
  //  them unsorted.

  while (B.eof() == false) {
    while (whitespaceSymbol[x] && (B.eof() == false))
      x = B.read();

    if (B.eof() == false) {
      idLen = 0;

      while (!whitespaceSymbol[x] && (B.eof() == false)) {
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

      findSequenceAndPrint(id);
    }
  }

  delete [] id;
}


md5_s *
computeMD5ForEachSequence(seqCache *F) {
  u32bit   numSeqs = F->getNumberOfSequences();
  md5_s   *result  = new md5_s [numSeqs];

  for (u32bit idx=0; idx < numSeqs; idx++) {
    seqInCore *s1 = F->getSequenceInCore(idx);
    md5_string(result+idx, s1->sequence(), s1->sequenceLength());
    result[idx].i = s1->getIID();
    delete s1;
  }

  return(result);
}

void
findDuplicates(char *filename) {
  seqInCore  *s1 = 0L;
  seqInCore  *s2 = 0L;
  seqCache   *A = new seqCache(filename);


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

      s1 = A->getSequenceInCore(result[idx-1].i);
      s2 = A->getSequenceInCore(result[idx].i);

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
mapDuplicates_Print(char *filea, seqInCore *sa,
                    char *fileb, seqInCore *sb) {
  if (strcmp(sa->sequence(), sb->sequence()) == 0) {
    fprintf(stdout, u32bitFMT" <-> "u32bitFMT"\n", sa->getIID(), sb->getIID());
  } else {
    fprintf(stderr, "COLLISION DETECTED BETWEEN %s:"u32bitFMT" AND %s:"u32bitFMT"!\nPLEASE REPORT THIS TO bri@walenz.org!\n",
            filea, sa->getIID(), fileb, sb->getIID());
  }
}



void
mapDuplicates(char *filea, char *fileb) {
  fprintf(stderr, "Computing MD5's for each sequence in '%s'.\n", filea);
  seqCache  *A = new seqCache(filea);
  md5_s     *resultA = computeMD5ForEachSequence(A);

  fprintf(stderr, "Computing MD5's for each sequence in '%s'.\n", fileb);
  seqCache  *B = new seqCache(fileb);
  md5_s     *resultB = computeMD5ForEachSequence(B);

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
      seqInCore *sa = A->getSequenceInCore(resultA[idxA].i);
      seqInCore *sb = B->getSequenceInCore(resultB[idxB].i);

      mapDuplicates_Print(filea, sa, fileb, sb);

      //  While the B sequence matches the current A sequence, output a match
      //
      u32bit idxBb = idxB+1;
      int resb = md5_compare(resultA+idxA, resultB+idxBb);
      while (resb == 0) {
        seqInCore *sbb = B->getSequenceInCore(resultB[idxBb].i);

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
        seqInCore *saa = A->getSequenceInCore(resultA[idxAa].i);

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
computeGCcontent(char *name) {
  seqCache   *A = new seqCache(name);

  for (u32bit idx=0; idx < A->getNumberOfSequences(); idx++) {
    seqInCore *S = A->getSequenceInCore(idx);
    char      *s = S->sequence();
    u32bit     genomeLength = S->sequenceLength();

    fprintf(stdout, ">%s\n", S->header());

    int gc[256] = {0};
    gc['c'] = 1;
    gc['C'] = 1;
    gc['g'] = 1;
    gc['G'] = 1;

    //  Replace the sequence with "g or c".  We can't do this inline,
    //  since output reports the sequence too.  The extra 1000 at the
    //  end is important, since we do not bother checking for the end
    //  of the valid data, just assume that it's zero.
    //
    char                *g = new char [S->sequenceLength() + 1000];
    for (u32bit i=0; i<genomeLength+1000; i++)
      g[i] = 0;
    for (u32bit i=0; i<genomeLength; i++)
      g[i] = gc[s[i]];

    //  This stolen from depthOfPolishes.C

    u32bit  ave3    = 0;
    u32bit  ave5    = 0;
    u32bit  ave11   = 0;
    u32bit  ave51   = 0;
    u32bit  ave101  = 0;
    u32bit  ave201  = 0;
    u32bit  ave501  = 0;
    u32bit  ave1001 = 0;
    u32bit  ave2001 = 0;

    //  Preload the averages
    ave3   += g[0];
    ave5   += g[0] + g[1];

    for (u32bit i=0; i<5; i++)     ave11   += g[i];
    for (u32bit i=0; i<25; i++)    ave51   += g[i];
    for (u32bit i=0; i<50; i++)    ave101  += g[i];
    for (u32bit i=0; i<100; i++)   ave201  += g[i];
    for (u32bit i=0; i<250; i++)   ave501  += g[i];
    for (u32bit i=0; i<500; i++)   ave1001 += g[i];
    for (u32bit i=0; i<1000; i++)  ave2001 += g[i];

    for (u32bit i=0; i<genomeLength; i++) {
      ave3    += g[i+1]    - ((i >    1) ? g[i-2]    : 0);
      ave5    += g[i+2]    - ((i >    2) ? g[i-3]    : 0);
      ave11   += g[i+5]    - ((i >    5) ? g[i-6]    : 0);
      ave51   += g[i+25]   - ((i >   25) ? g[i-25]   : 0);
      ave101  += g[i+50]   - ((i >   50) ? g[i-51]   : 0);
      ave201  += g[i+100]  - ((i >  100) ? g[i-101]  : 0);
      ave501  += g[i+250]  - ((i >  250) ? g[i-251]  : 0);
      ave1001 += g[i+500]  - ((i >  500) ? g[i-501]  : 0);
      ave2001 += g[i+1000] - ((i > 1000) ? g[i-1001] : 0);

      fprintf(stdout, u32bitFMT"\t"u32bitFMT"\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
              i,
              s[i],
              ave3    / (double)((i >=   1)  ? 3    - ((i < genomeLength -   1) ? 0 : i +    2 - genomeLength) : i+2),
              ave5    / (double)((i >=   2)  ? 5    - ((i < genomeLength -   2) ? 0 : i +    3 - genomeLength) : i+3),
              ave11   / (double)((i >=   5)  ? 11   - ((i < genomeLength -   4) ? 0 : i +    5 - genomeLength) : i+6),
              ave51   / (double)((i >=  25)  ? 51   - ((i < genomeLength -  24) ? 0 : i +   25 - genomeLength) : i+26),
              ave101  / (double)((i >=  50)  ? 101  - ((i < genomeLength -  49) ? 0 : i +   50 - genomeLength) : i+51),
              ave201  / (double)((i >= 100)  ? 201  - ((i < genomeLength -  99) ? 0 : i +  100 - genomeLength) : i+101),
              ave501  / (double)((i >= 250)  ? 501  - ((i < genomeLength - 249) ? 0 : i +  250 - genomeLength) : i+251),
              ave1001 / (double)((i >= 500)  ? 1001 - ((i < genomeLength - 499) ? 0 : i +  500 - genomeLength) : i+501),
              ave2001 / (double)((i >= 1000) ? 2001 - ((i < genomeLength - 999) ? 0 : i + 1000 - genomeLength) : i+1001));
    }

    delete [] g;
    delete    S;
  }
}
  





struct partition_s {
  u32bit  length;
  u32bit  index;
  u32bit  partition;
};

int
partition_s_compare(const void *A, const void *B) {
  const partition_s *a = (const partition_s *)A;
  const partition_s *b = (const partition_s *)B;
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
    p[i].length    = fasta->getSequenceLength(i);
    p[i].index     = i;
    p[i].partition = 0;
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
    if (p[i].partition == 0)
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
        if (p[i].partition == o) {
          seqInCore *S = fasta->getSequenceInCore(p[i].index);
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
        if (p[i].partition == o)
          sizeP += p[i].length;
      fprintf(stdout, u32bitFMT"]("u32bitFMT")", o, sizeP);
      for (u32bit i=0; i<n; i++)
        if (p[i].partition == o)
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
      p[i].partition = openP++;
      seqsP--;
    }
  }

  //  For the remaining, iterate through the list,
  //  greedily placing the longest sequence that fits
  //  into the open partition
  //
  while (seqsP > 0) {
    for (u32bit i=0; i<n; i++) {
      if ((p[i].partition == 0) &&
          (p[i].length + sizeP < partitionSize)) {
        p[i].partition = openP;
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
    p[nextS].partition = openP+1;
  }

  outputPartition(prefix, p, (u32bit)partitionSize, n);
  delete [] p;
}


void
partitionBySegment(char *prefix, u64bit numSegments) {
  u32bit        n = fasta->getNumberOfSequences();
  partition_s  *p = new partition_s [n];
  u32bit        numSeqPerPart = (u32bit)ceil(n / (double)numSegments);

  for (u32bit i=0; i<n; i++) {
    p[i].length    = fasta->getSequenceLength(i);
    p[i].index     = i;
    p[i].partition = i / numSeqPerPart + 1;
  }

  outputPartition(prefix, p, numSegments, n);
  delete [] p;
}



void
dumpBlocks(void) {
  seqInCore   *S     = 0L;
  u32bit       seqno = 0;

  failIfNoSource();

  bool                  V[256];
  for (u32bit i=0; i<256; i++)
    V[i] = false;
  V['n'] = true;
  V['N'] = true;

#warning old style seq iteration
  S = fasta->getSequenceInCore(seqno);
  while (S) {
    u32bit  len    = S->sequenceLength();
    char    begseq = S->sequence()[0];
    bool    nnn    = V[begseq];
    u32bit  begpos = 0;
    u32bit  pos    = 0;

    for (pos=0; pos<len; pos++) {
      char seq = S->sequence()[pos];

      if (nnn != V[seq]) {
        fprintf(stdout, "%c "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
                begseq, seqno, begpos, pos, pos - begpos);
        nnn = V[seq];
        begpos = pos;
        begseq = seq;
      }
    }

    fprintf(stdout, "%c "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
            begseq, seqno, begpos, pos, pos - begpos);
    fprintf(stdout, ". "u32bitFMT" "u32bitFMT" "u32bitFMT"\n", seqno, pos, u32bitZERO);

    delete S;
    S = fasta->getSequenceInCore(++seqno);
  }
}



void
stats(void) {
  seqInCore   *S     = 0L;

  failIfNoSource();

  bool                  V[256];
  for (u32bit i=0; i<256; i++)
    V[i] = false;
  V['n'] = true;
  V['N'] = true;

  u32bit  numSeq = fasta->getNumberOfSequences();

  //  sum of spans
  //  lengths of spans (for N50 computation)
  //  
  u64bit    Ss = 0;
  u32bit   *Ls = new u32bit [numSeq];

  //  sum of bases
  //  lengths of bases (for N50 computation)
  //  
  u64bit    Sb = 0;
  u32bit   *Lb = new u32bit [numSeq];

  for (u32bit i=0; i<numSeq; i++)
    Ls[i] = Lb[i] = 0;

  u32bit seqno = 0;

#warning old style seq iteration
  S = fasta->getSequenceInCore(seqno);
  while (S) {
    u32bit  len    = S->sequenceLength();
    u32bit  span   = len;
    u32bit  base   = len;

    for (u32bit pos=1; pos<len; pos++) {
      if (V[S->sequence()[pos]])
        base--;
    }

    Ss += span;
    Sb += base;

    Ls[S->getIID()] = span;
    Lb[S->getIID()] = base;

    delete S;
    S = fasta->getSequenceInCore(++seqno);
  }

  qsort(Ls, numSeq, sizeof(u32bit), u32bit_compare);
  qsort(Lb, numSeq, sizeof(u32bit), u32bit_compare);

  u32bit  n50s=0, n50b=0;

  for (u32bit i=0, sum=0; sum < Ss/2; i++) {
    n50s = Ls[i];
    sum += Ls[i];
  }

  for (u32bit i=0, sum=0; sum < Sb/2; i++) {
    n50b = Lb[i];
    sum += Lb[i];
  }

  fprintf(stderr, "%s\n", fasta->getSourceName());
  fprintf(stderr, "\n");
  fprintf(stderr, "numSeqs  "u32bitFMT"\n", fasta->getNumberOfSequences());
  fprintf(stderr, "\n");
  fprintf(stderr, "SPAN\n");
  fprintf(stderr, "n50      "u32bitFMTW(10)"\n", n50s);
  fprintf(stderr, "smallest "u32bitFMTW(10)"\n", Ls[0]);
  fprintf(stderr, "largest  "u32bitFMTW(10)"\n", Ls[numSeq-1]);
  fprintf(stderr, "totLen   "u64bitFMTW(10)"\n", Ss);
  fprintf(stderr, "\n");
  fprintf(stderr, "BASES\n");
  fprintf(stderr, "n50      "u32bitFMTW(10)"\n", n50b);
  fprintf(stderr, "smallest "u32bitFMTW(10)"\n", Lb[0]);
  fprintf(stderr, "largest  "u32bitFMTW(10)"\n", Lb[numSeq-1]);
  fprintf(stderr, "totLen   "u64bitFMTW(10)"\n", Sb);

  delete [] Ls;
  delete [] Lb;
}




void
processArray(int argc, char **argv) {

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "--findduplicates", 3) == 0) {
      arg++;
      findDuplicates(argv[arg]);
    } else if (strncmp(argv[arg], "--mapduplicates", 3) == 0) {
      arg += 2;
      mapDuplicates(argv[arg-1], argv[arg]);
    } else if ((strncmp(argv[arg], "--partition", 3) == 0) ||
               (strncmp(argv[arg], "--partitionmap", 3) == 0)) {
      failIfNoSource();

      char *prefix = 0L;
      if (strcmp(argv[arg], "--partition") == 0) {
        arg++;
        prefix = argv[arg];
      }

      arg++;
      //  does the next arg end with gbp, mbp, kbp or bp?  If so,
      //  partition by length, else partition into buckets.
      //
      int     al = strlen(argv[arg]);
      u64bit  ps = strtou64bit(argv[arg], 0L);

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
    } else if (strncmp(argv[arg], "--segment", 5) == 0) {
      failIfNoSource();
      partitionBySegment(argv[arg+1], strtou32bit(argv[arg+2], 0L));
      arg += 2;
    } else if (strncmp(argv[arg], "--gccontent", 4) == 0) {
      arg++;
      computeGCcontent(argv[arg]);
    } else if (strncmp(argv[arg], "--dumpblocks", 3) == 0) {
      dumpBlocks();
    } else if (strncmp(argv[arg], "--stats", 4) == 0) {
      stats();
    } else if (strncmp(argv[arg], "--errors", 3) == 0) {
      int    L = strtou32bit(argv[++arg], 0L);  //  Desired length
      int    l = 0;                             //  min of desired length, length of sequence
      int    N = strtou32bit(argv[++arg], 0L);  //  number of copies per sequence
      int    C = strtou32bit(argv[++arg], 0L);  //  number of mutations per copy
      double P = atof(argv[++arg]);             //  probability of mutation
      u32bit i = 0;

      seqInCore *S = fasta->getSequenceInCore(i++);
      while (S) {
        char   *seq = S->sequence();
        char   *hdr = S->header();
        int     len = S->sequenceLength();

        l = len;
        if ((L > 0) && (L < len))
          l = L;

        simseq(seq, hdr, len, N, l, C, P);

        delete S;
        S = fasta->getSequenceInCore(i++);
      }
    } else if (strncmp(argv[arg], "--seqstore", 5) == 0) {
      constructSeqStore(argv[++arg], fasta);
      exit(0);
    } else {
      switch(argv[arg][1]) {
        case 'f':
        case 'F':
          sourceFile = argv[++arg];
          delete fasta;
          fasta = new seqCache(sourceFile);
          break;
        case 'i':
          failIfNoSource();

          //  Check that the next arg is something reasonable -- not NULL and not another option.
          //
          ++arg;
          if ((argv[arg] == 0L) || (argv[arg][0] == '-'))
            fprintf(stderr, "ERROR: next arg to -i should be ATA label, I got '%s'\n",
                    (argv[arg] == 0L) ? "(nullpointer)" : argv[arg]), exit(1);

          printDescription(argv[arg]);
          break;
        case 'd':
          failIfNoSource();
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
          findSequenceAndPrint(argv[++arg]);
          break;
        case 'S':
          failIfNoSource();
          printRangeOfSequences(argv[arg+1], argv[arg+2]);
          arg += 2;
          break;
        case 'r':
          failIfNoSource();
          findAndPrintRandomSequences(strtou32bit(argv[++arg], 0L));
          break;
        case 'q':
          failIfNoSource();
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
          switch (argv[arg][2]) {
            case 0:
              begPos = strtou32bit(argv[arg+1], 0L);
              endPos = strtou32bit(argv[arg+2], 0L);
              arg += 2;
              break;
            case 'n':
              endExtract = strtou32bit(argv[++arg], 0L);
              break;
          }
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
