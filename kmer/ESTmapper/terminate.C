#include "util++.H"
#include "bio++.H"
#include "sim4.H"

#include "seqCache.H"

//  Terminates an ESTmapper run.
//
//  Splits a fasta file into multiple fasta files based on the first
//  occurrence of the iid.  So, if the iid is in polishes and
//  list-of-iid, the sequence is written to fasta1.  If the iid isn't
//  in the input (polishes or list-of-iid), put it into fasta3.
//  Any number of -p and -i can be specified.
//
//  -P polishes    fasta1
//  -I list-of-iid fasta2
//  -O fasta3
//  -i input.fasta
//
//  -P polishes MUST be sorted by cDNA iid.  Relatively easy to fix this,
//  just read all the polishes when building an iidReaderWriter, storing the
//  iid's we see into an array.

class iidReaderWriter {
public:
  iidReaderWriter(char *infile, char *otfile, bool ispolishes) {
    isPolishes = ispolishes;
    inPolishes = 0L;
    inFile     = 0L;

    if (isPolishes) {
      inPolishes = new sim4polishReader(infile);
    } else {
      errno = 0;
      inFile = fopen(infile, "r");
      if (errno)
        fprintf(stderr, "iidReaderWriter-- can't open '%s': %s\n", infile, strerror(errno)), exit(1);
    }

    errno = 0;
    otFile = fopen(otfile, "w");
    if (errno)
      fprintf(stderr, "iidReaderWriter-- can't open '%s': %s\n", otfile, strerror(errno)), exit(1);

    iids = 0L;
  };

  ~iidReaderWriter() {

    delete [] iids;

    if (isPolishes)
      delete inPolishes;
    else
      fclose(inFile);

    fclose(otFile);
  };

  bool     thisIID(u32bit targetiid) {
    return(iids[targetiid]);
  };

  void     writeSequence(seqInCore *S) {
    fprintf(otFile, ">%s\n%s\n", S->header(), S->sequence());
  };

  void     load(u32bit maxiid) {
    iids = new bool [maxiid];

    for (u32bit i=0; i<maxiid; i++)
      iids[i] = false;

    if (isPolishes) {
      sim4polish *p = inPolishes->nextAlignment();
      while (p) {
        iids[p->_estID] = true;
        delete p;
        p = inPolishes->nextAlignment();
      }
    } else {
      fscanf(inFile, u32bitFMT, &iid);
      while (!feof(inFile)) {
        iids[iid] = true;
        fscanf(inFile, u32bitFMT, &iid);
      }
    }
  };

private:
  bool               isPolishes;
  sim4polishReader  *inPolishes;
  FILE              *inFile;
  FILE              *otFile;
  u32bit             iid;
  bool              *iids;
};


int
main(int argc, char **argv) {
  u32bit                iidRWlen = 0;
  u32bit                iidRWmax = 128;
  iidReaderWriter     **iidRW    = new iidReaderWriter* [iidRWmax];

  FILE                 *defaultOut = 0L;

  seqCache             *F = 0L;
  seqInCore            *S = 0L;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-P") == 0) {
      iidRW[iidRWlen++] = new iidReaderWriter(argv[arg+1], argv[arg+2], true);
      arg+=2;
    } else if (strcmp(argv[arg], "-I") == 0) {
      iidRW[iidRWlen++] = new iidReaderWriter(argv[arg+1], argv[arg+2], false);
      arg+=2;
    } else if (strcmp(argv[arg], "-O") == 0) {
      errno = 0;
      defaultOut = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "Can't open '%s': %s\n", argv[arg], strerror(errno)), exit(1);
    } else if (strcmp(argv[arg], "-i") == 0) {
      F = new seqCache(argv[++arg]);
    } else {
      fprintf(stderr, "ESTmapper utility function -- not for human use.\n");
      exit(1);
    }
    arg++;
  }
  if ((iidRWlen == 0) || (defaultOut == 0L) || (F == 0L)) {
    fprintf(stderr, "spec error.\n");
    exit(1);
  }

  for (u32bit i=0; i<iidRWlen; i++)
    iidRW[i]->load(F->getNumberOfSequences());

  for (u32bit sid=0; ((S = F->getSequenceInCore(sid)) != 0L); sid++) {
    bool     found = false;
    u32bit   iid   = S->getIID();

    for (u32bit i=0; i<iidRWlen; i++) {
      if (iidRW[i]->thisIID(iid)) {
        found = true;
        iidRW[i]->writeSequence(S);
        break;
      }
    }

    if (found == false)
      fprintf(defaultOut, "%s\n%s\n", S->header(), S->sequence());

    delete S;
  }

  return(0);
}

