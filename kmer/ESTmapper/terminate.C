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

    errno = 0;
    inFile = fopen(infile, "r");
    if (errno)
      fprintf(stderr, "iidReaderWriter-- can't open '%s': %s\n", infile, strerror(errno)), exit(1);

    errno = 0;
    otFile = fopen(otfile, "w");
    if (errno)
      fprintf(stderr, "iidReaderWriter-- can't open '%s': %s\n", otfile, strerror(errno)), exit(1);

    iids = 0L;

    //nextIID();
  };

  ~iidReaderWriter() {
    fclose(inFile);
    fclose(otFile);
  };

  void     nextIID(void) {
    if (isPolishes) {
      sim4polish *p = s4p_readPolish(inFile);
      if (p) {
        if (p->estID < iid) {
          fprintf(stderr, "ERROR!  Polishes not sorted by cDNA id!\n");
          exit(1);
        }
        iid = p->estID;
        s4p_destroyPolish(p);
      }
    } else {
      fscanf(inFile, u32bitFMT, &iid);
    }
  };

  bool     thisIID(u32bit targetiid) {
    bool  result = false;

    if (iids) {
      result = iids[targetiid];
    } else {
      while ((iid < targetiid) && (!feof(inFile))) {
        nextIID();
        if (feof(inFile))
          iid = ~u32bitZERO;
      }
      result = iid == targetiid;
    }
    return(result);
  };

  void     writeSequence(seqInCore *S) {
    fprintf(otFile, "%s\n%s\n", S->header(), S->sequence());
  };

  void     load(u32bit maxiid) {
    iids = new bool [maxiid];

    for (u32bit i=0; i<maxiid; i++)
      iids[i] = false;

    //  Mostly, nextIID(), where we stuff the value of iid into an array
    //
    if (isPolishes) {
      sim4polish *p = s4p_readPolish(inFile);
      while (p) {
        iids[p->estID] = true;
        s4p_destroyPolish(p);
        p = s4p_readPolish(inFile);
      }
    } else {
      fscanf(inFile, u32bitFMT, &iid);
      while (!feof(inFile)) {
        iids[iid] = true;
        fscanf(inFile, u32bitFMT, &iid);
      }
    }

#ifdef OLDLOAD
    //  A very simple loader, but it still assumes the IIDs are sorted.
    //
    for (u32bit i=0; i<maxiid; i++) {
      if ((i & 0xfff) == 0)
        fprintf(stderr, "loading "u32bitFMT"\r", i);
      iidsload[i] = thisIID(i);
    }
    iids = iidsload;
    fprintf(stderr, "\n");
#endif


  };

private:
  bool           isPolishes;
  FILE          *inFile;
  FILE          *otFile;
  u32bit         iid;
  bool          *iids;
};


int
main(int argc, char **argv) {

  if (argc < 2)
    fprintf(stderr, "ESTmapper utility function -- not for human use.\n"), exit(1);

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

  if (true)
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

