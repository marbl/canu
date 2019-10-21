
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
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "kmers.H"
#include "system.H"
#include "sequence.H"
#include "bits.H"


#define OP_NONE       0
#define OP_DUMP       1
#define OP_EXISTENCE  2
#define OP_INCLUDE    3
#define OP_EXCLUDE    4



void
dumpExistence(dnaSeqFile           *sf,
              kmerCountExactLookup *kl) {
  uint32   nameMax = 0;
  char    *name    = NULL;
  uint64   seqLen  = 0;
  uint64   seqMax  = 0;
  char    *seq     = NULL;
  uint8   *qlt     = NULL;

  char     fString[64];
  char     rString[64];

  bool     fExists = false, rExists = false;
  uint64   fValue  = 0,     rValue  = 0;

  while (sf->loadSequence(name, nameMax, seq, qlt, seqMax, seqLen)) {
    kmerIterator  kiter(seq, seqLen);

    while (kiter.nextMer()) {
      fExists = kl->exists(kiter.fmer(), fValue);
      rExists = kl->exists(kiter.rmer(), rValue);

      fprintf(stdout, "%s\t%lu\t%c\t%s\t%lu\t%s\t%lu\n",
              name,
              kiter.position(),
              (fExists || rExists) ? 'T' : 'F',
              kiter.fmer().toString(fString), fValue,
              kiter.rmer().toString(rString), rValue);
    }
  }

  delete [] name;
  delete [] seq;
  delete [] qlt;
}



void
reportExistence(dnaSeqFile           *sf,
                kmerCountExactLookup *kl) {
  uint32   nameMax = 0;
  char    *name    = NULL;
  uint64   seqLen  = 0;
  uint64   seqMax  = 0;
  char    *seq     = NULL;
  uint8   *qlt     = NULL;

  while (sf->loadSequence(name, nameMax, seq, qlt, seqMax, seqLen)) {
    kmerIterator  kiter(seq, seqLen);

    uint64   nKmer      = 0;
    uint64   nKmerFound = 0;

    while (kiter.nextMer()) {
      nKmer++;

      if ((kl->value(kiter.fmer()) > 0) ||
          (kl->value(kiter.rmer()) > 0))
        nKmerFound++;
    }
      
    fprintf(stdout, "%s\t%lu\t%lu\t%lu\n", name, nKmer, kl->nKmers(), nKmerFound);
  }

  delete [] name;
  delete [] seq;
  delete [] qlt;
}

void
include( dnaSeqFile           *sf,
         dnaSeqFile           *sf2,
         kmerCountExactLookup *kl,
         const char           *r2name) {
  uint32   nameMax  = 0;
  uint32   nameMax2 = 0;
  char    *name     = NULL;
  char    *name2    = NULL;
  uint64   seqLen   = 0;
  uint64   seqLen2  = 0;
  uint64   seqMax   = 0;
  uint64   seqMax2  = 0;
  char    *seq      = NULL;
  char    *seq2     = NULL;
  uint8   *qlt      = NULL;
  uint8   *qlt2     = NULL;

  uint64   nReads        = 0;
  uint64   nReadsFound   = 0;
  uint64   nKmer         = 0;
  uint64   nKmerFound    = 0;

  // output file for R2
  compressedFileWriter W(r2name);
  FILE   *r2_file = NULL;
  if (sf2 != NULL) {
    r2_file = *W;
  }

  while (sf->loadSequence(name, nameMax, seq, qlt, seqMax, seqLen)) {
    kmerIterator  kiter(seq, seqLen);

    nReads++;
    nKmerFound = 0;

    while (kiter.nextMer()) {
      nKmer++;

      if ((kl->value(kiter.fmer()) > 0) ||
          (kl->value(kiter.rmer()) > 0))
        nKmerFound++;
    }

    if (sf2 != NULL) {
      sf2->loadSequence(name2, nameMax2, seq2, qlt2, seqMax2, seqLen2);
      kmerIterator  kiter2(seq2, seqLen2);

      while (kiter2.nextMer()) {
        nKmer++;

        if ((kl->value(kiter2.fmer()) > 0) ||
            (kl->value(kiter2.rmer()) > 0))
          nKmerFound++;
      }
    }

    if (nKmerFound > 0) {
      if (qlt[0] == 0) {	// is fasta
        fprintf(stdout , ">%s %lu\n%s\n", name , nKmerFound, seq);
        if (sf2 != NULL) fprintf(r2_file, ">%s %lu\n%s\n", name2, nKmerFound, seq2);
      }
      else			// is fastq
      {
        fprintf(stdout , "@%s %lu\n%s\n+\n%s\n", name , nKmerFound, seq , qlt);
        if (sf2 != NULL) fprintf(r2_file, "@%s %lu\n%s\n+\n%s\n", name2, nKmerFound, seq2, qlt2);
      }
      nReadsFound++;
    }
  }

  fprintf(stderr, "\nIncluding %lu reads (or read pairs) out of %lu.\n", nReadsFound, nReads);

  delete [] name;
  delete [] seq;
  delete [] qlt;

  delete [] name2;
  delete [] seq2;
  delete [] qlt2;
}

void
exclude(dnaSeqFile           *sf,
        dnaSeqFile           *sf2,
        kmerCountExactLookup *kl,
        const char*           r2name) {

  uint32   nameMax = 0;
  uint32   nameMax2= 0;
  char    *name    = NULL;
  char    *name2   = NULL;
  uint64   seqLen  = 0;
  uint64   seqLen2 = 0;
  uint64   seqMax  = 0;
  uint64   seqMax2 = 0;
  char    *seq     = NULL;
  char    *seq2    = NULL;
  uint8   *qlt     = NULL;
  uint8   *qlt2    = NULL;

  uint64   nReads = 0;
  uint64   nReadsFound = 0;
  uint64   nKmerFound = 0;

  // output file for R2
  compressedFileWriter W(r2name);
  FILE    *r2_file  = NULL;

  if (sf2 != NULL) {
    r2_file    = *W;
  }


  while (sf->loadSequence(name, nameMax, seq, qlt, seqMax, seqLen)) {
    kmerIterator  kiter(seq, seqLen);
    nReads++;
    nKmerFound = 0;

    while (kiter.nextMer()) {
      if ((kl->value(kiter.fmer()) > 0) ||
          (kl->value(kiter.rmer()) > 0)) {
        nKmerFound++;
      }
    }

    if (sf2 != NULL) {
      sf2->loadSequence(name2, nameMax2, seq2, qlt2, seqMax2, seqLen2);
      kmerIterator  kiter2(seq2, seqLen2);

      while (kiter2.nextMer()) {
        if ((kl->value(kiter2.fmer()) > 0) ||
            (kl->value(kiter2.rmer()) > 0)) {
          nKmerFound++;
        }
      }
    }

    if (nKmerFound == 0) {
       if (qlt[0] == 0) {	// is fasta
          fprintf(stdout, ">%s\n%s\n", name, seq);
          if (sf2 != NULL) fprintf(r2_file, ">%s\n%s\n", name2, seq2);
       }
       else			// is fastq
       {
          fprintf(stdout, "@%s\n%s\n+\n%s\n", name, seq, qlt);
          if (sf2 != NULL) fprintf(r2_file, "@%s\n%s\n+\n%s\n", name2, seq2, qlt2);
       }
       nReadsFound++;
    }
  }

  fprintf(stderr, "\nIncluding %lu reads (or read pairs) out of %lu.\n", nReadsFound, nReads);

  delete [] name;
  delete [] seq;
  delete [] qlt;
  delete [] name2;
  delete [] seq2;
  delete [] qlt2;
}

int
main(int argc, char **argv) {
  char   *inputSeqName = NULL;
  char   *inputSeqName2 = NULL;
  char   *inputDBname  = NULL;
  char   *r2name       = NULL;
  uint64  minV         = 0;
  uint64  maxV         = UINT64_MAX;
  uint32  threads      = omp_get_max_threads();
  uint32  memory       = 0;
  uint32  reportType   = OP_NONE;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-sequence") == 0) {   //  INPUT READS and RANGE TO PROCESS
      inputSeqName = argv[++arg];

    } else if (strcmp(argv[arg], "-sequence2") == 0) {
      inputSeqName2 = argv[++arg];

    } else if (strcmp(argv[arg], "-mers") == 0) {
      inputDBname = argv[++arg];

    } else if (strcmp(argv[arg], "-min") == 0) {
      minV = strtouint64(argv[++arg]);

    } else if (strcmp(argv[arg], "-max") == 0) {
      maxV = strtouint64(argv[++arg]);

    } else if (strcmp(argv[arg], "-threads") == 0) {
      threads = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-memory") == 0) {
      memory = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-dump") == 0) {
      reportType = OP_DUMP;

    } else if (strcmp(argv[arg], "-existence") == 0) {
      reportType = OP_EXISTENCE;

    } else if (strcmp(argv[arg], "-include") == 0) {
      reportType = OP_INCLUDE;

    } else if (strcmp(argv[arg], "-exclude") == 0) {
      reportType = OP_EXCLUDE;

    } else if (strcmp(argv[arg], "-r2") == 0) {
      r2name     = argv[++arg];

    }else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (inputSeqName == NULL)
    err.push_back("No input sequences (-sequence) supplied.\n");
  if (inputDBname == NULL)
    err.push_back("No query meryl database (-mers) supplied.\n");
  if (reportType == OP_NONE)
    err.push_back("No report-type (-existence, etc) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s <report-type> -sequence <input.fasta> -mers <input.meryl> [-sequence2 <input.fasta> -r2 <output.r2>]\n", argv[0]);
    fprintf(stderr, "  Query the kmers in meryl database <input.meryl> with the sequences\n");
    fprintf(stderr, "  in <input.fasta> (both FASTA and FASTQ supported, file can be compressed).\n");
    fprintf(stderr, "  -sequence2 and -r2 is only required for -include and -exclude\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  The meryl database can be filtered by value.  More advanced filtering\n");
    fprintf(stderr, "  requires a new database to be constructed using meryl.\n");
    fprintf(stderr, "    -min   m    Ignore kmers with value below m\n");
    fprintf(stderr, "    -max   m    Ignore kmers with value above m\n");
    fprintf(stderr, "    -threads t  Number of threads to use when constructing lookup table.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Memory usage can be limited, within reason, by sacrificing kmer lookup\n");
    fprintf(stderr, "  speed.  If the lookup table requires more memory than allowed, the program\n");
    fprintf(stderr, "  exits with an error.\n");
    fprintf(stderr, "    -memory m   Don't use more than m GB memory\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Exactly one report type must be specified.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -existence     Report a tab-delimited line for each sequence showing\n");
    fprintf(stderr, "                 the number of kmers in the sequence, in the database,\n");
    fprintf(stderr, "                 and in both.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "     output:  seqName <tab> mersInSeq <tab> mersInDB <tab> mersInBoth\n");
    fprintf(stderr, "         seqName    - name of the sequence\n");
    fprintf(stderr, "         mersInSeq  - number of mers in the sequence\n");
    fprintf(stderr, "         mersInDB   - number of mers in the meryl database\n");
    fprintf(stderr, "         mersInBoth - number of mers in the sequence that are\n");
    fprintf(stderr, "                      also in the database\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -dump          Report a tab-delimited line reporting each kmer in the input\n");
    fprintf(stderr, "                 sequences, in order, annotated with the value of the kmer in\n");
    fprintf(stderr, "                 the input database.  If the kmer does not exist in the database\n");
    fprintf(stderr, "                 its value will be reported as zero.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "     output:  seqName <tab> seqPos <tab> exists <tab> fwd-mer <tab> fwd-val <tab> rev-mer <tab> rev-val\n");
    fprintf(stderr, "         seqName    - name of the sequence this kmer is from\n");
    fprintf(stderr, "         seqPos     - start position (0-based) of the kmer in the sequence\n");
    fprintf(stderr, "         exists     - 'T' if the kmer exists in the database, 'F' if it does not\n");
    fprintf(stderr, "         fwd-mer    - forward mer sequence\n");
    fprintf(stderr, "         fwd-val    - value of the forward mer in the database\n");
    fprintf(stderr, "         rev-mer    - reverse mer sequence\n");
    fprintf(stderr, "         rev-val    - value of the reverse mer in the database\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -include       Extract sequences *containing* kmers in <input.meryl>.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "     output:  sequence given format (fasta or fastq) with the number of overlapping kmers appended\n");
    fprintf(stderr, "              if pairs of sequences are given, R1 will be stdout and R2 be named as <output.r2>\n");
    fprintf(stderr, "              <output.r2> will be automatically compressed if ends with .gz, .bz2, or xs\n");
    fprintf(stderr, "         seqName    - name of the sequence this kmer is from\n");
    fprintf(stderr, "         mersInBoth - number of mers in both sequence and in the database\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -exclude       Extract sequences *NOT containing* kmers in <input.meryl>.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "     output:  sequence given format (fasta or fastq) without reads containing kmers\n");
    fprintf(stderr, "              if pairs of sequences are given, R1 will be stdout and R2 be named as <output.r2>\n");
    fprintf(stderr, "              <output.r2> will be automatically compressed if ends with .gz, .bz2, or xs\n");
    fprintf(stderr, "         seqName    - name of the sequence this kmer is from\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  omp_set_num_threads(threads);

  //  Open the kmers, build a lookup table.

  fprintf(stderr, "-- Loading kmers from '%s' into lookup table.\n", inputDBname);

  kmerCountFileReader   *merylDB    = new kmerCountFileReader(inputDBname);
  kmerCountExactLookup  *kmerLookup = new kmerCountExactLookup(merylDB, memory, minV, maxV);

  if (kmerLookup->configure() == false) {
    exit(1);
  }

  kmerLookup->load();

  delete merylDB;   //  Not needed anymore.

  //  Open sequences.

  fprintf(stderr, "-- Opening sequences in '%s'.\n", inputSeqName);

  dnaSeqFile            *seqFile    = new dnaSeqFile(inputSeqName);
  dnaSeqFile            *seqFile2   = NULL;

  // only keep seqFile2 when provided
  if (inputSeqName2 != NULL) {
    fprintf(stderr, "-- Opening sequences in '%s'.\n", inputSeqName2);
    seqFile2   = new dnaSeqFile(inputSeqName2);
  }

  //  Do something.

  if (reportType == OP_DUMP) {
    dumpExistence(seqFile, kmerLookup);
  }

  if (reportType == OP_EXISTENCE) {
    reportExistence(seqFile, kmerLookup);
  }

  if (reportType == OP_INCLUDE) {
    include(seqFile, seqFile2, kmerLookup, r2name);
  }

  if (reportType == OP_EXCLUDE) {
    exclude(seqFile, seqFile2, kmerLookup, r2name);
  }

  //  Done!
  fprintf(stderr, "Bye!\n");

  delete seqFile;
  delete seqFile2;
  delete kmerLookup;

  exit(0);
}



