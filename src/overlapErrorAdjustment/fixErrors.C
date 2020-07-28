
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "runtime.H"
#include "types.H"
#include "arrays.H"

#include "sqStore.H"

#include "correctOverlaps.H"
#include "correctionOutput.H"


//  Technically, private to correctOverlaps-Correct_Frags.C
void
correctRead(uint32 curID,
            char *fseq, uint32 &fseqLen, Adjust_t *fadj, uint32 &fadjLen,
            const char *oseq, uint32  oseqLen,
            Correction_Output_t  *C,
            uint64               &Cpos,
            uint64                Clen,
            uint64               *changes);


int
main(int argc, char **argv) {
  char const *seqStoreIn = nullptr;
  char const *redIn      = nullptr;

  char const *outName    = nullptr;

  uint32      bgnID      = 1;
  uint32      endID      = UINT32_MAX;


  argc = AS_configure(argc, argv);

  vector<char const *>  err;
  int                   arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqStoreIn = argv[++arg];

    } else if (strcmp(argv[arg], "-red") == 0) {
      redIn = argv[++arg];

    } else if (strcmp(argv[arg], "-r") == 0) {
      decodeRange(argv[++arg], bgnID, endID);

    } else if (strcmp(argv[arg], "-O") == 0) {
      outName = argv[++arg];

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (seqStoreIn == nullptr)   err.push_back("No sequence store (-S option) supplied.\n");
  if (outName    == nullptr)   err.push_back("No output file (-O option) supplied.\n");
  if (bgnID      >  endID)     err.push_back("Inconsistent read range (-r option): begin ID must be no larger than end ID.");
  if (redIn      == nullptr)   err.push_back("No findErrors .red file (-red option) supplied.");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -S orig.seqStore -red f.red -O out.fasta [-r bgnID-endID]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S   seqStore           path to a sequence store\n");
    fprintf(stderr, "  -red f.red              path to the output .red file from findErrors\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -r   bgnID-endID        output only reads bgnID through endID, inclusive\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -O   out.fasta[.gz]     write FASTA to file 'out.fasta, optionally compressing\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  //  Open inputs check paramters.

  sqStore              *seqStore = new sqStore(seqStoreIn);

  if (seqStore->sqStore_lastReadID() < bgnID)
    bgnID = seqStore->sqStore_lastReadID();

  if (seqStore->sqStore_lastReadID() < endID)
    endID = seqStore->sqStore_lastReadID();

  //  Open outputs.

  compressedFileWriter *outFASTA = new compressedFileWriter(outName);

  //  Load corrections.

  memoryMappedFile     *Cfile = new memoryMappedFile(redIn);
  Correction_Output_t  *C     = (Correction_Output_t *)Cfile->get();
  uint64                Cpos  = 0;
  uint64                Clen  = Cfile->length() / sizeof(Correction_Output_t);

  sqRead                read;

  uint32                corBasesLen = 0;
  uint32                corBasesMax = 128 * 1024;
  char                 *corBases    = new char [corBasesMax];

  uint32                adjustsLen  = 0;
  Adjust_t             *adjusts     = new Adjust_t [Clen];

  uint64                totalBases  =  0;
  uint64                changes[12] = {0};

  //  Scan reads and corrections to figure out how much space to allocate.

  //  Iterate over reads.

  for (uint32 curID=bgnID; curID <= endID; curID++) {
    if (seqStore->sqStore_getReadLength(curID) == 0)
      continue;

    seqStore->sqStore_getRead(curID, &read);

    totalBases += read.sqRead_length();

    if (read.sqRead_length() > corBasesMax)
      resizeArray(corBases, 0, corBasesMax, read.sqRead_length() + 1024, resizeArray_doNothing);

    correctRead(curID,
                corBases,
                corBasesLen,
                adjusts,
                adjustsLen,
                read.sqRead_sequence(),
                read.sqRead_length(),
                C,
                Cpos,
                Clen,
                changes);

    AS_UTL_writeFastA(outFASTA->file(), corBases, corBasesLen, 0, ">%d\n", curID);

    corBasesLen = 0;   //  correctRead() appends to the corBases array, which we don't want.
  }


  delete [] corBases;
  delete [] adjusts;

  delete    Cfile;

  delete    outFASTA;

  delete    seqStore;

  fprintf(stderr, "Corrected " F_U64 " bases with " F_U64 " substitutions, " F_U64 " deletions and " F_U64 " insertions.\n",
          totalBases,
          changes[A_SUBST] + changes[C_SUBST] + changes[G_SUBST] + changes[T_SUBST],
          changes[DELETE],
          changes[A_INSERT] + changes[C_INSERT] + changes[G_INSERT] + changes[T_INSERT]);

  return(0);
}


