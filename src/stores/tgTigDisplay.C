
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


#include "sqStore.H"
#include "tgStore.H"


int
main(int argc, char **argv) {
  char const   *seqStorName = nullptr;
  stringList    tigFileNames;

  //  Delta Encoded Multi Align Display - DEMAD!
  bool          DEMADenable     = false;
  uint32        DEMADwidth      = 250;
  uint32        DEMADspace      = 10;
  bool          DEMADshowDiffs  = true;

  //  CIGAR Encoded Multi Align Display - CEMAD (or CIGAREMAD).
  bool          CEMADenable     = false;

  //  utgcns 'results' output -> BAM output
  bool          BAMenable       = false;
  char const   *BAMoutputName   = nullptr;
  bool          BAMignoreContigs = false;
  bool          BAMignoreRepeats = false;
  bool          BAMignoreBubbles = false;


  argc = AS_configure(argc, argv, 1);

  std::vector<char const *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    if      (strcmp(argv[arg], "-S") == 0)   { seqStorName = argv[++arg]; }
    //se if (strcmp(argv[arg], "-t") == 0)   { tigFileName = argv[++arg]; }

    else if (strcmp(argv[arg], "-d") == 0)   { DEMADenable = true; }

    //se if (strcmp(argv[arg], "-c") == 0)   { CEMADenable = true; }

    else if (strcmp(argv[arg], "-b") == 0)   { BAMenable = true;   }
    else if (strcmp(argv[arg], "-o") == 0)   { BAMoutputName = argv[++arg]; }

    else if (strcmp(argv[arg], "-L") == 0)
      tigFileNames.load(argv[++arg]);

    else if ((strcmp(argv[arg], "-") == 0) ||                //  A single '-' (for stdin) or
             (merylutil::fileExists(argv[arg]) == true))     //  an existing file
      tigFileNames.add(argv[arg]);

    else {                                                   //  No idea what this crap is.
      char *s = new char [1024];
      snprintf(s, 1024, "ERROR: Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }
  }

  if (seqStorName == nullptr)     err.push_back("ERROR: No seqStore (-S) supplied.\n");
  if (tigFileNames.size() == 0)   err.push_back("ERROR: No utgcns 'results' output files supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -S seqStore [options] utgcnsResultFile ...\n", argv[0]);
    fprintf(stderr, "  Display multialignments or creates a single BAM output file\n");
    fprintf(stderr, "   for tigs in utgcnsResultFile.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S s           input seqStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -d             display multialign from delta-encoded alignments\n");
    fprintf(stderr, "\n");
    //rintf(stderr, "  -c             display multialign from cigar-encoded alignments\n");
    //rintf(stderr, "                   (not implemented)\n");
    //rintf(stderr, "\n");
    fprintf(stderr, "  -b             convert tigs in input files to a single bam output\n");
    fprintf(stderr, "  -o             bam output file (stdout if not supplied)\n");
    fprintf(stderr, "\n");
    return 1;
  }


  //  Force use of corrected non-homopolymer compressed reads.
  sqRead_setDefaultVersion(sqRead_corrected | sqRead_trimmed | sqRead_normal);

  sqStore   *seqStore  = new sqStore(seqStorName);


  //
  //  Show a the classic Celera Assembler multi-align display with
  //  alignments regenerated from the delta-encoded alignments.
  //
  if (DEMADenable == true) {
    for (uint64 tf=0; tf<tigFileNames.size(); tf++) {
      readBuffer  *tigBuffer = new readBuffer(tigFileNames[tf]);
      tgTig        tig;

      while (tig.loadFromBuffer(tigBuffer))
        tig.display(stdout, seqStore, DEMADwidth, DEMADspace, DEMADshowDiffs);

      delete tigBuffer;
    }
  }


  //
  //  Show ... nothing.
  //
  if (CEMADenable == true) {
  }


  //
  //  Convert the set of input tig files into a single BAM output file.
  //    (see also tgTig::dumpBAM() in tgTig.C)
  //
  if (BAMenable   == true) {
    sam_hdr_t *outBAMhp = sam_hdr_init();
    samFile   *outBAMfp = hts_open(BAMoutputName, "wb");
    if (outBAMfp == NULL) {
      fprintf(stderr, "Failed to open BAM output file '%s': %s\n", BAMoutputName, strerror(errno));
      exit(1);
    }

    sam_hdr_add_line(outBAMhp, "HD", "VN", SAM_FORMAT_VERSION, nullptr);
    sam_hdr_add_pg  (outBAMhp, "utgcns", "VN", MERYL_UTILITY_VERSION, nullptr);

    //  The header needs to know how many 'reference' sequences (tigs) are
    //  in this file, their length and their name.  It annoyingly wants this
    //  as old-schoold malloc()d data, but for our sanity, the lists are built
    //  using C++ allocated memory, then converted to malloc() storage at the end.

    std::vector<uint32>   tLen;
    std::vector<char *>   tName;

    fprintf(stderr, "-- READING TIG NAMES\n");

    for (uint64 tf=0; tf<tigFileNames.size(); tf++) {              //  Load tig names and lengths
      readBuffer  *tigBuffer = new readBuffer(tigFileNames[tf]);   //  into C++-friendly storage.
      tgTig        tig;
      char         nam[16] = {0};

      fprintf(stderr, "--   %s\n", tigFileNames[tf]);

      while (tig.loadFromBuffer(tigBuffer)) {
        if ((BAMignoreRepeats == true) && (tig._suggestRepeat == true))    continue;
        if ((BAMignoreBubbles == true) && (tig._suggestBubble == true))    continue;

        sprintf(nam, "tig%08u", tig.tigID());
        tLen.push_back(tig.length());
        tName.push_back(strdup(nam));
      }

      delete tigBuffer;
    }

    outBAMhp->n_targets      = tLen.size();                   //  Setup header for N ref seqs.
    outBAMhp->target_len     = (uint32_t *)malloc(outBAMhp->n_targets * sizeof(uint32_t));
    outBAMhp->target_name    = (char    **)malloc(outBAMhp->n_targets * sizeof(char *));

    for (uint32 ii=0; ii<outBAMhp->n_targets; ii++) {         //  Copy names to the header.
      outBAMhp->target_len[ii]  = tLen[ii];
      outBAMhp->target_name[ii] = tName[ii];
    } 

    int ret = sam_hdr_write(outBAMfp, outBAMhp);   //  Only works for BAM; sam needs some extra field
    if (ret < 0) {                                 //  set so that hdr_write will output SQ records.
      fprintf(stderr, "Failed to write header to BAM file!\n");
      exit(1);
    }

    //  With a valid header, iterate over all the tigs again and add alignments records for
    //  each read.

    fprintf(stderr, "-- COPYING READ-to-TIG ALIGNMENTS\n");

    uint32   bamTarget = uint32max;   //  Rolls over to zero on the first tig.

    for (uint64 tf=0; tf<tigFileNames.size(); tf++) {
      readBuffer  *tigBuffer = new readBuffer(tigFileNames[tf]);
      sqRead       read;
      tgTig        tig;

      while (tig.loadFromBuffer(tigBuffer)) {
        if ((BAMignoreRepeats == true) && (tig._suggestRepeat == true))    continue;
        if ((BAMignoreBubbles == true) && (tig._suggestBubble == true))    continue;

        bamTarget++;

        if (tig.numberOfChildren() > 1)
          fprintf(stderr, "--   %s: tig%08u %8u reads -> bam target #%u\n", tigFileNames[tf], tig.tigID(), tig.numberOfChildren(), bamTarget);

        for (uint32 rr=0; rr<tig.numberOfChildren(); rr++) {
          tgPosition  *child = tig.getChild(rr);

          const char  *cigar         = tig.getChildCIGAR(rr);
          size_t       cigarArrayLen = (cigar == nullptr) ? 0 : strlen(cigar);

          uint32_t    *cigarArray    = (cigarArrayLen == 0) ? nullptr : new uint32_t [cigarArrayLen];
          ssize_t      cigarLenS     = (cigarArrayLen == 0) ? 0       : sam_parse_cigar(cigar, nullptr, &cigarArray, &cigarArrayLen);

          seqStore->sqStore_getRead(child->ident(), &read);

          char const  *readName      = read.sqRead_name();

          uint32       readlen       = read.sqRead_length() - child->_askip - child->_bskip;
          char        *readseq       = read.sqRead_sequence();

          if (child->isReverse() == true)                                       //  If reverse, get a copy of
            readseq = reverseComplementCopy(readseq + child->_bskip, readlen);  //  the RC of the read, otherwise
          else                                                                  //  get a copy of the forward seq.
            readseq = duplicateString(readseq + child->_askip);                 //  This duplicates unitigConsensus's
          readseq[readlen] = 0;                                                 //  addRead / abSequence.

          bam1_t      *bamRecord     = bam_init1();
          int          flags         = 0;

          flags |= (cigarArrayLen == 0) ? BAM_FUNMAP : 0;
          flags |= (child->isForward()) ? 0 : BAM_FREVERSE;

          int ret = bam_set1(bamRecord,                           //  Record to add to
                             strlen(readName), readName,          //  Name of entry to add
                             flags,                               //  Flags.
                             bamTarget,                           //  Target ID of the ref in the header
                             child->min(),                        //  Start position on target, 0-based
                             255,                                 //  Mapping Quality not available
                             cigarLenS, cigarArray,               //  Number of CIGAR operations, and operations
                             -1, -1,                              //  Position (target, begin) of next read in template
                             0, /*_tig->length(),*/               //  Length of template
                             readlen, readseq, nullptr,           //  Read length, sequence and quality values
                             0);                                  //  Space to reserve for auxiliary data
          if (ret < 0) { 
            fprintf(stderr, "Failed to create bam record:\n");
            fprintf(stderr, "  read %u %s\n", child->ident(), readName);
            fprintf(stderr, "  length=%u  askip=%u  bskip=%u\n", read.sqRead_length(), child->_askip, child->_bskip);
            fprintf(stderr, "  %s\n", cigar);
            exit(1);
          }

          ret = sam_write1(outBAMfp, outBAMhp, bamRecord);
          if (ret < 0) {
            fprintf(stderr, "Failed to write sam record! %s\n", strerror(errno));
            exit(1);
          }
          bam_destroy1(bamRecord);

          delete [] readseq;   //  A copy of either the rev-comp or forward sequence.
          delete [] cigarArray;
        }  //  Over all reads in a tig.
      }    //  Over all tigs in the file

      delete tigBuffer;
    }      //  Over all tigfiles.

    sam_close(outBAMfp);
    sam_hdr_destroy(outBAMhp);

    fprintf(stderr, "-- Success!  Bye.\n");
  }


  delete seqStore;

  return 0;
}
