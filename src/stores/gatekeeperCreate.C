
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
 *  This file is derived from:
 *
 *    src/AS_GKP/gkpStoreCreate.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2012-FEB-06 to 2013-AUG-01
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-OCT-09 to 2015-AUG-10
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-NOV-08
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "gkStore.H"
#include "findKeyAndValue.H"
#include "AS_UTL_fileIO.H"


#undef  UPCASE  //  Don't convert lowercase to uppercase, special case for testing alignments.
#define UPCASE  //  Convert lowercase to uppercase.  Probably needed.



//  Canu doesn't really use QVs as of late 2015.  Only consensus is aware of them, but since read
//  correction doesn't output QVs, the QVs that consensus sees are totally bogus.  So, disable
//  storage of them - gatekeeperCreate will store a default QV for every read, and consensus will
//  process as usual.
//
#define  DO_NOT_STORE_QVs



//  Support fastq of fasta, even in the same file.
//  Eventually want to support bax.h5 natively.


uint32  validSeq[256] = {0};


uint32
loadFASTA(char                 *L,
          char                 *H,
          char                 *S,
          uint32               &Slen,
          char                 *Q,
          compressedFileReader *F,
          FILE                 *errorLog,
          uint32               &nWARNS) {
  uint32  nLines = 0;     //  Lines read from the input
  uint32  nBases = 0;     //  Bases read from the input, used for reporting errors
  bool    valid  = true;

  //  We've already read the header.  It's in L.  But we want to use L to load the sequence, so the
  //  header is copied to H.  We need to do a copy anyway - we need to return the next header in L.
  //  We can either copy L to H now, or use H to read lines and copy H to L at the end.  We're stuck
  //  either way.

  strcpy(H, L + 1);

  //  Clear the sequence.

  S[0] = 0;
  Q[0] = 0;  //  Sentinel to tell gatekeeper to use the fixed QV value

  Slen = 0;

  //  Load sequence.  This is a bit tricky, since we need to peek ahead
  //  and stop reading before the next header is loaded.  Instead, we read the
  //  next line into what we'd read the header into outside here.

  fgets(L, AS_MAX_READLEN+1, F->file());  nLines++;
  chomp(L);

  //  Catch empty reads - reads with no sequence line at all.

  if (L[0] == '>') {
    fprintf(errorLog, "read '%s' is empty.\n", H);
    nWARNS++;
    return(nLines);
  }

  //  Copy in the sequence, as long as it is valid sequence.  If any invalid letters
  //  are found, set the base to 'N'.

  uint32  baseErrors = 0;

  while ((!feof(F->file())) && (L[0] != '>')) {
    nBases += strlen(L);  //  Could do this in the loop below, but it makes it ugly.

    for (uint32 i=0; (Slen < AS_MAX_READLEN) && (L[i] != 0); i++) {
      switch (L[i]) {
#ifdef UPCASE
        case 'a':   S[Slen] = 'A';  break;
        case 'c':   S[Slen] = 'C';  break;
        case 'g':   S[Slen] = 'G';  break;
        case 't':   S[Slen] = 'T';  break;
#else
        case 'a':   S[Slen] = 'a';  break;
        case 'c':   S[Slen] = 'c';  break;
        case 'g':   S[Slen] = 'g';  break;
        case 't':   S[Slen] = 't';  break;
#endif
        case 'A':   S[Slen] = 'A';  break;
        case 'C':   S[Slen] = 'C';  break;
        case 'G':   S[Slen] = 'G';  break;
        case 'T':   S[Slen] = 'T';  break;
        case 'n':   S[Slen] = 'N';  break;
        case 'N':   S[Slen] = 'N';  break;
        default:
          baseErrors++;
          S[Slen]   = 'N';
          break;
      }

      Slen++;
    }

    //  Grab the next line.  It should be more sequence, or the next header, or eof.
    //  The last two are stop conditions for the while loop.

    L[0] = 0;

    fgets(L, AS_MAX_READLEN+1, F->file());  nLines++;
    chomp(L);
  }

  //  Terminate the sequence.

  S[Slen] = 0;

  //  Report errors.

  if (baseErrors > 0) {
    fprintf(errorLog, "read '%s' has "F_U32" invalid base%s.  Converted to 'N'.\n",
            H, baseErrors, (baseErrors > 1) ? "s" : "");
    nWARNS++;
  }

  if (Slen == 0) {
    fprintf(errorLog, "read '%s' is empty.\n", H);
    nWARNS++;
  }

  if (Slen != nBases) {
    fprintf(errorLog, "read '%s' is too long; contains %u bases, but we can only handle %u.\n", H, nBases, AS_MAX_READLEN);
    nWARNS++;
  }

  //  Do NOT clear L, it contains the next header.

  return(nLines);
}



uint32
loadFASTQ(char                 *L,
          char                 *H,
          char                 *S,
          uint32               &Slen,
          char                 *Q,
          compressedFileReader *F,
          FILE                 *errorLog,
          uint32               &nWARNS) {

  //  We've already read the header.  It's in L.

  strcpy(H, L + 1);

  //  Load sequence.

  S[0] = 0;
  Slen = 0;

  S[AS_MAX_READLEN+1-2] = 0;  //  If this is ever set, the read is probably longer than we can support.
  S[AS_MAX_READLEN+1-1] = 0;  //  This will always be zero; fgets() sets it.

  Q[AS_MAX_READLEN+1-2] = 0;  //  This too.
  Q[AS_MAX_READLEN+1-1] = 0;

  fgets(S, AS_MAX_READLEN+1, F->file());
  chomp(S);

  //  Check for long reads.  If found, read the rest of the line, and report an error.  The -1 (in
  //  the print) is because fgets() and strlen() will count the newline, which isn't a base.

  if ((S[AS_MAX_READLEN+1-2] != 0) && (S[AS_MAX_READLEN+1-2] != '\n')) {
    char    *overflow = new char [1048576];
    uint32   nBases   = AS_MAX_READLEN;

    do {
      overflow[1048576-2] = 0;
      overflow[1048576-1] = 0;
      fgets(overflow, 1048576, F->file());
      nBases += strlen(overflow);
    } while (overflow[1048576-2] != 0);

    fprintf(errorLog, "read '%s' is too long; contains %u bases, but we can only handle %u.\n", H, nBases-1, AS_MAX_READLEN);

    nWARNS++;

    delete [] overflow;
  }

  //  Check for and correct invalid bases.

  uint32 baseErrors = 0;

  for (uint32 i=0; (Slen < AS_MAX_READLEN) && (S[i] != 0); i++) {
    switch (S[i]) {
#ifdef UPCASE
      case 'a':   S[i] = 'A';  break;
      case 'c':   S[i] = 'C';  break;
      case 'g':   S[i] = 'G';  break;
      case 't':   S[i] = 'T';  break;
#else
      case 'a':                break;
      case 'c':                break;
      case 'g':                break;
      case 't':                break;
#endif
      case 'A':                break;
      case 'C':                break;
      case 'G':                break;
      case 'T':                break;
      case 'n':   S[i] = 'N';  break;
      case 'N':                break;
      default:
        S[i] = 'N';
        Q[i] = '!';  //  QV=0, ASCII=33
        baseErrors++;
        break;
    }

    Slen++;
  }

  if (baseErrors > 0) {
        fprintf(errorLog, "read '%s' has "F_U32" invalid base%s.  Converted to 'N'.\n",
                L, baseErrors, (baseErrors > 1) ? "s" : "");
    nWARNS++;
  }

  //  Load the qv header, and then load the qvs themselves over the header.

  Q[0] = 0;
  fgets(Q, AS_MAX_READLEN+1, F->file());
  fgets(Q, AS_MAX_READLEN+1, F->file());
  chomp(Q);

  //  As with the base, we need to suck in the rest of the longer-than-allowed QV string.  But we don't need to report it
  //  or do anything fancy, just advance the file pointer.

  if ((Q[AS_MAX_READLEN-1] != 0) && (Q[AS_MAX_READLEN-1] != '\n')) {
    char    *overflow = new char [1048576];

    do {
      overflow[1048576-2] = 0;
      overflow[1048576-1] = 0;
      fgets(overflow, 1048576, F->file());
    } while (overflow[1048576-2] != 0);

    delete [] overflow;
  }

  //  Convert from the (assumed to be) Sanger QVs to plain ol' integers.

  uint32 QVerrors = 0;

#ifndef DO_NOT_STORE_QVs

  for (uint32 i=0; Q[i]; i++) {
    if (Q[i] < '!') {  //  QV=0, ASCII=33
      Q[i] = '!';
      QVerrors++;
    }

    if (Q[i] > '!' + 60) {  //  QV=60, ASCII=93=']'
      Q[i] = '!' + 60;
      QVerrors++;
    }
  }

  if (QVerrors > 0) {
    fprintf(errorLog, "read '%s' has "F_U32" invalid QV%s.  Converted to min or max value.\n",
            L, QVerrors, (QVerrors > 1) ? "s" : "");
    nWARNS++;
  }

#else

  //  If we're not using QVs, just reset the first value to -1.  This is the sentinel that FASTA sequences set,
  //  causing the encoding later to use a fixed QV for all bases.

  Q[0] = 0;

#endif

  //  Clear the lines, so we can load the next one.

  L[0] = 0;

  return(4);  //  FASTQ always reads exactly four lines
}





void
loadReads(gkStore    *gkpStore,
          gkLibrary  *gkpLibrary,
          uint32      gkpFileID,
          uint32      minReadLength,
          FILE       *nameMap,
          FILE       *htmlLog,
          FILE       *errorLog,
          char       *fileName,
          uint32     &nWARNS,
          uint32     &nLOADED,
          uint64     &bLOADED,
          uint32     &nSKIPPED,
          uint64     &bSKIPPED) {
  char    *L = new char [AS_MAX_READLEN + 1];  //  +1.  One for the newline, and one for the terminating nul.
  char    *H = new char [AS_MAX_READLEN + 1];
  char    *S = new char [AS_MAX_READLEN + 1];
  char    *Q = new char [AS_MAX_READLEN + 1];

  uint32   Slen = 0;
  uint32   Qlen = 0;

  uint64   lineNumber = 1;

  fprintf(stderr, "\n");
  fprintf(stderr, "  Loading reads from '%s'\n", fileName);

#if 0
  fprintf(htmlLog, "<tr id='gkpload%u'><td colspan='2'>%s</td></tr>\n", gkpFileID, fileName);
  fprintf(htmlLog, "<tr class='details'><td rowspan='9'>Parameters</td><td>preset=N/A</td></tr>\n");
  fprintf(htmlLog, "<tr class='details'><td>defaultQV=%u</td></tr>\n",            gkpLibrary->gkLibrary_defaultQV());
  fprintf(htmlLog, "<tr class='details'><td>isNonRandom=%s</td></tr>\n",          gkpLibrary->gkLibrary_isNonRandom()          ? "true" : "false");
  fprintf(htmlLog, "<tr class='details'><td>removeDuplicateReads=%s</td></tr>\n", gkpLibrary->gkLibrary_removeDuplicateReads() ? "true" : "false");
  fprintf(htmlLog, "<tr class='details'><td>finalTrim=%s</td></tr>\n",            gkpLibrary->gkLibrary_finalTrim()            ? "true" : "false");
  fprintf(htmlLog, "<tr class='details'><td>removeSpurReads=%s</td></tr>\n",      gkpLibrary->gkLibrary_removeSpurReads()      ? "true" : "false");
  fprintf(htmlLog, "<tr class='details'><td>removeChimericReads=%s</td></tr>\n",  gkpLibrary->gkLibrary_removeChimericReads()  ? "true" : "false");
  fprintf(htmlLog, "<tr class='details'><td>checkForSubReads=%s</td></tr>\n",     gkpLibrary->gkLibrary_checkForSubReads()     ? "true" : "false");
#else
  fprintf(htmlLog, "nam "F_U32" %s\n", gkpFileID, fileName);

  fprintf(htmlLog, "lib preset=N/A");
  fprintf(htmlLog,    " defaultQV=%u",            gkpLibrary->gkLibrary_defaultQV());
  fprintf(htmlLog,    " isNonRandom=%s",          gkpLibrary->gkLibrary_isNonRandom()          ? "true" : "false");
  fprintf(htmlLog,    " removeDuplicateReads=%s", gkpLibrary->gkLibrary_removeDuplicateReads() ? "true" : "false");
  fprintf(htmlLog,    " finalTrim=%s",            gkpLibrary->gkLibrary_finalTrim()            ? "true" : "false");
  fprintf(htmlLog,    " removeSpurReads=%s",      gkpLibrary->gkLibrary_removeSpurReads()      ? "true" : "false");
  fprintf(htmlLog,    " removeChimericReads=%s",  gkpLibrary->gkLibrary_removeChimericReads()  ? "true" : "false");
  fprintf(htmlLog,    " checkForSubReads=%s\n",   gkpLibrary->gkLibrary_checkForSubReads()     ? "true" : "false");
#endif


  compressedFileReader *F = new compressedFileReader(fileName);

  uint32   nFASTAlocal    = 0;  //  number of sequences read from disk
  uint32   nFASTQlocal    = 0;
  uint32   nWARNSlocal    = 0;

  uint32   nLOADEDAlocal  = 0;  //  Sequences actaully loaded into the store
  uint32   nLOADEDQlocal  = 0;

  uint64   bLOADEDAlocal  = 0;
  uint64   bLOADEDQlocal  = 0;

  uint32   nSKIPPEDAlocal = 0;  //  Sequences skipped because they are too short
  uint32   nSKIPPEDQlocal = 0;

  uint64   bSKIPPEDAlocal = 0;
  uint64   bSKIPPEDQlocal = 0;

  fgets(L, AS_MAX_READLEN+1, F->file());
  chomp(L);

  while (!feof(F->file())) {
    bool  isFASTA = false;
    bool  isFASTQ = false;

    if      (L[0] == '>') {
      lineNumber += loadFASTA(L, H, S, Slen, Q, F, errorLog, nWARNSlocal);
      isFASTA = true;
      nFASTAlocal++;
    }

    else if (L[0] == '@') {
      lineNumber += loadFASTQ(L, H, S, Slen, Q, F, errorLog, nWARNSlocal);
      isFASTQ = true;
      nFASTQlocal++;
    }

    else {
      fprintf(errorLog, "invalid read header '%.40s%s' in file '%s' at line "F_U64", skipping.\n",
              L, (strlen(L) > 80) ? "..." : "", fileName, lineNumber);
      L[0] = 0;
      nWARNSlocal++;
    }

    //  If S[0] isn't nul, we loaded a sequence and need to store it.

    if (Slen < minReadLength) {
      fprintf(errorLog, "read '%s' of length "F_U32" in file '%s' at line "F_U64" is too short, skipping.\n",
              H, Slen, fileName, lineNumber);

      if (isFASTA) {
        nSKIPPEDAlocal += 1;
        bSKIPPEDAlocal += Slen;
      }

      if (isFASTQ) {
        nSKIPPEDQlocal += 1;
        bSKIPPEDQlocal += Slen;
      }

      S[0] = 0;
      Q[0] = 0;
    }

    if (S[0] != 0) {
      gkRead     *nr = gkpStore->gkStore_addEmptyRead(gkpLibrary);
      gkReadData *nd = nr->gkRead_encodeSeqQlt(H, S, Q, gkpLibrary->gkLibrary_defaultQV());

      gkpStore->gkStore_stashReadData(nr, nd);

      delete nd;

      if (isFASTA) {
        nLOADEDAlocal += 1;
        bLOADEDAlocal += Slen;
      }

      if (isFASTQ) {
        nLOADEDQlocal += 1;
        bLOADEDQlocal += Slen;
      }

      fprintf(nameMap, F_U32"\t%s\n", gkpStore->gkStore_getNumReads(), H);
    }

    //  If L[0] is nul, we need to load the next line.  If not, the next line is the header (from
    //  the fasta loader).

    if (L[0] == 0) {
      fgets(L, AS_MAX_READLEN+1, F->file());  lineNumber++;
      chomp(L);
    }
  }

  delete    F;

  delete [] Q;
  delete [] S;
  delete [] H;
  delete [] L;

  lineNumber--;  //  The last fgets() returns EOF, but we still count the line.

  //  Write status to the screen

  fprintf(stderr, "    Processed "F_U64" lines.\n", lineNumber);

  fprintf(stderr, "    Loaded "F_U64" bp from:\n", bLOADEDAlocal + bLOADEDQlocal);
  if (nFASTAlocal > 0)
    fprintf(stderr, "      "F_U32" FASTA format reads ("F_U64" bp).\n", nFASTAlocal, bLOADEDAlocal);
  if (nFASTQlocal > 0)
    fprintf(stderr, "      "F_U32" FASTQ format reads ("F_U64" bp).\n", nFASTQlocal, bLOADEDQlocal);

  if (nWARNSlocal > 0)
    fprintf(stderr, "    WARNING: "F_U32" reads issued a warning.\n", nWARNSlocal);

  if (nSKIPPEDAlocal > 0)
    fprintf(stderr, "    WARNING: "F_U32" reads (%0.4f%%) with "F_U64" bp (%0.4f%%) were too short (< "F_U32"bp) and were ignored.\n",
            nSKIPPEDAlocal, 100.0 * nSKIPPEDAlocal / (nSKIPPEDAlocal + nLOADEDAlocal),
            bSKIPPEDAlocal, 100.0 * bSKIPPEDAlocal / (bSKIPPEDAlocal + bLOADEDAlocal),
            minReadLength);

  if (nSKIPPEDQlocal > 0)
    fprintf(stderr, "    WARNING: "F_U32" reads (%0.4f%%) with "F_U64" bp (%0.4f%%) were too short (< "F_U32"bp) and were ignored.\n",
            nSKIPPEDQlocal, 100.0 * nSKIPPEDQlocal / (nSKIPPEDQlocal + nLOADEDQlocal),
            bSKIPPEDQlocal, 100.0 * bSKIPPEDQlocal / (bSKIPPEDQlocal + bLOADEDQlocal),
            minReadLength);

  //  Write status to HTML

#if 0
  fprintf(htmlLog, "<tr class='details'><td rowspan='2'>FASTA</td><td>"F_U32" reads ("F_U64" bp)</td></tr>\n", nLOADEDAlocal, bLOADEDAlocal);
  fprintf(htmlLog, "<tr class='details'><td>"F_U32" reads ("F_U64" bp) were short and not loaded</td></tr>\n", nSKIPPEDAlocal, bSKIPPEDAlocal);

  fprintf(htmlLog, "<tr class='details'><td rowspan='2'>FASTQ</td><td>"F_U32" reads ("F_U64" bp)</td></tr>\n", nLOADEDQlocal, bLOADEDQlocal);
  fprintf(htmlLog, "<tr class='details'><td>"F_U32" reads ("F_U64" bp) were short and not loaded</td></tr>\n", nSKIPPEDQlocal, bSKIPPEDQlocal);

  fprintf(htmlLog, "<tr><td colspan='2'>"F_U32" reads ("F_U64" bp) loaded, "F_U32" reads ("F_U64" bp) skipped, "F_U32" warnings</td></tr>\n",
          nLOADEDAlocal + nLOADEDQlocal, bLOADEDAlocal + bLOADEDQlocal,
          nSKIPPEDAlocal + nSKIPPEDQlocal, bSKIPPEDAlocal + bSKIPPEDQlocal,
          nWARNSlocal);
#else
  fprintf(htmlLog, "dat "F_U32" "F_U64" "F_U32" "F_U64" "F_U32" "F_U64" "F_U32" "F_U64" "F_U32"\n",
          nLOADEDAlocal, bLOADEDAlocal,
          nSKIPPEDAlocal, bSKIPPEDAlocal,
          nLOADEDQlocal, bLOADEDQlocal,
          nSKIPPEDQlocal, bSKIPPEDQlocal,
          nWARNSlocal);
#endif

  //  Add the just loaded numbers to the global numbers

  nWARNS   += nWARNSlocal;

  nLOADED  += nLOADEDAlocal + nLOADEDQlocal;
  bLOADED  += bLOADEDAlocal + bLOADEDQlocal;

  nSKIPPED += nSKIPPEDAlocal + nSKIPPEDQlocal;
  bSKIPPED += bSKIPPEDAlocal + bSKIPPEDQlocal;
};




int
main(int argc, char **argv) {
  char            *gkpStoreName      = NULL;
  char            *outPrefix         = NULL;

  uint32           minReadLength     = 0;

  uint32           firstFileArg      = 0;

  char             errorLogName[FILENAME_MAX];
  char             htmlLogName[FILENAME_MAX];
  char             nameMapName[FILENAME_MAX];

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-minlength") == 0) {
      minReadLength = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "--") == 0) {
      firstFileArg = arg++;
      break;

    } else if (argv[arg][0] == '-') {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      err++;

    } else {
      firstFileArg = arg;
      break;
    }
    arg++;
  }

  if (gkpStoreName == NULL)
    err++;
  if (firstFileArg == 0)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s [...] -o gkpStore\n", argv[0]);
    fprintf(stderr, "  -o gkpStore         create this gkpStore\n");
    fprintf(stderr, "  \n");
    fprintf(stderr, "  -minlength L        discard reads shorter than L\n");
    fprintf(stderr, "  \n");
    fprintf(stderr, "  \n");

    if (gkpStoreName == NULL)
      fprintf(stderr, "ERROR: no gkpStore (-o) supplied.\n");
    if (firstFileArg == 0)
      fprintf(stderr, "ERROR: no input files supplied.\n");

    exit(1);
  }


  gkStore     *gkpStore     = gkStore::gkStore_open(gkpStoreName, gkStore_extend);
  gkRead      *gkpRead      = NULL;
  gkLibrary   *gkpLibrary   = NULL;
  uint32       gkpFileID    = 0;      //  Used for HTML output, an ID for each file loaded.

  uint32       inLineLen    = 1024;
  char         inLine[1024] = { 0 };

  validSeq['a'] = validSeq['c'] = validSeq['g'] = validSeq['t'] = validSeq['n'] = 1;
  validSeq['A'] = validSeq['C'] = validSeq['G'] = validSeq['T'] = validSeq['N'] = 1;

  errno = 0;

  sprintf(errorLogName, "%s/errorLog",    gkpStoreName);
  FILE    *errorLog = fopen(errorLogName, "w");
  if (errno)
    fprintf(stderr, "ERROR:  cannot open error file '%s': %s\n", errorLogName, strerror(errno)), exit(1);

  sprintf(htmlLogName,   "%s/load.dat", gkpStoreName);
  FILE    *htmlLog   = fopen(htmlLogName,   "w");
  if (errno)
    fprintf(stderr, "ERROR:  cannot open uid map file '%s': %s\n", htmlLogName, strerror(errno)), exit(1);

  sprintf(nameMapName,   "%s/readNames.txt", gkpStoreName);
  FILE    *nameMap   = fopen(nameMapName,   "w");
  if (errno)
    fprintf(stderr, "ERROR:  cannot open uid map file '%s': %s\n", nameMapName, strerror(errno)), exit(1);

  uint32  nERROR   = 0;  //  There aren't any errors, we just exit fatally if encountered.
  uint32  nWARNS   = 0;

  uint32  nLOADED  = 0;  //  Reads loaded
  uint64  bLOADED  = 0;  //  Bases loaded

  uint32  nSKIPPED = 0;
  uint64  bSKIPPED = 0;  //  Bases not loaded, too short

#if 0
  fprintf(htmlLog, "<!DOCTYPE html>\n");
  fprintf(htmlLog, "<html>\n");
  fprintf(htmlLog, "<head>\n");
  fprintf(htmlLog, "<title>gatekeeper load statistics</title>\n");
  fprintf(htmlLog, "<style type='text/css'>\n");
  fprintf(htmlLog, "body       { font-family: Helvetica, Verdana, sans-serif; }\n");
  fprintf(htmlLog, "h1, h2     { color: #ee3e80; }\n");
  fprintf(htmlLog, "p          { color: #665544; }\n");
  fprintf(htmlLog, "th, td     { border: 1px solid #111111; padding: 2px 2px 2px 2px; }\n");
  fprintf(htmlLog, "td:hover   { background-color: #e4e4e4; }\n");
  fprintf(htmlLog, "th:hover   { background-color: #d4d4d4; }\n");
  fprintf(htmlLog, "tr.details { visibility: collapse; }\n");
  fprintf(htmlLog, "</style>\n");
  fprintf(htmlLog, "</head>\n");
  fprintf(htmlLog, "<body>\n");
  fprintf(htmlLog, "<h2>Input Files</h2>\n");
  fprintf(htmlLog, "<table>\n");
#endif

  for (; firstFileArg < argc; firstFileArg++) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Starting file '%s'.\n", argv[firstFileArg]);

    compressedFileReader *inFile = new compressedFileReader(argv[firstFileArg]);
    char                 *line   = new char [10240];
    char                 *linekv = new char [10240];
    KeyAndValue           keyval;

    while (fgets(line, 10240, inFile->file()) != NULL) {
      chomp(line);
      strcpy(linekv, line);  //  keyval.find() modifies the input line, adding a nul byte to split the key and value.
      keyval.find(linekv);

      if (keyval.key() == NULL) {
        //  No key, so must be a comment or blank line
        continue;
      }

      if (strcasecmp(keyval.key(), "name") == 0) {
        gkpLibrary = gkpStore->gkStore_addEmptyLibrary(keyval.value());
        continue;
      }

      //  We'd better have a gkpLibrary defined, if not, the .gkp input file is incorrect.
      if (gkpLibrary == NULL) {
        fprintf(stderr, "WARNING: no 'name' tag in gkp input; creating library with name 'DEFAULT'.\n");
        gkpLibrary = gkpStore->gkStore_addEmptyLibrary(keyval.value());
        nWARNS++;
      }

      if        (strcasecmp(keyval.key(), "preset") == 0) {
        gkpLibrary->gkLibrary_parsePreset(keyval.value());

      } else if (strcasecmp(keyval.key(), "qv") == 0) {
        gkpLibrary->gkLibrary_setDefaultQV(keyval.value_double());

      } else if (strcasecmp(keyval.key(), "isNonRandom") == 0) {
        gkpLibrary->gkLibrary_setIsNonRandom(keyval.value_bool());

      } else if (strcasecmp(keyval.key(), "readType") == 0) {
        gkpLibrary->gkLibrary_setReadType(keyval.value());

      } else if (strcasecmp(keyval.key(), "removeDuplicateReads") == 0) {
        gkpLibrary->gkLibrary_setRemoveDuplicateReads(keyval.value_bool());

      } else if (strcasecmp(keyval.key(), "finalTrim") == 0) {
        gkpLibrary->gkLibrary_setFinalTrim(keyval.value());

      } else if (strcasecmp(keyval.key(), "removeSpurReads") == 0) {
        gkpLibrary->gkLibrary_setRemoveSpurReads(keyval.value_bool());

      } else if (strcasecmp(keyval.key(), "removeChimericReads") == 0) {
        gkpLibrary->gkLibrary_setRemoveChimericReads(keyval.value_bool());

      } else if (strcasecmp(keyval.key(), "checkForSubReads") == 0) {
        gkpLibrary->gkLibrary_setCheckForSubReads(keyval.value_bool());

      } else if (AS_UTL_fileExists(line, false, false)) {
        loadReads(gkpStore,
                  gkpLibrary,
                  gkpFileID++,
                  minReadLength,
                  nameMap,
                  htmlLog,
                  errorLog,
                  line,
                  nWARNS, nLOADED, bLOADED, nSKIPPED, bSKIPPED);

      } else {
        fprintf(stderr, "ERROR:  option '%s' not recognized, and not a file of reads.\n", line);
        exit(1);
      }
    }

    delete    inFile;
    delete [] line;
    delete [] linekv;
  }

#if 0
  fprintf(htmlLog, "</table>\n");
#endif

  gkpStore->gkStore_close();

  fclose(nameMap);
  fclose(errorLog);

  fprintf(stderr, "\n");
  fprintf(stderr, "Finished with:\n");
  fprintf(stderr, "  "F_U32" warnings (bad base or qv, too short, too long)\n", nWARNS);
  fprintf(stderr, "\n");
#if 0
  fprintf(stderr, "Read from inputs:\n");
  fprintf(stderr, "  "F_U64" bp.\n",    bLOADED);
  fprintf(stderr, "  "F_U32" reads.\n", nLOADED);
  fprintf(stderr, "\n");
#endif
  fprintf(stderr, "Loaded into store:\n");
  fprintf(stderr, "  "F_U64" bp.\n",    bLOADED);
  fprintf(stderr, "  "F_U32" reads.\n", nLOADED);
  fprintf(stderr, "\n");
  fprintf(stderr, "Skipped (too short):\n");
  fprintf(stderr, "  "F_U64" bp (%.4f%%).\n",    bSKIPPED, 100.0 * bSKIPPED / (bSKIPPED + bLOADED));
  fprintf(stderr, "  "F_U32" reads (%.4f%%).\n", nSKIPPED, 100.0 * nSKIPPED / (nSKIPPED + nLOADED));
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");

#if 0
  fprintf(htmlLog, "\n");
  fprintf(htmlLog, "<h2>Final Store</h2>\n");
  fprintf(htmlLog, "<table>\n");
  fprintf(htmlLog, "<tr><td colspan='2'>%s</td></tr>\n", gkpStoreName);
  fprintf(htmlLog, "<tr><td>readsLoaded</td><td>"F_U32" reads ("F_U64" bp)</td></tr>\n", nLOADED, bLOADED);
  fprintf(htmlLog, "<tr><td>readsSkipped</td><td>"F_U32" reads ("F_U64" bp) (read was too short)</td></tr>\n", nSKIPPED, bSKIPPED);
  fprintf(htmlLog, "<tr><td>warnings</td><td>"F_U32" warnings (invalid base or quality value)</td></tr>\n", nWARNS);
  fprintf(htmlLog, "</table>\n");
  fprintf(htmlLog, "\n");

  fprintf(htmlLog, "<script type='text/javascript'>\n");
  fprintf(htmlLog, "var toggleOne = function() {\n");
  fprintf(htmlLog, "  var table = this.closest('table');\n");
  fprintf(htmlLog, "  var elts  = table.querySelectorAll('.details');\n");
  fprintf(htmlLog, "\n");
  fprintf(htmlLog, "  for (var i=0; i<elts.length; i++) {\n");
  fprintf(htmlLog, "    if (!elts[i].enabled) {\n");
  fprintf(htmlLog, "      elts[i].enabled = true;\n");
  fprintf(htmlLog, "      elts[i].style.visibility = 'visible';\n");
  fprintf(htmlLog, "    } else {\n");
  fprintf(htmlLog, "      elts[i].enabled = false;\n");
  fprintf(htmlLog, "      elts[i].style.visibility = 'collapse';\n");
  fprintf(htmlLog, "    }\n");
  fprintf(htmlLog, "  }\n");
  fprintf(htmlLog, "}\n");
  fprintf(htmlLog, "\n");
  for (uint32 ii=0; ii<gkpFileID; ii++) {
    fprintf(htmlLog, "document.getElementById('gkpload%u').onclick = toggleOne;\n", ii);
    fprintf(htmlLog, "document.getElementById('gkpload%u').style   = 'cursor: pointer;';\n", ii);
  }
  fprintf(htmlLog, "</script>\n");
  fprintf(htmlLog, "\n");
  fprintf(htmlLog, "</body>\n");
  fprintf(htmlLog, "</html>\n");
#else
  fprintf(htmlLog, "sum "F_U32" "F_U64" "F_U32" "F_U64" "F_U32"\n", nLOADED, bLOADED, nSKIPPED, bSKIPPED, nWARNS);
#endif

  fclose(htmlLog);



  if (nERROR > 0)
    fprintf(stderr, "gatekeeperCreate did NOT finish successfully; too many errors.\n");

  if (bSKIPPED > 0.25 * (bSKIPPED + bLOADED))
    fprintf(stderr, "gatekeeperCreate did NOT finish successfully; too many bases skipped.  Check your reads.\n");

  if (nWARNS > 0.25 * (nLOADED))
    fprintf(stderr, "gatekeeperCreate did NOT finish successfully; too many warnings.  Check your reads.\n");

  if (nSKIPPED > 0.50 * (nLOADED))
    fprintf(stderr, "gatekeeperCreate did NOT finish successfully; too many short reads.  Check your reads!\n");

  if ((nERROR > 0) ||
      (bSKIPPED > 0.25 * (bSKIPPED + bLOADED)) ||
      (nWARNS   > 0.25 * (nSKIPPED + nLOADED)) ||
      (nSKIPPED > 0.50 * (nSKIPPED + nLOADED)))
    exit(1);

  fprintf(stderr, "gatekeeperCreate finished successfully.\n");

  exit(0);
}
