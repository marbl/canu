
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
 *    src/stores/gatekeeperCreate.C
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
 *    Sergey Koren beginning on 2019-JUL-17
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "sqStore.H"
#include "files.H"
#include "strings.H"

#include "mt19937ar.H"

#include <algorithm>

#undef  UPCASE  //  Don't convert lowercase to uppercase, special case for testing alignments.
#define UPCASE  //  Convert lowercase to uppercase.  Probably needed.



//  Canu doesn't really use QVs as of late 2015.  Only consensus is aware of them, but since read
//  correction doesn't output QVs, the QVs that consensus sees are totally bogus.  So, disable
//  storage of them - seqStoreCreate will store a default QV for every read, and consensus will
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
          int32                &Slen,
          uint8                *Q,
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
  Q[0] = 255;  //  Sentinel to tell sqStore to use the fixed QV value

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
        case 'u':   S[Slen] = 'T';  break;
#else
        case 'a':   S[Slen] = 'a';  break;
        case 'c':   S[Slen] = 'c';  break;
        case 'g':   S[Slen] = 'g';  break;
        case 't':   S[Slen] = 't';  break;
        case 'u':   S[Slen] = 't';  break;
#endif
        case 'A':   S[Slen] = 'A';  break;
        case 'C':   S[Slen] = 'C';  break;
        case 'G':   S[Slen] = 'G';  break;
        case 'T':   S[Slen] = 'T';  break;
        case 'U':   S[Slen] = 'T';  break;
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
    fprintf(errorLog, "read '%s' has " F_U32 " invalid base%s.  Converted to 'N'.\n",
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
          int32                &Slen,
          uint8                *Q,
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

  L[AS_MAX_READLEN+1-2] = 0;  //  This too.
  L[AS_MAX_READLEN+1-1] = 0;

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
      case 'u':   S[i] = 'T';  break;
#else
      case 'a':                break;
      case 'c':                break;
      case 'g':                break;
      case 't':                break;
      case 'u':                break;
#endif
      case 'A':                break;
      case 'C':                break;
      case 'G':                break;
      case 'T':                break;
      case 'U':   S[i] = 'T';  break;
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
        fprintf(errorLog, "read '%s' has " F_U32 " invalid base%s.  Converted to 'N'.\n",
                L, baseErrors, (baseErrors > 1) ? "s" : "");
    nWARNS++;
  }

  //  Load the qv header, and then load the qvs themselves over the header.

  L[0] = 0;
  fgets(L, AS_MAX_READLEN+1, F->file());   //  NOTE: Using L, not Q, since L is char and Q is uint8
  fgets(L, AS_MAX_READLEN+1, F->file());
  chomp(L);

  //  As with the base, we need to suck in the rest of the longer-than-allowed QV string.  But we don't need to report it
  //  or do anything fancy, just advance the file pointer.

  if ((L[AS_MAX_READLEN+1-2] != 0) && (L[AS_MAX_READLEN+1-2] != '\n')) {
    char    *overflow = new char [1048576];

    do {
      overflow[1048576-2] = 0;
      overflow[1048576-1] = 0;
      fgets(overflow, 1048576, F->file());
    } while (overflow[1048576-2] != 0);

    delete [] overflow;
  }

  //  If we're not using QVs, just terminate the sequence.

  Q[0] = 255;  //  Sentinel to tell sqStore to use the fixed QV value

  //  But if we are storing QVs, check lengths and convert from letters to integers

#ifndef DO_NOT_STORE_QVs
  int32    sLen = strlen(S);
  int32    qLen = strlen(L);

  if (sLen < qLen) {
    fprintf(errorLog, "read '%s' sequence length %u quality length %u; quality values trimmed.\n",
            H, sLen, qLen);
    nWARNS++;
    L[sLen] = 0;
  }

  if (sLen > qLen) {
    fprintf(errorLog, "read '%s' sequence length %u quality length %u; sequence trimmed.\n",
            H, sLen, qLen);
    nWARNS++;
    S[qLen] = 0;
  }

  uint32 QVerrors = 0;

  for (uint32 i=0; L[i]; i++) {
    if (L[i] < '!') {  //  QV=0, ASCII=33
      L[i] = '!';
      QVerrors++;
    }

    if (L[i] > '!' + 60) {  //  QV=60, ASCII=93=']'
      L[i] = '!' + 60;
      QVerrors++;
    }

    Q[i] = L[i] - '!';
  }

  if (QVerrors > 0) {
    fprintf(errorLog, "read '%s' has " F_U32 " invalid QV%s.  Converted to min or max value.\n",
            L, QVerrors, (QVerrors > 1) ? "s" : "");
    nWARNS++;
  }
#endif

  //  Clear the lines, so we can load the next one.

  L[0] = 0;

  return(4);  //  FASTQ always reads exactly four lines
}





void
loadReads(sqStore    *seqStore,
          sqLibrary  *seqLibrary,
          uint32      seqFileID,
          uint32      minReadLength,
          FILE       *nameMap,
          FILE       *loadLog,
          FILE       *errorLog,
          char       *fileName,
          uint32     &nWARNS,
          uint32     &nLOADED,
          uint64     &bLOADED,
          uint32     &nSKIPPED,
          uint64     &bSKIPPED) {
  char    *L = new char  [AS_MAX_READLEN + 1];  //  +1.  One for the newline, and one for the terminating nul.
  char    *H = new char  [AS_MAX_READLEN + 1];
  char    *S = new char  [AS_MAX_READLEN + 1];
  uint8   *Q = new uint8 [AS_MAX_READLEN + 1];

  int32    Sbgn = 0;
  int32    Send = 0;
  int32    Slen = 0;

  uint64   lineNumber = 1;

  fprintf(stderr, "\n");
  fprintf(stderr, "  Loading reads from '%s'\n", fileName);

  fprintf(loadLog, "nam " F_U32 " %s\n", seqFileID, fileName);

  fprintf(loadLog, "lib preset=N/A");
  fprintf(loadLog,    " defaultQV=%u",            seqLibrary->sqLibrary_defaultQV());
  fprintf(loadLog,    " isNonRandom=%s",          seqLibrary->sqLibrary_isNonRandom()          ? "true" : "false");
  fprintf(loadLog,    " removeDuplicateReads=%s", seqLibrary->sqLibrary_removeDuplicateReads() ? "true" : "false");
  fprintf(loadLog,    " finalTrim=%s",            seqLibrary->sqLibrary_finalTrim()            ? "true" : "false");
  fprintf(loadLog,    " removeSpurReads=%s",      seqLibrary->sqLibrary_removeSpurReads()      ? "true" : "false");
  fprintf(loadLog,    " removeChimericReads=%s",  seqLibrary->sqLibrary_removeChimericReads()  ? "true" : "false");
  fprintf(loadLog,    " checkForSubReads=%s\n",   seqLibrary->sqLibrary_checkForSubReads()     ? "true" : "false");

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
      fprintf(errorLog, "invalid read header '%.40s%s' in file '%s' at line " F_U64 ", skipping.\n",
              L, (strlen(L) > 80) ? "..." : "", fileName, lineNumber);
      L[0] = 0;
      nWARNSlocal++;
    }

    //  Trim N from the ends.

    Sbgn = 0;
    Send = Slen - 1;

    while ((Sbgn <= Send) && ((S[Sbgn] == 'N') ||
                              (S[Sbgn] == 'n'))) {
      S[Sbgn] = 0;
      Q[Sbgn] = 0;
      Sbgn++;
    }

    while ((Sbgn <= Send) && ((S[Send] == 'N') ||
                              (S[Send] == 'n'))) {
      S[Send] = 0;
      Q[Send] = 0;
      Send--;
    }

    Send++;

    if ((Sbgn > 0) && (Send < Slen))
      fprintf(errorLog, "read '%s' of length " F_U32 " in file '%s' at line " F_U64 " - trimmed " F_S32 " non-ACGT bases from the 5' and " F_S32 " non-ACGT bases from the 3' end.\n",
              H, Slen, fileName, lineNumber, Sbgn, Slen - Send);

    else if (Sbgn > 0)
      fprintf(errorLog, "read '%s' of length " F_U32 " in file '%s' at line " F_U64 " - trimmed " F_S32 " non-ACGT bases from the 5' end.\n",
              H, Slen, fileName, lineNumber, Sbgn);

    else if (Send < Slen)
      fprintf(errorLog, "read '%s' of length " F_U32 " in file '%s' at line " F_U64 " - trimmed " F_S32 " non-ACGT bases from the 3' end.\n",
              H, Slen, fileName, lineNumber, Slen - Send);

    Slen = Send - Sbgn;

    //  Drop short reads.  "Rick Wakeman, eat your heart out. Here we go!"

    if (Slen < minReadLength) {
      fprintf(errorLog, "read '%s' of length " F_U32 " in file '%s' at line " F_U64 " - too short, skipping.\n",
              H, Slen, fileName, lineNumber);

      if (isFASTA) {
        nSKIPPEDAlocal += 1;
        bSKIPPEDAlocal += Slen;
      }

      if (isFASTQ) {
        nSKIPPEDQlocal += 1;
        bSKIPPEDQlocal += Slen;
      }

      Slen = 0;

      S[0] = 0;
      Q[0] = 0;
    }

    //  Otherwise, load it!

    else {
      sqReadData *readData = seqStore->sqStore_addEmptyRead(seqLibrary);

      readData->sqReadData_setName(H);
      readData->sqReadData_setBasesQuals(S + Sbgn, Q + Sbgn);

      seqStore->sqStore_stashReadData(readData);

      delete readData;

      if (isFASTA) {
        nLOADEDAlocal += 1;
        bLOADEDAlocal += Slen;
      }

      if (isFASTQ) {
        nLOADEDQlocal += 1;
        bLOADEDQlocal += Slen;
      }

      fprintf(nameMap, F_U32"\t%s\n", seqStore->sqStore_getNumReads(), H);
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

  fprintf(stderr, "    Processed " F_U64 " lines.\n", lineNumber);

  fprintf(stderr, "    Loaded " F_U64 " bp from:\n", bLOADEDAlocal + bLOADEDQlocal);
  if (nFASTAlocal > 0)
    fprintf(stderr, "      " F_U32 " FASTA format reads (" F_U64 " bp).\n", nFASTAlocal, bLOADEDAlocal);
  if (nFASTQlocal > 0)
    fprintf(stderr, "      " F_U32 " FASTQ format reads (" F_U64 " bp).\n", nFASTQlocal, bLOADEDQlocal);

  if (nWARNSlocal > 0)
    fprintf(stderr, "    WARNING: " F_U32 " reads issued a warning.\n", nWARNSlocal);

  if (nSKIPPEDAlocal > 0)
    fprintf(stderr, "    WARNING: " F_U32 " reads (%0.4f%%) with " F_U64 " bp (%0.4f%%) were too short (< " F_U32 "bp) and were ignored.\n",
            nSKIPPEDAlocal, 100.0 * nSKIPPEDAlocal / (nSKIPPEDAlocal + nLOADEDAlocal),
            bSKIPPEDAlocal, 100.0 * bSKIPPEDAlocal / (bSKIPPEDAlocal + bLOADEDAlocal),
            minReadLength);

  if (nSKIPPEDQlocal > 0)
    fprintf(stderr, "    WARNING: " F_U32 " reads (%0.4f%%) with " F_U64 " bp (%0.4f%%) were too short (< " F_U32 "bp) and were ignored.\n",
            nSKIPPEDQlocal, 100.0 * nSKIPPEDQlocal / (nSKIPPEDQlocal + nLOADEDQlocal),
            bSKIPPEDQlocal, 100.0 * bSKIPPEDQlocal / (bSKIPPEDQlocal + bLOADEDQlocal),
            minReadLength);

  //  Write status to HTML

  fprintf(loadLog, "dat " F_U32 " " F_U64 " " F_U32 " " F_U64 " " F_U32 " " F_U64 " " F_U32 " " F_U64 " " F_U32 "\n",
          nLOADEDAlocal, bLOADEDAlocal,
          nSKIPPEDAlocal, bSKIPPEDAlocal,
          nLOADEDQlocal, bLOADEDQlocal,
          nSKIPPEDQlocal, bSKIPPEDQlocal,
          nWARNSlocal);

  //  Add the just loaded numbers to the global numbers

  nWARNS   += nWARNSlocal;

  nLOADED  += nLOADEDAlocal + nLOADEDQlocal;
  bLOADED  += bLOADEDAlocal + bLOADEDQlocal;

  nSKIPPED += nSKIPPEDAlocal + nSKIPPEDQlocal;
  bSKIPPED += bSKIPPEDAlocal + bSKIPPEDQlocal;
};





bool
createStore(const char *seqStoreName,
            uint32      firstFileArg,
            char      **argv,
            uint32      argc,
            uint32      minReadLength) {

  sqStore     *seqStore     = sqStore::sqStore_open(seqStoreName, sqStore_create);   //  sqStore_extend MIGHT work
  sqRead      *seqRead      = NULL;
  sqLibrary   *seqLibrary   = NULL;
  uint32       seqFileID    = 0;      //  Used for HTML output, an ID for each file loaded.

  uint32       inLineLen    = 1024;
  char         inLine[1024] = { 0 };

  FILE        *errorLog = AS_UTL_openOutputFile(seqStoreName, '/', "errorLog");
  FILE        *loadLog  = AS_UTL_openOutputFile(seqStoreName, '/', "load.dat");
  FILE        *nameMap  = AS_UTL_openOutputFile(seqStoreName, '/', "readNames.txt");

  uint32       nWARNS   = 0;

  uint32       nLOADED  = 0;  //  Reads loaded
  uint64       bLOADED  = 0;  //  Bases loaded

  uint32       nSKIPPED = 0;
  uint64       bSKIPPED = 0;  //  Bases not loaded, too short


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
        seqLibrary = seqStore->sqStore_addEmptyLibrary(keyval.value());
        continue;
      }

      //  We'd better have a seqLibrary defined, if not, the .seq input file is incorrect.
      if (seqLibrary == NULL) {
        fprintf(stderr, "WARNING: no 'name' tag in seq input; creating library with name 'DEFAULT'.\n");
        seqLibrary = seqStore->sqStore_addEmptyLibrary(keyval.value());
        nWARNS++;
      }

      if        (strcasecmp(keyval.key(), "preset") == 0) {
        seqLibrary->sqLibrary_parsePreset(keyval.value());

      } else if (strcasecmp(keyval.key(), "qv") == 0) {
        seqLibrary->sqLibrary_setDefaultQV(keyval.value_double());

      } else if (strcasecmp(keyval.key(), "isNonRandom") == 0) {
        seqLibrary->sqLibrary_setIsNonRandom(keyval.value_bool());

      } else if (strcasecmp(keyval.key(), "readType") == 0) {
        seqLibrary->sqLibrary_setReadType(keyval.value());

      } else if (strcasecmp(keyval.key(), "removeDuplicateReads") == 0) {
        seqLibrary->sqLibrary_setRemoveDuplicateReads(keyval.value_bool());

      } else if (strcasecmp(keyval.key(), "finalTrim") == 0) {
        seqLibrary->sqLibrary_setFinalTrim(keyval.value());

      } else if (strcasecmp(keyval.key(), "removeSpurReads") == 0) {
        seqLibrary->sqLibrary_setRemoveSpurReads(keyval.value_bool());

      } else if (strcasecmp(keyval.key(), "removeChimericReads") == 0) {
        seqLibrary->sqLibrary_setRemoveChimericReads(keyval.value_bool());

      } else if (strcasecmp(keyval.key(), "checkForSubReads") == 0) {
        seqLibrary->sqLibrary_setCheckForSubReads(keyval.value_bool());

      } else if (fileExists(line)) {
        loadReads(seqStore,
                  seqLibrary,
                  seqFileID++,
                  minReadLength,
                  nameMap,
                  loadLog,
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

  fprintf(loadLog, "sum " F_U32 " " F_U64 " " F_U32 " " F_U64 " " F_U32 "\n", nLOADED, bLOADED, nSKIPPED, bSKIPPED, nWARNS);

  seqStore->sqStore_close();

  AS_UTL_closeFile(nameMap,  seqStoreName, '/', "readNames.txt");
  AS_UTL_closeFile(loadLog,  seqStoreName, '/', "load.dat");
  AS_UTL_closeFile(errorLog, seqStoreName, '/', "errorLog");

  fprintf(stderr, "\n");
  fprintf(stderr, "Finished with:\n");
  fprintf(stderr, "  " F_U32 " warnings (bad base or qv, too short, too long)\n", nWARNS);
  fprintf(stderr, "\n");
  fprintf(stderr, "Loaded into store:\n");
  fprintf(stderr, "  " F_U64 " bp.\n",    bLOADED);
  fprintf(stderr, "  " F_U32 " reads.\n", nLOADED);
  fprintf(stderr, "\n");
  fprintf(stderr, "Skipped (too short):\n");
  fprintf(stderr, "  " F_U64 " bp (%.4f%%).\n",    bSKIPPED, (bSKIPPED + bLOADED > 0) ? (100.0 * bSKIPPED / (bSKIPPED + bLOADED)) : 0);
  fprintf(stderr, "  " F_U32 " reads (%.4f%%).\n", nSKIPPED, (nSKIPPED + nLOADED > 0) ? (100.0 * nSKIPPED / (nSKIPPED + nLOADED)) : 0);
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");

  if (bSKIPPED > 0.25 * (bSKIPPED + bLOADED))
    fprintf(stderr, "sqStoreCreate did NOT finish successfully; too many bases skipped.  Check your reads.\n");

  if (nWARNS > 0.25 * (nLOADED))
    fprintf(stderr, "sqStoreCreate did NOT finish successfully; too many warnings.  Check your reads.\n");

  // do not output an error for too many short reads

  if ((bSKIPPED > 0.25 * (bSKIPPED + bLOADED)) ||
      (nWARNS   > 0.25 * (nSKIPPED + nLOADED)) ||
      (nSKIPPED > 0.50 * (nSKIPPED + nLOADED)))
    exit(0);

  return(true);
}



struct rl_t {
  uint32   readID;
  uint32   length;
  double   score;
};

bool  byScore(const rl_t &a, const rl_t &b)   { return(a.score > b.score); }


bool
deleteShortReads(const char *seqStoreName,
                 uint64      genomeSize,
                 double      desiredCoverage,
                 double      lengthBias) {
  mtRandom     mtctx;
  uint64       desiredBases = (uint64)floor(genomeSize * desiredCoverage);

  if (desiredBases == 0)
    return(true);

  //
  //  Open the store for modification.
  //

  sqStore     *seqStore     = sqStore::sqStore_open(seqStoreName, sqStore_extend);

  uint32       nReads       = seqStore->sqStore_getNumReads();
  rl_t        *readLen      = new rl_t [nReads + 1];
  uint64       readLenSum   = 0;

  //
  //  Initialize our list of read scores, and report a summary.
  //

  fprintf(stderr, "Analyzing read lengths in store '%s'\n", seqStoreName);
  fprintf(stderr, "\n");
  fprintf(stderr, "For genome size of " F_U64 " bases, at coverage %.2f, want to keep " F_U64 " bases.\n",
          genomeSize, desiredCoverage, desiredBases);
  fprintf(stderr, "\n");

  for (uint32 ii=0; ii<nReads+1; ii++) {
    uint32  len = seqStore->sqStore_getRead(ii)->sqRead_sequenceLength();

    readLen[ii].readID = ii;
    readLen[ii].length = len;
    readLen[ii].score  = mtctx.mtRandomRealOpen() * pow(len, lengthBias);

    readLenSum  += len;
  }

  fprintf(stderr, "Found %12" F_U32P " reads of length\n", nReads);
  fprintf(stderr, "      %12" F_U64P " bases.\n", readLenSum);
  fprintf(stderr, "\n");

  //
  //  Sort the array by length, keep only the longest reads up to coverage * genomeSize.
  //

#ifdef _GLIBCXX_PARALLEL
  __gnu_sequential::
#endif
  sort(readLen, readLen + nReads+1, byScore);

  readLenSum = 0;

  uint32   readsKept = 0;
  uint32   readsLost = 0;

  for (uint32 ii=0; ii<nReads+1; ii++) {
    fprintf(stdout, "%7u %8u %12.4f%s\n",
            readLen[ii].readID,
            readLen[ii].length,
            readLen[ii].score,
            (readLenSum < desiredBases) ? "" : " REMOVED");

    if (readLenSum < desiredBases) {
      readLenSum  += readLen[ii].length;

      readsKept++;
    }

    else {
      seqStore->sqStore_setIgnore(readLen[ii].readID);

      readsLost++;
    }
  }

  delete [] readLen;

  seqStore->sqStore_close();

  fprintf(stderr, "Kept  %12" F_U32P " reads of length\n", readsKept);
  fprintf(stderr, "      %12" F_U64P " bases.\n", readLenSum);
  fprintf(stderr, "\n");

  return(true);
}





int
main(int argc, char **argv) {
  char            *seqStoreName      = NULL;

  uint32           minReadLength     = 0;
  uint64           genomeSize        = 0;
  double           desiredCoverage   = 0;
  double           lengthBias        = 1.0;

  uint32           firstFileArg      = 0;

  //  Initialize the global.

  validSeq['a'] = validSeq['c'] = validSeq['g'] = validSeq['t'] = validSeq['n'] = 1;
  validSeq['A'] = validSeq['C'] = validSeq['G'] = validSeq['T'] = validSeq['N'] = 1;

  //  Parse options.

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      seqStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-minlength") == 0) {
      minReadLength = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-genomesize") == 0) {
      genomeSize = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-coverage") == 0) {
      desiredCoverage = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-bias") == 0) {
      lengthBias = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "--") == 0) {
      firstFileArg = arg++;
      break;

    } else if (argv[arg][0] == '-') {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);

    } else {
      firstFileArg = arg;
      break;
    }
    arg++;
  }

  if (seqStoreName == NULL)
    err.push_back("ERROR: no seqStore (-o) supplied.\n");

  if (firstFileArg == 0)
    err.push_back("ERROR: no input files supplied.\n");

  if ((desiredCoverage > 0) && (genomeSize == 0))
    err.push_back("ERROR: no genome size (-genomesize) set, needed for coverage filtering (-coverage) to work.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -o seqStore [-minlength L] [-genomesize G -coverage C] input.ssi\n", argv[0]);
    fprintf(stderr, "  -o seqStore            load raw reads into new seqStore\n");
    fprintf(stderr, "  \n");
    fprintf(stderr, "  -minlength L           discard reads shorter than L\n");
    fprintf(stderr, "  \n");
    fprintf(stderr, "  -genomesize G          expected genome size, for keeping only the longest reads\n");
    fprintf(stderr, "  -coverage C            desired coverage in long reads\n");
    fprintf(stderr, "  \n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }


  if (createStore(seqStoreName, firstFileArg, argv, argc, minReadLength) &&
      deleteShortReads(seqStoreName, genomeSize, desiredCoverage, lengthBias)) {
    fprintf(stderr, "sqStoreCreate finished successfully.\n");
    exit(0);
  } else {
    fprintf(stderr, "sqStoreCreate terminated abnormally.\n");
    exit(1);
  }
}
