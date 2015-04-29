
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2012, J. Craig Venter Institute.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

const char *mainid = "$Id$";

#include "AS_global.H"
#include "gkStore.H"
#include "findKeyAndValue.H"
#include "AS_UTL_fileIO.H"



//  Support fastq of fasta, even in the same file.
//  Eventually want to support bax.h5 natively.


uint32  validSeq[256] = {0};


uint32
loadFASTA(gkStore              *gkpStore,
          gkLibrary            *gkpLibrary,
          char                 *L,
          char                 *H,
          char                 *S,
          uint32               &Slen,
          char                 *Q,
          compressedFileReader *F,
          FILE                 *errorLog,
          uint32               &nWARNS) {
  uint32  nLines = 0;
  bool    valid  = true;

  //  We've already read the header.  It's in L.  But we want to use L to load the sequence, so the
  //  header is copied to H.  We need to do a copy anyway - we need to return the next header in L.
  //  We can either copy L to H now, or use H to read lines and copy H to L at the end.  We're stuck
  //  either way.

  strcpy(H, L + 1);

  //  Load sequence.  This is a bit tricky, since we need to peek ahead
  //  and stop reading before the next header is loaded.


  S[0] = 0;
  Slen = 0;

  fgets(L, AS_MAX_READLEN, F->file());  nLines++;
  chomp(L);

  Q[0] = 0;
  
  while ((!feof(F->file())) && (L[0] != '>') && (valid == true)) {

    //  Copy in the sequence, as long as it is valid sequence.  If any invalid letters
    //  are found, stop copying, and stop reading sequence.

    bool  bogusBase = false;

    for (uint32 i=0; L[i]; i++) {
      switch (L[i]) {
        case 'a':   S[Slen] = 'A';  break;
        case 'c':   S[Slen] = 'C';  break;
        case 'g':   S[Slen] = 'G';  break;
        case 't':   S[Slen] = 'T';  break;
        case 'A':   S[Slen] = 'A';  break;
        case 'C':   S[Slen] = 'C';  break;
        case 'G':   S[Slen] = 'G';  break;
        case 'T':   S[Slen] = 'T';  break;
        case 'n':   S[Slen] = 'N';  break;
        case 'N':   S[Slen] = 'N';  break;
        default:
          fprintf(errorLog, "read '%s' has invalid base '%c' (0x%02x) at position %u.  Converted to 'N'.\n",
                 H, L[i], L[i], i);
          S[Slen]   = 'N';
          valid     = false;
          bogusBase = true;
          break;
      }

      Slen++;
    }

    if (bogusBase)
      nWARNS++;

    //  If we're still valid, grab the next line; it should be a header.  If not valid, pass the
    //  last line back so we can fail the header check on the next iteration.

    if (valid) {
      fgets(L, AS_MAX_READLEN, F->file());  nLines++;
      chomp(L);
    }
  }

  //  Do NOT clear L, it contains the next header.

  return(nLines);
}



uint32          
loadFASTQ(gkStore              *gkpStore,
          gkLibrary            *gkpLibrary,
          char                 *L,
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

  fgets(S, AS_MAX_READLEN, F->file());
  chomp(S);

  //  Check for and correct invalid bases.

  uint32 baseErrors = 0;

  for (uint32 i=0; S[i]; i++) {
    switch (S[i]) {
      case 'a':   S[i] = 'A';  break;
      case 'c':   S[i] = 'C';  break;
      case 'g':   S[i] = 'G';  break;
      case 't':   S[i] = 'T';  break;
      case 'A':                break;
      case 'C':                break;
      case 'G':                break;
      case 'T':                break;
      case 'n':   S[i] = 'N';  break;
      case 'N':                break;
      default:
        S[i]       = 'N';
        Q[i]       = '!';
        baseErrors++;
        break;
    }

    Slen++;
  }

  if (baseErrors > 0) {
        fprintf(errorLog, "read '%s' has %u invalid base%s.  Converted to 'N'.\n",
                L, baseErrors, (baseErrors > 1) ? "s" : "");
    nWARNS++;
  }

  //  Load the qv header, and then load the qvs themselves over the header.
  Q[0] = 0;
  fgets(Q, AS_MAX_READLEN, F->file());
  fgets(Q, AS_MAX_READLEN, F->file());
  chomp(Q);

  //  Convert from the (assumed to be) Sanger QVs to the CA offset '0' QVs.

  uint32 QVerrors = 0;

  for (uint32 i=0; Q[i]; i++) {
    if (Q[i] < '!') {
      Q[i] = '!';
      QVerrors++;
    }

    if (Q[i] > 'Z') {
      Q[i] = 'Z';
      QVerrors++;
    }

    Q[i] -= '!';
    Q[i] += '0';
  }

  if (QVerrors > 0) {
    fprintf(errorLog, "read '%s' has invalid %u QV%s.  Converted to min or max value.\n",
            L, QVerrors, (QVerrors > 1) ? "s" : "");
    nWARNS++;
  }

  //  Clear the lines, so we can load the next one.

  L[0] = 0;

  return(4);  //  FASTQ always reads exactly four lines
}





void
loadReads(gkStore    *gkpStore,
          gkLibrary  *gkpLibrary,
          uint32      minReadLength,
          FILE       *nameMap,
          FILE       *errorLog,
          char       *fileName,
          uint32     &nFASTA,
          uint32     &nFASTQ,
          uint32     &nWARNS,
          uint32     &nSHORT) {
  char    *L = new char [AS_MAX_READLEN + 1];
  char    *H = new char [AS_MAX_READLEN + 1];

  uint32   Slen = 0;
  char    *S = new char [AS_MAX_READLEN + 1];

  uint32   Qlen = 0;
  char    *Q = new char [AS_MAX_READLEN + 1];

  uint64   lineNumber = 1;

  fprintf(stderr, "\n");
  fprintf(stderr, "  Loading reads from '%s'\n", fileName);

  compressedFileReader *F = new compressedFileReader(fileName);

  uint32   nFASTAlocal = 0;
  uint32   nFASTQlocal = 0;
  uint32   nWARNSlocal = 0;
  uint32   nSHORTlocal = 0;

  fgets(L, AS_MAX_READLEN, F->file());
  chomp(L);

  while (!feof(F->file())) {

    if      (L[0] == '>') {
      lineNumber += loadFASTA(gkpStore, gkpLibrary, L, H, S, Slen, Q, F, errorLog, nWARNS);
      nFASTAlocal++;
    }

    else if (L[0] == '@') {
      lineNumber += loadFASTQ(gkpStore, gkpLibrary, L, H, S, Slen, Q, F, errorLog, nWARNS);
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
      fprintf(errorLog, "read '%s' of length %u in file '%s' at line "F_U64" is too short, skipping.\n",
              H, Slen, fileName, lineNumber);
      nSHORTlocal++;
      S[0] = 0;
      Q[0] = 0;
    }

    if (S[0] != 0) {
      gkRead     *nr = gkpStore->gkStore_addEmptyRead(gkpLibrary);
      gkReadData *nd = (Q[0] == 0) ? nr->gkRead_encodeSeqQlt(H, S, gkpLibrary->gkLibrary_defaultQV()) :
                                     nr->gkRead_encodeSeqQlt(H, S, Q);

      gkpStore->gkStore_stashReadData(nr, nd);

      delete nd;

      fprintf(nameMap, F_U32"\t%s\n", gkpStore->gkStore_getNumReads(), H);
    }

    //  If L[0] is nul, we need to load the next line.  If not, the next line is the header (from
    //  the fasta loader).

    if (L[0] == 0) {
      fgets(L, AS_MAX_READLEN, F->file());  lineNumber++;
      chomp(L);
    }
  }

  delete    F;

  delete [] Q;
  delete [] S;
  delete [] H;
  delete [] L;

  lineNumber--;  //  The last fgets() returns EOF, but we still count the line.

  fprintf(stderr, "    Processed "F_U32" lines.\n", lineNumber);

  if (nFASTAlocal > 0)
    fprintf(stderr, "    Loaded "F_U32" FASTA format reads.\n", nFASTAlocal);
  if (nFASTQlocal > 0)
    fprintf(stderr, "    Loaded "F_U32" FASTQ format reads.\n", nFASTQlocal);
  if (nWARNSlocal > 0)
    fprintf(stderr, "      WARNING: "F_U32" reads issued a warning.\n", nWARNSlocal);
  if (nSHORTlocal > 0)
    fprintf(stderr, "      WARNING: "F_U32" reads (%0.4f%%) were too short (< %ubp) and were ignored.\n",
            nSHORTlocal, 100.0 * nSHORTlocal / (nFASTAlocal + nFASTQlocal), minReadLength);

  nFASTA += nFASTAlocal;
  nFASTQ += nFASTQlocal;
  nWARNS += nWARNSlocal;
  nSHORT += nSHORTlocal;
};




int
main(int argc, char **argv) {
  char            *gkpStoreName      = NULL;
  char            *outPrefix         = NULL;

  uint32           minReadLength     = 0;

  uint32           firstFileArg      = 0;

  char             errorLogName[FILENAME_MAX];
  char             nameMapName[FILENAME_MAX];


  //argc = AS_configure(argc, argv);

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
      fprintf(stderr, "ERROR: no gkpStore (-g) supplied.\n");
    if (firstFileArg == 0)
      fprintf(stderr, "ERROR: no input files supplied.\n");

    exit(1);
  }


  gkStore     *gkpStore     = new gkStore(gkpStoreName, gkStore_extend);
  gkRead      *gkpRead      = NULL;
  gkReadData   gkpReadData;
  gkLibrary   *gkpLibrary   = NULL;

  uint32       inLineLen    = 1024;
  char         inLine[1024] = { 0 };

  validSeq['a'] = validSeq['c'] = validSeq['g'] = validSeq['t'] = validSeq['n'] = 1;
  validSeq['A'] = validSeq['C'] = validSeq['G'] = validSeq['T'] = validSeq['N'] = 1;

  errno = 0;

  sprintf(errorLogName, "%s.errorLog",    gkpStoreName);
  FILE    *errorLog = fopen(errorLogName, "w");
  if (errno)
    fprintf(stderr, "ERROR:  cannot open error file '%s': %s\n", errorLogName, strerror(errno)), exit(1);

  sprintf(nameMapName,   "%s/readNames.txt", gkpStoreName);
  FILE    *nameMap   = fopen(nameMapName,   "w");
  if (errno)
    fprintf(stderr, "ERROR:  cannot open uid map file '%s': %s\n", nameMapName, strerror(errno)), exit(1);

  uint32  nERROR = 0;  //  There aren't any errors, we just exit fatally if encountered.
  uint32  nWARNS = 0;
  uint32  nSHORT = 0;
  uint32  nFASTA = 0;
  uint32  nFASTQ = 0;

  for (; firstFileArg < argc; firstFileArg++) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Starting file '%s'.\n", argv[firstFileArg]);

    compressedFileReader *inFile = new compressedFileReader(argv[firstFileArg]);
    char                 *line   = new char [10240];
    KeyAndValue           keyval;

    while (fgets(line, 10240, inFile->file()) != NULL) {
      chomp(line);
      keyval.find(line);

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

      } else if (strcasecmp(keyval.key(), "trustHomopolymerRuns") == 0) {
        gkpLibrary->gkLibrary_setTrustHomopolymerRuns(keyval.value_bool());

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

      } else if (AS_UTL_fileExists(keyval.key(), false, false)) {
        loadReads(gkpStore,
                  gkpLibrary,
                  minReadLength,
                  nameMap,
                  errorLog,
                  keyval.key(),
                  nFASTA, nFASTQ, nWARNS, nSHORT);

      } else {
        fprintf(stderr, "ERROR:  option '%s' not recognized, and not a file of reads.\n", line);
        exit(1);
      }
    }

    delete    inFile;
    delete [] line;
  }

  delete gkpStore;

  fclose(nameMap);
  fclose(errorLog);

  fprintf(stderr, "\n");
  fprintf(stderr, "Finished with:\n");
  fprintf(stderr, "  %u errors\n", nERROR);
  fprintf(stderr, "  %u warnings.\n", nERROR, nWARNS);
  fprintf(stderr, "\n");
  fprintf(stderr, "Skipped:\n");
  fprintf(stderr, "  %u short reads (%.4f%%).\n", nSHORT, 100.0 * nSHORT / (nSHORT + nFASTA + nFASTQ));
  fprintf(stderr, "\n");
  fprintf(stderr, "Loaded:\n");
  fprintf(stderr, "  %u FASTA reads.\n", nFASTA);
  fprintf(stderr, "  %u FASTQ reads.\n", nFASTQ);
  fprintf(stderr, "\n");

  if (nERROR > 0)
    fprintf(stderr, "gatekeeperCreate did NOT finish successfully; too many errors.\n");

  if (nWARNS > 0.10 * (nFASTA + nFASTQ))
    fprintf(stderr, "gatekeeperCreate did NOT finish successfully; too many warnings.  Check your reads.\n");

  if (nSHORT > 0.25 * (nFASTA + nFASTQ))
    fprintf(stderr, "gatekeeperCreate did NOT finish successfully; too many short reads.  Check your reads!\n");

  if ((nERROR > 0) ||
      (nWARNS > 0.10 * (nSHORT + nFASTA + nFASTQ)) ||
      (nSHORT > 0.25 * (nSHORT + nFASTA + nFASTQ)))
    exit(1);

  fprintf(stderr, "gatekeeperCreate finished successfully.\n");
  exit(0);
}
