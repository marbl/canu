
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




uint32
loadFASTA(gkStore *gkpStore,
          gkLibrary *gkpLibrary,
          char *L,
          char *H,
          char *S,
          compressedFileReader *F) {

  //  We've already read the header.  It's in L.  But we want to use L to load the sequence, so the
  //  header is copied to H.  We need to do a copy anyway - we need to return the next header in L.
  //  We can either copy L to H now, or use H to read lines and copy H to L at the end.  We're stuck
  //  either way.

  strcpy(H, L + 1);

  //  Load sequence.  This is a bit tricky, since we need to peek ahead
  //  and stop reading before the next header is loaded.

  //  Load sequence.

  uint32  Slen = 0;

  fgets(L, AS_MAX_READLEN, F->file());
  chomp(L);
  
  while ((!feof(F->file())) && (L[0] != '>')) {
    strcat(S + Slen, L);
    Slen += strlen(L);

    fgets(L, AS_MAX_READLEN, F->file());
    chomp(L);
  }

  //  Add a new read to the store.

  gkRead     *nr = gkpStore->gkStore_addEmptyRead(gkpLibrary);

  //  Populate the read with the appropriate gunk, based on the library type.

  //assert(gkpLibrary->gkLibrary_encoding() == GKREAD_TYPE_SEQ_QLT);

  gkReadData *nd = nr->gkRead_encodeSeqQlt(H, S, gkpLibrary->gkLibrary_defaultQV());

  //  Stash the encoded data in the store

  gkpStore->gkStore_stashReadData(nr, nd);

  //  And then toss it out.

  delete nd;

  //  Do NOT clear L, it contains the next header.

  //H[0] = 0;
  //S[0] = 0;

  return(0);
}



uint32
loadFASTQ(gkStore *gkpStore,
          gkLibrary *gkpLibrary,
          char *L,
          char *H,
          char *S,
          char *Q,
          compressedFileReader *F) {

  //  We've already read the header.  It's in L.

  strcpy(H, L + 1);

  //  Load sequence.
  fgets(S, AS_MAX_READLEN, F->file());
  chomp(S);

  //  Check for and correct invalid bases.
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
        fprintf(stderr, "-- WARNING:  read '%s' has invalid base '%c' (0x%02x) at position %u.  Converted to 'N'.\n",
                L, S[i], S[i], i);
        S[i] = 'N';
        Q[i] = '!';
        break;
    }
  }

  //  Load the qv header, and then load the qvs themselves over the header.
  fgets(Q, AS_MAX_READLEN, F->file());
  fgets(Q, AS_MAX_READLEN, F->file());
  chomp(Q);

  //  Convert from whatever QV encoding they have to the CA encoding
#warning ASSUMING READS ARE SANGER QV ENCODED

  for (uint32 i=0; Q[i]; i++) {
    Q[i] -= '!';
    Q[i] += '0';
  }

  //  Add a new read to the store.

  gkRead     *nr = gkpStore->gkStore_addEmptyRead(gkpLibrary);

  //  Populate the read with the appropriate gunk, based on the library type.

  gkReadData *nd = nr->gkRead_encodeSeqQlt(H, S, Q);

  //  Stash the encoded data in the store

  gkpStore->gkStore_stashReadData(nr, nd);

  //  And then toss it out.

  delete nd;

  //  Clear the lines, so we can load the next one.

  L[0] = 0;
  //H[0] = 0;
  //S[0] = 0;
  //Q[0] = 0;

  return(4);  //  FASTQ always reads exactly four lines
}





uint32
loadReads(gkStore *gkpStore, gkLibrary *gkpLibrary, FILE *nameMap, char *fileName) {
  char    *L = new char [AS_MAX_READLEN + 1];
  char    *H = new char [AS_MAX_READLEN + 1];

  uint32   Slen = 0;
  char    *S = new char [AS_MAX_READLEN + 1];

  uint32   Qlen = 0;
  char    *Q = new char [AS_MAX_READLEN + 1];

  uint64   lineNumber = 0;

  fprintf(stderr, "-- Loading reads from '%s'\n", fileName);

  compressedFileReader *F = new compressedFileReader(fileName);

  uint32   nFASTA = 0;
  uint32   nFASTQ = 0;
  uint32   nERROR = 0;

  fgets(L, AS_MAX_READLEN, F->file());
  chomp(L);

  while (!feof(F->file())) {

    if      (L[0] == '>') {
      lineNumber += loadFASTA(gkpStore, gkpLibrary, L, H, S, F);
      nFASTA++;

      fprintf(nameMap, F_U32"\t%s\n", gkpStore->gkStore_getNumReads(), H);
    }

    else if (L[0] == '@') {
      lineNumber += loadFASTQ(gkpStore, gkpLibrary, L, H, S, Q, F);
      nFASTQ++;

      fprintf(nameMap, F_U32"\t%s\n", gkpStore->gkStore_getNumReads(), H);
    }

    else {
      fprintf(stderr, "-- WARNING:  invalid read header '%s' in file '%s' at line "F_U64", skipping.\n",
              L, fileName, lineNumber);
      L[0] = 0;
      nERROR++;
    }

    //  If L[0] is nul, we need to load the next line.  If not, the next line is the header (from
    //  the fasta loader).

    if (L[0] == 0) {
      fgets(L, AS_MAX_READLEN, F->file());
      chomp(L);
    }
  }

  delete    F;

  delete [] Q;
  delete [] S;
  delete [] H;
  delete [] L;

  if (nFASTA > 0)
    fprintf(stderr, "-- Loaded "F_U32" FASTA format reads from '%s'.\n", nFASTA, fileName);
  if (nFASTQ > 0)
    fprintf(stderr, "-- Loaded "F_U32" FASTQ format reads from '%s'.\n", nFASTQ, fileName);

  return(0);
};




int
main(int argc, char **argv) {
  char            *gkpStoreName      = NULL;
  char            *outPrefix         = NULL;

  bool             ignoreClear       = true;

  uint32           firstFileArg      = 0;

  char             errorLogName[FILENAME_MAX];
  char             nameMapName[FILENAME_MAX];


  //argc = AS_configure(argc, argv);

  //fprintf(stderr, "gkLibrary: %u\n", sizeof(gkLibrary));
  //fprintf(stderr, "gkRead:    %u\n", sizeof(gkRead));

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-useclear") == 0) {
      ignoreClear = false;


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
    fprintf(stderr, "  -useclear           enforce clear ranges on the name lines\n");
    fprintf(stderr, "  \n");
    fprintf(stderr, "  \n");

    if (gkpStoreName == NULL)
      fprintf(stderr, "ERROR: no gkpStore (-g) supplied.\n");
    if (firstFileArg == 0)
      fprintf(stderr, "ERROR: no input files supplied.\n");

    exit(1);
  }


  gkStore     *gkpStore = new gkStore(gkpStoreName, gkStore_extend);
  gkRead      *gkpRead;
  gkReadData   gkpReadData;
  gkLibrary   *gkpLibrary;

  uint32       inLineLen = 1024;
  char         inLine[1024];

  errno = 0;

  sprintf(errorLogName, "%s.errorLog",    gkpStoreName);
  FILE    *errorLog = fopen(errorLogName, "w");
  if (errno)
    fprintf(stderr, "ERROR:  cannot open error file '%s': %s\n", errorLogName, strerror(errno)), exit(1);

  sprintf(nameMapName,   "%s/readNames.txt", gkpStoreName);
  FILE    *nameMap   = fopen(nameMapName,   "w");
  if (errno)
    fprintf(stderr, "ERROR:  cannot open uid map file '%s': %s\n", nameMapName, strerror(errno)), exit(1);

  uint32  nErrs  = 0;
  uint32  nWarns = 0;

  for (; firstFileArg < argc; firstFileArg++) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Starting file '%s'.\n", argv[firstFileArg]);

    compressedFileReader *inFile = new compressedFileReader(argv[firstFileArg]);
    char                 *line   = new char [10240];
    KeyAndValue           keyval;

    fgets(line, 10240, inFile->file());
    chomp(line);

    while (!feof(inFile->file())) {
      //fprintf(stderr, "LINE '%s'\n", line);

      keyval.find(line);

      if        (keyval.key() == NULL) {
        //  Just a damn comment line to ignore.

      } else if (strcasecmp(keyval.key(), "name") == 0) {
        gkpLibrary = gkpStore->gkStore_addEmptyLibrary(keyval.value());

      } else if (strcasecmp(keyval.key(), "preset") == 0) {
        gkpLibrary->gkLibrary_parsePreset(keyval.value());

      } else if (strcasecmp(keyval.key(), "qv") == 0) {
        gkpLibrary->gkLibrary_setDefaultQV(keyval.value_double());

      } else if (strcasecmp(keyval.key(), "isNonRandom") == 0) {
        gkpLibrary->gkLibrary_setIsNonRandom(keyval.value_bool());

      } else if (strcasecmp(keyval.key(), "trustHomopolymerRuns") == 0) {
        gkpLibrary->gkLibrary_setTrustHomopolymerRuns(keyval.value_bool());

      } else if (strcasecmp(keyval.key(), "initialTrim") == 0) {
        gkpLibrary->gkLibrary_setInitialTrim(keyval.value());

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
        nWarns += loadReads(gkpStore, gkpLibrary, nameMap, keyval.key());

      } else {
        fprintf(stderr, "ERROR:  option '%s' not recognized, and not a file of reads.\n", line);
        nErrs++;
      }

      fgets(line, 10240, inFile->file());
      chomp(line);
    }

    delete    inFile;
    delete [] line;
  }

  delete gkpStore;

  fclose(nameMap);
  fclose(errorLog);

  fprintf(stderr, "\n");

  if (nErrs > 0) {
    fprintf(stderr, "gatekeeperCreate did NOT finish successfully; %u errors and %u warnings.\n", nErrs, nWarns);
    exit(1);
  }

  //return(AS_GKP_summarizeErrors());
  fprintf(stderr, "gatekeeperCreate finished successfully; %u errors and %u warnings.\n", nErrs, nWarns);

  exit(0);
}
