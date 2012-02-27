
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2009, J. Craig Venter Institute
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

static char const *rcsid = "$Id: AS_GKP_illumina.C,v 1.30 2012-02-27 22:45:57 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "AS_global.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_encodeSequenceQuality.h"
#include "AS_UTL_reverseComplement.h"
#include "AS_UTL_fasta.h"

#define FASTQ_SANGER    0
#define FASTQ_SOLEXA    1
#define FASTQ_ILLUMINA  2

#define FASTQ_INNIE     0
#define FASTQ_OUTTIE    1

//  This seems to ALWAYS be a bad idea.
#undef FASTQ_TRIM_JUNK

static int *isValidACGTN = NULL;

class ilFragment {
 public:
  char        snam[AS_READ_MAX_NORMAL_LEN];
  char        sstr[AS_READ_MAX_NORMAL_LEN];
  char        qnam[AS_READ_MAX_NORMAL_LEN];
  char        qstr[AS_READ_MAX_NORMAL_LEN];
  gkFragment  fr;
};


static
uint64
processSeq(char       *N,
           ilFragment *fr,
           char        end,
           uint32      fastqType,
           uint32      fastqOrient,
           bool        forcePacked,
           uint32      packedLength) {

  chomp(fr->snam);
  chomp(fr->sstr);
  chomp(fr->qnam);
  chomp(fr->qstr);

  uint32   slen = strlen(fr->sstr);
  uint32   qlen = strlen(fr->qstr);

  //  Determine if this should be a PACKED or a NORMAL fragment.  PACKED is definitely
  //  preferred for assemblies that are using most reads of a single length.

  fr->fr.gkFragment_clear();

  if (forcePacked == true)
    fr->fr.gkFragment_setType(GKFRAGMENT_PACKED);
  else
    fr->fr.gkFragment_setType(GKFRAGMENT_NORMAL);

  fr->fr.gkFragment_setIsDeleted(1);

  //  If we're packed, make sure the length is appropriate

  if ((forcePacked) &&
      (slen > packedLength)) {
    AS_GKP_reportError(AS_GKP_ILL_SEQ_TOO_LONG, fr->snam, slen, packedLength);
    fr->sstr[packedLength] = 0;
    fr->qstr[packedLength] = 0;
    slen = packedLength;
    qlen = packedLength;
  }
      
  //  Complicated, but fast, parsing of the 'snam' to find clear ranges.

  char   *tok = strtok(fr->snam, " \t");

  uint32   clrL=0, clrR=slen;  //  Defined range, whole read
  uint32   clvL=1, clvR=0;     //  Undefined range
  uint32   clmL=1, clmR=0;
  uint32   tntL=1, tntR=0;

  while (tok != NULL) {
    if (((tok[0] == 'c') || (tok[0] == 'C')) &&
        ((tok[1] == 'l') || (tok[1] == 'L')) &&
        ((tok[2] == 'r') || (tok[2] == 'R'))) {
      clrL = atoi(tok + 4);
      while ((*tok) && (*tok != ','))
        tok++;
      clrR = atoi(tok + 1);
    }
    if (((tok[0] == 'c') || (tok[0] == 'C')) &&
        ((tok[1] == 'l') || (tok[1] == 'L')) &&
        ((tok[2] == 'v') || (tok[2] == 'V'))) {
      clvL = atoi(tok + 4);
      while ((*tok) && (*tok != ','))
        tok++;
      clvR = atoi(tok + 1);
    }
    if (((tok[0] == 'm') || (tok[0] == 'M')) &&
        ((tok[1] == 'a') || (tok[1] == 'A')) &&
        ((tok[2] == 'x') || (tok[2] == 'X'))) {
      clmL = atoi(tok + 4);
      while ((*tok) && (*tok != ','))
        tok++;
      clmR = atoi(tok + 1);
    }
    if (((tok[0] == 't') || (tok[0] == 'T')) &&
        ((tok[1] == 'n') || (tok[1] == 'N')) &&
        ((tok[2] == 't') || (tok[2] == 'T'))) {
      tntL = atoi(tok + 4);
      while ((*tok) && (*tok != ','))
        tok++;
      tntR = atoi(tok + 1);
    }
    if (((tok[0] == 'r') || (tok[0] == 'R')) &&
        ((tok[1] == 'n') || (tok[1] == 'N')) &&
        ((tok[2] == 'd') || (tok[2] == 'D')) &&
        ((tok[4] == 'f') || (tok[4] == 'F'))) {
      fr->fr.gkFragment_setIsNonRandom(1);
    }

    tok = strtok(NULL, " \t");
  }

  if (clrR > slen) clrR = slen;    if (clrL > clrR) clrL = clrR;
  if (clvR > slen) clvR = slen;    if (clvL > clvR) clvL = clvR;
  if (clmR > slen) clmR = slen;    if (clmL > clmR) clmL = clmR;
  if (tntR > slen) tntR = slen;    if (tntL > tntR) tntL = tntR;


  //  Clean up what we read.  Remove trailing newline (whoops, already done), truncate read names to
  //  the first word.

  for (uint32 i=0; fr->snam[i]; i++)
    if (isspace(fr->snam[i])) {
      fr->snam[i] = 0;
      break;
    }

  for (uint32 i=0; fr->qnam[i]; i++)
    if (isspace(fr->qnam[i])) {
      fr->qnam[i] = 0;
      break;
    }

  //  Check that things are consistent.  Same names, same lengths, etc.

  if (fr->snam[0] != '@') {
    AS_GKP_reportError(AS_GKP_ILL_NOT_SEQ_START_LINE, N, fr->snam);
    return(0);
  }

  if (fr->qnam[0] != '+') {
    AS_GKP_reportError(AS_GKP_ILL_NOT_QLT_START_LINE, N, fr->qnam);
    return(0);
  }

  if ((fr->qnam[1] != 0) && (strcmp(fr->snam+1, fr->qnam+1) != 0)) {
    AS_GKP_reportError(AS_GKP_ILL_SEQ_QLT_NAME_DIFFER, N, fr->snam, fr->qnam);
    return(0);
  }

  if (slen != qlen) {
    AS_GKP_reportError(AS_GKP_ILL_SEQ_QLT_LEN_DIFFER, N, fr->snam, slen, qlen);
    return(0);
  }

  //  Convert QVs and check for errors
  //
  //  Sanger       range from ! to I
  //  Solexa       range from ; to h (with '@' == QV 0)
  //  Illumina 1.3 range from @ to h
  //  Illumina 1.5 range from B to h (B = Read Segment QC Indicator)
  //  Illumina 1.8 range from ! to J

  uint32 QVerrors = 0;

  if (fastqType == FASTQ_SANGER) {
    for (uint32 i=0; fr->qstr[i]; i++) {
      if (fr->qstr[i] < '!') {
        QVerrors++;
        fr->qstr[i] = '!';
      }
      //if ('I' < fr->qstr[i]) {
      //  QVerrors++;
      //  fr->qstr[i] = 'I';
      //}

      fr->qstr[i] -= '!';
      if (fr->qstr[i] > QUALITY_MAX)
        fr->qstr[i] = QUALITY_MAX;
      fr->qstr[i] += '0';
    }
  }

  if (fastqType == FASTQ_SOLEXA) {
    double qs;
    for (uint32 i=0; fr->qstr[i]; i++) {
      if (fr->qstr[i] < '@') {
        QVerrors++;
        fr->qstr[i] = '@';
      }
      //if ('h' < fr->qstr[i]) {
      //  QVerrors++;
      //  fr->qstr[i] = ';';
      //}
      qs  = fr->qstr[i];
      qs -= '@';
      qs /= 10.0;
      qs  = 10.0 * log10(pow(10.0, qs) + 1);
      if (qs > QUALITY_MAX)
        qs = QUALITY_MAX;
      fr->qstr[i] = lround(qs) + '0';
    }
  }

  if (fastqType == FASTQ_ILLUMINA) {
    for (uint32 i=0; fr->qstr[i]; i++) {
      if (fr->qstr[i] < '@') {
        QVerrors++;
        fr->qstr[i] = '@';
      }
      //if ('h' < fr->qstr[i]) {
      //  QVerrors++;
      //  fr->qstr[i] = 'h';
      //}
      fr->qstr[i] -= '@';
      if (fr->qstr[i] > QUALITY_MAX)
        fr->qstr[i] = QUALITY_MAX;
      fr->qstr[i] += '0';
    }
  }

  if (QVerrors > 0)
    AS_GKP_reportError(AS_GKP_ILL_BAD_QV, fr->snam, QVerrors);


  //  Reverse the read if it is from an outtie pair.  This ONLY works if the reads are the same
  //  length throughout the library.  WE DO NOT CHECK THAT IS SO.

  if (fastqOrient == FASTQ_OUTTIE)
    reverseComplement(fr->sstr, fr->qstr, slen);

  //  Make sure there aren't any bogus letters

  for (uint32 i=0; i<slen; i++)
    if (isValidACGTN[fr->sstr[i]] == 0) {
      fr->sstr[i] = 'N';
      fr->qstr[i] = '0';
    }

  //  Trim crap off the ends.

#ifdef FASTQ_TRIM_JUNK
  while ((clrR > 0) && (fr->qstr[clrR-1] < '6'))
    clrR--;

  while ((clrL < clrR) && (fr->qstr[clrL] < '6'))
    clrL++;
#endif

  assert(clrL <= clrR);

  if (clrR - clrL < AS_READ_MIN_LEN)
    return(0);

  //  Construct a UID for this read
  //
  //  A 64-bit unsigned holds 2^64 = 4611686018427387904
  //
  //  Our store holds at most 2^31 fragments = 2 billion.
  //
  //  The UIDs are constructed as:
  //  4611686018427387904
  //           2147483648
  //  LLLLLLLLR##########
  //

  uint64  libraryIID = gkpStore->gkStore_getNumLibraries();
  uint64  readUID    = 0;

  //  The plus one is to make the UID and IID match up when these are the first fragments in the
  //  store.  Pointless otherwise.

  switch (end) {
    case 'l':
      readUID = (libraryIID * 10 + 1) * 10000000000LLU + gkpStore->gkStore_getNumFragments() + 1;
      break;
    case 'r':
      readUID = (libraryIID * 10 + 2) * 10000000000LLU + gkpStore->gkStore_getNumFragments() + 1;
      break;
    case 'u':
      readUID = (libraryIID * 10 + 0) * 10000000000LLU + gkpStore->gkStore_getNumFragments() + 1;
      break;
    default:
      readUID = 0;
      break;
  }

  //  Got a good read, make it.

  fr->fr.gkFragment_setReadUID(AS_UID_fromInteger(readUID));
  fr->fr.gkFragment_setIsDeleted(0);

  fr->fr.gkFragment_setLibraryIID(libraryIID);
  fr->fr.gkFragment_setOrientation(AS_READ_ORIENT_UNKNOWN);

  memcpy(fr->fr.gkFragment_getSequence(), fr->sstr, sizeof(char) * slen);
  memcpy(fr->fr.gkFragment_getQuality(),  fr->qstr, sizeof(char) * qlen);

  fr->fr.gkFragment_setLength(slen);

  fr->fr.gkFragment_getSequence()[slen] = 0;
  fr->fr.gkFragment_getQuality() [qlen] = 0;

  fr->fr.clrBgn = clrL;
  fr->fr.clrEnd = clrR;

  fr->fr.maxBgn = clmL;
  fr->fr.maxEnd = clmR;

  fr->fr.vecBgn = clvL;
  fr->fr.vecEnd = clvR;

  fr->fr.tntBgn = tntL;
  fr->fr.tntEnd = tntR;

  return(readUID);
}



static
uint64
readSeq(FILE       *F,
        char       *N,
        ilFragment *fr,
        char        end,
        uint32      fastqType,
        uint32      fastqOrient,
        bool        forcePacked,
        uint32      packedLength) {

  fr->fr.gkFragment_setType(GKFRAGMENT_PACKED);
  fr->fr.gkFragment_setIsDeleted(1);

  fr->fr.gkFragment_setMateIID(0);
  fr->fr.gkFragment_setLibraryIID(0);

  fgets(fr->snam, AS_READ_MAX_NORMAL_LEN, F);  chomp(fr->snam);
  fgets(fr->sstr, AS_READ_MAX_NORMAL_LEN, F);  chomp(fr->sstr);
  fgets(fr->qnam, AS_READ_MAX_NORMAL_LEN, F);  chomp(fr->qnam);
  fgets(fr->qstr, AS_READ_MAX_NORMAL_LEN, F);  chomp(fr->qstr);

  if (feof(F))
    return(0);

  return(processSeq(N, fr, end, fastqType, fastqOrient, forcePacked, packedLength));
}


static
bool
openFile(char *name, FILE *&file) {
  char *cmd  = new char [strlen(name) + 256];
  bool  pipe = false;

  errno = 0;

  if (AS_UTL_fileExists(name, FALSE, FALSE) == FALSE) {
    fprintf(stderr, "ERROR:  Failed to open fastq input file '%s': %s\n", name, strerror(errno));
    AS_GKP_reportError(AS_GKP_ILL_CANT_OPEN_INPUT, name, strerror(errno));
    exit(1);
  }

  if        (name == NULL) {
    fprintf(stderr, "ERROR:  Failed to open fastq input file: no name supplied.\n");
    AS_GKP_reportError(AS_GKP_ILL_CANT_OPEN_INPUT, "(no-name-supplied)", "no name supplied");
    exit(1);

  } else if (strcasecmp(name+strlen(name)-3, ".gz") == 0) {
    sprintf(cmd, "gzip -dc %s", name);
    file = popen(cmd, "r");
    pipe = true;

  } else if (strcasecmp(name+strlen(name)-4, ".bz2") == 0) {
    sprintf(cmd, "bzip2 -dc %s", name);
    file = popen(cmd, "r");
    pipe = true;

  } else {
    file = fopen(name, "r");
    pipe = false;
  }

  if (errno) {
    fprintf(stderr, "ERROR:  Failed to open fastq input file '%s': %s\n", name, strerror(errno));
    AS_GKP_reportError(AS_GKP_ILL_CANT_OPEN_INPUT, name, strerror(errno));
    exit(1);
  }

  delete [] cmd;

  return(pipe);
}


static
void
loadFastQReads(char    *lname,
               char    *rname,
               uint32   fastqType,
               uint32   fastqOrient,
               bool     forcePacked,
               uint32   packedLength) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Processing %s %s QV encoding reads from:\n",
          (fastqOrient == FASTQ_INNIE) ? "INNIE" : "OUTTIE",
          (fastqType   == FASTQ_ILLUMINA) ? "ILLUMINA 1.3+" : ((fastqType == FASTQ_SANGER) ? "SANGER" : "SOLEXA pre-1.3"));
  if (lname == rname) {
    fprintf(stderr, "      '%s' (INTERLACED)\n", lname);
  } else {
    fprintf(stderr, "      '%s'\n", lname);
    fprintf(stderr, "  and '%s'\n", rname);
  }

  if (fastqUIDmap == NULL) {
    errno = 0;
    fastqUIDmap = fopen(fastqUIDmapName, "w");
    if (errno) {
      fprintf(stderr, "cannot open fastq UID map file '%s': %s\n", fastqUIDmapName, strerror(errno));
      exit(1);
    }
  }

  FILE *lfile = NULL, *rfile = NULL;
  bool  lpipe = false, rpipe = false;

  if (strcmp(lname, rname) == 0) {
    lpipe = openFile(lname, lfile);
    rfile = lfile;
  } else {
    lpipe = openFile(lname, lfile);
    rpipe = openFile(rname, rfile);
  }

  ilFragment  *lfrg = new ilFragment;
  ilFragment  *rfrg = new ilFragment;

  lfrg->fr.gkFragment_enableGatekeeperMode(gkpStore);
  rfrg->fr.gkFragment_enableGatekeeperMode(gkpStore);

  while (!feof(lfile) && !feof(rfile)) {
    uint32 nfrg = gkpStore->gkStore_getNumFragments();

    uint64 lUID = readSeq(lfile, lname, lfrg, 'l', fastqType, fastqOrient, forcePacked, packedLength);
    uint64 rUID = readSeq(rfile, rname, rfrg, 'r', fastqType, fastqOrient, forcePacked, packedLength);

    if       ((lfrg->fr.gkFragment_getIsDeleted() == 0) &&
              (rfrg->fr.gkFragment_getIsDeleted() == 0)) {
      //  Both OK, add a mated read.
      lfrg->fr.gkFragment_setMateIID(nfrg + 2);
      rfrg->fr.gkFragment_setMateIID(nfrg + 1);

      lfrg->fr.gkFragment_setOrientation(AS_READ_ORIENT_INNIE);
      rfrg->fr.gkFragment_setOrientation(AS_READ_ORIENT_INNIE);

      gkpStore->gkStore_addFragment(&lfrg->fr);
      gkpStore->gkStore_addFragment(&rfrg->fr);

      fprintf(fastqUIDmap, F_U64"\t"F_U32"\t%s\t"F_U64"\t"F_U32"\t%s\n",
              lUID, nfrg + 1, lfrg->snam+1,
              rUID, nfrg + 2, rfrg->snam+1);

    } else if (lfrg->fr.gkFragment_getIsDeleted() == 0) {
      //  Only add the left fragment.
      gkpStore->gkStore_addFragment(&lfrg->fr);

      fprintf(fastqUIDmap, F_U64"\t"F_U32"\t%s\n",
              lUID, nfrg + 1, lfrg->snam+1);


    } else if (rfrg->fr.gkFragment_getIsDeleted() == 0) {
      //  Only add the right fragment.
      gkpStore->gkStore_addFragment(&rfrg->fr);

      fprintf(fastqUIDmap, F_U64"\t"F_U32"\t%s\n",
              rUID, nfrg + 1, rfrg->snam+1);


    } else {
      //  Both deleted, do nothing.
    }
  }

  delete lfrg;
  delete rfrg;

  if (strcmp(lname, rname) == 0) {
    if (lpipe)  pclose(lfile);  else  fclose(lfile);
  } else {
    if (lpipe)  pclose(lfile);  else  fclose(lfile);
    if (rpipe)  pclose(rfile);  else  fclose(rfile);
  }
}



static
void
loadFastQReads(char    *uname,
               uint32   fastqType,
               uint32   fastqOrient,
               bool     forcePacked,
               uint32   packedLength) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Processing SINGLE-ENDED %s QV encoding reads from:\n",
          (fastqType   == FASTQ_ILLUMINA) ? "ILLUMINA 1.3+" : ((fastqType == FASTQ_SANGER) ? "SANGER" : "SOLEXA pre-1.3"));
  fprintf(stderr, "      '%s'\n", uname);

  if (fastqUIDmap == NULL) {
    errno = 0;
    fastqUIDmap = fopen(fastqUIDmapName, "w");
    if (errno) {
      fprintf(stderr, "cannot open fastq UID map file '%s': %s\n", fastqUIDmapName, strerror(errno));
      exit(1);
    }
  }

  FILE *ufile = NULL;
  bool  upipe = openFile(uname, ufile);

  ilFragment  *ufrg = new ilFragment;

  ufrg->fr.gkFragment_enableGatekeeperMode(gkpStore);

  while (!feof(ufile)) {
    uint32 nfrg = gkpStore->gkStore_getNumFragments();

    uint64 uUID = readSeq(ufile, uname, ufrg, 'u', fastqType, fastqOrient, forcePacked, packedLength);

    if (ufrg->fr.gkFragment_getIsDeleted() == 0) {
      //  Add a fragment.
      gkpStore->gkStore_addFragment(&ufrg->fr);

      fprintf(fastqUIDmap, F_U64"\t"F_U32"\t%s\n",
              uUID, nfrg + 1, ufrg->snam+1);

    } else {
      //  Junk read, do nothing.
    }
  }

  delete ufrg;

  if (upipe)  pclose(ufile);  else  fclose(ufile);
}




void
checkLibraryForFastQPointers(LibraryMesg *lib_mesg, uint32 packedLength) {
  uint32  fastqType   = FASTQ_SOLEXA;
  uint32  fastqOrient = FASTQ_INNIE;
  bool    forcePacked = false;

  isValidACGTN = AS_UTL_getValidACGTN();

  for (uint32 i=0; i<lib_mesg->num_features; i++) {
    if (strncasecmp(lib_mesg->features[i], "illumina", 8) == 0) {
      fprintf(stderr, "ERROR:  Obsolete LIB feature '%s' detected.\n",
              lib_mesg->features[i]);
      fprintf(stderr, "The presence of an obsolete flag indicates that the FRG file was generated ");
      fprintf(stderr, "by an older version of Celera Assembler.\n");
      fprintf(stderr, "Please re-generate your FRG file with this version of Celera Assembler.\n");
      exit(1);
    }
  }

  //  Search for forced normal - we technically should grab the library from gatekeeper
  //  and use the flag there, but this is easier.
  for (uint32 i=0; i<lib_mesg->num_features; i++) {
    if (strcasecmp(lib_mesg->features[i], "forceShortReadFormat") == 0) {
      if ((lib_mesg->values[i][0] == '1') ||
          (lib_mesg->values[i][0] == 't') ||
          (lib_mesg->values[i][0] == 'T'))
        forcePacked = true;
    }
  }

  //  Search for the type of the reads.
  for (uint32 i=0; i<lib_mesg->num_features; i++) {
    if (strcasecmp(lib_mesg->features[i], "fastqQualityValues") == 0) {
      if (strcasecmp(lib_mesg->values[i], "sanger") == 0)
        fastqType = FASTQ_SANGER;
      if (strcasecmp(lib_mesg->values[i], "solexa") == 0)
        fastqType = FASTQ_SOLEXA;
      if (strcasecmp(lib_mesg->values[i], "illumina") == 0)
        fastqType = FASTQ_ILLUMINA;
    }
  }

  //  Search for the orientation of the reads.
  for (uint32 i=0; i<lib_mesg->num_features; i++) {
    if (strcasecmp(lib_mesg->features[i], "fastqOrientation") == 0) {
      if (strcasecmp(lib_mesg->values[i], "innie") == 0)
        fastqOrient = FASTQ_INNIE;
      if (strcasecmp(lib_mesg->values[i], "outtie") == 0)
        fastqOrient = FASTQ_OUTTIE;
    }
  }

  //  Search for and load the reads.
  for (uint32 i=0; i<lib_mesg->num_features; i++) {
    if (strcasecmp(lib_mesg->features[i], "fastqMates") == 0) {
      char *sl = lib_mesg->values[i];
      char *sr = lib_mesg->values[i];

      while ((*sr) && (*sr != ','))
        sr++;

      if (*sr)
        *sr++ = 0;

      if (*sr == 0)
        loadFastQReads(sl, sl, fastqType, fastqOrient, forcePacked, packedLength);
      else
        loadFastQReads(sl, sr, fastqType, fastqOrient, forcePacked, packedLength);
    }

    if (strcasecmp(lib_mesg->features[i], "fastqReads") == 0) {
      char *sl = lib_mesg->values[i];
      char *sr = lib_mesg->values[i];

      while ((*sr) && (*sr != ','))
        sr++;

      if (*sr)
        fprintf(stderr, "ERROR:  Multipe FastQ files given to fastqReads.\n"), exit(1);

      loadFastQReads(sl, fastqType, fastqOrient, forcePacked, packedLength);
    }
  }
}
