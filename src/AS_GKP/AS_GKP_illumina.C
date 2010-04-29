
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

static char const *rcsid = "$Id: AS_GKP_illumina.C,v 1.11 2010-04-29 17:57:00 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "AS_global.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_encodeSequenceQuality.h"
#include "AS_UTL_reverseComplement.h"

#define FASTQ_SANGER    0
#define FASTQ_SOLEXA    1
#define FASTQ_ILLUMINA  2

#define FASTQ_INNIE     0
#define FASTQ_OUTTIE    1

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
processSeq(char *N, ilFragment *fr, char end, uint32 fastqType, uint32 fastqOrient) {

  fr->fr.gkFragment_setType(GKFRAGMENT_PACKED);
  fr->fr.gkFragment_setIsDeleted(1);

  fr->fr.gkFragment_setMateIID(0);
  fr->fr.gkFragment_setLibraryIID(0);

  //  Clean up what we read.  Remove trailing newline, truncate read names to the first word.

  chomp(fr->snam);
  chomp(fr->sstr);
  chomp(fr->qnam);
  chomp(fr->qstr);

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

  uint32   slen = strlen(fr->sstr);
  uint32   qlen = strlen(fr->qstr);

  uint32   clrL=0, clrR=slen;

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

  //  Convert QVs

  if (fastqType == FASTQ_SANGER) {
    for (uint32 i=0; fr->qstr[i]; i++) {
      if (fr->qstr[i] < '!') {
        AS_GKP_reportError(AS_GKP_ILL_BAD_QV, fr->snam, fr->qstr[i], "sanger");
        return(0);
      }
      fr->qstr[i] -= '!';
      if (fr->qstr[i] > QUALITY_MAX)
        fr->qstr[i] = QUALITY_MAX;
      fr->qstr[i] += '0';
    }
  }

  if (fastqType == FASTQ_SOLEXA) {
    double qs;
    for (uint32 i=0; fr->qstr[i]; i++) {
      if (fr->qstr[i] < ';') {
        AS_GKP_reportError(AS_GKP_ILL_BAD_QV, fr->snam, fr->qstr[i], "solexa");
        return(0);
      }
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
        AS_GKP_reportError(AS_GKP_ILL_BAD_QV, fr->snam, fr->qstr[i], "illumina");
        return(0);
      }
      fr->qstr[i] -= '@';
      if (fr->qstr[i] > QUALITY_MAX)
        fr->qstr[i] = QUALITY_MAX;
      fr->qstr[i] += '0';
    }
  }

  //  Reverse the read if it is from an outtie pair.  This ONLY works if the reads are the same
  //  length throughout the library.  WE DO NOT CHECK THAT IS SO.

  if (fastqOrient == FASTQ_OUTTIE)
    reverseComplement(fr->sstr, fr->qstr, slen);

  //  Trim crap off the ends.

  while ((fr->qstr[clrR-1] < '6') && (clrR > 1))
    clrR--;

  while ((fr->qstr[clrL] < '6') && (clrL < clrR))
    clrL++;

  assert(clrL <= clrR);

  if (clrR - clrL < AS_READ_MIN_LEN)
    return(0);

  //  Make sure there aren't any bogus letters

  for (uint32 i=0; i<slen; i++)
    if (fr->sstr[i] == '.')
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

  fr->fr.maxBgn = 1;  //  No max yet.
  fr->fr.maxEnd = 0;

  fr->fr.vecBgn = 1;  //  There is no vector clear defined for Illumina reads.
  fr->fr.vecEnd = 0;

  fr->fr.tntBgn = 1;  //  Nothing contaminated.
  fr->fr.tntEnd = 0;

  return(readUID);
}


static
uint64
readQSeq(FILE *F, char *N, ilFragment *fr, char end, uint32 fastqType, uint32 fastqOrient) {

  fr->fr.gkFragment_setType(GKFRAGMENT_PACKED);
  fr->fr.gkFragment_setIsDeleted(1);

  fr->fr.gkFragment_setMateIID(0);
  fr->fr.gkFragment_setLibraryIID(0);

  fgets(fr->qstr, AS_READ_MAX_NORMAL_LEN, F);  chomp(fr->qstr);

  if (feof(F))
    return(0);

  char *v[32];
  int   nv = 0;
  char *p;
  int   is = 1;

  for (char *p=fr->qstr; *p; p++) {
    if (isspace(*p)) {
      *p = 0;
      is = 1;
    } else {
      if (is == 1) {
        v[nv++] = p;
        is = 0;
      }
    }
  }

  if (nv != 11)
    fprintf(stderr, "ERROR:  qseq not in expected format.  Please convert to standard FastQ.\n");
  assert(nv == 11);

  sprintf(fr->snam, "@%s:%s:%s:%s:%s#%s/%s", v[0], v[2], v[3], v[4], v[5], v[6], v[7]);
  sprintf(fr->sstr, "%s", v[8]);
  sprintf(fr->qnam, "+%s:%s:%s:%s:%s#%s/%s", v[0], v[2], v[3], v[4], v[5], v[6], v[7]);
  sprintf(fr->qstr, "%s", v[9]);

  return(processSeq(N, fr, end, fastqType, fastqType));
}


static
uint64
readSeq(FILE *F, char *N, ilFragment *fr, char end, uint32 fastqType, uint32 fastqOrient) {

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

  return(processSeq(N, fr, end, fastqType, fastqType));
}


static
bool
openFile(char *name, FILE *&file) {
  char *cmd  = new char [strlen(name) + 256];
  bool  pipe = false;

  errno = 0;
  if        (name == NULL) {
    fprintf(stderr, "ERROR:  Failed to open illumina input file: no name supplied.\n");
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
    fprintf(stderr, "ERROR:  Failed to open illumina input file '%s': %s\n", name, strerror(errno));
    AS_GKP_reportError(AS_GKP_ILL_CANT_OPEN_INPUT, name, strerror(errno));
    exit(1);
  }

  delete [] cmd;

  return(pipe);
}


static
void
loadIlluminaReads(char *lname, char *rname, bool isSeq, uint32 fastqType, uint32 fastqOrient) {
  fprintf(stderr, "Processing illumina reads from '%s'\n", lname);
  fprintf(stderr, "                           and '%s'\n", rname);

  FILE *lfile = NULL;
  bool  lpipe = openFile(lname, lfile);
  FILE *rfile = NULL;
  bool  rpipe = openFile(rname, rfile);

  ilFragment  *lfrg = new ilFragment;
  ilFragment  *rfrg = new ilFragment;

  lfrg->fr.gkFragment_enableGatekeeperMode(gkpStore);
  rfrg->fr.gkFragment_enableGatekeeperMode(gkpStore);

  while (!feof(lfile) && !feof(rfile)) {
    uint32 nfrg = gkpStore->gkStore_getNumFragments();

    uint64 lUID = (isSeq) ? readSeq(lfile, lname, lfrg, 'l', fastqType, fastqOrient) : readQSeq(lfile, lname, lfrg, 'l', fastqType, fastqOrient);
    uint64 rUID = (isSeq) ? readSeq(rfile, rname, rfrg, 'r', fastqType, fastqOrient) : readQSeq(rfile, rname, rfrg, 'r', fastqType, fastqOrient);

    if       ((lfrg->fr.gkFragment_getIsDeleted() == 0) &&
              (rfrg->fr.gkFragment_getIsDeleted() == 0)) {
      //  Both OK, add a mated read.
      lfrg->fr.gkFragment_setMateIID(nfrg + 2);
      rfrg->fr.gkFragment_setMateIID(nfrg + 1);

      lfrg->fr.gkFragment_setOrientation(AS_READ_ORIENT_INNIE);
      rfrg->fr.gkFragment_setOrientation(AS_READ_ORIENT_INNIE);

      gkpStore->gkStore_addFragment(&lfrg->fr);
      gkpStore->gkStore_addFragment(&rfrg->fr);


    } else if (lfrg->fr.gkFragment_getIsDeleted() == 0) {
      //  Only add the left fragment.
      gkpStore->gkStore_addFragment(&lfrg->fr);


    } else if (rfrg->fr.gkFragment_getIsDeleted() == 0) {
      //  Only add the right fragment.
      gkpStore->gkStore_addFragment(&rfrg->fr);


    } else {
      //  Both deleted, do nothing.
    }
  }

  delete lfrg;
  delete rfrg;

  if (lpipe)  pclose(lfile);  else  fclose(lfile);
  if (rpipe)  pclose(rfile);  else  fclose(rfile);
}



static
void
loadIlluminaReads(char *uname, bool isSeq, uint32 fastqType, uint32 fastqOrient) {
  fprintf(stderr, "Processing illumina reads from '%s'.\n", uname);

  FILE *ufile = NULL;
  bool  upipe = openFile(uname, ufile);

  ilFragment  *ufrg = new ilFragment;

  ufrg->fr.gkFragment_enableGatekeeperMode(gkpStore);

  while (!feof(ufile)) {
    uint32 nfrg = gkpStore->gkStore_getNumFragments();

    uint64 uUID = (isSeq) ? readSeq(ufile, uname, ufrg, 'u', fastqType, fastqOrient) : readQSeq(ufile, uname, ufrg, 'u', fastqType, fastqOrient);

    if (ufrg->fr.gkFragment_getIsDeleted() == 0) {
      //  Add a fragment.
      gkpStore->gkStore_addFragment(&ufrg->fr);

    } else {
      //  Junk read, do nothing.
    }
  }

  delete ufrg;

  if (upipe)  pclose(ufile);  else  fclose(ufile);
}




void
checkLibraryForIlluminaPointers(LibraryMesg *lib_mesg) {
  uint32  fastqType   = FASTQ_SOLEXA;
  uint32  fastqOrient = FASTQ_INNIE;

  //  Search for the type of the reads.
  for (uint32 i=0; i<lib_mesg->num_features; i++) {
    if (strcasecmp(lib_mesg->features[i], "illuminaFastQType") == 0) {
      if (strcasecmp(lib_mesg->values[i], "sanger") == 0) {
        fprintf(stderr, "Type set to SANGER.\n");
        fastqType = FASTQ_SANGER;
      }
      if (strcasecmp(lib_mesg->values[i], "solexa") == 0) {
        fprintf(stderr, "Type set to SOLEXA, pre-1.3.\n");
        fastqType = FASTQ_SOLEXA;
      }
      if (strcasecmp(lib_mesg->values[i], "illumina") == 0) {
        fprintf(stderr, "Type set to ILLUMINA 1.3+.\n");
        fastqType = FASTQ_ILLUMINA;
      }
    }
  }

  //  Search for the orientation of the reads.
  for (uint32 i=0; i<lib_mesg->num_features; i++) {
    if (strcasecmp(lib_mesg->features[i], "illuminaOrientation") == 0) {
      if (strcasecmp(lib_mesg->values[i], "innie") == 0) {
        fprintf(stderr, "Orientation set to INNIE.\n");
        fastqOrient = FASTQ_INNIE;
      }
      if (strcasecmp(lib_mesg->values[i], "outtie") == 0) {
        fprintf(stderr, "Orientation set to OUTTIE.\n");
        fastqOrient = FASTQ_OUTTIE;
      }
    }
  }

  //  Search for and load the reads.
  for (uint32 i=0; i<lib_mesg->num_features; i++) {
    if (strcasecmp(lib_mesg->features[i], "illuminaQSequence") == 0) {
      char *sl = lib_mesg->values[i];
      char *sr = lib_mesg->values[i];

      while ((*sr) && (*sr != ','))
        sr++;

      if (*sr) {
        *sr = 0;
        sr++;
        loadIlluminaReads(sl, sr, false, fastqType, fastqOrient);
      } else {
        loadIlluminaReads(sl, false, fastqType, fastqOrient);
      }
    }

    if (strcasecmp(lib_mesg->features[i], "illuminaSequence") == 0) {
      char *sl = lib_mesg->values[i];
      char *sr = lib_mesg->values[i];

      while ((*sr) && (*sr != ','))
        sr++;

      if (*sr) {
        *sr = 0;
        sr++;
        loadIlluminaReads(sl, sr, true, fastqType, fastqOrient);
      } else {
        loadIlluminaReads(sl, true, fastqType, fastqOrient);
      }
    }
  }
}
