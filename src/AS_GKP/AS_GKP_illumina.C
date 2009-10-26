
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

static char const *rcsid = "$Id: AS_GKP_illumina.C,v 1.2 2009-10-26 13:20:26 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>

#include "AS_global.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"


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
processSeq(char *N, ilFragment *fr, char end) {

  fr->fr.gkFragment_setType(GKFRAGMENT_PACKED);
  fr->fr.gkFragment_setIsDeleted(1);

  fr->fr.gkFragment_setMateIID(0);
  fr->fr.gkFragment_setLibraryIID(0);

  chomp(fr->snam);
  chomp(fr->sstr);
  chomp(fr->qnam);
  chomp(fr->qstr);

  uint32   slen = strlen(fr->sstr);
  uint32   qlen = strlen(fr->qstr);

  uint32   clrL=0, clrR=slen;

  if (fr->snam[0] != '@')
    fprintf(stderr, "ERROR:  file '%s': seq name '%s' is not a sequence start line.\n", N, fr->snam);
  if (fr->qnam[0] != '+')
    fprintf(stderr, "ERROR:  file '%s': qlt name '%s' is not a quality start line.\n", N, fr->qnam);

  if (strcmp(fr->snam+1, fr->qnam+1) != 0)
    fprintf(stderr, "ERROR:  file: '%s': seq/qlt names differ; seq='%s' qlt='%s'\n", N, fr->snam, fr->qnam);

  if (slen != qlen)
    fprintf(stderr, "ERROR:  file '%s': seq/qlt lengths differ for read '%s'; seq=%d qlt=%d\n", N, fr->snam, slen, qlen);

  //  Convert QVs

  for (int32 i=0; i<qlen; i++)
    fr->qstr[i] -= 'B' - '0';

  //  Trim crap off the ends.

  while ((fr->qstr[clrR-1] <= '0') && (clrR > 1))
    clrR--;

  while ((fr->qstr[clrL] <= '0') && (clrL < clrR))
    clrL++;

  assert(clrL <= clrR);

  if (clrR - clrL < AS_READ_MIN_LEN)
    return(0);

  //  Make sure there aren't any bogus letters

  for (int32 i=0; i<slen; i++)
    if (fr->sstr[i] == '.') {
      //fprintf(stderr, "%s %s %d %d\n", fr->sstr, fr->qstr, clrL, clrR);
      return(0);
    }


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
readQSeq(FILE *F, char *N, ilFragment *fr, char end) {

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

  assert(nv == 11);

  sprintf(fr->snam, "@%s:%s:%s:%s:%s#%s/%s", v[0], v[2], v[3], v[4], v[5], v[6], v[7]);
  sprintf(fr->sstr, "%s", v[8]);
  sprintf(fr->qnam, "+%s:%s:%s:%s:%s#%s/%s", v[0], v[2], v[3], v[4], v[5], v[6], v[7]);
  sprintf(fr->qstr, "%s", v[9]);

  return(processSeq(N, fr, end));
}


static
uint64
readSeq(FILE *F, char *N, ilFragment *fr, char end) {

  fgets(fr->snam, AS_READ_MAX_NORMAL_LEN, F);  chomp(fr->snam);
  fgets(fr->sstr, AS_READ_MAX_NORMAL_LEN, F);  chomp(fr->sstr);
  fgets(fr->qnam, AS_READ_MAX_NORMAL_LEN, F);  chomp(fr->qnam);
  fgets(fr->qstr, AS_READ_MAX_NORMAL_LEN, F);  chomp(fr->qstr);

  if (feof(F))
    return(0);

  return(processSeq(N, fr, end));
}


static
void
loadIlluminaReads(char *lname, char *rname, bool isSeq) {
  fprintf(stderr, "Processing illumina reads from '%s'\n", lname);
  fprintf(stderr, "                           and '%s'\n", rname);

  errno = 0;
  FILE  *lfile = fopen(lname, "r");
  if (errno)
    fprintf(stderr, "ERROR:  Coulnd't open Illumina file '%s' for reading: %s\n", lname, strerror(errno)), exit(1);

  FILE  *rfile = fopen(rname, "r");
  if (errno)
    fprintf(stderr, "ERROR:  Coulnd't open Illumina file '%s' for reading: %s\n", rname, strerror(errno)), exit(1);

  ilFragment  *lfrg = new ilFragment;
  ilFragment  *rfrg = new ilFragment;

  lfrg->fr.gkFragment_enableGatekeeperMode(gkpStore);
  rfrg->fr.gkFragment_enableGatekeeperMode(gkpStore);

  while (!feof(lfile) && !feof(rfile)) {
    uint32 nfrg = gkpStore->gkStore_getNumFragments();

    uint64 lUID = (isSeq) ? readSeq(lfile, lname, lfrg, 'l') : readQSeq(lfile, lname, lfrg, 'l');
    uint64 rUID = (isSeq) ? readSeq(rfile, rname, rfrg, 'r') : readQSeq(rfile, rname, rfrg, 'r');

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
}



static
void
loadIlluminaReads(char *fname, bool isSeq) {
  fprintf(stderr, "Processing illumina reads from '%s'.\n", fname);


}



void
checkLibraryForIlluminaPointers(LibraryMesg *lib_mesg) {

  for (uint32 i=0; i<lib_mesg->num_features; i++) {
    if (strcasecmp(lib_mesg->features[i], "illuminaQSequence") == 0) {
      char *sl = lib_mesg->values[i];
      char *sr = lib_mesg->values[i];

      while ((*sr) && (*sr != ','))
        sr++;

      if (*sr) {
        *sr = 0;
        sr++;
        loadIlluminaReads(sl, sr, false);
      } else {
        loadIlluminaReads(sl, false);
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
        loadIlluminaReads(sl, sr, true);
      } else {
        loadIlluminaReads(sl, true);
      }
    }
  }
}
