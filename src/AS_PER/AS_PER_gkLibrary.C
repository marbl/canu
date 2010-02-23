
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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

static const char *rcsid = "$Id: AS_PER_gkLibrary.C,v 1.7 2010-02-23 03:49:28 brianwalenz Exp $";

#include "AS_PER_gkpStore.h"


static
int
decodeBoolean(char *feature, char *value) {
  int  ret = 0;

  //  Decodes a string with 0/1, false/true, no/yes into an integer flag.

  switch (value[0]) {
    case '0':
    case 'f':
    case 'F':
    case 'n':
    case 'N':
      ret = 0;
      break;
    case '1':
    case 't':
    case 'T':
    case 'y':
    case 'Y':
      ret = 1;
      break;
    default:
      fprintf(stderr, "AS_PER_decodeLibraryFeatures()-- Found feature '%s' but has unknown boolean value '%s'\n",
              feature, value);
      break;
  }

  return(ret);
}



void
gkLibrary::gkLibrary_decodeFeatures(LibraryMesg *lmesg) {
  int f;
  for (f=0; f<lmesg->num_features; f++) {
    char *fea = lmesg->features[f];
    char *val = lmesg->values[f];

    //  Unitigger options
    if      (strcasecmp(fea, "forceBOGunitigger") == 0)
      forceBOGunitigger = decodeBoolean("forceBOGunitigger", val);

    else if (strcasecmp(fea, "isNotRandom") == 0)
      isNotRandom = decodeBoolean("isNotRandom", val);

    //  Alignment options
    else if (strcasecmp(fea, "doNotTrustHomopolymerRuns") == 0)
      doNotTrustHomopolymerRuns = decodeBoolean("doNotTrustHomopolymerRuns", val);

    //  OBT options
    else if (strcasecmp(fea, "doRemoveDuplicateReads") == 0)
      doRemoveDuplicateReads = decodeBoolean("doRemoveDuplicateReads", val);

    else if (strcasecmp(fea, "doNotQVTrim") == 0)
      doNotQVTrim = decodeBoolean("doNotQVTrim", val);

    else if (strcasecmp(fea, "goodBadQVThreshold") == 0)
      goodBadQVThreshold = strtoul(val, NULL, 10);

    else if (strcasecmp(fea, "doNotOverlapTrim") == 0)
      doNotOverlapTrim = decodeBoolean("doNotOverlapTrim", val);

    //  Gatekeeper options
    else if (strcasecmp(fea, "usePackedFragments") == 0)
      usePackedFragments = decodeBoolean("usePackedFragments", val);

    //  Illumina options, just to make it not complain about unknown features

    else if (strcasecmp(fea, "illuminaFastQType") == 0)
      ;
    else if (strcasecmp(fea, "illuminaQSequence") == 0)
      ;
    else if (strcasecmp(fea, "illuminaSequence") == 0)
      ;
    //  Library options (orientation is not a feature, it's part of the library)

    else
      fprintf(stderr, "gkLibrary_decodeFeatures()-- found feature '%s' but don't understand it.\n",
              fea);
  }
}


void
gkLibrary::gkLibrary_encodeFeaturesCleanup(LibraryMesg *lmesg) {
  while (lmesg->num_features > 0) {
    lmesg->num_features--;
    safe_free(lmesg->features[lmesg->num_features]);
    safe_free(lmesg->values  [lmesg->num_features]);
  }
  safe_free(lmesg->features);
  safe_free(lmesg->values);
}


void
gkLibrary::gkLibrary_encodeFeatures(LibraryMesg *lmesg) {

  //  Examine the gkpl, allocate space to encode the features into
  //  features/values, return the number of features encoded.
  //
  //  Be sure to call gkLibrary_encodeFeaturesCleanup to properly
  //  cleanup the LibraryMesg after it is written!

  //  We can hardcode the maximum number of features we expect to be
  //  writing.  Otherwise, we should count the number of features we
  //  want to encode, allocate....but what a pain.  We'll just assert
  //  if there are too many.
  //
  lmesg->num_features = 0;
  lmesg->features     = (char **)safe_malloc(64 * sizeof(char*));
  lmesg->values       = (char **)safe_malloc(64 * sizeof(char*));

  int    nf  = 0;
  char **fea = lmesg->features;
  char **val = lmesg->values;

  //  Mostly for debugging, but just might be generally a
  //  GoodThing(tm) to always specify optional features.
  int    alwaysEncode = 1;

  //  Unitigger options
  if (forceBOGunitigger || alwaysEncode) {
    fea[nf] = (char *)safe_malloc(32 * sizeof(char));
    val[nf] = (char *)safe_malloc(32 * sizeof(char));
    sprintf(fea[nf], "forceBOGunitigger");
    sprintf(val[nf], "%d", forceBOGunitigger);
    nf++;
  }

  if (isNotRandom || alwaysEncode) {
    fea[nf] = (char *)safe_malloc(32 * sizeof(char));
    val[nf] = (char *)safe_malloc(32 * sizeof(char));
    sprintf(fea[nf], "isNotRandom");
    sprintf(val[nf], "%d", isNotRandom);
    nf++;
  }

  //  Alignment options
  if (doNotTrustHomopolymerRuns || alwaysEncode) {
    fea[nf] = (char *)safe_malloc(32 * sizeof(char));
    val[nf] = (char *)safe_malloc(32 * sizeof(char));
    sprintf(fea[nf], "doNotTrustHomopolymerRuns");
    sprintf(val[nf], "%d", doNotTrustHomopolymerRuns);
    nf++;
  }

  //  OBT options
  if (doRemoveDuplicateReads || alwaysEncode) {
    fea[nf] = (char *)safe_malloc(32 * sizeof(char));
    val[nf] = (char *)safe_malloc(32 * sizeof(char));
    sprintf(fea[nf], "doRemoveDuplicateReads");
    sprintf(val[nf], "%d", doRemoveDuplicateReads);
    nf++;
  }

  if (doNotQVTrim || alwaysEncode) {
    fea[nf] = (char *)safe_malloc(32 * sizeof(char));
    val[nf] = (char *)safe_malloc(32 * sizeof(char));
    sprintf(fea[nf], "doNotQVTrim");
    sprintf(val[nf], "%d", doNotQVTrim);
    nf++;
  }

  if (goodBadQVThreshold || alwaysEncode) {
    fea[nf] = (char *)safe_malloc(32 * sizeof(char));
    val[nf] = (char *)safe_malloc(32 * sizeof(char));
    sprintf(fea[nf], "goodBadQVThreshold");
    sprintf(val[nf], "%d", goodBadQVThreshold);
    nf++;
  }

  if (doNotOverlapTrim || alwaysEncode) {
    fea[nf] = (char *)safe_malloc(32 * sizeof(char));
    val[nf] = (char *)safe_malloc(32 * sizeof(char));
    sprintf(fea[nf], "doNotOverlapTrim");
    sprintf(val[nf], "%d", doNotOverlapTrim);
    nf++;
  }

  //  Gatekeeper options
  if (usePackedFragments || alwaysEncode) {
    fea[nf] = (char *)safe_malloc(32 * sizeof(char));
    val[nf] = (char *)safe_malloc(32 * sizeof(char));
    sprintf(fea[nf], "usePackedFragments");
    sprintf(val[nf], "%d", usePackedFragments);
    nf++;
  }

  //  Library options (orientation is not a feature, it's part of the library)

  lmesg->num_features = nf;

  assert(nf < 64);
}
