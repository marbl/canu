
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2010, J. CRaig Venter Institute.
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

const char *mainid = "$Id: fastqToCA.C,v 1.1 2010-01-23 04:05:27 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
//#include <assert.h>
//#include <ctype.h>
//#include <sys/stat.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_PER_gkpStore.h"
#include "AS_MSG_pmesg.h"

int
main(int argc, char **argv) {
  int       insertSize       = 0;
  int       insertStdDev     = 0;
  char     *libraryName      = 0L;

  bool      isMated          = false;

  char    **fastq            = new char * [argc];
  int32     fastqLen         = 0;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-insertsize") == 0) {
      insertSize   = atoi(argv[++arg]);
      insertStdDev = atoi(argv[++arg]);
      isMated      = true;

    } else if (strcmp(argv[arg], "-libraryname") == 0) {
      libraryName = argv[++arg];

    } else if (strcmp(argv[arg], "-fastq") == 0) {
      fastq[fastqLen++] = argv[++arg];

    } else {
      err++;
    }

    arg++;
  }

  if ((insertSize == 0) && (insertStdDev >  0))
    err++;
  if ((insertSize >  0) && (insertStdDev == 0))
    err++;
  if (libraryName == 0L)
    err++;
  if (fastq == 0L)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s [-insertsize <mean> <stddev>] [-libraryname <name>] [-output <name>]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -insertsize i d        Mates are on average i +- d bp apart.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -libraryname n         The UID of the library these reads are added to.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -fastq A               Single ended reads, in fastq format.\n");
    fprintf(stderr, "  -fastq A,B             Paired end reads, in fastq format.\n");
    fprintf(stderr, "\n");

    if ((insertSize == 0) && (insertStdDev >  0))
      fprintf(stderr, "ERROR:  Invalid -insertsize.  Both mean and std.dev must be greater than zero.\n");
    if ((insertSize >  0) && (insertStdDev == 0))
      fprintf(stderr, "ERROR:  Invalid -insertsize.  Both mean and std.dev must be greater than zero.\n");
    if (libraryName == 0L)
      fprintf(stderr, "ERROR:  No library name supplied with -libraryname.\n");
    if (fastq == 0L)
      fprintf(stderr, "ERROR:  No reads supplied with -fastq.\n");

    exit(1);
  }

  //  Check that all the reads are mated or unmated, as required by the insert size description.

  for (int32 i=0; i<fastqLen; i++) {
    int32 ncomma = 0;

    for (int32 j=0; fastq[i][j]; j++)
      if (fastq[i][j] == ',')
        ncomma++;

    if ((isMated == true)  && (ncomma != 1))
      fprintf(stderr, "ERROR:  Library is mated, but -fastq '%s' doesn't supply exactly two files.\n", fastq[i]), err++;
    if ((isMated == false) && (ncomma != 0))
      fprintf(stderr, "ERROR:  Library is unmated, but -fastq '%s' doesn't supply exactly one file.\n", fastq[i]), err++;
  }

  if (err)
    exit(1);

  //  Construct the library

  gkLibrary         gkl;

  gkl.libraryUID                 = AS_UID_load(libraryName);

  gkl.mean                       = insertSize;
  gkl.stddev                     = insertStdDev;

  gkl.forceBOGunitigger          = 1;
  gkl.isNotRandom                = 0;

  gkl.doNotTrustHomopolymerRuns  = 0;

  gkl.doRemoveDuplicateReads     = 0;
  gkl.doNotQVTrim                = 0;
  gkl.goodBadQVThreshold         = 0;
  gkl.doNotOverlapTrim           = 0;

  gkl.usePackedFragments         = 1;

  gkl.orientation                = (isMated) ? AS_READ_ORIENT_INNIE : AS_READ_ORIENT_UNKNOWN;

  //  Construct the messages.

  VersionMesg       vr2Mesg;
  VersionMesg       vr1Mesg;
  LibraryMesg       libMesg;

  vr2Mesg.version      = 2;

  vr1Mesg.version      = 1;

  libMesg.action       = AS_ADD;
  libMesg.eaccession   = gkl.libraryUID;
  libMesg.mean         = gkl.mean;
  libMesg.stddev       = gkl.stddev;
  libMesg.source       = NULL;
#warning unsafe conversion of orient
  libMesg.link_orient = (OrientType)AS_READ_ORIENT_NAMES[gkl.orientation][0];

  gkl.gkLibrary_encodeFeatures(&libMesg);

  //  Add in the non-standard Illumina features.

  libMesg.features = (char **)safe_realloc(libMesg.features, sizeof(char *) * (libMesg.num_features + fastqLen));
  libMesg.values   = (char **)safe_realloc(libMesg.values,   sizeof(char *) * (libMesg.num_features + fastqLen));

  for (int32 i=0; i<fastqLen; i++) {
    int32 nf = libMesg.num_features;

    libMesg.features[nf] = (char *)safe_malloc(sizeof(char) * 32);
    libMesg.values  [nf] = (char *)safe_malloc(sizeof(char) * (strlen(fastq[i]) + 1));

    sprintf(libMesg.features[nf], "illuminaSequence");
    sprintf(libMesg.values[nf],   "%s", fastq[i]);

    libMesg.num_features++;
  }

  //  Emit the VER and LIB messages.  Enable version 2, write the LIB, switch back to version 1.

  GenericMesg       pmesg;

  pmesg.m = &vr2Mesg;
  pmesg.t = MESG_VER;
  WriteProtoMesg_AS(stdout, &pmesg);

  pmesg.m = &libMesg;
  pmesg.t = MESG_LIB;
  WriteProtoMesg_AS(stdout, &pmesg);

  pmesg.m = &vr1Mesg;
  pmesg.t = MESG_VER;
  WriteProtoMesg_AS(stdout, &pmesg);

  //  Clean up.

  gkl.gkLibrary_encodeFeaturesCleanup(&libMesg);

  delete [] fastq;

  exit(0);
}

