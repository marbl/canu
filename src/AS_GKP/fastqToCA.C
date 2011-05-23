
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

const char *mainid = "$Id: fastqToCA.C,v 1.13 2011-05-23 19:04:55 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_PER_gkpStore.h"
#include "AS_MSG_pmesg.h"



void
addFeature(LibraryMesg *libMesg, char *feature, char *value) {
  int32 nf = libMesg->num_features;

  libMesg->features[nf] = (char *)safe_malloc(sizeof(char) * (strlen(feature) + 1));
  libMesg->values  [nf] = (char *)safe_malloc(sizeof(char) * (strlen(value) + 1));

  strcpy(libMesg->features[nf], feature);
  strcpy(libMesg->values[nf],   value);

  libMesg->num_features++;
}



int
main(int argc, char **argv) {
  int       insertSize       = 0;
  int       insertStdDev     = 0;
  char     *libraryName      = 0L;

  bool      isMated          = false;

  char     *type             = "illumina";

  char     *orientInnie      = "innie";
  char     *orientOuttie     = "outtie";
  char     *orient           = orientInnie;

  bool      interlaced       = false;

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

    } else if (strcmp(argv[arg], "-type") == 0) {
      type = argv[++arg];

    } else if (strcmp(argv[arg], "-innie") == 0) {
      orient = orientInnie;

    } else if (strcmp(argv[arg], "-outtie") == 0) {
      orient = orientOuttie;

    } else if (strcmp(argv[arg], "-interlaced") == 0) {
      interlaced = true;

    } else if (strcmp(argv[arg], "-fastq") == 0) {
      fastq[fastqLen++] = argv[++arg];

    } else {
      fprintf(stderr, "ERROR:  Unknown option '%s'\n", argv[arg]);
      exit(1);
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
  if ((strcasecmp(type, "sanger") != 0) && (strcasecmp(type, "solexa") != 0) && (strcasecmp(type, "illumina") != 0))
    err++;
  if (fastqLen == 0)
    err++;
  if ((isMated == false) && (interlaced == true))
    err++;

  if (err) {
    fprintf(stderr, "usage: %s [-insertsize <mean> <stddev>] [-libraryname <name>]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Writes a CA FRG file consisting of a single LIB message to stdout.  When this file is fed to\n");
    fprintf(stderr, "gatekeeper, gatekeeper will read fragments from the original fastq files supplied with the\n");
    fprintf(stderr, "-fastq option.  The full absolute path name must be used for fastq files.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -insertsize i d    Mates are on average i +- d bp apart.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -libraryname n     The UID of the library these reads are added to.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -type t            What type of fastq ('illumina' is the default):\n");
    fprintf(stderr, "                       'sanger'   -- QV's are PHRED, offset=33 '!', NCBI SRA data.\n");
    fprintf(stderr, "                       'solexa'   -- QV's are Solexa, early Solexa data.\n");
    fprintf(stderr, "                       'illumina' -- QV's are PHRED, offset=64 '@', Illumina reads from version 1.3 on.\n");
    fprintf(stderr, "                     See Cock, et al., 'The Sanger FASTQ file format for sequences with quality scores, and\n");
    fprintf(stderr, "                     the Solexa/Illumina FASTQ variants', doi:10.1093/nar/gkp1137\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -innie             The paired end reads are 5'-3' <-> 3'-5' (usually for paired-end reads) (default)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -outtie            The paired end reads are 3'-5' <-> 5'-3' (usually for mate-pair reads)\n");
    fprintf(stderr, "                     This switch will reverse-complement every read, transforming outtie-oriented\n");
    fprintf(stderr, "                     mates into innie-oriented mates.  This trick only works if all reads are the\n");
    fprintf(stderr, "                     same length.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -interlaced        The paired reads come one after another in a single fastq file.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -fastq A           Single ended or interlaced paired end reads, in fastq format.\n");
    fprintf(stderr, "  -fastq A,B         Paired end reads, in fastq format.\n");
    fprintf(stderr, "\n");

    if ((insertSize == 0) && (insertStdDev >  0))
      fprintf(stderr, "ERROR:  Invalid -insertsize.  Both mean and std.dev must be greater than zero.\n");
    if ((insertSize >  0) && (insertStdDev == 0))
      fprintf(stderr, "ERROR:  Invalid -insertsize.  Both mean and std.dev must be greater than zero.\n");
    if (libraryName == 0L)
      fprintf(stderr, "ERROR:  No library name supplied with -libraryname.\n");
    if ((strcasecmp(type, "sanger") != 0) && (strcasecmp(type, "solexa") != 0) && (strcasecmp(type, "illumina") != 0))
      fprintf(stderr, "ERROR:  Invalid type '%s' supplied with -type.\n", type);
    if (fastqLen == 0)
      fprintf(stderr, "ERROR:  No reads supplied with -fastq.\n");
    if ((isMated == false) && (interlaced == true))
      fprintf(stderr, "ERROR:  Interlaced reads (-interlaced) must have am insert size (-insertsize).\n");

    exit(1);
  }

  //  Check that all the reads are mated or unmated, as required by the insert size description.

  for (int32 i=0; i<fastqLen; i++) {
    int32 ncomma = 0;

    for (int32 j=0; fastq[i][j]; j++)
      if (fastq[i][j] == ',')
        ncomma++;

    if ((isMated == true)  && (interlaced == false) && (ncomma != 1))
      fprintf(stderr, "ERROR:  Library is mated, but -fastq '%s' doesn't supply exactly two files.\n", fastq[i]), err++;
    if ((isMated == false) && (ncomma != 0))
      fprintf(stderr, "ERROR:  Library is unmated, but -fastq '%s' doesn't supply exactly one file.\n", fastq[i]), err++;
    if ((interlaced == true) && (ncomma != 0))
      fprintf(stderr, "ERROR:  Library is interlaced, but -fastq '%s' doesn't supply exactly one file.\n", fastq[i]), err++;
  }

  if (err)
    exit(1);

  //  Check that all the read files exist, and that the paths are absolute.

  int32   fastqPathErrors = 0;

  for (int32 i=0; i<fastqLen; i++) {
    char   *f1  = fastq[i];
    char   *f2  = strrchr(fastq[i], ',');
    char    cwd[FILENAME_MAX];

    getcwd(cwd, FILENAME_MAX);

    if (f2) {
      *f2 = 0;
      f2++;
    }

    if ((f1) && (AS_UTL_fileExists(f1, FALSE, FALSE) == false))
      fprintf(stderr, "ERROR: fastq file '%s' doesn't exist.\n", f1), fastqPathErrors++;

    if ((f2) && (AS_UTL_fileExists(f2, FALSE, FALSE) == false))
      fprintf(stderr, "ERROR: fastq file '%s' doesn't exist.\n", f2), fastqPathErrors++;

    if ((f1) && (f1[0] != '/')) {
      char *n1 = new char [FILENAME_MAX];
      sprintf(n1, "%s/%s", cwd, f1);
      f1 = n1;

      if (AS_UTL_fileExists(f1, FALSE, FALSE) == false)
        fprintf(stderr, "ERROR: absolute-path fastq file '%s' doesn't exist.\n", f1), fastqPathErrors++;
    }

    if ((f2) && (f2[0] != '/')) {
      char *n2 = new char [FILENAME_MAX];
      sprintf(n2, "%s/%s", cwd, f2);
      f2 = n2;

      if (AS_UTL_fileExists(f2, FALSE, FALSE) == false)
        fprintf(stderr, "ERROR: absolute-path fastq file '%s' doesn't exist.\n", f2), fastqPathErrors++;
    }

    char *n = new char [FILENAME_MAX + FILENAME_MAX];

    if (f2)
      sprintf(n, "%s,%s", f1, f2);
    else
      sprintf(n, "%s", f1);

    fastq[i] = n;
  }

  if (fastqPathErrors)
    fprintf(stderr, "ERROR: some fastq files not found.\n"), exit(1);

  //  Construct the library

  gkLibrary         gkl;

  gkl.libraryUID                 = AS_UID_load(libraryName);

  gkl.mean                       = insertSize;
  gkl.stddev                     = insertStdDev;

  gkl.forceBOGunitigger          = 1;
  gkl.isNotRandom                = 0;

  gkl.doNotTrustHomopolymerRuns  = 0;

  gkl.doMerBasedTrimming         = 1;
  gkl.doRemoveDuplicateReads     = 1;
  gkl.doNotQVTrim                = 0;
  gkl.goodBadQVThreshold         = 0;
  gkl.doNotOverlapTrim           = 0;

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

  libMesg.link_orient.setIsUnknown();

  switch(gkl.orientation) {
    case AS_READ_ORIENT_INNIE:
      libMesg.link_orient.setIsInnie();
      break;
    case AS_READ_ORIENT_OUTTIE:
      libMesg.link_orient.setIsOuttie();
      break;
    case AS_READ_ORIENT_NORMAL:
      libMesg.link_orient.setIsNormal();
      break;
    case AS_READ_ORIENT_ANTINORMAL:
      libMesg.link_orient.setIsAnti();
      break;
    case AS_READ_ORIENT_UNKNOWN:
      libMesg.link_orient.setIsUnknown();
      break;
    default:
      //  Cannot happen, unless someone adds a new orientation to gkFragment.
      assert(0);
      break;
  }

  gkl.gkLibrary_encodeFeatures(&libMesg);

  //  Add in the non-standard Illumina features.

  libMesg.features = (char **)safe_realloc(libMesg.features, sizeof(char *) * (libMesg.num_features + 2 + fastqLen));
  libMesg.values   = (char **)safe_realloc(libMesg.values,   sizeof(char *) * (libMesg.num_features + 2 + fastqLen));

  addFeature(&libMesg, "illuminaFastQType", type);
  addFeature(&libMesg, "illuminaOrientation", orient);

  for (int32 i=0; i<fastqLen; i++)
    addFeature(&libMesg, "illuminaSequence", fastq[i]);

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

