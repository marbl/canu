
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

const char *mainid = "$Id: gkpStoreDumpFASTQ.C,v 1.1 2012-02-06 08:30:19 brianwalenz Exp $";

#include "AS_global.H"
#include "AS_PER_gkpStore.H"
#include "AS_UTL_fasta.H"

  //  Currently doesn't support a single library, or figuring out which
  //  libraries are not in the dump range.


class libInfo {
public:
  libInfo(char *outPrefix, char *libName) {
    sprintf(aname, "%s.%s.1.fastq",       outPrefix, libName);
    sprintf(bname, "%s.%s.2.fastq",       outPrefix, libName);
    sprintf(pname, "%s.%s.paired.fastq",  outPrefix, libName);
    sprintf(uname, "%s.%s.unmated.fastq", outPrefix, libName);

    errno = 0;

    a = fopen(aname, "w");
    if (errno)
      fprintf(stderr, "Failed to open output file '%s': %s\n", aname, strerror(errno)), exit(1);

    b = fopen(bname, "w");
    if (errno)
      fprintf(stderr, "Failed to open output file '%s': %s\n", bname, strerror(errno)), exit(1);

    p = fopen(pname, "w");
    if (errno)
      fprintf(stderr, "Failed to open output file '%s': %s\n", pname, strerror(errno)), exit(1);

    u = fopen(uname, "w");
    if (errno)
      fprintf(stderr, "Failed to open output file '%s': %s\n", uname, strerror(errno)), exit(1);
  };
  ~libInfo() {
    fclose(a);
    fclose(b);
    fclose(p);
    fclose(u);
  };

  char aname[FILENAME_MAX];
  char bname[FILENAME_MAX];
  char pname[FILENAME_MAX];
  char uname[FILENAME_MAX];

  FILE *a;
  FILE *b;
  FILE *p;
  FILE *u;
};    


int
main(int argc, char **argv) {
  char            *gkpStoreName      = NULL;
  char            *outPrefix         = NULL;

  AS_IID           libToDump         = 0;
  uint32           clrToDump         = AS_READ_CLEAR_LATEST;

  AS_IID           bgnIID            = 1;
  AS_IID           endIID            = AS_IID_MAX;

  bool             dumpAllBases      = true;
  bool             dumpAllReads      = false;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-l") == 0) {
      libToDump = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-b") == 0) {
      bgnIID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      endIID  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-c") == 0) {
      clrToDump = gkStore_decodeClearRegionLabel(argv[++arg]);

    } else if (strcmp(argv[arg], "-g") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      outPrefix = argv[++arg];

    } else {
      err++;
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if (gkpStoreName == NULL)
    err++;
  if (outPrefix == NULL)
    err++;
  if (clrToDump == AS_READ_CLEAR_ERROR)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s [...] -o fastq-prefix -g gkpStore\n", argv[0]);
    fprintf(stderr, "  -g gkpStore\n");
    fprintf(stderr, "  -o fastq-prefix     write files fastq-prefix.1.fastq, fastq-prefix.2.fastq, fastq-prefix.paired.fastq, fastq-prefix.unmated.fastq\n");
    fprintf(stderr, "  \n");
    fprintf(stderr, "  -l libToDump        output only fragments in library number libToDump (NOT IMPLEMENTED)\n");
    fprintf(stderr, "  -b iid              output starting at fragment iid\n");
    fprintf(stderr, "  -e iid              output stopping after fragment iid\n");
    fprintf(stderr, "  -c clrName          output clear range 'clrName'\n");
    fprintf(stderr, "  \n");

    if (gkpStoreName == NULL)
      fprintf(stderr, "ERROR: no gkpStore (-g) supplied.\n");
    if (outPrefix == NULL)
      fprintf(stderr, "ERROR: no output prefix (-o) supplied.\n");
    if (clrToDump == AS_READ_CLEAR_ERROR)
      fprintf(stderr, "ERROR: clear range (-c) is not a valid clear range.\n");
    exit(1);
  }

  gkStore    *gkp       = new gkStore(gkpStoreName, FALSE, FALSE);

  AS_IID    numFrags    = gkp->gkStore_getNumFragments();
  AS_IID    numLibs     = gkp->gkStore_getNumLibraries();

  libInfo **lib         = new libInfo * [numLibs];

  lib[0] = new libInfo(outPrefix, "legacy");

  for (uint32 i=1; i<numLibs; i++)
    lib[i] = new libInfo(outPrefix, gkp->gkStore_getLibrary(i)->libraryName);

  if (bgnIID < 1)
    bgnIID = 1;
  if (numFrags < endIID)
    endIID = numFrags;

  //AS_IID    streamBgn = AS_IID_MIN;
  //AS_IID    streamEnd = AS_IID_MAX;

  gkStream   *fs        = new gkStream(gkp, bgnIID, endIID, GKFRAGMENT_QLT);
  gkFragment  fr;



  while (fs->next(&fr)) {
    int32   lclr   = fr.gkFragment_getClearRegionBegin(clrToDump);
    int32   rclr   = fr.gkFragment_getClearRegionEnd  (clrToDump);

    AS_IID  id1    = fr.gkFragment_getReadIID();
    AS_IID  id2    = fr.gkFragment_getMateIID();

    AS_IID  libIID = fr.gkFragment_getLibraryIID();

    if ((dumpAllReads == false) && (fr.gkFragment_getIsDeleted() == true))
      //  Fragment is deleted, don't dump.
      continue;

    if ((libToDump != 0) && (fr.gkFragment_getLibraryIID() == libToDump))
      //  Fragment isn't marked for dumping, don't dump.
      continue;

    if ((dumpAllBases == false) && (lclr >= rclr))
      //  Fragment has null or invalid clear range, don't dump.
      continue;

    if ((id2 != 0) && (id2 < id1))
      //  Mated, and the mate is the first frag.  We've already reported this one.
      continue;

    char *seq = fr.gkFragment_getSequence() + ((dumpAllBases == false) ? fr.gkFragment_getClearRegionBegin(clrToDump) : 0);
    char *qlt = fr.gkFragment_getQuality()  + ((dumpAllBases == false) ? fr.gkFragment_getClearRegionBegin(clrToDump) : 0);

    int32 len = (dumpAllBases == false) ? fr.gkFragment_getClearRegionLength(clrToDump) : fr.gkFragment_getSequenceLength();

    seq[len] = 0;
    qlt[len] = 0;

    if (dumpAllBases == true) {
      for (uint32 i=0; i<lclr; i++)
        seq[i] += (seq[i] >= 'A') ? 'a' - 'A' : 0;

      for (uint32 i=lclr; i<rclr; i++)
        seq[i] += (seq[i] >= 'A') ? 0 : 'A' - 'a';

      for (uint32 i=rclr; seq[i]; i++)
        seq[i] += (seq[i] >= 'A') ? 'a' - 'A' : 0;
    }

    if (id2 == 0) {
      //  Unmated read, dump to the unmated reads file.
      AS_UTL_writeFastQ(lib[libIID]->u, seq, len, qlt, len,
                        "@%s clr="F_U32","F_U32" clv="F_U32","F_U32" max="F_U32","F_U32" tnt="F_U32","F_U32" rnd=%c\n",
                        AS_UID_toString(fr.gkFragment_getReadUID()),
                        fr.gkFragment_getClearRegionBegin(clrToDump),
                        fr.gkFragment_getClearRegionEnd  (clrToDump),
                        fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_VEC),
                        fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_VEC),
                        fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_MAX),
                        fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_MAX),
                        fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_TNT),
                        fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_TNT),
                        fr.gkFragment_getIsNonRandom() ? 'f' : 't');
      continue;
    }

    //  Write the first fragment (twice).
    AS_UTL_writeFastQ(lib[libIID]->a, seq, len, qlt, len,
                      "@%s clr="F_U32","F_U32" clv="F_U32","F_U32" max="F_U32","F_U32" tnt="F_U32","F_U32" rnd=%c\n",
                      AS_UID_toString(fr.gkFragment_getReadUID()),
                      fr.gkFragment_getClearRegionBegin(clrToDump),
                      fr.gkFragment_getClearRegionEnd  (clrToDump),
                      fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_VEC),
                      fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_VEC),
                      fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_MAX),
                      fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_MAX),
                      fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_TNT),
                      fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_TNT),
                      fr.gkFragment_getIsNonRandom() ? 'f' : 't');

    AS_UTL_writeFastQ(lib[libIID]->p, seq, len, qlt, len,
                      "@%s clr="F_U32","F_U32" clv="F_U32","F_U32" max="F_U32","F_U32" tnt="F_U32","F_U32" rnd=%c\n",
                      AS_UID_toString(fr.gkFragment_getReadUID()),
                      fr.gkFragment_getClearRegionBegin(clrToDump),
                      fr.gkFragment_getClearRegionEnd  (clrToDump),
                      fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_VEC),
                      fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_VEC),
                      fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_MAX),
                      fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_MAX),
                      fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_TNT),
                      fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_TNT),
                      fr.gkFragment_getIsNonRandom() ? 'f' : 't');

    //  Grab the second fragment.

    gkp->gkStore_getFragment(id2, &fr, GKFRAGMENT_QLT);

    lclr = fr.gkFragment_getClearRegionBegin(clrToDump) + 1;
    rclr = fr.gkFragment_getClearRegionEnd  (clrToDump);

    seq = fr.gkFragment_getSequence() + ((dumpAllBases == false) ? fr.gkFragment_getClearRegionBegin(clrToDump) : 0);
    qlt = fr.gkFragment_getQuality()  + ((dumpAllBases == false) ? fr.gkFragment_getClearRegionBegin(clrToDump) : 0);
    len = (dumpAllBases == false) ? fr.gkFragment_getClearRegionLength(clrToDump) : fr.gkFragment_getSequenceLength();

    seq[len] = 0;
    qlt[len] = 0;

    //  Write the second fragment (twice).
    AS_UTL_writeFastQ(lib[libIID]->b, seq, len, qlt, len,
                      "@%s clr="F_U32","F_U32" clv="F_U32","F_U32" max="F_U32","F_U32" tnt="F_U32","F_U32" rnd=%c\n",
                      AS_UID_toString(fr.gkFragment_getReadUID()),
                      fr.gkFragment_getClearRegionBegin(clrToDump),
                      fr.gkFragment_getClearRegionEnd  (clrToDump),
                      fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_VEC),
                      fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_VEC),
                      fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_MAX),
                      fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_MAX),
                      fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_TNT),
                      fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_TNT),
                      fr.gkFragment_getIsNonRandom() ? 'f' : 't');

    AS_UTL_writeFastQ(lib[libIID]->p, seq, len, qlt, len,
                      "@%s clr="F_U32","F_U32" clv="F_U32","F_U32" max="F_U32","F_U32" tnt="F_U32","F_U32" rnd=%c\n",
                      AS_UID_toString(fr.gkFragment_getReadUID()),
                      fr.gkFragment_getClearRegionBegin(clrToDump),
                      fr.gkFragment_getClearRegionEnd  (clrToDump),
                      fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_VEC),
                      fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_VEC),
                      fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_MAX),
                      fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_MAX),
                      fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_TNT),
                      fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_TNT),
                      fr.gkFragment_getIsNonRandom() ? 'f' : 't');
  }

  delete fs;
  delete gkp;

  exit(0);
}
