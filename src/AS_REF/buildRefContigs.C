
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 20008, J. Craig Venter Institute.
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

const char *mainid = "$Id: buildRefContigs.C,v 1.2 2008-10-30 04:49:44 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

extern "C" {
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_fileIO.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_gkpStore.h"
}

GateKeeperStore   *gkpStore = 0L;

#include "splitToWords.H"
#include "refAlignment.H"

fragmentData    *frg    = 0L;
uint32           frgLen = 0;

mappingData     *ali    = 0L;
uint32           aliLen = 0;
uint32           aliMax = 0;


static
void
readMapping(char *filename) {
  char  L[1024];
  FILE *F;

  errno = 0;
  F = fopen(filename, "r");
  if (errno)
    fprintf(stderr, "Failed to open input mapping '%s': %s\n", filename, strerror(errno)), exit(1);

  fgets(L, 1024, F);

  while (!feof(F)) {
    if (aliLen >= aliMax) {
      if (aliMax == 0)
        aliMax = 1048576 / 2;
      aliMax *= 2;
      mappingData *A = new mappingData [aliMax];
      memcpy(A, ali, sizeof(mappingData) * aliLen);
      delete [] ali;
      ali = A;
    }

    mappingDataParseLine(ali + aliLen++, L);

    fgets(L, 1024, F);
  }

  fclose(F);

  if (ali == 0L)
    fprintf(stderr, "Failed to read alignments from '%s'.\n", filename);

  qsort(ali, aliLen, sizeof(mappingData), mappingDataSortByRefPosition);

  //  Add sentinel to the end of ali, so that the main processing loop
  //  will find the last contig naturally.

  ali[aliLen].frgIID  = 0;
  ali[aliLen].frgUID  = AS_UID_undefined();
  ali[aliLen].frgBgn  = 0;
  ali[aliLen].frgEnd  = 0;
  ali[aliLen].fwd     = 0;
  ali[aliLen].refUID  = AS_UID_undefined();
  ali[aliLen].refBgn  = 0;
  ali[aliLen].refEnd  = 0;

  aliLen++;

  fprintf(stderr, "Read "F_U32" alignments.\n", aliLen);
}


static
void
loadFragments(void) {
  fragRecord  fr;
  uint32      err=0, alive=0, dead=0, unmapped=0;

  frgLen = getNumGateKeeperFragments(gkpStore) + 1;
  frg    = (fragmentData *)safe_calloc(frgLen, sizeof(fragmentData));

  memset(frg, 0, sizeof(fragmentData) * frgLen);

  //  Load the fragment data.  Set the mappingDataIndex to an invalid
  //  value, we'll set the real value in the loop after this.

  for (uint32 i=1; i<frgLen; i++) {
    getFrag(gkpStore, i, &fr, FRAG_S_INF);

    if (getFragRecordIsDeleted(&fr)) {
      dead++;
      frg[i].fragUID    = AS_UID_undefined();
      frg[i].mateIID    = 0;
      frg[i].libraryIID = 0;
    } else {
      alive++;
      frg[i].fragUID    = getFragRecordUID(&fr);
      frg[i].mateIID    = getFragRecordMateIID(&fr);
      frg[i].libraryIID = getFragRecordLibraryIID(&fr);
    }

    frg[i].mappingDataIndex = aliLen;
  }

  //  For each mapping, get the read IID, and update the
  //  mappingDataIndex pointer from the fragment data to the mapping
  //  data.
  //
  //  aliLen-1 since the last thing is a sentinel end of list marker,
  //  not a real alignment.

  unmapped = alive;

  for (uint32 i=0; i<aliLen-1; i++) {
    uint32  iid = ali[i].frgIID;

    if (AS_UID_compare(frg[iid].fragUID, AS_UID_undefined()) == 0)
      fprintf(stderr, "ERROR:  Fragment %s,%d is deleted in gkpStore, but used in the mapping.\n",
              AS_UID_toString(frg[iid].fragUID), iid), err++;

    if (frg[iid].mappingDataIndex != aliLen)
      fprintf(stderr, "ERROR:  Fragment %s,%d appears more than once in the mapping.\n",
              AS_UID_toString(frg[iid].fragUID), iid), err++;

    unmapped--;

    frg[iid].mappingDataIndex = i;
  }

  if (err)
    fprintf(stderr, "There were errors in the mapping.  Fail.\n"), exit(1);

  fprintf(stderr, "Found "F_U32" alive, "F_U32" dead, and "F_U32" unmapped fragments.\n", alive, dead, unmapped);
}


static
void
outputDistances(void) {
  int               imdLen = getNumGateKeeperLibraries(gkpStore)+1;
  IntMateDistMesg  *imd = (IntMateDistMesg *)safe_calloc(imdLen, sizeof(IntMateDistMesg));

  for (int i=1; i<imdLen; i++) {
    imd[i].refines     = i;
    imd[i].mean        = 0.0;
    imd[i].stddev      = 0.0;
    imd[i].min         = 0;
    imd[i].max         = 0;
    imd[i].num_buckets = 0;
    imd[i].histogram   = NULL;
  }

  //  Count how many samples we'll have
  //
  for (int i=1; i<frgLen; i++)
    imd[frg[i].libraryIID].num_buckets++;

  //  Allocate space for the histograms
  //
  for (int i=1; i<imdLen; i++) {
    imd[i].histogram   = (int32 *)safe_malloc(sizeof(int32) * imd[i].num_buckets);
    imd[i].num_buckets = 0;
  }

  //  Populate the histogram
  //
  for (int i=1; i<frgLen; i++) {
    int  l  = frg[i].libraryIID;
    int  m  = frg[i].mateIID;
    int  ii = frg[i].mappingDataIndex;
    int  mm = frg[m].mappingDataIndex;

    if ((m == 0) || (i < m))
      //  No mate, or mate is later in the list
      continue;

    if (ali[ii].fwd == ali[mm].fwd)
      //  Same orientation.
      continue;

    if (AS_UID_compare(ali[ii].refUID, ali[mm].refUID) != 0)
      //  On different reference sequences
      continue;

    if ((ali[ii].refBgn < ali[mm].refBgn) && (ali[ii].fwd == 0))
      //  Not innie.  Frag i is before frag m, but pointing the wrong way.
      continue;

    if ((ali[ii].refBgn > ali[mm].refBgn) && (ali[ii].fwd == 1))
      //  Not innie.  Frag i is after frag m, but pointing the wrong way.
      continue;

    //  So, frag i has an opposite oriented mate on the same scaffold,
    //  and the pair is innie.
    //
    //  ii--->   <---mm  OR  mm--->   <---ii
    //
    imd[l].histogram[imd[l].num_buckets++] = ((ali[ii].fwd) ?
                                              (ali[mm].refEnd - ali[ii].refBgn) :
                                              (ali[ii].refEnd - ali[mm].refBgn));

  }

  //  Compute statistics
  //
  for (int i=1; i<imdLen; i++) {
    double  mu = 0;
    double  si = 0;

    imd[i].min = imd[i].histogram[0];
    imd[i].max = imd[i].histogram[0];

    for (int j=0; j<imd[i].num_buckets; j++) {
      if (imd[i].histogram[i] < imd[i].min)
        imd[i].min = imd[i].histogram[i];
      if (imd[i].max < imd[i].histogram[i])
        imd[i].max = imd[i].histogram[i];

      mu += imd[i].histogram[j];
      si += imd[i].histogram[j] * imd[i].histogram[j];
    }

    fprintf(stderr, "mu=%f si=%f\n", mu, si);

#warning stddev computation is bogus
    imd[i].mean   = mu / imd[i].num_buckets;
    imd[i].stddev = sqrt((si - mu * mu / imd[i].num_buckets) / (imd[i].num_buckets - 1));

    safe_free(imd[i].histogram);
    imd[i].num_buckets = 0;
    imd[i].histogram   = NULL;
  }

  //  Output
  //
  for (int i=1; i<imdLen; i++) {
    GenericMesg  pmesg;
    pmesg.m = imd + i;
    pmesg.t = MESG_IMD;
    WriteProtoMesg_AS(stdout, &pmesg);
  }

  safe_free(imd);
}


static
void
outputFragments(void) {
  IntAugFragMesg     iaf;
  GenericMesg        pmesg;
  uint32             num = 0;

  int i, m;

  for (i=1; i<frgLen; i++) {
    if (AS_UID_compare(frg[i].fragUID, AS_UID_undefined()) != 0) {
      m = frg[i].mateIID;

      iaf.iaccession    = i;
      iaf.type          = AS_READ;
      iaf.chaff         = 0;
      iaf.clear_rng.bgn = -1;
      iaf.clear_rng.end = -1;
      iaf.mate_status   = UNASSIGNED_MATE;

      pmesg.m = &iaf;
      pmesg.t = MESG_IAF;
      WriteProtoMesg_AS(stdout, &pmesg);

      num++;
    }
  }

  fprintf(stderr, "Wrote "F_U32" fragment messages.\n", num);
}


static
void
outputMates(void) {
  IntAugMatePairMesg iam;
  GenericMesg        pmesg;
  uint32             num = 0;

  int i, m;

  for (i=1; i<frgLen; i++) {
    m = frg[i].mateIID;

    if ((m != 0) && (m < i)) {
      iam.fragment1     = i;
      iam.fragment2     = m;
      iam.mate_status   = UNASSIGNED_MATE;

      assert(AS_UID_compare(frg[i].fragUID, AS_UID_undefined()) != 0);
      assert(AS_UID_compare(frg[m].fragUID, AS_UID_undefined()) != 0);

      pmesg.m = &iam;
      pmesg.t = MESG_IAM;
      WriteProtoMesg_AS(stdout, &pmesg);

      num++;
    }
  }

  fprintf(stderr, "Wrote "F_U32" mate pair messages.\n", num);
}


static
void
outputUnitigs(void) {
  GenericMesg     pmesg;
  IntUnitigMesg   ium;
  int             impMax = 128 * 1024;
  int             impLen = 0;
  IntMultiPos    *imp    = (IntMultiPos *)safe_malloc(impMax * sizeof(IntMultiPos));
  uint32          num    = 0;

  uint32 iumId = 1;

  uint32 i = 0;
  uint32 b = ali[0].refBgn;
  uint32 e = ali[0].refEnd;
  AS_UID r = ali[0].refUID;

  for (i=0; i<aliLen; i++) {
    //fprintf(stdout, "%s "F_U32"-"F_U32"\n", AS_UID_toString(ali[i].refID), ali[i].refBgn, ali[i].refEnd);

    if ((AS_UID_compare(r, ali[i].refUID) != 0) || (e < ali[i].refBgn)) {
      //fprintf(stderr, "Unitig: "F_U32" %s "F_U32"-"F_U32"\n", i, AS_UID_toString(r), b, e);

      ium.iaccession    = iumId++;
      ium.coverage_stat = 10.0;
      ium.status        = AS_UNIQUE;
      ium.unique_rept   = AS_FORCED_UNIQUE;
      ium.length        = e;
      ium.consensus     = 0L;
      ium.quality       = 0L;
      ium.forced        = 0;
      ium.num_frags     = impLen;
      ium.f_list        = imp;

      pmesg.m = &ium;
      pmesg.t = MESG_IUM;
      WriteProtoMesg_AS(stdout, &pmesg);

      num++;

      b = ali[i].refBgn;
      e = ali[i].refEnd;
      r = ali[i].refUID;

      impLen = 0;
    }

    if (e < ali[i].refEnd)
      e = ali[i].refEnd;

    imp[impLen].type           = AS_READ;
    imp[impLen].ident          = ali[i].frgIID;
    imp[impLen].sourceInt      = 0;
    if (ali[i].fwd) {
      imp[impLen].position.bgn   = ali[i].refBgn - b;
      imp[impLen].position.end   = ali[i].refEnd - b;
    } else {
      imp[impLen].position.bgn   = ali[i].refEnd - b;
      imp[impLen].position.end   = ali[i].refBgn - b;
    }
    imp[impLen].contained      = 0;
    imp[impLen].delta_length   = 0;
    imp[impLen].delta          = 0L;

    impLen++;
  }

  fprintf(stderr, "Wrote "F_U32" unitig messages.\n", num);
}


static
void
outputUnitigLinks(void) {
  //  No unitigs here.
}


static
void
outputContigs(void) {
  GenericMesg     pmesg;
  IntConConMesg   icm;
  int             impMax = 128 * 1024;
  int             impLen = 0;
  IntMultiPos    *imp    = (IntMultiPos *)safe_malloc(impMax * sizeof(IntMultiPos));
  int             iupMax = 128 * 1024;
  int             iupLen = 0;
  IntUnitigPos   *iup    = (IntUnitigPos *)safe_malloc(iupMax * sizeof(IntUnitigPos));
  uint32          num    = 0;

  uint32 icmId = 1;

  uint32 i = 0;
  uint32 b = ali[0].refBgn;
  uint32 e = ali[0].refEnd;
  AS_UID r = ali[0].refUID;

  for (i=0; i<aliLen; i++) {
    //fprintf(stdout, "%s "F_U32"-"F_U32"\n", AS_UID_toString(ali[i].refID), ali[i].refBgn, ali[i].refEnd);

    if ((AS_UID_compare(r, ali[i].refUID) != 0) || (e < ali[i].refBgn)) {
      //fprintf(stderr, "Contig: "F_U32" %s "F_U32"-"F_U32"\n", i, AS_UID_toString(r), b, e);

      iup[iupLen].type           = AS_UNIQUE_UNITIG;
      iup[iupLen].ident          = icmId;  //  One unitig per contig
      iup[iupLen].position.bgn   = 0;
      iup[iupLen].position.end   = e;
      iup[iupLen].delta_length   = 0;
      iup[iupLen].delta          = 0L;

      iupLen++;

      icm.iaccession  = icmId++;
      icm.placed      = AS_PLACED;
      icm.length      = e;
      icm.consensus   = 0L;
      icm.quality     = 0L;
      icm.forced      = 0;
      icm.num_pieces  = impLen;
      icm.num_unitigs = iupLen;
      icm.num_vars    = 0;
      icm.pieces      = imp;
      icm.unitigs     = iup;
      icm.v_list      = 0L;

      pmesg.m = &icm;
      pmesg.t = MESG_ICM;
      WriteProtoMesg_AS(stdout, &pmesg);

      num++;

      b = ali[i].refBgn;
      e = ali[i].refEnd;
      r = ali[i].refUID;

      impLen = 0;
      iupLen = 0;
    }

    if (e < ali[i].refEnd)
      e = ali[i].refEnd;

    imp[impLen].type           = AS_READ;
    imp[impLen].ident          = ali[i].frgIID;
    imp[impLen].sourceInt      = 0;
    if (ali[i].fwd) {
      imp[impLen].position.bgn   = ali[i].refBgn - b;
      imp[impLen].position.end   = ali[i].refEnd - b;
    } else {
      imp[impLen].position.bgn   = ali[i].refEnd - b;
      imp[impLen].position.end   = ali[i].refBgn - b;
    }
    imp[impLen].contained      = 0;
    imp[impLen].delta_length   = 0;
    imp[impLen].delta          = 0L;

    impLen++;
  }

  fprintf(stderr, "Wrote "F_U32" contig messages.\n", num);
}


static
void
outputContigLinks(void) {
  //  No contig links, yet.
}


static
void
outputScaffolds(void) {
  GenericMesg     pmesg;
  IntScaffoldMesg isf;
  int             icpMax = 128 * 1024;
  int             icpLen = 0;
  IntContigPairs *icp    = (IntContigPairs *)safe_malloc(icpMax * sizeof(IntContigPairs));
  uint32          num    = 0;

  uint32 icmId = 1;
  uint32 isfId = 1;

  uint32 i = 0;
  uint32 b = ali[0].refBgn;
  uint32 e = ali[0].refEnd;
  AS_UID r = ali[0].refUID;

  uint32 mean = 0;

  for (i=0; i<aliLen; i++) {
    //fprintf(stdout, "%s "F_U32"-"F_U32"\n", AS_UID_toString(ali[i].refID), ali[i].refBgn, ali[i].refEnd);

    if ((AS_UID_compare(r, ali[i].refUID) != 0) || (e < ali[i].refBgn)) {
      //fprintf(stderr, "Contig: "F_U32" %s "F_U32"-"F_U32"\n", i, AS_UID_toString(r), b, e);

      if (icpLen == 0) {
        //  First contig, assume it's the only contig in the scaffold.
        icp[0].contig1        = icmId;
        icp[0].contig2        = icmId;
        icp[0].mean           = mean;
        icp[0].stddev         = 0;
        icp[0].orient         = AS_NORMAL;
        icpLen++;
      } else if ((icpLen == 1) && (icp[0].contig1 == icp[0].contig2)) {
        //  Second contig.  Make the first icp be a real contig pair.
        icp[0].contig1        = icmId - 1;
        icp[0].contig2        = icmId;
        icp[0].mean           = mean;
        icp[0].stddev         = 0;
        icp[0].orient         = AS_NORMAL;
      } else {
        //  Later contigs, as expected.
        icp[icpLen].contig1        = icmId - 1;
        icp[icpLen].contig2        = icmId;
        icp[icpLen].mean           = mean;
        icp[icpLen].stddev         = 0;
        icp[icpLen].orient         = AS_NORMAL;
        icpLen++;
      }

      if (AS_UID_compare(r, ali[i].refUID) != 0) {
        isf.iaccession       = isfId++;
        isf.num_contig_pairs = (icp[0].contig1 == icp[0].contig2) ? 0 : icpLen;
        isf.contig_pairs     = icp;

        pmesg.m = &isf;
        pmesg.t = MESG_ISF;
        WriteProtoMesg_AS(stdout, &pmesg);

        num++;

        icpLen = 0;
      }

      icmId++;

      mean = ali[i].refBgn - e;

      b = ali[i].refBgn;
      e = ali[i].refEnd;
      r = ali[i].refUID;
    }

    if (e < ali[i].refEnd)
      e = ali[i].refEnd;
  }

  fprintf(stderr, "Wrote "F_U32" scaffold messages.\n", num);
}


static
void
outputScaffoldLinks(void) {
  //  No scaffold links, yet.
}


int
main(int argc, char **argv) {
  char    *mappingFileName   = 0L;
  char    *gkpStoreName      = 0L;
  uint32   buildOnlyUnitigs  = 0;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-m") == 0) {
      mappingFileName = argv[++arg];

    } else if (strcmp(argv[arg], "-U") == 0) {
      buildOnlyUnitigs = 1;

    } else if (strcmp(argv[arg], "-S") == 0) {
      buildOnlyUnitigs = 0;

    } else {
      fprintf(stderr, "Invalid option: '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if ((err) || (gkpStoreName == 0L) || (mappingFileName == 0L)) {
    fprintf(stderr, "usage: %s [-U | -S] -g gkpStore -m mapping\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -g gkpStore\n");
    fprintf(stderr, "  -m mapping\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -U             build unitigs, for input to cgw\n");
    fprintf(stderr, "  -S             build scaffolds, for input to terminator\n");
    exit(1);
  }

  gkpStore = openGateKeeperStore(gkpStoreName, false);
  if (gkpStore == 0L)
    fprintf(stderr, "Failed to open gkpStore '%s'.\n", gkpStoreName);

  readMapping(mappingFileName);
  loadFragments();

  if (buildOnlyUnitigs) {
    outputUnitigs();
  } else {
    outputDistances();
    outputFragments();
    outputMates();
    outputUnitigs();
    outputUnitigLinks();
    outputContigs();
    outputContigLinks();
    outputScaffolds();
    outputScaffoldLinks();
  }

  exit(0);
}
