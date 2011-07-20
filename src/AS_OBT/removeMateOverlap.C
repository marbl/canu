
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005-2007, J. Craig Venter Institute.
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

const char *mainid = "$Id: removeMateOverlap.C,v 1.4 2011-07-20 20:01:37 mkotelbajcvi Exp $";

//  Remove mate relationships for any fragments that overlap

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "readOverlap.H"

#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "AS_OVS_overlapStore.h"



class fragT {
public:
  uint32   doDelete:1;
  uint32   isDeleted:1;
  uint32   mateIID;
};



fragT *
loadFragments(gkStore *gkp) {
  gkStream       *fs = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);
  gkFragment      fr;
  fragT          *frag = new fragT [gkp->gkStore_getNumFragments() + 1];

  while (fs->next(&fr)) {
    AS_IID  iid  = fr.gkFragment_getReadIID();

    frag[iid].doDelete         = 0;
    frag[iid].isDeleted        = fr.gkFragment_getIsDeleted() ? 1 : 0;
    frag[iid].mateIID          = fr.gkFragment_getMateIID();

    //  HACK!  Ignore non-illumina reads
    if (fr.gkFragment_getLibraryIID() != 1)
      frag[iid].mateIID = 0;
  }

  delete fs;

  return(frag);
}






int
main(int argc, char **argv) {
  gkStore           *gkpstore     = 0L;
  OverlapStore      *ovlstore     = 0L;
  OVSoverlap         ovl;

  uint64             totalOverlaps = 0;
  uint64             matesRemoved  = 0;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-gkp", 2) == 0) {
      gkpstore = new gkStore(argv[++arg], FALSE, TRUE);

    } else if (strncmp(argv[arg], "-ovl", 2) == 0) {
      ovlstore = AS_OVS_openOverlapStore(argv[++arg]);

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((gkpstore == 0L) || (ovlstore == 0L) || (err)) {
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    exit(1);
  }

  fragT  *frag = loadFragments(gkpstore);



  while (AS_OVS_readOverlapFromStore(ovlstore, &ovl, AS_OVS_TYPE_OVL)) {
    totalOverlaps++;

    if ((frag[ovl.a_iid].mateIID == 0) ||
        (frag[ovl.b_iid].mateIID == 0) ||
        (frag[ovl.a_iid].mateIID != ovl.b_iid) ||
        (frag[ovl.b_iid].mateIID != ovl.a_iid))
      //  Frags must be mated to each other.
      continue;

    if (ovl.dat.ovl.flipped == 0)
      //  Innie mate overlaps must be flipped
      continue;

    if ((ovl.dat.ovl.a_hang < 0) || (ovl.dat.ovl.b_hang < 0))
      //  ...And on the innie pointing end
      continue;

    //  ...Before we delete.

    matesRemoved++;

    frag[ovl.a_iid].doDelete = 1;
    frag[ovl.b_iid].doDelete = 1;

    //  Print the overlapping sequence
    {
      uint32  iid = ovl.a_iid;
      uint32  mid = ovl.b_iid;

      gkFragment  fr;
      gkFragment  mr;

      gkpstore->gkStore_getFragment(iid, &fr, GKFRAGMENT_SEQ);
      gkpstore->gkStore_getFragment(mid, &mr, GKFRAGMENT_SEQ);

      fprintf(stdout, ">%d,%d\n%s\n",
              ovl.dat.ovl.a_hang,
              ovl.dat.ovl.b_hang,
              fr.gkFragment_getSequence() + ovl.dat.ovl.a_hang);
    }

    //fprintf(stderr, "DELETE %d <-> %d\n", ovl.a_iid, ovl.b_iid);
  }


  for (uint32 iid=0; iid<=gkpstore->gkStore_getNumFragments(); iid++) {
    if ((frag[iid].doDelete) && (frag[iid].isDeleted == 0)) {
      uint32      mid = frag[iid].mateIID;
      gkFragment  fr;
      gkFragment  mr;

      gkpstore->gkStore_getFragment(iid, &fr, GKFRAGMENT_INF);
      gkpstore->gkStore_getFragment(mid, &mr, GKFRAGMENT_INF);

      fr.gkFragment_setMateIID(0);
      mr.gkFragment_setMateIID(0);

      //gkpstore->gkStore_setFragment(&fr);
      //gkpstore->gkStore_setFragment(&mr);

      frag[iid].isDeleted = 1;
      frag[mid].isDeleted = 1;

      //fprintf(stdout, ">%s\n%s\n", AS_UID_toString(fr.gkFragment_getReadUID()), fr.gkFragment_getSequence());
      //fprintf(stdout, ">%s\n%s\n", AS_UID_toString(mr.gkFragment_getReadUID()), mr.gkFragment_getSequence());
    }
  }

  delete [] frag;
  delete    gkpstore;

  fprintf(stderr, "Total overlaps   "F_U64"\n", totalOverlaps);
  fprintf(stderr, "Mates removed    "F_U64"\n", matesRemoved);

  exit(0);
}
