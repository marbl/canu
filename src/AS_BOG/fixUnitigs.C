
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2009, The Venter Institute. All rights reserved.
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

const char *mainid = "$Id: fixUnitigs.c,v 1.7 2011-09-06 02:15:18 mkotelbajcvi Exp $";

#include "AS_global.H"
#include "AS_MSG_pmesg.H"
#include "AS_OVS_overlapStore.H"
#include "AS_PER_gkpStore.H"


int
updateFragmentWithParent(IntUnitigMesg *iunitig, int thisFrag, OverlapStore *ovs) {
  uint32         ovlMax = 0;
  uint32         ovlLen = 0;
  OVSoverlap    *ovl    = NULL;

  int     testFrag = thisFrag - 1;
  int     testOvl  = 0;

  int     oldParent = iunitig->f_list[thisFrag].parent;
  int     oldAHang  = iunitig->f_list[thisFrag].ahang;
  int     oldBHang  = iunitig->f_list[thisFrag].bhang;

  uint32  consensusCutoff = AS_OVS_encodeQuality(AS_CNS_ERROR_RATE);

  int     contained       = 0;
  int     fragment        = -1;
  int     overlap         = -1;
  int     overlapIdentity = consensusCutoff;
  int     overlapBHang    = AS_READ_MAX_NORMAL_LEN;

  HashTable_AS  *ovlBefore = CreateScalarHashTable_AS();
  HashTable_AS  *ovlAfter  = CreateScalarHashTable_AS();
  HashTable_AS  *iidIndex  = CreateScalarHashTable_AS();

  int     hangSlop = 0;

  int     failed   = -1;

  fprintf(stderr, "\n");
  fprintf(stderr, "WORKING on fragment %d == %d\n", thisFrag, iunitig->f_list[thisFrag].ident);

  //  Save in the hash table the fragments before/after this one.
  //
  for (testFrag=0; testFrag<iunitig->num_frags; testFrag++) {
    InsertInHashTable_AS(iidIndex,
                         (uint64)iunitig->f_list[testFrag].ident, sizeof(uint64),
                         (uint64)testFrag, 0);

    if (testFrag < thisFrag)
      InsertInHashTable_AS(ovlBefore,
                           (uint64)iunitig->f_list[testFrag].ident, sizeof(uint64),
                           ~(uint64)0, 0);
    if (testFrag > thisFrag)
      InsertInHashTable_AS(ovlAfter,
                           (uint64)iunitig->f_list[testFrag].ident, sizeof(uint64),
                           ~(uint64)0, 0);
  }

  //  Get the overlaps for this fragment.
  //
  AS_OVS_setRangeOverlapStore(ovs, iunitig->f_list[thisFrag].ident, iunitig->f_list[thisFrag].ident);

  if (ovlMax < AS_OVS_numOverlapsInRange(ovs)) {
    ovlMax = AS_OVS_numOverlapsInRange(ovs) * 2;
    ovl    = (OVSoverlap *)safe_realloc(ovl, sizeof(OVSoverlap) * ovlMax);
  }
  ovlLen = 0;
  while (AS_OVS_readOverlapFromStore(ovs, ovl+ovlLen, AS_OVS_TYPE_OVL)) {
    int  aid=0,  bid=0;
    int  afwd=0, bfwd=0;
    int  correct=0;

    //  Reorient the overlap so the b_iid is thisFrag.
    //
    {
      AS_IID x = ovl[ovlLen].a_iid;
      ovl[ovlLen].a_iid = ovl[ovlLen].b_iid;
      ovl[ovlLen].b_iid = x;

      if (ovl[ovlLen].dat.ovl.flipped) {
        int x = ovl[ovlLen].dat.ovl.a_hang;
        ovl[ovlLen].dat.ovl.a_hang = ovl[ovlLen].dat.ovl.b_hang;
        ovl[ovlLen].dat.ovl.b_hang = x;
      } else {
        ovl[ovlLen].dat.ovl.a_hang = -ovl[ovlLen].dat.ovl.a_hang;
        ovl[ovlLen].dat.ovl.b_hang = -ovl[ovlLen].dat.ovl.b_hang;
      }
    }

    //  Make sure we get the correct overlap.  We seem to be allowed
    //  to have both an I and an N overlap for a given pair of
    //  fragments.  At least, I hope that's all we're allowed.
    //
    aid = LookupValueInHashTable_AS(iidIndex, (uint64)ovl[ovlLen].a_iid, sizeof(uint64));
    bid = LookupValueInHashTable_AS(iidIndex, (uint64)ovl[ovlLen].b_iid, sizeof(uint64));

    afwd = (iunitig->f_list[aid].position.bgn < iunitig->f_list[aid].position.end);
    bfwd = (iunitig->f_list[bid].position.bgn < iunitig->f_list[bid].position.end);

    if ((afwd == bfwd) && (ovl[ovlLen].dat.ovl.flipped == 0))
      correct = 1;
    if ((afwd != bfwd) && (ovl[ovlLen].dat.ovl.flipped == 1))
      correct = 1;


    if (ExistsInHashTable_AS(ovlBefore, (uint64)ovl[ovlLen].a_iid, sizeof(uint64))) {
      if (correct)
        ReplaceInHashTable_AS(ovlBefore,
                              (uint64)ovl[ovlLen].a_iid, sizeof(uint64),
                              (uint64)ovlLen, 0);
      fprintf(stderr, "%s before overlap for %d (%c) to %d (%c) ("F_S64","F_S64",%c) at ovl position %d\n",
              correct ? "save" : "skip",
              ovl[ovlLen].a_iid, afwd ? 'F' : 'R',
              ovl[ovlLen].b_iid, bfwd ? 'F' : 'R',
              ovl[ovlLen].dat.ovl.a_hang,
              ovl[ovlLen].dat.ovl.b_hang,
              ovl[ovlLen].dat.ovl.flipped ? 'I' : 'N',
              ovlLen);
    }

    if (ExistsInHashTable_AS(ovlAfter, (uint64)ovl[ovlLen].a_iid, sizeof(uint64))) {
      if (correct)
        ReplaceInHashTable_AS(ovlAfter,
                              (uint64)ovl[ovlLen].a_iid, sizeof(uint64),
                              (uint64)ovlLen, 0);
      fprintf(stderr, "%s after  overlap for %d (%c) to %d (%c) ("F_S64","F_S64",%c) at ovl position %d\n",
              correct ? "save" : "skip",
              ovl[ovlLen].a_iid, afwd ? 'F' : 'R',
              ovl[ovlLen].b_iid, bfwd ? 'F' : 'R',
              ovl[ovlLen].dat.ovl.a_hang,
              ovl[ovlLen].dat.ovl.b_hang,
              ovl[ovlLen].dat.ovl.flipped ? 'I' : 'N',
              ovlLen);
    }

    ovlLen++;
  }

 tryAgain:

  //  See if we're contained in any of these overlaps.
  if (overlap == -1) {
    for (testFrag=thisFrag-1; testFrag>=0; testFrag--) {
      if (ExistsInHashTable_AS(ovlBefore, (uint64)iunitig->f_list[testFrag].ident, sizeof(uint64))) {
        testOvl = LookupValueInHashTable_AS(ovlBefore, (uint64)iunitig->f_list[testFrag].ident, sizeof(uint64));

        //  Fragment has no overlap
        if (testOvl == -1)
          continue;

        fprintf(stderr, "found testFrag = %d testOvl = %d erates "F_U64" %u hang "F_S64" "F_S64" (CONTAIN) slop=%d\n",
                testFrag, testOvl,
                ovl[testOvl].dat.ovl.orig_erate, consensusCutoff,
                ovl[testOvl].dat.ovl.a_hang,
                ovl[testOvl].dat.ovl.b_hang,
                hangSlop);

        //  Three if's for documentation:
        //  1)  If we're an overlap we care about
        //  2)  If we're a contained overlap
        //  3)  If we're better than what we've seen so far
        //  Then save the overlap
        //
        if (ovl[testOvl].dat.ovl.orig_erate < consensusCutoff) {
          if ((ovl[testOvl].dat.ovl.a_hang >= -hangSlop) &&
              (ovl[testOvl].dat.ovl.b_hang <= hangSlop)) {
            if (ovl[testOvl].dat.ovl.orig_erate < overlapIdentity) {
              contained       = 1;
              fragment        = testFrag;
              overlap         = testOvl;
              overlapBHang    = 0;
              overlapIdentity = ovl[testOvl].dat.ovl.orig_erate;
            }
          }
        }
      }
    }
  }


  //  If not contained, scan the overlaps again, looking for the
  //  thickest/bestest.  This will be the overlap with the smallest a
  //  or b hang -- depending on the orientation of the parent
  //  fragment.
  //
  //  Instead of working through overlaps, we work through fragments.
  //
  if (overlap == -1) {
    for (testFrag=thisFrag-1; testFrag>=0; testFrag--) {
      if (ExistsInHashTable_AS(ovlBefore, (uint64)iunitig->f_list[testFrag].ident, sizeof(uint64))) {
        int ahang = 0;
        int bhang = 0;

        testOvl = LookupValueInHashTable_AS(ovlBefore, (uint64)iunitig->f_list[testFrag].ident, sizeof(uint64));

        //  Fragment has no overlap
        if (testOvl == -1)
          continue;

        //  Overlap is too noisy
        if (ovl[testOvl].dat.ovl.orig_erate >= consensusCutoff)
          continue;

        if (iunitig->f_list[testFrag].position.bgn < iunitig->f_list[testFrag].position.end) {
          ahang = ovl[testOvl].dat.ovl.a_hang;
          bhang = ovl[testOvl].dat.ovl.b_hang;
        } else {
          ahang = -ovl[testOvl].dat.ovl.b_hang;
          bhang = -ovl[testOvl].dat.ovl.a_hang;
        }

        //  Overlap isn't dovetail -- negative ahang
        if (ahang < 0)
          continue;

        //  Overlap isn't dovetail -- containment
        if (bhang < 0)
          continue;

        fprintf(stderr, "found testFrag = %d testOvl = %d erates "F_U64" %u hang "F_S64" "F_S64" (DOVETAIL) slop=%d\n",
                testFrag, testOvl,
                ovl[testOvl].dat.ovl.orig_erate, consensusCutoff,
                ovl[testOvl].dat.ovl.a_hang,
                ovl[testOvl].dat.ovl.b_hang,
                hangSlop);

        if (bhang < overlapBHang) {
          contained       = 0;
          fragment        = testFrag;
          overlap         = testOvl;
          overlapIdentity = ovl[testOvl].dat.ovl.orig_erate;
          overlapBHang    = bhang;
        }
      }
    }
  }


  //  Now, if we have found the parent fragment, update.
  //
  if (overlap >= 0) {
    testOvl  = overlap;
    testFrag = fragment;

    iunitig->f_list[thisFrag].parent = ovl[testOvl].a_iid;

    if (contained)
      iunitig->f_list[thisFrag].contained = iunitig->f_list[thisFrag].parent;
    else
      iunitig->f_list[thisFrag].contained = 0;

    //  Reorient again based on the orientation of the testFrag.
    //
    if (iunitig->f_list[testFrag].position.bgn < iunitig->f_list[testFrag].position.end) {
      //  testFrag is forward
      iunitig->f_list[thisFrag].ahang  = ovl[testOvl].dat.ovl.a_hang;
      iunitig->f_list[thisFrag].bhang  = ovl[testOvl].dat.ovl.b_hang;
    } else {
      //  testFrag is reverse
      iunitig->f_list[thisFrag].ahang  = -ovl[testOvl].dat.ovl.b_hang;
      iunitig->f_list[thisFrag].bhang  = -ovl[testOvl].dat.ovl.a_hang;
    }

    //  Report we did something.
    //
    fprintf(stderr, "Updated fragment "F_IID" from "F_IID",%d,%d to "F_IID",%d,%d\n",
            iunitig->f_list[thisFrag].ident,
            oldParent,
            oldAHang,
            oldBHang,
            iunitig->f_list[thisFrag].parent,
            iunitig->f_list[thisFrag].ahang,
            iunitig->f_list[thisFrag].bhang);

    goto successfullyUpdated;
  }


  //  Otherwise, try to find an overlap again, this time allowing a
  //  bit of slop in the hangs.
  //
  if (hangSlop == 0) {
    hangSlop = 10;
    goto tryAgain;
  }


  //  Now, we're convinced there is no decent overlap between this
  //  fragment and any fragment before it.
  //
  //  Scan forward for the first thing we overlap.

  for (testFrag=thisFrag+1; testFrag < iunitig->num_frags; testFrag++) {
    int ahang = 0;
    int bhang = 0;

    testOvl = LookupValueInHashTable_AS(ovlAfter, (uint64)iunitig->f_list[testFrag].ident, sizeof(uint64));

    //  Fragment has no overlap
    if (testOvl == -1)
      continue;

    //  Overlap is too noisy
    if (ovl[testOvl].dat.ovl.orig_erate >= consensusCutoff)
      continue;

    if (iunitig->f_list[testFrag].position.bgn < iunitig->f_list[testFrag].position.end) {
      ahang = ovl[testOvl].dat.ovl.a_hang;
      bhang = ovl[testOvl].dat.ovl.b_hang;
    } else {
      ahang = -ovl[testOvl].dat.ovl.b_hang;
      bhang = -ovl[testOvl].dat.ovl.a_hang;
    }

    //  Don't allow negative ahangs.  At all.  This catches the case
    //  where the parent might be contained in us, and generally makes
    //  consensus happier.
    //
    //  Don't allow empty hangs - this can lead to infinite loops
    //  where we keep swapping the same two fragments.  OK, not
    //  infinite, since we eventually run out of stack space and
    //  crash.
    //
    if (ahang <= 0)
      continue;

    fprintf(stderr, "shifttest ovl=%d testFrag="F_IID" pos %d-%d  thisFrag="F_IID" pos %d-%d  hangs %d,%d\n",
            testOvl,
            iunitig->f_list[testFrag].ident,
            iunitig->f_list[testFrag].position.bgn,
            iunitig->f_list[testFrag].position.end,
            iunitig->f_list[thisFrag].ident,
            iunitig->f_list[thisFrag].position.bgn,
            iunitig->f_list[thisFrag].position.end,
            ahang, bhang);

    IntMultiPos  fragCopy = iunitig->f_list[thisFrag];

    memmove(iunitig->f_list + thisFrag,
            iunitig->f_list + thisFrag + 1,
            sizeof(IntMultiPos) * (testFrag - thisFrag));

    iunitig->f_list[testFrag] = fragCopy;

    fprintf(stderr, "Shifted fragment "F_IID" from position %d to position %d\n",
            iunitig->f_list[testFrag].ident,
            thisFrag, testFrag);

    //  Since we moved things around, we must process the new fragment
    //  at 'thisFrag's location.
    //
    failed = updateFragmentWithParent(iunitig, thisFrag, ovs);

    if (failed == -1)
      goto successfullyUpdated;

    break;
  }


  //  And we failed.  Good luck with this one.
  //
  fprintf(stderr, "Failed to update fragment "F_IID" from "F_IID",%d,%d.\n",
          iunitig->f_list[thisFrag].ident,
          oldParent,
          oldAHang,
          oldBHang);

  failed = thisFrag;

 successfullyUpdated:
  DeleteHashTable_AS(ovlBefore);
  DeleteHashTable_AS(ovlAfter);
  safe_free(ovl);

  return(failed);
}




void
fixUnitig(IntUnitigMesg *iunitig, OverlapStore *ovs) {
  int            thisFrag;
  int            thatFrag;

  for (thisFrag=1; thisFrag<iunitig->num_frags; thisFrag++) {
    int failed = updateFragmentWithParent(iunitig, thisFrag, ovs);

    //  If that failed, the iunitig is guaranteed good up until the
    //  'failed' fragment.  It'll get written out back in main; all we
    //  need to do is fix up the rest of the fragments, possibly into
    //  multiple unitigs.

    if (failed != -1) {
      IntUnitigMesg  junitig = *iunitig;

      assert(failed == thisFrag);

      //  Make the iacc big, just to label this as needing a new iacc.
      junitig.iaccession += 1000000000;

      junitig.num_frags  = iunitig->num_frags - failed;
      junitig.f_list     = iunitig->f_list    + failed;

      junitig.f_list[0].parent    = 0;
      junitig.f_list[0].ahang     = 0;
      junitig.f_list[0].bhang     = 0;
      junitig.f_list[0].contained = 0;

      iunitig->num_frags = failed;

      fixUnitig(&junitig, ovs);

      GenericMesg   pmesg;

      pmesg.t = MESG_IUM;
      pmesg.m = &junitig;

      WriteProtoMesg_AS(stdout, &pmesg);
    }
  }
}



int
main(int argc, char **argv) {
  OverlapStore  *ovs        = NULL;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-O") == 0) {
      ovs = AS_OVS_openOverlapStore(argv[++arg]);
    } else {
      err++;
    }

    arg++;
  }
  if ((ovs == NULL) || (err)) {
    fprintf(stderr, "usage: %s -O ovlStore < unitigs.cgb > fixedUnitigs.cgb\n", argv[0]);
    exit(1);
  }

  GenericMesg   *pmesg = NULL;
  while ((ReadProtoMesg_AS(stdin, &pmesg) != EOF)) {
    if (pmesg->t == MESG_IUM)
      fixUnitig((IntUnitigMesg *)(pmesg->m), ovs);

    WriteProtoMesg_AS(stdout, pmesg);
  }

  exit(0);
}
