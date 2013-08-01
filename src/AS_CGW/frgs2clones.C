
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

const char *mainid = "$Id$";

#include "AS_global.H"
#include "AS_PER_gkpStore.H"
#include "AS_UTL_reverseComplement.H"
#include "SYS_UIDclient.H"
#include "MultiAlignment_CNS.H"

int
main(int argc, char **argv) {
  char   *gkpStorePath   = NULL;
  int     runConsensus   = 0;
  int     Ngaps          = 0;
  uint64  UIDstart       = 1230000;

  argc = AS_configure(argc, argv);

  int  arg=1;
  int  err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpStorePath = argv[++arg];
      
    } else if (strcmp(argv[arg], "-C") == 0) {
      runConsensus = 1;

    } else if (strcmp(argv[arg], "-N") == 0) {
      Ngaps = 1;

    } else if (strcmp(argv[arg], "-U") == 0) {
      UIDstart = 0;

    } else {
      err++;
    }

    arg++;
  }
  if ((err) || (gkpStorePath == NULL)) {
    fprintf(stderr, "usage: %s -g gkpStore [opts]\n", argv[0]);
    fprintf(stderr, "  -g gkpStore    read fragments from this gatekeeper store\n");
    fprintf(stderr, "  -C             build an alignment for overlapping fragments\n");
    fprintf(stderr, "  -N             build a scaffold for non-overlapping mated fragments\n");
    fprintf(stderr, "  -U             use real UIDs from the UID server\n");
    exit(1);
  }

  UIDserver *uidServer = UIDserverInitialize(256, UIDstart);

  gkStore    *gs = new gkStore(gkpStorePath, FALSE, FALSE);
  gkStream   *fs = new gkStream(gs, 0, 0, GKFRAGMENT_SEQ);

  gkFragment fr;
  gkFragment fm;

  char       *frseq = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_NORMAL_LEN);
  char       *fmseq = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_NORMAL_LEN);

  while (fs->next(&fr)) {
    if (fr.gkFragment_getIsDeleted())
      continue;

    //  Get the sequence of the fragment
    //
    uint32   frbgn = fr.gkFragment_getClearRegionBegin();
    uint32   frend = fr.gkFragment_getClearRegionEnd();
    uint32   frlen = frend - frbgn;

    strcpy(frseq, fr.gkFragment_getSequence() + frbgn);    frseq[frlen] = 0;

    //  No mate pair, just emit the sequence.
    //
    if (fr.gkFragment_getMateIID() == 0){
      fprintf(stdout, ">%s\n%s\n", AS_UID_toString(fr.gkFragment_getReadUID()), frseq);
      continue;
    }

    //  Only do the pair once.
    //
    if (fr.gkFragment_getMateIID() < fr.gkFragment_getReadIID())
      continue;

    //  Get the sequence of the mate
    //
    gs->gkStore_getFragment(fr.gkFragment_getMateIID(), &fm, GKFRAGMENT_SEQ);

    uint32   fmbgn = fm.gkFragment_getClearRegionBegin();
    uint32   fmend = fm.gkFragment_getClearRegionEnd();
    uint32   fmlen = fmend - fmbgn;

    strcpy(fmseq, fm.gkFragment_getSequence() + fmbgn);    fmseq[fmlen] = 0;

    uint64  mergeUID = getUID(uidServer);

    //fprintf(stderr, ">A %d len=%d\n%s\n", fr.gkFragment_getReadIID(), frlen, frseq);
    //fprintf(stderr, ">B %d len=%d\n%s\n", fr.gkFragment_getMateIID(), fmlen, fmseq);

    ALNoverlap *ovl = Local_Overlap_AS_forCNS(frseq, fmseq,
                                              -fmlen, frlen,
                                              fmlen,  frlen,
                                              1,
                                              AS_CNS_ERROR_RATE,
                                              1e-6,
                                              40,
                                              AS_FIND_LOCAL_ALIGN_NO_TRACE);

    //  If they don't overlap reasonably, join together with a gap, or as two sequences.
    //
    if ((ovl == NULL) ||
        (((ovl->begpos < 0) || (ovl->endpos < 0)) &&
         (frlen + fmlen - abs(ovl->begpos) - abs(ovl->endpos)) / 2 < 100)) {

      if(Ngaps)
        fprintf(stdout, ">"F_S64" from mated fragments %s and %s\n%sNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN%s\n",
                mergeUID,
                AS_UID_toString(fr.gkFragment_getReadUID()),
                AS_UID_toString(fm.gkFragment_getReadUID()),
                frseq, fmseq);
      else
        fprintf(stdout, ">"F_S64"a (fragment %s)\n%s\n>"F_S64"b (fragment %s)\n%s\n",
                mergeUID, AS_UID_toString(fr.gkFragment_getReadUID()), frseq,
                mergeUID, AS_UID_toString(fm.gkFragment_getReadUID()), fmseq);

      continue;
    }

    //  Otherwise, we have an overlap.

    //fprintf(stderr, ">A %d len=%d\n%s\n", fr.gkFragment_getReadIID(), frlen, frseq);
    //fprintf(stderr, ">B %d len=%d\n%s\n", fr.gkFragment_getMateIID(), fmlen, fmseq);
    //fprintf(stderr, "OVERLAP: begpos=%d endpos=%d  length=%d  diffs=%d  comp=%d\n",
    //        ovl->begpos, ovl->endpos, ovl->length, ovl->diffs, ovl->comp);

    //  Emit a concatenation of the sequences.
    //
    if (runConsensus == 0) {
      reverseComplementSequence(fmseq, fmlen);

      fprintf(stdout, ">"F_S64" merged sequence of mated fragments %s and %s\n%s%s\n",
              mergeUID,
              AS_UID_toString(fr.gkFragment_getReadUID()),
              AS_UID_toString(fm.gkFragment_getReadUID()),
              frseq, fmseq + fmlen - ovl->endpos);

      continue;
    }

    //  Emit an alignment of the sequences.
    //
    MultiAlignT  *ma = CreateEmptyMultiAlignT();
    IntMultiPos   imps[2];

    //  Based on the overlap, decide where to put the fragment and its mate so
    //  that the fragments are ordered correctly.
    //
    int f = (ovl->begpos >= 0) ? 0 : 1;
    int s = (ovl->begpos >= 0) ? 1 : 0;

    //  Add the first fragment, forward.
    //
    imps[f].type         = AS_READ;
    imps[f].ident        = fr.gkFragment_getReadIID();
    imps[f].contained    = 0;
    imps[f].position.bgn = (ovl->begpos >= 0) ? (0)     : (0     - ovl->begpos);
    imps[f].position.end = (ovl->begpos >= 0) ? (frlen) : (frlen - ovl->begpos);
    imps[f].delta_length = 0;
    imps[f].delta        = NULL;

    //  Add the mate, reversed.
    //
    imps[s].type         = AS_READ;
    imps[s].ident        = fm.gkFragment_getReadIID();
    imps[s].contained    = 0;
    imps[s].position.bgn = (ovl->begpos >= 0) ? (ovl->begpos + fmlen) : (fmlen);
    imps[s].position.end = (ovl->begpos >= 0) ? (ovl->begpos)         : (0);
    imps[s].delta_length = 0;
    imps[s].delta        = NULL;

    //  Build the MultiAlign

    ma->maID                      = 1;
    ma->data.unitig_coverage_stat = 0.0;
    ma->data.unitig_microhet_prob = 1.0;

    ma->data.unitig_status        = AS_UNASSIGNED;
    ma->data.unitig_unique_rept   = AS_FORCED_NONE;

    ma->data.contig_status        = AS_UNPLACED;

    //  Add the fragments

    ResetVA_IntMultiPos(ma->f_list);
    SetRangeVA_IntMultiPos(ma->f_list, 0, 2, imps);

    //  Call consensus

    //VERBOSE_MULTIALIGN_OUTPUT = 1;

    if (MultiAlignUnitig(ma, gs, NULL, NULL) == 0)
      fprintf(stderr, "MultiAlignUnitig() failed for overlap of fragments %s and %s\n",
              AS_UID_toString(fr.gkFragment_getReadUID()),
              AS_UID_toString(fm.gkFragment_getReadUID())), exit(1);

    VA_TYPE(char)  *cns = CreateVA_char(2 * AS_READ_MAX_NORMAL_LEN);
    VA_TYPE(char)  *qlt = CreateVA_char(2 * AS_READ_MAX_NORMAL_LEN);

    GetMultiAlignUngappedConsensus(ma, cns, qlt);

    fprintf(stdout, ">"F_S64" merged sequence of mated fragments %s and %s\n%s\n",
            mergeUID,
            AS_UID_toString(fr.gkFragment_getReadUID()),
            AS_UID_toString(fm.gkFragment_getReadUID()),
            Getchar(cns, 0));

    Delete_VA(cns);
    Delete_VA(qlt);

    DeleteMultiAlignT(ma);
  }

  exit(0);
}
