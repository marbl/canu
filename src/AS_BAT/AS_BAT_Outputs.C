
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
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

static const char *rcsid = "$Id: AS_BAT_Outputs.C,v 1.3 2011-12-29 09:26:03 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Instrumentation.H"

#include "AS_BAT_Outputs.H"

#include "AS_CGB_histo.h"
#include "MultiAlignStore.h"




//  Massage the Unitig into a MultiAlignT (also used in SplitChunks_CGW.c)
void
unitigToMA(MultiAlignT *ma,
           uint32       iumiid,
           Unitig      *utg) {

  ma->maID                      = iumiid;
  ma->data.unitig_coverage_stat = utg->getCovStat();
  ma->data.unitig_microhet_prob = 1.0;  //  Default to 100% probability of unique

  ma->data.unitig_status        = AS_UNASSIGNED;
  ma->data.unitig_unique_rept   = AS_FORCED_NONE;

  ma->data.contig_status        = AS_UNPLACED;

  //  Add the fragments

  ResetVA_IntMultiPos(ma->f_list);

  for (uint32 fi=0; fi<utg->ufpath.size(); fi++) {
    ufNode        *frg = &utg->ufpath[fi];
    IntMultiPos    imp;

    imp.type         = AS_READ;
    imp.ident        = frg->ident;
    imp.contained    = frg->contained;
    imp.parent       = frg->parent;
    imp.ahang        = frg->ahang;
    imp.bhang        = frg->bhang;
    imp.position.bgn = frg->position.bgn;
    imp.position.end = frg->position.end;
    imp.delta_length = 0;
    imp.delta        = NULL;

    AppendVA_IntMultiPos(ma->f_list, &imp);
  }
}



void
writeIUMtoFile(UnitigVector  &unitigs,
               char          *fileprefix,
               char          *tigStorePath,
               uint32         frg_count_target,
               bool           isFinal) {
  uint32      utg_count              = 0;
  uint32      frg_count              = 0;
  uint32      prt_count              = 1;
  char        filename[FILENAME_MAX] = {0};
  uint32     *partmap                = new uint32 [unitigs.size()];

  //  This code closely follows that in AS_CGB_unitigger.c::output_the_chunks()

  if (isFinal)
    checkUnitigMembership(unitigs);

  // Open up the initial output file

  sprintf(filename, "%s.iidmap", fileprefix);
  FILE *iidm = fopen(filename, "w");
  assert(NULL != iidm);

  sprintf(filename, "%s.partitioning", fileprefix);
  FILE *part = fopen(filename, "w");
  assert(NULL != part);

  sprintf(filename, "%s.partitioningInfo", fileprefix);
  FILE *pari = fopen(filename, "w");
  assert(NULL != pari);

  //  Step through all the unitigs once to build the partition mapping and IID mapping.

  memset(partmap, 0xff, sizeof(uint32) * unitigs.size());

  for (uint32 iumiid=0, ti=0; ti<unitigs.size(); ti++) {
    Unitig  *utg = unitigs[ti];
    uint32   nf  = (utg) ? utg->getNumFrags() : 0;

    if ((utg == NULL) || (nf == 0))
      continue;

    assert(utg->getLength() > 0);
    assert(nf == utg->ufpath.size());

    if ((frg_count + nf >= frg_count_target) &&
        (frg_count      >  0)) {
      fprintf(pari, "Partition %d has %d unitigs and %d fragments.\n",
              prt_count, utg_count, frg_count);

      prt_count++;
      utg_count = 0;
      frg_count = 0;
    }

    uint32 tigid = (isFinal) ? iumiid : ti;

    assert(tigid < unitigs.size());
    partmap[tigid] = prt_count;

    fprintf(iidm, "Unitig "F_U32" == IUM "F_U32" (in partition "F_U32" with "F_U32" frags)\n",
            utg->id(),
            (tigid),
            partmap[(tigid)],
            nf);

    for (uint32 fragIdx=0; fragIdx<nf; fragIdx++) {
      ufNode  *f = &utg->ufpath[fragIdx];

      fprintf(part, "%d\t%d\n", prt_count, f->ident);
    }

    utg_count += 1;
    frg_count += nf;

    iumiid++;
  }

  fprintf(pari, "Partition %d has %d unitigs and %d fragments.\n",
          prt_count, utg_count, frg_count);

  fclose(pari);
  fclose(part);
  fclose(iidm);

  //  Step through all the unitigs once to build the partition mapping and IID mapping.

  MultiAlignStore  *MAS = new MultiAlignStore(tigStorePath);
  MultiAlignT      *ma  = CreateEmptyMultiAlignT();

  MAS->writeToPartitioned(partmap, unitigs.size(), NULL, 0);

  for (uint32 iumiid=0, ti=0; ti<unitigs.size(); ti++) {
    Unitig  *utg = unitigs[ti];
    uint32   nf  = (utg) ? utg->getNumFrags() : 0;

    if ((utg == NULL) || (nf == 0))
      continue;

    unitigToMA(ma, (isFinal) ? iumiid : ti, utg);

    //  NOTE!  This is not currently a valid multialign as it has NO IntUnitigPos.  That is
    //  added during consensus.  CGW will correctly assert that it reads in unitigs with
    //  exactly one IUP.

    //  Stash the unitig in the store

    MAS->insertMultiAlign(ma, TRUE, FALSE);

    iumiid++;
  }

  DeleteMultiAlignT(ma);

  delete    MAS;
  delete [] partmap;
}


//  For every unitig, report the best overlaps contained in the
//  unitig, and all overlaps contained in the unitig.
void
writeOVLtoFile(UnitigVector &unitigs,
               char         *fileprefix) {
  char         filename[FILENAME_MAX] = {0};
  GenericMesg  pmesg;
  OverlapMesg  omesg;

  sprintf(filename, "%s.unused.ovl", fileprefix);
  FILE *file = fopen(filename, "w");
  assert(file != NULL);

  for (uint32  ti=0; ti<unitigs.size(); ti++) {
    Unitig  *utg = unitigs[ti];

    if (utg == NULL)
      continue;

    for (uint32 fi=0; fi<utg->ufpath.size(); fi++) {
      ufNode  *frg = &utg->ufpath[fi];

      //  Where is our best overlap?  Contained or dovetail?

      BestEdgeOverlap *bestedge5 = OG->getBestEdgeOverlap(frg->ident, false);
      BestEdgeOverlap *bestedge3 = OG->getBestEdgeOverlap(frg->ident, true);

      int              bestident5 = 0;
      int              bestident3 = 0;

      if (bestedge5) {
        bestident5 = bestedge5->fragId();

        if ((bestident5 > 0) && (utg->fragIn(bestident5) != utg->id())) {
          omesg.aifrag          = frg->ident;
          omesg.bifrag          = bestident5;
          omesg.ahg             = bestedge5->ahang();
          omesg.bhg             = bestedge5->bhang();
          omesg.orientation.setIsUnknown();
          omesg.overlap_type    = AS_DOVETAIL;
          omesg.quality         = 0.0;
          omesg.min_offset      = 0;
          omesg.max_offset      = 0;
          omesg.polymorph_ct    = 0;
          omesg.alignment_trace = NULL;
#ifdef AS_MSG_USE_OVL_DELTA
          omesg.alignment_delta = NULL;
#endif

          //  This overlap is off of the 5' end of this fragment.
          if (bestedge5->frag3p() == false)
            omesg.orientation.setIsOuttie();
          if (bestedge5->frag3p() == true)
            omesg.orientation.setIsAnti();

          pmesg.t = MESG_OVL;
          pmesg.m = &omesg;

          WriteProtoMesg_AS(file, &pmesg);
        }
      }

      if (bestedge3) {
        bestident3 = bestedge3->fragId();

        if ((bestident3 > 0) && (utg->fragIn(bestident3) != utg->id())) {
          omesg.aifrag          = frg->ident;
          omesg.bifrag          = bestident3;
          omesg.ahg             = bestedge3->ahang();
          omesg.bhg             = bestedge3->bhang();
          omesg.orientation.setIsUnknown();
          omesg.overlap_type    = AS_DOVETAIL;
          omesg.quality         = 0.0;
          omesg.min_offset      = 0;
          omesg.max_offset      = 0;
          omesg.polymorph_ct    = 0;
          omesg.alignment_trace = NULL;
#ifdef AS_MSG_USE_OVL_DELTA
          omesg.alignment_delta = NULL;
#endif

          //  This overlap is off of the 3' end of this fragment.
          if (bestedge3->frag3p() == false)
            omesg.orientation.setIsNormal();
          if (bestedge3->frag3p() == true)
            omesg.orientation.setIsInnie();

          pmesg.t = MESG_OVL;
          pmesg.m = &omesg;

          WriteProtoMesg_AS(file, &pmesg);
        }
      }
    }
  }

  fclose(file);
}



void
writeCGAtoFile(UnitigVector &unitigs,
               char         *outputprefix,
               double        globalARate) {
  char  filename[FILENAME_MAX] = {0};
  sprintf(filename, "%s.cga.0", outputprefix);

  errno = 0;
  FILE *stats = fopen(filename,"w");
  if (errno) {
    fprintf(logFile, "Failed to open '%s' for writing: %s\n", filename, strerror(errno));
    fprintf(logFile, "Writing CGA to the log instead.\n");
    stats = logFile;
  }

  fprintf(stats, "Global Arrival Rate: %f\n", globalARate);

  const int nsample = 500;
  const int nbucket = 500;

  MyHistoDataType  z;

  memset(&z, 0, sizeof(MyHistoDataType));

  Histogram_t *len = create_histogram(nsample, nbucket, 0, TRUE);
  Histogram_t *cvg = create_histogram(nsample, nbucket, 0, TRUE);
  Histogram_t *arv = create_histogram(nsample, nbucket, 0, TRUE);

  extend_histogram(len, sizeof(MyHistoDataType), myindexdata, mysetdata, myaggregate, myprintdata);
  extend_histogram(cvg, sizeof(MyHistoDataType), myindexdata, mysetdata, myaggregate, myprintdata);
  extend_histogram(arv, sizeof(MyHistoDataType), myindexdata, mysetdata, myaggregate, myprintdata);

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig          *u = unitigs[ti];

    if (u == NULL)
      continue;

    int  l = u->getLength();
    int  c = (int)rintf(u->getCovStat());
    int  a = (int)rintf(u->getLocalArrivalRate() * 10000);

    if (c < -200000) {
      fprintf(stderr, "ti=%d c=%d\n", ti, c);
    }

    z.nsamples     = 1;
    z.sum_frags    = z.min_frags    = z.max_frags    = u->getNumFrags();
    z.sum_bp       = z.min_bp       = z.max_bp       = l;;
    z.sum_rho      = z.min_rho      = z.max_rho      = (int)rintf(u->getAvgRho());
    z.sum_discr    = z.min_discr    = z.max_discr    = c;
    z.sum_arrival  = z.min_arrival  = z.max_arrival  = a;
    z.sum_rs_frags = z.min_rs_frags = z.max_rs_frags = 0;
    z.sum_nr_frags = z.min_nr_frags = z.max_nr_frags = 0;

    add_to_histogram(len, l, &z);
    add_to_histogram(cvg, c, &z);
    add_to_histogram(arv, a, &z);
  }

  fprintf(stats, "\n");
  fprintf(stats, "\n");
  fprintf(stats, "Unitig Length\n");
  fprintf(stats, "label\tsum\tcummulative\tcummulative   min  average  max\n");
  fprintf(stats, "     \t   \t sum       \t fraction\n");
  print_histogram(stats,len, 0, 1);

  fprintf(stats, "\n");
  fprintf(stats, "\n");
  fprintf(stats, "Unitig Coverage Stat\n");
  fprintf(stats, "label\tsum\tcummulative\tcummulative   min  average  max\n");
  fprintf(stats, "     \t   \t sum       \t fraction\n");
  print_histogram(stats,cvg, 0, 1);

  fprintf(stats, "\n");
  fprintf(stats, "\n");
  fprintf(stats, "Unitig Arrival Rate\n");
  fprintf(stats, "label\tsum\tcummulative\tcummulative   min  average  max\n");
  fprintf(stats, "     \t   \t sum       \t fraction\n");
  print_histogram(stats,arv, 0, 1);

  free_histogram(len);
  free_histogram(cvg);
  free_histogram(arv);

  if (stats != logFile)
    fclose(stats);
}
