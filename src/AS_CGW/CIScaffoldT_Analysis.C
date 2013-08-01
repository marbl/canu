
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * Copyright (C) 2007, J. Craig Venter Institute.
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

static const char *rcsid = "$Id$";

#include "CIScaffoldT_Analysis.H"

instrumentLIB::instrumentLIB(AS_IID iid_, double mean_, double stddev_, bool innie_) {
  iid    = iid_;
  mean   = mean_;
  stddev = stddev_;
  innie  = innie_;

  minDist = mean - 3.0 * stddev;
  maxDist = mean + 3.0 * stddev;

  pdf.resize(6 * stddev + 1);

  double  c = 1 / (stddev * sqrt(2 * M_PI));
  double  d = 2 * stddev * stddev;

  double  s = 0;

  //  Compute a fudge scaling factor because we aren't computing this over infinity.
  for (int32 i=0; i<=6 * stddev; i++)
    s += c * exp(-((i - 3 * stddev) * (i - 3 * stddev) / d));

  c = c / s;

  for (int32 i=0; i<=6 * stddev; i++)
    pdf[i] = c * exp(-((i - 3 * stddev) * (i - 3 * stddev) / d));

  s = 0.0;

  for (int32 i=0; i<=6 * stddev; i++)
    s += pdf[i];

  fprintf(stderr, "instrumentLIB %d %.2f+-%.2f [%.2f,%.2f] innie=%d sum=%f\n",
          iid, mean, stddev, minDist, maxDist, innie, s);
}



void
instrumentSCF::clearStats(void) {
  numMateInternal = 0;
  numMateExternal = 0;
  numMateFree     = 0;

  numEcstatic   = 0.0;
  numDejected   = 0.0;
  numApathetic  = 0.0;
  numHappy      = 0.0;
  numGap        = 0.0;
  numMisClose   = 0.0;
  numMis        = 0.0;
  numMisFar     = 0.0;
  numTooClose   = 0.0;
  numTooFar     = 0.0;
  numMissing    = 0.0;
  numExternal   = 0.0;
}

void
instrumentSCF::init(CIScaffoldT *scaffold, bool forward, double offMean, double offVariance) {

  iid = 0;

  clearStats();

  if (scaffold == NULL)
    return;

  assert(scaffold->type == REAL_SCAFFOLD);
  assert(scaffold->info.Scaffold.numElements > 0);

  iid    = scaffold->id;

  CIScaffoldTIterator      contigs;
  NodeCGW_T               *contig;

  InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &contigs);

  bool      fwd = false;
  LengthT   bgn = { 0, 0 };
  LengthT   end = { 0, 0 };

  while ((contig = NextCIScaffoldTIterator(&contigs)) != NULL) {
    assert(contig->scaffoldID == scaffold->id);

    bgn = contig->offsetAEnd;
    end = contig->offsetBEnd;

    //  This is not the normal reverse complement.  We're encoding the orientation
    //  in coordinates, in effect, this is not a span on the scaffold, it's two
    //  independent points.
    //    bgn end                  bgn end
    //    0   10 (fwd)   becomes   100 90  (rev)
    //    10  0  (rev)   becomes   90  100 (fwd)
    //
    if (forward == false) {
      bgn.mean     = scaffold->bpLength.mean     - contig->offsetAEnd.mean;
      bgn.variance = scaffold->bpLength.variance - contig->offsetAEnd.variance;
      end.mean     = scaffold->bpLength.mean     - contig->offsetBEnd.mean;
      end.variance = scaffold->bpLength.variance - contig->offsetBEnd.variance;
#ifdef CIA_SUPER_VERBOSE
      fprintf(stderr, "CONTIG %d rv %.0f,%.0f\n", contig->id, bgn.mean, end.mean);
#endif
    }

    bgn.mean     += offMean;
    bgn.variance += offVariance;

    end.mean     += offMean;
    end.variance += offVariance;

    instrumentCTG  ctg(contig->id, scaffold->id, bgn, end);

#ifdef CIA_SUPER_VERBOSE
    fprintf(stderr, "CONTIG %d at %.0f,%.0f\n", contig->id, contig->offsetAEnd.mean, contig->offsetBEnd.mean);
#endif

#if 0
    NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph, contig->info.Contig.AEndCI);

    if ((ScaffoldGraph->tigStore->getNumUnitigs(contig->id, FALSE) == 1) &&
        (contig->scaffoldID == NULLINDEX) &&
        (unitig->info.CI.numInstances > 0))
      //  Contig is a surrogate instance
      continue;
#endif

    MultiAlignT *ma = ScaffoldGraph->tigStore->loadMultiAlign(contig->id, FALSE);

    for(int32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
      IntMultiPos *imp = GetIntMultiPos(ma->f_list, i);
      CIFragT     *cif = GetCIFragT(ScaffoldGraph->CIFrags, imp->ident);

      FRGmap[imp->ident] = FRG.size();

      FRG.push_back(instrumentFRG(imp->ident, contig->id, ctg.fwd, ctg.bgn, ctg.end, cif->contigOffset5p, cif->contigOffset3p));
    }

    CTG.push_back(ctg);
  }
}



void
instrumentSCF::init(CIScaffoldT *scaffoldA, SEdgeT *edge, CIScaffoldT *scaffoldB) {
  bool      fwdA = (edge->orient.isAB_AB() || edge->orient.isAB_BA());
  bool      fwdB = (edge->orient.isAB_AB() || edge->orient.isBA_AB());

#ifdef CIA_SUPER_VERBOSE
  fprintf(stderr, "EDGE: %d %d %.2f +- %.2f orient %c (fwd %d %d) scfLen %.2f %.2f\n",
          edge->idA, edge->idB,
          edge->distance.mean, sqrt(edge->distance.variance), edge->orient.toLetter(),
          fwdA, fwdB,
          scaffoldA->bpLength.mean, scaffoldB->bpLength.mean);
#endif

  LengthT   offA = { 0, 0 };
  LengthT   offB = { 0, 0 };

  if        (edge->distance.mean >= 0) {
    offB.mean      = scaffoldA->bpLength.mean     + edge->distance.mean;
    offB.variance  = scaffoldA->bpLength.variance + edge->distance.variance;

  } else if (scaffoldA->bpLength.mean < -edge->distance.mean) {
    offA.mean      = -edge->distance.mean - scaffoldA->bpLength.mean;
    offA.variance  =  edge->distance.variance;

  } else {
    offB.mean      = scaffoldA->bpLength.mean     - edge->distance.mean;
    offB.variance  = scaffoldA->bpLength.variance - edge->distance.variance;
  }

#ifdef CIA_SUPER_VERBOSE
  fprintf(stderr, "A fwd %d offset %f,%f\n", fwdA, offA.mean, offA.variance);
  fprintf(stderr, "B fwd %d offset %f,%f\n", fwdB, offB.mean, offB.variance);
#endif

  init(scaffoldA, fwdA, offA.mean, offA.variance);
  init(scaffoldB, fwdB, offB.mean, offB.variance);

  iid = 0;
}




void
instrumentSCF::analyze(vector<instrumentLIB> &libs) {
  AS_IID  fiid;
  AS_IID  miid;

  char   *scfCov = NULL;
  uint32  scfLen = 0;

  clearStats();

  for (uint32 fi=0; fi<FRG.size(); fi++) {
    CIFragT  *cif = GetCIFragT(ScaffoldGraph->CIFrags, FRG[fi].iid);
    CIFragT  *mif = GetCIFragT(ScaffoldGraph->CIFrags, cif->mate_iid);

    if (cif->mate_iid == 0) {
      numMateFree++;
      continue;
    }

    fiid = cif->read_iid;
    miid = mif->read_iid;

    assert(cif->mate_iid == mif->read_iid);
    assert(cif->read_iid == mif->mate_iid);

    assert(cif->flags.bits.isDeleted == false);
    assert(mif->flags.bits.isDeleted == false);
      
    assert(cif->read_iid == fiid);
    assert(mif->read_iid == miid);

    assert(cif->dist == mif->dist);

    instrumentLIB   &lib = libs[cif->dist];

    //  If external...
    if (FRGmap.find(miid) == FRGmap.end()) {
      numMateExternal++;

      //  If needed, build up the coverage map for the scaffold.  If set, this position is covered
      //  by a contig.
      //
      //  Why the +1?  Because we're comparing int vs float, and the rules for MAX() seem to differ
      //  from less than.  MAX(2, 2.5) is 2 in integer space.  But 2 < 2.5 in float space.  PITA.
      //
      if (scfLen == 0) {
        for (uint32 i=0; i<CTG.size(); i++) {
          assert(CTG[i].bgn.mean < CTG[i].end.mean);
          scfLen = MAX(scfLen, CTG[i].end.mean + 1);
        }

        scfCov = new char [scfLen];

        memset(scfCov, 0, sizeof(char) * scfLen);

        for (uint32 i=0; i<CTG.size(); i++)
          for (uint32 p=CTG[i].bgn.mean; p<CTG[i].end.mean; p++) {
            assert(p < scfLen);
            scfCov[p] = true;
          }
      }

      //  Now sum over the PDF.

      instrumentFRG   &frg  = FRG[fi];

      int32   scfPos = 0;  //  MUST be signed

      if        (((lib.innie == true)  && (frg.fwd == true)) ||
                 ((lib.innie == false) && (frg.fwd == false))) {
        scfPos = (frg.bgn + frg.end) / 2 + lib.mean - lib.pdf.size() / 2;

      } else if (((lib.innie == true)  && (frg.fwd == false)) ||
                 ((lib.innie == false) && (frg.fwd == true))) {
        scfPos = (frg.bgn + frg.end) / 2 - lib.mean - lib.pdf.size() / 2;

      } else {
        assert(0);
      }

      for (uint32 ii=0; ii<lib.pdf.size(); ii++, scfPos++) {
        if      (scfPos <  0)             numExternal += lib.pdf[ii];
        else if (scfLen <= scfPos)        numExternal += lib.pdf[ii];
        else if (scfCov[scfPos] == false) numGap      += lib.pdf[ii];
        else                              numMissing  += lib.pdf[ii];
      }

#ifdef CIA_SUPER_VERBOSE
      fprintf(stderr, "EXTERNAL frg %d (%u,%u fwd %d) mate range %d,%d - external %.4f gap %.4f missing %.4f\n",
              fiid,
              frg.bgn, frg.end, frg.fwd,
              scfPos - (int32)lib.pdf.size(), scfPos,
              numExternal, numGap, numMissing);
#endif

      continue;
    }

    //  Mated, and both fragments are in the scaffold.

    if (miid < fiid)
      continue;

    numMateInternal++;

    instrumentFRG   &frg  = FRG[fi];
    instrumentFRG   &mrg  = FRG[FRGmap[miid]];

    assert(frg.iid == fiid);
    assert(mrg.iid == miid);

    int32  frgmin = (frg.bgn < frg.end) ? frg.bgn : frg.end;
    int32  frgmax = (frg.bgn < frg.end) ? frg.end : frg.bgn;

    int32  mrgmin = (mrg.bgn < mrg.end) ? mrg.bgn : mrg.end;
    int32  mrgmax = (mrg.bgn < mrg.end) ? mrg.end : mrg.bgn;

    int32  posmin = (frgmin < mrgmin) ? frgmin : mrgmin;
    int32  posmax = (frgmax < mrgmax) ? mrgmax : frgmax;

    int32  dist      = posmax - posmin;
    bool   misOrient = false;

    if        ((frgmin <= mrgmin) && (mrgmax <= frgmax)) {
      //  mrg contained in frg
      misOrient = true;

    } else if ((mrgmin <= frgmin) && (frgmax <= mrgmax)) {
      //  frg contained in mrg
      misOrient = true;
        
    } else if (frgmin < mrgmin) {
      //  frg before mrg.
      assert(frgmax < mrgmax);

      if ((lib.innie == true)  && ((frg.fwd == false) || (mrg.fwd == true)))   misOrient = true;
      if ((lib.innie == false) && ((frg.fwd == true)  || (mrg.fwd == false)))  misOrient = true;

    } else {
      //  mrg before frg
      assert(mrgmax < frgmax);

      if ((lib.innie == true)  && ((frg.fwd == true)  || (mrg.fwd == false)))  misOrient = true;
      if ((lib.innie == false) && ((frg.fwd == false) || (mrg.fwd == true)))   misOrient = true;
    }


    if (misOrient == false) {
      if      (dist < lib.minDist) {
#ifdef CIA_SUPER_VERBOSE
        fprintf(stderr, "ORICLOSE frg %d (%u,%u fwd %d) mrg %d (%u,%u fwd %d)\n",
                fiid, frg.bgn, frg.end, frg.fwd,
                miid, mrg.bgn, mrg.end, mrg.fwd);
#endif
        numTooClose++;
      } else if (lib.maxDist < dist) {
#ifdef CIA_SUPER_VERBOSE
        fprintf(stderr, "ORIFAR   frg %d (%u,%u fwd %d) mrg %d (%u,%u fwd %d)\n",
                fiid, frg.bgn, frg.end, frg.fwd,
                miid, mrg.bgn, mrg.end, mrg.fwd);
#endif
        numTooFar++;
      } else {
#ifdef CIA_SUPER_VERBOSE
        fprintf(stderr, "ORIHAPPY frg %d (%u,%u fwd %d) mrg %d (%u,%u fwd %d)\n",
                fiid, frg.bgn, frg.end, frg.fwd,
                miid, mrg.bgn, mrg.end, mrg.fwd);
#endif
        numHappy++;
      }

      continue;
    }

    if      (dist < lib.minDist) {
#ifdef CIA_SUPER_VERBOSE
      fprintf(stderr, "MISCLOSE frg %d (%u,%u fwd %d) mrg %d (%u,%u fwd %d)\n",
              fiid, frg.bgn, frg.end, frg.fwd,
                miid, mrg.bgn, mrg.end, mrg.fwd);
#endif
      numMisClose++;
    } else if (lib.maxDist < dist) {
#ifdef CIA_SUPER_VERBOSE
      fprintf(stderr, "MISFAR   frg %d (%u,%u fwd %d) mrg %d (%u,%u fwd %d)\n",
              fiid, frg.bgn, frg.end, frg.fwd,
              miid, mrg.bgn, mrg.end, mrg.fwd);
#endif
      numMisFar++;
    } else {
#ifdef CIA_SUPER_VERBOSE
      fprintf(stderr, "MIS      frg %d (%u,%u fwd %d) mrg %d (%u,%u fwd %d)\n",
                fiid, frg.bgn, frg.end, frg.fwd,
                miid, mrg.bgn, mrg.end, mrg.fwd);
#endif
      numMis++;
    }
  }

  delete [] scfCov;

  numEcstatic  = numHappy + numGap;
  numDejected  = numMisClose + numMis + numMisFar + numTooClose + numTooFar + numMissing;
  numApathetic = numExternal;
}



void
instrumentSCF::report(void) {
  fprintf(stderr, "Instrumented scaffold "F_IID" with "F_SIZE_T" contigs and "F_SIZE_T" fragments.\n", iid, CTG.size(), FRG.size());
  fprintf(stderr, "  mates:     "F_U32"\n",          numMateInternal + numMateExternal);
  fprintf(stderr, "             "F_U32" internal\n", numMateInternal);
  fprintf(stderr, "             "F_U32" external\n", numMateExternal);
  fprintf(stderr, "             "F_U32" unmated\n",  numMateFree);
  fprintf(stderr, "\n");
  fprintf(stderr, "  ecstatic:  %.2f\n",       numEcstatic);
  fprintf(stderr, "             %.2f happy\n", numHappy);
  fprintf(stderr, "             %.2f gap\n",   numGap);
  fprintf(stderr, "\n");
  fprintf(stderr, "  dejected:  %.2f\n", numDejected);
  fprintf(stderr, "             %.2f misoriented, close\n", numMisClose);
  fprintf(stderr, "             %.2f misoriented\n",        numMis);
  fprintf(stderr, "             %.2f misoriented, far\n",   numMisFar);
  fprintf(stderr, "             %.2f oriented, close\n",    numTooClose);
  fprintf(stderr, "             %.2f oriented, far\n",      numTooFar);
  fprintf(stderr, "             %.2f missing\n",            numMissing);
  fprintf(stderr, "\n");
  fprintf(stderr, "  apathetic: %.2f\n",          numApathetic);
  fprintf(stderr, "             %.2f external\n", numExternal);
}




//  If we're just after mate happiness during merging (NOT ENABLED) we are allowed to rearrange
//  contigs willy nilly.
//
//  If we're after a replacement for least squares, we cannot rearrange.  Worse, we need someway of
//  ensuring that the gaps we compute don't implicitly rearrange for us.
//
//  This DOES NOT WORK, especially in heavily interleaved scaffolds.  It tries to compute gap by
//  gap, using mates that span a specific gap.  If those mates span other gaps, and those gap sizes
//  are incorrect, then the estimate for this gap size is incorrect.

void
instrumentSCF::estimateGaps(vector<instrumentLIB> &libs, bool allowReorder) {

  //  Initialize the GAP size vector.  Set variance to zero; it will be computed later.

  if (allowReorder)
    sort(CTG.begin(), CTG.end());

#ifdef CIA_SUPER_VERBOSE
  for (uint32 i=0; i<CTG.size(); i++)
    fprintf(stderr, "CTG %d %7.0f-%7.0f\n", CTG[i].iid, CTG[i].bgn.mean, CTG[i].end.mean);
#endif

  for (uint32 i=1; i<CTG.size(); i++) {
    GAPmean.push_back(CTG[i].bgn.mean - CTG[i-1].end.mean);
    GAPvari.push_back(0.0);
  }

  //  Build a map from CTG iid to CTG vector element.

  for (uint32 i=0; i<CTG.size(); i++)
    CTGmap[CTG[i].iid] = i;

  //  Without trying to be clever or efficient, for each gap, iterate over all fragments.
  //  if the pair spans the gap, save it for the estimate.

  for (uint32 it=0; it<1; it++) {
    for (uint32 gi=0; gi<GAPmean.size(); gi++) {
      vector<int32>   actualM;  //  Current actual mean distance
      vector<int32>   expectM;  //  Expected mean distance
      vector<int32>   expectS;  //  Expected std.dev

      for (uint32 fi=0; fi<FRG.size(); fi++) {
        CIFragT  *cif = GetCIFragT(ScaffoldGraph->CIFrags, FRG[fi].iid);
        CIFragT  *mif = GetCIFragT(ScaffoldGraph->CIFrags, cif->mate_iid);

        if (cif->mate_iid == 0) {
          numMateFree++;
          continue;
        }

        AS_IID fiid = cif->read_iid;
        AS_IID miid = mif->read_iid;

        assert(cif->mate_iid == mif->read_iid);
        assert(cif->read_iid == mif->mate_iid);

        assert(cif->flags.bits.isDeleted == false);
        assert(mif->flags.bits.isDeleted == false);
      
        assert(cif->read_iid == fiid);
        assert(mif->read_iid == miid);

        assert(cif->dist == mif->dist);

        instrumentLIB   &lib = libs[cif->dist];

        if (FRGmap.find(miid) == FRGmap.end())
          //  Externally mated.
          continue;

        instrumentFRG   &frg  = FRG[fi];
        instrumentFRG   &mrg  = FRG[FRGmap[miid]];

        if (mrg.bgn < frg.bgn)
          //  Process it later in canonical form
          continue;

        if ((frg.fwd != true) || (mrg.fwd != false))
          //  Bad orientation.
          continue;

        if (frg.end >= CTG[gi].end.mean) 
          //  Left most fragment is after the gap begins.
          continue;

        if (mrg.bgn <= CTG[gi+1].bgn.mean)
          //  Right most fragment is before the gap ends.
          continue;

        int32  dist = mrg.end - frg.bgn;

        if (dist > lib.mean + 3 * lib.stddev)
          //  Pair is excessively stretched.
          continue;

        if (dist < lib.mean - 3 * lib.stddev)
          //  Pair is excessively compressed.
          continue;

        //  Pair should be spanning the gap.  Negative gaps screw up this test.  All we can assert
        //  is the negation of what we test above.
        //
        assert(CTG[gi].end.mean >= frg.bgn);
        assert(mrg.end          >= CTG[gi+1].bgn.mean);

        int32   lb = CTG[gi].end.mean - frg.bgn;
        int32   rb = mrg.end          - CTG[gi+1].bgn.mean;

#ifdef CIA_SUPER_VERBOSE
        fprintf(stderr, "GAP %d frg %d %d-%d mrg %d %d-%d actual %d / %d expected %.2f +- %2.f / %f\n",
                gi,
                frg.iid, frg.bgn, frg.end,
                mrg.iid, mrg.bgn, mrg.end,
                mrg.end - frg.bgn,
                mrg.end - frg.bgn - lb - rb,
                lib.mean, lib.stddev,
                lib.mean - lb - rb);
#endif

        actualM.push_back(mrg.end - frg.bgn - lb - rb);  //  Yeah, would have been easier as ctg[1].bgn - ctg[0].end

        expectM.push_back(lib.mean - lb - rb);
        expectS.push_back(lib.stddev);
      } //  Over all frags

      //  Compute some estimate of the gap

      double gapM = 0.0;
      double gapV = 0.0;

      if        (actualM.size() == 0) {
        gapM = CTG[gi+1].bgn.mean     - CTG[gi].end.mean;
        gapV = CTG[gi+1].bgn.variance - CTG[gi].end.variance;

      } else if (actualM.size() == 1) {
        gapM = expectM[0];
        gapV = expectS[0] * expectS[0];

      } else {
        double   sumDists   = 0.0;
        double   sumSquares = 0.0;

        for (uint32 i=0; i<expectM.size(); i++)
          sumDists += expectM[i];
    
        gapM = sumDists / expectM.size();

        for (uint32 i=0; i<expectM.size(); i++)
          sumSquares += (expectM[i] - gapM) * (expectM[i] - gapM);

        gapV = sumSquares / (expectM.size() - 1);
      }

      if ((allowReorder == false) && (gapM < -20))
        gapM = -20;

      GAPmean[gi] = gapM;
      GAPvari[gi] = gapV;
    }  //  Over all gaps

    //  All gap sizes are estimated.  Update CTG and FRG positions.

#warning this is inefficient

    for (uint32 fi=0; fi<FRG.size(); fi++) {
      uint32 ci = CTGmap[FRG[fi].cid];

      FRG[fi].bgn -= CTG[ci].bgn.mean;
      FRG[fi].end -= CTG[ci].bgn.mean;
    }

#ifdef CIA_SUPER_VERBOSE
    fprintf(stderr, "ITER %u\n", it);
    fprintf(stderr, "CTG %d from %7.0f +- %7.0f -- %7.0f +- %7.0f GAP %7.0f +- %7.0f --TO-- %7.0f +- %7.0f -- %7.0f +- %7.0f -- GAP %7.0f +- %7.0f\n",
            0,
            CTG[0].bgn.mean, sqrt(CTG[0].bgn.variance), CTG[0].end.mean, sqrt(CTG[0].end.variance),
            0.0, 0.0,
            CTG[0].bgn.mean, sqrt(CTG[0].bgn.variance), CTG[0].end.mean, sqrt(CTG[0].end.variance),
            0.0, 0.0);
#endif
    
    for (uint32 gi=0; gi<GAPmean.size(); gi++) {
      uint32  ca = gi;
      uint32  cb = gi + 1;
      int32   cm = CTG[cb].end.mean - CTG[cb].bgn.mean;
      int32   cv = CTG[cb].end.mean - CTG[cb].bgn.mean;

      LengthT  olda = CTG[ca].end;
      LengthT  oldb = CTG[cb].bgn;
      LengthT  olde = CTG[cb].end;

      CTG[cb].bgn.mean     = CTG[ca].end.mean     + GAPmean[gi];
      CTG[cb].bgn.variance = CTG[ca].end.variance + GAPvari[gi];

      CTG[cb].end.mean     = CTG[cb].bgn.mean     + cm;
      CTG[cb].end.variance = CTG[cb].bgn.variance + cv;

#ifdef CIA_SUPER_VERBOSE
      fprintf(stderr, "CTG %d from %7.0f +- %7.0f -- %7.0f +- %7.0f GAP %7.0f +- %7.0f --TO-- %7.0f +- %7.0f -- %7.0f +- %7.0f -- GAP %7.0f +- %7.0f\n",
              cb,
              oldb.mean, sqrt(oldb.variance), olde.mean, sqrt(olde.variance),
              oldb.mean     - olda.mean,
              sqrt(oldb.variance - olda.variance),
              CTG[cb].bgn.mean, sqrt(CTG[cb].bgn.variance), CTG[cb].end.mean, sqrt(CTG[cb].end.variance),
              CTG[cb].bgn.mean     - CTG[ca].end.mean,
              sqrt(CTG[cb].bgn.variance - CTG[ca].end.variance));
#endif
    }

    for (uint32 fi=0; fi<FRG.size(); fi++) {
      uint32 ci = CTGmap[FRG[fi].cid];

      FRG[fi].bgn += CTG[ci].bgn.mean;
      FRG[fi].end += CTG[ci].bgn.mean;
    }
  }
}

