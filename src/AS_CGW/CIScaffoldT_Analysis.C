
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

static const char *rcsid = "$Id: CIScaffoldT_Analysis.C,v 1.1 2012-07-16 09:13:17 brianwalenz Exp $";

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

  fprintf(stderr, "LIB %d %.2f+-%.2f [%.2f,%.2f] innie=%d sum=%f\n",
          iid, mean, stddev, minDist, maxDist, innie, s);
}




void
instrumentSCF::init(CIScaffoldT *scaffold, bool forward, double offMean, double offVariance) {

  iid = 0;

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
      //fprintf(stderr, "REVERSE to %.2f,%.2f\n", bgn.mean, end.mean);
    }

    bgn.mean     += offMean;
    bgn.variance += offVariance;

    end.mean     += offMean;
    end.variance += offVariance;

    instrumentCTG  ctg(contig->id, scaffold->id, bgn, end);

    //fprintf(stderr, "CONTIG %d at %.0f,%.0f\n", contig->id, contig->offsetAEnd.mean, contig->offsetBEnd.mean);

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

      FRGmap[imp->ident] = FRG.size();

      FRG.push_back(instrumentFRG(imp->ident, contig->id, ctg.fwd, ctg.bgn, ctg.end, imp->position));
    }

    CTG.push_back(ctg);
  }
}



void
instrumentSCF::init(CIScaffoldT *scaffoldA, SEdgeT *edge, CIScaffoldT *scaffoldB) {
  bool      fwdA = (edge->orient.isAB_AB() || edge->orient.isAB_BA());
  bool      fwdB = (edge->orient.isAB_AB() || edge->orient.isBA_AB());

  //fprintf(stderr, "EDGE: %d %d %.2f +- %.2f orient %c (fwd %d %d) scfLen %.2f %.2f\n",
  //        edge->idA, edge->idB,
  //        edge->distance.mean, sqrt(edge->distance.variance), edge->orient.toLetter(),
  //        fwdA, fwdB,
  //        scaffoldA->bpLength.mean, scaffoldB->bpLength.mean);

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

  //fprintf(stderr, "A fwd %d offset %f,%f\n", fwdA, offA.mean, offA.variance);
  //fprintf(stderr, "B fwd %d offset %f,%f\n", fwdB, offB.mean, offB.variance);

  init(scaffoldA, fwdA, offA.mean, offA.variance);
  init(scaffoldB, fwdB, offB.mean, offB.variance);

  iid = 0;
}




void
instrumentSCF::analyze(vector<instrumentLIB> &libs) {
  AS_IID  fiid;
  AS_IID  miid;

  char   *scfGap = NULL;
  uint32  scfLen = 0;

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

      //  If needed, build up the gap map for the scaffold.  Why the +1?  Because we're comparing
      //  int vs float, and the rules for MAX() seem to differ from less than.  MAX(2, 2.5) is 2 in
      //  integer space.  But 2 < 2.5 in float space.  PITA.
      //
      if (scfLen == 0) {
        for (uint32 i=0; i<CTG.size(); i++) {
          assert(CTG[i].bgn.mean < CTG[i].end.mean);
          scfLen = MAX(scfLen, CTG[i].end.mean + 1);
        }

        scfGap = new char [scfLen];

        memset(scfGap, 0, sizeof(char) * scfLen);

        for (uint32 i=0; i<CTG.size(); i++)
          for (uint32 p=CTG[i].bgn.mean; p<CTG[i].end.mean; p++) {
            assert(p < scfLen);
            scfGap[p] = true;
          }
      }

      //  Now sum over the PDF.

      instrumentFRG   &frg  = FRG[fi];

      int32   scfPos = 0;  //  MUST be signed

      if        (((lib.innie == true)  && (frg.fwd == true)) ||
                 ((lib.innie == false) && (frg.fwd == false))) {
        scfPos = frg.bgn + lib.minDist;

      } else if (((lib.innie == true)  && (frg.fwd == false)) ||
                 ((lib.innie == false) && (frg.fwd == true))) {
        scfPos = frg.end - lib.maxDist;

      } else {
        assert(0);
      }

      for (uint32 ii=0; ii<lib.pdf.size(); ii++, scfPos++) {
        if      (scfPos <  0)            numExternal += lib.pdf[ii];
        else if (scfLen <= scfPos)       numExternal += lib.pdf[ii];
        else if (scfGap[scfPos] == true) numGap      += lib.pdf[ii];
        else                             numMissing  += lib.pdf[ii];
      }

      //fprintf(stderr, "EXTERNAL frg %d (%u,%u ori %d) at scfPos %d,%d external %.4f gap %.4f missing %.4f\n",
      //        fiid,
      //        frg.bgn, frg.end, frg.fwd,
      //        scfPos, scfPos + lib.pdf.size(),
      //        numExternal, numGap, numMissing);

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

    uint32           dist      = 0;
    bool             misOrient = false;

    if (frg.bgn < mrg.bgn) {
      //  frg before mrg.
      assert(mrg.end > frg.bgn);
      dist = mrg.end - frg.bgn;

      if ((lib.innie == true)  && ((frg.fwd == false) || (mrg.fwd == true)))   misOrient = true;
      if ((lib.innie == false) && ((frg.fwd == true)  || (mrg.fwd == false)))  misOrient = true;

    } else {
      //  mrg before frg
      assert(frg.end > mrg.bgn);
      dist = frg.end - mrg.bgn;

      if ((lib.innie == true)  && ((frg.fwd == true)  || (mrg.fwd == false)))  misOrient = true;
      if ((lib.innie == false) && ((frg.fwd == false) || (mrg.fwd == true)))   misOrient = true;
    }

    if (misOrient == false) {
      if      (dist < lib.minDist)
        numTooClose++;
      else if (lib.maxDist < dist)
        numTooFar++;
      else
        numHappy++;

      continue;
    }

    if      (dist < lib.minDist) {
      numMisClose++;
    } else if (lib.maxDist < dist) {
      numMisFar++;
    } else {
      numMis++;
    }
  }

  delete [] scfGap;

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



