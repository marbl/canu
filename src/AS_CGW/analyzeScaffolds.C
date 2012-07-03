
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

const char *mainid = "$Id: analyzeScaffolds.C,v 1.1 2012-07-03 01:32:56 brianwalenz Exp $";

#include "AS_global.h"

#include "Globals_CGW.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"

#include <vector>
#include <map>

using namespace std;



//  This came from analyzePosMap.C, class libEntry.

class instrumentLIB {
public:
  instrumentLIB() {
    iid    = 0;
    mean   = 0;
    stddev = 0;
    innie  = true;
  };

  instrumentLIB(AS_IID iid_, double mean_, double stddev_, bool innie_) {
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

    for (int32 i=0; i<=6 * stddev; i++) {
      double x = i - 3 * stddev;

      pdf[i] = c * exp(-(x * x / d));
      s += pdf[i];
    }

    fprintf(stderr, "LIB %d %.2f+-%.2f [%.2f,%.2f] innie=%d sum=%f\n",
            iid, mean, stddev, minDist, maxDist, innie, s);
  };

  ~instrumentLIB() {
  };

  AS_IID           iid;

  double           mean;
  double           stddev;
  bool             innie;

  double           minDist;
  double           maxDist;

  vector<double>   pdf;
};




class instrumentFRG {
public:
  instrumentFRG();
  instrumentFRG(AS_IID       iid_,
                AS_IID       cid_,
                bool         ctgfwd_,
                LengthT      ctgbgn_,
                LengthT      ctgend_,
                SeqInterval  pos_) {

    uint32   rbgn = pos_.bgn;
    uint32   rend = pos_.end;

    //  If the contig is flipped, reverse complement the fragment positions.

    if (ctgfwd_ == false) {
      double   ctgLen = ctgend_.mean - ctgbgn_.mean;

      rbgn = ctgLen - pos_.end;
      rend = ctgLen - pos_.bgn;
    }

    //  Offset to scaffold positions.

    rbgn += ctgbgn_.mean;
    rend += ctgbgn_.mean;

    //  And set our variables, noting if the fragment is flipped.

    iid  = iid_;
    cid  = cid_;

    if (rbgn < rend) {
      fwd  = (ctgfwd_) ? true : false;
      bgn  = rbgn;
      end  = rend;
    } else {
      fwd  = (ctgfwd_) ? false : true;
      bgn  = rend;
      end  = rbgn;
    }

    //fprintf(stderr, "FRG %d at orig %u,%u ctg %u,%u scf %u,%u fwd=%d\n",
    //        iid, pos_.bgn, pos_.end, rbgn, rend, bgn, end, fwd);
  };

  AS_IID                   iid;  //  My IID
  AS_IID                   cid;  //  My contig IID
  bool                     fwd;  //  Orientation in scaffold
  uint32                   bgn;  //  Position in scaffold
  uint32                   end;
};


class instrumentCTG {
public:
  instrumentCTG();
  instrumentCTG(AS_IID       iid_,
                AS_IID       sid_,
                LengthT      bgn_,
                LengthT      end_) {
    iid = iid_;
    sid = sid_;

    if (bgn_.mean < end_.mean) {
      fwd = true;
      bgn = bgn_;
      end = end_;
    } else {
      fwd = false;
      bgn = end_;
      end = bgn_;
    }

    fprintf(stderr, "CONTIG %d at %.0f,%.0f fwd=%d\n",
            iid, bgn.mean, end.mean, fwd);
  };

  AS_IID                   iid;  //  My IID
  AS_IID                   sid;  //  My scaffold IID
  bool                     fwd;  //  Orientation in scaffold
  LengthT                  bgn;  //  Position in scaffold
  LengthT                  end;
};


class instrumentSCF {
public:
  instrumentSCF() {
    init(NULL);
  };
  instrumentSCF(CIScaffoldT *scaffold) {
    init(scaffold);
  };

  void    init(CIScaffoldT *scaffold);

  void    analyze(vector<instrumentLIB> &lib);
  void    report(void);

  AS_IID                   iid;

  vector<instrumentCTG>    CTG;
  vector<instrumentFRG>    FRG;

  map<AS_IID,size_t>       FRGmap;        //  Maps fragment IID to location in FRG vector

  //  Should keep this per library?

  //  Summary
  uint32                   numMateInternal;
  uint32                   numMateExternal;
  uint32                   numMateFree;

  double                   numEcstatic;
  double                   numDejected;
  double                   numApathetic;

  //  Ecstatic
  double                   numHappy;      //  good orient, good distance
  double                   numGap;        //  not present in scaffold, but in a gap

  //  Dejected
  double                   numMisClose;   //  misoriented and too close
  double                   numMis;        //  misoriented
  double                   numMisFar;     //  misoriented and too far

  double                   numTooClose;   //  oriented, but too close
  //                       numHappy       //  oriented, correct distance == numHappy
  double                   numTooFar;     //  oriented, but too far

  double                   numMissing;    //  not present in scaffold, but enough space for it

  //  Apathetic
  double                   numExternal;   //  not present in scaffold, off the end
};





void
instrumentSCF::init(CIScaffoldT *scaffold) {

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

    instrumentCTG  ctg(contig->id, scaffold->id, contig->offsetAEnd, contig->offsetBEnd);

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

      instrumentFRG   &frg  = FRG[fi];

      double  fracGap      = 0.0;
      double  fracMissing  = 0.0;
      double  fracExternal = 0.0;

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

      if (scfLen == 0) {
        scfLen = CTG[CTG.size()-1].end.mean;
        scfGap = new char [scfLen];

        memset(scfGap, 0, sizeof(char) * scfLen);

        for (uint32 i=0; i<CTG.size(); i++)
          for (uint32 p=CTG[i].bgn.mean; p<CTG[i].end.mean; p++)
            scfGap[p] = true;
      }

      for (uint32 ii=0; ii<lib.pdf.size(); ii++, scfPos++) {
        if        (scfPos < 0)
          fracExternal += lib.pdf[ii];

        else if (scfLen < scfPos)
          fracExternal += lib.pdf[ii];

        else if (scfGap[scfPos] == true)
          fracGap += lib.pdf[ii];

        else
          fracMissing += lib.pdf[ii];
      }

      numGap      += fracGap;       //  Could be placed in a gap
      numMissing  += fracMissing;   //  Should be in the scaffold
      numExternal += fracExternal;  //  Off the end of the scaffold
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

    if      (dist < lib.minDist)
      numMisClose++;
    else if (lib.maxDist < dist)
      numMisFar++;
    else
      numMis++;
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







int main (int argc, char *argv[]) {
  int32       checkpointVers           = 0;
  int32       tigStoreVers             = 0;

  GlobalData = new Globals_CGW();

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      strcpy(GlobalData->gkpStoreName, argv[++arg]);

    } else if (strcmp(argv[arg], "-t") == 0) {
      strcpy(GlobalData->tigStoreName, argv[++arg]);
      tigStoreVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-c") == 0) {
      strcpy(GlobalData->outputPrefix, argv[++arg]);
      checkpointVers = atoi(argv[++arg]);

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((GlobalData->gkpStoreName[0] == 0) ||
      (GlobalData->tigStoreName[0] == 0) ||
      (err)) {
    fprintf(stderr, "usage: %s -g gkpStore [-o prefix] [-s firstUID] [-n namespace] [-E server] [-h]\n", argv[0]);
    fprintf(stderr, "  -g gkpStore             mandatory path to the gkpStore\n");
    fprintf(stderr, "  -t tigStore version     mandatory path to the tigStore and version\n");
    fprintf(stderr, "  -c checkpoint version   optional path to a checkpoint and version\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  LoadScaffoldGraphFromCheckpoint(GlobalData->outputPrefix, checkpointVers, FALSE);


  vector<instrumentLIB>   lib;

  for (int32 i=0; i<GetNumDistTs(ScaffoldGraph->Dists); i++) {
    DistT *dptr = GetDistT(ScaffoldGraph->Dists, i);

    lib.push_back(instrumentLIB(i, dptr->mu, dptr->sigma, true));
  }

  GraphNodeIterator   scaffolds;
  CIScaffoldT        *scaffold;

  InitGraphNodeIterator(&scaffolds, ScaffoldGraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while ((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL) {
    if(scaffold->type != REAL_SCAFFOLD)
      continue;

    //if (scaffold->id != 14)
    //  continue;

    instrumentSCF  scf(scaffold);

    scf.analyze(lib);
    scf.report();
  }


  DestroyScaffoldGraph(ScaffoldGraph);

  return(0);
}
