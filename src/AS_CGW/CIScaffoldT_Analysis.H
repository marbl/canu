
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

#ifndef CISCAFFOLDT_ANALYSIS_H
#define CISCAFFOLDT_ANALYSIS_H

static const char *rcsid_CISCAFFOLDT_ANALYSIS_H = "$Id: CIScaffoldT_Analysis.H,v 1.1 2012-07-16 09:13:17 brianwalenz Exp $";

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
  instrumentLIB(AS_IID iid_, double mean_, double stddev_, bool innie_);
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
  instrumentFRG() {
    iid = 0;
    cid = 0;
    fwd = false;
    bgn = 0;
    end = 0;
  };
  instrumentFRG(AS_IID       iid_,
                AS_IID       cid_,
                bool         ctgfwd_,
                LengthT      ctgbgn_,
                LengthT      ctgend_,
                SeqInterval  pos_) {

    uint32   rbgn = pos_.bgn;
    uint32   rend = pos_.end;

    //  If the contig is flipped, reverse complement the fragment positions.
    //  Yes, we really do bodge contig length.  This can cause differences
    //  if you compare an instrument scaffold against the same scaffold with
    //  forward==false.

    if (ctgfwd_ == false) {
      double   ctgLen = ctgend_.mean - ctgbgn_.mean;

      rbgn = (ctgLen >= pos_.end) ? ctgLen - pos_.end : 0;
      rend = (ctgLen >= pos_.bgn) ? ctgLen - pos_.bgn : 0;
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

    //fprintf(stderr, "instrumentFRG %d at orig %u,%u ctg %u,%u scf %u,%u fwd=%d\n",
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
  instrumentCTG() {
    iid = 0;
    sid = 0;
    fwd = false;
    bgn.mean = 0;
    bgn.variance = 0;
    end.mean = 0;
    end.variance = 0;
  };
  instrumentCTG(AS_IID iid_, AS_IID sid_, LengthT bgn_, LengthT end_) {
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

    //fprintf(stderr, "instrumentCTG %d orig %.2f,%.2f now at %.0f,%.0f fwd=%d\n",
    //        iid, bgn_.mean, end_.mean, bgn.mean, end.mean, fwd);
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
    init(NULL, false, 0, 0);
  };
  instrumentSCF(CIScaffoldT *scaffold) {
    init(scaffold, true, 0, 0);
  };
  instrumentSCF(CIScaffoldT *scaffold, bool forward) {
    init(scaffold, forward, 0, 0);
  };
  instrumentSCF(CIScaffoldT *scaffoldA, SEdgeT *edge, CIScaffoldT *scaffoldB) {
    init(scaffoldA, edge, scaffoldB);
  };

  void    init(CIScaffoldT *scaffold, bool forward, double offMean, double offVariance);
  void    init(CIScaffoldT *scaffoldA, SEdgeT *edge, CIScaffoldT *scaffoldB);

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


#endif  //  CISCAFFOLDT_ANALYSIS_H
