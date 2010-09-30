
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

static const char *rcsid = "$Id: AS_BOG_MateLocation.cc,v 1.1 2010-09-30 05:40:21 brianwalenz Exp $";

#include "AS_BOG_MateLocation.hh"


MateLocation::MateLocation(Unitig *utg) {
  MateLocationEntry   mle;

  mle.mlePos1.bgn = mle.mlePos1.end = 0;
  mle.mlePos2.bgn = mle.mlePos2.end = 0;

  mle.mleFrgID1 = 0;
  mle.mleFrgID2 = 0;

  mle.mleUtgID1 = 0;
  mle.mleUtgID2 = 0;

  mle.isGrumpy = false;

  _table.clear();
  _table.push_back(mle);

  _tigLen = utg->getLength();

  goodGraph   = new int32 [_tigLen + 1];
  badFwdGraph = new int32 [_tigLen + 1];
  badRevGraph = new int32 [_tigLen + 1];

  badExternalFwd = new int32 [_tigLen + 1];
  badExternalRev = new int32 [_tigLen + 1];

  badCompressed = new int32 [_tigLen + 1];
  badStretched  = new int32 [_tigLen + 1];
  badNormal     = new int32 [_tigLen + 1];
  badAnti       = new int32 [_tigLen + 1];
  badOuttie     = new int32 [_tigLen + 1];

  memset(goodGraph,   0, sizeof(int32) * (_tigLen + 1));
  memset(badFwdGraph, 0, sizeof(int32) * (_tigLen + 1));
  memset(badRevGraph, 0, sizeof(int32) * (_tigLen + 1));

  memset(badExternalFwd, 0, sizeof(int32) * (_tigLen + 1));
  memset(badExternalRev, 0, sizeof(int32) * (_tigLen + 1));

  memset(badCompressed, 0, sizeof(int32) * (_tigLen + 1));
  memset(badStretched,  0, sizeof(int32) * (_tigLen + 1));
  memset(badNormal,     0, sizeof(int32) * (_tigLen + 1));
  memset(badAnti,       0, sizeof(int32) * (_tigLen + 1));
  memset(badOuttie,     0, sizeof(int32) * (_tigLen + 1));

  buildTable(utg);
  buildHappinessGraphs(utg);
}



MateLocation::~MateLocation() {
  delete [] goodGraph;
  delete [] badFwdGraph;
  delete [] badRevGraph;

  delete [] badExternalFwd;
  delete [] badExternalRev;

  delete [] badCompressed;
  delete [] badStretched;
  delete [] badNormal;
  delete [] badAnti;
  delete [] badOuttie;
}
            


void
MateLocation::buildTable(Unitig *utg) {

#if 0
  fprintf(logFile, "buildTable()-- unitig %d\n", utg->id());
#endif

  for (uint32 fi=0; fi<utg->dovetail_path_ptr->size(); fi++) {
    DoveTailNode  *frag = &(*utg->dovetail_path_ptr)[fi];

    if (FI->mateIID(frag->ident) == 0)
      //  Not mated.
      continue;

    uint32  mid = FI->mateIID(frag->ident);

    if (_iidToTableEntry.find(mid) == _iidToTableEntry.end()) {
      //  We didn't find the mate in the _table, so we know that we haven't seen
      //  either pair, and can create a new entry.
      //
      MateLocationEntry  mle;

      mle.mleFrgID1 = frag->ident;
      mle.mlePos1   = frag->position;
      mle.mleUtgID1 = utg->id();

      mle.mleFrgID2 = 0;
      mle.mlePos2   = NULL_SEQ_LOC;
      mle.mleUtgID2 = 0;

      mle.isGrumpy  = false;

      _iidToTableEntry[frag->ident] = _table.size();
      _table.push_back(mle);

#if 0
      fprintf(logFile, "buildTable()-- unitig %d frag %d at %d,%d\n",
              mle.mleUtgID1, mle.mleFrgID1, mle.mlePos1.bgn, mle.mlePos1.end);
#endif

    } else {
      //  Found the mate in the table.  Use that entry.
      //
      uint32  tid = _iidToTableEntry[mid];

      _table[tid].mleFrgID2 = frag->ident;
      _table[tid].mlePos2   = frag->position;
      _table[tid].mleUtgID2 = utg->id();

      _iidToTableEntry[frag->ident] = tid;

#if 0
      fprintf(logFile, "buildTable()-- unitig %d frag %d at %d,%d AND unitig %d frag %d at %d,%d\n",
              _table[tid].mleUtgID1, _table[tid].mleFrgID1, _table[tid].mlePos1.bgn, _table[tid].mlePos1.end,
              _table[tid].mleUtgID2, _table[tid].mleFrgID2, _table[tid].mlePos2.bgn, _table[tid].mlePos2.end);
#endif
    }
  }

  std::sort(_table.begin(), _table.end());

  for (uint32 i=0; i<_table.size(); i++) {
    _iidToTableEntry[_table[i].mleFrgID1] = i;
    _iidToTableEntry[_table[i].mleFrgID2] = i;
  }
}



void
MateLocation::buildHappinessGraphs(Unitig *utg) {

  for (uint32 mleidx=0; mleidx<_table.size(); mleidx++) {
    MateLocationEntry &loc = _table[mleidx];

    uint32 lib =  FI->libraryIID(loc.mleFrgID1);

    if (IS->valid(lib) == false)
      // Don't check libs that we didn't generate good stats for
      continue;

    int32 badMaxInter = static_cast<int32>(IS->mean(lib) + BADMATE_INTER_STDDEV * IS->stddev(lib));
    int32 badMinInter = static_cast<int32>(IS->mean(lib) - BADMATE_INTER_STDDEV * IS->stddev(lib));

    int32 badMaxIntra = static_cast<int32>(IS->mean(lib) + BADMATE_INTRA_STDDEV * IS->stddev(lib));
    int32 badMinIntra = static_cast<int32>(IS->mean(lib) - BADMATE_INTRA_STDDEV * IS->stddev(lib));

    //  To keep the results the same as the previous version (1.89)
    badMaxIntra = badMaxInter;
    badMinIntra = badMinInter;

    int32 dist = 0;
    int32 bgn  = 0;
    int32 end  = 0;

    //  Bgn and End MUST be signed.

    int32 matBgn = loc.mlePos2.bgn;
    int32 matEnd = loc.mlePos2.end;
    int32 matLen = (matBgn < matEnd) ? (matEnd - matBgn) : (matBgn - matEnd);

    int32 frgBgn = loc.mlePos1.bgn;
    int32 frgEnd = loc.mlePos1.end;
    int32 frgLen = (frgBgn < frgEnd) ? (frgEnd - frgBgn) : (frgBgn - frgEnd);

    if ((matLen >= MIN(badMaxInter, badMaxIntra)) ||
        (frgLen >= MIN(badMaxInter, badMaxIntra)))
      //  Yikes, fragment longer than insert size!
      continue;


    //  Until reset, assume this is a bad mate pair.
    loc.isGrumpy = true;


    //  If the mate is in another unitig, mark bad only if there is enough space to fit the mate in
    //  this unitig.
    if (loc.mleUtgID1 != loc.mleUtgID2) {
      if ((isReverse(loc.mlePos1) == true)  && (badMaxInter < frgBgn)) {
        incrRange(badExternalRev, -1, frgBgn - badMaxInter, frgEnd);
        incrRange(badExternalRev, -1, frgBgn, frgEnd);
        if (logFileFlagSet(LOG_HAPPINESS))
          fprintf(logFile, "buildHappinessGraph()--  pos %d,%d (len %d) and pos %d,%d (len %d) in tig of length %u -- bad external reverse\n",
                  frgBgn, frgEnd, frgLen,
                  matBgn, matEnd, matLen,
                  _tigLen);
        goto markBad;
      }

      if ((isReverse(loc.mlePos1) == false) && (badMaxInter < _tigLen - frgBgn)) {
        incrRange(badExternalFwd, -1, frgEnd, frgBgn + badMaxInter);
        incrRange(badExternalFwd, -1, frgEnd, frgBgn);
        if (logFileFlagSet(LOG_HAPPINESS))
          fprintf(logFile, "buildHappinessGraph()--  pos %d,%d (len %d) and pos %d,%d (len %d) in tig of length %u -- bad external forward\n",
                  frgBgn, frgEnd, frgLen,
                  matBgn, matEnd, matLen,
                  _tigLen);
        goto markBad;
      }

      //  Not enough space.  Not a grumpy mate pair.
      loc.isGrumpy = false;
      if (logFileFlagSet(LOG_HAPPINESS))
        fprintf(logFile, "buildHappinessGraph()--  pos %d,%d (len %d) and pos %d,%d (len %d) in tig of length %u -- not bad, not enough space\n",
                frgBgn, frgEnd, frgLen,
                matBgn, matEnd, matLen,
                _tigLen);
      continue;
    }


    //  Both mates are in this unitig.


    //  Same orientation?
    if ((isReverse(loc.mlePos1) == false) &&
        (isReverse(loc.mlePos2) == false)) {
      incrRange(badNormal, -1, MIN(frgBgn, matBgn), MAX(frgEnd, matEnd));
      if (logFileFlagSet(LOG_HAPPINESS))
        fprintf(logFile, "buildHappinessGraph()--  pos %d,%d (len %d) and pos %d,%d (len %d) in tig of length %u -- bad normal\n",
                frgBgn, frgEnd, frgLen,
                matBgn, matEnd, matLen,
                _tigLen);
      goto markBad;
    }

    if ((isReverse(loc.mlePos1) == true) &&
        (isReverse(loc.mlePos2) == true)) {
      incrRange(badAnti, -1, MIN(frgEnd, matEnd), MAX(frgBgn, matBgn));
      if (logFileFlagSet(LOG_HAPPINESS))
        fprintf(logFile, "buildHappinessGraph()--  pos %d,%d (len %d) and pos %d,%d (len %d) in tig of length %u -- bad anti\n",
                frgBgn, frgEnd, frgLen,
                matBgn, matEnd, matLen,
                _tigLen);
      goto markBad;
    }


    //  Check a special case for a circular unitig, outtie mates, but close enough to the end to
    //  plausibly be linking the ends together.
    //
    //   <---              --->
    //  ========unitig==========
    //
    if ((isReverse(loc.mlePos1) == true) &&
        (badMinIntra               <= frgBgn + _tigLen - matBgn) &&
        (frgBgn + _tigLen - matBgn <= badMaxIntra)) {
      loc.isGrumpy = false;  //  IT'S GOOD, kind of.
      if (logFileFlagSet(LOG_HAPPINESS))
        fprintf(logFile, "buildHappinessGraph()--  pos %d,%d (len %d) and pos %d,%d (len %d) in tig of length %u -- good because circular\n",
                frgBgn, frgEnd, frgLen,
                matBgn, matEnd, matLen,
                _tigLen);
      continue;
    }


    //  Outties?  True if pos1.end < pos2.bgn.  (For the second case, swap pos1 and pos2)
    //
    //  (pos1.end) <------   (pos1.bgn)
    //  (pos2.bgn)  -------> (pos2.end)
    //
    if ((isReverse(loc.mlePos1) == true)  && (loc.mlePos1.end < loc.mlePos2.bgn)) {
      incrRange(badOuttie, -1, MIN(frgBgn, frgEnd), MAX(matBgn, matEnd));
      if (logFileFlagSet(LOG_HAPPINESS))
        fprintf(logFile, "buildHappinessGraph()--  pos %d,%d (len %d) and pos %d,%d (len %d) in tig of length %u -- bad outtie (case 1)\n",
                frgBgn, frgEnd, frgLen,
                matBgn, matEnd, matLen,
                _tigLen);
      goto markBad;
    }
    if ((isReverse(loc.mlePos1) == false) && (loc.mlePos2.end < loc.mlePos1.bgn)) {
      incrRange(badOuttie, -1, MIN(frgBgn, frgEnd), MAX(matBgn, matEnd));
      if (logFileFlagSet(LOG_HAPPINESS))
        fprintf(logFile, "buildHappinessGraph()--  pos %d,%d (len %d) and pos %d,%d (len %d) in tig of length %u -- bad outtie (case 2)\n",
                frgBgn, frgEnd, frgLen,
                matBgn, matEnd, matLen,
                _tigLen);
      goto markBad;
    }

    //  So, now not NORMAL or ANTI or OUTTIE.  We must be left with innies.

    if (isReverse(loc.mlePos1) == false)
      //  First fragment is on the left, second is on the right.
      dist = loc.mlePos2.bgn - loc.mlePos1.bgn;
    else
      //  First fragment is on the right, second is on the left.
      dist = loc.mlePos1.bgn - loc.mlePos2.bgn;

    assert(dist >= 0);

    if (dist < badMinIntra) {
      incrRange(badCompressed, -1, MIN(frgBgn, matBgn), MAX(frgBgn, matBgn));
      if (logFileFlagSet(LOG_HAPPINESS))
        fprintf(logFile, "buildHappinessGraph()--  pos %d,%d (len %d) and pos %d,%d (len %d) in tig of length %u -- bad compressed\n",
                frgBgn, frgEnd, frgLen,
                matBgn, matEnd, matLen,
                _tigLen);
      goto markBad;
    }

    if (badMaxIntra < dist) {
      incrRange(badStretched, -1, MIN(frgBgn, matBgn), MAX(frgBgn, matBgn));
      if (logFileFlagSet(LOG_HAPPINESS))
        fprintf(logFile, "buildHappinessGraph()--  pos %d,%d (len %d) and pos %d,%d (len %d) in tig of length %u -- bad stretched\n",
                frgBgn, frgEnd, frgLen,
                matBgn, matEnd, matLen,
                _tigLen);
      goto markBad;
    }

    assert(badMinIntra <= dist);
    assert(dist        <= badMaxIntra);

    incrRange(goodGraph, 1, MIN(frgBgn, matBgn), MAX(frgBgn, matBgn));
    loc.isGrumpy = false;  //  IT'S GOOD!
    if (logFileFlagSet(LOG_HAPPINESS))
      fprintf(logFile, "buildHappinessGraph()--  pos %d,%d (len %d) and pos %d,%d (len %d) in tig of length %u -- GOOD!\n",
              frgBgn, frgEnd, frgLen,
              matBgn, matEnd, matLen,
              _tigLen);
    continue;

  markBad:

    //  Mark bad from the 3' end of the fragment till the upper limit where the mate should go.

    if (loc.mleUtgID1 == utg->id()) {
      assert(loc.mleFrgID1 != 0);
      if (isReverse(loc.mlePos1) == false) {
        //  Mark bad for forward fagment 1
        assert(frgBgn < frgEnd);
        bgn = frgEnd;
        end = frgBgn + badMaxIntra;
        incrRange(badFwdGraph, -1, bgn, end);
      } else {
        //  Mark bad for reverse fragment 1
        assert(frgEnd < frgBgn);
        bgn = frgBgn - badMaxIntra;
        end = frgEnd;
        incrRange(badRevGraph, -1, bgn, end);
      }
    }

    if (loc.mleUtgID2 == utg->id()) {
      assert(loc.mleFrgID2 != 0);
      if (isReverse(loc.mlePos2) == false) {
        //  Mark bad for forward fragment 2
        assert(matBgn < matEnd);
        bgn = matEnd;
        end = matBgn + badMaxIntra;
        incrRange(badFwdGraph, -1, bgn, end);
      } else {
        //  Mark bad for reverse fragment 2
        assert(matEnd < matBgn);
        bgn = matBgn - badMaxIntra;
        end = matEnd;
        incrRange(badRevGraph, -1, bgn, end);
      }
    }
  }  //  Over all MateLocationEntries in the table
}  //  buildHappinessGraph()
