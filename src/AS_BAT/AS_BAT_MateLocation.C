
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

static const char *rcsid = "$Id: AS_BAT_MateLocation.C,v 1.5 2012-07-30 01:21:01 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"

#include "AS_BAT_MateLocation.H"


MateLocation::MateLocation(UnitigVector &unitigs, Unitig *utg) {
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

  _tig    = utg;
  _tigLen = utg->getLength();

  _numMates = 0;

  good   = new int32 [_tigLen + 1];
  badFwd = new int32 [_tigLen + 1];
  badRev = new int32 [_tigLen + 1];

  badExternalFwd = new int32 [_tigLen + 1];
  badExternalRev = new int32 [_tigLen + 1];

  badCompressed = new int32 [_tigLen + 1];
  badStretched  = new int32 [_tigLen + 1];
  badNormal     = new int32 [_tigLen + 1];
  badAnti       = new int32 [_tigLen + 1];
  badOuttie     = new int32 [_tigLen + 1];

  memset(good,   0, sizeof(int32) * (_tigLen + 1));
  memset(badFwd, 0, sizeof(int32) * (_tigLen + 1));
  memset(badRev, 0, sizeof(int32) * (_tigLen + 1));

  memset(badExternalFwd, 0, sizeof(int32) * (_tigLen + 1));
  memset(badExternalRev, 0, sizeof(int32) * (_tigLen + 1));

  memset(badCompressed, 0, sizeof(int32) * (_tigLen + 1));
  memset(badStretched,  0, sizeof(int32) * (_tigLen + 1));
  memset(badNormal,     0, sizeof(int32) * (_tigLen + 1));
  memset(badAnti,       0, sizeof(int32) * (_tigLen + 1));
  memset(badOuttie,     0, sizeof(int32) * (_tigLen + 1));

  nunmated = 0;

  ngood[0]            = ngood[1]            = ngood[2]            = 0;
  nbadFwd[0]          = nbadFwd[1]          = nbadFwd[2]          = 0;
  nbadRev[0]          = nbadRev[1]          = nbadRev[2]          = 0;

  ngoodExternalFwd[0] = ngoodExternalFwd[1] = ngoodExternalFwd[2] = 0;
  ngoodExternalRev[0] = ngoodExternalRev[1] = ngoodExternalRev[2] = 0;

  nbadExternalFwd[0]  = nbadExternalFwd[1]  = nbadExternalFwd[2]  = 0;
  nbadExternalRev[0]  = nbadExternalRev[1]  = nbadExternalRev[2]  = 0;

  nbadCompressed[0]   = nbadCompressed[1]   = nbadCompressed[2]   = 0;
  nbadStretched[0]    = nbadStretched[1]    = nbadStretched[2]    = 0;
  nbadNormal[0]       = nbadNormal[1]       = nbadNormal[2]       = 0;
  nbadAnti[0]         = nbadAnti[1]         = nbadAnti[2]         = 0;
  nbadOuttie[0]       = nbadOuttie[1]       = nbadOuttie[2]       = 0;

  buildTable();
  buildHappinessGraphs(unitigs);
}



MateLocation::~MateLocation() {
  delete [] good;
  delete [] badFwd;
  delete [] badRev;

  delete [] badExternalFwd;
  delete [] badExternalRev;

  delete [] badCompressed;
  delete [] badStretched;
  delete [] badNormal;
  delete [] badAnti;
  delete [] badOuttie;
}
            



void
MateLocation::dumpHappiness(const char *prefix, const char *name) {
  char  dirname[FILENAME_MAX] = {0};
  char  outname[FILENAME_MAX] = {0};

  sprintf(dirname, "%s.%03u.%s.mateHappiness",
          prefix, logFileOrder, name);
  sprintf(outname, "%s.%03u.%s.mateHappiness/utg%09u.mateHappiness",
          prefix, logFileOrder, name, _tig->id());

  if (AS_UTL_fileExists(dirname, TRUE, TRUE) == 0)
    AS_UTL_mkdir(dirname);

  FILE *F = fopen(outname, "w");

  for (int32 i=0; i<_tigLen; i++)
    fprintf(F, "%u\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
            i,
            good[i],
            badFwd[i],
            badRev[i],
            badExternalFwd[i],
            badExternalRev[i],
            badCompressed[i],
            badStretched[i],
            badNormal[i],
            badAnti[i],
            badOuttie[i]);

  fclose(F);
}




void
MateLocation::buildTable(void) {

#if 0
  writeLog("buildTable()-- unitig %d\n", _tig->id());
#endif

  for (uint32 fi=0; fi<_tig->ufpath.size(); fi++) {
    ufNode  *frag = &_tig->ufpath[fi];

    if (FI->mateIID(frag->ident) == 0) {
      //  Not mated.
      nunmated++;
      continue;
    }

    uint32  mid = FI->mateIID(frag->ident);

    if (_iidToTableEntry.find(mid) == _iidToTableEntry.end()) {
      //  We didn't find the mate in the _table, so we know that we haven't seen
      //  either pair, and can create a new entry.
      //
      MateLocationEntry  mle;

      mle.mleFrgID1 = frag->ident;
      mle.mlePos1   = frag->position;
      mle.mleUtgID1 = _tig->id();

      mle.mleFrgID2 = 0;
      mle.mlePos2   = NULL_SEQ_LOC;
      mle.mleUtgID2 = 0;

      mle.isGrumpy  = false;

      assert(_table.size() > 0);

      _iidToTableEntry[frag->ident] = _table.size();
      _table.push_back(mle);

#if 0
      writeLog("buildTable()-- unitig %d frag %d at %d,%d\n",
              mle.mleUtgID1, mle.mleFrgID1, mle.mlePos1.bgn, mle.mlePos1.end);
#endif

    } else {
      //  Found the mate in the table.  Use that entry.
      //
      uint32  tid = _iidToTableEntry[mid];

      assert(tid != 0);

      _table[tid].mleFrgID2 = frag->ident;
      _table[tid].mlePos2   = frag->position;
      _table[tid].mleUtgID2 = _tig->id();

      _iidToTableEntry[frag->ident] = tid;

      _numMates++;

#if 0
      writeLog("buildTable()-- unitig %d frag %d at %d,%d AND unitig %d frag %d at %d,%d\n",
              _table[tid].mleUtgID1, _table[tid].mleFrgID1, _table[tid].mlePos1.bgn, _table[tid].mlePos1.end,
              _table[tid].mleUtgID2, _table[tid].mleFrgID2, _table[tid].mlePos2.bgn, _table[tid].mlePos2.end);
#endif
    }
  }

  std::sort(_table.begin(), _table.end());

  assert(_table[0].mleFrgID1 == 0);
  assert(_table[0].mleFrgID2 == 0);

  assert(_table[0].mleUtgID1 == 0);
  assert(_table[0].mleUtgID2 == 0);

  for (uint32 i=0; i<_table.size(); i++) {
    _iidToTableEntry[_table[i].mleFrgID1] = i;
    _iidToTableEntry[_table[i].mleFrgID2] = i;
  }
}



void
MateLocation::buildHappinessGraphs(UnitigVector &unitigs) {

  //  First entry is always zero.  Needed for the getById() accessor.

  assert(_table[0].mleFrgID1 == 0);
  assert(_table[0].mleFrgID2 == 0);

  assert(IS != NULL);

  for (uint32 mleidx=1; mleidx<_table.size(); mleidx++) {
    MateLocationEntry &loc = _table[mleidx];

    //  We MUST have mleFrgID1 defined.  If mleFrgID2 is not defined, then the mate is external.
    assert(loc.mleFrgID1 != 0);

    //  Well, this bit of ugly is here to fill out the location of the mate (when it is in a
    //  different unitig, the location for this fragment is not set)....EXCEPT that the mate might
    //  not even be placed in a unitig yet (if it is a contain, and we're called before contains are
    //  placed).
    //
    //  This is currently only used for logging -- AND statistics on mate happiness.  Later we
    //  should use it to determine if the mate is buried in the other unitig, which would make the
    //  fragment in this unitig bad.
    //
    if (loc.mleFrgID2 == 0) {
      assert(loc.mleUtgID2 == 0);
      loc.mleFrgID2 = FI->mateIID(loc.mleFrgID1);
      loc.mleUtgID2 = _tig->fragIn(loc.mleFrgID2);

      if (loc.mleUtgID2 != 0) {
        Unitig  *mt = unitigs[loc.mleUtgID2];
        uint32   fi = _tig->pathPosition(loc.mleFrgID2);
        loc.mlePos2 = mt->ufpath[fi].position;
      }
    }

    uint32 lib =  FI->libraryIID(loc.mleFrgID1);

    if (lib == 0)
      //  Shouldn't occur, but just in case, ignore fragments in the legacy library.
      continue;

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

    int32 nContained  = 0;

    if (OG->getBestContainer(loc.mleFrgID1)->isContained == true)
      nContained++;

    if (OG->getBestContainer(loc.mleFrgID2)->isContained == true)
      nContained++;

    if ((matLen >= MIN(badMaxInter, badMaxIntra)) ||
        (frgLen >= MIN(badMaxInter, badMaxIntra)))
      //  Yikes, fragment longer than insert size!
      continue;


    //  Until reset, assume this is a bad mate pair.
    loc.isGrumpy = true;

    int32 ULEN1 = 0;  //unitigs[loc.mleUtgID1]->getLength()
    int32 ULEN2 = 0;  //unitigs[loc.mleUtgID2]->getLength()


    //  If the mate is in another unitig, mark bad only if there is enough space to fit the mate in
    //  this unitig.
    //
    //  TODO: OR if there isn't enough space on the end of the OTHER unitig
    //
    if (loc.mleUtgID1 != loc.mleUtgID2) {
      if ((isReverse(loc.mlePos1) == true)  && (badMaxInter < frgBgn)) {
        incrRange(badExternalRev, -1, frgBgn - badMaxInter, frgEnd);
        incrRange(badExternalRev, -1, frgBgn, frgEnd);
        nbadExternalRev[nContained]++;
        if (logFileFlagSet(LOG_HAPPINESS))
          writeLog("buildHappinessGraph()--  unitig %d (len %d) frag %d pos %d,%d (len %d) and unitig %d (len %d) frag %d pos %d,%d (len %d) -- bad external reverse\n",
                  loc.mleUtgID1, ULEN1, loc.mleFrgID1, frgBgn, frgEnd, frgLen,
                  loc.mleUtgID2, ULEN2, loc.mleFrgID2, matBgn, matEnd, matLen);
        goto markBad;
      }

      if ((isReverse(loc.mlePos1) == false) && (badMaxInter < _tigLen - frgBgn)) {
        incrRange(badExternalFwd, -1, frgEnd, frgBgn + badMaxInter);
        incrRange(badExternalFwd, -1, frgEnd, frgBgn);
        nbadExternalFwd[nContained]++;
        if (logFileFlagSet(LOG_HAPPINESS))
          writeLog("buildHappinessGraph()--  unitig %d (len %d) frag %d pos %d,%d (len %d) and unitig %d (len %d) frag %d pos %d,%d (len %d) -- bad external forward\n",
                  loc.mleUtgID1, ULEN1, loc.mleFrgID1, frgBgn, frgEnd, frgLen,
                  loc.mleUtgID2, ULEN2, loc.mleFrgID2, matBgn, matEnd, matLen);
        goto markBad;
      }

      //  Not enough space.  Not a grumpy mate pair.
      loc.isGrumpy = false;

      if (isReverse(loc.mlePos1) == true)
        ngoodExternalRev[nContained]++;
      else
        ngoodExternalFwd[nContained]++;

      if (logFileFlagSet(LOG_HAPPINESS))
        writeLog("buildHappinessGraph()--  unitig %d (len %d) frag %d pos %d,%d (len %d) and unitig %d (len %d) frag %d pos %d,%d (len %d) -- not bad, not enough space\n",
                loc.mleUtgID1, ULEN1, loc.mleFrgID1, frgBgn, frgEnd, frgLen,
                loc.mleUtgID2, ULEN2, loc.mleFrgID2, matBgn, matEnd, matLen);
      continue;
    }


    //  Both mates are in this unitig.


    //  Same orientation?
    if ((isReverse(loc.mlePos1) == false) &&
        (isReverse(loc.mlePos2) == false)) {
      incrRange(badNormal, -1, MIN(frgBgn, matBgn), MAX(frgEnd, matEnd));
      nbadNormal[nContained]++;
      if (logFileFlagSet(LOG_HAPPINESS))
        writeLog("buildHappinessGraph()--  unitig %d (len %d) frag %d pos %d,%d (len %d) and unitig %d (len %d) frag %d pos %d,%d (len %d) -- bad normal\n",
                loc.mleUtgID1, ULEN1, loc.mleFrgID1, frgBgn, frgEnd, frgLen,
                loc.mleUtgID2, ULEN2, loc.mleFrgID2, matBgn, matEnd, matLen);
      goto markBad;
    }

    if ((isReverse(loc.mlePos1) == true) &&
        (isReverse(loc.mlePos2) == true)) {
      incrRange(badAnti, -1, MIN(frgEnd, matEnd), MAX(frgBgn, matBgn));
      nbadAnti[nContained]++;
      if (logFileFlagSet(LOG_HAPPINESS))
        writeLog("buildHappinessGraph()--  unitig %d (len %d) frag %d pos %d,%d (len %d) and unitig %d (len %d) frag %d pos %d,%d (len %d) -- bad anti\n",
                loc.mleUtgID1, ULEN1, loc.mleFrgID1, frgBgn, frgEnd, frgLen,
                loc.mleUtgID2, ULEN2, loc.mleFrgID2, matBgn, matEnd, matLen);
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
      ngood[nContained]++;
      if (logFileFlagSet(LOG_HAPPINESS))
        writeLog("buildHappinessGraph()--  unitig %d (len %d) frag %d pos %d,%d (len %d) and unitig %d (len %d) frag %d pos %d,%d (len %d) -- good because circular\n",
                loc.mleUtgID1, ULEN1, loc.mleFrgID1, frgBgn, frgEnd, frgLen,
                loc.mleUtgID2, ULEN2, loc.mleFrgID2, matBgn, matEnd, matLen);
      continue;
    }


    //  Outties?  True if pos1.end < pos2.bgn.  (For the second case, swap pos1 and pos2)
    //
    //  (pos1.end) <------   (pos1.bgn)
    //  (pos2.bgn)  -------> (pos2.end)
    //
    if ((isReverse(loc.mlePos1) == true)  && (loc.mlePos1.end < loc.mlePos2.bgn)) {
      incrRange(badOuttie, -1, MIN(frgBgn, frgEnd), MAX(matBgn, matEnd));
      nbadOuttie[nContained]++;
      if (logFileFlagSet(LOG_HAPPINESS))
        writeLog("buildHappinessGraph()--  unitig %d (len %d) frag %d pos %d,%d (len %d) and unitig %d (len %d) frag %d pos %d,%d (len %d) -- bad outtie (case 1)\n",
                loc.mleUtgID1, ULEN1, loc.mleFrgID1, frgBgn, frgEnd, frgLen,
                loc.mleUtgID2, ULEN2, loc.mleFrgID2, matBgn, matEnd, matLen);
      goto markBad;
    }
    if ((isReverse(loc.mlePos1) == false) && (loc.mlePos2.end < loc.mlePos1.bgn)) {
      incrRange(badOuttie, -1, MIN(frgBgn, frgEnd), MAX(matBgn, matEnd));
      nbadOuttie[nContained]++;
      if (logFileFlagSet(LOG_HAPPINESS))
        writeLog("buildHappinessGraph()--  unitig %d (len %d) frag %d pos %d,%d (len %d) and unitig %d (len %d) frag %d pos %d,%d (len %d) -- bad outtie (case 2)\n",
                loc.mleUtgID1, ULEN1, loc.mleFrgID1, frgBgn, frgEnd, frgLen,
                loc.mleUtgID2, ULEN2, loc.mleFrgID2, matBgn, matEnd, matLen);
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
      nbadCompressed[nContained]++;
      if (logFileFlagSet(LOG_HAPPINESS))
        writeLog("buildHappinessGraph()--  unitig %d (len %d) frag %d pos %d,%d (len %d) and unitig %d (len %d) frag %d pos %d,%d (len %d) -- bad compressed\n",
                loc.mleUtgID1, ULEN1, loc.mleFrgID1, frgBgn, frgEnd, frgLen,
                loc.mleUtgID2, ULEN2, loc.mleFrgID2, matBgn, matEnd, matLen);
      goto markBad;
    }

    if (badMaxIntra < dist) {
      incrRange(badStretched, -1, MIN(frgBgn, matBgn), MAX(frgBgn, matBgn));
      nbadStretched[nContained]++;
      if (logFileFlagSet(LOG_HAPPINESS))
        writeLog("buildHappinessGraph()--  unitig %d (len %d) frag %d pos %d,%d (len %d) and unitig %d (len %d) frag %d pos %d,%d (len %d) -- bad stretched\n",
                loc.mleUtgID1, ULEN1, loc.mleFrgID1, frgBgn, frgEnd, frgLen,
                loc.mleUtgID2, ULEN2, loc.mleFrgID2, matBgn, matEnd, matLen);
      goto markBad;
    }

    assert(badMinIntra <= dist);
    assert(dist        <= badMaxIntra);

    incrRange(good, 1, MIN(frgBgn, matBgn), MAX(frgBgn, matBgn));
    loc.isGrumpy = false;  //  IT'S GOOD!
    ngood[nContained]++;
    if (logFileFlagSet(LOG_HAPPINESS))
      writeLog("buildHappinessGraph()--  unitig %d (len %d) frag %d pos %d,%d (len %d) and unitig %d (len %d) frag %d pos %d,%d (len %d) -- GOOD!\n",
              loc.mleUtgID1, ULEN1, loc.mleFrgID1, frgBgn, frgEnd, frgLen,
              loc.mleUtgID2, ULEN2, loc.mleFrgID2, matBgn, matEnd, matLen);
    continue;

  markBad:

    //  Mark bad from the 3' end of the fragment till the upper limit where the mate should go.

    if (loc.mleUtgID1 == _tig->id()) {
      assert(loc.mleFrgID1 != 0);
      if (isReverse(loc.mlePos1) == false) {
        //  Mark bad for forward fagment 1
        assert(frgBgn < frgEnd);
        bgn = frgEnd;
        end = frgBgn + badMaxIntra;
        incrRange(badFwd, -1, bgn, end);
      } else {
        //  Mark bad for reverse fragment 1
        assert(frgEnd < frgBgn);
        bgn = frgBgn - badMaxIntra;
        end = frgEnd;
        incrRange(badRev, -1, bgn, end);
      }
    }

    if (loc.mleUtgID2 == _tig->id()) {
      assert(loc.mleFrgID2 != 0);
      if (isReverse(loc.mlePos2) == false) {
        //  Mark bad for forward fragment 2
        assert(matBgn < matEnd);
        bgn = matEnd;
        end = matBgn + badMaxIntra;
        incrRange(badFwd, -1, bgn, end);
      } else {
        //  Mark bad for reverse fragment 2
        assert(matEnd < matBgn);
        bgn = matBgn - badMaxIntra;
        end = matEnd;
        incrRange(badRev, -1, bgn, end);
      }
    }
  }  //  Over all MateLocationEntries in the table
}  //  buildHappinessGraph()
