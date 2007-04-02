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

/* RCS info
 * $Id: AS_BOG_MateChecker.hh,v 1.5 2007-04-02 21:00:16 eliv Exp $
 * $Revision: 1.5 $
*/

#ifndef INCLUDE_AS_BOG_MATECHEKER
#define INCLUDE_AS_BOG_MATECHEKER

#include <map>

#include "AS_BOG_UnitigGraph.hh"

extern "C" {
#include "AS_PER_gkpStore.h"
}

namespace AS_BOG{
    typedef std::map<iuid,iuid> IdMap;
    typedef IdMap::const_iterator IdMapConstIter;
    typedef std::vector<int> DistanceList;
    typedef DistanceList::const_iterator DistanceListCIter;
    typedef std::map<iuid,DistanceList> LibraryDistances;
    typedef LibraryDistances::const_iterator LibDistsConstIter;

    struct MateInfo {
        iuid mate;
        iuid lib;
    };
    typedef std::map<iuid,MateInfo> MateMap;
    static const MateInfo NULL_MATE_INFO = {0,0};

    struct MateLocation {
        SeqInterval pos1;
        SeqInterval pos2;
    };
    /*
    inline bool operator==(SeqInterval a, SeqInterval b) {
        if (a.bgn == b.bgn && a.end == b.end ||
            a.bgn == b.end && a.end == b.bgn)
            return true;
        else
            return false;
    };
    inline bool operator<(SeqInterval a, SeqInterval b)
    {
        if ( isReverse(a) ) {
            if ( isReverse(b) ) return a.end < b.end;
            else                return a.end < b.bgn;
        } else {
            if ( isReverse(b) ) return a.bgn < b.end;
            else                return a.bgn < b.bgn; 
        }
    };
    inline bool operator==(MateLocation a, MateLocation b) {
        if (a.pos1 == b.pos1 && a.pos2 == b.pos2)
            return true;
        else
            return false;
    };
    inline bool operator<(MateLocation a, MateLocation b) {
        if (a.pos1 < b.pos1) return true;
        if (a.pos2 < b.pos2) return true;
        else                 return false;
    };
    */
    static const SeqInterval NULL_MATE_LOC = {0,0};
    typedef std::map<iuid,MateLocation> MateLocMap;

    struct DistanceCompute {
        double stddev;
        double mean;
        double sumSquares;
        double sumDists;
        int numPairs;
    };
    typedef std::map<iuid,DistanceCompute> LibraryStats;

    ///////////////////////////////////////////////////////////////////////////

    struct MateChecker{
        ~MateChecker();

        void readStore(const char *);   // reads the gkpStore mate info into memory
        // returns the mate iid for the given iid, or zero if none
        MateInfo getMateInfo(iuid);

        // Checks size of mates internal to unitig
        LibraryStats* checkUnitig( Unitig* );
        // Compute good and bad coverage graphs for a unitig
        void computeMateCoverage( Unitig*, LibraryStats &);
        // Computes stddev and mate coverage over all unitigs
        void checkUnitigGraph( UnitigGraph& );

        private:
            MateMap _mates;
            LibraryDistances _dists; // all distances 
    };
}
#endif
