
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

#ifndef INCLUDE_AS_BAT_PLACEFRAGUSINGOVERLAPS
#define INCLUDE_AS_BAT_PLACEFRAGUSINGOVERLAPS

static const char *rcsid_INCLUDE_AS_BAT_PLACEFRAGUSINGOVERLAPS = "$Id: AS_BAT_PlaceFragUsingOverlaps.H,v 1.1 2010-12-06 08:03:48 brianwalenz Exp $";

class overlapPlacement {
public:
  overlapPlacement() {
    frgID                 = 0;

    tigID                 = 0;
    position.bgn          = 0;
    position.end          = 0;

    nReverse              = 0;
    nForward              = 0;

    fCoverage             = 0.0;

    errors                = 0.0;
    aligned               = 0;
    covered.bgn           = 0;
    covered.end           = 0;
  };
  ~overlapPlacement() {
  };

  bool operator<(const overlapPlacement &that) const {
    if (tigID != that.tigID)
      return(tigID < that.tigID);
    return(position < that.position);
  };

public:
  uint32            frgID;      //  ID of the fragment this position is for.

  uint32            tigID;      //  Unitig ID of this placement
  SeqInterval       position;   //  Unitig position of this placement

  uint32            nReverse;   //  Num overlaps aligned reverse
  uint32            nForward;   //  Num overlaps aligned forward

  double            fCoverage;  //  Coverage of the fragment

  double            errors;     //  number of errors in alignments
  uint32            aligned;    //  number of bases in alignments
  SeqInterval       covered;    //  Position of the overlap on the read
};


bool
placeFragUsingOverlaps(UnitigVector &unitigs,
                       uint32 fid,
                       OverlapStore *ovlStoreUniq,
                       OverlapStore *ovlStoreRept,
                       vector<overlapPlacement> &placements);

#endif  //  INCLUDE_AS_BAT_PLACEFRAGUSINGOVERLAPS
