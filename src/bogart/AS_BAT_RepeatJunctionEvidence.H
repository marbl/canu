

class repeatJunctionEvidence {
public:
  repeatJunctionEvidence(Unitig *target, overlapPlacement &op);

  bool operator<(repeatJunctionEvidence const that) const {
    return(point < that.point);
  };

  //  Placed position in the unitig.  The fragment has overlap evidence that covers
  //  'coveredbgn-coveredend' and the unaligned hangs of the fragment cover the remaining.
  //  Orientation of the fragment is not tracked in this position.

  int32               uncovered5bgn;
  int32               uncovered5end;

  int32               coveredbgn;
  int32               coveredend;

  int32               uncovered3bgn;
  int32               uncovered3end;

  bool                is3;

  int32               point;

  AS_IID              eviFrag;  //  Evidence fragment (orientation is meaningless)
  FragmentEnd         tigFrag;  //  Unitig fragment
  
  uint32              numOvl3p;  //  Weight/score of the incoming fragment
  uint32              numOvl5p;

  overlapPlacement    place;  //  Possibly unused
};


class repeatUniqueBreakPoint {
public:
  repeatUniqueBreakPoint() {
    point       = 0;
    breakFrag   = FragmentEnd();
    rptLeft     = false;
  };
  repeatUniqueBreakPoint(uint32 point_, FragmentEnd breakFrag_, bool rptLeft_) {
    point       = point_;
    breakFrag   = breakFrag_;
    rptLeft     = rptLeft_;
  };

  bool operator<(repeatUniqueBreakPoint const that) const {
    return(point < that.point);
  };

  uint32        point;       //  Position in the unitig where we want to break
  FragmentEnd   breakFrag;   //  Fragment end that position is
  bool          rptLeft;     //  Repeat is to the left of the point
};



class repeatRegion {
public:
  repeatRegion(uint32 bgn_, uint32 end_) {
    bgn             = bgn_;
    end             = end_;
    ejectUnanchored = false;
  };

  uint32         bgn;
  uint32         end;

  bool           ejectUnanchored;

  repeatUniqueBreakPoint  rujBgn;
  repeatUniqueBreakPoint  rujEnd;
};








repeatJunctionEvidence::repeatJunctionEvidence(Unitig *target, overlapPlacement &op) {
  FragmentEnd  *end3;
  FragmentEnd  *end5;

  if (op.position.bgn < op.position.end) {
    //  Fragment is placed forward.
    uncovered5bgn = op.position.bgn;
    uncovered5end = op.verified.bgn;

    coveredbgn    = op.verified.bgn;
    coveredend    = op.verified.end;

    uncovered3bgn = op.verified.end;
    uncovered3end = op.position.end;

    end5 = &op.frag5p;
    end3 = &op.frag3p;

  } else {
    //  Fragment is placed reverse.
    uncovered5bgn = op.position.end;
    uncovered5end = op.verified.end;

    coveredbgn    = op.verified.end;
    coveredend    = op.verified.bgn;

    uncovered3bgn = op.verified.bgn;
    uncovered3end = op.position.bgn;

    end5 = &op.frag3p;
    end3 = &op.frag5p;
  }

  if (uncovered5bgn < 0) {
    uncovered5bgn = 0;   //  Ignore ends that extend past the edge of the unitig
    uncovered5end = 0;
  }

  if (target->getLength() <= uncovered3end) {
    uncovered3bgn = 0;
    uncovered3end = 0;
  }

  //  Set flags indicating (true) if the unaligned portion on that end is significant enough to
  //  be trusted as a real unaligned piece.
  //
  //  If both are true, this looks like a short repeat in a long fragment.  The only way (I can
  //  think of anyway) this could occur is from an alignment to a short contained fragment.
  //
  //  If both are false, the fragment has aligned fully.

#define UNCOVERED_NOISE_FILTER 10

  bool  save5 = (uncovered5bgn + UNCOVERED_NOISE_FILTER < uncovered5end);
  bool  save3 = (uncovered3bgn + UNCOVERED_NOISE_FILTER < uncovered3end);

  eviFrag = op.frgID;

  numOvl3p = 0;
  numOvl5p = 0;

  if (save5 == save3) {
    tigFrag       = FragmentEnd();
    uncovered5bgn = 0;
    uncovered5end = 0;
    uncovered3bgn = 0;
    uncovered3end = 0;

    return;
  }

  place = op;

  if (save5) {
    tigFrag       = *end5;
    point         = coveredbgn;
    is3           = false;
    uncovered3bgn = 0;
    uncovered3end = 0;

  } else {
    tigFrag       = *end3;
    point         = coveredend;
    is3           = true;
    uncovered5bgn = 0;
    uncovered5end = 0;
  }

  //  When the fragment is placed on the unitig, it returns the first/last fragment aligning
  //  to each end, and the (dovetail) end it aligns to.  If the fragment aligns with a spur
  //  (as here) the end computation is backwards.
  //
  //        ------------       aligned fragment
  //    5--------3             unitig fragment -- placement returns 3' end
  //
  //        \-----------       aligned fragment with spur
  //         5--------------3  unitig fragment -- placement return 3' end also
  //
  //  We reverse these ends to show the end with the break point.
  //
  tigFrag = FragmentEnd(tigFrag.fragId(), !tigFrag.frag3p());

#if 0
  writeLog("markRepeats()-- tig frag %8d ovl frag %8d %6d-%6d %5.2f%% tig pos %8d/%c' (%4d) %6d-%6d (%4d) %8d/%c'\n",
          op[pl].refID, op[pl].frgID,
          op[pl].covered.bgn,
          op[pl].covered.end,
          op[pl].fCoverage * 100.0,
          end5->fragId(), (end5->frag5p() ? '5' : '3'),
          uncovered5end - uncovered5bgn,
          coveredbgn, coveredend,
          uncovered3end - uncovered3bgn,
          end3->fragId(), (end3->frag5p() ? '5' : '3'));
#endif
}
