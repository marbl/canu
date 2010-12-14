
/**************************************************************************
 * Copyright (C) 2010, J Craig Venter Institute. All rights reserved.
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

const char *mainid = "$Id: bogusness.C,v 1.4 2010-12-14 01:46:35 brianwalenz Exp $";

#include "AS_BAT_bogusUtil.H"

//  Reads snapper/nucmer output, figures out what the minimal set of alignments that cover each
//  aligned sequence, then compares those alignments against the ideal unitigs.
//
//  It reports:
//    bubbles that should have been popped.
//    unitigs that span a unique-repeat junction (one end in unique, one end in repeat).
//    unitigs that include a repeat.
//    unitigs that include a repeat, and are chimeric at the repeat.
//    unitigs that should be joined.
//
//

#define IDEAL_REPT        0  //  A repeat unitig
#define IDEAL_UNIQ        1  //  A unique unitig
#define IDEAL_REPTUNIQ    2  //  A portion of this repeat unitig might be separable
#define IDEAL_UNIQWEAK    3  //  A portion of this unique unitig might be unjoinable
#define IDEAL_MIXED       4  //  Used for classifying bubbles

char *types[4] = { "REPT", "UNIQ", "SEPR", "WEAK" };

#define STATUS_BEGINSin   0
#define STATUS_ENDSin     1
#define STATUS_CONTAINS   2
#define STATUS_CONTAINED  3

char *statuses[4] = { "BEGINSin", "ENDSin", "CONTAINS", "CONTAINED" };



class bogusResult {
public:
  bogusResult(const char  *_utgID,
              int32 _utgIID,
              int32 _alignNum, int32 _alignTotal,
              int32 _utgBgn, int32 _utgEnd,
              int32 _genBgn, int32 _genEnd,
              int32 _status,
              int32 _type,
              int32 _idlNum,
              int32 _idlBgn, int32 _idlEnd,
              int32 _alignLen,
              double _utgCov,
              double _idlCov) {

    if (_utgID)
      strcpy(utgID, _utgID);
    else
      utgID[0] = 0;

    utgIID     = _utgIID;

    alignNum   = _alignNum;
    alignTotal = _alignTotal;

    utgBgn     = _utgBgn;
    utgEnd     = _utgEnd;

    genBgn     = _genBgn;
    genEnd     = _genEnd;

    status     = _status;
    type       = _type;

    idlNum     = _idlNum;
    idlBgn     = _idlBgn;
    idlEnd     = _idlEnd;

    alignLen   = _alignLen;
    utgCov     = _utgCov;
    idlCov     = _idlCov;
  };
  ~bogusResult() {
  };

  bool operator<(const bogusResult &that) const {
    if (genBgn < that.genBgn)
      return(true);
    if (genBgn > that.genBgn)
      return(false);
    return(idlBgn < that.idlBgn);
  };

  char    utgID[32];
  int32   utgIID;

  int32   alignNum;
  int32   alignTotal;

  int32   utgBgn;
  int32   utgEnd;

  int32   genBgn;
  int32   genEnd;

  int32   status;
  int32   type;

  int32   idlNum;
  int32   idlBgn;
  int32   idlEnd;

  int32   alignLen;
  double  utgCov;
  double  idlCov;
};



class idealUnitig {
public:
  idealUnitig() {
    bgn  = 0;
    end  = 0;
    type = 0;
  };

  idealUnitig(int32 b, int32 e, char t) {
    bgn  = b;
    end  = e;
    type = t;
  };

  ~idealUnitig() {
  };

  int32  bgn;
  int32  end;
  char   type;
};


#if 0
class actualUnitig {
public:
  actualUnitig() {
    bgn  = 0;
    end  = 0;
    type = 0;
  };

  actualUnitig(int32 b, int32 e, char t) {
    bgn  = b;
    end  = e;
    type = t;
  };

  ~actualUnitig() {
  };

  int32  bgn;
  int32  end;
  int32  unitigID;
  int32  alignNum;
  int32  alignMax;
  char   type;
};
#endif


void
loadIdealUnitigs(char *idealName,
                 vector<idealUnitig>  &ideal) {

  errno = 0;
  FILE *F = fopen(idealName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", idealName, strerror(errno)), exit(1);

  char          L[1024];
  splitToWords  S;
  char          t = 0;

  fgets(L, 1024, F);

  while (!feof(F)) {
    chomp(L);

    S.split(L);

    if        (strcmp(S[2], "REPT") == 0)
      t = (S.numWords() == 4) ? IDEAL_REPT : IDEAL_UNIQWEAK;

    else if (strcmp(S[2], "UNIQ") == 0)
      t = (S.numWords() == 4) ? IDEAL_UNIQ : IDEAL_REPTUNIQ;

    else
      fprintf(stderr, "Unknown type in '%s'\n", L), exit(1);

    ideal.push_back(idealUnitig(S(0), S(1), t));

    fgets(L, 1024, F);
  }

  fclose(F);

  fprintf(stderr, "Loaded %lu ideal unitigs.\n", ideal.size());
}





bool
isUnitigContained(int32  frgIID,
                  idealUnitig  &ideal,
                  genomeAlignment  &cover,
                  int32  &alen, int32 &ulen, int32 &ilen,
                  double &ufrac,
                  double &ifrac) {

  if ((ideal.bgn    <= cover.genBgn) &&
      (cover.genEnd <= ideal.end)) {
    alen  = cover.genEnd - cover.genBgn;
    ulen  = cover.genEnd - cover.genBgn;
    ilen  = ideal.end - ideal.bgn;
    ufrac = 100.0 * alen / ulen;
    ifrac = 100.0 * alen / ilen;
    return(true);
  }

  return(false);
}


bool
isUnitigContaining(int32  frgIID,
                   idealUnitig  &ideal,
                   genomeAlignment  &cover,
                   int32  &alen, int32 &ulen, int32 &ilen,
                   double &ufrac,
                   double &ifrac) {

  if ((cover.genBgn <= ideal.bgn) &&
      (ideal.end    <= cover.genEnd)) {
    alen  = ideal.end - ideal.bgn;
    ulen  = cover.genEnd - cover.genBgn;
    ilen  = ideal.end - ideal.bgn;
    ufrac = 100.0 * alen / ulen;
    ifrac = 100.0 * alen / ilen;
    return(true);
  }

  return(false);
}


bool
isUnitigEnding(int32  frgIID,
               idealUnitig  &ideal,
               genomeAlignment  &cover,
               int32  &alen, int32 &ulen, int32 &ilen,
               double &ufrac,
               double &ifrac) {

  if ((cover.genBgn <= ideal.bgn) &&
      (ideal.bgn    <= cover.genEnd)) {
    alen  = cover.genEnd - ideal.bgn;
    ulen  = cover.genEnd - cover.genBgn;
    ilen  = ideal.end - ideal.bgn;
    ufrac = 100.0 * alen / ulen;
    ifrac = 100.0 * alen / ilen;
    return(true);
  }

  return(false);
}


bool
isUnitigBeginning(int32  frgIID,
                  idealUnitig  &ideal,
                  genomeAlignment  &cover,
                  int32  &alen, int32 &ulen, int32 &ilen,
                  double &ufrac,
                  double &ifrac) {

  if ((cover.genBgn <= ideal.end) &&
      (ideal.end    <= cover.genEnd)) {
    alen  = ideal.end - cover.genBgn;
    ulen  = cover.genEnd - cover.genBgn;
    ilen  = ideal.end - ideal.bgn;
    ufrac = 100.0 * alen / ulen;
    ifrac = 100.0 * alen / ilen;
    return(true);
  }

  return(false);
}





int
main(int argc, char **argv) {
  char                      *nucmerName   = 0L;
  char                      *snapperName  = 0L;
  char                      *idealName = 0L;

  vector<idealUnitig>        ideal;

  vector<genomeAlignment>    genome;
  map<string, int32>         IIDmap;       //  Maps an ID string to an IID.
  vector<string>             IIDname;      //  Maps an IID to an ID string.
  vector<uint32>             IIDcount;     //  Maps an IID to the number of alignments

  vector<bogusResult>        results;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-nucmer") == 0) {
      nucmerName = argv[++arg];

    } else if (strcmp(argv[arg], "-snapper") == 0) {
      snapperName = argv[++arg];

    } else if (strcmp(argv[arg], "-ideal") == 0) {
      idealName = argv[++arg];

    } else {
      err++;
    }
    arg++;
  }
  if ((nucmerName == 0L) && (snapperName == 0L))
    fprintf(stderr, "ERROR: No input matches supplied (either -nucmer or -snapper).\n"), err++;
  if (idealName == 0L)
    fprintf(stderr, "ERROR: No ideal unitigs supplied (-ideal).\n"), err++;
  if (err) {
    exit(1);
  }

  loadIdealUnitigs(idealName, ideal);

  if (nucmerName)
    loadNucmer(nucmerName, genome, IIDmap, IIDname, IIDcount, NULL);

  if (snapperName)
    loadSnapper(snapperName, genome, IIDmap, IIDname, IIDcount, NULL);


  sort(genome.begin(), genome.end(), byFragmentID);



  for (uint32 bgn=0, lim=genome.size(); bgn<lim; ) {
    uint32 frgIID = genome[bgn].frgIID;

    vector<genomeAlignment>  cover;

    //  Find the range of alignments for a single fragment.

    uint32 end = bgn + 1;
    while ((end < genome.size()) &&
           (genome[bgn].frgIID == genome[end].frgIID))
      end++;

    //  Basically, discard any alignment that is contained in some other alignment.  We're not
    //  trying to explain the genome using the unitigs, we're trying to explain (evaluate) the
    //  unitigs using the genome.

    for (uint32 i=bgn; i<end; i++) {
      genomeAlignment  &I = genome[i];

      if (I.isDeleted)
        //  [i] already deleted
        continue;

      for (uint32 j=bgn; j<end; j++) {
        genomeAlignment  &J = genome[j];

        if (i == j)
          //  We are not contained in ourself
          continue;

        if (J.isDeleted)
          //  [j] already deleted
          continue;

        if ((I.frgBgn <= J.frgBgn) &&
            (J.frgEnd <= I.frgEnd))
          //  J contained in I.
          J.isDeleted = true;
      }
    }

    //  Copy whatever is left over to a list of 'covering' alignments

    for (uint32 i=bgn; i<end; i++) {
      genomeAlignment  &I = genome[i];

      if (I.isDeleted)
        //  [i] contained
        continue;
    
      cover.push_back(I);
    }

    //  Attempt to classify

    for (uint32 c=0; c<cover.size(); c++) {
      int32  alen  = 0;
      int32  ulen  = 0;
      int32  ilen  = 0;
      double ufrac = 0.0;
      double ifrac = 0.0;

      //  Pass one, is our unitig (region) completely 100% explained by any single ideal unitig?
      //  What usually happens is that our unitig is 100% contained in a unique ideal unitig, but
      //  there is (less than a fragment) of overlap with the adjacent repeat unitig.
      //
      //  This was annotated as either CORRECT or EXPLAINED, but that hasn't been useful so far.
      //    (ideal.size() == 1) ? " (CORRECT)" : " (EXPLAINED)");

      bool explained = false;

      for (uint32 i=0; i<ideal.size(); i++) {
        if (isUnitigContained(frgIID, ideal[i], cover[c], alen, ulen, ilen, ufrac, ifrac)) {
          results.push_back(bogusResult(IIDname[frgIID].c_str(), frgIID,
                                        c+1, cover.size(),
                                        cover[c].frgBgn, cover[c].frgEnd,
                                        cover[c].genBgn, cover[c].genEnd,
                                        STATUS_CONTAINED,
                                        ideal[i].type,
                                        i,
                                        ideal[i].bgn, ideal[i].end,
                                        alen,
                                        ufrac, ifrac));
          explained = true;
          break;
        }
      }

      if (explained)
        continue;

      for (uint32 i=0; i<ideal.size(); i++) {
        if (isUnitigContained(frgIID, ideal[i], cover[c], alen, ulen, ilen, ufrac, ifrac)) {
          results.push_back(bogusResult(IIDname[frgIID].c_str(), frgIID,
                                        c+1, cover.size(),
                                        cover[c].frgBgn, cover[c].frgEnd,
                                        cover[c].genBgn, cover[c].genEnd,
                                        STATUS_CONTAINED,
                                        ideal[i].type,
                                        i,
                                        ideal[i].bgn, ideal[i].end,
                                        alen,
                                        ufrac, ifrac));
          continue;
        }

        if (isUnitigContaining(frgIID, ideal[i], cover[c], alen, ulen, ilen, ufrac, ifrac)) {
          results.push_back(bogusResult(IIDname[frgIID].c_str(), frgIID,
                                        c+1, cover.size(),
                                        cover[c].frgBgn, cover[c].frgEnd,
                                        cover[c].genBgn, cover[c].genEnd,
                                        STATUS_CONTAINS,
                                        ideal[i].type,
                                        i,
                                        ideal[i].bgn, ideal[i].end,
                                        alen,
                                        ufrac, ifrac));
          continue;
        }

        if (isUnitigBeginning(frgIID, ideal[i], cover[c], alen, ulen, ilen, ufrac, ifrac)) {
          results.push_back(bogusResult(IIDname[frgIID].c_str(), frgIID,
                                        c+1, cover.size(),
                                        cover[c].frgBgn, cover[c].frgEnd,
                                        cover[c].genBgn, cover[c].genEnd,
                                        STATUS_BEGINSin,
                                        ideal[i].type,
                                        i,
                                        ideal[i].bgn, ideal[i].end,
                                        alen,
                                        ufrac, ifrac));

        }

        if (isUnitigEnding(frgIID, ideal[i], cover[c], alen, ulen, ilen, ufrac, ifrac)) {
          results.push_back(bogusResult(IIDname[frgIID].c_str(), frgIID,
                                        c+1, cover.size(),
                                        cover[c].frgBgn, cover[c].frgEnd,
                                        cover[c].genBgn, cover[c].genEnd,
                                        STATUS_ENDSin,
                                        ideal[i].type,
                                        i,
                                        ideal[i].bgn, ideal[i].end,
                                        alen,
                                        ufrac, ifrac));
        }
      }
    }

    bgn = end;
  }
  
  sort(results.begin(), results.end());

  for (uint32 i=0; i<results.size(); i++) {
    bogusResult *bi = &results[i];

    fprintf(stdout, "| %s || %d of %d || %d-%d || %d-%d || %s || %s || %05d || %d-%d || %dbp || %.2f%% || %.2f%%\n",
            bi->utgID,
            bi->alignNum, bi->alignTotal,
            bi->utgBgn, bi->utgEnd,
            bi->genBgn, bi->genEnd,
            statuses[bi->status],
            types[bi->type],
            bi->idlNum,
            bi->idlBgn, bi->idlEnd,
            bi->alignLen,
            bi->utgCov, bi->idlCov);
  }

  ////////////////////////////////////////
  //
  //  Attempt to analyze the result.
  //
  ////////////////////////////////////////
  //
  //  The array variables count over various tolerances at the edges.
  //
  //                         --------     spans the REPT, but includes a bit of the UNIQ on
  //        For example:  uuuuu    uuuu   either end.  If the UNIQ is less than 250bp, we'd
  //                          rrrrrr      numURU{0000}=0  numURU{0250}=1  numURU{0500}=1

  int32  tolerances[4] = { 0, 250, 500, 1000 };

  uint32  numChimeraTotal = 0;
  uint32  numChimera[1024] = {0};      //  Chimeric unitigs by number of pieces

  uint32  numBubbleUNIQ = 0;           //  Number of bubbles in unique regions (contained in a UNIQ and another unitig)
  uint32  numBubbleREPT = 0;           //  Number of bubbles in repeat regions (contained in a REPT and another unitig)
  uint32  numBubbleOTHR = 0;           //  Number of bubbles that span a region

  uint32  numIncompleteUNIQ[4] = {0};  //  Number of unitigs that do not completely span a UNIQ region (**)
  uint32  numIncompleteREPT[4] = {0};  //  Number of unitigs that do not completely span a REPT region (**)

  uint32  numURU[4] = {0};             //  Number of unitigs that correctly span a R, and might have a bit in U on either end.
  uint32  numRUR[4] = {0};             //  Number of unitigs that correctly span a U, and might have a bit in R on either end.

  uint32  numIllegalBorder[4] = {0};   //  Number of unitigs that cross a UNIQ/REPT border, or multiple borders

  ////////////////////////////////////////////////////////////////////////////////
  //  Chimera

  int32    *chm = new int32 [IIDname.size()];

  memset(chm, 0, sizeof(int32) * IIDname.size());

  for (uint32 i=0; i<results.size(); i++) {
    bogusResult *bi = &results[i];

    assert(bi->alignTotal < 1024);

    if ((bi->alignNum == 1) && (bi->alignTotal > 1))
      chm[bi->utgIID]++;
  }

  for (uint32 i=0; i<IIDname.size(); i++) {
    if (chm[i] == 0)
      continue;

    fprintf(stderr, "%s ", IIDname[i].c_str());
    numChimeraTotal++;
  }
  fprintf(stderr, "\n");

  delete [] chm;

  ////////////////////////////////////////////////////////////////////////////////
  //  Bubbles
  //
  //  Bubbles are painful to detect in this data.  Distinguishing between a bubble in a repeat
  //  (which has multiple alignments) and a chimera (which also has multiple alignments) is non
  //  trivial.
  //
  //  We can't use the bogusResult data, which breaks alignments across ideal unitig boundaries - we
  //  have to use the raw genomeAlignment data.
  //
  //  So, we do n^2 comparisons on genomeAlignment.  We want to find, well see the comment below.
  //
  //  NOTE THAT WE COUNT EACH BUBBLE INSTANCE SEPARATELY.  We are NOT counting bubble unitigs, but
  //  bubble placements.  This is probably a bug.

  int32    *bgn = new int32 [IIDname.size()];
  int32    *end = new int32 [IIDname.size()];
  int32    *bub = new int32 [IIDname.size()];
  int32    *typ = new int32 [IIDname.size()];

  memset(bgn, 0, sizeof(int32) * IIDname.size());
  memset(end, 0, sizeof(int32) * IIDname.size());
  memset(bub, 0, sizeof(int32) * IIDname.size());
  memset(typ, 0, sizeof(int32) * IIDname.size());

  for (uint32 i=0; i<genome.size(); i++) {
    int32  id = genome[i].frgIID;

    if (end[id] == 0) {
      bgn[id] = genome[i].frgBgn;
      end[id] = genome[i].frgEnd;
    } else {
      bgn[id] = MIN(bgn[id], genome[i].frgBgn);
      end[id] = MAX(end[id], genome[i].frgEnd);
    }
  }

  for (uint32 i=0; i<genome.size(); i++) {
    genomeAlignment *gi = &genome[i];

    for (uint32 j=0; j<genome.size(); j++) {
      genomeAlignment *gj = &genome[j];

      //  j contained in i AND j the maximal alignment -- A BUBBLE!
      //
      if ((i != j) &&
          (gi->genBgn <= gj->genBgn) &&
          (gj->genEnd <= gi->genEnd) &&
          (gj->frgBgn == bgn[gj->frgIID]) &&
          (gj->frgEnd == end[gj->frgIID])) {
        bub[gj->frgIID]++;
        //typ[gj->frgIID] = IDEAL_MIXED;
      }
    }
  }

  //  Now, for any unitig marked as a bubble, try to determine the type.  The rule is
  //  simple.  Three cases:
  //    matches are all UNIQ
  //    matches are all REPT
  //    matches are of mixed type (bubble spans a junction)
  //  
  for (uint32 b=0; b<results.size(); b++) {
    int32  id = results[b].utgIID;

    if (bub[id] == 0)
      continue;

    if (typ[id] == 0)
      typ[id] = results[b].type;

    if (typ[id] != results[b].type)
      typ[id] = IDEAL_MIXED;
  }

  for (uint32 i=0; i<IIDname.size(); i++) {
    if (bub[i] == 0)
      continue;

    switch (typ[i]) {
      case IDEAL_UNIQ:  numBubbleUNIQ++;  break;
      case IDEAL_REPT:  numBubbleREPT++;  break;
      case IDEAL_MIXED: numBubbleOTHR++;  break;
      default:                            break;
    }
  }

  delete [] bgn;
  delete [] end;
  delete [] bub;
  delete [] typ;

  ////////////////////////////////////////////////////////////////////////////////
  //  












  //  Write a 'nice' report.

  fprintf(stderr, "numUnitigs:      %lu\n", IIDname.size());
  fprintf(stderr, "\n");
  fprintf(stderr, "numChimera       %u\n", numChimeraTotal);
  fprintf(stderr, "\n");
  fprintf(stderr, "numBubbleUNIQ:   %u\n", numBubbleUNIQ);
  fprintf(stderr, "numBubbleREPT:   %u\n", numBubbleREPT);
  fprintf(stderr, "numBubbleOTHR:   %u\n", numBubbleOTHR);
  fprintf(stderr, "\n");

#if 0
  for (uint32 t=0; t<4; t++) {
    fprintf(stderr, "TOLERANCE %d\n", tolerances[t]);
    fprintf(stderr, "numIncompleteUNIQ   \n");
    fprintf(stderr, "numIncompleteREPT   \n");
    fprintf(stderr, "numUNIQ-REPT-UNIQ   \n");
    fprintf(stderr, "numREPT-UNIQ-REPT   \n");
    fprintf(stderr, "numIllegal\n");
    fprintf(stderr, "\n");
  }
#endif

  return(0);
}
