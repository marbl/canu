
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_BAT/bogusness.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010-2011,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-19 to 2014-DEC-23
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "bogusUtil.H"

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
              int32 _chnBgn, int32 _chnEnd,
              int32 _genIID,
              int32 _genBgn, int32 _genEnd,
              bool  _isReverse,
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

    chnBgn     = _chnBgn;
    chnEnd     = _chnEnd;

    genIID     = _genIID;
    genBgn     = _genBgn;
    genEnd     = _genEnd;
    isReverse  = _isReverse;

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
    if (chnBgn < that.chnBgn)
      return(true);
    if (chnBgn > that.chnBgn)
      return(false);
    return(idlBgn < that.idlBgn);
  };

  char    utgID[32];
  int32   utgIID;

  int32   alignNum;
  int32   alignTotal;

  int32   utgBgn;
  int32   utgEnd;

  int32   chnBgn;    //  Position in the chained reference
  int32   chnEnd;

  int32   genIID;      //  Position in the actual reference (for output)
  int32   genBgn;
  int32   genEnd;
  bool    isReverse;

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



vector<referenceSequence>  refList;
map<string,uint32>         refMap;

vector<idealUnitig>        ideal;

vector<genomeAlignment>    genome;
map<string, int32>         IIDmap;       //  Maps an ID string to an IID.
vector<string>             IIDname;      //  Maps an IID to an ID string.

vector<bogusResult>        results;



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

    if        (strcmp(S[3], "REPT") == 0)
      t = (S.numWords() == 4) ? IDEAL_REPT : IDEAL_UNIQWEAK;

    else if (strcmp(S[3], "UNIQ") == 0)
      t = (S.numWords() == 4) ? IDEAL_UNIQ : IDEAL_REPTUNIQ;

    else
      fprintf(stderr, "Unknown type in '%s'\n", L), exit(1);

#warning ONLY WORKS ON ONE REFERENCE
    ideal.push_back(idealUnitig(S(1), S(2), t));

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

  if ((ideal.bgn    <= cover.chnBgn) &&
      (cover.genEnd <= ideal.end)) {
    alen  = cover.genEnd - cover.chnBgn;
    ulen  = cover.genEnd - cover.chnBgn;
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

  if ((cover.chnBgn <= ideal.bgn) &&
      (ideal.end      <= cover.chnEnd)) {
    alen  = ideal.end - ideal.bgn;
    ulen  = cover.chnEnd - cover.chnBgn;
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

  if ((cover.chnBgn <= ideal.bgn) &&
      (ideal.bgn      <= cover.chnEnd)) {
    alen  = cover.chnEnd - ideal.bgn;
    ulen  = cover.chnEnd - cover.chnBgn;
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

  if ((cover.chnBgn <= ideal.end) &&
      (ideal.end      <= cover.chnEnd)) {
    alen  = ideal.end - cover.chnBgn;
    ulen  = cover.chnEnd - cover.chnBgn;
    ilen  = ideal.end - ideal.bgn;
    ufrac = 100.0 * alen / ulen;
    ifrac = 100.0 * alen / ilen;
    return(true);
  }

  return(false);
}





int
main(int argc, char **argv) {
  uint32   nucmerNamesLen  = 0;
  uint32   snapperNamesLen = 0;
  char    *nucmerNames[1024];
  char    *snapperNames[1024];
  char    *idealName = 0L;
  char    *refName = 0L;
  char    *outputPrefix = 0L;

  FILE    *gffOutput     = 0L;
  FILE    *resultsOutput = 0L;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-nucmer") == 0) {
      nucmerNames[nucmerNamesLen++] = argv[++arg];

    } else if (strcmp(argv[arg], "-snapper") == 0) {
      snapperNames[snapperNamesLen++] = argv[++arg];

    } else if (strcmp(argv[arg], "-ideal") == 0) {
      idealName = argv[++arg];

    } else if (strcmp(argv[arg], "-reference") == 0) {
      refName = argv[++arg];

    } else if (strcmp(argv[arg], "-output") == 0) {
      outputPrefix = argv[++arg];

    } else {
      err++;
    }
    arg++;
  }
  if ((nucmerNamesLen == 0) && (snapperNamesLen == 0L))
    fprintf(stderr, "ERROR: No input matches supplied (either -nucmer or -snapper).\n"), err++;
  if (refName == 0L)
    fprintf(stderr, "ERROR: No referece supplied (-reference).\n"), err++;
  if (idealName == 0L)
    fprintf(stderr, "ERROR: No ideal unitigs supplied (-ideal).\n"), err++;
  if (outputPrefix == 0L)
    fprintf(stderr, "ERROR: No output prefix supplied (-output).\n"), err++;
  if (err) {
    exit(1);
  }


  {
    char   outputName[FILENAME_MAX];

    errno = 0;

    sprintf(outputName, "%s.bogusness", outputPrefix);
    resultsOutput = fopen(outputName, "w");

    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n",
              outputName, strerror(errno)), exit(1);

    sprintf(outputName, "%s.gff3", outputPrefix);
    gffOutput = fopen(outputName, "w");

    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n",
              outputName, strerror(errno)), exit(1);

    fprintf(gffOutput, "##gff-version 3\n");
  }


  loadIdealUnitigs(idealName, ideal);
  loadReferenceSequence(refName, refList, refMap);


  for (uint32 nn=0; nn<nucmerNamesLen; nn++)
    loadNucmer(nucmerNames[nn], genome, IIDmap, IIDname, refList, refMap, 0.0);

  for (uint32 nn=0; nn<snapperNamesLen; nn++)
    loadSnapper(snapperNames[nn], genome, IIDmap, IIDname, refList, refMap, 0.0);


  sort(genome.begin(), genome.end(), byFragmentID);


  for (uint32 bgn=0, lim=genome.size(); bgn<lim; ) {
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
      genome[i].isDeleted = false;
      genome[i].isRepeat  = false;
    }

    for (uint32 i=bgn; i<end; i++) {
      if (genome[i].isDeleted)
        continue;

      for (uint32 j=bgn; j<end; j++) {
        if ((i == j) ||
            (genome[j].isDeleted))
          continue;

        //  If J is contained in I, mark J as deleted.

        assert(genome[i].frgBgn < genome[i].frgEnd);
        assert(genome[j].frgBgn < genome[j].frgEnd);

        if (((genome[i].frgBgn <  genome[j].frgBgn) && (genome[j].frgEnd <  genome[i].frgEnd)) ||
            ((genome[i].frgBgn <= genome[j].frgBgn) && (genome[j].frgEnd <  genome[i].frgEnd)) ||
            ((genome[i].frgBgn <  genome[j].frgBgn) && (genome[j].frgEnd <= genome[i].frgEnd)) ||
            ((genome[i].frgBgn <= genome[j].frgBgn) && (genome[j].frgEnd <= genome[i].frgEnd) && (genome[j].identity <  genome[i].identity)))
          genome[j].isDeleted = true;
      }
    }

    //  Mark repeats -- mapped at exactly the same unitig coords with exactly the same
    //  identitiy.  We'll draw these differently; no parent (not exon-like) and in a different
    //  color.

    for (uint32 i=bgn; i<end; i++) {
      if (genome[i].isDeleted)
        continue;

      for (uint32 j=bgn; j<end; j++) {
        if ((i == j) ||
            (genome[j].isDeleted))
          continue;

        if ((genome[i].frgBgn == genome[j].frgBgn) && (genome[j].frgEnd == genome[i].frgEnd)) {
          //fprintf(stderr, "MARK REPEAT %d %s %d %d -- %d %s %d %d\n",
          //        i, IIDname[genome[i].frgIID].c_str(), genome[i].frgBgn, genome[i].frgEnd,
          //        j, IIDname[genome[j].frgIID].c_str(), genome[j].frgBgn, genome[j].frgEnd);
          genome[i].isRepeat = true;
          genome[j].isRepeat = true;
        }
      }
    }


    //  Write to the gff3 output.

    {
      bool    *foundData = new bool   [refList.size()];
      int32   *spanBgn   = new int32  [refList.size()];
      int32   *spanEnd   = new int32  [refList.size()];
      char   **spanHdr   = new char * [refList.size()];

      for (uint32 i=0; i<refList.size(); i++) {
        foundData[i] = false;
        spanBgn[i]   = 0;
        spanEnd[i]   = 0;
        spanHdr[i]   = NULL;
      }

      //  Find the span of the non-repeat pieces - only used if we aren't 100% repeat.

      for (uint32 i=bgn; i<end; i++) {
        if ((genome[i].isDeleted == false) &&
            (genome[i].isRepeat  == false)) {
          int32 iid = genome[i].genIID;

          if (foundData[iid] == false) {
            foundData[iid] = true;
            spanBgn[iid]   = genome[i].genBgn;  //  OUTPUT: find the span on this reference
            spanEnd[iid]   = genome[i].genEnd;
            spanHdr[iid]   = refList[genome[i].genIID].rsrefName;
          } else {
            spanBgn[iid]   = MIN(spanBgn[iid], genome[i].genBgn);
            spanEnd[iid]   = MAX(spanEnd[iid], genome[i].genEnd);
          }
        }
      }

      //  Are we 100% repeat?

      bool allRepeat = true;

      for (uint32 i=bgn; i<end; i++) {
        //fprintf(stderr, "%s del=%d rep=%d %d-%d %d-%d\n",
        //        IIDname[genome[i].frgIID].c_str(),
        //        genome[i].isDeleted,
        //        genome[i].isRepeat,
        //        genome[i].genBgn, genome[i].genEnd,
        //        genome[i].frgBgn, genome[i].frgEnd);
        if ((genome[i].isDeleted == false) &&
            (genome[i].isRepeat  == false))
          allRepeat = false;
      }

      //  Write unitigs that are entirely repeats without the span, or write
      //  them with a span, but still labeled as repeat.
      //
      if (allRepeat) {
        //  Ignore the span and draw each match individually.
        for (uint32 i=bgn; i<end; i++) {
          if (genome[i].isDeleted)
            continue;

          fprintf(gffOutput, "%s\t.\tbogusness_repeat\t%d\t%d\t.\t%c\t.\tID=REPEAT-%s-%d-%d-%d-%d;Name=%s-%d-%d-REPEAT\n",
                  refList[genome[i].genIID].rsrefName,
                  genome[i].genBgn,
                  genome[i].genEnd,
                  (genome[i].isReverse) ? '-' : '+',
                  IIDname[genome[i].frgIID].c_str(),
                  genome[i].genBgn, genome[i].genEnd,
                  genome[i].frgBgn, genome[i].frgEnd,
                  IIDname[genome[i].frgIID].c_str(), genome[i].frgBgn, genome[i].frgEnd);
        }

      } else {
        //  Draw repeats outside the span as repeats.  Draw those in the span as....what

        //  Dump the span.
        for (uint32 iid=0; iid<refList.size(); iid++)
          if (foundData[iid])
            fprintf(gffOutput, "%s\t.\tbogusness_span\t%d\t%d\t.\t.\t.\tID=%s;Name=SPAN-%s-%d-%d\n",
                    spanHdr[iid], spanBgn[iid], spanEnd[iid],
                    IIDname[genome[bgn].frgIID].c_str(),
                    IIDname[genome[bgn].frgIID].c_str(), spanBgn[iid], spanEnd[iid]);

        //  And the unitigs.
        for (uint32 i=bgn; i<end; i++) {
          int32 iid = genome[i].genIID;

          if (genome[i].isDeleted)
            continue;

          if ((spanBgn[iid] <= genome[i].genBgn) && (genome[i].genEnd <= spanEnd[iid])) {
            //  Alignment within the 'unique' region, draw as a flaw in the unique span if it is labeled repeat
            //
            //  XXX   IMPLEMENT ME!
            //
            if (genome[i].isRepeat) {
              fprintf(stderr, "REPEAT at%d %d\n", genome[i].genBgn, genome[i].genEnd);
            }

            fprintf(gffOutput, "%s\t.\tbogusness_match\t%d\t%d\t.\t%c\t.\tParent=%s;Name=%s-%d-%d\n",
                    refList[genome[i].genIID].rsrefName,
                    genome[i].genBgn,
                    genome[i].genEnd,
                    (genome[i].isReverse) ? '-' : '+',
                    IIDname[genome[i].frgIID].c_str(),
                    IIDname[genome[i].frgIID].c_str(), genome[i].frgBgn, genome[i].frgEnd);
          } else {
            //  Alignment outside the 'unique' region, draw as a repeat.
            fprintf(gffOutput, "%s\t.\tbogusness_repeat\t%d\t%d\t.\t%c\t.\tParent=%s;Name=%s-%d-%d-REPEAT\n",
                    refList[genome[i].genIID].rsrefName,
                    genome[i].genBgn,
                    genome[i].genEnd,
                    (genome[i].isReverse) ? '-' : '+',
                    IIDname[genome[i].frgIID].c_str(),
                    IIDname[genome[i].frgIID].c_str(), genome[i].frgBgn, genome[i].frgEnd);
          }
        }
      }

      delete [] foundData;
      delete [] spanBgn;
      delete [] spanEnd;
      delete [] spanHdr;
    }


    //  Copy whatever is left over to a list of 'covering' alignments.

    for (uint32 i=bgn; i<end; i++) {
      if (genome[i].isDeleted == false)
        cover.push_back(genome[i]);
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

      uint32  frgIID = genome[bgn].frgIID;

      bool explained = false;

      for (uint32 i=0; i<ideal.size(); i++) {
        if (isUnitigContained(frgIID, ideal[i], cover[c], alen, ulen, ilen, ufrac, ifrac)) {
          results.push_back(bogusResult(IIDname[frgIID].c_str(), frgIID,
                                        c+1, cover.size(),
                                        cover[c].frgBgn, cover[c].frgEnd,
                                        cover[c].chnBgn, cover[c].chnEnd,
                                        cover[c].genIID,
                                        cover[c].genBgn, cover[c].genEnd,
                                        cover[c].isReverse,
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
                                        cover[c].chnBgn, cover[c].chnEnd,
                                        cover[c].genIID,
                                        cover[c].genBgn, cover[c].genEnd,
                                        cover[c].isReverse,
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
                                        cover[c].chnBgn, cover[c].chnEnd,
                                        cover[c].genIID,
                                        cover[c].genBgn, cover[c].genEnd,
                                        cover[c].isReverse,
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
                                        cover[c].chnBgn, cover[c].chnEnd,
                                        cover[c].genIID,
                                        cover[c].genBgn, cover[c].genEnd,
                                        cover[c].isReverse,
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
                                        cover[c].chnBgn, cover[c].chnEnd,
                                        cover[c].genIID,
                                        cover[c].genBgn, cover[c].genEnd,
                                        cover[c].isReverse,
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

    fprintf(resultsOutput, "| %s || %d of %d || %d-%d || %s || %d-%d || %s || %s || %05d || %d-%d || %dbp || %.2f%% || %.2f%%\n",
            bi->utgID,
            bi->alignNum, bi->alignTotal,
            bi->utgBgn, bi->utgEnd,
            refList[bi->genIID].rsrefName,
            bi->genBgn, bi->genEnd,
            statuses[bi->status],
            types[bi->type],
            bi->idlNum,
            bi->idlBgn, bi->idlEnd,
            bi->alignLen,
            bi->utgCov, bi->idlCov);
  }

  fclose(gffOutput);
  fclose(resultsOutput);



#if 0 // USELESS BROKEN STUFF


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
          (gi->chnBgn <= gj->chnBgn) &&
          (gj->chnEnd <= gi->chnEnd) &&
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

#endif // USELESS BROKEN STUFF

  return(0);
}
