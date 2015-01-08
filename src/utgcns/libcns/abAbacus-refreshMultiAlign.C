
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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

static char *rcsid = "$Id$";

#include "abAbacus.H"

//#include "AS_UTL_reverseComplement.H"


//#define MIN_SIZE_OF_MANODE 10000

//  Next ID to use for a VAR record
//int32 vreg_id = 0;


//  Stats we need to stash somewhere better
uint32   NumRunsOfGaps = 0;
uint32   NumGaps       = 0;
uint32   NumColumns    = 0;






void
abVarRegion::showConfirmedReads(void) {

  fprintf(stderr, "Confirmed reads=\n");

  for (uint32 j=0; j<nr; j++) {
    if (reads[j].allele_id >= nca)
      continue;

    for (uint32 l=0; l< end-beg; l++)
      fprintf(stderr, "%c", reads[j].bases[l]);

    fprintf(stderr, " allele= %d iid= %d\n", reads[j].allele_id, reads[j].iid);
  }

  fprintf(stderr, "\n");
}

void
abVarRegion::showReadsOfConfirmedAlleles(void) {

  fprintf(stderr, "Confirmed reads=\n");

  for (uint32 j=0; j<nca; j++)
    for (uint32 l=0; l< alleles[j].num_reads; l++)
      fprintf(stderr, "%d  %d \n", alleles[j].read_ids[l], alleles[j].read_iids[l]);

  fprintf(stderr, "\n");
}


void
abVarRegion::showReads(void) {

  fprintf(stderr, "Num_reads= %d vreg.beg= %d vreg.end= %d\n", nr, beg, end);
  fprintf(stderr, "Reads=\n");

  for (uint32 j=0; j<nr; j++) {
    for (uint32 l=0; l< end-beg; l++)
      fprintf(stderr, "%c", reads[j].bases[l]);
    fprintf(stderr, "\n");
  }

  fprintf(stderr, "Ave_qvs= \n");

  for (uint32 j=0; j<nr; j++)
    fprintf(stderr, "%f ", reads[j].ave_qv);

  fprintf(stderr, "\n");
}






void
abAbacus::getReadsForVAR(abVarRegion &region, abColID *cids) {

  //  Allocate memory for reads.  Sadly, no non-default constructor is possible
  //  for array construction.

  //AllocateMemoryForReads(&vreg.reads, vreg.nr, vreg.end - vreg.beg, 0);
  //vreg.reads = new abVarRead [vreg.nr] (vreg.end - vreg.beg, 0);

  region.reads = new abVarRead [region.nr];

  for (uint32 ii=0; ii<region.nr; ii++)
    region.reads[ii].initialize(region.end - region.beg, 0);

  //  Over all columns we care about

  for (uint32 k=region.beg; k<region.end; k++) {
    abColumn          *column = getColumn(cids[k]);
    abColBeadIterator *ci     = createColBeadIterator(cids[k]);

    //  Collect bases and IDs in the coluimn

    for (abBeadID bid = ci->next(); bid.isValid(); bid = ci->next()) {
      abBead      *bead = getBead(bid);
      abSequence  *seq  = getSequence(bead->frag_index);

      if (seq->isRead() == false)
        //  Not a sequencing read, skip it.
        continue;

      uint32  readLoc = 0;

      for (readLoc=0; readLoc<region.nr; readLoc++)
        if (seq->iid == region.iids[readLoc])
          break;

      if (readLoc >= region.nr)
        //  Not a read we care about, skip it.
        continue;

      //  Grab the base, and deterime the QV for it.

      char    base = getBase(bead->baseIdx());
      uint32  qv   = UINT32_MAX;

      if (base != '-') {
        //  An actual base, grab the QV
        qv = getQual(bead->baseIdx()) - '0';

      } else {
        //  An actual gap, grab surrounding QVs.  (But why?  The QV should be set to this already!)
        //  Original comment indicated that this was only applicable to 'boundary gaps', I think,
        //  those with real bases on at least one side.  The code below doesn't capture that.

        abBead *prevBead = getBead(bead->prev);
        abBead *nextBead = getBead(bead->next);

        //  if there is a previous bead, and that bead is not a gap, get the qv.
        //  original comment also said: 'otherwise, it stays QV_FOR_MULTI_GAP'

        if ((prevBead != NULL) &&
            (getBase(prevBead->baseIdx()) != '-'))
          qv = getQual(prevBead->baseIdx()) - '0';

        //  If there is a next bead, and that bead is not a gap, and the qv is worse
        //  than the qv we've found so far, use that.

        if ((nextBead != NULL) &&
            (getBase(nextBead->baseIdx()) != '-') &&
            (getQual(nextBead->baseIdx()) < qv))
          qv = getQual(nextBead->baseIdx()) - '0';
      }

      assert(qv != UINT32_MAX);

      uint32 p = k - region.beg;

      region.reads[readLoc].bases[p] = base;
      region.reads[readLoc].qvs[p]   = qv;
      region.reads[readLoc].iid      = seq->iid;
    }
  }

  //  Reset qvs of internal gaps to MIN(qv_first_gap, qv_last_gap), and compute average qvs.

  for (uint32 k=0; k<region.nr; k++) {
    uint32 m = region.end - region.beg;

    region.reads[k].uglen  = 0;
    region.reads[k].ave_qv = 0;

    for (uint32 i=0; i<m; i++) {

      //  If an actual base, sum the QVs.

      if (region.reads[k].bases[i] != '-') {
        region.reads[k].ave_qv += region.reads[k].qvs[i];
        region.reads[k].uglen++;
        continue;
      }

      //  Otherwise, ... do someting complicated.

      uint32 first_gap = i;
      uint32 first_qv  = region.reads[k].qvs[first_gap];
      uint32 last_gap  = i;

      if ((first_qv == 0) && (i > 0))
        first_qv = region.reads[k].qvs[i-1];

      while ((last_gap < m) && (region.reads[k].bases[last_gap] == '-'))
        last_gap++;

      if ((last_gap == m) || (region.reads[k].bases[last_gap] != '-'))
        last_gap--;

      uint32  last_qv = region.reads[k].qvs[last_gap];

      if ((last_qv == 0) && (last_gap < m-1))
        last_qv = region.reads[k].qvs[last_gap+1];

      uint32  min_qv  = QV_FOR_MULTI_GAP;  //  When both are zero, the final else below

      if ((first_qv != 0) && (last_qv != 0))
        min_qv = (first_qv < last_qv) ? first_qv : last_qv;

      else if ((first_qv == 0) && (last_qv != 0))
        min_qv = last_qv;

      else if ((first_qv != 0) && (last_qv == 0))
        min_qv = first_qv;

      //  Set all QVs in the gap to the minimum, and sum.

      for (uint32 j=first_gap; j <= last_gap; j++) {
        region.reads[k].qvs[j]  = min_qv;
        region.reads[k].ave_qv += min_qv;
      }

      i = last_gap;
    }

    //  Done with the loop, convert the sum to an aveage.

    region.reads[k].ave_qv /= m;
  }

#if 0
  fprintf(stderr, "In GetReads: ave_qvs= ");
  for (uint32 k=0; k<region.nr; k++)
    fprintf(stderr, "%3.1f ", region.reads[k].ave_qv);
  fprintf(stderr, "\n");
#endif
}





static
void
smoothenVariation(double *var, uint32 len, uint32 window) {

  if (window <= 0)
    return;

  double *y = new double [len];

  for (uint32 i=0; i<len; i++) {
    int32  left_win  = 0;
    int32  right_win = 0;

    double sum_var = (var[i] > ZERO_MINUS) ? var[i] : ((var[i] < -1) ? 0 : -var[i]);

    int32 max_left_win  = window / 2;
    int32 max_right_win = window - max_left_win;

    for (int32 j = i-1; ((j >= 0) && (left_win <= max_left_win)); j--) {
      if (var[j] > ZERO_MINUS) {
        // consensus is not gap
        left_win++;
        sum_var += var[j];

      } else if (var[j] > -1 + ZERO_MINUS)
        // consensus is gap, var != 0
        sum_var -= var[j];
    }

    for (int32 j = i+1; ((j < len) && (right_win <= max_right_win)); j++) {
      if (var[j] > ZERO_MINUS) {
        // consensus is not gap
        right_win++;
        sum_var += var[j];

      } else if (var[j] > -1 + ZERO_MINUS)
        // consensus is gap, var != 0
        sum_var -= var[j];
    }

    y[i] = (left_win+right_win > 0) ? sum_var / (left_win + right_win) : var[i];
  }

  for (uint32 i=0; i<len; i++)
    var[i] = y[i];

  delete [] y;
}



                                                                  
bool
abVarRegion::phaseWithPreviousRegion(uint32* &allele_map,
                                     uint32   prev_nca,   uint32  *prev_nca_iid,   uint32  prev_nca_iid_max,
                                     uint32   prev_ncr,   uint32  *prev_ncr_iid,   uint32  prev_ncr_iid_max) {
  bool     is_phased = 0;
  uint32 **allele_matrix;
  uint32  *check;
  uint32   num_reads = 0;

  if (prev_nca == nca && nca >= 2) {
    allele_map       = new uint32   [nca];
    check            = new uint32   [nca];
    allele_matrix    = new uint32 * [prev_nca];
    allele_matrix[0] = new uint32   [prev_nca * nca];

    for (int32 i=1; i<prev_nca; i++)
      allele_matrix[i] = allele_matrix[i-1] + nca;

    assert(prev_ncr < prev_ncr_iid_max);
    assert(prev_nca < prev_nca_iid_max);

    //  Populate the allele_matrix:
    //   columns  = confirmed alleles in the previous VAR record
    //   rows     = confirmed alleles in the current  VAR record
    //   elements = # of reads shared by the prototype and image alleles
    //
    for (uint32 i=0; i<nca; i++) {                      // i = allele id in current VAR record
      for (uint32 j=0; j<alleles[i].num_reads; j++) {
        uint32 read_id=0;

        uint32 l = 0; // allele id in the prev VAR record

        for (uint32 k=0; k<prev_ncr; k++) {
          assert (l < prev_nca);  //  data in prev_nca_iid only valid up to prev_nca!

          if (read_id == prev_nca_iid[l]) { // start of a new allele
            l++;
            read_id = 0;
          }

          read_id++;

          if (prev_ncr_iid[k] == alleles[i].read_iids[j]) {
            assert(l < prev_nca);  //  matrix only allocated up to prev_nca

            allele_matrix[l][i]++;
            num_reads++;
          }
        }
      }
    }

    // Check if    alleles of previous VAR record
    // map well on alleles of current  VAR record. They will do, if:
    // - maximal element in each row in allele_matrix is greater than half
    //   the sum of all elements, and
    // - maximal elements of different rows are located in different
    //   columns
    //
    for (int32 i=0; i<nca; i++)
      allele_map[i] = check[i] = UINT32_MAX;
 
    for (int32 i=0; i<nca; i++) {
      int32 sum = 0;
      int32 max = -1;
      int32 j_best = -1;

      for (int32 j=0; j<prev_nca; j++) {
        sum += allele_matrix[i][j];

        if (max < allele_matrix[i][j]) {
          max = allele_matrix[i][j];
          j_best = j;
        }
      }

      if (2 * max > sum) {
        allele_map[i] = j_best;
        check[j_best] = 1;
      }
    }

    {
      int32 product = 1;

      for (uint32 i=0; i<nca; i++)
        product *= (check[i] + 1);

      if (product > 0)
        is_phased = true;
    }

    // Check if    alleles of current  VAR record
    // map well on alleles of previous VAR record. They will do, if:
    // - maximal element in each col in allele_matrix is greater than half
    //   the sum of all elements, and
    // - maximal elements of different cols are located in different
    //   rows
    //
    if (!is_phased) {
      for (int32 i=0; i<nca; i++) {
        allele_map[i] = check[i] = -1;
      }
      for (int32 j=0; j<nca; j++) { // loop through all columns
        int32 sum = 0, max = -1, i_best = -1;
        for (uint32 i=0; i<nca; i++) {
          sum += allele_matrix[i][j];
          if (max < allele_matrix[i][j]) {
            max = allele_matrix[i][j];
            i_best = i;
          }
        }
        if (2*max > sum) {
          allele_map[i_best] = j;
          check[i_best] = 1;
        }
      }
      {
        int32 product = 1;
        for (int32 i=0; i< nca; i++)
          product *= (check[i]+1);
        if (product > 0)
          is_phased = 1;
      }
    }

    // If still not phased, try another approach for the particular case of diploid genome
    if (!is_phased && nca == 2) {
      int32 max = -1, max2 = -1, i_best = -1, j_best = -1;

      // Find the max and 2nd max elements
      for (uint32 i=0; i<nca; i++) {
        for (uint32 j=0; j<prev_nca; j++) {
          if (max < allele_matrix[i][j]) {
            max2 = max;
            max = allele_matrix[i][j];
            i_best = i;
            j_best = j;
          } else if (max >= allele_matrix[i][j] && max2 < allele_matrix[i][j]) {
            max2 = allele_matrix[i][j];
          }
        }
      }
      if (max > max2) {
        is_phased = 1;
        allele_map[j_best] = i_best;
        allele_map[j_best == 0 ? 1 : 0] = i_best == 0 ? 1 : 0;
      }
    }

    delete [] allele_matrix[0];
    delete [] allele_matrix;
    delete [] check;
  } // if (prev_nca == nca)

  return(is_phased);
}






void
abVarRegion::setConsensusToMajorAllele(abAbacus *abacus, abColID *cids) {

#ifdef DEBUG_VAR_RECORDS
  char *bases = new char [end - beg + 1];

  bases[end-beg]='\0';

  fprintf(stderr, "VAR beg= %d end= %d\n", beg, end);

  OutputReads(stderr, reads, nr, end-beg);
  OutputDistMatrix(stderr, &vreg);
  OutputAlleles(stderr, &vreg);
#endif

  for (uint32 m=0; m<end-beg; m++) {
    abColumn *column = abacus->getColumn(cids[beg+m]);
    abBead   *call   = abacus->getBead(column->callID());

    abacus->setBase(call->baseIdx(), reads[ alleles[0].read_ids[0] ].bases[m]);

#ifdef DEBUG_VAR_RECORDS
    bases[m] = cbase;
#endif
  }

#ifdef DEBUG_VAR_RECORDS
  fprintf(stderr, "beg= %d end= %d Num_reads= %d \nConsensus= \n%s\n", beg, end, nr, bases);
  fprintf(stderr, "Reads=\n");
  for (uint32 read_id=0; read_id<nr; read_id++) {
    for (uint32 m=0; m<end-beg; m++)
      fprintf(stderr, "%c", reads[read_id].bases[m]);
    fprintf(stderr, "   allele=%d \n", reads[read_id].allele_id);
  }

  delete [] bases;
#endif
}







//  Compute the gapped and ungapped distance, and save the smaller of the two
void
abVarRegion::populateDistanceMatrix(void) {
  char  *ugread1 = new char [end - beg];
  char  *ugread2 = new char [end - beg];

  for (uint32 i=0; i<nr; i++) {
    for (uint32 j=i; j<nr; j++) {
      uint32  gapped_dist = 0;

      uint32  uglen1 = 0;
      uint32  uglen2 = 0;

      for (uint32 k=0; k<end - beg; k++) {
        if (reads[i].bases[k] != reads[j].bases[k])
          gapped_dist++;

        if (reads[i].bases[k] != '-')
          ugread1[uglen1++] = reads[i].bases[k];

        if (reads[j].bases[k] != '-')
          ugread2[uglen2++] = reads[j].bases[k];
      }

      uint32  uglen = (uglen1 < uglen2) ? uglen2 : uglen1;

      uint32  ungapped_dist = 0;

      for (uint32 k=0; k<uglen; k++)
        if (((k < uglen1) && (k < uglen2) && (ugread1[k] != ugread2[k])) ||
            ((k < uglen1) && (k >= uglen2)) ||
            ((k < uglen2) && (k >= uglen1)))
          ungapped_dist++;

      dist_matrix[i][j] =
        dist_matrix[j][i] = (gapped_dist < ungapped_dist) ? gapped_dist : ungapped_dist;
    }
  }

  delete [] ugread1;
  delete [] ugread2;
}







//  detect allele and split reads between the alleles

void
abVarRegion::clusterReads(void) {

  // Allocate memory for alleles

  na      = 0;
  alleles = new abVarAllele [nr];

  for (uint32 ii=0; ii<nr; ii++)
    alleles[ii].initialize(nr);

  // Initialize alleles

  // Process zero elements first
  for (uint32 row=0; row<nr; row++) {
    for (uint32 col=row+1; col<nr; col++) {
      if (dist_matrix[row][col] != 0)
        continue;

      if ((reads[row].allele_id < 0) && (reads[col].allele_id < 0)) {

        // New allele
        reads[row].allele_id     = na;
        reads[col].allele_id     = na;
        alleles[na].weight       = ROUND(reads[row].ave_qv) + ROUND(reads[col].ave_qv);
        alleles[na].uglen        = reads[row].uglen;
        alleles[na].read_ids[0]  = row;
        alleles[na].read_ids[1]  = col;
        alleles[na].read_iids[0] = reads[row].iid;
        alleles[na].read_iids[1] = reads[col].iid;
        alleles[na].num_reads    = 2;
        alleles[na].id           = na;

        na++;

      } else if ((reads[row].allele_id < 0) && (reads[col].allele_id >= 0)) {
        // Already existing allele
        uint32 aid = reads[col].allele_id;
        reads[row].allele_id     = aid;
        alleles[aid].weight     += ROUND(reads[row].ave_qv);

        uint32 anr = alleles[aid].num_reads;
        alleles[aid].read_ids[anr]  = row;
        alleles[aid].read_iids[anr] = reads[row].iid;

        alleles[aid].num_reads++;

      } else if ((reads[row].allele_id >= 0) && (reads[col].allele_id < 0)) {
        // Already existing allele
        uint32 aid = reads[row].allele_id;
        reads[col].allele_id    = aid;
        alleles[aid].weight    += ROUND(reads[col].ave_qv);

        uint32 anr = alleles[aid].num_reads;
        alleles[aid].read_ids[anr]  = col;
        alleles[aid].read_iids[anr] = reads[col].iid;

        alleles[aid].num_reads++;
      }
    }
  }

  nca = na;

  //Now process the remaining reads; assign each to its "own" allele

  for (uint32 row=0; row<nr; row++) {
    if (reads[row].allele_id < 0) {
      // New allele
      reads[row].allele_id     = na;
      alleles[na].weight       = ROUND(reads[row].ave_qv);
      alleles[na].uglen        = reads[row].uglen;
      alleles[na].read_ids[0]  = row;
      alleles[na].read_iids[0] = reads[row].iid;
      alleles[na].num_reads    = 1;
      alleles[na].id           = na;

      na++;
    }
  }
}







// Reverse sort confirmed alleles by ungapped length
static
void
sortAllelesByLength(abVarAllele *alleles, uint32 num_alleles, abVarRead *reads) {

  for (uint32 i=0; i<num_alleles; i++) {
    uint32 best_uglen = alleles[i].uglen;
    uint32 best_id    = UINT32_MAX;

    for (uint32 j=i+1; j<num_alleles; j++) {
      if (best_uglen  < alleles[j].uglen ) {
        best_uglen = alleles[j].uglen;
        best_id    = j;
      }
    }

    if (best_id != UINT32_MAX) {
      abVarAllele temp     = alleles[i];  //  probably calls destructor and kills arrays
      alleles[i]           = alleles[best_id];
      alleles[best_id]     = temp;
    }
  }

  //  Update allele_id of reads

  for (uint32 i=0; i<num_alleles; i++)
    for (uint32 j=0; j<alleles[i].num_reads; j++)
      reads[ alleles[i].read_ids[j] ].allele_id = i;
}


// Reverse sort by weight
static
void
sortAllelesByWeight(abVarAllele *alleles, uint32 num_alleles, abVarRead *reads) {

  for (uint32 i=0; i<num_alleles; i++) {
    uint32 best_weight = alleles[i].weight;
    uint32 best_id     = UINT32_MAX;

    for (uint32 j=i+1; j<num_alleles; j++) {
      if (best_weight < alleles[j].weight) {
        best_weight = alleles[j].weight;
        best_id     = j;
      }
    }

    if (best_id != UINT32_MAX) {
      abVarAllele temp = alleles[i];  //  probably calls destructor and kills arrays
      alleles[i]       = alleles[best_id];
      alleles[best_id] = temp;
    }
  }

  //  Update allele_id of reads

  for (uint32 i=0; i<num_alleles; i++)
    for (uint32 j=0; j<alleles[i].num_reads; j++)
      reads[ alleles[i].read_ids[j] ].allele_id = i;
}


// Sort confirmed alleles according to their mapping
// between two "phased" VAR records
static
void
sortAllelesByMapping(abVarAllele *alleles, uint32 nca, abVarRead *reads, uint32 *allele_map) {

  for (uint32 i=0; i<nca; i++) {
    uint32  j = 0;      // j is id of the allele that should be at i-th place

    for (j=0; j<nca; j++)
      if (allele_map[j] == i)
        break;

    for (uint32 k=i; k<nca; k++) {
      if (alleles[k].id == j) {
        abVarAllele temp = alleles[i];
        alleles[i]       = alleles[k];
        alleles[k]       = temp;
        break;
      }
    }
  }

  // Update allele_ids
  for (uint32 i=0; i<nca; i++)
    alleles[i].id = i;

  for (uint32 i=0; i<nca; i++)
    for (uint32 j=0; j<alleles[i].num_reads; j++)
      reads[ alleles[i].read_ids[j] ].allele_id = i;
}







void
abAbacus::refreshMultiAlign(abMultiAlignID  mid,
                            uint32          quality,
                            uint32          splitAlleles,
                            uint32          smoothWindow,
                            uint32          doPhasing,
                            uint32         *nvars,
                            abVarRead     **v_list,
                            uint32          make_v_list,   //  0, 1 or 2
                            bool            getScores) {  //  0, 1 or 2

  // refresh columns from cid to end
  // if quality == -1, don't recall the consensus base

  if (nvars != NULL)
    *nvars = 0;

  abMultiAlign *ma = getMultiAlign(mid);

  ma->columnList.clear();

  if (ma->first.isValid() == false);
    return;

  abVarRegion  vreg;

  uint32       varLen = 1024;
  double      *varf   = new double  [varLen];
  abColID     *cids   = new abColID [varLen];

  uint32   maxPrev   = 1024;
  char    *prevBases = (getScores > 0) ? new char   [maxPrev] : NULL;
  uint32  *prevIDs   = (getScores > 0) ? new uint32 [maxPrev] : NULL;

  //
  // Calculate variation as a function of position in MANode.
  //

  // Call consensus using all alleles; the goal is to detect a variation at a given position

  uint32  prev_nr           = 0;

  // Variables used to phase VAR records
  uint32  prev_nca          = 0;     // valid size of array prev_nca_iid
  uint32  prev_ncr          = 0;     // valid size of array prev_ncr_iid
  uint32  prev_nca_iid_max  = 4096;  // allocated size of arrays; err on the large side
  uint32  prev_ncr_iid_max  = 4096;  //   hopefully avoiding reallocs.
  uint32 *prev_nca_iid      = NULL;  // number of reads in 10 first confirmed alleles
  uint32 *prev_ncr_iid      = NULL;  // iids of the first 100 reads, rev. sorted by allele

  uint32   index = 0;
  abColID  cid   = ma->first;

  while (cid.isValid()) {
    abColumn  *column = getColumn(cid);

    if (quality != -2) {
      if (index >= varLen)
        resizeArrayPair(varf, cids, varLen, varLen, varLen + 1024);

      //BaseCall(cid, quality, varf[index], &vreg, -1, cbase, getScores, opp);
      char cbase = baseCall(vreg, cid, quality, varf[index], -1, getScores, splitAlleles, smoothWindow);

      if (cbase == '-')
        NumGaps++;

      cids[index] = cid;
    }

    column->ma_position = index;

    ma->columnList.push_back(cid);  //  Does this need to be stored?

    // sanity check
    if (index > 0) {
      abColID   prev = ma->columnList[index-1];
      abColumn *pcol = getColumn(prev);

      if ((prev != column->prev) || (pcol->next != column->lid))
        fprintf(stderr, "RefreshMANode column relationships violated");

      assert(prev == column->prev);
      assert(pcol->next != column->lid);
    }

    if (getScores == true) {
#if 0
      fprintf(stderr, "vreg.nb=%d vreg.curr_bases=", vreg.nb);
      for (i=0; i<vreg.nb; i++)
        fprintf(stderr, "%c", vreg.curr_bases[i]);
      fprintf(stderr, " prev_nr=%d prevBases=", prev_nr);
      for (i=0; i<prev_nr; i++)
        fprintf(stderr, "%c", prevBases[i]);
      fprintf(stderr, " NumRunsOfGaps=%d \nvreg.iids= ", NumRunsOfGaps);
      for (i=0; i<vreg.nb; i++)
        fprintf(stderr, "%d ", vreg.iids[i]);
      fprintf(stderr, "\n");
      fprintf(stderr, "prevIDs= ");
      for (i=0; i<prev_nr; i++)
        fprintf(stderr, "%d ", prevIDs[i]);
      fprintf(stderr, "\n");
#endif

      //  Update stats: number of runs of gaps and number of gaps.

      for (int32 i=0; i<prev_nr; i++) {
        if (prevBases[i] == '-')
          continue;

        for (int32 j=0; j<vreg.nb; j++) {
          if (vreg.curr_bases[j] != '-')
            continue;

          if (prevIDs[i] == vreg.iids[j])
            NumRunsOfGaps++;
        }
      }

      //

      if (vreg.nb > maxPrev) {
        resizeArrayPair(prevBases, prevIDs, maxPrev, maxPrev, vreg.nb + 128);
        //maxPrev =  vreg.nb;
        //prevBases = (char  *)safe_realloc(prevBases, maxPrev*sizeof(char));
        //prevIDs   = (int32 *)safe_realloc(prevIDs, maxPrev*sizeof(int32));
      }

      prev_nr = vreg.nb;

      for (uint32 i=0; i<vreg.nb; i++) {
        prevBases[i] = vreg.curr_bases[i];
        prevIDs[i]  = vreg.iids[i];
      }
    }  //  getScore == true

    cid = column->next;
    index++;
  }

  if (getScores == true)
    NumColumns += index;

  if ((splitAlleles == false) ||
      (quality     <= 0) ||
      (make_v_list == false)) {
    vreg.curr_bases.clear();
    vreg.iids.clear();

    delete [] varf;
    delete [] cids;

    delete [] prevBases;
    delete [] prevIDs;

    return;
  }

  assert(make_v_list == 1 || nvars  != NULL);
  assert(make_v_list == 1 || v_list != NULL);

  // Proceed further only if accurate base calls are needed

  varLen = index -1;

  double *svarf  = new double [varLen];

  for (uint32 i=0; i<varLen; i++) {
    svarf[i] = varf[i];

    if (varf[i] < ZERO_MINUS)
      varf[i] = (varf[i] < -1.) ? 0. : -varf[i];

#if 0
    if (varf[i] > 0)
      fprintf(stderr, "i= %d varf= %f\n", i, varf[i]);
#endif
  }

  smoothenVariation(svarf, varLen, smoothWindow);

  prev_nca_iid = new uint32 [prev_nca_iid_max];
  prev_ncr_iid = new uint32 [prev_ncr_iid_max];

  for (uint32 i=0; i<varLen; i++) {
    if (svarf[i] == 0)
      continue;

    // Process a region of variation
    //
    // smoothWindow <  0 means no columns with variation will be grouped
    //                         into regions of variation
    // smoothWindow == 0 means only immediately adjacent columns with
    //                         variation will be grouped into regions of variation
    uint32  *conf_read_iids = NULL;
    double   fict_var;

    vreg.beg = i;
    vreg.end = i;

    uint32 svbeg = i;
    uint32 svend = i;

    // Set the beginning of unsmoothed VAR region to the 1st position with varf != 0

    if (smoothWindow > 0)
      while ((vreg.beg < varLen - 1) && (DBL_EQ_DBL(varf[vreg.beg], 0.0)))
        vreg.beg++;

    // Set the end of smoothed VAR region to the 1st position with svarf == 0

    if (smoothWindow > 0)
      while ((svend < varLen) && (svarf[svend] > ZERO_PLUS))
        svend++;

    else if (svend < varLen)
      svend++;

    // Set the end of unsmoothed VAR region to the 1st position with varf == 0

    vreg.end = vreg.beg;

    if (smoothWindow > 0)
      while (vreg.end < varLen && varf[vreg.end] > ZERO_PLUS)
        vreg.end++;

    else if (vreg.end < varLen)
      vreg.end++;

    assert(vreg.beg < vreg.end);

    // Store iids of all the reads in current region

    vreg.iids.clear();
      
    // Get all the read iids
    // Calculate the total number of reads, vreg.nr (corresponding to any allele)

    for (uint32 j=vreg.beg; j<vreg.end; j++) {
      abColumn          *column = getColumn(cids[j]);
      abBead            *call   = getBead(column->call);
      abColBeadIterator *ci     = createColBeadIterator(cids[j]);
        
      for (abBeadID bid = ci->next(); bid.isValid(); bid=ci->next()) {
        abBead      *bead = getBead(bid);
        char         base = getBase(bead->baseIdx());
        abSequence  *seq  = getSequence(bead->frag_index);

        if (base == 'N')
          //  Not a base we care about.
          continue;

        if (seq->isRead() == false)
          //  Not a sequence we care about.
          continue;
          
        bool  newRead = true;

        for (uint32 i=0; i<vreg.iids.size(); i++)
          if (vreg.iids[i] == seq->iid)
            newRead = false;

        if (newRead == false)
          //  We've already seen this read
          continue;

        //  This used to set the allocated but unused iids to -1.

        vreg.iids.push_back(seq->iid);
      }
    }

    vreg.nr = vreg.iids.size();

    //  BPW:  Hmmm.  How do we have no reads here?
    if (vreg.nr == 0) {
      fprintf(stderr, "no reads for vreg beg=%d end=%d SKIPPING IT (bug?)\n", vreg.beg, vreg.end);
      continue;
    }

    getReadsForVAR(vreg, cids);

#ifdef DEBUG_VAR_RECORDS
    show_reads(&vreg);
#endif

    // Calculate a sum of qvs for each read within a variation region
    // Populate the distance matrix

    vreg.allocateDistanceMatrix(-1);
    vreg.populateDistanceMatrix();

#ifdef DEBUG_VAR_RECORDS
    OutputDistMatrix(stderr, &vreg);
#endif

    // Populate vreg.alleles array
    // Determine the best allele and the number of reads in this allele
    vreg.clusterReads();

#if 0
    if (vreg.nca < 2) {
      fprintf(stderr, "nca= %d\n", vreg.nca);
      show_reads(&vreg);
    }
#endif

    bool     is_phased  = false;
    uint32  *allele_map = NULL;
    
    if ((make_v_list == 2) && (doPhasing))
      is_phased = vreg.phaseWithPreviousRegion(allele_map,
                                               prev_nca, prev_nca_iid, prev_nca_iid_max,
                                               prev_ncr, prev_ncr_iid, prev_ncr_iid_max);

    if (is_phased)
      sortAllelesByMapping(vreg.alleles, vreg.nca, vreg.reads, allele_map);
    else
      sortAllelesByWeight(vreg.alleles, vreg.na, vreg.reads);

    delete [] allele_map;

    // Create a list of iids of confirmed reads
    {
      uint32 start_num = 0;
      uint32 num_conf_reads = 0;

      for (uint32 i1=0; i1<vreg.nca; i1++)
        num_conf_reads += vreg.alleles[i1].num_reads;

      conf_read_iids = new uint32 [num_conf_reads];

      for (uint32 i1=0; i1<vreg.nca; i1++) {
        for (uint32 j=0; j<vreg.alleles[i1].num_reads; j++)
          conf_read_iids[start_num+j] = vreg.alleles[i1].read_iids[j];
        start_num += vreg.alleles[i1].num_reads;
      }
    }

#ifdef DEBUG_VAR_RECORDS
    show_confirmed_reads(&vreg);
    show_reads_of_confirmed_alleles(&vreg);
#endif

    // Update the info about the previous VAR record :
    // prev_nca, prev_nca_iid and prev_ncr_iid
    {

      //  we should have space enough for these, right?  If
      //  not, count and reallocate new space -- be sure to
      //  clear it!

      prev_nca = vreg.nca;
      prev_ncr = 0;

      for (uint32 iv = 0; iv < vreg.nca; iv++)
        prev_ncr += vreg.alleles[iv].num_reads;

      if (prev_nca_iid_max < prev_nca) {
        while (prev_nca_iid_max < prev_nca)
          prev_nca_iid_max *= 2;

        delete [] prev_nca_iid;
        prev_nca_iid = new uint32 [prev_nca_iid_max];
      }

      if (prev_ncr_iid_max < prev_ncr) {
        while (prev_ncr_iid_max < prev_ncr)
          prev_ncr_iid_max *= 2;

        delete [] prev_ncr_iid;
        prev_ncr_iid = new uint32 [prev_ncr_iid_max];
      }

      memset(prev_nca_iid, 0, sizeof(int32) * prev_nca_iid_max);
      memset(prev_ncr_iid, 0, sizeof(int32) * prev_ncr_iid_max);

      for (uint32 iv=0, kv=0; iv < vreg.nca; iv++) {
        prev_nca_iid[iv] = vreg.alleles[iv].num_reads;

        for (uint32 jv=0; jv < vreg.alleles[iv].num_reads; jv++)
          prev_ncr_iid[kv++] = vreg.alleles[iv].read_iids[jv];
      }

#ifdef DEBUG_VAR_RECORDS
      OutputAlleles(stderr, &vreg);
#endif

      /* Store variations in a v_list */
      if (make_v_list == 2)
#warning NOT POPULATING OUTPUT
        //vreg.populateVARRecord(is_phased, cids, *nvars, min_len_vlist, *v_list, vreg, opp, getScores, conf_read_iids);
        ;
      else
        vreg.setConsensusToMajorAllele(this, cids);
    }

    if (smoothWindow > 0)
      i = svend;

#if 0
    for (uint32 j=0; j<vreg.nr; j++) {
      delete [] vreg.dist_matrix[j];
      delete [] vreg.reads[j].bases;
      delete [] vreg.reads[j].qvs;
      delete [] vreg.alleles[j].read_ids;
      delete [] vreg.alleles[j].read_iids;
    }

    vreg.nr = 0;

    delete [] vreg.dist_matrix;
    delete [] vreg.reads;
    delete [] vreg.alleles;
#endif

    delete [] conf_read_iids;
  }

  delete [] varf;
  delete [] svarf;
  delete [] cids;

  delete [] prevBases;
  delete [] prevIDs;

  delete [] prev_nca_iid;
  delete [] prev_ncr_iid;
}


