
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

static char *rcsid = "$Id: RefreshMANode.c,v 1.2 2009-06-10 18:05:13 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"
#include "MicroHetREZ.h"
#include "AS_UTL_reverseComplement.h"


#define MIN_SIZE_OF_MANODE 10000
#define MIN_ALLOCATED_DEPTH 100

//  Next ID to use for a VAR record
int vreg_id = 0;

static
void
UpdateScores(VarRegion vreg, char *cbase, int nca) {
  int i, j;

  for (i=0; i<nca; i++)
      for (j=i+1; j<nca; j++)
          if (cbase[i] != cbase[j])
            NumAAMismatches++;
}

static
void
UpdateScoreNumRunsOfGaps(VarRegion vreg, int prev_nr, char *prev_bases,
                         int32 *prev_iids, int get_scores) {
  int i, j;

  // Updating count of stretches of gaps
  for (i=0; i<prev_nr; i++) {
    if (prev_bases[i] == '-')
      continue;

    for (j=0; j<vreg.nb; j++) {
      if (vreg.curr_bases[j] != '-')
        continue;

      if (prev_iids[i] == vreg.iids[j]) {
          if (get_scores == 1)
            NumRunsOfGapsInUnitigReads++;
          else if (get_scores == 2)
            NumRunsOfGapsInContigReads++;
        }
    }
  }
}

static
void
UpdateScoreNumGaps(char cbase, int get_scores) {
  if (cbase == '-') {
      if (get_scores == 1)
        NumGapsInUnitigs++;
      else if (get_scores == 2)
        NumGapsInContigs++;
    }
}

static
int
IsNewRead(int32 iid, int32 *iid_list, int nr) {
  int i;

  for (i=0; i<nr; i++)
      if (iid_list[i] == iid)
        return 0;

  return 1;
}

static
void
GetReadIidsAndNumReads(int cid, VarRegion  *vreg) {
  int      cind;
  int16    bi;
  Column  *column=GetColumn(columnStore,cid);
  Bead    *call = GetBead(beadStore, column->call);
  Bead    *bead;
  int32    bid;
  int32    iid;
  FragType type;
  ColumnBeadIterator ci;
  int nr=0, max_nr=100;

  CreateColumnBeadIterator(cid, &ci);

  while ( (bid = NextColumnBead(&ci)) != -1 ) {
      char base;

      bead =  GetBead(beadStore,bid);
      base = *Getchar(sequenceStore,bead->soffset);
      if ( base == 'N' )
        continue;
      type = GetFragment(fragmentStore,bead->frag_index)->type;
      iid  = GetFragment(fragmentStore,bead->frag_index)->iid;

      if ((type == AS_READ) ||
          (type == AS_EXTR) ||
          (type == AS_TRNR)) {
          if (IsNewRead(iid, vreg->iids, vreg->nr)) {
            if (vreg->nr == vreg->max_nr) {
              int l;
              vreg->max_nr += MIN_ALLOCATED_DEPTH;
              vreg->iids = (int32 *)safe_realloc(vreg->iids, vreg->max_nr*sizeof(int32));

              for (l=vreg->nr; l<vreg->max_nr; l++)
                vreg->iids[l] = -1;
            }
            vreg->iids[vreg->nr] = iid;
            vreg->nr++;
          }
        }
    }
}


static
void
GetReadsForVARRecord(Read *reads, int32 *iids, int32 nvr,
                     int beg, int end, int32 *cids) {
  int k;

  for (k=beg; k<end; k++) {
      Column *column=GetColumn(columnStore,cids[k]);
      Bead  *bead;
      int32  bid;
      int32  iid;
      FragType type;
      ColumnBeadIterator ci;
      char  base, qv;
      int i, j;
      int    nr=0, max_nr=100;

      CreateColumnBeadIterator(cids[k], &ci);

      // Collect bases and usids in the coluimn
      while ( (bid = NextColumnBead(&ci)) != -1 ) {
          bead = GetBead(beadStore,bid);
          type = GetFragment(fragmentStore,bead->frag_index)->type;
          iid  = GetFragment(fragmentStore,bead->frag_index)->iid;

          if ((type == AS_READ)
              //              || (type == AS_EXTR)
              //              || (type == AS_TRNR)
              )
            {
              base = *Getchar(sequenceStore,bead->soffset);
              if (base != '-') {
                  qv = (int)(*Getchar(qualityStore,bead->soffset)-'0');
                } else {
 // set qvs of boundary gaps to qvs of adjacent bases
                  Bead *prev_bead = GetBead(beadStore, bead->prev);
                  Bead *next_bead = GetBead(beadStore, bead->next);
                  qv = 0;
                  if (prev_bead != NULL) {
                      char prev_base=*Getchar(sequenceStore, prev_bead->soffset);
                      int  prev_qv  =(int)(*Getchar(qualityStore,
                                                    prev_bead->soffset)-'0');
                      if (prev_base != '-') { qv = prev_qv; }
                      // otherwise, it stays QV_FOR_MULTI_GAP
                    }
                  if (next_bead != NULL) {
                      char next_base=*Getchar(sequenceStore, next_bead->soffset);
                      int  next_qv  =(int)(*Getchar(qualityStore,
                                                    next_bead->soffset)-'0');
                      if (next_base != '-' &&
                          (qv == 0 || qv > next_qv))
                        qv = next_qv;
                    }
                }
              iid  =  GetFragment(fragmentStore,bead->frag_index)->iid;

              for (i=0; i<nvr; i++)
                if (iid == iids[i])
                  break;

              if (i >= nvr)
                continue;

              reads[i].bases[k-beg] = base;
              reads[i].qvs[k-beg] = qv;
              reads[i].iid = iid;
            }
        }
    }

  // Reset qvs of internal gaps to MIN(qv_first_gap, qv_last_gap);
  // Compute ave_qvs
  for (k=0; k<nvr; k++) {
      int i, j, m = end-beg;
      reads[k].uglen = 0;
      reads[k].ave_qv = 0.;
      for (i=0; i<m; i++) {
          if (reads[k].bases[i] != '-') {
              reads[k].ave_qv += (double)reads[k].qvs[i];
              reads[k].uglen++;
            } else {
              // gap
              int first_gap = i;
              int first_qv  = reads[k].qvs[first_gap];
              int last_gap  = i;
              int last_qv;
              int min_qv;
              if (first_qv == 0 && i>0) first_qv = reads[k].qvs[i-1];
              while (last_gap<m && reads[k].bases[last_gap] == '-')
                last_gap++;
              if (last_gap == m || reads[k].bases[last_gap] != '-')
                last_gap--;
              last_qv = reads[k].qvs[last_gap];
              if (last_qv == 0 && last_gap<m-1)
                last_qv = reads[k].qvs[last_gap+1];
              if (first_qv != 0 && last_qv  != 0)
                min_qv = (first_qv < last_qv) ? first_qv : last_qv;
              else if (first_qv == 0 && last_qv != 0)
                min_qv = last_qv;
              else if (first_qv != 0 && last_qv == 0)
                min_qv = first_qv;
              else // both == 0
                min_qv = QV_FOR_MULTI_GAP;
              for (j=first_gap; j<=last_gap; j++)
                {
                  reads[k].qvs[j] = min_qv;
                  reads[k].ave_qv += (double)min_qv;
                }
              i = last_gap;
            }
        }
      reads[k].ave_qv /= (double)m;
    }

#if 0
  fprintf(stderr, "In GetReads: ave_qvs= ");
  for (k=0; k<nvr; k++)
    fprintf(stderr, "%3.1f ", reads[k].ave_qv);
  fprintf(stderr, "\n");
#endif
}












static
void
SmoothenVariation(double *var, int len, int window) {
  int i;
  double *y = (double *)safe_malloc(len * sizeof(double));

  if (window <= 0)
    {
      safe_free(y);
      return;
    }
  for (i=0; i<len; i++)
    {
      int j, left_win=0, right_win=0;
      double sum_var = (var[i]>ZERO_MINUS)?var[i]:((var[i]<-1.)?0.:-var[i]);
      int max_left_win = window/2;
      int max_right_win = window - max_left_win;

      j = i-1;
      while (j>=0 && left_win<=max_left_win)
        {
          if (var[j] > ZERO_MINUS)   // consensus is not gap
            {
              left_win++;
              sum_var += var[j];
            }
          else if (var[j]>-1+ZERO_MINUS) // consensus is gap, var != 0
            sum_var -= var[j];

          j--;
        }
      j = i+1;
      while (j<len && right_win<=max_right_win)
        {
          if (var[j] > ZERO_MINUS)  // consensus is not gap
            {
              right_win++;
              sum_var += var[j];
            }
          else if (var[j]>-1+ZERO_MINUS) // consensus is gap, var != 0
            sum_var -= var[j];
          j++;
        }
      y[i] = (left_win+right_win > 0) ?
        sum_var/(double)(left_win+right_win) : var[i];
    }
  for (i=0; i<len; i++)
    {
      var[i] = y[i];
    }
  safe_free(y);
}





static
int
PhaseWithPrevVreg(int32 nca, Allele *alleles, Read *reads, int32 **allele_map,
                  int32 prev_nca, int32 *prev_nca_iid, int32 prev_nca_iid_max,
                  int32 prev_ncr, int32 *prev_ncr_iid, int32 prev_ncr_iid_max) {
  int   i, j, k, l;
  int   is_phased = 0;
  int32 **allele_matrix;
  int32  *check;
  int   num_reads = 0;

  if (prev_nca == nca && nca >= 2)
    {
      *allele_map = (int32 *)safe_calloc(nca, sizeof(int32));
      check      = (int32 *)safe_calloc(nca, sizeof(int32));
      allele_matrix    = (int32 **)safe_calloc(prev_nca, sizeof(int32 *));
      allele_matrix[0] = (int32  *)safe_calloc(prev_nca * nca, sizeof(int32));
      for (i=1; i<prev_nca; i++)
        allele_matrix[i] = allele_matrix[i-1] + nca;

      assert(prev_ncr < prev_ncr_iid_max);
      assert(prev_nca < prev_nca_iid_max);

      /* Populate the allele_matrix:
       *   columns  = confirmed alleles in the previous VAR record
       *   rows     = confirmed alleles in the current  VAR record
       *   elements = # of reads shared by the prototype and image alleles
       */
      for (i=0; i<nca; i++)   // i = allele id in current VAR record
        {
          for (j=0; j<alleles[i].num_reads; j++)
            {
              int read_id=0;
              l = 0; // allele id in the prev VAR record

              for (k=0; k<prev_ncr; k++)
                {
                  assert (l < prev_nca);  //  data in prev_nca_iid only valid up to prev_nca!

                  if (read_id == prev_nca_iid[l])  // start of a new allele
                    {
                      l++;
                      read_id = 0;
                    }
                  read_id++;

                  if (prev_ncr_iid[k] == alleles[i].read_iids[j])
                    {
                      assert(l < prev_nca);  //  matrix only allocated up to prev_nca
                      allele_matrix[l][i]++;
                      num_reads++;
                    }
                }
            }
        }

      /* Check if    alleles of previous VAR record
       * map well on alleles of current  VAR record. They will do, if:
       * - maximal element in each row in allele_matrix is greater than half
       *   the sum of all elements, and
       * - maximal elements of different rows are located in different
       *   columns
       */
      for (i=0; i<nca; i++)
        {
          (*allele_map)[i] = check[i] = -1;
        }
      for (i=0; i<nca; i++)
        {
          int sum = 0, max = -1, j_best = -1;
          for (j=0; j<prev_nca; j++)
            {
              sum += allele_matrix[i][j];
              if (max < allele_matrix[i][j])
                {
                  max = allele_matrix[i][j];
                  j_best = j;
                }
            }
          if (2*max > sum)
            {
              (*allele_map)[i] = j_best;
              check[j_best] = 1;
            }
        }

      {
        int product = 1;
        for (i=0; i< nca; i++)
          product *= (check[i]+1);
        if (product > 0)
          is_phased = 1;
      }

      /* Check if    alleles of current  VAR record
       * map well on alleles of previous VAR record. They will do, if:
       * - maximal element in each col in allele_matrix is greater than half
       *   the sum of all elements, and
       * - maximal elements of different cols are located in different
       *   rows
       */
      if (!is_phased)
        {
          for (i=0; i<nca; i++)
            {
              (*allele_map)[i] = check[i] = -1;
            }
          for (j=0; j<nca; j++) // loop through all columns
            {
              int sum = 0, max = -1, i_best = -1;
              for (i=0; i<nca; i++)
                {
                  sum += allele_matrix[i][j];
                  if (max < allele_matrix[i][j])
                    {
                      max = allele_matrix[i][j];
                      i_best = i;
                    }
                }
              if (2*max > sum)
                {
                  (*allele_map)[i_best] = j;
                  check[i_best] = 1;
                }
            }
          {
            int product = 1;
            for (i=0; i< nca; i++)
              product *= (check[i]+1);
            if (product > 0)
              is_phased = 1;
          }
        }

      /* If still not phased, try another approach for the particular case of diploid genome */
      if (!is_phased && nca == 2)
        {
          int max = -1, max2 = -1, i_best = -1, j_best = -1;

          /* Find the max and 2nd max elements */
          for (i=0; i<nca; i++)
            {
              for (j=0; j<prev_nca; j++)
                {
                  if (max < allele_matrix[i][j])
                    {
                      max2 = max;
                      max = allele_matrix[i][j];
                      i_best = i;
                      j_best = j;
                    }
                  else if (max >= allele_matrix[i][j] && max2 < allele_matrix[i][j])
                    {
                      max2 = allele_matrix[i][j];
                    }
                }
            }
          if (max > max2)
            {
              is_phased = 1;
              (*allele_map)[j_best] = i_best;
              (*allele_map)[j_best == 0 ? 1 : 0] = i_best == 0 ? 1 : 0;
            }
        }

      safe_free(allele_matrix[0]);
      safe_free(allele_matrix);
      safe_free(check);
    } /* if (prev_nca == nca) */

  return is_phased;
}







#ifdef DEBUG_VAR_REGIONS
static
void
show_confirmed_reads(VarRegion *vreg) {
  int j, l;
  fprintf(stderr, "Confirmed reads=\n");
  for (j=0; j<vreg->nr; j++) {
      if (vreg->reads[j].allele_id >= vreg->nca)
        continue;

      for (l=0; l< vreg->end-vreg->beg; l++)
          fprintf(stderr, "%c", vreg->reads[j].bases[l]);

      fprintf(stderr, " allele= %d iid= %d\n",
              vreg->reads[j].allele_id, vreg->reads[j].iid);
    }
  fprintf(stderr, "\n");
}

static
void
show_reads_of_confirmed_alleles(VarRegion *vreg) {
  int j, l;
  fprintf(stderr, "Confirmed reads=\n");
  for (j=0; j<vreg->nca; j++)
      for (l=0; l< vreg->alleles[j].num_reads; l++)
          fprintf(stderr, "%d  %d \n", vreg->alleles[j].read_ids[l],
                  vreg->alleles[j].read_iids[l]);
  fprintf(stderr, "\n");
}


static
void
show_reads(VarRegion *vreg) {
  int j, l;

  fprintf(stderr, "Num_reads= %d vreg.beg= %d vreg.end= %d\n",
          vreg->nr, vreg->beg, vreg->end);

  fprintf(stderr, "Reads=\n");

  for (j=0; j<vreg->nr; j++)
      for (l=0; l< vreg->end-vreg->beg; l++)
          fprintf(stderr, "%c", vreg->reads[j].bases[l]);
      fprintf(stderr, "\n");
    }

  fprintf(stderr, "Ave_qvs= \n");

  for (j=0; j<vreg->nr; j++)
    fprintf(stderr, "%f ", vreg->reads[j].ave_qv);

  fprintf(stderr, "\n");
}
#endif







static
int
GetTheMostDistantRead(int curr_read_id, int32 nr, int32 **dist_matrix) {
  int i, dist_read_id = -1;
  int max_dist = -1;

  for (i=0; i<nr; i++) {
      if (i == curr_read_id)
        continue;

      if (max_dist < dist_matrix[curr_read_id][i]) {
          max_dist = dist_matrix[curr_read_id][i];
          dist_read_id = i;
        }
    }

  return dist_read_id;
}







static
void
PopulateVARRecord(int is_phased, int32 *cids, int32 *nvars, int32 *min_len_vlist,
                  IntMultiVar **v_list, VarRegion vreg, CNS_Options *opp, int get_scores,
                  int32 *conf_read_iids) {
  double fict_var;
  int   m;
  int   num_reported_alleles= (vreg.nca < 2) ? 2 : vreg.nca;
  char  cbase;
  char *base = (char*)safe_calloc(num_reported_alleles,sizeof(char));
  char  buf[10000];
  int   tot_num_conf_reads = 0;

  for (m=0; m<vreg.nca; m++)
    tot_num_conf_reads += vreg.alleles[m].num_reads;

  if (!(*v_list)) {
    *v_list = (IntMultiVar *)safe_malloc(*min_len_vlist*
                                         sizeof(IntMultiVar));
  }
  if (*nvars == *min_len_vlist) {
    *min_len_vlist += 10;
    *v_list = (IntMultiVar *)safe_realloc(*v_list, *min_len_vlist*
                                          sizeof(IntMultiVar));
  }
  vreg_id++;
  (*v_list)[*nvars].position.bgn = vreg.beg;
  (*v_list)[*nvars].position.end = vreg.end;
  (*v_list)[*nvars].num_reads = (int32)vreg.nr;
  (*v_list)[*nvars].num_conf_alleles = vreg.nca;
  (*v_list)[*nvars].min_anchor_size = opp->smooth_win;
  (*v_list)[*nvars].var_length = vreg.end-vreg.beg;
  (*v_list)[*nvars].curr_var_id = vreg_id;
  (*v_list)[*nvars].phased_var_id = is_phased ? vreg_id-1 : -1;
  (*v_list)[*nvars].weights = (char*)safe_calloc(num_reported_alleles,
                                                 (5+2)* sizeof(char));
  (*v_list)[*nvars].nr_conf_alleles = (char*)safe_calloc(num_reported_alleles,
                                                         (2+2)* sizeof(char));
  (*v_list)[*nvars].var_seq = (char*)safe_malloc(num_reported_alleles*
                                                 (vreg.end-vreg.beg+1)* sizeof(char));
  (*v_list)[*nvars].conf_read_iids=
    (char*)safe_calloc(vreg.nr, (15+2)*sizeof(char));
  NumVARRecords++;
  {
    int al, rd;
    int32 shift   = vreg.end-vreg.beg+1;
    int distant_read_id = -42, distant_allele_id = -42;
    if (vreg.nca < 2)
      {
        distant_read_id = GetTheMostDistantRead(vreg.alleles[0].read_ids[0],
                                                vreg.nr, vreg.dist_matrix);
        distant_allele_id = vreg.reads[distant_read_id].allele_id;
      }
#ifdef DEBUG_VAR_RECORDS
    fprintf(stderr, "VAR beg= %d end= %d\n", vreg.beg, vreg.end);
    OutputReads(stderr, vreg.reads, vreg.nr, vreg.end-vreg.beg);
    OutputDistMatrix(stderr, &vreg);
    OutputAlleles(stderr, &vreg);
#endif

    for (m=0; m<vreg.end-vreg.beg; m++)
      {
        for (al=num_reported_alleles-1; al >=0; al--)
          {
            if (al == 0 || al < vreg.nca)
              {
                int read_id = vreg.alleles[al].read_ids[0];
                base[al] = vreg.reads[read_id].bases[m];
                if (al == 0)
                  {
                    int32 cid = cids[vreg.beg+m];
                    Column *column=GetColumn(columnStore,cid);
                    Bead *call = GetBead(beadStore, column->call);

#if 0
                    // Check the consistency of a consensus call
                    BaseCall(cids[vreg.beg+m], 1, &fict_var, &vreg,
                             vreg.alleles[al].id, &cbase, 0, 0, opp);
                    if (cbase != base[al])
                      {
                        fprintf(stderr, "Error setting the consensus base #%d  %c/%c\n", m, cbase, base[al]);
                        OutputReads(stderr, vreg.reads, vreg.nr, vreg.end-vreg.beg);
                        OutputDistMatrix(stderr, &vreg);
                        OutputAlleles(stderr, &vreg);
                        exit(1);
                      }
#endif

                    // Set the consensus quality
                    BaseCall(cids[vreg.beg+m], 1, &fict_var, &vreg,
                             -1, &cbase, 0, 0, opp);

                    // Set the consensus base
                    Setchar(sequenceStore, call->soffset, &base[al]);

                  }
              }
            else // vreg.nca < 2 and al == 1
              {
                base[al] = vreg.reads[distant_read_id].bases[m];
              }
            (*v_list)[*nvars].var_seq[m+al*shift] = base[al];
          }
        if (get_scores > 0)
          UpdateScores(vreg, base, num_reported_alleles);
      }

    sprintf((*v_list)[*nvars].weights, "");
    sprintf((*v_list)[*nvars].nr_conf_alleles, "");
    for (al=0; al < num_reported_alleles; al++)
      {
        int    weight = vreg.alleles[al].weight;
        int    num_reads = vreg.alleles[al].num_reads;
        const char *format_weight = (al < num_reported_alleles-1) ? "%d/" : "%d\0";
        const char *format_num_reads = (al < num_reported_alleles-1) ? "%d/" : "%d\0";

        (*v_list)[*nvars].var_seq[-1+(al+1)*shift] =
          (al < num_reported_alleles-1) ? '/' : '\0';

        sprintf(buf, format_weight, ROUND(weight));
        (*v_list)[*nvars].weights = strcat((*v_list)[*nvars].weights, buf);

        sprintf(buf, format_num_reads, num_reads);
        (*v_list)[*nvars].nr_conf_alleles = strcat((*v_list)[*nvars].nr_conf_alleles, buf);
      }
#if 0
    fprintf(stderr, "len= %d var_seq = %s\n", vreg.end-vreg.beg, (*v_list)[*nvars].var_seq);

#endif
    sprintf((*v_list)[*nvars].conf_read_iids, "");
    for (rd=0; rd < tot_num_conf_reads; rd++)
      {
        const char *format_iids = (rd < tot_num_conf_reads-1) ? "%d/" : "%d\0";
        sprintf(buf, format_iids, conf_read_iids[rd]);
        (*v_list)[*nvars].conf_read_iids = strcat((*v_list)[*nvars].conf_read_iids, buf);
      }

    for (al=0; al < num_reported_alleles; al++)
      {
        if ((*v_list)[*nvars].var_seq[al*shift]             == '-' &&
            (*v_list)[*nvars].var_seq[al*shift + shift - 2] == '-')
          {
            NumVARStringsWithFlankingGaps++;
          }
      }
  }
  safe_free(base);
#ifdef DEBUG_VAR_RECORDS
  fprintf(stderr, "VARiation= %s\n", (*v_list)[*nvars].var_seq);
#endif
  (*nvars)++;
}




static
void
SetConsensusToMajorAllele(int32 *cids, VarRegion vreg, CNS_Options *opp, int get_scores,
                          int32 *conf_read_iids) {
  int i,m;
  char  cbase;
  int read_id = vreg.alleles[0].read_ids[0];
#ifdef DEBUG_VAR_RECORDS
  char *bases = (char*)safe_calloc((vreg.end-vreg.beg+1),sizeof(char));
  bases[vreg.end-vreg.beg]='\0';
  fprintf(stderr, "VAR beg= %d end= %d\n", vreg.beg, vreg.end);
  OutputReads(stderr, vreg.reads, vreg.nr, vreg.end-vreg.beg);
  OutputDistMatrix(stderr, &vreg);
  OutputAlleles(stderr, &vreg);
#endif
  for (m=0; m<vreg.end-vreg.beg; m++)
    {
      int32   cid = cids[vreg.beg+m];
      Column *column=GetColumn(columnStore,cid);
      Bead   *call = GetBead(beadStore, column->call);

      // Set the consensus base
      cbase = vreg.reads[read_id].bases[m];
      Setchar(sequenceStore, call->soffset, &cbase);
#ifdef DEBUG_VAR_RECORDS
      bases[m]=cbase;
    }
  fprintf(stderr, "beg= %d end= %d Num_reads= %d \nConsensus= \n%s\n", vreg.beg, vreg.end, vreg.nr, bases);
  fprintf(stderr, "Reads=\n");
  for (read_id=0; read_id<vreg.nr; read_id++)
    {
      for (m=0; m<vreg.end-vreg.beg; m++)
        {
          fprintf(stderr, "%c", vreg.reads[read_id].bases[m]);
        }
      fprintf(stderr, "   allele=%d \n", vreg.reads[read_id].allele_id);
    }
  safe_free(bases);
#else
  }
#endif
}





int
RefreshMANode(int32 mid, int quality, CNS_Options *opp, int32 *nvars,
              IntMultiVar **v_list, int make_v_list, int get_scores) {
  // refresh columns from cid to end
  // if quality == -1, don't recall the consensus base
  int     i=0, i1, j=0, l=0, index=0, len_manode = MIN_SIZE_OF_MANODE;
  int32   cid=0, *cids=NULL, *prev_iids=NULL;
  int     window=0, svbeg=0, svend=0, max_prev_nr=INITIAL_NR, prev_nr=0;
  char    cbase=0, abase=0, *prev_bases = NULL;
  char   *var_seq=NULL;
  double *varf=NULL, *svarf=NULL;
  Column *column = NULL;
  VarRegion  vreg = {0};
  char  **reads = NULL;
  MANode *ma = GetMANode(manodeStore,mid);
  int32   min_len_vlist = 10;

  // Variables used to phase VAR records
  int32  prev_nca          = 0;     // valid size of array prev_nca_iid
  int32  prev_ncr          = 0;     // valid size of array prev_ncr_iid
  int32  prev_nca_iid_max  = 4096;  // allocated size of arrays; err on the large side
  int32  prev_ncr_iid_max  = 4096;  //   hopefully avoiding reallocs.
  int32 *prev_nca_iid      = NULL;  // number of reads in 10 first confirmed alleles
  int32 *prev_ncr_iid      = NULL;  // iids of the first 100 reads, rev. sorted by allele

  //  Make sure that we have valid options here, we then reset the
  //  pointer to the freshly copied options, so that we can always
  //  assume opp is a valid pointer
  //
  CNS_Options  opp_private;
  if (opp == NULL) {
    opp_private.split_alleles   = CNS_OPTIONS_SPLIT_ALLELES_DEFAULT;
    opp_private.smooth_win      = CNS_OPTIONS_MIN_ANCHOR_DEFAULT;
    opp = &opp_private;
  }

  window = opp->smooth_win;

  if (nvars)
    *nvars = 0;

  assert(ma != NULL);

  if ( ma->first == -1 )
    return 1;

  vreg.nr     = 0;
  vreg.max_nr = MIN_ALLOCATED_DEPTH;
  vreg.iids  = (int32 *)safe_calloc(vreg.max_nr, sizeof(int32));
  vreg.curr_bases =  (char *)safe_calloc(vreg.max_nr, sizeof(char));

  varf     = (double *)safe_calloc(len_manode, sizeof(double));
  cids     =  (int32 *)safe_calloc(len_manode, sizeof(int32));
  Resetint32(ma->columns);
  cid = ma->first;
  vreg.nr = -1;

  if (get_scores > 0) {
    prev_bases = (char  *)safe_malloc(max_prev_nr*sizeof(char ));
    prev_iids  = (int32 *)safe_malloc(max_prev_nr*sizeof(int32));
  }

  // Calculate variation as a function of position in MANode.
  while ( cid  > -1 )
    {
      column = GetColumn(columnStore, cid);
      assert(column != NULL);
      if ( quality != -2 )
        {
          if (index >= len_manode)
            {
              len_manode += MIN_SIZE_OF_MANODE;
              varf  = (double *)safe_realloc(varf,  len_manode*sizeof(double));
              cids = (int32 *)safe_realloc(cids, len_manode*sizeof(int32));
            }
          // Call consensus using all alleles
          // The goal is to detect a variation at a given position
          BaseCall(cid, quality, &(varf[index]), &vreg, -1, &cbase, 0, get_scores,
                   opp);
          cids[index] = cid;
        }
      column->ma_index = index;
      AppendVA_int32(ma->columns, &cid);
      // sanity check
      if (index>0) {
        int32 prev= *Getint32(ma->columns, index-1);
        Column *pcol= GetColumn(columnStore, prev);
        if( prev != column->prev ||  pcol->next != column->lid)
          {
            fprintf(stderr, "RefreshMANode column relationships violated");
            assert(0);
          }
      }

      if (get_scores> 0)
        {
#if 0
          fprintf(stderr, "vreg.nb=%d vreg.curr_bases=", vreg.nb);
          for (i=0; i<vreg.nb; i++)
            fprintf(stderr, "%c", vreg.curr_bases[i]);
          fprintf(stderr, " prev_nr=%d prev_bases=", prev_nr);
          for (i=0; i<prev_nr; i++)
            fprintf(stderr, "%c", prev_bases[i]);
          fprintf(stderr, " NumRunsOfGaps=%d \nvreg.iids= ", NumRunsOfGaps);
          for (i=0; i<vreg.nb; i++)
            fprintf(stderr, "%d ", vreg.iids[i]);
          fprintf(stderr, "\n");
          fprintf(stderr, "prev_iids= ");
          for (i=0; i<prev_nr; i++)
            fprintf(stderr, "%d ", prev_iids[i]);
          fprintf(stderr, "\n");
#endif
          UpdateScoreNumRunsOfGaps(vreg, prev_nr, prev_bases, prev_iids,
                                   get_scores);
          UpdateScoreNumGaps(cbase, get_scores);
          if (vreg.nb > max_prev_nr) {
            max_prev_nr =  vreg.nb;
            prev_bases = (char *)safe_realloc(prev_bases,
                                              max_prev_nr*sizeof(char));
            prev_iids  = (int32 *)safe_realloc(prev_iids,
                                               max_prev_nr*sizeof(int32));
          }
          prev_nr = vreg.nb;
          for (i=0; i<vreg.nb; i++) {
            prev_bases[i] = vreg.curr_bases[i];
            prev_iids[i]  = vreg.iids[i];
          }
        }

      cid = column->next;
      index++;
    }

  if (get_scores == 1) {
    NumColumnsInUnitigs += index;
  }
  else if (get_scores == 2) {
    NumColumnsInContigs += index;
  }

  if (opp->split_alleles == 0 || quality <= 0 || make_v_list == 0)
    {
      safe_free(vreg.curr_bases);
      safe_free(vreg.iids);
      safe_free(varf);
      safe_free(cids);
      if (get_scores > 0) {
        safe_free(prev_bases);
        safe_free(prev_iids);
      }
      return 1;
    }

  assert(make_v_list == 1 || nvars  != NULL);
  assert(make_v_list == 1 || v_list != NULL);

  // Proceed further only if accurate base calls are needed
  // Smoothen variation
  len_manode = index -1;
  svarf= (double *)safe_calloc(len_manode, sizeof(double));
  for (i=0; i<len_manode; i++) {
    svarf[i] = varf[i];
    if (varf[i] < ZERO_MINUS)
      varf[i] = (varf[i] < -1.) ? 0. : -varf[i];
#if 0
    if (varf[i] > 0)
      fprintf(stderr, "i= %d varf= %f\n", i, varf[i]);
#endif
  }
  SmoothenVariation(svarf, len_manode, window);

  prev_nca_iid = (int32 *)safe_calloc(prev_nca_iid_max, sizeof(int32));
  prev_ncr_iid = (int32 *)safe_calloc(prev_ncr_iid_max, sizeof(int32));

  for (i=0; i<len_manode; i++)
    {
      if (svarf[i] == 0) {
        continue;
      }
      else
        {
          // Process a region of variation
          //
          // opp->smooth_win <  0 means no columns with variation will be grouped
          //                      into regions of variation
          // opp->smooth_win == 0 means only immediately adjacent columns with
          //                      variation will be grouped into regions of variation
          int32  is_phased = 0;
          int32 *conf_read_iids = NULL;
          double fict_var;
          int32 *allele_map;

          vreg.beg = vreg.end = svbeg = svend = i;

          /* Set the beginning of unsmoothed VAR region to the 1st position with varf != 0 */
          if (opp->smooth_win > 0)
            {
              while (vreg.beg < len_manode - 1 &&
                     DBL_EQ_DBL(varf[vreg.beg], (double)0.0))
                vreg.beg++;
            }

          /* Set the end of smoothed VAR region to the 1st position with svarf == 0 */
          if (opp->smooth_win > 0)
            {
              while ((svend < len_manode-1) && (svarf[svend] > ZERO_PLUS))
                svend++;
            }
          else if (svend < len_manode-1)
            svend++;

          /* Set the end of unsmoothed VAR region to the 1st position with varf == 0 */
          vreg.end = vreg.beg;
          if (opp->smooth_win > 0)
            {
              while (vreg.end < len_manode-1 && varf[vreg.end] > ZERO_PLUS)
                vreg.end++;
            }
          else if (svend < len_manode-1)
            vreg.end++;

          // Store iids of all the reads in current region
          vreg.nr = 0;
          for(l=0; l<vreg.max_nr; l++)
            vreg.iids[l] = -1;

          // Get all the read iids
          // Calculate the total number of reads, vreg.nr (corresponding to any allele)
          for (j=vreg.beg; j<vreg.end; j++)
            GetReadIidsAndNumReads(cids[j], &vreg);

          //  BPW:  Hmmm.  How do we have no reads here?
          if (vreg.nr == 0) {
            fprintf(stderr, "no reads for vreg beg=%d end=%d SKIPPING IT (bug?)\n", vreg.beg, vreg.end);
            continue;
          }

          // Allocate memrory for reads
          AllocateMemoryForReads(&vreg.reads, vreg.nr, vreg.end - vreg.beg,
                                 0);

          GetReadsForVARRecord(vreg.reads, vreg.iids, vreg.nr, vreg.beg,
                               vreg.end, cids);
          // Calculate a sum of qvs for each read within a variation region
          // Populate the distance matrix

#ifdef DEBUG_VAR_RECORDS
          show_reads(&vreg);
#endif

          AllocateDistMatrix(&vreg, -1);
          PopulateDistMatrix(vreg.reads, vreg.end-vreg.beg, &vreg);

#ifdef DEBUG_VAR_RECORDS
          OutputDistMatrix(stderr, &vreg);
#endif

          // Allocate memory for alleles
          AllocateMemoryForAlleles(&vreg.alleles, vreg.nr, &vreg.na);

          // Populate vreg.alleles array
          // Determine the best allele and the number of reads in this allele
          ClusterReads(vreg.reads, vreg.nr, vreg.alleles, &(vreg.na),
                       &(vreg.nca), vreg.dist_matrix);
#if 0
          if (vreg.nca < 2) {
            fprintf(stderr, "nca= %d\n", vreg.nca);
            show_reads(&vreg);
          }
#endif
          if (make_v_list == 2)
            is_phased = PhaseWithPrevVreg(vreg.nca, vreg.alleles, vreg.reads,
                                          &allele_map,
                                          prev_nca, prev_nca_iid, prev_nca_iid_max,
                                          prev_ncr, prev_ncr_iid, prev_ncr_iid_max);

          if (is_phased)
            {
              SortAllelesByMapping(vreg.alleles, vreg.nca, vreg.reads, allele_map);
              safe_free(allele_map);
            }
          else
            SortAllelesByWeight(vreg.alleles, vreg.na, vreg.reads);

          // Create a list of iids of confirmed reads
          {
            int start_num = 0, num_conf_reads = 0;

            for (i1=0; i1<vreg.nca; i1++)
              num_conf_reads += vreg.alleles[i1].num_reads;

            conf_read_iids = (int32 *)safe_calloc(num_conf_reads, sizeof(int32));

            for (i1=0; i1<vreg.nca; i1++)
              {
                for (j=0; j<vreg.alleles[i1].num_reads; j++)
                  {
                    conf_read_iids[start_num+j] = vreg.alleles[i1].read_iids[j];
                  }
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
            int iv, jv, kv = 0;

            //  we should have space enough for these, right?  If
            //  not, count and reallocate new space -- be sure to
            //  clear it!

            prev_nca = vreg.nca;
            prev_ncr = 0;
            for (iv = 0; iv < vreg.nca; iv++)
              prev_ncr += vreg.alleles[iv].num_reads;

            if (prev_nca_iid_max < prev_nca) {
              while (prev_nca_iid_max < prev_nca)
                prev_nca_iid_max *= 2;
              safe_free(prev_nca_iid);
              prev_nca_iid = (int32 *)safe_calloc(prev_nca_iid_max, sizeof(int32));
            }
            if (prev_ncr_iid_max < prev_ncr) {
              while (prev_ncr_iid_max < prev_ncr)
                prev_ncr_iid_max *= 2;
              safe_free(prev_ncr_iid);
              prev_ncr_iid = (int32 *)safe_calloc(prev_ncr_iid_max, sizeof(int32));
            }

            memset(prev_nca_iid, 0, sizeof(int32) * prev_nca_iid_max);
            memset(prev_ncr_iid, 0, sizeof(int32) * prev_ncr_iid_max);

            for (iv = 0; iv < vreg.nca; iv++) {
              prev_nca_iid[iv] = vreg.alleles[iv].num_reads;
              for (jv = 0; jv < vreg.alleles[iv].num_reads; jv++)
                prev_ncr_iid[kv++] = vreg.alleles[iv].read_iids[jv];
            }

#ifdef DEBUG_VAR_RECORDS
            OutputAlleles(stderr, &vreg);
#endif

            /* Store variations in a v_list */
            if (make_v_list == 2)
              PopulateVARRecord(is_phased, cids, nvars, &min_len_vlist, v_list,
                                vreg, opp, get_scores, conf_read_iids);
            else
              SetConsensusToMajorAllele(cids, vreg, opp, get_scores,
                                        conf_read_iids);
          }

          if (opp->smooth_win > 0)
            i = svend;

          for (j=0; j<vreg.nr; j++)
            {
              safe_free(vreg.dist_matrix[j]);
              safe_free(vreg.reads[j].bases);
              safe_free(vreg.reads[j].qvs);
              safe_free(vreg.alleles[j].read_ids);
              safe_free(vreg.alleles[j].read_iids);
            }
          safe_free(vreg.dist_matrix);
          safe_free(vreg.reads);
          safe_free(vreg.alleles);
          vreg.nr = 0;
          safe_free(conf_read_iids);
        }
    }
  safe_free(vreg.curr_bases);
  safe_free(vreg.iids);
  safe_free(varf);
  safe_free(svarf);
  safe_free(cids);
  if (get_scores > 0) {
    safe_free(prev_bases);
    safe_free(prev_iids);
  }
  safe_free(prev_nca_iid);
  safe_free(prev_ncr_iid);
  return 1;
}


