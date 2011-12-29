
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

static const char *rcsid = "$Id: AS_BOG_ArrivalRate.cc,v 1.5 2011-12-29 09:26:03 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_BestOverlapGraph.hh"

#include "MultiAlignStore.h"




double UnitigGraph::getGlobalArrivalRate(long total_random_frags_in_genome, long genome_size){

  if (genome_size != 0)
    return((double)total_random_frags_in_genome / (double)genome_size);

  double globalArrivalRate;
  double total_rho=0, avg_rho;
  double total_arrival_frags=0;
  size_t rho_gt_10000 = 0;

  // Go through all the unitigs to sum rho and unitig arrival frags
  UnitigVector::const_iterator iter;
  for(iter=unitigs.begin();
      iter!=unitigs.end();
      iter++){

    if (*iter == NULL)
      continue;

    avg_rho = (*iter)->getAvgRho();
    total_rho += avg_rho;
    if (avg_rho > 10000.0)
      rho_gt_10000 += (size_t)avg_rho / 10000;

    double unitig_random_frags = (*iter)->getNumRandomFrags();
    if (--unitig_random_frags < 0)
      unitig_random_frags = 0;

    total_arrival_frags += unitig_random_frags;
    (*iter)->setLocalArrivalRate( unitig_random_frags / avg_rho );
  }
  // Estimate GAR
  globalArrivalRate = (total_rho > 0) ? (total_arrival_frags / total_rho): 0;

  fprintf(logFile, "Calculated Global Arrival rate %f\n", globalArrivalRate);

  // Now recalculate based on big unitigs, copied from AS_CGB/AS_CGB_cgb.c
  if (rho_gt_10000 * 20000 > total_rho) {
    double min_10_local_arrival_rate          = globalArrivalRate;
    double median_local_arrival_rate          = globalArrivalRate;
    double max_local_arrival_rate             = globalArrivalRate;
    double recalibrated_fragment_arrival_rate = globalArrivalRate;
    size_t num_arrival_rates=0;
    int median_index;
    std::vector<double> arrival_rate_array(rho_gt_10000);

    for( iter=unitigs.begin(); iter!=unitigs.end(); iter++) {
      if (*iter == NULL)
        continue;
      avg_rho = (*iter)->getAvgRho();
      if (avg_rho > 10000.0) {
        const int num_10000 = (size_t)avg_rho / 10000;
        const double local_arrival_rate =
          (*iter)->getNumRandomFrags() / avg_rho;
        assert(num_10000 > 0);
        for(uint32 i=0;i<num_10000;i++){
          assert(i < rho_gt_10000);
          arrival_rate_array[i] = local_arrival_rate;
          num_arrival_rates++;
        }
        if(num_arrival_rates > 0){
          double tmp_fragment_arrival_rate, max_diff_arrival_rate;
          double prev_arrival_rate, cur_arrival_rate, diff_arrival_rate;
          int max_diff_index;
          std::sort(arrival_rate_array.begin(),arrival_rate_array.end());
          min_10_local_arrival_rate = arrival_rate_array[num_arrival_rates / 10];
          median_index = (num_arrival_rates * 5) / 10;
          median_local_arrival_rate = arrival_rate_array[median_index];
          max_local_arrival_rate = arrival_rate_array[num_arrival_rates-1];
          recalibrated_fragment_arrival_rate = arrival_rate_array[(num_arrival_rates * 19) / 20];
          prev_arrival_rate = min_10_local_arrival_rate;
          max_diff_arrival_rate = 0.0;
          for(uint32 i=num_arrival_rates / 10;i<median_index;i++){
            cur_arrival_rate = arrival_rate_array[i];
            diff_arrival_rate = cur_arrival_rate - prev_arrival_rate;
            prev_arrival_rate = cur_arrival_rate;
            if(diff_arrival_rate > max_diff_arrival_rate){
              max_diff_arrival_rate = diff_arrival_rate;
            }
          }
          max_diff_arrival_rate *= 2.0;
          max_diff_index = num_arrival_rates - 1;
          for(uint32 i=median_index;i<num_arrival_rates;i++){
            cur_arrival_rate = arrival_rate_array[i];
            diff_arrival_rate = cur_arrival_rate - prev_arrival_rate;
            prev_arrival_rate = cur_arrival_rate;
            if(diff_arrival_rate > max_diff_arrival_rate){
              max_diff_arrival_rate = diff_arrival_rate;
              max_diff_index = i - 1;
              break;
            }
          }
          max_diff_arrival_rate = arrival_rate_array[max_diff_index];
          tmp_fragment_arrival_rate =MIN(min_10_local_arrival_rate * 2.0, median_local_arrival_rate * 1.25);
          if(tmp_fragment_arrival_rate < recalibrated_fragment_arrival_rate){
            recalibrated_fragment_arrival_rate = tmp_fragment_arrival_rate;
          }
          if(max_diff_arrival_rate < recalibrated_fragment_arrival_rate){
            recalibrated_fragment_arrival_rate = max_diff_arrival_rate;
          }
        }
        if(recalibrated_fragment_arrival_rate > globalArrivalRate){
          globalArrivalRate = recalibrated_fragment_arrival_rate;
        }
      }
    }
  }

  fprintf(logFile, "Computed genome_size: %f\n",
          ((globalArrivalRate > 0.0) ? total_random_frags_in_genome / globalArrivalRate : 0.0));

  return(globalArrivalRate);

}

