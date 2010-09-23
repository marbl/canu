
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

static const char *rcsid = "$Id: AS_BOG_ArrivalRate.cc,v 1.1 2010-09-23 09:34:50 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_BestOverlapGraph.hh"

#include "MultiAlignStore.h"

#undef max



float UnitigGraph::getGlobalArrivalRate(long total_random_frags_in_genome, long genome_size){

  float _globalArrivalRate;

  //  If the genome size has not been specified, estimate the GAR.
  if(genome_size == 0){

    float total_rho=0, avg_rho;
    float total_arrival_frags=0;
    size_t rho_gt_10000 = 0;

    // Go through all the unitigs to sum rho and unitig arrival frags
    UnitigVector::const_iterator iter;
    for(
        iter=unitigs->begin();
        iter!=unitigs->end();
        iter++){

      if (*iter == NULL)
        continue;

      avg_rho = (*iter)->getAvgRho(_fi);
      total_rho += avg_rho;
      if (avg_rho > 10000.0)
        rho_gt_10000 += (size_t)avg_rho / 10000;

      float unitig_random_frags = (*iter)->getNumRandomFrags();
      if (--unitig_random_frags < 0)
        unitig_random_frags = 0;

      total_arrival_frags += unitig_random_frags;
      (*iter)->setLocalArrivalRate( unitig_random_frags / avg_rho );
    }
    // Estimate GAR
    _globalArrivalRate = (total_rho > 0) ? (total_arrival_frags / total_rho): 0;

    std::cerr << "Calculated Global Arrival rate " << _globalArrivalRate <<
      std::endl;
    // Now recalculate based on big unitigs, copied from AS_CGB/AS_CGB_cgb.c
    if (rho_gt_10000 * 20000 > total_rho) {
      float min_10_local_arrival_rate          = _globalArrivalRate;
      float median_local_arrival_rate          = _globalArrivalRate;
      float max_local_arrival_rate             = _globalArrivalRate;
      float recalibrated_fragment_arrival_rate = _globalArrivalRate;
      size_t num_arrival_rates=0;
      int median_index;
      std::vector<float> arrival_rate_array(rho_gt_10000);

      for( iter=unitigs->begin(); iter!=unitigs->end(); iter++) {
        if (*iter == NULL)
          continue;
        avg_rho = (*iter)->getAvgRho(_fi);
        if (avg_rho > 10000.0) {
          const int num_10000 = (size_t)avg_rho / 10000;
          const float local_arrival_rate =
            (*iter)->getNumRandomFrags() / avg_rho;
          assert(num_10000 > 0);
          int i;
          for(i=0;i<num_10000;i++){
            assert(i < rho_gt_10000);
            arrival_rate_array[i] = local_arrival_rate;
            num_arrival_rates++;
          }
          if(num_arrival_rates > 0){
            float tmp_fragment_arrival_rate, max_diff_arrival_rate;
            float prev_arrival_rate, cur_arrival_rate, diff_arrival_rate;
            int max_diff_index;
            std::sort(arrival_rate_array.begin(),arrival_rate_array.end());
            min_10_local_arrival_rate = arrival_rate_array[num_arrival_rates / 10];
            median_index = (num_arrival_rates * 5) / 10;
            median_local_arrival_rate = arrival_rate_array[median_index];
            max_local_arrival_rate = arrival_rate_array[num_arrival_rates-1];
            recalibrated_fragment_arrival_rate =
              arrival_rate_array[(num_arrival_rates * 19) / 20];
            prev_arrival_rate = min_10_local_arrival_rate;
            max_diff_arrival_rate = 0.0;
            for(i=num_arrival_rates / 10;i<median_index;i++){
              cur_arrival_rate = arrival_rate_array[i];
              diff_arrival_rate = cur_arrival_rate - prev_arrival_rate;
              prev_arrival_rate = cur_arrival_rate;
              if(diff_arrival_rate > max_diff_arrival_rate){
                max_diff_arrival_rate = diff_arrival_rate;
              }
            }
            max_diff_arrival_rate *= 2.0;
            max_diff_index = num_arrival_rates - 1;
            for(i=median_index;i<num_arrival_rates;i++){
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
            tmp_fragment_arrival_rate =MIN(min_10_local_arrival_rate * 2.0,
                                           median_local_arrival_rate * 1.25);
            if(tmp_fragment_arrival_rate < recalibrated_fragment_arrival_rate){
              recalibrated_fragment_arrival_rate = tmp_fragment_arrival_rate;
            }
            if(max_diff_arrival_rate < recalibrated_fragment_arrival_rate){
              recalibrated_fragment_arrival_rate = max_diff_arrival_rate;
            }
          }
          if(recalibrated_fragment_arrival_rate > _globalArrivalRate){
            _globalArrivalRate = recalibrated_fragment_arrival_rate;
#if 0
            std::cerr <<
              "Used recalibrated global_fragment_arrival_rate="
                      << _globalArrivalRate << std::endl
                      << "Used recalibrated global_fragment_arrival_distance="
                      <<((_globalArrivalRate > 0.) ? 1./(_globalArrivalRate) : 0.)
                      << std::endl
                      << "Chunk arrival rates sorted at 1/100s"
                      << std::endl ;
            for(i=0;i<100;i++) {
              std::cerr<<arrival_rate_array[((num_arrival_rates * i) / 100)]
                       << std::endl;
            }
            std::cerr << max_local_arrival_rate << std::endl;
#endif
          }
        }
      }
    }

    std::cerr << 
      "Computed genome_size="
              << (_globalArrivalRate > 0.f ? ((float)total_random_frags_in_genome / _globalArrivalRate) : 0.f) 
              << std::endl;
  }else{
    // Compute actual GAR
    _globalArrivalRate = (float)total_random_frags_in_genome / (float)genome_size;
  }

  return(_globalArrivalRate);

}

