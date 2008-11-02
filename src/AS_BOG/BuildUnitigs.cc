
/**************************************************************************
 * Copyright (C) 2005, J Craig Venter Institute. All rights reserved.
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

const char *mainid = "$Id: BuildUnitigs.cc,v 1.54 2008-11-02 06:27:53 brianwalenz Exp $";

#include<vector>
#include<cmath>
#include<cstdlib>

#include "AS_BOG_BestOverlapGraph.hh"
#include "AS_BOG_ChunkGraph.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_MateChecker.hh"

using std::vector;

extern "C" {
#include "getopt.h"
#include "AS_OVS_overlapStore.h"
#include "AS_CGB_histo.h"
}


void outputHistograms(UnitigGraph *utg, FragmentInfo *fi, FILE *stats) {
  const int nsample=500;
  const int nbucket=500;
  MyHistoDataType zork;

  Histogram_t *length_of_unitigs_histogram = create_histogram(nsample,nbucket,0,TRUE);
  Histogram_t *covg_histogram = create_histogram(nsample,nbucket,0,TRUE);
  Histogram_t *arate_histogram = create_histogram(nsample,nbucket,0,TRUE);

  extend_histogram(length_of_unitigs_histogram, sizeof(MyHistoDataType),
                   myindexdata,mysetdata,myaggregate,myprintdata);
  extend_histogram(covg_histogram, sizeof(MyHistoDataType),
                   myindexdata,mysetdata,myaggregate,myprintdata);
  extend_histogram(arate_histogram, sizeof(MyHistoDataType),
                   myindexdata,mysetdata,myaggregate,myprintdata);

  UnitigVector::const_iterator uiter = utg->unitigs->begin();
  for(;uiter != utg->unitigs->end(); uiter++) {

    Unitig *u = *uiter;
    if (u == NULL)
      continue;
    zork.nsamples = 1;
    int numFrags = u->getNumFrags();
    zork.sum_frags = zork.min_frags = zork.max_frags = numFrags;

    int bases = u->getLength();
    zork.sum_bp = zork.min_bp = zork.max_bp = bases;

    int rho = static_cast<int>(u->getAvgRho(fi));
    zork.sum_rho = zork.min_rho = zork.max_rho = rho;

    //long covg = static_cast<long>(rint(u->getCovStat() * 1000));
    int covg = static_cast<int>(rintf(u->getCovStat(fi)));
    zork.sum_discr = zork.min_discr = zork.max_discr = covg;

    float arateF = u->getLocalArrivalRate(fi) * 10000;
    int arate = static_cast<int>(rintf(arateF));
    if (arate < 0)
      fprintf(stats, "Negative Local ArrivalRate %f id %d arate %d\n", arateF, u->id(), arate);
    zork.sum_arrival = zork.min_arrival = zork.max_arrival = arate;

    zork.sum_rs_frags=zork.min_rs_frags=zork.max_rs_frags=0;
    zork.sum_nr_frags=zork.min_nr_frags=zork.max_nr_frags=0;

    add_to_histogram(length_of_unitigs_histogram, u->getLength(), &zork);
    add_to_histogram(covg_histogram, covg, &zork);
    add_to_histogram(arate_histogram, arate, &zork);
  }

  fprintf(stats, "Length of Unitigs histogram\n");
  print_histogram(stats,length_of_unitigs_histogram, 0, 1);

  fprintf(stats, "Unitig Coverage Stat histogram\n");
  print_histogram(stats,covg_histogram, 0, 1);

  fprintf(stats, "Unitig Arrival Rate histogram\n");
  print_histogram(stats,arate_histogram, 0, 1);

  free_histogram( length_of_unitigs_histogram );
  free_histogram( covg_histogram );
  free_histogram( arate_histogram );
}


int
main (int argc, char * argv []) {

  // Get path/names of olap and frg stores from command line
  const char* OVL_Store_Path;
  const char* GKP_Store_Path;

  double        erate = .015;

  long genome_size=0;
  int ch;
  bool argsDone=false;

  int   fragment_count_target = 0;
  char *output_prefix         = "bog";

  bool  unitigIntersectBreaking = false;
  bool  ejectUnhappyContained   = false;
  int   badMateBreakThreshold   = -7;

  argc = AS_configure(argc, argv);

  optarg = NULL;
  while(!argsDone && (ch = getopt(argc, argv,"B:o:O:G:e:m:s:bk"))) {
    switch(ch) {
      case -1:
        argsDone=true;
        break;

      case 'B':
        fragment_count_target = atoi(optarg);
        break;
      case 'o':
        output_prefix = optarg;
        break;

      case 'G':
        GKP_Store_Path = strdup(optarg);
        assert( GKP_Store_Path != NULL ); break;
      case 'O':
        OVL_Store_Path = strdup(optarg);
        assert( OVL_Store_Path != NULL ); break;
      case 'b':
        unitigIntersectBreaking = true; break;
      case 'e':
        {
          erate = atof(optarg);
          if ((erate < 0.0) || (AS_MAX_ERROR_RATE < erate))
            fprintf(stderr, "Invalid overlap error threshold %s; must be between 0.00 and %.2f.\n",
                    optarg, AS_MAX_ERROR_RATE), exit(1);

          fprintf(stderr, "The overlap error threshold = %.3f = %.3f%%\n",
                  erate, AS_OVS_decodeQuality(erate) * 100.0);
        }
        break;
      case 'k':
        ejectUnhappyContained = true; break;
      case 'm':
        badMateBreakThreshold = -atoi(optarg); break;
      case 's':
        genome_size = atol(optarg); break;
      default:
        fprintf(stderr,"Unrecognized option -%c optarg %s\n\n",optopt, optarg);
        fprintf(stderr, "usage: %s -O <OVL Store Path> -G <GKP Store Path>\n", argv[0]);
        fprintf(stderr, "\n");
        fprintf(stderr, "[-B b]       Target number of fragments per IUM batch.\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "[-o prefix]  Output prefix name\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "[-s <genome size>]\n");
        fprintf(stderr, "  If the genome size is set to 0, this will cause the unitigger\n");
        fprintf(stderr, "  to try to estimate the genome size based on the constructed\n");
        fprintf(stderr, "  unitig lengths.\n");
        fprintf(stderr, "[-b] Break promisciuous unitigs at unitig intersection points\n");
        fprintf(stderr, "[-e] Fraction error to generate unitigs for; default is 0.015\n");
        fprintf(stderr, "[-k] Kick out unhappy contained mated reads into singleton unitigs\n");
        fprintf(stderr, "[-m] Number of bad mates in a region required to break a unitig\n");
        fprintf(stderr, "     default is 7\n");
        fprintf(stderr, " \n");
        exit(1);
    }
  }

  fprintf(stderr, "Genome Size: "F_S64"\n", genome_size);

  GateKeeperStore  *gkpStore = openGateKeeperStore(GKP_Store_Path, FALSE);
  OverlapStore     *ovlStore = AS_OVS_openOverlapStore(OVL_Store_Path);

  FragmentInfo     *fragInfo = new FragmentInfo(gkpStore);

  BestOverlapGraph      *BOG = new BestOverlapGraph(fragInfo, ovlStore, erate);

  ChunkGraph *cg = new ChunkGraph(fragInfo, BOG);
  UnitigGraph utg(fragInfo, BOG);
  utg.build(cg, unitigIntersectBreaking, output_prefix);

  MateChecker  mateChecker(fragInfo);
  mateChecker.checkUnitigGraph(utg, badMateBreakThreshold);

  float globalARate = utg.getGlobalArrivalRate(getNumGateKeeperRandomFragments(gkpStore), genome_size);
  Unitig::setGlobalArrivalRate(globalARate);

  utg.writeIUMtoFile(output_prefix, fragment_count_target);

  {
    char  filename[FILENAME_MAX] = {0};
    sprintf(filename, "%s.cga.0", output_prefix);

    FILE *stats = fopen(filename,"w");
    assert(NULL != stats);

    fprintf(stats, "Global Arrival Rate: %f\n", globalARate);
    fprintf(stats, "There were %d unitigs generated.\n", utg.unitigs->size());

    outputHistograms( &utg, fragInfo, stats );

    fclose(stats);
  }

  delete BOG;
  delete cg;
  delete fragInfo;

  AS_OVS_closeOverlapStore(ovlStore);
  closeGateKeeperStore(gkpStore);

  fprintf(stderr, "Bye.\n");

  return 0;
}


