
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

const char *mainid = "$Id: BuildUnitigs.cc,v 1.67 2010-01-26 02:27:04 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_ChunkGraph.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_BestOverlapGraph.hh"
#include "AS_BOG_MateChecker.hh"

#include "getopt.h"
#include "AS_CGB_histo.h"

//  Used in AS_BOG_Unitig.cc
FragmentInfo     *debugfi = 0L;
BestOverlapGraph *bog     = 0L;

FILE             *logFile      = NULL;
uint32            logFileOrder = 1;

//  Closes the current logFile, opens a new one called 'prefix.logFileOrder.name'.  If 'name' is
//  NULL, the logFile is reset to stderr.
void
setLogFile(char *prefix, char *name) {
  char  logFileName[FILENAME_MAX];

  if (logFile != stderr)
    fclose(logFile);

  if (name == NULL) {
    logFile = stderr;
    return;
  }

  sprintf(logFileName, "%s.%03u.%s", prefix, logFileOrder++, name);

  errno = 0;
  logFile = fopen(logFileName, "w");
  if (errno) {
    fprintf(stderr, "setLogFile()-- Failed to open logFile '%s': %s.\n", logFileName, strerror(errno));
    fprintf(stderr, "setLogFile()-- Will now log to stderr instead.\n");
    logFile = stderr;
  }
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

  for (uint32 ti=0; ti<utg->unitigs->size(); ti++) {
    Unitig  *u = (*utg->unitigs)[ti];

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
      fprintf(stats, "Negative Local ArrivalRate %f id %d arate %d\n", arateF, ti, arate);
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
  char      *gkpStorePath           = NULL;
  char      *ovlStoreUniqPath       = NULL;
  char      *ovlStoreReptPath       = NULL;
  char      *tigStorePath           = NULL;

  double    erate                   = 0.015;
  long      genome_size             = 0;

  int       fragment_count_target   = 0;
  char     *output_prefix           = NULL;

  bool      popBubbles              = false;
  bool      breakIntersections      = false;
  bool      joinUnitigs             = false;
  int       badMateBreakThreshold   = -7;

  argc = AS_configure(argc, argv);

  int err = 0;
  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-B") == 0) {
      fragment_count_target = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      output_prefix = argv[++arg];

    } else if (strcmp(argv[arg], "-G") == 0) {
      gkpStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      if      (ovlStoreUniqPath == NULL)
        ovlStoreUniqPath = argv[++arg];
      else if (ovlStoreReptPath == NULL)
        ovlStoreReptPath = argv[++arg];
      else
        err++;

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-U") == 0) {
      popBubbles = true;

    } else if (strcmp(argv[arg], "-b") == 0) {
      breakIntersections = true;

    } else if (strcmp(argv[arg], "-J") == 0) {
      joinUnitigs = true;

    } else if (strcmp(argv[arg], "-e") == 0) {
      erate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-m") == 0) {
      badMateBreakThreshold = -atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-s") == 0) {
      genome_size = atol(argv[++arg]);

    } else {
      err++;
    }

    arg++;
  }

  if ((erate < 0.0) || (AS_MAX_ERROR_RATE < erate))
    err++;
  if (output_prefix == NULL)
    err++;
  if (gkpStorePath == NULL)
    err++;
  if (ovlStoreUniqPath == NULL)
    err++;
  if (tigStorePath == NULL)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s -o outputName -O ovlStore -G gkpStore -T tigStore\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -O         Mandatory path to an ovlStore.\n");
    fprintf(stderr, "  -G         Mandatory path to a gkpStore.\n");
    fprintf(stderr, "  -T         Mandatory path to a tigStore (can exist or not).\n");
    fprintf(stderr, "  -o prefix  Mandatory name for the output files\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -B b       Target number of fragments per IUM batch.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -s size    If the genome size is set to 0, this will cause the unitigger\n");
    fprintf(stderr, "             to try to estimate the genome size based on the constructed\n");
    fprintf(stderr, "             unitig lengths.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -U         Enable EXPERIMENTAL short unitig merging (aka bubble popping).\n");
    fprintf(stderr, "  -J         Enable EXPERIMENTAL long unitig joining.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -b         Break promisciuous unitigs at unitig intersection points\n");
    fprintf(stderr, "  -e         Fraction error to generate unitigs for; default is 0.015\n");
    fprintf(stderr, "  -m num     Number of bad mates in a region required to break a unitig; default is 7\n");
    fprintf(stderr, " \n");

    if ((erate < 0.0) || (AS_MAX_ERROR_RATE < erate))
      fprintf(stderr, "Invalid overlap error threshold (-e option); must be between 0.00 and %.2f.\n",
              AS_MAX_ERROR_RATE);

    if (output_prefix == NULL)
      fprintf(stderr, "No output prefix name (-o option) supplied.\n");

    if (gkpStorePath == NULL)
      fprintf(stderr, "No gatekeeper store (-G option) supplied.\n");

    if (ovlStoreUniqPath == NULL)
      fprintf(stderr, "No overlap store (-O option) supplied.\n");

    if ((ovlStoreUniqPath != NULL) && (ovlStoreUniqPath == ovlStoreReptPath))
      fprintf(stderr, "Too many overlap stores (-O option) supplied.\n");

    if (tigStorePath == NULL)
      fprintf(stderr, "No output tigStore (-T option) supplied.\n");

    exit(1);
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "Bubble popping        = %s\n", (popBubbles) ? "on" : "off");
  fprintf(stderr, "Intersection breaking = %s\n", (breakIntersections) ? "on" : "off");
  fprintf(stderr, "Bad mate threshold    = %d\n", badMateBreakThreshold);
  fprintf(stderr, "Error threshold       = %.3f%%\n", erate * 100);
  fprintf(stderr, "Genome Size           = "F_S64"\n", genome_size);
  fprintf(stderr, "\n");

  gkStore          *gkpStore     = new gkStore(gkpStorePath, FALSE, FALSE);
  OverlapStore     *ovlStoreUniq = AS_OVS_openOverlapStore(ovlStoreUniqPath);
  OverlapStore     *ovlStoreRept = ovlStoreReptPath ? AS_OVS_openOverlapStore(ovlStoreReptPath) : NULL;

  FragmentInfo     *fragInfo = new FragmentInfo(gkpStore);

  debugfi = fragInfo;

  BestOverlapGraph      *BOG = new BestOverlapGraph(fragInfo, ovlStoreUniq, ovlStoreRept, erate);

  bog = BOG;

  ChunkGraph *cg = new ChunkGraph(fragInfo, BOG);
  UnitigGraph utg(fragInfo, BOG);
  utg.build(cg, breakIntersections, joinUnitigs, popBubbles, output_prefix);

  MateChecker  mateChecker(fragInfo);
  mateChecker.checkUnitigGraph(utg, badMateBreakThreshold);

  float globalARate = utg.getGlobalArrivalRate(gkpStore->gkStore_getNumRandomFragments(), genome_size);
  Unitig::setGlobalArrivalRate(globalARate);

  utg.setParentAndHang(cg);

  utg.writeIUMtoFile(output_prefix, tigStorePath, fragment_count_target);
  utg.writeOVLtoFile(output_prefix);

  {
    char  filename[FILENAME_MAX] = {0};
    sprintf(filename, "%s.cga.0", output_prefix);

    FILE *stats = fopen(filename,"w");
    assert(NULL != stats);

    fprintf(stats, "Global Arrival Rate: %f\n", globalARate);
    fprintf(stats, "There were %d unitigs generated.\n", (int)utg.unitigs->size());

    outputHistograms( &utg, fragInfo, stats );

    fclose(stats);
  }

  delete BOG;
  delete cg;
  delete fragInfo;

  delete gkpStore;
  AS_OVS_closeOverlapStore(ovlStoreUniq);
  AS_OVS_closeOverlapStore(ovlStoreRept);

  fprintf(stderr, "Bye.\n");

  return 0;
}


