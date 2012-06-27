
/**************************************************************************
 * Copyright (C) 2011, J Craig Venter Institute. All rights reserved.
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

const char *mainid = "$Id: classifyMates.C,v 1.29 2012-06-27 20:11:51 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_UTL_decodeRange.H"
#include "AS_OVS_overlapStore.h"
#include "AS_PER_gkpStore.h"

#include "classifyMates.H"
#include "classifyMates-globalData.H"

#include "sweatShop.H"


void
cmWorker(void *G, void *T, void *S) {
  cmGlobalData    *g = (cmGlobalData  *)G;
  cmThreadData    *t = (cmThreadData  *)T;
  cmComputation   *s = (cmComputation *)S;

  t->clear();

  if (g->testSpur(s, t))
    return;

  if (g->testChimers(s, t))
    return;

  if ((g->nodesMax > 0) && (s->result.classified == false))
    g->doSearchBFS(s, t);

  if ((g->depthMax > 0) && (s->result.classified == false))
    g->doSearchDFS(s, t);

  if ((g->pathsMax > 0) && (s->result.classified == false))
    g->doSearchRFS(s, t);

  if (g->suspiciousSearch) {
    g->doSearchSuspicious(s->fragIID, true,  s->mateIID, t, t->distInnie, t->distNormal);
    g->doSearchSuspicious(s->fragIID, false, s->mateIID, t, t->distAnti,  t->distOuttie);

    s->nInnie  = t->distInnie.size();
    s->nNormal = t->distNormal.size();
    s->nAnti   = t->distAnti.size();
    s->nOuttie = t->distOuttie.size();

    //  We're expecting no Normal/Anti oriented results, and, depending on PE vs MP,
    //  no results for one of Innie/Outtie.

    if ((g->innie == true) && (s->nOuttie + s->nNormal + s->nAnti > s->nInnie))
      //  Searching for innie-oriented MP.
      s->result.suspicious = true;

    if ((g->innie == false) && (s->nInnie + s->nNormal + s->nAnti > s->nOuttie))
      //  Searching for outtie-oriented PE.
      s->result.suspicious = true;
  }
}



void *
cmReader(void *G) {
  cmGlobalData    *g = (cmGlobalData  *)G;
  cmComputation   *s = NULL;

  for (; g->curFragIID < g->numFrags; g->curFragIID++) {
    if (g->fi[g->curFragIID].doSearch == false)
      continue;

    if (g->fi[g->curFragIID].mateIID == 0)
      continue;

    if (g->fi[g->curFragIID].mateIID < g->curFragIID)
      continue;

    s = new cmComputation(g->curFragIID, g->fi[g->curFragIID].mateIID, g->innie);

    g->curFragIID++;

    break;
  }

  return(s);
}



void
cmWriter(void *G, void *S) {
  cmGlobalData    *g = (cmGlobalData  *)G;
  cmComputation   *c = (cmComputation *)S;

  g->resultOutput->write(c->result);

  if ((c->result.suspicious) &&
      (g->resultSuspicious)) {
    fprintf(g->resultSuspicious, "IID %u/%u classified %c limited %c exhausted %c suspicious %c ni %u nn %u na %u no %u\n",
            c->fragIID,
            c->mateIID,
            c->result.classified ? 'T' : 'F',
            c->result.limited    ? 'T' : 'F',
            c->result.exhausted  ? 'T' : 'F',
            c->result.suspicious ? 'T' : 'F',
            c->nInnie,
            c->nNormal,
            c->nAnti,
            c->nOuttie);
    //  Don't have the values here, only in the thread data
    //fprintf(stderr, "ni %u values: ", );
    //fprintf(stderr, "nn %u values: ");
    //fprintf(stderr, "na %u values: ");
    //fprintf(stderr, "no %u values: ");
  }


  //  Only use 'classified' and 'exhausted' for adjusting the search time, since those are the two
  //  categorites where we found a result (either the other mate pair correctly, or a gap in
  //  coverage).
  //
  //  EDIT - don't use 'exhausted'.  It's a non-result.
  //
  if (c->result.classified == true) {
    g->runTime.addDataPoint(c->result.iteration);
  } else {
    delete c;
    return;
  }

  delete c;

  if (g->runTime.numData() < 1000)
    //  Not enough data to recompute statistics
    return;


  if ((g->runTime.numData() % 1000) != 0)
    //  Too soon to recompute statistics
    return;

  g->runTime.recompute();

  uint32  ci = 0;

  if (g->nodesMax > 0)  ci = g->nodesMax;
  if (g->depthMax > 0)  ci = g->depthMax;
  if (g->pathsMax > 0)  ci = g->pathsMax;

  //  If this stuff is normally distributed, 4 stddev will include 99.993666% of the data points.
  //  The new limit is the average of the old limit and the current performance.
  uint32  ni = (uint32)floor((ci + g->runTime.mean() + 4 * g->runTime.stddev()) / 2);

  if ((g->runTime.numData() % 10000) == 0)
    fprintf(stderr, "\nRUNTIME: %f +- %f  min/max %u/%u  RESET iteration limit to %u\n",
            g->runTime.mean(), g->runTime.stddev(), g->runTime.min(), g->runTime.max(), ni);

  if (g->nodesMax > 0)  g->nodesMax = ni;
  if (g->depthMax > 0)  g->depthMax = ni;
  if (g->pathsMax > 0)  g->pathsMax = ni;
}




int
main(int argc, char **argv) {
  char       *gkpStoreName      = NULL;
  char       *ovlStoreName      = NULL;
  char       *resultsName       = NULL;

  bool        doCache           = false;

  double      maxErrorFraction  = 0.045;

  uint32      distMin           = 0;
  uint32      distMax           = 0;
  bool        innie             = false;  //  require mates to be innie

  uint32      nodesMax          = 0;
  uint32      depthMax          = 0;
  uint32      pathsMax          = 0;

  bool        searchSuspicious  = true;

  set<AS_IID> searchLibs;
  set<AS_IID> backboneLibs;

  uint32      numThreads        = 4;

  uint64      memoryLimit       = UINT64_MAX;

  argc = AS_configure(argc, argv);

  int err = 0;
  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-e") == 0) {
      maxErrorFraction = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      resultsName = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      numThreads = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-m") == 0) {
      memoryLimit  = atoi(argv[++arg]);
      memoryLimit *= 1024;
      memoryLimit *= 1024;
      memoryLimit *= 1024;

    } else if (strcmp(argv[arg], "-sl") == 0) {
      AS_UTL_decodeRange(argv[++arg], searchLibs);
      
    } else if (strcmp(argv[arg], "-bl") == 0) {
      AS_UTL_decodeRange(argv[++arg], backboneLibs);

    } else if (strcmp(argv[arg], "-cache") == 0) {
      doCache = true;

    } else if (strcmp(argv[arg], "-min") == 0) {
      distMin = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-max") == 0) {
      distMax = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-innie") == 0) {
      innie = true;

    } else if (strcmp(argv[arg], "-outtie") == 0) {
      innie = false;

    } else if (strcmp(argv[arg], "-bfs") == 0) {
      nodesMax = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-dfs") == 0) {
      depthMax = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-rfs") == 0) {
      pathsMax = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-nosuspicious") == 0) {
      searchSuspicious = false;

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if (gkpStoreName == 0L)
    err++;
  if (ovlStoreName == 0L)
    err++;
  if (resultsName == 0L)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -G gkpStore -O ovlStore -o resultsFile ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -G gkpStore      Read fragments from here\n");
    fprintf(stderr, "  -O ovlStore      Read overlaps from here\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o results       Write results here\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -e maxError      Use overlaps with less than 'maxError' fraction error (default 0.045)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t n             Use 'n' compute threads\n");
    fprintf(stderr, "  -m m             Use at most 'm' GB memory (default: unlimited)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -sl l[-m]        Search for mates in libraries l-m\n");
    fprintf(stderr, "  -bl l[-m]        Use libraries l-m for searching\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -cache           Write loaded data to cache files, useful for restarting\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "STANDARD CONFIGURATION:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Will search for outtie-oriented PE mates, from 100 to 1500bp apart, using the\n");
    fprintf(stderr, "BFS algorithm.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "ADVANCED CONFIGURATION:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -min m           Mates must be at least m bases apart\n");
    fprintf(stderr, "  -max m           Mates must be at most m bases apart\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -innie           Mates must be innie\n");
    fprintf(stderr, "  -outtie          Mates must be outtie\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -bfs N           Use 'breadth-first search'; search at most N fragments\n");
    fprintf(stderr, "  -dfs N           Use 'depth-first search'; search to depth N overlaps\n");
    fprintf(stderr, "  -rfs N           Use 'random-path search'; search at most N paths\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -nosuspicious    Don't search for suspicious pairs.\n");
    fprintf(stderr, "\n");
    if (gkpStoreName == 0L)
      fprintf(stderr, "No gatekeeper store (-G) supplied.\n");
    if (ovlStoreName == 0L)
      fprintf(stderr, "No overlap store (-O) supplied.\n");
    if (resultsName == 0L)
      fprintf(stderr, "No results output (-o) supplied.\n");
    exit(1);
  }

  cmGlobalData  *g = new cmGlobalData(resultsName,
                                      distMin, distMax, innie,
                                      nodesMax,
                                      depthMax,
                                      pathsMax,
                                      memoryLimit,
                                      searchSuspicious);

  if (g->load(searchLibs,
              backboneLibs,
              maxErrorFraction) == false) {
    g->saveCache = doCache;  //  No better place...

    g->loadFragments(gkpStoreName, searchLibs, backboneLibs);
    g->loadOverlaps(ovlStoreName, maxErrorFraction);
    g->save();
  }

  sweatShop *ss = new sweatShop(cmReader, cmWorker, cmWriter);

  ss->setLoaderQueueSize(1048576);
  ss->setWriterQueueSize(65536);

  ss->setLoaderBatchSize(128);
  ss->setWorkerBatchSize(128);

  ss->setNumberOfWorkers(numThreads);

  for (u32bit w=0; w<numThreads; w++)
    ss->setThreadData(w, new cmThreadData());  //  these leak

  ss->run(g, true);

  delete g;
}
