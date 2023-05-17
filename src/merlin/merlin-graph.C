
/******************************************************************************
 *
 *  This is a tool for linking haplotype specific k-mers in long-reads.
 *
 *  This software is based on:
 *    'Meryl'                  (https://github.com/marbl/meryl)
 *
 *  This is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "merlin-globals.H"

/****
 * This stub builds the graph
 * Input: global G and input S
 * Output: print the graph structure as-seen
 */ 
void
processBuild(void *G, void *T, void *S) {
  merlinGlobal  *g = (merlinGlobal  *)G;
  merlinThrData *t = (merlinThrData *)T;
  merlinInput   *s = (merlinInput   *)S;

  kmer  mer;
  char  orientMer     = node::MER_REV;
  char kmerstring[65]; // temp for retrieving kmer sequence in debug

  //  Print debug
  if (g->debug == true && t->oDebug == nullptr) {
    char  name[FILENAME_MAX+1];
    snprintf(name, FILENAME_MAX, "%s.%02d.debug.gz", g->outName, t->threadID);
    t->oDebug = new compressedFileWriter(name);
  }

  while (s->kiter.nextBase()) {
    if (s->kiter.isValid() == false)
      continue;

    s->totalK++;

    //  Look up the mer in s exists in marker database.
    int val = g->lookupMarker(s->kiter.fmer(), s->kiter.rmer());
    if (val > 0) {
      if (s->kiter.fmer() <= s->kiter.rmer()) {
        mer = s->kiter.fmer();
        orientMer = node::MER_FWD;
      } else {
        mer = s->kiter.rmer();
        orientMer = node::MER_REV;
      }

      //  for debug: add kmers found
      s->mersPos.push_back(s->kiter.bgnPosition());
      s->mersFound.push_back(mer);
      s->mersOrient.push_back(orientMer);
      s->mersValue.push_back(val);
    }
  }

  //  Output debug.
  if (g->debug) {
    char kmerstring[65];  // for temporarily storing the string form of the kmer
    fprintf(t->oDebug->file(), "%s\t", s->seq.ident());

    //  Only output pos and mers when there are >= 2 mers
    if (s->mersPos.size() > 1) {
      fprintf(t->oDebug->file(), "%d", s->mersPos[0]);
      for ( int ii = 1; ii < s->mersPos.size(); ii++)
        fprintf(t->oDebug->file(), ":%d", s->mersPos.at(ii));

      fprintf(t->oDebug->file(), "\t");

      fprintf(t->oDebug->file(), "%s", s->mersFound[0].toString(kmerstring));
      for ( int ii = 1; ii < s->mersFound.size(); ii++) 
        fprintf(t->oDebug->file(), ":%s", s->mersFound[ii].toString(kmerstring));

    }

    fprintf(t->oDebug->file(), "\n");
  }

  fprintf(stderr, "%s has %d markers out of %d total mers\n",
      s->seq.ident(), s->mersFound.size(), s->totalK);

}


//  Update graph
//  This is called after each process thread, one at a time.
//  Now it's time to define edges and placing them to g.
//
void
updateGraph(void *G, void *S) {
  
  merlinGlobal  *g = (merlinGlobal *)G;
  merlinInput   *s = (merlinInput *)S;

  bool isFirst = true;  // is this the first marker found?
  kmer  mer;
  kmer  prevMer; // pointer to the prev. mer node
  char  orientMer     = node::MER_REV;
  char  orientPrevMer = node::MER_REV;
  node *curNode;
  node *prevNode;

  for (int ii = 0; ii < s->mersFound.size(); ii++) {
    mer = s->mersFound.at(ii);
    orientMer = s->mersOrient.at(ii);

    if (g->nodes.find(mer) == g->nodes.end()) {
      curNode = new node();
      curNode->setNode(mer, s->mersValue.at(ii));
      g->nodes.insert(pair<kmer, node*>(mer, curNode));
      //  fprintf(stderr, "[ DEBUG ] :: Insert new node %s\n", curNode->nodeId.toString(kmerstring));
    } else {
      curNode = g->nodes.at(mer);
      //  fprintf(stderr, "[ DEBUG ] :: Node %s found.\n", curNode->nodeId.toString(kmerstring));
    }

    if (isFirst) {
      isFirst = false; // check this off, we will start with the second one
    } else {
      //  we can make a new edge, or add count
      //  fprintf(stderr, "[ DEBUG ] :: adding new edge between nodes: begin\n");
      // if (prevNode->nodeId <= curNode->nodeId) {
      prevNode->addEdge(curNode, orientPrevMer, orientMer);
      // something is wrong here. let's keep the order as it was
      //} else {
      //  curNode->addEdge(prevNode, orientMer, orientPrevMer);
      //}
      //  fprintf(stderr, "[ DEBUG ] :: adding new edge between nodes: done\n");
    }

    prevMer  = mer;
    prevNode = curNode;
    orientPrevMer = orientMer;
  }

}
