
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
 *  Modifications by:
 *
 *    Sergey Koren beginning on 2018-APR-18
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Brian P. Walenz beginning on 2018-APR-18
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "runtime.H"
#include "ovStore.H"
#include "strings.H"
#include "tgStore.H"

#include <vector>
#include <map>

using namespace std;

double MIN_READ_FRACTION = 0.5;
double MAX_READ_STRETCH  = 1.2;

void save_tig(sqStore *seqStore, tgStore *tigStore, tgTig *tig,
              map<uint32, map<uint32, int32> > &readToStart,
              map<uint32, map<uint32, int32> > &readToEnd,
              map<uint32, map<uint32, bool> >   &readToOri,
              map<uint32, bool>                 &readUsed,
              map<uint32, map<uint32, double> > &readFraction,
              map<uint32, map<uint32, uint32> > &readPieces) {
  if (tig->_layoutLen > 0) {
    resizeArray(tig->_children, tig->_childrenLen, tig->_childrenMax, readToOri.size(), resizeArray_doNothing);

    // first loop finds min value from those we placed here to make new 0
    int32 minOffset = 0;

   /*

    for (map<uint32, map<uint32, bool> >::iterator it=readToOri.begin(); it != readToOri.end(); ++it) {
      // now we pick the best location
      uint32 best = 0;
      uint32 bestCount = 0;
      for (map<uint32, uint32>::iterator chunks=readPieces[it->first].begin(); chunks != readPieces[it->first].end(); ++chunks) {
        if (bestCount < chunks->second) {
          // sanity check the placement
          if (readToEnd[it->first][chunks->first] - readToStart[it->first][chunks->first] < MAX_READ_STRETCH*seqStore->sqStore_getRead(it->first)->sqRead_length() && readFraction[it->first][chunks->first] > MIN_READ_FRACTION) {
            bestCount = chunks->second;
            best = chunks->first;
          }
        }
      }
      if (bestCount == 0) continue; // we couldn't find a good match covering the read in this tig, skip it, worst case it ends up as chaff
      if (readToStart[it->first][best] < minOffset) minOffset = readToStart[it->first][best];
    }
    */
    minOffset *= -1;
    fprintf(stderr, "Set min offset to be %d\n", minOffset);

    for (map<uint32, map<uint32, bool> >::iterator it=readToOri.begin(); it != readToOri.end(); ++it) {
      // now we pick the best location
      uint32 best = 0;
      uint32 bestCount = 0;
      for (map<uint32, uint32>::iterator chunks=readPieces[it->first].begin(); chunks != readPieces[it->first].end(); ++chunks) {
        if (bestCount < chunks->second) {
          // sanity check the placement
          if (readToEnd[it->first][chunks->first] - readToStart[it->first][chunks->first] < MAX_READ_STRETCH*seqStore->sqStore_getReadLength(it->first)
              && readToEnd[it->first][chunks->first] - readToStart[it->first][chunks->first] > seqStore->sqStore_getReadLength(it->first) / MAX_READ_STRETCH
              && readFraction[it->first][chunks->first] > MIN_READ_FRACTION) {
            bestCount = chunks->second;
            best = chunks->first;
          }
        }
      }
      if (bestCount == 0) continue; // we couldn't find a good match covering the read in this tig, skip it, worst case it ends up as chaff

      fprintf(stderr, "For read %d picked best index %d with %d chunks which is at positions %d-%d\n", it->first, best, bestCount, readToStart[it->first][best], readToEnd[it->first][best]);
      readUsed[it->first] = true;

      if ((it->second)[best] == true)
        tig->addChild()->set(it->first,
                      0, 0, 0,
                      minOffset + readToStart[it->first][best], minOffset + readToEnd[it->first][best]);
      else
        tig->addChild()->set(it->first,
                      0, 0, 0,
                      minOffset + readToEnd[it->first][best], minOffset + readToStart[it->first][best]);
    }
    if (tig->_childrenLen > 0) tigStore->insertTig(tig, false);
  }
  readToStart.clear();
  readToEnd.clear();
  readToOri.clear();
  readFraction.clear();
  readPieces.clear();
}

int
main(int argc, char **argv) {
  char           *outName  = NULL;
  char           *seqName  = NULL;
  uint32          minOverlapLength = 0;

  vector<char *>  files;

  int32     arg = 1;
  int32     err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      outName = argv[++arg];

    } else if (strcmp(argv[arg], "-S") == 0) {
      seqName = argv[++arg];

    } else if (fileExists(argv[arg])) {
      files.push_back(argv[arg]);

    } else {
      fprintf(stderr, "ERROR:  invalid arg '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if ((err) || (seqName == NULL) || (outName == NULL) || (files.size() == 0)) {
    fprintf(stderr, "usage: %s [options] file.dbg.lay[.gz]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  Converts wtdbg layout to tigStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o out     output prefix\n");
    fprintf(stderr, "\n");

    if (seqName == NULL)
      fprintf(stderr, "ERROR:  no seqStore (-S) supplied\n");
    if (files.size() == 0)
      fprintf(stderr, "ERROR:  no overlap files supplied\n");

    exit(1);
  }

  char        *ovStr = new char [1024*1024];

  sqStore    *seqStore = new sqStore(seqName);
  char        filename[FILENAME_MAX] = {0};
  snprintf(filename, FILENAME_MAX, "%s.%sStore", outName, "ctg");
  tgStore     *tigStore = new tgStore(filename);
  tgTig       *tig      = new tgTig;

  // the wtdbg layout breaks read into 2kb pieces
  // each read may be split across multiple contigs and multiple times in a contig
  // we only use the read in the contig where it makes sense (large fraction covered and not too big span)
  // within a contig we select the position of the read with the most 2kbp pieces putting it there
  map<uint32, map<uint32, int32> > readToStart;
  map<uint32, map<uint32, int32> > readToEnd;
  map<uint32, map<uint32, bool> >   readToOri;
  map<uint32, bool>   readUsed;
  map<uint32, map<uint32, double> > readFraction;
  map<uint32, map<uint32, uint32> > readPieces;

  double      offset = 0;

  tig->clear();
  for (uint32 ff=0; ff<files.size(); ff++) {
    compressedFileReader  *in = new compressedFileReader(files[ff]);

    while (fgets(ovStr, 1024*1024, in->file()) != NULL) {
       splitToWords  W(ovStr);

       if (ovStr[0] == '>') {
          save_tig(seqStore, tigStore, tig, readToStart, readToEnd, readToOri, readUsed, readFraction, readPieces);

          offset = 0;
          readToStart.clear();
          readToEnd.clear();
          tig->clear();
          tig->_tigID = atoi(W[0]+4);

          //  Set the class and some flags.

          tig->_class           = tgTig_contig;
          tig->_suggestRepeat   = false;
          tig->_suggestCircular = false;

          tig->_layoutLen       = atoi(W[2]+4);
       } if (ovStr[0] == 'E') {
          offset = W.toint32(1) * 1.10;
          fprintf(stderr, "The offset is updated to be %f\n", offset);
       } if (ovStr[0] == 'S' || ovStr[0] == 's') {
          uint32 rid = atoi(W[1]+4);
          uint32 rLen = strlen(W[1]);
          uint32 index = 0;
          char *c = W[1]+rLen-3;
          if (*c == '_') {
             index = atoi(W[1]+rLen-1);
          }
          fprintf(stderr, "The char is %c for string %s which made index %d\n", *c, W[1], index);

          if (readUsed.find(rid) != readUsed.end()) {
             continue;
          }

          int32 bgn = 0;
          int32 end = 0;

          if (W[2][0] == '+') {
             readToOri[rid][index] = true;
             bgn            = (int)(offset) - W.toint32(3);
             end            = int(offset) + W.toint32(4) + seqStore->sqStore_getReadLength(rid) - (W.toint32(3) + W.toint32(4));
          } else if (W[2][0] == '-') {
             readToOri[rid][index] = false;
             bgn            = int(offset) - (seqStore->sqStore_getReadLength(rid) - (W.toint32(3) + W.toint32(4)));
             end            = int(offset) + W.toint32(4);
          }
         if (readToStart.find(rid) == readToStart.end() || readToStart[rid].find(index) == readToStart[rid].end()) {
             readToStart[rid][index] = max(0, bgn);
             readToEnd[rid][index] = end;
             readPieces[rid][index] = 1;
             readFraction[rid][index] = W.todouble(4) / seqStore->sqStore_getReadLength(rid);
             fprintf(stderr, "Initialized read %d at index %d of length %d at offset %f to %d-%d\n", rid, index, seqStore->sqStore_getReadLength(rid), offset, bgn, end);
          }
          if (readToEnd[rid][index] < end) {
             readToEnd[rid][index] = end;
             ++readPieces[rid][index];
             readFraction[rid][index] += W.todouble(4) / seqStore->sqStore_getReadLength(rid);
             fprintf(stderr, "Updated read %d at index %d to %d-%d based on %d and %d of length %d\n", rid, index, readToStart[rid][index], end, W.toint32(3), W.toint32(4), seqStore->sqStore_getReadLength(rid));
          }
       }
    }
  }
  save_tig(seqStore, tigStore, tig, readToStart, readToEnd, readToOri, readUsed, readFraction, readPieces);

  delete tig;
  delete tigStore;

  delete seqStore;

  exit(0);
}
