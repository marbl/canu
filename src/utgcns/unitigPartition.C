
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "unitigPartition.H"



void
tigPartitioning::loadTigInfo(tgStore *tigStore) {

  _nTigs  = tigStore->numTigs();
  _nReads = 0;

  _tigInfo.reserve(_nTigs);

  for (uint32 ti=0; ti<_nTigs; ti++) {
    uint64  len = 0;   //  64-bit so we don't overflow the various
    uint64  nc  = 0;   //  multiplications below.

    //  If there's a tig here, load it and get the info.

    if (tigStore->isDeleted(ti) == false) {
      tgTig *tig = tigStore->loadTig(ti);

      len = tig->length();
      nc  = tig->numberOfChildren();

      _nReads += nc;

      tigStore->unloadTig(ti);
    }

    //  Initialize the tigInfo.  If no tig is here, all the fields will end
    //  up zero, and we'll not put it into a partition.

    _tigInfo[ti].tigID           = ti;
    _tigInfo[ti].tigLength       = len;
    _tigInfo[ti].tigChildren     = nc;

    _tigInfo[ti].consensusArea   = 0;
    _tigInfo[ti].consensusMemory = 0;

    _tigInfo[ti].partition       = 0;
  }
}


void
tigPartitioning::loadTigInfo(std::vector<tgTig *> &tigList) {

  _nTigs  = tigList.size();
  _nReads = 0;

  _tigInfo.resize(_nTigs);

  for (uint32 ti=0; ti<_nTigs; ti++) {
    uint64  len = 0;   //  64-bit so we don't overflow the various
    uint64  nc  = 0;   //  multiplications below.

    //  If there's a tig here, load it and get the info.

    if (tigList[ti]) {
      len = tigList[ti]->length();
      nc  = tigList[ti]->numberOfChildren();

      _nReads += nc;
    }

    //  Initialize the tigInfo.  If no tig is here, all the fields will end
    //  up zero, and we'll not put it into a partition.

    if (tigList[ti])
      assert(ti == tigList[ti]->tigID());

    _tigInfo[ti].tigID           = ti;
    _tigInfo[ti].tigLength       = len;
    _tigInfo[ti].tigChildren     = nc;

    _tigInfo[ti].consensusArea   = 0;
    _tigInfo[ti].consensusMemory = 0;

    _tigInfo[ti].partition       = 0;
  }
}




void
tigPartitioning::greedilyPartition(double   partitionSizeScale,
                                   double   tigLengthScale,
                                   double   maxReadsPer,
                                   bool     verbose) {

  //  Compute the effort 'area' and estimated memory for each tig.

  for (uint32 ti=0; ti<_tigInfo.size(); ti++) {
    _tigInfo[ti].consensusArea   = _tigInfo[ti].tigLength * tigLengthScale * _tigInfo[ti].tigChildren;
    _tigInfo[ti].consensusMemory = _tigInfo[ti].tigLength * tigLengthScale * 1024;
  }

  //  Sort the tigInfo by decreasing area.

  auto byArea = [](tigInfo const &A, tigInfo const &B) {
    return(A.consensusArea > B.consensusArea);
  };

  std::sort(_tigInfo.begin(), _tigInfo.end(), byArea);

  //  Grab the biggest tig (it's number 0) and compute a maximum area per partition.

  uint64   maxArea         = _tigInfo[0].consensusArea * partitionSizeScale;
  uint32   maxReads        = (uint32)ceil(_nReads * maxReadsPer);
  uint32   currentPart     = 1;
  uint64   currentArea     = 0;
  uint32   currentTigs     = 0;
  uint32   currentChildren = 0;
  bool     stillMore       = true;

  if (maxArea == 0)
    maxArea = uint64max;

  if (verbose) {
    //fprintf(stderr, "\n");
    //fprintf(stderr, "maxArea  = %s\n", 
    //fprintf(stderr, "maxReads = %u\n",  maxReads);
    //fprintf(stderr, "\n");
    fprintf(stderr, "--------------------- TIG -------------------  ------- PARTITION --------  maxArea:  %s\n", (maxArea != uint64max) ? toDec(maxArea) : "infinite");
    fprintf(stderr, "    ID   Reads    Length         Area Mem(GB)    ID   Total Area TotReads  maxReads: %u\n", maxReads);
    fprintf(stderr, "------ ------- --------- ------------ -------  ---- ------------ --------\n");
  }

  while (stillMore) {
    stillMore = false;

    for (uint32 ti=0; ti<_tigInfo.size(); ti++) {

      //  Do nothing, it's already in a partition or if the area is zero.
      if      ((_tigInfo[ti].partition     != 0) ||
               (_tigInfo[ti].consensusArea == 0)) {
      }

      //  If nothing in the current partition, or still space in the
      //  partition, add this tig to the current partition.
      //
      //  In particular, this allows partitions of a single tig to be larger
      //  than the maximum (e.g., a partitionSize < 1.0).
      //
      //  It also ensures we don't assign too many (singleton) reads to a
      //  partition.

      else if ((currentTigs == 0) ||
               ((currentArea     + _tigInfo[ti].consensusArea < maxArea) &&
                (currentChildren + _tigInfo[ti].tigChildren   < maxReads))) {
        _tigInfo[ti].partition = currentPart;

        currentArea     += _tigInfo[ti].consensusArea;
        currentTigs     += 1;
        currentChildren += _tigInfo[ti].tigChildren;

        if (verbose)
          fprintf(stderr, "%6u %7lu %9lu %12lu %7.3f  %4u %12lu %8u\n",
                  _tigInfo[ti].tigID,
                  _tigInfo[ti].tigChildren,
                  _tigInfo[ti].tigLength,
                  _tigInfo[ti].consensusArea,
                  _tigInfo[ti].consensusMemory / 1024.0 / 1024.0 / 1024.0,
                  _tigInfo[ti].partition,
                  currentArea,
                  currentChildren);
      }

      //  Do nothing.  This tig is too large for the current partition.
      else {
        stillMore = true;
      }
    }

    //  Nothing else will fit in this partition.  Move to the next.

    currentPart    += 1;
    currentArea     = 0;
    currentTigs     = 0;
    currentChildren = 0;
  }

  if (verbose)
    fprintf(stderr, "------ ------- --------- ------------ -------  ---- ------------ --------\n");

  _nPartitions = currentPart;
}



void
tigPartitioning::outputPartitions(sqStore *seqStore, tgStore *tigStore, char const *storeName) {
  std::map<uint32, uint32>   readToPart;

  //  Sort by tigID.

  auto byTigID = [](tigInfo const &A, tigInfo const &B) {
    return(A.tigID < B.tigID);
  };

  std::sort(_tigInfo.begin(), _tigInfo.end(), byTigID);

  //  Build a mapping from readID to partitionID.

  for (uint32 ti=0; ti<_tigInfo.size(); ti++) {
    if (_tigInfo[ti].partition > 0) {
      tgTig  *tig = tigStore->loadTig(_tigInfo[ti].tigID);

      for (uint32 fi=0; fi<tig->numberOfChildren(); fi++)
        readToPart[tig->getChild(fi)->ident()] = _tigInfo[ti].partition;

      tigStore->unloadTig(_tigInfo[ti].tigID);
    }
  }

  //  Create output files for each partition and write a small header.

  writeBuffer              **parts = new writeBuffer * [_nPartitions];
  char                       partName[FILENAME_MAX+1];

  for (uint32 pi=1; pi<_nPartitions; pi++) {
    snprintf(partName, FILENAME_MAX, "%s/partition.%04u", storeName, pi);

    parts[pi] = new writeBuffer(partName, "w");

    uint64  magc = (uint64)0x5f5f656c69467173llu;   //  'sqFile__'
    uint64  vers = (uint64)0x0000000000000001llu;
    uint64  defv = (uint64)sqRead_defaultVersion;

    parts[pi]->writeIFFobject("MAGC", magc);
    parts[pi]->writeIFFobject("VERS", vers);
    parts[pi]->writeIFFobject("DEFV", defv);
  }

  //  Scan the store, copying read data to partition files.

  sqRead                    *rd    = new sqRead;
  sqReadDataWriter          *wr    = new sqReadDataWriter;

  for (uint32 fi=1; fi<seqStore->sqStore_lastReadID()+1; fi++)
    if (readToPart.count(fi) > 0)
      seqStore->sqStore_saveReadToBuffer(parts[readToPart[fi]], fi, rd, wr);

  delete    wr;
  delete    rd;

  //  All done!  Cleanup.

  for (uint32 pi=1; pi<_nPartitions; pi++)
    delete parts[pi];

  delete [] parts;
}


void
tigPartitioning::reportPartitioning(FILE *partFile) {

  fprintf(partFile, "      Tig     Reads    Length         Area  Memory GB  Partition\n");
  fprintf(partFile, "--------- --------- --------- ------------  ---------  ---------\n");

  for (uint32 ti=0; ti<_tigInfo.size(); ti++)
    if (_tigInfo[ti].partition != 0)
      fprintf(partFile, "%9u %9lu %9lu %12lu  %9.3f  %9u\n",
              _tigInfo[ti].tigID,
              _tigInfo[ti].tigChildren,
              _tigInfo[ti].tigLength,
              _tigInfo[ti].consensusArea,
              _tigInfo[ti].consensusMemory / 1024.0 / 1024.0 / 1024.0,
              _tigInfo[ti].partition);
}
