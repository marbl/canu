
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * Copyright (C) 2007, J. Craig Venter Institute.
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

const char *mainid = "$Id: cgwDump.C,v 1.1 2012-10-25 17:15:54 brianwalenz Exp $";

#include "AS_global.h"

#include "Globals_CGW.h"
#include "AS_CGW_dataTypes.h"
#include "Output_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"

#include "AS_UTL_fasta.h"
#include "AS_UTL_reverseComplement.h"

#define DUMP_NOTHING       0
#define DUMP_READS         1
#define DUMP_UNITIGS       2
#define DUMP_CONTIGS       3
#define DUMP_SCAFFOLDS     4

#define DUMP_CONSENSUS     0x01
#define DUMP_LAYOUT        0x02
#define DUMP_EDGES         0x04


class outputObj {
public:
  outputObj() {
    id = 0;
    fwd = false;
    gapBgn = gapEnd = posBgn = posEnd = 0;
  };

  uint32             id;
  bool               fwd;

  int32              gapBgn;  //  Gapped coordinates
  int32              gapEnd;

  int32              posBgn;  //  Ungapped output coordinates
  int32              posEnd;
};


class outputTig {
public:
  outputTig() {};
  outputTig(NodeCGW_T *tig, MultiAlignT *ma, bool isUnitig) {
    init(tig, ma, isUnitig);
  };

  void    init(NodeCGW_T *tig, MultiAlignT *ma, bool isUnitig);

  int32   setPosition(int32 bgn);

  outputObj          obj;  //  This unitig or contig
  vector<outputObj>  frg;  //  Component fragments
  vector<outputObj>  utg;  //  Component unitigs

  vector<char>       seq;
};


class outputScf {
public:
  outputScf() {};
  outputScf(NodeCGW_T *scf) {
    init(scf);
  };

  void    init(NodeCGW_T *scaffold);

  vector<outputObj>  frg;  //  Component fragments
  vector<outputTig>  utg;  //  Component unitigs
  vector<outputTig>  ctg;  //  Component contigs

  vector<int32>      gap;  //  Gap from contig [i] to contig [i+1].

  vector<char>       seq;
};






void
outputTig::init(NodeCGW_T    *tig,
                MultiAlignT  *ma,
                bool          isUnitig) {
  outputObj   tmpObj;

  obj.id     = tig->id;
  obj.fwd    = (tig->offsetAEnd.mean < tig->offsetBEnd.mean);

  //  Build a map from gapped to ungapped coordinates.

  char    *gapseq     = Getchar(ma->consensus, 0);
  uint32   seqLen     = GetMultiAlignLength(ma);
  uint32  *gapToUngap = new uint32 [seqLen];

  if (obj.fwd == false)
    reverseComplementSequence(gapseq, seqLen);

  uint32   gp = 0;  //  Gapped length, after the loop
  uint32   up = 0;  //  Ungapped length

  for (; gp<seqLen; gp++) {
    gapToUngap[gp] = up;

    if (gapseq[gp] != '-') {
      seq.push_back(gapseq[gp]);
      up++;
    }
  }

  //  Set the position of the tig in the parent.  For now, set the ungapped
  //  position to the start of the parent; the parent (scaffold) must update this.

  obj.gapBgn = (obj.fwd) ? tig->offsetAEnd.mean : tig->offsetBEnd.mean;
  obj.gapEnd = (obj.fwd) ? tig->offsetBEnd.mean : tig->offsetAEnd.mean;

  obj.posBgn = 0;
  obj.posEnd = up;

  //  Populate a list of reads.

  for (uint32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
    IntMultiPos *imp = GetIntMultiPos(ma->f_list, i);

    tmpObj.id     = imp->ident;
    tmpObj.fwd    = (imp->position.bgn < imp->position.end) ? true : false;
    tmpObj.gapBgn = imp->position.bgn;
    tmpObj.gapEnd = imp->position.end;
    tmpObj.posBgn = gapToUngap[tmpObj.gapBgn];
    tmpObj.posEnd = gapToUngap[tmpObj.gapEnd];

    frg.push_back(tmpObj);
  }

  //  Populate a list of unitigs.

  for (uint32 i=0; i<GetNumIntUnitigPoss(ma->u_list); i++) {
    IntUnitigPos *iup = GetIntUnitigPos(ma->u_list, i);

    tmpObj.id     = iup->ident;
    tmpObj.fwd    = (iup->position.bgn < iup->position.end) ? true : false;
    tmpObj.gapBgn = iup->position.bgn;
    tmpObj.gapEnd = iup->position.end;
    tmpObj.posBgn = gapToUngap[tmpObj.gapBgn];
    tmpObj.posEnd = gapToUngap[tmpObj.gapEnd];

    frg.push_back(tmpObj);
  }

  //  Cleanup.

  delete [] gapToUngap;
}


int32
outputTig::setPosition(int32 bgn) {

  obj.posBgn += bgn;
  obj.posEnd += bgn;

  for (uint32 fi=0; fi<frg.size(); fi++) {
    frg[fi].posBgn += bgn;
    frg[fi].posEnd += bgn;
  }

  for (uint32 ui=0; ui<utg.size(); ui++) {
    utg[ui].posBgn += bgn;
    utg[ui].posEnd += bgn;
  }

  return(obj.posEnd);
}



void
outputScf::init(NodeCGW_T *scaffold) {
  CIScaffoldTIterator    contigs;
  ChunkInstanceT        *contig;

  uint32  nFrg = 0;
  uint32  nUtg = 0;
  uint32  nCtg = 0;

  //  Pass 1: find gap sizes between contigs.

  double   lastEnd = 0;

  InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &contigs);
  while ((contig = NextCIScaffoldTIterator(&contigs)) != NULL) {
    bool          isFwd   = (contig->offsetAEnd.mean < contig->offsetBEnd.mean);
    double        bgn     = (isFwd) ? contig->offsetAEnd.mean : contig->offsetBEnd.mean;
    double        end     = (isFwd) ? contig->offsetBEnd.mean : contig->offsetAEnd.mean;

    if (lastEnd > 0) {
      if (bgn < lastEnd + 20)
        gap.push_back(20);
      else
        gap.push_back((int32)(bgn - lastEnd));

      fprintf(stderr, "GAP %u %d from last %.2f this %.2f %.2f\n", gap.size(), gap.back(), lastEnd, bgn, end);
    }

    lastEnd = end;
  }

  gap.push_back(0);  //  For convenience, a zero base gap after the last contig

  //  Pass 2: Populate frg, utg and ctg vectors.

  InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &contigs);
  while ((contig = NextCIScaffoldTIterator(&contigs)) != NULL) {
    ctg.push_back(outputTig());
    ctg.back().init(contig,
                    ScaffoldGraph->tigStore->loadMultiAlign(contig->id, FALSE),
                    FALSE);
  }

  //  Pass 3: Place each contig in the correct location.  This also adjusts the component
  //  reads and unitigs.

  int32 scfLen = 0;

  for (uint32 ci=0; ci<ctg.size(); ci++)
    scfLen = ctg[ci].setPosition(scfLen) + gap[ci];

  //  Pass 3: Construct scaffold sequence

  for (uint32 ci=0; ci<ctg.size(); ci++) {
    for (uint32 ss=0; ss<ctg[ci].seq.size(); ss++)
      seq.push_back(ctg[ci].seq[ss]);

    for (uint32 gg=0; gg<gap[ci]; gg++)
      seq.push_back('n');
  }

  seq.push_back(0);

  //AS_UTL_writeFastA(stderr, &seq[0], seq.size()-1, seq.size()-1, ">scf"F_U32"\n", scaffold->id);
  //fprintf(stderr, "ctg%08d\tscf%08d\t%d\t%d\t%c\n", contig->id, scaffold->id, bgn, end, (isFwd) ? 'f' : 'r');
}













FILE *
openPrefixFile(char *prefix, char *label) {
  char     N[FILENAME_MAX];
  FILE    *F;

  sprintf(N, "%s%s", prefix, label);

  errno = 0;
  F = fopen(N, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for write: %s\n", N, strerror(errno)), exit(1);

  return(F);
}


void
dumpReads(uint32          bgnIID,
          uint32          endIID,
          vector<uint32>  objIID,
          uint32          dumpType,
          char           *outputPrefix) {
}



void
dumpUnitigs(uint32          bgnIID,
            uint32          endIID,
            vector<uint32>  objIID,
            uint32          dumpType,
            char           *outputPrefix) {

  FILE    *outputUnique  = openPrefixFile(outputPrefix, ".utg.unique.fasta");
  FILE    *outputSingle  = openPrefixFile(outputPrefix, ".utg.single.fasta");
  FILE    *outputRock    = openPrefixFile(outputPrefix, ".utg.rock.fasta");
  FILE    *outputStone   = openPrefixFile(outputPrefix, ".utg.stone.fasta");
  FILE    *outputPebble  = openPrefixFile(outputPrefix, ".utg.pebble.fasta");

  FILE    *outputLayout  = openPrefixFile(outputPrefix, ".utg.layout");
  FILE    *outputEdges   = openPrefixFile(outputPrefix, ".utg.edges");

  bool     isVec = (objIID.size() > 0);

  if (endIID > GetNumGraphNodes(ScaffoldGraph->CIGraph))
    endIID = GetNumGraphNodes(ScaffoldGraph->CIGraph);

  if (isVec) {
    bgnIID = 0;
    endIID = objIID.size();
  }

  for (uint32 ii=bgnIID; ii<endIID; ii++) {
    NodeCGW_T    *unitig = GetGraphNode(ScaffoldGraph->CIGraph, (isVec) ? objIID[ii] : ii);

    if (unitig->flags.bits.isDead)
      continue;

    MultiAlignT  *ma     = ScaffoldGraph->tigStore->loadMultiAlign(unitig->id, TRUE);
    UnitigType    type   = finalUnitigType(unitig);

    uint32  seqLen = GetMultiAlignUngappedLength(ma);
    char   *seq    = new char [seqLen + 1];
    char   *qlt    = new char [seqLen + 1];

    GetMultiAlignUngappedConsensus(ma, seq, qlt);

    if (dumpType & DUMP_LAYOUT) {
    }

    if (dumpType & DUMP_EDGES) {
    }

    if (dumpType & DUMP_CONSENSUS) {
      if (type == AS_UNIQUE_UNITIG)
        AS_UTL_writeFastA(outputUnique, seq, seqLen, seqLen, ">utg"F_U32"\n", unitig->id);

      if (type == AS_SINGLE_UNITIG)
        AS_UTL_writeFastA(outputSingle, seq, seqLen, seqLen, ">utg"F_U32"\n", unitig->id);

      if (type == AS_ROCK_UNITIG)
        AS_UTL_writeFastA(outputRock,   seq, seqLen, seqLen, ">utg"F_U32"\n", unitig->id);

      if (type == AS_STONE_UNITIG)
        AS_UTL_writeFastA(outputStone,  seq, seqLen, seqLen, ">utg"F_U32"\n", unitig->id);

      if (type == AS_PEBBLE_UNITIG)
        AS_UTL_writeFastA(outputPebble, seq, seqLen, seqLen, ">utg"F_U32"\n", unitig->id);
    }
  }

  if (outputUnique)  fclose(outputUnique);
  if (outputSingle)  fclose(outputSingle);
  if (outputRock)    fclose(outputRock);
  if (outputStone)   fclose(outputStone);
  if (outputPebble)  fclose(outputPebble);

  if (outputLayout)  fclose(outputLayout);
  if (outputEdges)   fclose(outputEdges);
}



void
dumpContigs(uint32          bgnIID,
            uint32          endIID,
            vector<uint32>  objIID,
            uint32          dumpType,
            char           *outputPrefix) {
}



void
dumpScaffolds(uint32          bgnIID,
              uint32          endIID,
              vector<uint32>  objIID,
              uint32          dumpType,
              char           *outputPrefix) {
  FILE    *outputMulti   = openPrefixFile(outputPrefix, ".scf.multi.fasta");
  FILE    *outputSingle  = openPrefixFile(outputPrefix, ".scf.single.fasta");

  FILE    *outputLayout  = openPrefixFile(outputPrefix, ".scf.layout");
  FILE    *outputEdges   = openPrefixFile(outputPrefix, ".scf.edges");

  bool     isVec = (objIID.size() > 0);

  if (endIID > GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph))
    endIID = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);

  if (isVec) {
    bgnIID = 0;
    endIID = objIID.size();
  }

  for (uint32 ii=bgnIID; ii<endIID; ii++) {
    NodeCGW_T             *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, (isVec) ? objIID[ii] : ii);

    if ((scaffold->type              != REAL_SCAFFOLD) ||
        (scaffold->flags.bits.isDead == true))
      continue;

    assert(scaffold->info.Scaffold.numElements > 0);

    outputScf              scf(scaffold);

    if (dumpType & DUMP_CONSENSUS) {
      if (scf.ctg.size() == 1)
        AS_UTL_writeFastA(outputSingle, &scf.seq[0], scf.seq.size()-1, scf.seq.size()-1, ">scf"F_U32"\n", scaffold->id);

      if (scf.ctg.size() > 1)
        AS_UTL_writeFastA(outputMulti, &scf.seq[0], scf.seq.size()-1, scf.seq.size()-1, ">scf"F_U32"\n", scaffold->id);
    }
  }

  if (outputMulti)   fclose(outputMulti);
  if (outputSingle)  fclose(outputSingle);

  if (outputLayout)  fclose(outputLayout);
  if (outputEdges)   fclose(outputEdges);
}



int main (int argc, char *argv[]) {
  int32           checkpointVers           = 0;
  int32           tigStoreVers             = 0;

  char           *outputPrefix             = NULL;

  uint32          objectType               = DUMP_NOTHING;

  uint32          bgnIID                   = 0;
  uint32          endIID                   = UINT32_MAX;
  vector<uint32>  objIID;

  uint32          dumpType                 = DUMP_NOTHING;


  GlobalData = new Globals_CGW();

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      strcpy(GlobalData->gkpStoreName, argv[++arg]);

    } else if (strcmp(argv[arg], "-t") == 0) {
      strcpy(GlobalData->tigStoreName, argv[++arg]);
      tigStoreVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-c") == 0) {
      strcpy(GlobalData->outputPrefix, argv[++arg]);
      checkpointVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      outputPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-reads") == 0) {
      objectType = DUMP_READS;
    } else if (strcmp(argv[arg], "-unitigs") == 0) {
      objectType = DUMP_UNITIGS;
    } else if (strcmp(argv[arg], "-contigs") == 0) {
      objectType = DUMP_CONTIGS;
    } else if (strcmp(argv[arg], "-scaffolds") == 0) {
      objectType = DUMP_SCAFFOLDS;

    } else if (strcmp(argv[arg], "-b") == 0) {
      bgnIID = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-e") == 0) {
      endIID = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-i") == 0) {
      objIID.push_back(atoi(argv[++arg]));

    } else if (strcmp(argv[arg], "-consensus") == 0) {
      dumpType |= DUMP_CONSENSUS;
    } else if (strcmp(argv[arg], "-layout") == 0) {
      dumpType |= DUMP_LAYOUT;
    } else if (strcmp(argv[arg], "-edges") == 0) {
      dumpType |= DUMP_EDGES;

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }

  if (GlobalData->gkpStoreName[0] == 0)
    err++;
  if (GlobalData->tigStoreName[0] == 0)
    err++;
  if (GlobalData->outputPrefix[0] == 0)
    err++;
  if (outputPrefix == NULL)
    err++;

  if (err > 0) {
    fprintf(stderr, "usage: %s -g gkpStore -t tigStore version -c checkpoint version [action] -o prefix\n", argv[0]);
    fprintf(stderr, "  -g gkpStore             mandatory path to the gkpStore\n");
    fprintf(stderr, "  -t tigStore version     mandatory path to the tigStore and version\n");
    fprintf(stderr, "  -c checkpoint version   mandatory path to a checkpoint and version\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o prefix               output is written to files starting with 'prefix'\n");
    fprintf(stderr, "                          (e.g., prefix.ctg.fasta)\n");
    fprintf(stderr, "                          (e.g., prefix.posmap.frgscf)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "action - what object types to dump\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -reads                  dumps reads\n");
    fprintf(stderr, "  -unitigs                dumps unitigs\n");
    fprintf(stderr, "  -contigs                dumps contigs\n");
    fprintf(stderr, "  -scaffolds              dumps scaffolds\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "selection - what objects to dump\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -b bgnIID               dumps objects bgnIID <= IID <= endIID\n");
    fprintf(stderr, "  -e endIID               \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -i singleIID            dumps a specific object (multiple -i allowed)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "selection - what data to dump\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -consensus              dumps consensus sequence\n");
    fprintf(stderr, "  -layout                 dumps posmap layout of component objects\n");
    fprintf(stderr, "  -edges                  dumps unused mate edges involving selected objects\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "processing - cleanup / label before dumping\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -labeltigs              assign labels to unitigs/contigs\n");
    fprintf(stderr, "  -labelreads             assign labels to mate pairs\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    if (GlobalData->gkpStoreName[0] == 0)
      fprintf(stderr, "ERROR: No gatekeeper store (-g) supplied.\n");
    if (GlobalData->tigStoreName[0] == 0)
      fprintf(stderr, "ERROR: No tigStore store (-t) supplied.\n");
    if (GlobalData->outputPrefix[0] == 0)
      fprintf(stderr, "ERROR: No checkpoint (-c) supplied.\n");
    if (outputPrefix == NULL)
      fprintf(stderr, "ERROR: No output prefix (-o) supplied.\n");
    exit(1);
  }

  LoadScaffoldGraphFromCheckpoint(GlobalData->outputPrefix, checkpointVers, FALSE);

  //  If the tigStore version is different than the checkpoint version, reopen the tigStore
  if (tigStoreVers != checkpointVers) {
    delete ScaffoldGraph->tigStore;
    ScaffoldGraph->tigStore = new MultiAlignStore(GlobalData->tigStoreName, tigStoreVers, 0, 0, FALSE, FALSE);
  }

  bool labelTigs  = true;
  bool labelMates = true;

  if (labelTigs) {
  }

  if (labelMates)
    MarkContigEdges();

  if (objectType == DUMP_READS)
    dumpReads(bgnIID, endIID, objIID, dumpType, outputPrefix);

  if (objectType == DUMP_UNITIGS)
    dumpUnitigs(bgnIID, endIID, objIID, dumpType, outputPrefix);

  if (objectType == DUMP_CONTIGS)
    dumpContigs(bgnIID, endIID, objIID, dumpType, outputPrefix);

  if (objectType == DUMP_SCAFFOLDS)
    dumpScaffolds(bgnIID, endIID, objIID, dumpType, outputPrefix);

  DestroyScaffoldGraph(ScaffoldGraph);

  return(0);
}
