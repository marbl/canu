
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
#ifndef INSTRUMENT_CGW_H
#define INSTRUMENT_CGW_H

#include <math.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_Hash.h"
#include "InputDataTypes_CGW.h"
#include "ScaffoldGraph_CGW.h"

VA_DEF(cds_float32);
VA_DEF(IntContigPairs);



/*
  Bit vector switches for instrumenting options - field in instrumenter
*/
#define INST_OPT_INTRA_MATES    1
#define INST_OPT_INTER_MATES    2
#define INST_OPT_ALL_MATES      (INST_OPT_INTRA_MATES | \
                                 INST_OPT_INTER_MATES)
#define INST_OPT_BREAKPOINTS    4
#define INST_OPT_FRAGMENTS      8
#define INST_OPT_CONTIG_PAIRS  16
#define INST_OPT_ALL       (INST_OPT_INTRA_MATES | \
                            INST_OPT_INTER_MATES | \
                            INST_OPT_BREAKPOINTS | \
                            INST_OPT_FRAGMENTS | \
                            INST_OPT_CONTIG_PAIRS)
#define INST_ONE_CELAMY        32
#define INST_ALL_CELAMY        64

#define TRACK_3P  /*store fragment 3-prime end locations*/

/*
  levels of verbosity
*/
typedef enum
  {
    InstrumenterSilent,
    InstrumenterVerbose1,
    InstrumenterVerbose2,
    InstrumenterVerbose3,
    InstrumenterVerbose4,
    InstrumenterVerbose5
  } InstrumenterVerbosity;

/* what to print for mate details for scaffold instrumenter */
typedef enum { PRINTCELAMY, PRINTTABLE }  printtypes;

/*
  types of results
*/
typedef enum
  {
    InstrumenterIndeterminate,
    InstrumenterSame,
    InstrumenterBetter,
    InstrumenterWorse
  } InstrumenterResult;

/*
  types of instrumenters
*/
typedef enum
  {
    InstrumenterUnitigLevel,
    InstrumenterContigLevel,
    InstrumenterScaffoldLevel,
    InstrumenterScaffoldGraphLevel
  } InstrumenterLevel;


/*
  InstrumenterStatistics
*/
typedef struct
{
  cds_int32  numNegatives;
  cds_int32  numPositives;
  cds_int32  num;
  cds_float32 sumOfSquares;
  cds_float32 min;
  cds_float32 minNegative;
  cds_float32 max;
  cds_float32 mean;
  cds_float32 stddev;
} InstrumenterStatistics;


/*
  Set of instrumenter orientations. Used for tracking orientations and
  for array indices
*/
typedef enum
  {
    INNIE_INSTR = 0,
    OUTTIE_INSTR,
    NORMAL_INSTR,
    ANTINORMAL_INSTR,
    NUM_ORIENTATIONS_INSTR
  } InstrumentOrientations;


/*
  Set of distance status categories
*/
typedef enum
  {
    DISTANCE_OKAY_INSTR,
    DISTANCE_TOO_CLOSE_INSTR,
    DISTANCE_TOO_FAR_INSTR
  } InstrumentDistStatus;


/*
  Structure to hold mate pair fragment indices & positions in contig/scaffold
*/
typedef struct
{
  CDS_CID_t   fragIID;
  cds_float32 fragOffset5p;
  CDS_CID_t   fragChunkIID;
  CDS_CID_t   mateIID;
  cds_float32 mateOffset5p;
  CDS_CID_t   mateChunkIID;
  CDS_CID_t   libIID;
  FragType    type;
#ifdef TRACK_3P
  cds_float32 fragOffset3p;
  cds_float32 mateOffset3p;
#endif
} MateDetail;
VA_DEF(MateDetail);


/*
  Structure to hold fragment & position
*/
typedef struct
{
  CDS_CID_t  iid;
  FragType    type;
  cds_float32 offset5p;
} FragDetail;
VA_DEF(FragDetail);

typedef struct
{
  VA_TYPE(MateDetail) * happy[NUM_ORIENTATIONS_INSTR];
  VA_TYPE(MateDetail) * misoriented[NUM_ORIENTATIONS_INSTR][NUM_ORIENTATIONS_INSTR];
  VA_TYPE(MateDetail) * misseparatedClose[NUM_ORIENTATIONS_INSTR];
  VA_TYPE(MateDetail) * misseparatedFar[NUM_ORIENTATIONS_INSTR];
  VA_TYPE(FragDetail) * inter;
} MateStatusPositions;

/*
  Structure to hold variable arrays of MateDetails to keep track of
  how many & where good/bad mate pairs are
*/
typedef struct
{
  MateStatusPositions * intra;
  MateStatusPositions * inter;
} MateStatusPositionsSet;


typedef struct
{
  cds_int32 reads;
  cds_int32 externalReads;
  cds_int32 externalFrags;
} MateStats;


typedef struct
{
  MateStats happy;
  MateStats misoriented;
  MateStats misseparatedClose;
  MateStats misseparatedFar;
  MateStats inter;
} MateStatsSet;

/*
  Structure to hold mate instrumentation
*/
typedef struct
{
  cds_uint32 options;
  MateStatusPositionsSet * mateStatus;
  MateStatsSet intra;
  MateStatsSet inter;
  MateStats mateless;
  VA_TYPE(FragDetail) * noMate;
} MateInstrumenter;

/*
  Structure to keep track of a position of a fragment in a surrogate in
  a contig or scaffold
*/
typedef struct
{
  CDS_CID_t contig;
  cds_float32 offset5p;
  cds_float32 offset3p;
  CDS_CID_t nextIndex;
} SurrogateFragLocation;


/*
  Types of breakpoints
*/
typedef enum
  {
    INST_BP_TOO_CLOSE,
    INST_BP_TOO_FAR,
    INST_BP_NORMAL,
    INST_BP_ANTINORMAL,
    INST_BP_OUTTIE,
    INST_BP_INNIE,
    INST_BP_UNKNOWN
  } InstBreakpointType;


/*
  Sections of breakpoint intervals
*/
typedef enum
  {
    BP_ALL,
    BP_LEFT,
    BP_MIDDLE,
    BP_RIGHT
  } BreakpointSection;

/*
  Structure to identify a breakpoint
  NOTE: orientations are deducible, but kept explicitly for convenience
*/
typedef struct
{
  CDS_CID_t iid;
  BreakpointSection section;
  CDS_CID_t contig1;
  CDS_COORD_t end1;
  //  FragOrient orient1;
  CDS_CID_t contig2;
  CDS_COORD_t end2;
  // because end2 may be to the left of the rightmost contig
  CDS_COORD_t contigEnd;
  //  FragOrient orient2;
  InstBreakpointType type;
  VA_TYPE(MateDetail *) mates;
  cds_int32 pairs;
} InstrumenterBreakpoint;
VA_DEF(InstrumenterBreakpoint);

#define INST_MIN_BREAK_MATES   2

/*
  Structure to consolidate mate pairs between a contig pair
  to indicate how things are & how they would be preferred
*/
typedef struct
{
  cds_int32   numPairs;
  cds_float32 distPref;
  ChunkOrientationType orientPref;
} InstrumenterContigPair;
VA_DEF(InstrumenterContigPair);

typedef struct
{
  CDS_CID_t id;
  cds_float32 offset;
  cds_float32 length;
  FragOrient orient;
} ContigPlacement;
VA_DEF(ContigPlacement);

typedef struct
{
  CDS_CID_t        contig1;
  CDS_COORD_t        size1;
  CDS_CID_t        contig2;
  cds_float32      dist;
  ChunkOrientationType orient;
} CP_Index;
VA_DEF(CP_Index);
  
/*
  Structure to keep track of many fragments in surrogates
*/
typedef struct
{
  // bookkeeping for the multiple positions of fragments in surrogate unitigs
  HashTable_AS * surrogateFragHT;
  // but DO NOT reallocate, because ht points to elements of locs
  cds_int32 numAllocatedLocs;
  // number of loc elements in use - to keep track while adding on
  cds_int32 numUsedLocs;
  // fixed size array of surrogate fragment locations
  // estimate (number of fragments in surrogates) * (instances of each surrogate)
  SurrogateFragLocation * surrogateFragLocs;
} SurrogateTracker;


/*
  Structure to keep track of fragments & locales in contig/scaffold
*/
typedef struct
{
  // bookkeeping items for mate relationships
  // for looking up & adding mates in this entity - place all frags here
  HashTable_AS * fragHT;

  // since iterating through the above hash table gets slow,
  // also keep the index in an array
  VA_TYPE(CDS_CID_t) * fragArray;
  
  // for those frags whose mates aren't in this entity
  // to be populated after populating & looping over fragHT
  VA_TYPE(MateDetail) * wExtMates;

} InstrumenterBookkeeping;


// structure for instrumenting a unitig
typedef struct
{
  // IID of unitig
  CDS_CID_t id;

  // is it a surrogate?
  int isSurrogate;
  
  // bit vector of instrumenting option switches
  cds_uint32 options;
  
  // unitig & above level of aggregation
  CDS_COORD_t leftEnd;
  CDS_COORD_t rightEnd;
  int orientation;

  // simple counts of (some) fragment types
  cds_int32 numReads;
  cds_int32 numExtReads;
  cds_int32 numExtFrags;

  // mate status & positions for each library & no-mates list
  MateInstrumenter mates;

  // bookkeeping of fragments, fragment locations, mates, & locales
  InstrumenterBookkeeping bookkeeping;
  
  // intra-contig breakpoints
  // NOTE: just a hook for now
  VA_TYPE(InstrumenterBreakpoint) * breakpoints;

  // for use of this structure in functions, also need surrogate frag data
  // see ScaffoldGraphInstrumenter structure
} UnitigInstrumenter;


// structure for instrumenting a contig
typedef struct
{
  // IID of contig
  CDS_CID_t id;

  // bit vector of instrumenting option switches
  cds_uint32 options;
  
  // contig & above level of aggregation
  CDS_COORD_t leftEnd;
  CDS_COORD_t rightEnd;
  int orientation;

  // unitig counting/sizing
  VA_TYPE(cds_float32) * unitigSizes;
  InstrumenterStatistics unitigSizeStats;

  // surrogate counting/sizing
  VA_TYPE(cds_float32) * surrogateSizes;
  InstrumenterStatistics surrogateSizeStats;

  // simple counts of (some) fragment types
  cds_int32 numReads;
  cds_int32 numExtReads;
  cds_int32 numExtFrags;

  // for instrumenting each unitig
  UnitigInstrumenter reusableUI;

  // for aggregating unitig-level data
  UnitigInstrumenter unitig;
  
  // mate status & positions for each library & no-mates list
  MateInstrumenter mates;

  // bookkeeping of fragments, fragment locations, mates, & locales
  InstrumenterBookkeeping bookkeeping;
  
  // intra-contig breakpoints
  // NOTE: just a hook for now
  VA_TYPE(InstrumenterBreakpoint) * breakpoints;

  // for use of this structure in functions, also need surrogate frag data
  // see ScaffoldGraphInstrumenter structure
} ContigInstrumenter;


// structure for scaffold instrumentation
typedef struct
{
  // IID of scaffold
  CDS_CID_t id;
  
  // bit vector of instrumenting option switches
  cds_uint32 options;
  
  // scaffold & above level of aggregation
  cds_float32 size;
  
  // gap counting & sizes
  VA_TYPE(cds_float32) * scaffoldGapSizes;
  InstrumenterStatistics scaffoldGapSizeStats;

  // inferred edge stddevs to detect scaffold zippering problem
  // these will tend to increase/decrease
  VA_TYPE(cds_float32) * inferredEdgeStddevs;

  // contig counting & sizes
  VA_TYPE(cds_float32) * contigSizes;
  InstrumenterStatistics contigSizeStats;

  /*
  // for aggregating unitig-level data
  UnitigInstrumenter unitig;
  */
  
  // for instrumenting each contig
  ContigInstrumenter reusableCI;
  
  // for aggregating contig-level data
  ContigInstrumenter contig;

  // mate status & positions for each library & no-mates list
  MateInstrumenter mates;

  // bookkeeping of fragments, fragment locations, mates, & locales
  InstrumenterBookkeeping bookkeeping;
  
  // inter-contig breakpoints
  VA_TYPE(InstrumenterBreakpoint) * breakpoints;

  // for finding mates in surrogates within a scaffold
  SurrogateTracker surrogateTracker;

  // to support InstrumentScaffoldMesg implementation
  // hashtable to retrieve contig indices by ID
  HashTable_AS * cpHT;
  // array of contigs with their orientation & location in scaffold
  VA_TYPE(ContigPlacement) * cpArray;

  // hashtable to keep track of contigs that are anchored
  HashTable_AS * anchoredHT;
  // variable array of IntContigPairs for use in IntScaffoldMesg
  // VA_TYPE(IntContigPairs) * icps;
  /*
  // for exploring contig pair issues
  VA_TYPE(InstrumenterContigPair) * contigPairs;
  VA_TYPE(CP_Index) * cp_index;
  */
} ScaffoldInstrumenter;


// structure for scaffold graph instrumentation
typedef struct
{
  // bit vector of instrumenting option switches
  cds_uint32 options;
  
  // scaffold graph level of aggregation
  cds_int32 numFragments;
  cds_int32 numNotInUnitigs;
  cds_int32 numNotInContigs;
  cds_int32 numNotPlaced;
  cds_int32 numChaff;
  cds_int32 numInUnresolvedChunks;
  
  // singleton scaffold counting & sizes (single contig, multiple unitig)
  VA_TYPE(cds_float32) * singletonScaffoldSizes;
  InstrumenterStatistics singletonScaffoldSizeStats;
  VA_TYPE(cds_int32) * unitigsPerSingletonScaffold;
  InstrumenterStatistics unitigsPerSingletonStats;

  // degenerate scaffold counting & sizes (single contig, single unitig)
  VA_TYPE(cds_float32) * degenerateScaffoldSizes;
  InstrumenterStatistics degenerateScaffoldSizeStats;
  cds_int32 numDegenerateScaffoldsWithoutReads;

  // scaffold counting & sizes
  VA_TYPE(cds_float32) * scaffoldSizes;
  InstrumenterStatistics scaffoldSizeStats;

  // for aggregating scaffold-level data
  ScaffoldInstrumenter scaffold;
  
  // bookkeeping of fragments, fragment locations, mates, & locales
  InstrumenterBookkeeping bookkeeping;

} ScaffoldGraphInstrumenter;

int DoSimpleScaffoldChecks(FILE * fp,
                           ScaffoldGraphT * graph,
                           CIScaffoldT * scaff);
int DoSimpleScaffoldGraphChecks(ScaffoldGraphT * graph, char * filename);

void DestroyMateInstrumenter(MateInstrumenter * mi);
int InitializeMateInstrumenter(ScaffoldGraphT * graph,
                               MateInstrumenter * mi);
void CreateMateInstrumenterFromScaffoldGraphInstrumenter(
                                                         MateInstrumenter * mi,
                                                         ScaffoldGraphInstrumenter *sgi);
void CreateMateInstrumenterFromScaffoldInstrumenter(MateInstrumenter * mi,
                                                    ScaffoldInstrumenter * si);
void AddMateInstrumenterCounts(MateInstrumenter * dest,
                               MateInstrumenter * src);
InstrumenterResult CompareMateInstrumenters(MateInstrumenter * miBefore,
                                            MateInstrumenter * miAfter,
                                            InstrumenterVerbosity verbose,
                                            FILE * printTo);
void ResetMateInstrumenterCounts(MateInstrumenter * mi);

void DestroyContigInstrumenter(ContigInstrumenter * ci);
ContigInstrumenter * CreateContigInstrumenter(ScaffoldGraphT * graph,
                                              cds_uint32 options);
int InitializeContigInstrumenter(ScaffoldGraphT * graph,
                                 ContigInstrumenter * ci);
int InstrumentContig(ScaffoldGraphT * graph,
		     HashTable_AS *cpHT,
                     SurrogateTracker * st,
                     ChunkInstanceT * contig,
                     ContigInstrumenter * ci,
                     cds_float32 aEnd,
                     cds_float32 bEnd);
void ComputeContigInstrumenterStats(ScaffoldGraphT * graph,
                                    ContigInstrumenter * ci);
void PrintContigInstrumenter(ScaffoldGraphT * graph,
                             ContigInstrumenter * ci,
                             InstrumenterVerbosity verbose,
                             char * prefix,
                             FILE * printTo);

void PrintScaffoldInstrumenterMateDetails(ScaffoldInstrumenter * si,
                                          FILE * printTo,
					  int printType);
void PrintExternalMateDetailsAndDists(ScaffoldGraphT * graph,
				      VA_TYPE(MateDetail) * mda,
				      char * prefix,
				      FILE * printTo,
				      int printtype);
void PrintUnmatedDetails(ScaffoldInstrumenter * si,
			 FILE * printTo,
			 int printType);


void DestroyScaffoldInstrumenter(ScaffoldInstrumenter * si);
ScaffoldInstrumenter * CreateScaffoldInstrumenter(ScaffoldGraphT * graph,
                                                  cds_uint32 options);
int InitializeScaffoldInstrumenter(ScaffoldGraphT * graph,
                                   ScaffoldInstrumenter * si);
int InstrumentScaffold(ScaffoldGraphT * graph,
                       CIScaffoldT * scaffold,
                       ScaffoldInstrumenter * si,
                       InstrumenterVerbosity verbose,
                       FILE * printTo);
int InstrumentScaffoldPair(ScaffoldGraphT * graph,
                           SEdgeT * sEdge,
                           ScaffoldInstrumenter * si,
                           InstrumenterVerbosity verbose,
                           FILE * printTo);
void ComputeScaffoldInstrumenterStats(ScaffoldGraphT * graph,
                                      ScaffoldInstrumenter * si);
void PrintScaffoldInstrumenter(ScaffoldGraphT * graph,
                               ScaffoldInstrumenter * si,
                               InstrumenterVerbosity verbose,
                               char * prefix,
                               FILE * printTo);

void DestroyScaffoldGraphInstrumenter(ScaffoldGraphInstrumenter * sgi);
ScaffoldGraphInstrumenter * CreateScaffoldGraphInstrumenter(
                                                            ScaffoldGraphT * graph,
                                                            cds_uint32 options);
int InitializeScaffoldGraphInstrumenter(ScaffoldGraphT * graph,
                                        ScaffoldGraphInstrumenter * sgi);
int InstrumentScaffoldGraph(ScaffoldGraphT * graph,
                            ScaffoldGraphInstrumenter * sgi,
                            CDS_COORD_t lowerLimit,
                            CDS_COORD_t upperLimit,
                            InstrumenterVerbosity verbose,
                            FILE * printTo);
void ComputeScaffoldGraphInstrumenterStats(ScaffoldGraphT * graph,
                                           ScaffoldGraphInstrumenter * sgi);
void PrintScaffoldGraphInstrumenter(ScaffoldGraphT * graph,
                                    ScaffoldGraphInstrumenter * sgi,
                                    InstrumenterVerbosity verbose,
                                    FILE * printTo);

int InstrumentContigEnd(ScaffoldGraphT * graph,
                        ScaffoldInstrumenter * si,
                        ChunkInstanceT * thisCI,
                        int end);
int InstrumentContigEndPartial(ScaffoldGraphT * graph,
                               ScaffoldInstrumenter * si,
                               ChunkInstanceT * thisCI,
                               int end,
                               cds_int32 numContigs);
void GetMateInstrumenterFromScaffoldGraphInstrumenter(
                                                      MateInstrumenter * mi,
                                                      ScaffoldGraphInstrumenter * sgi);
void GetMateInstrumenterFromScaffoldInstrumenter(MateInstrumenter * mi,
                                                 ScaffoldInstrumenter * si);
int InstrumentContigPath(ScaffoldGraphT * graph,
                         ScaffoldInstrumenter * si,
                         CDS_CID_t firstID,
                         int firstEnd,
                         CDS_CID_t lastID);

int AdjustCIScaffoldLabels(ScaffoldGraphT * graph,
                           cds_int32 * numScaffoldIDs);

cds_int32 GetMateStatsBad(MateStatsSet * mss);
cds_int32 GetMateStatsHappy(MateStatsSet * mss);
void PrintFragment(CIFragT * frag, CDS_CID_t index, FILE * printTo);
void PrintBreakpoint(InstrumenterBreakpoint * bp,
                     char * prefix,
                     FILE * printTo);
void FindRockStoneUnitigs(ScaffoldGraphT * graph);

#ifdef TRACK_3P
void safelyAppendInstInfo(char **locs,int32 utgIID, int *lenloc, int *lenUsed);
#endif

#endif
