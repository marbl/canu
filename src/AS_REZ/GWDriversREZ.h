
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
/**********************************************************************

        Module:  GWDriversREZ.h

   Description:  constants/protos for GWDriversREZ.h
                 notes about gap walking and bcc

    Programmer:  S. Lonardi (stelo@cs.purdue.edu)

       Written:   8 July 99
 **********************************************************************

  Notes: Filling Gaps by "walking" (a.k.a., phase 2 or repeat rez)
 (warning: I use the words "unitig", "chunk" and "CI" interchangeably)

Three strategies has been designed to fill the gaps after the first
phase of repeat resolution.

* Note1: before running -r 4 create a directory ./cam in the directory
* where the cgb file is located. The walker will create a .cam file for
* each gap walked with the name
* <prefix>.<begin_chunk_id>.<end_chunk_id>.ps.cam and another .cam file
* for each gap NOT walked with the same naming convention but with an
* underscore as a prefix.

* Note2 (08.10.99): The current "phase 1" repeat resolution is so
* aggressive that not much is left on small testcases like a004 and
* smaller (thanks, Art :-).  If you want to use them anyway for speed
* of testing you can comment the line "// if (level == 2)" in
* GapFillREZ.c The effect of this modification is that phase one will
* insert less chunks in the scaffolds, and as a consequence, it will
* leave the gap-walker more work to do

1) (cgw -r 4 ... )

The "simple gap walking" consists of selecting the subgraph induced by
a path walking using links (NOT the overlaps) between the two u-unitgs
that flank the gap.  The "walking" is a simple DFS with some
constraints.

The main constraint is that we do not go too far from the starting
u-unitigs. It is achieved by having an estimate of the gap size.  The
gap size can be computed from the offsets in the scaffold, if the gap
is internal, or from a scaffold edge, if the gap is in-between
scaffolds.  While traversing the graph we keep track of the mean and
variance of links and chunks visited.  When the
three-standard-deviation interval around the mean of the traveled
distance has empty intersection (and is over) the interval of the
estimated gap size, we break the traversal.

If we find a path that reaches the other u-unitig, we mark all the
chunks that lay on that path. A new DFS is started from the original
u-unitig that could possibly find a alternative paths (a path that
reaches a previously path-marked node is considered another good
path). Another set of DFSes are started from the other side of the
gap.

All the chunks that has been marked will become the subgraph of the
chunks to place in the gap.

This simple strategy gives good results on artificial data sets.  For
the test-cases a005, a004, a003, all the chunks selected actually
belongs to the gap.  However, some gaps cannot be walked because not
enough links exist. We plan to use other strategies for those cases.


                               USING ALL EDGES - artificial dataset
        +--------------+------------+-------------------------------+----------+------------
        | gaps > 1Kbps |   walked   |   subgraphs that contains     |  walked  | non walked
        |              |            | at least 1 CI outside the gap | gap size |  gap size
 -------+--------------+------------+-------------------------------+----------+------------
   a005 |       5      |       5    |               2               |   2.5K   |     N/A
   a004 |      44      |      44    |              25               |   3.2K   |     N/A
   a003 |     719      |     717#   |             552               |   2.6K   |    2.6K
 -------+--------------+------------+-------------------------------+----------+------------

 # a correction to get those two gaps that we are missing has been
   tried (set ALLOW_RECOMPUTE_SHORTER_PATHS to 1 and cross the
   fingers: for some reasons the dfs loops ... happy debugging :-) The
   problem comes from the fact that we mark nodes when we visit during
   the traversal but some traversals of the graph could end up pruned
   because of the gap estimate. However, it could happen that later we
   find another path that is "shorter", therefore compatible with the
   gap size, that wants to use that node already used.

                               USING ONLY LINKS - artificial dataset
        +--------------+------------+-------------------------------+----------+------------
        | gaps > 1Kbps |   walked   |   subgraphs that contains     |  walked  | non walked
        |              |            | at least 1 CI outside the gap | gap size |  gap size
 -------+--------------+------------+-------------------------------+----------+------------
   a005 |       5      |       2    |               0               |   2.5K   |    2.6K
   a004 |      44      |      27    |              1 ur             |   2.2K   |    5.6K
   a003 |     719      |     445    | 1 due to rr + 1 to other comb |   2.5K   |    2.6K
 -------+--------------+------------+-------------------------------+----------+------------


                               USING ALL EDGES
          +--------------+------------+-------------------------------+----------+------------
          | gaps > 1Kbps |   walked   |   subgraphs that contains     |  walked  | non walked
          |              |            | at least 1 CI outside the gap | gap size |  gap size
  --------+--------------+------------+-------------------------------+----------+------------
  hflu001 |       7      |       7    |              7                |    5K    |    N/A
ww001_1_1 |    1931      |    1919    |            831                |   3.2K   |    2K
  --------+--------------+------------+-------------------------------+----------+------------

                               USING ONLY LINKS
          +--------------+------------+-------------------------------+----------+------------
          | gaps > 1Kbps |   walked   |   subgraphs that contains     |  walked  | non walked
          |              |            | at least 1 CI outside the gap | gap size |  gap size
  --------+--------------+------------+-------------------------------+----------+------------
  hflu001 |       7      |       7    |             5 **              |    5K    |    4.9K
ww001_1_1 |    1931      |     765    |           100 ***             |   3.1K   |    3.1K  
  --------+--------------+------------+-------------------------------+----------+------------


   ** For 5 subgraphs of hflu001 out of 7 we observed at least one
   chunks from a different place in the genome. All the five times the
   chunks wrongly selected are only "rr".  Using the threshold -10.0
   on the coverage (see below) we can walk only two path on hflu001,
   and both are probably correct ("probably" because the path contains
   a "rr" chunk).

   *** With the test-case ww001_1_1 we made several mistakes in the
   *** selection. However, if we filter out chunks that are not "uu"
   *** before walking (using 99724 chunks) we are able to walk 653
   *** gaps with no (0) mistakes.  Unfortunately, we cannot assume to
   *** have that annotation for real data.

  Ideas to improve gap-walking:

  a - using path confirmed links (set ANNOTATE_CONFIRMING_PATH to 1 or re-use
      the code that Granger wrote in TransitiveReduction_CGW.c) My preliminary testings
      suggest that this kind of preprocessing is heavy in terms of computations and
      seems to help not too much, but maybe my code is not working ok

  b - use coverage statistics (implemented). Filtering out chunks that have "very low" coverage
      statistic is very easy and effective as the following table explains

                          USING ONLY LINKS on ww001_1_1
 ---------------------+--------------+-------------+-------------+-------------+-------------+
   repeat contents of |              |             |             |             |             |
   the bad subgraph   | no filtering | cov > -10.0 |  cov > -5.0 | cov > -2.0  |  cov > 0.0  |
 ---------------------+--------------+-------------+-------------+-------------+-------------+
       uu + rr        |      73(10)  |      27(15) |      17(?)  |       6(2)  |       0     |
       uu + ur        |       7(6)   |       7     |       8     |       5     |       0     |
       uu + ru        |       1(0)   |       1     |       1     |       2     |       0     |
 uu + any other combo |      19      |      12     |       9     |       3     |       0     |
 ---------------------+--------------+-------------+-------------+-------------+-------------+
                      |     100      |      47     |      35     |      16     |       0     |  total "bad" 
                      +--------------+-------------+-------------+-------------+-------------+
                      |     765      |     720     |     709     |     647     |     448     |  gaps walked
                      +--------------+-------------+-------------+-------------+-------------+
                      |  107176      |  106747     |  106213     |  103210     |   17152     |  total chunks 
                      +--------------+-------------+-------------+-------------+-------------+


  c - check that we did not use the same chunk in multiple paths. We detected 12 multiple
      uses of the same chunk (always a "rr") in the case of coverage > -10.0
      and 4 in the case of coverage > -2.0 (see the number in parenthesis is the number of
      "bad" gaps that in priciple we still have after detecting those multiple uses)
      (the detection has been implemented, not it's not clear what to do after we detect multiple uses:
      discard all the uses of chunk? discard all the paths that contain that chunk? The correction "phase"
      has not been implemented yet)

  d - start the walking not only from the two flanking chunks, but also on all the unique chunks
      that overlap the two flanking chunks; a preliminary result in a004 shows an increase of 25%
      in the number of gaps walked.

  e - as for 3), create a subgraph of all the chunks that has a direct link to 
      any unique in the same scaffold and try to walk (using overlap, too) in that subgraph. 
      This strategy is similar to approach 3) below




2) (cgw -r 5 ... )


A way to get some of the chunks in the gap, when there is not enough
edges to walk on the other side is to use biconnected components.  We
can build the biconnected components of the chunk graph using only
links (not overlaps), and then look at the components that contains at
least one u-unitig that flank the gap. The biconnectedness guarantees
that each chunk in the component has at least two paths to the
u-unitig. We think it should be enough evidence to select some of the
chunks on one side of the gap. This procedure should be iterated.

BTW, this strategy is part of a more general strategy that does not
necessarily involve gaps:

- build the bcc on the chunk graph  with edges of  weight 2 or more or
  with edge of any weight that are not pure overlap (done)
- scan the components that are not completely formed by u-unitigs or by
  non-u-unitigs (done Scan_Components() in BccREZ.c)
- place the non-u-unitigs with respect to the u-unitig in the same
  component (to be done, Compute_Tentative_Position_Bcc())





3) The following is not currently implemented.  We initially select
chunks that have links to the unique chunks in the scaffold. The
chunks selected constitute a sort of "stepping-stones" of the gap and
should be correctly placed with very high probability (to be
tested!). Then we try to find a path of overlapping chunks that is
consistent with the stones and with the estimate of the gap size.


TO DO:

- collect more data about the strategies (a,b,c,d,e) to improve
  walking performance (i.e., getting more chunks/coverage, walking
  more paths, making less errors). Some addition code could be
  required

- write routine(s) to compute a tentative position for the
  chunks selected by gap walking (Compute_Tentative_Position())
  (Mike?)

- write routine(s) to compute a tentative position for the
  chunks selected by the bcc (Compute_Tentative_Position_Bcc(),
  similar to the previous?)

- implement strategy 3) (Karin?)

- interscaffold walking? (see Test_Walker1() and Inter_Scaffold_Gap_Walker())
  maybe it will be not necessary: we hope to get one big scaffold, so every gap
  will be internal

 **********************************************************************/

/*********************************************************************
   CVS_ID: $Id: GWDriversREZ.h,v 1.3 2005-03-22 19:07:37 jason_miller Exp $
 *********************************************************************/

#ifndef GWDRIVERS_REZ_H
#define GWDRIVERS_REZ_H

//
// constants for testing the GapWalker (to be removed at some point)
//
#define WALK_BEGIN_POS     1090000   
// position for testing
#define WALK_END_POS       1170000   
// "
#define WALK_BEGIN_S_ID          8   
// left scaffold (to test the walker)
#define WALK_END_S_ID            6   
// right scaffold (to test the walker)
#define WALK_BEGIN_CID          50   
// begin chunk id
#define WALK_END_CID            60   
// end chunk id
#define ANNOTATE_CONFIRMING_PATH 0   
// set this to 1 to run the path confirmation on CIedges
#define CHECK_MULTIPLE_USE       0   
// set this to 1 to check if chunk happen to be used in multiple paths

#define LOW_BAYESIAN_THRESHOLD  0.3
#define HIGH_BAYESIAN_THRESHOLD 1.0
#define SWITCH_THRESHOLD 2000

#define BAC_WALKING_THRESHOLD  0.3

#define CHECKPOINT_DISTANCE 10000000    
// this is the number of basepairs we have inspected in scaffolds before
// gapwalking issues another checkpoint


typedef  struct
{
  int  scaff_id;
  int  scaff_length;
}  Scaffold_Size_t;

int Compare_Scaff_Size(const void *e1, const void *e2);
// 
//
// protos
//

//
// main Walker function. it calls the appropriate function given the
// <level> as input.  The debug information are stored in "*.gwlog"
//
// The levels > 5 are experimental
//

//
// the main biconnected components function
//
// it will returns the number of chunks inserted
// (when the implementation is finished)
//
int Biconnected_Components(float (*)(CIEdgeT *, CDS_CID_t));


//
// the main intrascaffold gap walking
//
// void Intra_Scaffold_Path_Finding(float (*)(CIEdgeT *, CDS_CID_t), float qt, int startWalkFrom, int bac_walking);
void Intra_Scaffold_Path_Finding( int startWalkFrom,
                                  double gapSizeStdDevs, 
                                  int doSmallContigRemoval,
                                  int trimScaffoldEnds);

//
// the main interscaffold gap walking
//
void Inter_Scaffold_Path_Finding(void);

void Inter_Scaffold_Analysis();

// a helper function for walking
void ComputeGapStatisitics(void);

//
// Compute the Dijkstra's shortest path: a_cid is the source, b_cid is the
// sink. The subgraph is obtained selecting the chunks in the
// interval [begin_pos,end_pos]
// 
void Test_Walker0(int32, int32, CDS_CID_t, CDS_CID_t,
		  int (*)(ChunkInstanceT *, chunk_subgraph *, int32, int32),
		  float (*)(CIEdgeT *, CDS_CID_t));

//
// Test the gap walker between two scaffold
//
// note: Inter_Scaffold_Gap_Walker() is not done yet
//
void Test_Walker1(CDS_CID_t, CDS_CID_t,
		  float (*)(CIEdgeT *, CDS_CID_t));

//
// Test the walker selecting all the chunks between <begin_pos> and
// <end_pos> and trying to find a path
//
void Test_Walker2(int32, int32, CDS_CID_t, CDS_CID_t, int,
		  int (*)(ChunkInstanceT *, chunk_subgraph *, int32, int32),
		  float (*)(CIEdgeT *, CDS_CID_t));

//
// collect various statistical measures on the chunks,
// like outdegree, number of unique chunk reachable
//
void Test_Walker3(int32, int32,
		  int (*)(ChunkInstanceT *, chunk_subgraph *, int32, int32),
		  float (*)(CIEdgeT *, CDS_CID_t));

//
// Compute coordinates for chunks found walking across gaps 
void Compute_Tentative_Position(chunk_subgraph *, CDS_CID_t, CDS_CID_t);

// 
// This function update the scaffolds by filling the appropriate
// gaps as indicated by the infos in the gapAssignment (only the
// ones which has a "keep" flag TRUE will be considered). if the flag
// RebuildSEdges is TRUE then we will rebuild SEdges. if the flag
// ChiSquare is TRUE then we will compute the ChiSquare test at
// each scaffold. returns the number of inserted chunks
int Update_Scaffold_Graph
    (ScaffoldGraphT *, Scaffold_Fill_t *, int, int, int,
     int copyAllOverlaps, int, Kind_Of_Fill_t);

LengthT Compute_Gap_Length(ChunkInstanceT *, ChunkInstanceT *, int);

// can this be merged with localeInfo struct in FbacREZ.h?
typedef struct ScaffoldEnd
{
	  int localeNumber;
	  int spanningBACs;
	  int leftScaffID;
	  int rightScaffID;
	  int leftScaffEnd;  // A or B end of scaffold
	  int rightScaffEnd;
	  CDS_CID_t leftExtremalFragID;
	  CDS_CID_t rightExtremalFragID;
	  CDS_CID_t leftExtremalFragIid;
	  CDS_CID_t rightExtremalFragIid;
	  int leftContigID;
	  int rightContigID;
	  int leftBasesFromEnd;
	  int rightBasesFromEnd;
	  int leftBacOffset;
	  int rightBacOffset;
	  int leftInSurrogate;   // currently not allowing ScaffoldEnds based on frags in surrogates
	  int rightInSurrogate;  // so these fields could be eliminated
	  int onWalkList;
	  int leftIidTrend;   // tracks whether the locale iids are increasing or decreasing
	  int rightIidTrend;   // tracks whether the locale iids are increasing or decreasing
} ScaffoldEndT;

CIScaffoldT* CreateNewScaffold();

#endif
