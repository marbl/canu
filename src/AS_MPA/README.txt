README file for /dev4/dewim/mps/bin

Instructions to generate assembly mate data from SNAPPER output

**********************************************************************
SNAPPER output has the form

0[150-0-0] 1[221670143-221670291] <148-0-98-forward-unknown>

where:
0   is the index of the fragment in the input multifasta file used by SNAPPER
150 is the length of the fragment
1   is the index of the chromosome in the other SNAPPER input multifasta file
221670143-221670291 is the interval of the match on the chromosome
148 is the number of basepairs of the fragment in the match
98  is the percent identity
forward is the orientation of the fragment wrt the chromosome
**********************************************************************


**********************************************************************
buildMapStore input has the form

17000178532974 0 194348501 194347841

where:
17000178532974 is the fragment UID
0              is the chromosome index in the multifasta file
194348501      is the coordinate of the 5' fragment end on the chromosome
194347841      is the coordinate of the 3' fragment end on the chromosome
**********************************************************************

To get from SNAPPER output to buildMapStore input you need:

1. a file listing fragment UIDs corresponding to indices in the
   SNAPPER input multifasta file. The format is one index & UID per
   line, like:

0 17000179060752
1 17000179060954
2 17000179060956
3 17000179060958

   see /dev4/dewim/mps/human/stores/matedFragUIDs.txt for an example

2. If there are multiple SNAPPER output files that need to be
   concatenated, modify concatRegions.c to identify filenames & index
   offsets of the fragments - in order. Run concatRegions on the set
   of files to produce a single file.

3. Generate a celagram file of matches per fragment, in order. Run
   matchCounts on the SNAPPER output file to produce this.

4. Generate a file listing uniquely mapped mated frags from the
   SNAPPER output, UIDS file, & celagram file by running
   getUniqueMatchesFromCGM.

5. Generate buildMapStore input by including just fields 1, 3, 4, and
   5 from this file.

Then, run buildMapStore (in cds/AS/src/AS_CVT) to create a new map
store for the assembly.

To get a file of intra-chromosome mate pairs, run dumpMappedMatePairs
(in cds/AS/src/AS_CVT).

To get a file of inter-chromosome mate pairs, run dumpMappedElsewheres
(in cds/AS/src/AS_CVT).
