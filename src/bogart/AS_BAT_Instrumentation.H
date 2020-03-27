
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

#ifndef INCLUDE_AS_BAT_INSTRUMENTATION
#define INCLUDE_AS_BAT_INSTRUMENTATION


void  checkUnitigMembership(TigVector &tigs);
void  reportOverlaps(TigVector &tigs, const char *prefix, const char *name);
void  reportTigs(TigVector &tigs, const char *prefix, const char *name, uint64 genomeSize);

void  classifyTigsAsUnassembled(TigVector    &tigs,
                                uint32        fewReadsNumber,
                                uint32        tooShortLength,
                                double        spanFraction,
                                double        lowcovFraction,   uint32  lowcovDepth);

#endif //  INCLUDE_AS_BAT_INSTRUMENTATION
