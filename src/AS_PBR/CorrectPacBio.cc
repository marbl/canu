/*
Copyright (C) 2011, Battelle National Biodefense Institute (BNBI);
all rights reserved. Authored by: Sergey Koren

This Software was prepared for the Department of Homeland Security
(DHS) by the Battelle National Biodefense Institute, LLC (BNBI) as
part of contract HSHQDC-07-C-00020 to manage and operate the National
Biodefense Analysis and Countermeasures Center (NBACC), a Federally
Funded Research and Development Center.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

 * Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

 * Neither the name of the Battelle National Biodefense Institute nor
  the names of its contributors may be used to endorse or promote
  products derived from this software without specific prior written
  permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

const char *mainid = "$Id: CorrectPacBio.cc,v 1.16 2012-08-20 13:10:37 skoren Exp $";

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_gkpStore.h"

#include "AS_PBR_util.hh"
#include "AS_PBR_correct.hh"
#include "AS_PBR_filter.hh"
#include "AS_PBR_mates.hh"
#include "AS_PBR_output.hh"

#include <map>
#include <pthread.h>

static const char* EXECUTABLE_NAME="correctPacBio";

using namespace std;

static uint32 loadFragments(gkStream *fs, uint32* includeLib, map<AS_IID, uint32>& frgToLen, map<AS_IID, uint8>& frgToLib, map<AS_IID, AS_IID>& frgToMate) {
    gkFragment  fr;
    uint32 counter = 0;

    fprintf(stderr, "Streaming fragments\n");
    // figure out which libraries we want to use and store fragment info
    while (fs->next(&fr)) {
        int32 len = fr.gkFragment_getClearRegionLength();
        AS_IID libID = fr.gkFragment_getLibraryIID();

        if (includeLib[libID] == TRUE) {
            counter++;
        }
        frgToLib[fr.gkFragment_getReadIID()] = libID;
        frgToLen[fr.gkFragment_getReadIID()] = len;
        frgToMate[fr.gkFragment_getReadIID()] = fr.gkFragment_getMateIID();
    }
    return counter;
}

// partition the work
static AS_IID partitionWork( uint32 counter, uint32 *libToInclude, map<AS_IID, uint8>& frgToLib, int &numThreads, int &partitions, uint32& perFile, PBRThreadWorkArea *wa) {
    uint32 lastEnd = 0;
    uint32 lastFile = 0;
    uint32 currThread = 0;
    AS_IID lastFrag = 0;
    AS_IID firstFrag = 0;
    PBRThreadGlobals* waGlobal = wa[0].globals;

    if (partitions > counter) {
        partitions = counter;
    }
    if (numThreads > counter) {
        numThreads = counter;
    }

    uint32 perThread = (uint32) floor((double)counter / numThreads);
    perFile   = (uint32) round((double)counter / partitions);
    fprintf(stderr, "Each thread responsible for %d (%d threads) and each file %d of %d fragments\n", perThread, numThreads, perFile, counter);

    counter = 0;

    for(map<AS_IID, uint8>::const_iterator iter = frgToLib.begin(); iter != frgToLib.end(); iter++) {
        if (libToInclude[iter->second] == TRUE) {
            if (counter == 0) { firstFrag = iter->first; lastEnd = iter->first; }
            counter++;
        }
        if (currThread < waGlobal->numThreads -1 ) {
            if ((counter > 0 && counter % perThread == 0)) {
                wa[currThread].start = lastEnd;
                wa[currThread].end = iter->first-1;
                wa[currThread].fileStart = (uint32) ceil((double)(wa[currThread].start-wa[0].start+1) / perFile);
                wa[currThread].fileEnd = (uint32) floor((double)(wa[currThread].end-wa[0].start+1) / perFile);
                if (wa[currThread].fileStart <= lastFile) { wa[currThread].fileStart = lastFile + 1; }
                if (wa[currThread].fileEnd < wa[currThread].fileStart) { wa[currThread].fileEnd = wa[currThread].fileStart; }
                lastFile = MAX(wa[currThread].fileStart, wa[currThread].fileEnd);
                wa[currThread].id = currThread;

                lastEnd = iter->first;
                currThread++;
            }
        }
        lastFrag = (iter->first > lastFrag ? iter->first : lastFrag);
    }
    // give the rest to the last thread
    wa[currThread].start = lastEnd;
    wa[currThread].end = lastFrag;
    wa[currThread].fileStart = (uint32) ceil((double)(wa[currThread].start-wa[0].start+1) / perFile);
    wa[currThread].fileEnd = MIN(partitions, (uint32) floor((double)(wa[currThread].end-wa[0].start+1) / perFile));
    if (wa[currThread].fileStart <= lastFile) { wa[currThread].fileStart = lastFile + 1; }
    if (wa[currThread].fileEnd < wa[currThread].fileStart) { wa[currThread].fileEnd = MIN(partitions, wa[currThread].fileStart); }

    wa[currThread].id = currThread;

    return firstFrag;
}

int
main (int argc, char * argv []) {
    fprintf(stderr, "Parsing arguments\n");

    PBRThreadGlobals  thread_globals;

    char      *gkpStorePath           = NULL;
    thread_globals.ovlStoreUniqPath   = NULL;

    // initialize default parameters
    thread_globals.numThreads        = 2;
    thread_globals.maxErate          = 0.25;
    thread_globals.erate             = 0.15;
    thread_globals.elimit            = 4.5;
    thread_globals.globalRepeats     = TRUE;
    thread_globals.repeatMultiplier  = 2.0;
    thread_globals.partitions        = 100;
    thread_globals.minLength         = 500;
    thread_globals.fixedMemory       = FALSE;
    thread_globals.maxUncorrectedGap = 0;
    thread_globals.hasMates			= false;
    thread_globals.percentToEstimateInserts = DEFAULT_SAMPLE_SIZE;
    thread_globals.percentShortReadsToStore = DEFAULT_SHORT_READ_STORE_SIZE;
    thread_globals.verboseLevel 			  = VERBOSE_DEBUG;

    strcpy(thread_globals.prefix, "asm");
    argv[0][strlen(argv[0])-strlen(EXECUTABLE_NAME)] = '\0';
    strcpy(thread_globals.exePath, argv[0]);
    fprintf(stderr, "Starting %s, running from %s\n", EXECUTABLE_NAME, thread_globals.exePath);

    argc = AS_configure(argc, argv);

    int err = 0;
    int arg = 1;
    while (arg < argc) {
        if (strcmp(argv[arg], "-G") == 0) {
            gkpStorePath = argv[++arg];

        } else if (strcmp(argv[arg], "-O") == 0) {
            if      (thread_globals.ovlStoreUniqPath == NULL)
                thread_globals.ovlStoreUniqPath = argv[++arg];
            else
                err++;

        } else if (strcmp(argv[arg], "-f") == 0) {
            thread_globals.fixedMemory = (atoi(argv[++arg]) != 0 ? TRUE : FALSE);

        } else if (strcmp(argv[arg], "-p") == 0) {
            thread_globals.partitions = atoi(argv[++arg]);
            if (thread_globals.partitions <= 0) { thread_globals.partitions = 1; }

        } else if (strcmp(argv[arg], "-o") == 0) {
            strncpy(thread_globals.prefix, argv[++arg], FILENAME_MAX);

        } else if (strcmp(argv[arg], "-l") == 0) {
            thread_globals.minLength = atoi(argv[++arg]);

        } else if (strcmp(argv[arg], "-e") == 0) {
            thread_globals.erate = atof(argv[++arg]);

        } else if (strcmp(argv[arg], "-E") == 0) {
            thread_globals.elimit = atof(argv[++arg]);

        } else if (strcmp(argv[arg], "-c") == 0) {
            thread_globals.maxErate = atof(argv[++arg]);

        } else if (strcmp(argv[arg], "-R") == 0) {
            thread_globals.globalRepeats = FALSE;
            thread_globals.repeatMultiplier = atof(argv[++arg]);

        } else if (strcmp(argv[arg], "-t") == 0) {
            thread_globals.numThreads = atoi(argv[++arg]);
            if (thread_globals.numThreads <= 0) { thread_globals.numThreads = 1; }

        } else if (strcmp(argv[arg], "-C") == 0){
            int32 cutoff = atoi(argv[++arg]);
            if (cutoff <= 0) { thread_globals.covCutoff = 0; }
            else if (cutoff >= MAX_COV_HIST) {thread_globals.covCutoff = MAX_COV_HIST - 1; }
            else { thread_globals.covCutoff = cutoff; }
            fprintf(stderr, "The cutoff is set to be %d\n", thread_globals.covCutoff);

        } else if (strcmp(argv[arg], "-M") == 0) {
            int32 maxGap = atoi(argv[++arg]);
            if (maxGap <= 0) { thread_globals.maxUncorrectedGap = 0; }
            else { thread_globals.maxUncorrectedGap = maxGap; }

        } else if (strcmp(argv[arg], "-m") == 0) {
            thread_globals.percentShortReadsToStore = atof(argv[++arg]);
            if (thread_globals.percentShortReadsToStore <= 0 || thread_globals.percentShortReadsToStore > 1) { thread_globals.percentShortReadsToStore = DEFAULT_SHORT_READ_STORE_SIZE; }
            fprintf (stderr, "Percent short reads will be %f\n", thread_globals.percentShortReadsToStore);

        } else if (strcmp(argv[arg], "-S") == 0) {
            thread_globals.percentToEstimateInserts = atof(argv[++arg]);
            if (thread_globals.percentToEstimateInserts <= 0 || thread_globals.percentToEstimateInserts > 1) { thread_globals.percentToEstimateInserts = DEFAULT_SAMPLE_SIZE; }

        } else if (strcmp(argv[arg], "-v") == 0) {
            thread_globals.verboseLevel = atoi(argv[++arg]);
            if (thread_globals.verboseLevel < VERBOSE_OFF) { thread_globals.verboseLevel = VERBOSE_OFF; }
            if (thread_globals.verboseLevel > VERBOSE_DEVELOPER) { thread_globals.verboseLevel = VERBOSE_DEVELOPER; }

        } else {
            err++;
        }

        arg++;
    }

    if ((thread_globals.maxErate < 0.0) || (AS_MAX_ERROR_RATE < thread_globals.maxErate))
        err++;
    if ((thread_globals.erate < 0.0) || (AS_MAX_ERROR_RATE < thread_globals.erate))
        err++;
    if (thread_globals.elimit < 0.0)
        err++;
    if (gkpStorePath == NULL)
        err++;
    if (thread_globals.ovlStoreUniqPath == NULL)
        err++;
    if (thread_globals.minLength <= 0)
        err++;

    if (err) {
        fprintf(stderr, "usage: %s -O ovlStore -G gkpStore [options]\n", argv[0]);
        fprintf(stderr, "\n");
        fprintf(stderr, "  -O         Mandatory path to an ovlStore.\n");
        fprintf(stderr, "  -G         Mandatory path to a gkpStore.\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "  -e 0.15   no more than 0.015 fraction (1.5%%) error\n");
        fprintf(stderr, "  -E 0      no more than 0 errors\n");
        fprintf(stderr, "  -c 0.25   ignore overlaps over this rate before correction\n");
        fprintf(stderr, "  -l %d     ignore corrected fragments less than %d bp\n", thread_globals.minLength, thread_globals.minLength);
        fprintf(stderr, "  -t %d     use %d threads to process correction in parallel\n", thread_globals.numThreads, thread_globals.numThreads);
        fprintf(stderr, "  -p %d     output %d results files, corresponds to #of parallel consensus jobs desired\n", thread_globals.partitions, thread_globals.partitions);
        fprintf(stderr, "  -p %s     output prefix of %s\n", thread_globals.prefix, thread_globals.prefix);
        fprintf(stderr, "\n");
        fprintf(stderr, " -C %d 	 Specify the pacBio coverage (integer) instead of automatically estimating.\n", thread_globals.covCutoff);
        fprintf(stderr, " -M %d	 	 The maximum uncorrected PacBio gap that will be allowed. When there is no short-read coverage for a region, by default the pipeline will split a PacBio sequence. This option allows a number of PacBio sequences without short-read coverage to remain. For example, specifying 50, will mean 50bp can have no short-read coverage without splitting the PacBio sequence. Warning: this will allow more sequences that went through the SMRTportal to not be fixed.\n", thread_globals.maxUncorrectedGap);
        fprintf(stderr, " -m %f	 	 The percentage of short reads to use to recruit other PacBio sequences to fill coverage gaps. Must be a decimal value between 0 and 1. Higher values lead to a longer runtime but more gaps being recovered. The default is %f\n", thread_globals.percentShortReadsToStore, DEFAULT_SHORT_READ_STORE_SIZE);
        fprintf(stderr, " -S %f	 	 The percentage of mates used to estimate the library insert size. Must be a decimal value between 0 and 1. For example, specifying 0.10 will use 10 percent of the mates to estimate the insert size. The default is %f\n", thread_globals.percentToEstimateInserts, DEFAULT_SAMPLE_SIZE);
        fprintf(stderr, "\n");
        fprintf(stderr, "  -v %d     level of verbosity. Higher values generate more debugging output\n", thread_globals.verboseLevel);
        fprintf(stderr, "\n");
        fprintf(stderr, "  -R %f     consider a pileup of %f times the mean for a single corrected read to be a repeat and distribute reads to their best locations (this is only useful for metagenomic or non-even coverage datasets. Otherwise the global repeat estimate is used instead)\n", thread_globals.repeatMultiplier, thread_globals.repeatMultiplier);
        fprintf(stderr, "\n");

        if ((thread_globals.maxErate < 0.0) || (AS_MAX_ERROR_RATE < thread_globals.maxErate))
            fprintf(stderr, "Invalid overlap error threshold (-c option); must be between 0.00 and %.2f.\n", AS_MAX_ERROR_RATE);

        if ((thread_globals.erate < 0.0) || (AS_MAX_ERROR_RATE < thread_globals.erate))
            fprintf(stderr, "Invalid overlap error threshold (-e option); must be between 0.00 and %.2f.\n",
                    AS_MAX_ERROR_RATE);

        if (thread_globals.elimit < 0.0)
            fprintf(stderr, "Invalid overlap error limit (-E option); must be above 0.\n");

        if (gkpStorePath == NULL)
            fprintf(stderr, "No gatekeeper store (-G option) supplied.\n");

        if (thread_globals.ovlStoreUniqPath == NULL)
            fprintf(stderr, "No overlap store (-O option) supplied.\n");

        if (thread_globals.minLength <= 0)
            fprintf(stderr, "Invalid min length, must be positive: %d\n", thread_globals.minLength);

        exit(1);
    }

    // set up global data structures
    fprintf(stderr, "Opening stores\n");
    thread_globals.partitionStarts = new pair<AS_IID, AS_IID>[thread_globals.partitions];
    thread_globals.gkp = new gkStore(gkpStorePath, FALSE, FALSE, TRUE);
    uint32 i;
    thread_globals.libToInclude = new uint32[thread_globals.gkp->gkStore_getNumLibraries() + 1];
    memset(thread_globals.libToInclude, 0, (thread_globals.gkp->gkStore_getNumLibraries() + 1) * sizeof(uint32));

    fprintf(stderr, "Loading library information\n");
    for (i=1; i<=thread_globals.gkp->gkStore_getNumLibraries(); i++) {
        gkLibrary      *gkpl = thread_globals.gkp->gkStore_getLibrary(i);
        thread_globals.libToInclude[i] |= gkpl->doConsensusCorrection;
        if (gkpl->doConsensusCorrection) {
            strncpy(thread_globals.libName, gkpl->libraryName, LIBRARY_NAME_SIZE);
        }

        // store the library size as well
        thread_globals.libToOrientation[i] = gkpl->orientation;
        if (gkpl->orientation == AS_READ_ORIENT_UNKNOWN) {
            thread_globals.libToSize[i] = pair<double, double>(0, 0);
        } else {
            thread_globals.libToSize[i] = pair<double, double>(gkpl->mean - CGW_CUTOFF * gkpl->stddev, gkpl->mean + CGW_CUTOFF * gkpl->stddev);
            thread_globals.hasMates = true;
        }
    }

    fprintf(stderr, "Correcting fragments\n");
    pthread_attr_t    attr;
    pthread_t  *      thread_id;
    PBRThreadWorkArea* thread_wa;
    thread_id = new pthread_t[thread_globals.numThreads];
    thread_wa = new PBRThreadWorkArea[thread_globals.numThreads];

    // partition data
    thread_wa[0].globals = &thread_globals;
    gkStream    *fs = new gkStream(thread_globals.gkp, 1, thread_globals.gkp->gkStore_getNumFragments(), GKFRAGMENT_INF);
    uint32 numFrags = loadFragments(fs, thread_globals.libToInclude, thread_globals.frgToLen, thread_globals.frgToLib, thread_globals.frgToMate);
    AS_IID firstFrag = partitionWork(numFrags, thread_globals.libToInclude, thread_globals.frgToLib, thread_globals.numThreads, thread_globals.partitions, thread_globals.perFile, thread_wa);

    // record info on which fragments surround gaps and can help recruit long reads for correction
    thread_globals.bitMax = 0;
    thread_globals.gappedReadSet = initGappedReadSet(&thread_globals, firstFrag);
    delete fs;

    // create parallel threads
    if  (thread_globals.numThreads > 1) {
        pthread_attr_init (& attr);
        pthread_attr_setstacksize (& attr, THREAD_STACKSIZE);
        pthread_mutex_init (& thread_globals.overlapMutex, NULL);
        pthread_mutex_init (& thread_globals.globalDataMutex, NULL);
        pthread_mutex_init (& thread_globals.gkpMutex, NULL);
        pthread_mutex_init (& thread_globals.countMutex, NULL);
    }

    for (int i = 1; i < thread_globals.numThreads; i++) {
        thread_wa[i].globals = &thread_globals;
        int status = pthread_create(&thread_id[i], &attr, correctFragments, &(thread_wa[i]));
        if (status != 0) {
            fprintf (stderr, "pthread_create error at line %d:  %s\n",  __LINE__, strerror (status));
            exit (-3);
        }
    }
    // create the 0 thread (this thread is always created since we have at least one thread always)
    correctFragments(&(thread_wa[0]));
    for (int i = 1; i < thread_globals.numThreads; i++) {
        void *ptr;
        int status = pthread_join(thread_id[i], &ptr);
        if (status != 0) {
            fprintf (stderr, "pthread_create error at line %d:  %s\n",  __LINE__, strerror (status));
            exit (-3);
        }
    }

    char command[FILENAME_MAX] = {0};
    if (thread_globals.hasMates == true) {
        // keep well-mated reads, first estimate the library insert sizes
        fprintf(stderr, "Re-estimating insert sizes\n");
        thread_globals.mpa = new matePairAnalysis((char *)thread_globals.gkp->gkStore_path());
        delete[] thread_wa;
        thread_wa = new PBRThreadWorkArea[thread_globals.numThreads];
        delete[] thread_id;
        thread_id = new pthread_t[thread_globals.numThreads];
        thread_wa[0].globals = &thread_globals;
        thread_wa[0].id = 0;

        for  (uint32 i = 0;  i < thread_globals.partitions; i++ ) {
            if (drand48() < thread_globals.percentToEstimateInserts) {
                thread_globals.toOutput.push(pair<AS_IID, AS_IID>(i+1, i+1));
                fprintf(stderr, "Using partition %d to re-estimate insert sizes, total partitions so far %d/%d\n", i+1, thread_globals.toOutput.size(), thread_globals.partitions);
            }
        }
        for (int i = 1; i < thread_globals.numThreads; i++) {
            thread_wa[i].globals = &thread_globals;
            thread_wa[i].id = i;
            int status = pthread_create(&thread_id[i], &attr, estimateInsertSizes, &(thread_wa[i]));
            if (status != 0) {
                fprintf (stderr, "pthread_create error at line %d:  %s\n",  __LINE__, strerror (status));
                exit (-3);
            }
        }
        // create the 0 thread
        estimateInsertSizes(&(thread_wa[0]));
        for  (int i = 1;  i < thread_globals.numThreads;  i ++) {
            void  * ptr;
            int status = pthread_join  (thread_id [i], & ptr);
        }

        // record the updated sizes
        thread_globals.mpa->finalize();
        thread_globals.mpa->printSummary(stderr);
        for (uint32 i=1; i<=thread_globals.gkp->gkStore_getNumLibraries(); i++) {
            if (thread_globals.mpa->mean(i) != 0) {
                thread_globals.libToSize[i] = pair<double, double>(thread_globals.mpa->mean(i) - CGW_CUTOFF * thread_globals.mpa->stddev(i), thread_globals.mpa->mean(i) + CGW_CUTOFF * thread_globals.mpa->stddev(i));
            }
        }
        delete thread_globals.mpa;

        fprintf(stderr, "Filtering mates\n");
        // filter using the mates
        delete[] thread_wa;
        thread_wa = new PBRThreadWorkArea[thread_globals.numThreads];
        delete[] thread_id;
        thread_id = new pthread_t[thread_globals.numThreads];
        thread_wa[0].globals = &thread_globals;
        thread_wa[0].id = 0;

        for  (uint32 i = 0;  i < thread_globals.partitions; i++ ) {
            thread_globals.toOutput.push(pair<AS_IID, AS_IID>(i+1, i+1));
        }
        for (int i = 1; i < thread_globals.numThreads; i++) {
            thread_wa[i].globals = &thread_globals;
            thread_wa[i].id = i;
            int status = pthread_create(&thread_id[i], &attr, screenBadMates, &(thread_wa[i]));
            if (status != 0) {
                fprintf (stderr, "pthread_create error at line %d:  %s\n",  __LINE__, strerror (status));
                exit (-3);
            }
        }
        // create the 0 thread
        screenBadMates(&(thread_wa[0]));
        for  (int i = 1;  i < thread_globals.numThreads;  i ++) {
            void  * ptr;
            int status = pthread_join  (thread_id [i], & ptr);
        }

        // create and open the filtered store
        sprintf(command, "rm -rf %s.paired.ovlStore", thread_globals.prefix);
        assert(system(command) == 0);
        sprintf(command, "find . \\( -name %s\\*ovb \\) -print > %s.paired.list", thread_globals.prefix, thread_globals.prefix);
        assert(system(command) == 0);
        sprintf(command, "%soverlapStoreBuild -o %s.paired.ovlStore -g %s -M 8192 -L %s.paired.list", thread_globals.exePath, thread_globals.prefix, gkpStorePath, thread_globals.prefix);
        assert(system(command) == 0);
        sprintf(thread_globals.ovlStoreUniqPath, "%s.paired.ovlStore", thread_globals.prefix);
        fprintf(stderr, "Updated overlap store path to be %s\n", thread_globals.ovlStoreUniqPath);

        // remove intermediate files
        sprintf(command, "rm -rf %s\\*ovb", thread_globals.prefix);
        assert(system(command) == 0);
    }

    // filter repeat reads out, currently based on coverage pre-mate filtering, should it be post-mate filtering?
    fprintf(stderr, "Filtering repeats\n");
    computeRepeatThreshold(thread_globals, firstFrag);
    if (thread_globals.fixedMemory == FALSE) {
        loadSequence(thread_globals.gkp, thread_globals.readsToPrint, thread_globals.frgToEnc);
    }
    thread_globals.readsToPrint.clear();

    // remove our paired store
    sprintf(command, "rm -rf %s.paired.ovlStore", thread_globals.prefix);
    assert(system(command) == 0);

    // output our tiling
    delete[] thread_wa;
    thread_wa = new PBRThreadWorkArea[thread_globals.numThreads];
    delete[] thread_id;
    thread_id = new pthread_t[thread_globals.numThreads];
    thread_wa[0].globals = &thread_globals;
    thread_wa[0].id = 0;
    thread_globals.readsToPrint.clear();

    for  (uint32 i = 0;  i < thread_globals.partitions; i++ ) {
        thread_globals.toOutput.push(pair<AS_IID, AS_IID>(i+1, i+1));
    }
    for  (i = 1;  i < thread_globals.numThreads ;  i ++) {
        thread_wa[i].globals = &thread_globals;
        thread_wa[i].id = i;
        int status = pthread_create(&thread_id[i], &attr, outputResults, &thread_wa[i]);
        if  (status != 0) {
            fprintf (stderr, "pthread_create error at line %d:  %s\n",  __LINE__, strerror (status));
            exit (-3);
        }
    }
    // create the 0 thread
    outputResults(&(thread_wa[0]));
    for  (i = 1;  i < thread_globals.numThreads;  i ++) {
        void  * ptr;
        int status = pthread_join  (thread_id [i], & ptr);
    }

    // clean up
    fprintf(stderr, "Cleaning up\n");
    if (thread_globals.numThreads > 1) {
        pthread_mutex_destroy(& thread_globals.overlapMutex);
        pthread_mutex_destroy (& thread_globals.globalDataMutex);
        pthread_mutex_destroy (& thread_globals.gkpMutex);
        pthread_mutex_destroy (& thread_globals.countMutex);
    }
    for (map<AS_IID, char*>::iterator iter = thread_globals.frgToEnc.begin(); iter != thread_globals.frgToEnc.end(); iter++) {
        delete[] iter->second;
    }
    thread_globals.frgToEnc.clear();
    delete thread_globals.gappedReadSet;
    delete[] thread_globals.partitionStarts;
    delete[] thread_globals.libToInclude;
    delete[] thread_id;
    delete[] thread_wa;
    delete thread_globals.gkp;
}
