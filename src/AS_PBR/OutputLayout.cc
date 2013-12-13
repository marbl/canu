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

const char *mainid = "$Id$";

#include "AS_global.H"
#include "AS_MSG_pmesg.H"
#include "AS_PER_gkpStore.H"

#include "AS_PBR_util.hh"
#include "AS_PBR_output.hh"

#include <map>

using namespace std;

static uint32 loadFragments(gkStream *fs, uint32* includeLib, map<AS_IID, uint32>& frgToLen) {
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
        frgToLen[fr.gkFragment_getReadIID()] = len;
    }
    return counter;
}

int
main (int argc, char * argv []) {
    fprintf(stderr, "Parsing arguments\n");

    PBRThreadGlobals  thread_globals;

    char      *gkpStorePath           = NULL;
    uint32     part		      = 0;
    bool       outputAsLay            = true;

    // initialize default parameters
    thread_globals.globalRepeats     = TRUE;
    thread_globals.minLength         = 500;
    thread_globals.maxUncorrectedGap = 0;
    thread_globals.hasMates                     = false;
    thread_globals.verboseLevel 			  = VERBOSE_OFF;
    thread_globals.allowLong                = FALSE;
    thread_globals.erate             = 0.15;

    strcpy(thread_globals.prefix, "asm");

    argc = AS_configure(argc, argv);

    int err = 0;
    int arg = 1;
    while (arg < argc) {
        if (strcmp(argv[arg], "-G") == 0) {
            gkpStorePath = argv[++arg];

        } else if (strcmp(argv[arg], "-o") == 0) {
            strncpy(thread_globals.prefix, argv[++arg], FILENAME_MAX);

        } else if (strcmp(argv[arg], "-l") == 0) {
            thread_globals.minLength = atoi(argv[++arg]);

        } else if (strcmp(argv[arg], "-L") == 0) {
            thread_globals.allowLong = TRUE;

        } else if (strcmp(argv[arg], "-R") == 0) {
            thread_globals.globalRepeats = FALSE;
            thread_globals.repeatMultiplier = atof(argv[++arg]);

        } else if (strcmp(argv[arg], "-M") == 0) {
            int32 maxGap = atoi(argv[++arg]);
            if (maxGap <= 0) { thread_globals.maxUncorrectedGap = 0; }
            else { thread_globals.maxUncorrectedGap = maxGap; }

       } else if (strcmp(argv[arg], "-p") == 0) {
            int32 partition = atoi(argv[++arg]);
            if (partition <= 0) { err++; }
            else { part = partition; }

        } else if (strcmp(argv[arg], "-v") == 0) {
            thread_globals.verboseLevel = atoi(argv[++arg]);
            if (thread_globals.verboseLevel < VERBOSE_OFF) { thread_globals.verboseLevel = VERBOSE_OFF; }
            if (thread_globals.verboseLevel > VERBOSE_DEVELOPER) { thread_globals.verboseLevel = VERBOSE_DEVELOPER; }

        } else if (strcmp(argv[arg], "-e") == 0) {
            thread_globals.erate = atof(argv[++arg]);

        } else if (strcmp(argv[arg], "-P") == 0) {
           outputAsLay = false;

        } else {
            err++;
        }

        arg++;
    }

    if (gkpStorePath == NULL)
        err++;
    if (thread_globals.minLength <= 0)
        err++;

    if (err) {
        fprintf(stderr, "usage: %s -G gkpStore [options]\n", argv[0]);
        fprintf(stderr, "\n");
        fprintf(stderr, "  -G         Mandatory path to a gkpStore.\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "  -l %d     ignore corrected fragments less than %d bp\n", thread_globals.minLength, thread_globals.minLength);
        fprintf(stderr, "  -o %s     output prefix of %s\n", thread_globals.prefix, thread_globals.prefix);
        fprintf(stderr, "\n");
        fprintf(stderr, " -M %d	 	 The maximum uncorrected PacBio gap that will be allowed. When there is no short-read coverage for a region, by default the pipeline will split a PacBio sequence. This option allows a number of PacBio sequences without short-read coverage to remain. For example, specifying 50, will mean 50bp can have no short-read coverage without splitting the PacBio sequence. Warning: this will allow more sequences that went through the SMRTportal to not be fixed.\n", thread_globals.maxUncorrectedGap);
        fprintf(stderr, "\n");
        fprintf(stderr, "  -v %d     level of verbosity. Higher values generate more debugging output\n", thread_globals.verboseLevel);
        fprintf(stderr, "\n");
        fprintf(stderr, "  -R %f     consider a pileup of %f times the mean for a single corrected read to be a repeat and distribute reads to their best locations (this is only useful for metagenomic or non-even coverage datasets. Otherwise the global repeat estimate is used instead)\n", thread_globals.repeatMultiplier, thread_globals.repeatMultiplier);
        fprintf(stderr, "\n");

        if (gkpStorePath == NULL)
            fprintf(stderr, "No gatekeeper store (-G option) supplied.\n");

        if (thread_globals.minLength <= 0)
            fprintf(stderr, "Invalid min length, must be positive: %d\n", thread_globals.minLength);

        exit(1);
    }

    // set up global data structures
    fprintf(stderr, "Opening stores\n");
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

    gkStream    *fs = new gkStream(thread_globals.gkp, 1, thread_globals.gkp->gkStore_getNumFragments(), GKFRAGMENT_INF);
    uint32 numFrags = loadFragments(fs, thread_globals.libToInclude, thread_globals.frgToLen);
    outputResults(&thread_globals, part, outputAsLay);
    delete[] thread_globals.libToInclude;
    delete thread_globals.gkp;
}
