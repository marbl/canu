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

using namespace std;

#include "AS_global.H"
#include "AS_PBR_store.hh"

#include <vector>

int
main(int argc, char **argv) {

    //  Options for everyone.  Everyone needs a filename!
    //
    char            *fileName    = NULL;

    char *progName = argv[0];
    LayRecordStore* store = NULL;

    argc = AS_configure(argc, argv);

    int arg = 1;
    int err = 0;
    int hlp = 0;
    while (arg < argc) {
        if        (strcmp(argv[arg++], "-l") == 0) {
            fileName = argv[arg++];
        } else {
            fprintf(stderr, "Invalid option: %s\n", argv[arg++]);
            err++;
        }
    }

    if (err || fileName == NULL) {
        fprintf(stderr, "Error: invalid options specified. Usage: pacBioDumpLayout -l <layout file>\n");
        exit(1);
    }

    store = openLayFile(fileName);
    if (store == NULL) {
        exit(1);
    }

    // read the layouts and dump them in a human-readable format
    LayRecord layout;
    while (readLayRecord(store, layout)) {
        fprintf(stdout, "LAY\t"F_IID"\t"F_U32"\n", layout.iid, layout.mp.size());
        for (vector<OverlapPos>::const_iterator iter = layout.mp.begin(); iter != layout.mp.end(); iter++) {
            fprintf(stdout, "TLE\t"F_IID"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_S32"\t"F_S32"\n", iter->ident, iter->position.bgn, iter->position.end, layout.bClrs[iter->ident].bgn, layout.bClrs[iter->ident].end, layout.bOvls[iter->ident].dat.ovl.a_hang, layout.bOvls[iter->ident].dat.ovl.b_hang);
        }
    }
    closeLayFile(store);
}
