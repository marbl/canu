#include "AS_UTL_fileIO.h"
#include <map>
#include <iostream>
#include <fstream>

#include "sam.h"

using namespace std;

/**
 * Check if read is properly mapped
 * @return true if read mapped, false otherwise
 */
static bool is_mapped(const bam1_core_t *core) {

    if (core->flag & BAM_FUNMAP) {
        return false;
    }

    return true;
}

samfile_t * open_alignment_file(std::string path) {
    samfile_t * fp = NULL;
    std::string flag = "r";
    if (path.substr(path.size() - 3).compare("bam") == 0) {
        //BAM file!
        flag += "b";
    }
    if ((fp = samopen(path.c_str(), flag.c_str(), 0)) == 0) {
        fprintf(stderr, "qaCompute: Failed to open file %s\n", path.c_str());
    }
    return fp;
}

int main(int argc, char * argv[]) {

    bool debug = false;
    if (argc < 2) {
        cerr << "Usage: " << argv[0]
                << " <sam file> <eid to iid mapping> [iid to length mapping]"
                << endl;
        exit(1);
    }
    if (argc > 4) {
        debug = true;
    }

    ifstream iids(argv[2], ifstream::in);
    map < string, string > eidToIID;
    map<string, uint32> iidToLen;

    while (iids.good()) {
        string illumina;
        string corrected;

        iids >> illumina >> corrected;
        eidToIID[illumina] = corrected;
    }

    iids.close();

    if (argc >= 4) {
        iids.open(argv[3], ifstream::in);

        while (iids.good()) {
            string illumina;
            string corrected;

            iids >> illumina >> corrected;
            uint32 correctedInt = atol(corrected.c_str());

            if (correctedInt == 0) {
                continue;
            }
            if (eidToIID.find(illumina) != eidToIID.end()) {
                iidToLen[eidToIID[illumina]] = correctedInt;
            } else {
                iidToLen[illumina] = correctedInt;
            }
        }
        iids.close();
    }

    bam1_t *b = bam_init1();
    samfile_t * fp = open_alignment_file(argv[1]);
    if (fp == NULL) {
        exit(1);
    }
    bam_header_t* head = fp->header; // sam header
    while (samread(fp, b) >= 0) {
        //Get bam core.
        const bam1_core_t *core = &b->core;
        if (core == NULL) {
            printf("Input file is corrupt!");
            exit(1);
        }
        if (!is_mapped(core)) {
            continue;
        }

        string id = bam1_qname(b);
        string biid;
        if (eidToIID.find(id) == eidToIID.end() && id.find_last_of(",")
                != string::npos) {
            id = id.substr(0, id.find_last_of(","));

            if (eidToIID.find(id) == eidToIID.end() && id.find_last_of("/")
                    != string::npos) {
                id = id.substr(0, id.find_last_of("/"));
                if (eidToIID.find(id) == eidToIID.end()) {
                    continue;
                }
                biid = eidToIID[id];
            } else {
                biid = eidToIID[id];
            }
        } else {
            if (eidToIID.find(id) == eidToIID.end()) {
                continue;
            }
            biid = eidToIID[id];
        }

        string ref = head->target_name[core->tid];
        string refiid;
        if (eidToIID.find(ref) == eidToIID.end() && ref.find_last_of(",")
                != string::npos) {
            ref = ref.substr(0, ref.find_last_of(","));
            if (eidToIID.find(ref) == eidToIID.end()) {
                continue;
            }
            refiid = eidToIID[ref];
        } else {
            if (eidToIID.find(ref) == eidToIID.end()) {
                continue;
            }
            refiid = eidToIID[ref];
        }

        uint32_t* cigar = bam1_cigar(b);
        uint32_t refLen = head->target_len[core->tid];
        uint32_t refLo = core->pos + 1;
        uint32_t refHigh = bam_calend(core, cigar);
        uint32 blen = (iidToLen.find(biid) != iidToLen.end() ? iidToLen[biid]
                : core->l_qseq);
        bool isFwd = !(core->flag & BAM_FREVERSE);
        int errors = 0;
        for (int k = 0; k < core->n_cigar; ++k) {
            int cop = cigar[k] & BAM_CIGAR_MASK; // operation
            int cl = cigar[k] >> BAM_CIGAR_SHIFT; // length
            switch (cop) {
            case BAM_CMATCH:
                break;
            default:
                errors += cl;
                break;
            }
        }

        double errRate = ((double) errors / (refHigh - refLo + 1)) * 100;
        if (refHigh - refLo + 1 >= (double) 0.9 * blen) {
            int ahang = refLo;
            int bhang = -1 * (refLen - refHigh);
            //"  OVL:   aIID bIID [I|N] aHang bHang error error_corrected\n");
            cout << refiid << "\t" << biid << "\t" << (isFwd == true ? "N"
                    : "I") << "\t" << ahang << "\t" << bhang << "\t" << errRate
                    << "\t" << errRate << endl;
        }

    }
    samclose(fp);
}
