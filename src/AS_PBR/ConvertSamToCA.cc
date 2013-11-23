#include "AS_UTL_fileIO.H"
#include "AS_UTL_IID.H"
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>

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

    map <string, AS_IID > eidToIID;
    map <AS_IID, uint32> iidToLen;
    uint32 counter = 0;


    stringstream inputIIDS(argv[2]);
    char inputFile[FILENAME_MAX];

    while (inputIIDS.getline(inputFile, FILENAME_MAX, ',')) {
       cerr << "Processing input file " << inputFile << " from " << argv[2] << endl;

       ifstream iids(inputFile, ifstream::in);
       while (iids.good()) {
          string eid;
          string iid;

          iids >> eid >> iid;
          eidToIID[eid] = atoi(iid.c_str());
          if (counter % 1000000 == 0) cerr << "Loaded " << counter << " eid entries" << endl;
          counter++;
       }
       iids.close();
    }

    if (argc >= 4) {
        stringstream inputLens(argv[3]);
        while (inputLens.getline(inputFile, FILENAME_MAX, ',')) {
           cerr << "Processing input file " << inputFile << " from " << argv[3] << endl;
           ifstream iids(inputFile, ifstream::in);

           counter = 0;
           while (iids.good()) {
              string iidStr;
              string lenStr;

              iids >> iidStr >> lenStr;
              uint32 len = atol(lenStr.c_str());

              if (len == 0) {
                continue;
              }
              if (eidToIID.find(iidStr) != eidToIID.end()) {
                iidToLen[eidToIID[iidStr]] = len;
              } else {
                AS_IID iid = atoi(iidStr.c_str());
                iidToLen[iid] = len;
              }
              if (counter % 1000000 == 0) cerr << "Loaded " << counter << " eid entries" << endl;
              counter++;
           }
           iids.close();
       }
    }
    cerr << "Loaded " << eidToIID.size() << " and " << iidToLen.size() << " eid/length translation" << endl;

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
        AS_IID biid;
        if (eidToIID.find(id) == eidToIID.end() && id.find_last_of(",")
                != string::npos) {
            id = id.substr(0, id.find_last_of(","));

            if (eidToIID.find(id) == eidToIID.end() && id.find_last_of("/")
                    != string::npos) {
                id = id.substr(0, id.find_last_of("/"));
                if (eidToIID.find(id) == eidToIID.end()) {
                    cerr << "Could not find iid for " << id << endl;
                    continue;
                }
                biid = eidToIID[id];
            } else {
                biid = eidToIID[id];
            }
        } else {
            if (eidToIID.find(id) == eidToIID.end()) {
                cerr << "Could not find iid for " << id << endl;
                continue;
            }
            biid = eidToIID[id];
        }

        string ref = head->target_name[core->tid];
        AS_IID refiid;
        if (eidToIID.find(ref) == eidToIID.end() && ref.find_last_of(",")
                != string::npos) {
            ref = ref.substr(0, ref.find_last_of(","));
            if (eidToIID.find(ref) == eidToIID.end()) {
                cerr << "Could not find iid for " << ref << endl;
                continue;
            }
            refiid = eidToIID[ref];
        } else {
            if (eidToIID.find(ref) == eidToIID.end()) {
                cerr << "Could not find iid for " << ref << endl;
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
