
#include <vector>

using namespace std;

class abVarRead {
public:
  abVarRead() {
    id        = 0;
    iid       = 0;
    bases     = NULL;
    qvs       = NULL;
    ave_qv    = 0;
    allele_id = 0;
    uglen     = 0;
  };
  ~abVarRead() {
    delete [] bases;
    delete [] qvs;
  };

  void   initialize(uint32 len, uint32 qv) {
    id        = UINT32_MAX;
    iid       = UINT32_MAX;
    bases     = new char   [len];
    qvs       = new uint32 [len];
    ave_qv    = 0;
    allele_id = 0;
    uglen     = 0;

    for (uint32 j=0; j<len; j++) {
      bases[j] = '-';
      qvs[j]   = qv;
    }
  };

  uint32         id;         // ?
  uint32         iid;        // gkp IID
  char          *bases;      // gapped sequence
  uint32        *qvs;        // quality values
  double         ave_qv;
  uint32         allele_id;
  uint32         uglen;      // ungapped length
};


class abVarAllele {
public:
  abVarAllele() {
    id        = 0;
    num_reads = 0;
    read_ids  = NULL;
    read_iids = NULL;
    weight    = 0;
    uglen     = 0;
  };
  ~abVarAllele() {
    delete [] read_ids;
    delete [] read_iids;
  };

  void   initialize(uint32 nr) {
    id        = UINT32_MAX;
    num_reads = 0;
    read_ids  = new uint32 [nr];  memset(read_ids,  0, sizeof(uint32) * nr);
    read_iids = new uint32 [nr];  memset(read_iids, 0, sizeof(uint32) * nr);
    weight    = 0;
    uglen     = 0;
  };

  uint32         id;
  uint32         num_reads;
  uint32        *read_ids;
  uint32        *read_iids;
  uint32         weight;
  uint32         uglen;      // ungapped length
};


//  This structure is used when recalling consensus bases
//  to use only one of two alleles
//
class abVarRegion {
public:
  abVarRegion() {
    beg         = 0;
    end         = 0;
    nb          = 0;
    na          = 0;
    nca         = 0;
    nr          = 0;

    reads       = NULL;
    alleles     = NULL;
    dist_matrix = NULL;
  }
  ~abVarRegion() {
    delete [] reads;
    delete [] alleles;
    delete [] dist_matrix;
  };

  void           allocateReads(void) {
  };

  //  Grab a chunk of memory (ok, two chunks), clear the values, and make a 2d array.
  void           allocateDistanceMatrix(int32 initDist) {
    dist_matrix      = new int32 * [nr];
    dist_matrix_vals = new int32   [nr * nr];

    for (uint32 ii=0; ii<nr * nr; ii++)
      dist_matrix_vals[ii] = initDist;
    
    for (uint32 ii=0; ii<nr; ii++)
      dist_matrix[ii] = dist_matrix_vals + ii * nr;
  };
  
  void           outputDistanceMatrix(FILE *F) {
    fprintf(F, "Distance matrix=\n");

    for (uint32 j=0; j<nr; j++) {
      for (uint32 k=0; k<nr; k++)
        fprintf(F, " %4d", dist_matrix[j][k]);

      fprintf(F, "\n");
    }
  };

  bool            phaseWithPreviousRegion(uint32* &allele_map,
                                          uint32   prev_nca,   uint32  *prev_nca_iid,   uint32  prev_nca_iid_max,
                                          uint32   prev_ncr,   uint32  *prev_ncr_iid,   uint32  prev_ncr_iid_max);

  void            setConsensusToMajorAllele(abAbacus *abacus, abColID *cids);

  void            populateDistanceMatrix(void);
  void            clusterReads(void);

  void            sortAllelesByLength(void);  //  nca
  void            sortAllelesByWeight(void);  //  na
  void            sortAllelesByMapping(uint32 *allele_map);  //  nca



  //  Debug interface
  void showConfirmedReads(void);
  void showReadsOfConfirmedAlleles(void);
  void showReads(void);


  uint32         beg;         // position of the left boundary
  uint32         end;         // position of the right boundary

  uint32         nb;          // number of "current" bases
  uint32         na;          // total number of detected alleles
  uint32         nca;         // number of confirmed alleles

  uint32         nr;          // number of reads in the region of variation
  //uint32         max_nr;

  //char          *curr_bases;  // dim = nr
  //uint32        *iids;        // iids of the reads

  vector<char>   curr_bases;  // dim = nr
  vector<uint32> iids;        // iids of the reads - gkp IID

  abVarRead     *reads;
  abVarAllele   *alleles;

  int32         *dist_matrix_vals;
  int32        **dist_matrix;       // nr x nr matrix of cross-distances between reads
};

