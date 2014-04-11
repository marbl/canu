#ifndef SFF_H
#define SFF_H

#include "util++.H"
#include "bio++.H"

#include "seqFile.H"

#define SFF_KEY_SEQUENCE_MAX         64

#define SFF_NAME_LENGTH_MAX         256
#define SFF_NUMBER_OF_FLOWS_MAX     512
#define SFF_NUMBER_OF_BASES_MAX    2048  //  The assembler itself cannot handle longer


struct sffHeader {
  //  The next block is read in one swoop from the sff file.  DO NOT MODIFY!
  uint32   magic_number;
  char     version[4];
  uint64   index_offset;
  uint32   index_length;
  uint32   number_of_reads;
  uint16   header_length;
  uint16   key_length;
  uint16   number_of_flows_per_read;
  uint8    flowgram_format_code;

  char     flow_chars[SFF_NUMBER_OF_FLOWS_MAX];     //  h->number_of_flows_per_read
  char     key_sequence[SFF_KEY_SEQUENCE_MAX];      //  h->key_length

  uint32   swap_endianess;
};


struct sffRead {
  //  The next block is read in one swoop from the sff file.  DO NOT MODIFY!
  uint16   read_header_length;
  uint16   name_length;
  uint32   number_of_bases;
  uint16   clip_quality_left;
  uint16   clip_quality_right;
  uint16   clip_adapter_left;
  uint16   clip_adapter_right;

  char     name[SFF_NAME_LENGTH_MAX];                     //  r->name_length

  uint16   flowgram_values[SFF_NUMBER_OF_FLOWS_MAX];      //  h->number_of_flows_per_read
  uint8    flow_index_per_base[SFF_NUMBER_OF_BASES_MAX];  //  r->number_of_bases
  char     bases[SFF_NUMBER_OF_BASES_MAX];                //  r->number_of_bases
  uint8    quality_scores[SFF_NUMBER_OF_BASES_MAX];       //  r->number_of_bases

  char     quality[SFF_NUMBER_OF_BASES_MAX];              //  quality_scores converted to CA-format qv
};

struct sffIndex {
  uint64  _seqPos;
  uint32  _seqLen;
  uint32  _namLen;
};


class sffFile : public seqFile {
protected:
  sffFile(const char *filename);
  sffFile();

public:
  ~sffFile();

protected:
  seqFile      *openFile(const char *name);

public:
  uint32        find(const char *sequencename) {
    assert(0);
    return(0);
  };

  uint32        getSequenceLength(uint32 iid) { return(_index[iid]._seqLen); };

  bool          getSequence(uint32 iid,
                            char *&h, uint32 &hLen, uint32 &hMax,
                            char *&s, uint32 &sLen, uint32 &sMax);
  bool          getSequence(uint32 iid,
                            uint32 bgn, uint32 end, char *s);

private:
  void          clear(void);

  readBuffer  *_rb;

  sffHeader    _header;
  sffRead      _read;

  sffIndex    *_index;

  uint64       _firstReadLocation;
  uint64       _readIID;

  friend class seqFactory;
};


#endif  //  SFF_H
