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
  u32bit   magic_number;
  char     version[4];
  u64bit   index_offset;
  u32bit   index_length;
  u32bit   number_of_reads;
  u16bit   header_length;
  u16bit   key_length;
  u16bit   number_of_flows_per_read;
  u8bit    flowgram_format_code;

  char     flow_chars[SFF_NUMBER_OF_FLOWS_MAX];     //  h->number_of_flows_per_read
  char     key_sequence[SFF_KEY_SEQUENCE_MAX];      //  h->key_length

  u32bit   swap_endianess;
};


struct sffRead {
  //  The next block is read in one swoop from the sff file.  DO NOT MODIFY!
  u16bit   read_header_length;
  u16bit   name_length;
  u32bit   number_of_bases;
  u16bit   clip_quality_left;
  u16bit   clip_quality_right;
  u16bit   clip_adapter_left;
  u16bit   clip_adapter_right;

  char     name[SFF_NAME_LENGTH_MAX];                     //  r->name_length

  u16bit   flowgram_values[SFF_NUMBER_OF_FLOWS_MAX];      //  h->number_of_flows_per_read
  u8bit    flow_index_per_base[SFF_NUMBER_OF_BASES_MAX];  //  r->number_of_bases
  char     bases[SFF_NUMBER_OF_BASES_MAX];                //  r->number_of_bases
  u8bit    quality_scores[SFF_NUMBER_OF_BASES_MAX];       //  r->number_of_bases

  char     quality[SFF_NUMBER_OF_BASES_MAX];              //  quality_scores converted to CA-format qv
};

struct sffIndex {
  u64bit  _seqPos;
  u32bit  _seqLen;
  u32bit  _namLen;
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
  u32bit        find(const char *sequencename) {
    assert(0);
    return(0);
  };

  u32bit        getSequenceLength(u32bit iid) { return(_index[iid]._seqLen); };

  bool          getSequence(u32bit iid,
                            char *&h, u32bit &hLen, u32bit &hMax,
                            char *&s, u32bit &sLen, u32bit &sMax);
  bool          getSequence(u32bit iid,
                            u32bit bgn, u32bit end, char *s);

private:
  void          clear(void);

  readBuffer  *_rb;

  sffHeader    _header;
  sffRead      _read;

  sffIndex    *_index;

  u64bit       _firstReadLocation;
  u64bit       _readIID;

  friend class seqFactory;
};


#endif  //  SFF_H
