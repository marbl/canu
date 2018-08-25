
//  Responsible for taking one kmer at a time and writing
//  it to 

class kmerCountFileWriter;

class kmerCountStreamWriter {
public:
  kmerCountStreamWriter(kmerCountFileWriter *writer,
                        uint32               fileNumber);

  ~kmerCountStreamWriter();

public:
  void    addMer(kmer k, uint32 c);

private:
  void    dumpBlock(uint64 nextPrefix=UINT64_MAX);

private:
  kmerCountFileWriter       *_writer;
  char                       _outName[FILENAME_MAX+1];

  //  Encoding data

  uint32                     _prefixSize;

  uint32                     _suffixSize;
  uint64                     _suffixMask;

  uint32                     _numFilesBits;
  uint32                     _numBlocksBits;
  uint64                     _numFiles;
  uint64                     _numBlocks;

  //  File data

  uint32                     _filePrefix;

  FILE                      *_datFile;
  kmerCountFileIndex        *_datFileIndex;

  //  Kmer data and etc for writing the stream.

  uint64                     _batchPrefix;
  uint64                     _batchNumKmers;
  uint64                     _batchMaxKmers;
  uint64                    *_batchSuffixes;
  uint32                    *_batchCounts;

  //kmerCountStatistics        _stats;
};



