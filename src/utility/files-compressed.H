
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/utility/AS_UTL_fasta.H
 *    src/utility/AS_UTL_fileIO.H
 *    src/utility/files.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-JUL-20
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef FILES_COMPRESSED_H
#define FILES_COMPRESSED_H

//  Do not include directly.  Use 'files.H' instead.

enum cftType {
  cftNONE  = 0,
  cftGZ    = 1,
  cftBZ2   = 2,
  cftXZ    = 3,
  cftSTDIN = 4
};

cftType  compressedFileType(char const *filename);



class compressedFileReader {
public:
  compressedFileReader(char const *filename);
  ~compressedFileReader();

  FILE *operator*(void)     {  return(_file);              };
  FILE *file(void)          {  return(_file);              };

  char *filename(void)      {  return(_filename);          };

  bool  isCompressed(void)  {  return(_pipe == true);      };
  bool  isNormal(void)      {  return((_pipe == false) &&
                                      (_stdi == false));   };

private:
  FILE  *_file;
  char  *_filename;
  bool   _pipe;
  bool   _stdi;
};



class compressedFileWriter {
public:
  compressedFileWriter(char const *filename, int32 level=1);
  ~compressedFileWriter();

  FILE *operator*(void)     {  return(_file);          };
  FILE *file(void)          {  return(_file);          };

  char *filename(void)      {  return(_filename);          };

  bool  isCompressed(void)  {  return(_pipe == true);  };

private:
  FILE  *_file;
  char  *_filename;
  bool   _pipe;
  bool   _stdi;
};



#endif  //  FILES_COMPRESSED_H
