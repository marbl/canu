
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#ifndef SQLIBRARY_H
#define SQLIBRARY_H

//  DO NOT INCLUDE THIS FILE DIRECTLY, include sqStore.H.

typedef enum {
  sqTechType_generic      = 0x0000,
  sqTechType_contig       = 0x0001,
  sqTechType_pacbio       = 0x0002,
  sqTechType_nanopore     = 0x0003,
  sqTechType_pacbio_hifi  = 0x0004,
} sqLibrary_tech;

static
char const *
toString(sqLibrary_tech t) {
  switch (t) {
    case sqTechType_generic:       return("generic");     break;
    case sqTechType_contig:        return("contig");      break;
    case sqTechType_pacbio:        return("PacBio");      break;
    case sqTechType_nanopore:      return("Nanopore");    break;
    case sqTechType_pacbio_hifi:   return("PacBioHiFi");  break;
    default:                       assert(0);             break;
  }

  return(NULL);
}



class sqLibrary {
public:
  sqLibrary(char *name=NULL, sqLibrary_tech tech=sqTechType_generic) {

    memset(_libraryName, 0, sizeof(char) * LIBRARY_NAME_SIZE);

    if (name)
      strncpy(_libraryName, name,        LIBRARY_NAME_SIZE-1);
    else
      strncpy(_libraryName, "UNDEFINED", LIBRARY_NAME_SIZE-1);

    _libraryID = UINT32_MAX;

    _techType  = tech;

    _keyLen    = 0;
    _keyAlloc  = 0;
    _keys      = NULL;
    _values    = NULL;

    _kvLen     = 0;
    _kvAlloc   = 0;
    _kvData    = NULL;
  }

  ~sqLibrary() {
    delete [] _keys;
    delete [] _values;
    delete [] _kvData;
  };


public:
  char const      *sqLibrary_libraryName(void)           { return(_libraryName);          };
  uint32           sqLibrary_libraryID(void)             { return(_libraryID);            };

  sqLibrary_tech   sqLibrary_techType(void)              { return(_techType);             };
  char const      *sqLibrary_techTypeString(void)        { return(toString(_techType));   };

  char            *sqLibrary_lookupKey(char const *key)  { return(NULL);                  };

private:
  char            _libraryName[LIBRARY_NAME_SIZE];
  uint32          _libraryID;

  sqLibrary_tech  _techType;

  uint32          _keyLen;     //  Number of key/value pairs we have.
  uint32          _keyAlloc;   //  Number of key/value pairs we have space for.
  uint32         *_keys;       //  Offset into _kvData for each NUL-terminated key string.
  uint32         *_values;     //  Offset into _kvData for each NUL-terminated value string.

  uint32          _kvLen;      //  Butes of kvData we have used.
  uint32          _kvAlloc;    //  Bytes of kvData we have allocated.
  char           *_kvData;     //  Bulk list of key and value strings.

  friend class sqStore;
};


#endif  //  SQLIBRARY_H
