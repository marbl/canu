
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



//  This object is allocated in
//    sqStore::sqStore                 (which calls sqLibrary_initialize to clear things)
//    sqStore::sqStore_loadMetadata    (which loads from disk and does not call initialize)
//    sqStore::sqStore_addEmptyLibrary (which calls sqLibrary_initialize)
//
class sqLibrary {
public:
  void             sqLibrary_initialize(char const    *name = nullptr,
                                        uint32         id   = 0,
                                        sqLibrary_tech tech = sqTechType_generic) {

    //  Bullet proof the library name - so we can make files with this prefix.

    memset(_libraryName, 0, sizeof(char) * LIBRARY_NAME_SIZE);

    if (name) {
      uint32  lp = 0;
      
      for (char const *orig=name; *orig; orig++) {
        if        (*orig == '/') {
          _libraryName[lp++] = '_';

        } else if (isspace(*orig) == 0) {
          _libraryName[lp++] = *orig;

        } else {
          _libraryName[lp++] = '_';
        }

        if (lp >= LIBRARY_NAME_SIZE) {
          _libraryName[LIBRARY_NAME_SIZE-1] = 0;
          break;
        }
      }
    }

    else {
      strncpy(_libraryName, "UNDEFINED", LIBRARY_NAME_SIZE-1);
    }

    //  Set the boring stuff.

    _libraryID = id;
    _techType  = tech;

    //  Clear the unused stuff.

    _U_keyLen    = 0;
    _U_keyAlloc  = 0;
    _U_keys      = NULL;
    _U_values    = NULL;

    _U_kvLen     = 0;
    _U_kvAlloc   = 0;
    _U_kvData    = NULL;
  };


public:
  char const      *sqLibrary_libraryName(void)           { return(_libraryName);          };
  uint32           sqLibrary_libraryID(void)             { return(_libraryID);            };

  sqLibrary_tech   sqLibrary_techType(void)              { return(_techType);             };
  char const      *sqLibrary_techTypeString(void)        { return(toString(_techType));   };

  //char          *sqLibrary_lookupKey(char const *key)  { return(NULL);                  };

private:
  char            _libraryName[LIBRARY_NAME_SIZE];
  uint32          _libraryID;

  sqLibrary_tech  _techType;

  //  The rest of these are unused.  They used to hold a set of key/value pairs
  //  for library features.  Canu never really used them, and they were removed.
  //  But to keep file compatibility, we still keep them around.
  //
  uint32          _U_keyLen;     //  Number of key/value pairs we have.
  uint32          _U_keyAlloc;   //  Number of key/value pairs we have space for.
  uint32         *_U_keys;       //  Offset into _kvData for each NUL-terminated key string.
  uint32         *_U_values;     //  Offset into _kvData for each NUL-terminated value string.

  uint32          _U_kvLen;      //  Butes of kvData we have used.
  uint32          _U_kvAlloc;    //  Bytes of kvData we have allocated.
  char           *_U_kvData;     //  Bulk list of key and value strings.

  friend class sqStore;
};


#endif  //  SQLIBRARY_H
