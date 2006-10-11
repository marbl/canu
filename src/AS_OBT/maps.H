#ifndef OBT_MAPS_H
#define OBT_MAPS_H

#include <stdio.h>
#include <stdlib.h>

using namespace std;
#include <map>

class vectorInfo {
public:
  vectorInfo() {
    vecL      = 0;
    vecR      = 0;
    hasVec    = 0;
    immutable = 0;
  };

  u32bit  vecL      : 12;
  u32bit  vecR      : 12;
  u32bit  hasVec    : 1;
  u32bit  immutable : 1;
};


class vectorMap {
public:
  vectorMap() {};
  ~vectorMap() {};

private:
  map<u64bit, vectorInfo>  m;

public:
  vectorInfo &operator[](u64bit uid) {
    return(m[uid]);
  };

  bool        exists(u64bit uid) {
    return(m.find(uid) != m.end());
  };

public:
  void readVectorMap(char *vectorFileName) {
    bool fatal = false;

    if (vectorFileName == 0L)
      return;

    errno = 0;
    FILE *intFile = fopen(vectorFileName, "r");
    if (errno)
      fprintf(stderr, "Can't open '%s': %s\n", vectorFileName, strerror(errno)), exit(1);

    char intLine[1024] = {0};
    fgets(intLine, 1024, intFile);
    while (!feof(intFile)) {
      chomp(intLine);
      splitToWords  W(intLine);
      if ((W[0] == 0L) || (W[1] == 0L) || (W[2] == 0L) || (W[3] != 0L)) {
        fprintf(stderr, "readVectorMap()-- Invalid line '%s'\n", intLine);
        fatal = true;
      } else {

        //  These look base-based!

        u64bit  uid  = strtou64bit(W[0], 0L);
        u32bit  intl = strtou32bit(W[1], 0L) - 1;
        u32bit  intr = strtou32bit(W[2], 0L);

        //  Silently swap, if needed

        if (intl > intr) {
          u32bit  s = intl;
          intl = intr;
          intr = s;
        }

        vectorInfo v;
        v.vecL      = intl;
        v.vecR      = intr;
        v.hasVec    = 1;
        v.immutable = 0;

        m[uid] = v;
      }

      fgets(intLine, 1024, intFile);
    }

    fclose(intFile);

    if (fatal)
      exit(1);
  };


  void readImmutableMap(char *immutableFileName) {
    bool fatal = false;

    if (immutableFileName == 0L)
      return;

    errno = 0;
    FILE *intFile = fopen(immutableFileName, "r");
    if (errno)
      fprintf(stderr, "Can't open '%s': %s\n", immutableFileName, strerror(errno)), exit(1);

    char intLine[1024] = {0};
    fgets(intLine, 1024, intFile);
    while (!feof(intFile)) {
      chomp(intLine);
      splitToWords  W(intLine);
      if ((W[0] == 0L) || (W[1] != 0L)) {
        fprintf(stderr, "readImmutableMap()-- Invalid line '%s'\n", intLine);
        fatal = true;
      } else {
        u64bit  uid  = strtou64bit(W[0], 0L);
        m[uid].immutable = 1;
      }

      fgets(intLine, 1024, intFile);
    }

    fclose(intFile);

    if (fatal)
      exit(1);
  }
};




#endif  //  OBT_MAPS_H

