#ifndef OBT_MAPS_H
#define OBT_MAPS_H

#include <stdio.h>
#include <stdlib.h>

#include "util++.H"

class vectorInfo {
public:
  vectorInfo() {
    vecL      = 0;
    vecR      = 0;
    hasVec    = 0;
    immutable = 0;
  };

  uint32  vecL      : 12;
  uint32  vecR      : 12;
  uint32  hasVec    : 1;
  uint32  immutable : 1;
};


class vectorMap {
public:
  vectorMap(GateKeeperStore *gkp) {
    _gkp = gkp;
    _inf = new vectorInfo [getLastElemFragStore(_gkp) + 1];
  };
  ~vectorMap() {
    delete [] _inf;
  };

private:
  GateKeeperStore          *_gkp;
  vectorInfo               *_inf;

public:
  vectorInfo &operator[](uint32 iid) {
    return(_inf[iid]);
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
        CDS_UID_t      uid  = STR_TO_UID(W[0], 0L, 10);
        PHashValue_AS  value;

        //  If the UID exists, set the values
        if (HASH_FAILURE != getGatekeeperUIDtoIID(_gkp, uid, &value)) {
          uint32 iid = value.IID;

          //  These look base-based!
          _inf[iid].vecL   = atoi(W[1]) - 1;
          _inf[iid].vecR   = atoi(W[2]);
          _inf[iid].hasVec = 1;

          //  Silently swap, if needed
          if (_inf[iid].vecL > _inf[iid].vecR) {
            uint32       s = _inf[iid].vecL;
            _inf[iid].vecL = _inf[iid].vecR;
            _inf[iid].vecR = s;
          }
        }
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
        CDS_UID_t      uid  = STR_TO_UID(W[0], 0L, 10);
        PHashValue_AS  value;

        if (HASH_FAILURE != getGatekeeperUIDtoIID(_gkp, uid, &value))
          _inf[value.IID].immutable = 1;
      }

      fgets(intLine, 1024, intFile);
    }

    fclose(intFile);

    if (fatal)
      exit(1);
  }
};




#endif  //  OBT_MAPS_H

