#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

extern "C" {

#include "AS_global.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_fragStore.h"
#include "AS_GKP_include.h"
}


//  Removes from the gatekeeperstore any links to deleted fragments.
//
//  XXX:  Does not delete the fragments from gatekeeper.  
GateKeeperStore gkpStore;



void
deleteLinksForFrag(CDS_IID_t fragIID) {
  GateKeeperFragmentRecord gkFrag1;
  GateKeeperFragmentRecord gkFrag2;

  getGateKeeperFragmentStore(gkpStore.frgStore, fragIID, &gkFrag1);

  GateKeeperLinkRecordIterator gkpLinks;
  GateKeeperLinkRecord         gkpLink;

  //  XXX  Once we fix the other CreateGateKeeperLinkRecordIterator()
  //  to listen to the return value, we can get rid of this test.

  if (gkFrag1.linkHead == 0)
    return;

  if (CreateGateKeeperLinkRecordIterator(gkpStore.lnkStore,
                                         gkFrag1.linkHead,
                                         fragIID,
                                         &gkpLinks))
    return;

  fprintf(stdout, "delete all links to frag iid " F_IID "\n", fragIID);
  
  while(NextGateKeeperLinkRecordIterator(&gkpLinks, &gkpLink)) {
    if (gkpLink.type == AS_MATE) {
      GateKeeperLinkRecord     newLink;

      newLink.frag1       = gkpLink.frag1;
      newLink.frag2       = gkpLink.frag2;
      newLink.deleted     = 0;
      newLink.frag1Next   = 0;
      newLink.frag2Next   = 0;
      newLink.type        = AS_MATE;
      newLink.distance    = 0;                      // Matches any
      newLink.orientation = AS_GKP_UNKNOWN;         // Matches any

      int deleteLink = findLink(gkpStore.lnkStore, fragIID, gkFrag1.linkHead, &newLink, &newLink);

      if (deleteLink == 0) {
        fprintf(stdout, "can't find link of type %c for iids "F_UID" and "F_UID"\n",
                newLink.type, gkpLink.frag1, gkpLink.frag2);
      } else {
        // idew 2006-12-01
        // fragIID may be gkpLink.frag1 or gkpLink.frag2
        // gkFrag1 may be data for gkpLink.frag1 or gkpLink.frag2
        // need to make it so gkFrag2 holds data for non-fragIID frag
        getGateKeeperFragmentStore(gkpStore.frgStore,
             (gkpLink.frag2 == fragIID ? gkpLink.frag1 : gkpLink.frag2),
                                   &gkFrag2);

        fprintf(stdout, "delete link %d of type %c for iids "F_UID" and "F_UID"\n",
                deleteLink, newLink.type, gkpLink.frag1, gkpLink.frag2);

        // idew 2006-12-01
        // see above comment
        unlinkLink_GKP(gkpStore.lnkStore,
                       gkpStore.frgStore,
                       gkpLink.frag1,
                       gkpLink.frag2,
                       (gkpLink.frag1 == fragIID ? &gkFrag1 : &gkFrag2),
                       (gkpLink.frag1 == fragIID ? &gkFrag2 : &gkFrag1),
                       &newLink,
                       deleteLink);

        UnRefPHashTable_AS(gkpStore.hashTable, UID_NAMESPACE_AS, newLink.distance);
        UnRefPHashTable_AS(gkpStore.hashTable, UID_NAMESPACE_AS, gkFrag2.readUID);
        UnRefPHashTable_AS(gkpStore.hashTable, UID_NAMESPACE_AS, gkFrag1.readUID);

#if 0
        //  This doesn't work.  It just asserts.
        //
        if(HASH_SUCCESS == DeleteFromPHashTable_AS(gkpStore.hashTable, UID_NAMESPACE_AS, gkFrag1.readUID))
          deleteAndMarkGateKeeperFragmentStore(gkpStore.frgStore, fragIID, gkFrag1.birthBatch);
        else
          assert(0);
#endif
      }
    }
  }
}


int
main(int argc, char **argv) {
  char *frgStoreName = 0L;
  char *gkpStoreName = 0L;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-f") == 0) {
      frgStoreName = argv[++arg];
    } else if (strcmp(argv[arg], "-g") == 0) {
      gkpStoreName = argv[++arg];
    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      exit(1);
    }

    arg++;
  }

  ////////////////////////////////////////

  strcpy(gkpStore.storePath, gkpStoreName);
  OpenGateKeeperStore(&gkpStore);

  FragStoreHandle frgStore = openFragStore(frgStoreName, "r");
  if (frgStore == NULLSTOREHANDLE)
    fprintf(stderr, "Failed to open fragStore %s!\n", frgStoreName), exit(1);

  uint32            fe      = getFirstElemFragStore(frgStore);
  uint32            le      = getLastElemFragStore(frgStore) + 1;
  ReadStructp       rd      = new_ReadStruct();
  bool             *deleted = new bool [le];

  for (uint32 iid=0; iid<le; iid++)
    deleted[iid] = false;

  int    fragStreamBufferSize = 1048576;
  char  *fragStreamBuffer     = new char [fragStreamBufferSize];

  FragStreamHandle frgStream = openFragStream(frgStore, fragStreamBuffer, fragStreamBufferSize);
  resetFragStream(frgStream, STREAM_FROMSTART, STREAM_UNTILEND);

  fprintf(stderr, "Reading frags "F_UL" to "F_UL", remembering which are deleted.\n", fe, le-1);

  while (nextFragStream(frgStream, rd, FRAG_S_SOURCE)) {
    uint32  isDeleted = 0;
    uint32  fragIID   = 0;

    getIsDeleted_ReadStruct(rd, &isDeleted);
    getReadIndex_ReadStruct(rd, &fragIID);

    if (isDeleted)
      deleteLinksForFrag(fragIID);
  }

  delete_ReadStruct(rd);
  closeFragStream(frgStream);
  closeFragStore(frgStore);
}
