
static char *rcsid = "$Id:  $";

#include "outputFalcon.H"

#include "AS_UTL_reverseComplement.H"


//  The falcon consensus format:
//
//  name sequence
//  read sequence
//  read sequence
//  read sequence
//  read sequence
//  + +            #  generate consensus for the 'name' sequence using 'read' sequences
//  ...
//  - -            #  To end processing
//


void
outputFalcon(gkStore      *gkpStore,
             tgTig        *tig,
             bool          trimToAlign,
             FILE         *F,
             gkReadData   *readData) {

  gkpStore->gkStore_loadReadData(tig->tigID(), readData);

  fprintf(F, "read"F_U32" %s\n", tig->tigID(), readData->gkReadData_getSequence());

  for (uint32 cc=0; cc<tig->numberOfChildren(); cc++) {
    tgPosition  *child = tig->getChild(cc);

    gkpStore->gkStore_loadReadData(child->ident(), readData);

    if (child->isReverse())
      reverseComplementSequence(readData->gkReadData_getSequence(),
                                readData->gkReadData_getRead()->gkRead_sequenceLength());

    //  For debugging/testing, skip one orientation of overlap.
    //
    //if (child->isReverse() == false)
    //  continue;
    //if (child->isReverse() == true)
    //  continue;

    //  Trim the read to the aligned bit
    char   *seq = readData->gkReadData_getSequence();

    if (trimToAlign) {
      seq += child->_askip;
      seq[ readData->gkReadData_getRead()->gkRead_sequenceLength() - child->_askip - child->_bskip ] = 0;
    }

    fprintf(F, "data"F_U32" %s\n", tig->getChild(cc)->ident(), seq);
  }

  fprintf(F, "+ +\n");
}

