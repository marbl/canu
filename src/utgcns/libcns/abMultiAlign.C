

uint32
abMultiAlign::getConsensus(abAbacus *ab,
                           vector<char> &bases,
                           vector<char> &quals) {
    uint32   length = length();

    bases.clear();
    quals.clear();

    abColumn   *col = GetColumn(columnStore, GetMANode(manodeStore, mid)->first);
    abBeadID   bid = col->call;

    for (int32 i=0; bid.isValid(); i++) {
      Bead *bead = ab->getBead(bid);

      bases.push_back(ab->getBase(bead->soffset));
      bases.push_back(ab->getQual(bead->soffset));

      bid = bead->next;
    }

    return(length);
  };
