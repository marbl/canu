
//  Contructs a search tree from a matchList
class matchTree {
public:
  matchTree(matchList *L, u32bit side);
  ~matchTree() {
    dict_free_nodes(_tree);
    dict_free(_tree);
  };

  dict_t        *_tree;
  dict_load_t    _load;
};





matchTree::matchTree(matchList *L, u32bit side) {

  //  Construct a list of pointers to the matchList data
  //
  //  kazlib was modified to be qsort() compatible and so it passes a
  //  pointer to whatever it is sorting.  Since kazlib operates on
  //  pointers anyway, this means that it passes the compare function
  //  a pointer to a pointer to the object.
  //
  //  Which really fails in this case.  We have a list of pointers to
  //  objects that we sort, then want to load.
  //
  //  Uhhh, no, this is correct.  We give kazlib a pointer to the
  //  object, it gives the compare function a pointer to that pointer.
  //
  //  qsort() below sorts pointers to objects, and does the same.

  match_t  **matchPointers = new match_t * [L->_matchesLen];
  for (u32bit i=0; i<L->_matchesLen; i++)
    matchPointers[i] = L->_matches + i;

  //  Choose a comparison function based on the side we want

  int (*sortMatches)(const void *, const void *) = sortMatches1;
  if (side == 1)
    sortMatches = sortMatches2;
  
  //  Sort

  qsort(matchPointers, L->_matchesLen, sizeof(match_t *), sortMatches);

  //  Load the tree (use DICTCOUNT_T_MAX for max nodes)

  _tree = dict_create(L->_matchesLen, sortMatches);
  dict_allow_dupes(_tree);

  dict_load_begin(&_load, _tree);

  for (u32bit i=0; i<L->_matchesLen; i++) {
#if 0
    fprintf(stderr, "Load "u32bitFMT" "u32bitFMT" "u32bitFMT" -- "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
            matchPointers[i]->iid1, matchPointers[i]->pos1, matchPointers[i]->len1,
            matchPointers[i]->iid2, matchPointers[i]->pos2, matchPointers[i]->len2);
#endif

    dnode_t   *node = (dnode_t *)malloc(sizeof(dnode_t));
    dnode_init(node, 0L);
    dict_load_next(&_load, node, matchPointers[i]);
  }

  dict_load_end(&_load);

  //  Clean up
  delete [] matchPointers;
}
