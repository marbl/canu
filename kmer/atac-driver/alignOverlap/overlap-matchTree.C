// This file is part of A2Amapper.
// Copyright (c) 2005 J. Craig Venter Institute
// Author: Brian Walenz
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "overlap.H"

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
    dnode_t   *node = (dnode_t *)malloc(sizeof(dnode_t));
    dnode_init(node, 0L);
    dict_load_next(&_load, node, matchPointers[i]);
  }

  dict_load_end(&_load);

  //  Clean up
  delete [] matchPointers;
}
