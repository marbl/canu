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

#ifndef OVERLAP_SPANTREE_H
#define OVERLAP_SPANTREE_H

#include "overlap.H"

class spanTree {
public:
  spanTree() {
    _tree = dict_create(DICTCOUNT_T_MAX, spanCompare);
  };
  ~spanTree() {
    dict_free_nodes(_tree);
    dict_free(_tree);
  };

  void addNewSpan(u32bit iid, u32bit len) {
    span_t       *span = new span_t(iid, 0, len);
    dict_alloc_insert(_tree, span, 0L);
  };

  u32bit size(void) {
    return((u32bit)dict_count(_tree));
  };

  void addMatch(match_t *match, u32bit side, u32bit color) {

    //  Query the tree for the first match before this position.
    //  We're guaranteed to find one before, since the tree was
    //  initialized with a span for the whole sequence.
    //
    span_t   *span = 0L;
    u32bit    beg  = 0;
    u32bit    end  = 0;

    if (side == 0) {
      span = new span_t(match->iid1, match->pos1, match->pos1 + match->len1);
      beg  = match->pos1;
      end  = match->pos1 + match->len1;
    } else {
      span = new span_t(match->iid2, match->pos2, match->pos2 + match->len2);
      beg  = match->pos2;
      end  = match->pos2 + match->len2;
    }

    dnode_t  *node = dict_upper_bound(_tree, span);
    delete span;
    span = (span_t *)dnode_getkey(node);

    //  We need to split the span pointed to by node, iterate through
    //  all the spans, and split the last one.

    if (span->_beg != beg) {
      span_t  *l = 0L;
      span_t  *r = 0L;
      span->split(beg, l, r);

      //  Kill this node, insert the new ones

      dict_delete(_tree, node);
      dnode_destroy(node);
      dict_alloc_insert(_tree, l, 0L);
      dict_alloc_insert(_tree, r, 0L);

      delete span;
      span = r;

      //  Argh!  Now find the node we just inserted...

      node = dict_lookup(_tree, r);
    }

    //  Until we hit the last span, add the match to the span

    while (span->_end < end) {
      span->addMatch(match->matchiid, color);
      node = dict_next(_tree, node);
      span = (span_t *)dnode_getkey(node);
    }

    //  We're at the last span, so split it like the beginning

    if (span->_end != end) {
      span_t  *l = 0L;
      span_t  *r = 0L;
      span->split(end, l, r);

      //  Kill this node, insert the new ones

      dict_delete(_tree, node);
      dnode_destroy(node);
      dict_alloc_insert(_tree, l, 0L);
      dict_alloc_insert(_tree, r, 0L);

      delete span;
      span = l;
    }

    //  FInally, add the match to the last span

    span->addMatch(match->matchiid, color);
  };

  dict_t        *_tree;
  dict_load_t    _load;
};

#endif // OVERLAP_SPANTREE_H
