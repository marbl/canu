// This file is part of A2Amapper.
// Copyright (c) 2004 Applera Corporation
// Author: Clark Mobarry
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


#ifndef __CYGWIN__
#include <libgen.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>

// using namespace std;
#include "heavyChainsMod.h"

// #include "libbritypes.h"

extern "C" {
  void    *construct(char *options);
  void     destruct(void *handle);
  void     addHit(void *handle,
                  char    orientation,
                  u32bit  id1,
                  u32bit  pos1,
                  u32bit  len1,
                  u32bit  id2,
                  u32bit  pos2,
                  u32bit  len2,
                  u32bit  filled);
  void     filter(void *handle);
  void     output(void *handle, FILE *file);


  void* constructStats(char *options);
  void  destructStats(void *handle);
  void  addStats(void  *handle, void *sp);
  void  showStats(void *handle, FILE *file);

}

///////////////////////////////////////////

void*
construct(char *options) {
#if 0
  int argc = 1;
  char *argv[2];
  argv[0] = NULL;
  argv[1] = opts;
  void * handle = (void *)(new StrandPair(argc, argv));
#else
  void * handle = (void *)(new StrandPair(10,100000,100.));
#endif
  fprintf(stderr, "construct strand pair %p\n",handle);
  return handle;
}

void
destruct(void *handle) {
  fprintf(stderr, "destruct strand pair %p\n", handle);
  delete (StrandPair *)handle;
}

void
addHit(void   *handle,
       char    orientation,
       u32bit  id1,
       u32bit  pos1,
       u32bit  len1,
       u32bit  id2,
       u32bit  pos2,
       u32bit  len2,
       u32bit  filled) {
  ((StrandPair *)handle)->addHit(orientation, id1, pos1, len1, id2, pos2, len2, filled);
}

void
filter(void *handle) {
  ((StrandPair *)handle)->process();
}

void
output(void *handle, FILE *file) {
  long matchid = 0;
  matchid = ((StrandPair *)handle)->print(file,1,"","",matchid);
}

//////////////////////////////////////////////////////////

void*
constructStats(char *options) {
#if 0
  int argc = 1;
  char *argv[2];
  argv[0] = NULL;
  argv[1] = opts;
  void * handle = (void *)(new TheStats(argc, argv));
#else
  std::map<std::string,std::string> globals;
  void * handle = (void *)(new TheStats(globals));
#endif
  fprintf(stderr, "construct stats object %p\n", handle);
  return handle;
}

void
destructStats(void *handle) {
  fprintf(stderr, "destruct stats object %p\n", handle);
  delete (TheStats *)handle;
}

void
addStats(void   *handle, void *sp) {
  ((TheStats *)handle)->add((StrandPair *)sp);
}

void
showStats(void *handle, FILE *file) {
  ((TheStats *)handle)->print(file);
}

//////////////////////////////////////////////////////
