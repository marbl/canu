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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <Python.h>

#include "GF_ALN_local.h"

static Local_Overlap *desc = NULL;


//  This should be in the "library" not in the client.  Sigh.

void syntenicSegments(char const * const Aseq, int const Astart, int const Astop,
                      char const * const Bseq, int const Bstart, int const Bstop,
                      double const erate) {

  desc = NULL; // In case an early exit happens!

  //  Key data types ("Local_Segment" and "Local_Overlap") are defined
  //  in "CA_ALN_local.h"


  assert(Astop >= Astart);
  assert(Bstop >= Bstart);

  // Step 1: get local segments:
  char const * const Ausable = Aseq + Astart;
  char const * const Busable = Bseq + Bstart;

  int const alen = Astop - Astart;
  int const blen = Bstop - Bstart;

  int NumSegs = 0; /* number of local matches returned */

  Local_Segment *local_results = Find_Local_Segments(Ausable,         /* sequence A */
                                                     alen,
                                                     Busable, 
                                                     blen,
                                                     LOCAL_FORW,      /* whether to compute a forward search , reverse, or both */
                                                     16,              /* minimum length of a reportable match */
                                                     erate,           /* maximum error for a match to be returned */
                                                     &NumSegs);       /* number of local matches returned */
  
  if(NumSegs==0)
    return;
  

  // Step 2: get a chain of local segments:

  Local_Overlap *Ov = Find_Local_Overlap(alen,          /* length of sequence A */
                                         blen,          /* length of sequence B */
                                         0,             /* comp==0 -> fwd orientation */
                                         0,             /* nextbest==0 -> find best overlap*/
                                         local_results, /* the input set of local segments */
                                         NumSegs,       /* number of input local segments */
                                         20 - 6,        /* shortest "overlap" to report" */
                                         1.0);          /* fraction of overlap not in a match -- needs to be large to allow substantial mismatches */

  if(Ov == NULL)
    return;

          
  // Step 3 (optional): 
  //
  // a) fix the chain of segments so that the segments don't overlap.
  // It must be a 1-1 mapping. (can either trim or delete segments--or
  // leave them completely alone)
  //
  // b) construct an alignment "trace" 
  //
  // The "trace" is the standard "AS" encoding of an alignment.

  if(Ov != NULL) {

    // coordinate munge between Gene's local aligner and
    // DP_Compare()-related routines coordinates from Find_Local
    // routines will be one off from those expected by the trace
    // routines, so adjust them!

    for(int i=0;i<=Ov->num_pieces;i++){
      if(i<Ov->num_pieces){
        Ov->chain[i].piece.abpos++;
        Ov->chain[i].piece.bbpos++;
        Ov->chain[i].piece.aepos++;
        Ov->chain[i].piece.bepos++;
      }
    }

    //  AS_Local_Trace assumes string pointer one before start of string!

    if(Ov != NULL) {
      int *trace = AS_Local_Trace(Ov,Ausable-1,Busable-1);

      if(trace == NULL)
        fprintf(stderr,"EXCEPTION Ov=%p trace=NULL\n", Ov);
    }

    for(int i=0;i<=Ov->num_pieces;i++){
      if(i<Ov->num_pieces){
        Ov->chain[i].piece.abpos--;
        Ov->chain[i].piece.bbpos--;
        Ov->chain[i].piece.aepos--;
        Ov->chain[i].piece.bepos--;
      }
    }
  }

  if(Ov != NULL)
    Ov->next = 0;

  desc = Ov;
}




int iterate_Local_Overlap(int &seg_abpos, int &seg_bbpos,
                          int &seg_alen,  int &seg_blen,
                          double &seg_error) {

  if (desc == NULL)
    return(0);

  Local_Chain *chain = desc->chain;

  assert(NULL != desc->chain);

  for(; 0 <= desc->next && desc->next < desc->num_pieces; ) {
    int the_piece = (desc->next)++;
      
    Local_Segment *seg = &(chain[the_piece].piece);

    assert(NULL != seg);
    assert(!chain[the_piece].reversed);
      
    // Set the return data

    seg_abpos = seg->abpos;
    seg_alen  = seg->aepos - seg->abpos;
    seg_bbpos = seg->bbpos;
    seg_blen  = seg->bepos - seg->bbpos;
    seg_error = seg->error;

    // Skip over the "deleted in-place" segments.
    if((seg->aepos <= seg->abpos)&&(seg->bepos <= seg->bbpos)) 
      continue;
      
    // the data is valid
    return(1);
  }

  //  Nothing left.
  return(0);
}




static PyObject *
spam_syntenicSegments(PyObject *self, PyObject *args) {
  char *Aseq   = "undefined";
  int   Astart = -1;
  int   Astop  = -1; // substring of Aseq
  char *Bseq   = "undefined";
  int   Bstart = -1;
  int   Bstop  = -1; // substring of Bseq
  double erate  = 1.0 / 3.0;

  PyObject *py_outfile = NULL;

  if (!PyArg_ParseTuple(args, "Osiisiif", &py_outfile, &Aseq, &Astart, &Astop, &Bseq, &Bstart, &Bstop, &erate))
    return NULL;

  try {
    syntenicSegments(Aseq, Astart, Astop, // substring of Aseq
                     Bseq, Bstart, Bstop, // substring of Bseq
                     erate);
  } catch (...) {
    PyErr_SetString(PyExc_RuntimeError,"sytenicSegments failed");
    return(Py_None);
  }

  Py_INCREF(Py_None);  // This is a module function returning void.
  return(Py_None);
}



static PyObject *
spam_iterateSegments(PyObject *self, PyObject *args) {
  int     seg_bgn1  = 0;
  int     seg_bgn2  = 0;
  int     seg_len1  = 0;
  int     seg_len2  = 0;
  double  seg_error = 0.0;
    
  if (iterate_Local_Overlap(seg_bgn1, seg_bgn2, seg_len1, seg_len2, seg_error))
    return(Py_BuildValue("(iiiif)", seg_bgn1, seg_bgn2, seg_len1, seg_len2, seg_error));

  Py_INCREF(Py_None);  // This is a module function returning void.
  return(Py_None);
}



static
PyMethodDef
registration_table[] = {
  {"syntenicSegments", spam_syntenicSegments, METH_VARARGS, "Compute syntenic segments"},
  {"iterateSegments",  spam_iterateSegments,  METH_VARARGS, "Iterator returning syntenic segments"},
  {NULL, NULL, 0, NULL}
};


extern "C"
void initlocalAlignerInterface() {
  Py_InitModule("localAlignerInterface", registration_table);
}

