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


//
// For Linux:
// 
// [cmobarry@localhost GF_SYN]$ gmake localAlignerInterfacemodule.so
// g++ -shared -fpic -I /usr/include/python2.2 -I /home/cmobarry/work/trunk/cds/SYS/inc -o localAlignerInterfacemodule.so localAlignerInterfacemodule.cc localAlignerInterface.cc  GF_ALN_overlap.cc GF_ALN_local.cc GF_ALN_dpaligner.cc GF_ALN_qvaligner.cc

//
// python
//>>> import localAlignerInterface
//>>> Aseq = "aaaaaaaaaaaaaaaaaaaaaaaacccctttttttttttttttttttttttt"
//>>> Bseq = "aaaaaaaaaaaaaaaaaaaaaaccaagggtttttttttttttttttttttttt"
//>>> Astart = 0
//>>> Astop = len(Aseq)
//>>> Bstart = 0
//>>> Bstop = len(Bseq)
//>>> dir(Aseq)
//>>> dir(Astop)
//>>> 
//>>> localAlignerInterface.syntenicSegments(
//>>>     sys.stdout,
//>>>     Aseq, Astart, Astop,
//>>>     Bseq, Bstart, Bstop,
//>>>     0, "match1", "run1", "seq1", "seq2", 0, 0, 0, 0.3, 10, 10, 10, 10,
//>>> )
//>>>     
//
////////////////////////////////////////////////////////////////////////////////

#include <Python.h>
#include "localAlignerInterface.h"

#undef DEBUGGING

static PyObject *
spam_syntenicSegments( PyObject *self, PyObject *args) {
  char * Aseq = "undefined"; int Astart = -1; int Astop = -1; // substring of Aseq
  char * Bseq = "undefined"; int Bstart = -1; int Bstop = -1; // substring of Bseq
  float erate = 1./3.;

#ifdef DEBUGGING
  fprintf( stderr, "Before PyArg_ParseTuple in spam_sytenicSegments\n");
#endif
  
  PyObject *py_outfile = NULL;
  if (!PyArg_ParseTuple(args,
			"Osiisiif"
                        ,&py_outfile
			,&Aseq
			,&Astart, &Astop // substring of Aseq
			,&Bseq
			,&Bstart, &Bstop // substring of Bseq
                        ,&erate
			)) {
    return NULL;
  }
#ifdef DEBUGGING
  fprintf( stderr, "After PyArg_ParseTuple\n");
#endif
  
  FILE * outfile = PyFile_AsFile(py_outfile);
  
#ifdef DEBUGGING
  fprintf( stderr, "After PyArg_ParseTuple\n");
  fprintf( stderr,
           "outfile=%p\n"
	   "Aseq=%s\n"
	   "Bseq=%s\n"
	   "Astart=%d Astop=%d\n"
	   "Bstart=%d Bstop=%d\n"
           "erate=%f\n"
           ,outfile
	   ,Aseq
	   ,Bseq
	   ,Astart, Astop // substring of Aseq
	   ,Bstart, Bstop // substring of Bseq
           ,erate
	   );
#endif

  try {
    syntenicSegments
      (
       outfile,
       Aseq, Astart, Astop, // substring of Aseq
       Bseq, Bstart, Bstop, // substring of Bseq
       erate
       );
#if 0
  } catch(somethingReallyBad){
    Py_Exit(1);
#endif
  } catch(...) {
    PyErr_SetString(PyExc_RuntimeError,"sytenicSegments failed");
    return Py_None;
    //return NULL;
  }

#ifdef DEBUGGING
  fprintf( stderr, "After syntenicSegments\n");
#endif

  Py_INCREF(Py_None);  // This is a module function returning void.
  return Py_None;
}

#if 0
static PyObject *
spam_iterateSegments( PyObject *self, PyObject *args) {
  char * selfid = ".";
  char * parentid = ".";
  char * a_uid = "undefined";
  char * b_uid = "undefined";
  int a_scaf_pos=-1; // Lowest position in string "A" for the alignment.
  int b_scaf_pos=-1; // Lowest position in string "B" for the alignment.
  int alen = -1;     // Length of the alignment box for string "A".
  int blen = -1;     // Length of the alignment box for string "B".
  int sequence_reversed_match=-1;  /* bool */
  int a_lft_seed_len=-1;
  int b_lft_seed_len=-1;
  int a_rht_seed_len=-1;
  int b_rht_seed_len=-1;

#ifdef DEBUGGING
  fprintf( stderr, "Before PyArg_ParseTuple in spam_iterateSegments\n");
#endif
  
  PyObject *py_outfile = NULL;
  if (!PyArg_ParseTuple(args,
			"Ossssiiiiiiiii"
                        ,&py_outfile
                        ,&selfid
                        ,&parentid
			,&a_uid
			,&b_uid
			,&a_scaf_pos // Lowest position in string "A" for the alignment.
			,&b_scaf_pos // Lowest position in string "B" for the alignment.
			,&alen       // Length of the alignment for string "A".
			,&blen       // Length of the alignment for string "B".
			,&sequence_reversed_match /* bool */
			,&a_lft_seed_len
			,&b_lft_seed_len
			,&a_rht_seed_len
			,&b_rht_seed_len
			)) {
    return NULL;
  }

  FILE * outfile = PyFile_AsFile(py_outfile);
  
#ifdef DEBUGGING
  fprintf( stderr, "After PyArg_ParseTuple\n");
  fprintf( stderr,
           "outfile=%p\n"
           "selfid=%s\n"
           "parentid=%s\n"
	   "a_uid=%s b_uid=%s\n"
	   "a_scaf_pos=%d b_scaf_pos=%d \n"
	   "a_len=%d b_len=%d \n"
           "sequence_reversed_match=%d\n"
	   "a_lft_seed_len=%d b_lft_seed_len=%d\n"
	   "a_rht_seed_len=%d b_lft_seed_len=%d\n"
           ,outfile
           ,selfid
           ,parentid
	   ,a_uid
	   ,b_uid
	   ,a_scaf_pos // Lowest position in string "A" for the alignment.
	   ,b_scaf_pos // Lowest position in string "B" for the alignment.
	   ,alen       // Length of the alignment for string "A".
	   ,blen        // Length of the alignment for string "B".
	   ,sequence_reversed_match /* bool */
	   ,a_lft_seed_len
	   ,b_lft_seed_len
	   ,a_rht_seed_len
	   ,b_rht_seed_len
	   );
#endif

  int iret = Print_Local_Overlap4(
                                  outfile,
                                  selfid,
                                  parentid,
                                  a_uid,
                                  b_uid,
                                  a_scaf_pos, // Lowest position in string "A" for the alignment.
                                  b_scaf_pos, // Lowest position in string "B" for the alignment.
                                  alen,       // Length of the alignment box for string "A".
                                  blen,       // Length of the alignment box for string "B".
                                  sequence_reversed_match, /* bool */
                                  a_lft_seed_len,
                                  b_lft_seed_len,
                                  a_rht_seed_len,
                                  b_rht_seed_len
                                  );
  
  Py_INCREF(Py_None);  // This is a module function returning void.
  return Py_None;
}
#endif

#if 1
static PyObject *
spam_iterateSegments( PyObject *self, PyObject *args) {
#ifdef DEBUGGING
  fprintf( stderr, "Before PyArg_ParseTuple\n");
#endif
#if 0
  if (!PyArg_ParseTuple(args,
			"i"
			,&first
			)) {
    return NULL;
  }
#endif
#ifdef DEBUGGING
  fprintf( stderr, "After PyArg_ParseTuple\n");
#endif

  {
    int seg_bgn1, seg_bgn2, seg_len1, seg_len2;
    float seg_error;
    int valid = iterate_Local_Overlap( seg_bgn1, seg_bgn2, seg_len1, seg_len2, seg_error);
    if(valid){ 
      // build a tuple here
      return Py_BuildValue("(iiiif)", seg_bgn1, seg_bgn2, seg_len1, seg_len2, seg_error);
    } else {
      Py_INCREF(Py_None);  // This is a module function returning void.
      return Py_None;
    }
  }
}
#endif


static PyObject *
spam_system( PyObject *self, PyObject *args) {
  char *command;
  int sts;
  
  if (!PyArg_ParseTuple(args, "s", &command))
    return NULL;
  sts = system(command);
  return Py_BuildValue("i", sts);
}


static PyObject*
spam_hello( PyObject *self, PyObject *args) {
  fprintf( stderr, "Before PyArg_ParseTuple\n");
  PyObject *py_outfile;
  if (!PyArg_ParseTuple(args,
			"O"
                        ,&py_outfile
                        )) {
        return NULL;
  }
  FILE * outfile = PyFile_AsFile(py_outfile);
  fprintf(stderr,"outfile = %p\n", outfile);
  fprintf(outfile,"Hello, world!\n");
  Py_INCREF(Py_None);  // This is a module function returning void.
  return Py_None;
}


static PyMethodDef registration_table[] = {
  {"hello", spam_hello, METH_VARARGS, "Hello Man"},
  {"system", spam_system, METH_VARARGS, "Execute a shell command"},
  {"syntenicSegments", spam_syntenicSegments, METH_VARARGS, "Compute syntenic segments"},
  {"iterateSegments", spam_iterateSegments, METH_VARARGS, "Iterator returning syntenic segments"},
  {NULL, NULL, 0, NULL} /* Sentinel... the end of the table */
};

#ifdef __cplusplus
extern "C" 
#endif
void initlocalAlignerInterface()
{
  (void) Py_InitModule("localAlignerInterface", registration_table);
}
