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

////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>

#undef DEBUGGING
#include "halign.h"

static H_Alignment_t* aln_ptr = NULL;

static PyObject *
spam_halignStart( PyObject *self, PyObject *args) {
  char *seq1, *seq2;
  int   len1=0, len2=0;
  int   offset1=0, offset2=0;

  if (!PyArg_ParseTuple(args,
			"ss"
			,&seq1
			,&seq2
			)) {
    return NULL;
  }

  len1 = strlen(seq1);
  len2 = strlen(seq2);
  // Sequence coordinates are base-based, starting from 0
  if(len1 > 0 && len2 >0 ) {
    halignStart(seq1+offset1, // This is the first base in the comparison.
                seq2+offset2,
                offset1, offset2,
                len1, len2,
                &aln_ptr);
  } else {
    if(aln_ptr != NULL){
      Free_align(aln_ptr); // Must call for each halign() but after printing output.
    }
    aln_ptr = NULL;
  }
  Py_INCREF(Py_None);  // This is a module function returning void.
  return Py_None;
}

  
static PyObject *
spam_halignDedash( PyObject *self, PyObject *args) {

  int bgn1, bgn2, len1, len2, nmat;
  int valid = iterateUngappedAlignSharpEnds(aln_ptr, bgn1, bgn2, len1, len2, nmat);
  if(valid){ 
    // build a tuple here
    return Py_BuildValue("(iiiii)", bgn1, bgn2, len1, len2, nmat);
  } else {
    Py_INCREF(Py_None);  // This is a module function returning void.
    return Py_None;
  }
}

static PyObject *
spam_Free_align( PyObject *self, PyObject *args) {
  if(aln_ptr != NULL){
    Free_align(aln_ptr); // Must call for each halign() but after printing output.
    aln_ptr = NULL;
  }
  Py_INCREF(Py_None);  // This is a module function returning void.
  return Py_None;
}

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
  {"halignStart", spam_halignStart, METH_VARARGS, "Initialize halign"},
  {"halignDedash", spam_halignDedash, METH_VARARGS, "dedashed subalignment"},
  {"halignFree", spam_Free_align, METH_VARARGS, "Free halign resources"},
  {NULL, NULL, 0, NULL} /* Sentinel... the end of the table */
};

#ifdef __cplusplus
extern "C" 
#endif
void inithalign()
{
  (void) Py_InitModule("halign", registration_table);
}
