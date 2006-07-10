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

#include <Python.h>
#include "halign.h"

H_Alignment_t aln = { 0, 0, 0, 0, 0, 0, 0, 0};

static PyObject *
spam_halignStart(PyObject *self, PyObject *args) {
  char *seq1 = 0L;
  char *seq2 = 0L;

  if (!PyArg_ParseTuple(args, "ss", &seq1, &seq2))
    return(NULL);

  halignStart(seq1, seq2, &aln);

  Py_INCREF(Py_None);
  return(Py_None);
}
  
static PyObject *
spam_halignDedash( PyObject *self, PyObject *args) {
  int bgn1=0, bgn2=0, len1=0, len2=0, nmat=0;

  if (iterateUngappedAlignSharpEnds(&aln, bgn1, bgn2, len1, len2, nmat))
    return(Py_BuildValue("(iiiii)", bgn1, bgn2, len1, len2, nmat));

  Py_INCREF(Py_None);
  return(Py_None);
}

static
PyMethodDef
registration_table[] = {
  {"halignStart",  spam_halignStart,  METH_VARARGS, "initialize halign"},
  {"halignDedash", spam_halignDedash, METH_VARARGS, "dedashed subalignment"},
  {NULL, NULL, 0, NULL}
};

extern "C" 
void inithalign(void) {
  Py_InitModule("halign", registration_table);
}
