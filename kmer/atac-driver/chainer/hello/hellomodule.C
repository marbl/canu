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
// g++ -shared -fpic -I /usr/include/python2.2 -o hellomodule.so hellomodule.cc

//// For Compaq
// cxx -shared -I/usr/local/ir/Python-2.2.2 -I/usr/local/ir/Python-2.2.2/Include -o hellomodule.so hellomodule.cc

//
// python
//>>> import hello
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
#include <malloc.h>

/////////////////////////////////////////////////////////////////////////////////////

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
  printf("Hello, world!" __TIME__ "\n");
  Py_INCREF(Py_None);  // This is a module function returning void.
  return Py_None;
}

static PyMethodDef registration_table[] = {
  {"hello", spam_hello, METH_VARARGS, "Hello Man"},
  {"system", spam_system, METH_VARARGS, "Execute a shell command"},
  {NULL, NULL, 0, NULL} /* Sentinel... the end of the table */
};

#ifdef __cplusplus
extern "C" 
#endif
void inithello()
{
  (void) Py_InitModule("hello", registration_table);
}
