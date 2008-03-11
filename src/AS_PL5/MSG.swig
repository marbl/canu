%module "AS::MSG"
%{
#include "AS_MSG/AS_MSG_pmesg.h"
#include "MSG.h"
%}

/* AS pmsg definitions
 */
%include pmesg.h

/* Parser definitions
 */
%include parser.h

/* 3 letters acronymes of struct members
 * TODO: This should really be moved to Parser.pm.
 */
%include acronyms.h
