
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received (LICENSE.txt) a copy of the GNU General Public 
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/
/*
 $Id: aug-cam.c,v 1.4 2005-03-22 19:48:32 jason_miller Exp $

 Usage: 

 Purpose:
*/

#include  "delcher.h"
#define  DEFAULT_COLOUR  "CF00CF"            // Magenta

int main(int argc, char * argv [])
{
  char  * colour = DEFAULT_COLOUR;
  int  ct, lo, hi;
  int  last_layout = 0;
  char  line [1000], * p;
  char  * unique_tag = "default";
#ifdef NEVER
  int num, last_category = -1;
  char ch;
  FILE  * fp = NULL;
#endif
  
  if  (argc < 1)
    {
      fprintf (stderr, "USAGE:  %s <colour-string> <unique-tag>\n", argv [0]);
      fprintf (stderr,
	       "  Adds repeat annotations (from celsim .cms file) onto end\n"
	       "  of <cam-file>.  Output is to  stdout .\n"
               "  unique-tag is used to create tag for celamy.\n");
      exit (-1);
    }
  
  if  (argc >= 2)
    colour = argv [1];
  if  (argc >= 3)
    unique_tag = argv [2];
  
  
  
#ifdef NEVER
  fp = File_Open (argv [1], "r");
  while  (fgets (line, 1000, fp) != NULL)
    {
      fputs (line, stdout);
      
      //  Assume Clark's format exactly
      ct = sscanf (line, "%d: %c", & num, & ch);
      if  (ct != 2)
	{
	  fprintf (stderr, "Bad cam line (skipping):  %s", line);
	  continue;
	}
      if  (isalpha (ch))
	last_category = num;
      else
	last_layout = num;
    }
  fclose (fp);
#endif
  

  // color definition line
  printf ("8RPT%s: C%s T7 # %s\n", unique_tag, colour, unique_tag);
  
  while  (fgets (line, 1000, stdin) != NULL)
    {
      fputs (line, stderr);
      p = strstr (line, " at ");
      if  (p == NULL)
	{
	  fprintf (stderr, "Bad cms line (skipping):  %s", line);
	  continue;
	}
      ct = sscanf (p + 4, "%d-%d", & lo, & hi);
      if  (ct != 2)
	{
	  fprintf (stderr, "Bad cms line (skipping):  %s", line);
	  continue;
	}

      p = strstr (line, "#");
      
      printf ("9RPT%s%04d: %d A8RPT%s %d R0 %s", unique_tag,
              ++ last_layout, lo, unique_tag, hi,
              (p == NULL) ? "\n" : p);
    }
  
  return 0;
}
