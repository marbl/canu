
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
#include  "AS_OVL_delcher.h"


const int  INCR_SIZE = 10000;
const int  MAX_LINE = 1000;
const int  FASTA_WIDTH = 60;


typedef struct
  {
   char  * header;
   char  * seq;
   int  len;
   int  keep;
  }  String_t;



static void  Fasta_Print
    (FILE * fp, char * hdr, char * seq);
static int  Read_String
    (FILE * fp, char * * T, int * Size, char header []);



int  main
    ()

  {
   String_t  * S = NULL;
   char  * buff = NULL;
   char  header [MAX_LINE];
   int  buff_size;
   int  i, j, ct = 0, kept = 0;

   buff_size = 10000;
   buff = (char *) Safe_malloc (buff_size);

   while  (Read_String (stdin, & buff, & buff_size, header))
     {
      S = (String_t *) Safe_realloc (S, (ct + 1) * sizeof (String_t));
      S [ct] . header = strdup (header);
      S [ct] . seq = strdup (buff);
      S [ct] . len = strlen (buff);
      S [ct] . keep = TRUE;
      ct ++;
     }

   for  (i = 0;  i < ct - 1;  i ++)
     for  (j = i + 1;  j < ct;  j ++)
       if  (S [i] . len <= S [j] . len)
           {
            if  (S [i] . keep && S [j] . keep
                   && strstr (S [j] . seq, S [i] . seq) != NULL)
                S [i] . keep = FALSE;
           }
         else
           {
            if  (S [i] . keep && S [j] . keep
                   && strstr (S [i] . seq, S [j] . seq) != NULL)
                S [j] . keep = FALSE;
           }

   for  (i = 0;  i < ct;  i ++)
     if  (S [i] . keep)
         {
          Fasta_Print (stdout, S [i] . header, S [i] . seq);
          kept ++;
         }

   fprintf (stderr, "Kept %d of %d sequences\n", kept, ct);

   return  0;
  }



static void  Fasta_Print
    (FILE * fp, char * hdr, char * seq)

//  Print to  fp  the sequence  seq  in FASTA format
//  with  hdr  as the label on the header line.

  {
   char  * p;
   int  ct = 0;

   fprintf (fp, ">%s\n", hdr);

   for  (p = seq;  * p != '\0';  p ++)
     {
      if  (ct == FASTA_WIDTH)
          {
           fputc ('\n', fp);
           ct = 0;
          }
      fputc (* p, fp);
      ct ++;
     }
   if  (ct > 0)
       fputc ('\n', fp);

   return;
  }



static int  Read_String
    (FILE * fp, char * * T, int * Size, char header [])

/* Read next string from  fp  (assuming FASTA format) into  (* T) [1 ..]
*  which has  (* Size)  characters.  Allocate extra memory if needed
*  and adjust  (* Size)  accordingly.  Return  TRUE  if successful,  FALSE
*  otherwise (e.g., EOF).  Set  header  to the contents of the FASTA
*  header line. */

  {
   char  Line [MAX_LINE];
   long int  Len;
   int  Ch, Ct;

   while  ((Ch = fgetc (fp)) != EOF && Ch != '>')
     ;

   if  (Ch == EOF)
       return  FALSE;

   fgets (Line, MAX_LINE, fp);
   Len = strlen (Line);
   assert (Len > 0 && Line [Len - 1] == '\n');
   Line [Len - 1] = '\0';
   strcpy (header, Line);


   Ct = 0;
   Len = 0;
   while  ((Ch = fgetc (fp)) != EOF && Ch != '>')
     {
      if  (isspace (Ch))
          continue;

      Ct ++;

      if  (Ct >= (* Size) - 1)
          {
           (* Size) += INCR_SIZE;
           (* T) = (char *) Safe_realloc ((* T), (* Size));
          }
      Ch = tolower (Ch);
      if  (! isalpha (Ch))
          {
           fprintf (stderr, "Unexpected character `%c\' in string %s\n",
                                 Ch, header);
           Ch = 'x';
          }
      (* T) [Len ++] = Ch;
     }

   (* T) [Len] = '\0';
   if  (Ch == '>')
       ungetc (Ch, fp);

   return  TRUE;
  }



