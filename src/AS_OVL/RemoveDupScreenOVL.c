
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
/*************************************************
* Module:  RemoveDupScreenOVL.c
* Description:
*   Read a multifasta file of screen elements; sort it into descending
*   order by size; and check for overlaps.  Any contained elements
*   are removed; dovetail overlapping elements are merged.  Output
*   is written as a multi-fasta file.  Input is from  stdin ;  output
*   is to  stdout .
* 
*    Programmer:  A. Delcher
*       Started:   28 Jan 2002  Adapted from simple removedups program
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: RemoveDupScreenOVL.c,v 1.1.1.1 2004-04-14 13:52:44 catmandew Exp $
 * $Revision: 1.1.1.1 $
*/

static char CM_ID[] = "$Id: RemoveDupScreenOVL.c,v 1.1.1.1 2004-04-14 13:52:44 catmandew Exp $";


#include  "AS_OVL_delcher.h"
#include  "AS_ALN_aligners.h"
#include  "AS_CGW_dataTypes.h"


const int  INCR_SIZE = 10000;
const int  MAX_LINE = 1000;
const int  FASTA_WIDTH = 60;
const int  MIN_OLAP_LEN = 40;
const double  OLAP_ERATE = 0.02;


typedef  struct
  {
   char  * header;
   char  * seq;
   int  len;
   int  keep;
  }  String_t;


int  Cmp
    (const void * A, const void * B)

//  Compare  String_t 's for qsort .

  {
   String_t  * X, * Y;

   X = (String_t *) A;
   Y = (String_t *) B;

   if  (X -> len > Y -> len)
       return  -1;
   else if  (X -> len < Y -> len)
       return  1;

   return  0;
  }


static char  * Check_File_Name = NULL;
    // Name of file to read strings from and see if
    // covered by the results of this program.
static int  Strictness_Level = 0;
    // Determines how good overlaps must be to be used
static int  Verbose = 0;



static void  Fasta_Print
    (FILE * fp, char * hdr, char * seq);
static Overlap *  Find_Overlap
    (char * seq1, char * seq2, ChunkOrientationType orientation, 
     int min_ahang, int max_ahang,
     double erate, double thresh, int minlen, CompareOptions what);
static void  Merge_Headers
    (char * * h1, char * h2, int len, int reverse);
static void  Parse_Command_Line
    (int argc, char * argv []);
static int  Read_String
    (FILE * fp, char * * T, int * Size, char header []);
static void  Show_Olap
    (const String_t * a, const String_t * b, Overlap * olap,
     ChunkOrientationType orient);
static void  Usage
    (char * command);



int  main
    (int argc, char * argv [])

  {
   String_t  * S = NULL;
   Overlap  * result, best_result;
   ChunkOrientationType  best_orient = AS_UNKNOWN;
   char  * buff = NULL;
   char  header [MAX_LINE];
   int  buff_size;
   int  best_i = 0, best_len = 0;
   int  i, j, ct = 0, kept = 0;

   Parse_Command_Line  (argc, argv);

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

   qsort (S, ct, sizeof (String_t), Cmp);

   for  (j = 1;  j < ct;  j ++)
     {
      int  min_ahang, max_ahang, min_olap;

      best_len = 0;
      
      for  (i = 0;  i < j;  i ++)
        {
         if  (! (S [i] . keep && S [j] . keep))
             continue;

         switch  (Strictness_Level)
           {
            case  0 :
              min_ahang = - (S [j] . len - MIN_OLAP_LEN);
              max_ahang = S [i] . len - MIN_OLAP_LEN;
              break;

            case  1 :
              if  (S [i] . len <= S [j] . len)
                  min_olap = S [i] . len / 2;
                else
                  min_olap = S [j] . len / 2;
              if  (min_olap < MIN_OLAP_LEN)
                  min_olap = MIN_OLAP_LEN;
              min_ahang = - (S [j] . len - min_olap);
              max_ahang = S [i] . len - min_olap;
              break;

            case  2 :
              if  (S [i] . len <= S [j] . len)
                  min_olap = S [i] . len - MIN_OLAP_LEN;
                else
                  min_olap = S [j] . len - MIN_OLAP_LEN;
              if  (min_olap < MIN_OLAP_LEN)
                  min_olap = MIN_OLAP_LEN;
              min_ahang = - (S [j] . len - min_olap);
              max_ahang = S [i] . len - min_olap;
              break;

            case  3 :
              if  (S [i] . len <= S [j] . len)
                  min_olap = S [i] . len;
                else
                  min_olap = S [j] . len;
              if  (min_olap < MIN_OLAP_LEN)
                  min_olap = MIN_OLAP_LEN;
              min_ahang = - (S [j] . len - min_olap);
              max_ahang = S [i] . len - min_olap;
              break;

            default :
              fprintf (stderr, "ERROR:  Bad strictness level = %d\n",
                       Strictness_Level);
              exit (EXIT_FAILURE);
           }

         result = Find_Overlap
                     (S [i] . seq, S [j] . seq, AB_AB,
                      min_ahang, max_ahang,
                      OLAP_ERATE, CGW_DP_THRESH, MIN_OLAP_LEN,
                      AS_FIND_ALIGN);

         if  (result != NULL && result -> length > best_len)
             {
              best_len = result -> length;
              best_i = i;
              best_result = (* result);
              best_orient = AB_AB;
             }

         result = Find_Overlap
                     (S [i] . seq, S [j] . seq, AB_BA,
                      min_ahang, max_ahang,
                      OLAP_ERATE, CGW_DP_THRESH, MIN_OLAP_LEN,
                      AS_FIND_ALIGN);

         if  (result != NULL && result -> length > best_len)
             {
              best_len = result -> length;
              best_i = i;
              best_result = (* result);
              best_orient = AB_BA;
             }
        }

      if  (best_len > 0)
          {
           if  (Verbose > 0)
               {
                Show_Olap (S + best_i, S + j, & best_result, best_orient);
               }
           if  (best_result . begpos >= 0 && best_result . endpos <= 0)
               {  // j is contained in best_i
                S [j] . keep = FALSE;
                continue;
               }
           if  (best_result . begpos <= 0 && best_result . endpos >= 0)
               {  // best_i is contained in j  should be rare since best_i is longer
                S [best_i] = S [j];
                S [j] . keep = FALSE;
                continue;
               }
           if  (best_result . endpos > 0)
               {
                S [best_i] . len += best_result . endpos;
                S [best_i] . seq = (char *) Safe_realloc (S [best_i] . seq,
                                       1 + S [best_i] . len);
                if  (best_orient == AB_BA)
                    Complement_Seq (S [j] . seq);
                strcat (S [best_i] . seq,
                        S [j] . seq + S [j] . len - best_result . endpos);
                S [j] . keep = FALSE;
                Merge_Headers (& (S [best_i] . header), S [j] . header,
                               S [best_i] . len, FALSE);
                continue;
               }
           if  (best_result . begpos < 0)
               {
                S [best_i] . len -= best_result . begpos;
                S [best_i] . seq = (char *) Safe_realloc (S [best_i] . seq,
                                       1 + S [best_i] . len);
                memmove (S [best_i] . seq - best_result . begpos, S [best_i] . seq,
                         1 + S [best_i] . len + best_result . begpos);
                if  (best_orient == AB_BA)
                    Complement_Seq (S [j] . seq);
                strncpy (S [best_i] . seq, S [j] . seq, - best_result . begpos);
                S [j] . keep = FALSE;
                Merge_Headers (& (S [best_i] . header),
                               S [j] . header, S [best_i] . len, TRUE);
                continue;
               }
          }
     }

   if  (Check_File_Name != NULL)
       {
        int  check_id = 0;
        FILE  * fp = fopen (Check_File_Name, "r");

        while  (Read_String (fp, & buff, & buff_size, header))
          {
           int  len, max_cover;

           check_id ++;
           len = strlen (buff);

           max_cover = 0;
           for  (i = 0;  i < ct;  i ++)
             if  (S [i] . keep)
                 {
                  result = Find_Overlap
                              (buff, S [i] . seq, AB_AB,
                               - (S [i] . len - len) - 2,
                               2,
                               OLAP_ERATE + 0.01, CGW_DP_THRESH, MIN_OLAP_LEN,
                               AS_FIND_ALIGN);
                  if  (result != NULL && result -> begpos <= 2
                         && result -> endpos >= -2)
                      {
                       printf ("Check #%3d covered by screen #%d %s\n",
                               check_id, i, S [i] . header);
                       break;
                      }
                  if  (result != NULL && result -> length > max_cover)
                      max_cover = result -> length;

                  result = Find_Overlap
                              (buff, S [i] . seq, AB_BA,
                               - (S [i] . len - len) - 2,
                               2,
                               OLAP_ERATE + 0.01, CGW_DP_THRESH, MIN_OLAP_LEN,
                               AS_FIND_ALIGN);
                  if  (result != NULL && result -> begpos <= 2
                         && result -> endpos >= -2)
                      {
                       printf ("Check #%3d covered by rev screen #%d %s\n",
                               check_id, i, S [i] . header);
                       break;
                      }
                  if  (result != NULL && result -> length > max_cover)
                      max_cover = result -> length;
                 }

           if  (i == ct)
               {
                printf ("Check #%3d not covered  len = %d  max_cover = %d\n",
                        check_id, len, max_cover);
                printf ("  %s\n", header);
                
               }
          }
       }

   // Sort again, lengths have changed

   qsort (S, ct, sizeof (String_t), Cmp);

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



static Overlap *  Find_Overlap
    (char * seq1, char * seq2, ChunkOrientationType orientation, 
     int min_ahang, int max_ahang,
     double erate, double thresh, int minlen, CompareOptions what)

//  Same as  OverlapSequences  in CGW

{
  Overlap *omesg;
  int flip = 0;
  
  // if the orientation is BA_AB or BA_BA, we need to reverse complement the first contig
  if (orientation == BA_AB || orientation == BA_BA)
	Complement_Seq( seq1 );
  
  // if the orientation is AB_BA or BA_BA, we need to set the flip variable for the second contig
  if (orientation == AB_BA || orientation == BA_BA)
	flip = 1;

  // min_ahang and end are essentially bounds on the a-hang
  omesg = DP_Compare(seq1, seq2,
					 min_ahang, max_ahang, (int) flip,
					 erate, thresh, minlen,
					 what);

  // return seq1 to its original state
  if (orientation == BA_AB || orientation == BA_BA)
	Complement_Seq( seq1 );
  
  // omesg->begpos is the a-hang, omesg->endpos is the b-hang
  return omesg;
}



static void  Merge_Headers
    (char * * h1, char * h2, int len, int reverse)

//  Merge string header  h2  into  (* h1)  allocating extra space
//  if needed.  Set the length field of the result to  len .
//  If  reverse  is true, put the  h2  info in front of the  h1  info;
//  otherwise, put  h1  in front of  h2 .

  {
   char  tag1 [MAX_LINE];
   char  tag2 [MAX_LINE];
   char  range1 [MAX_LINE];
   char  range2 [MAX_LINE];
   char  new_header [MAX_LINE];

   sscanf ((* h1), "%s %s", tag1, range1);
   sscanf (h2, "%s %s", tag2, range2);

   if  (reverse)
       sprintf (new_header, " iids %s,%s  len = %d",
                range2, range1, len);
     else
       sprintf (new_header, " iids %s,%s  len = %d",
                range1, range2, len);
   
   (* h1) = (char *) Safe_realloc ((* h1), 1 + strlen (new_header));
   strcpy ((* h1), new_header);

   return;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
    int  ch, errflg = FALSE;
   char  * p;

   optarg = NULL;

   while  (! errflg
             && ((ch = getopt (argc, argv, "C:L:v:")) != EOF))
     switch  (ch)
       {
        case  'C' :
          Check_File_Name = optarg;
          break;
         
        case  'L' :
          Strictness_Level = (int) strtol (optarg, & p, 10);
          if  (p == optarg)
              {
               fprintf (stderr, "ERROR:  Illegal strictness level \"%s\"\n",
                        optarg);
               errflg = TRUE;
              }
          break;

        case  'v' :
          Verbose = (int) strtol (optarg, & p, 10);
          if  (p == optarg)
              {
               fprintf (stderr, "ERROR:  Illegal verbose level \"%s\"\n",
                        optarg);
               errflg = TRUE;
              }
          break;

        case  '?' :
          fprintf (stderr, "Unrecognized option -%c\n", optopt);

        default :
          errflg = TRUE;
       }

   if  (errflg || optind != argc)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

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



static void  Show_Olap
    (const String_t * a, const String_t * b, Overlap * olap,
     ChunkOrientationType orient)

//  Show in  stdout a picture of the overlap between  a  and  b
//  that's indicated by  olap .   orient indicates if  b  has
//  been reverse-complemented

  {
   printf ("\n");
   if  (olap -> begpos >= 0 && olap -> endpos <= 0)
       {
        // a contains b
        printf ("------------------------------------> %s\n",
                a -> header);
        printf ("%5d %c-----------------------%c %5d %s\n",
                olap -> begpos,
                orient == AB_AB ? '-' : '<',
                orient == AB_AB ? '>' : '-',
                - olap -> endpos, b -> header);
       }
   else if  (olap -> begpos <= 0 && olap -> endpos >= 0)
       {
        // b contains 1
        printf ("%5d ------------------------> %5d %s\n",
                - olap -> begpos, olap -> endpos, a -> header);
        printf ("%c-----------------------------------%c %s\n",
                orient == AB_AB ? '-' : '<',
                orient == AB_AB ? '>' : '-',
                b -> header);
       }
   else if  (olap -> begpos > 0)
       {
        // a -> b dovetail
        printf ("------------------------------> %5d %s\n",
                olap -> endpos, a -> header);
        printf ("%5d %c-----------------------------%c %s\n",
                olap -> begpos,
                orient == AB_AB ? '-' : '<',
                orient == AB_AB ? '>' : '-',
                b -> header);
       }
     else
       {
        // b -> a dovetail
        printf ("%5d ------------------------------> %s\n",
                - olap -> begpos, a -> header);
        printf ("%c-----------------------------%c %5d %s\n",
                orient == AB_AB ? '-' : '<',
                orient == AB_AB ? '>' : '-',
                - olap -> endpos, b -> header);
       }

   return;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
           "USAGE:  %s [-v VerboseLevel]\n"
           "\n"
           "Combines overlapping and contained screen elements\n"
           "Input is from  stdin ; output is to  stdout \n"
           "\n"
           "Options:\n"
           "-L   set strictness level of overlaps\n"
           "     0  any overlap >= MIN_OLAP_LEN (default)\n"
           "     1  any overlap >= 1/2 of shorter string\n"
           "     2  overlaps within MIN_OLAP_LEN of being containments\n"
           "     3  strict containments\n"
           "-v   set level of verbose output\n",
           command);

   return;
  }



