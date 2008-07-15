
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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_MSG_pmesg_internal.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_fasta.h"
#include "AS_UTL_fileIO.h"
#include "AS_SDB_SequenceDB.h"

#define AS_SEQAN_MAX_HEADER_LENGTH     4
#define AS_SEQAN_MAX_LINE_LENGTH      80
#define AS_SEQAN_MAX_BUFFER_LENGTH  1024
#define AS_SEQAN_MAX_RESULT_LENGTH  AS_FRAG_MAX_LEN+AS_SEQAN_MAX_HEADER_LENGTH
#define AS_SEQAN_INPUT_NAME         "in.reads"
#define AS_SEQAN_RESULT             "temp.result"
#define AS_SEQAN_CNS                "temp.cns"

char *replace_str(char *str, char *orig, char *rep)
{
  static char buffer[AS_SEQAN_MAX_BUFFER_LENGTH];
  char *p;

  if(!(p = strstr(str, orig))) {
    return str;
  }

  int len = p-str;
  strncpy(buffer, str, len);
  buffer[p-str] = '\0';

  sprintf(buffer+len, "%s%s", rep, p+strlen(orig));

  return buffer;
}

void getFileName(char *prefix, char * wrkDir, char *fileName, char *result) {
   if (wrkDir != NULL) {
      // make sure we don't have wrkDir in both the prefix and the wrkDir
      char *wrkPrefix = replace_str(prefix, wrkDir, "");
      sprintf(result, "%s%s_%s", wrkDir, wrkPrefix, fileName);
   } else {
      sprintf(result, "%s_%s", prefix, fileName);
   }
}

char * readMultiLine(FILE* inFile) {
   char * result = (char *) safe_malloc(sizeof(char)*AS_SEQAN_MAX_BUFFER_LENGTH);
   int position = 0;
   int resultSize = AS_SEQAN_MAX_BUFFER_LENGTH;

   while (!feof(inFile)) {
      if ((position + AS_SEQAN_MAX_LINE_LENGTH) > resultSize) {
         resultSize += AS_SEQAN_MAX_BUFFER_LENGTH; // increase buffer size

         char *temp = safe_realloc(result, resultSize);
         assert(temp);
         result = temp;
      }
      fgets(result+position, AS_SEQAN_MAX_LINE_LENGTH, inFile);
      if (strcmp(result+position, "\n") == 0) {
fprintf(stderr, "Comparing %s to the end of line to see what I get\n", result+position);
         chomp(result);
         break;         
      }
      chomp(result);

      position = strlen(result);
   }

   return result;
}

char *getConsensus(char *inFile, char *seqAn, char * prefix, char * wrkDir) {
   char command[AS_SEQAN_MAX_BUFFER_LENGTH];
   char resultFile[AS_SEQAN_MAX_BUFFER_LENGTH];
   char cnsFile[AS_SEQAN_MAX_BUFFER_LENGTH];
   getFileName(prefix, wrkDir, AS_SEQAN_RESULT, resultFile);
   getFileName(prefix, wrkDir, AS_SEQAN_CNS, cnsFile);

   FILE *tempOut;
   sprintf(command, "%s -reads %s -outfile %s", seqAn, inFile, resultFile);
   assert(system(command) == 0);
   sprintf(command, "grep \"C:\" %s | awk '{print $2}' > %s", resultFile, cnsFile);
   assert(system(command) == 0);

   tempOut = fopen(cnsFile,"r");
   char * result = readMultiLine(tempOut); 
   fclose(tempOut);

   return result;
}

void updateRecord(IntUnitigMesg *ium_mesg, char * inFile, char *seqAn, char * prefix, char * wrkDir) {
   // update the consensus
   ium_mesg->consensus = getConsensus(inFile, seqAn, prefix, wrkDir);
   ium_mesg->length = strlen(ium_mesg->consensus);

   // update quality
   ium_mesg->quality = (char *) safe_malloc(sizeof(char) * ium_mesg->length+1);
   memset(ium_mesg->quality, '1', ium_mesg->length);
   ium_mesg->quality[ium_mesg->length] = '\0';

   //update read data
   int currRead = 0;
   char line[AS_SEQAN_MAX_RESULT_LENGTH];
   char resultFile[AS_SEQAN_MAX_BUFFER_LENGTH];
   getFileName(prefix, wrkDir, AS_SEQAN_RESULT, resultFile);

   FILE *tempOut;

   tempOut = fopen(resultFile,"r");
   // skip header of output
   while (!feof(tempOut)) {
      fgets(line, AS_SEQAN_MAX_RESULT_LENGTH, tempOut);
      if (strncasecmp(line, "typ", 3) == 0) {
         break;
      }
   }

   // now read alignments of each read
   for (currRead = 0; currRead < ium_mesg->num_frags; currRead++) {
      if (currRead != 0) {
         // read the typ: line
         fgets(line, AS_SEQAN_MAX_RESULT_LENGTH, tempOut);
      }
      // read the seq: line
      fgets(line, AS_SEQAN_MAX_RESULT_LENGTH, tempOut);
      // read the Pos line
      fgets(line, AS_SEQAN_MAX_RESULT_LENGTH, tempOut);
      chomp(line);
      CDS_COORD_t begin, end;
      sscanf(line,"Pos:"F_COORD","F_COORD,&begin,&end);
      ium_mesg->f_list[currRead].position.bgn = begin;
      ium_mesg->f_list[currRead].position.end = end;

      // read the dln line
      fgets(line, AS_SEQAN_MAX_RESULT_LENGTH, tempOut);
      chomp(line);
      sscanf(line,"dln:"F_S32, &ium_mesg->f_list[currRead].delta_length);

      fgets(line, AS_SEQAN_MAX_RESULT_LENGTH, tempOut);
      chomp(line);

      if (ium_mesg->f_list[currRead].delta_length > 0) {
         char *dlnStr = line+AS_SEQAN_MAX_HEADER_LENGTH;

         ium_mesg->f_list[currRead].delta = (int32 *)safe_malloc(sizeof(int32) * ium_mesg->f_list[currRead].delta_length);
         int i = 0;
         while (i < ium_mesg->f_list[currRead].delta_length) {
            ium_mesg->f_list[currRead].delta[i] = (int32) strtol(dlnStr,&dlnStr,10);
            i++;
         }
      }

      // read blank line
      fgets(line, AS_SEQAN_MAX_RESULT_LENGTH, tempOut);
   }
   fclose(tempOut);
}

void updateICMRecord(IntConConMesg *icm_mesg, char * inFile, char *seqAn, char *prefix, char * wrkDir) {
   // update the consensus
   icm_mesg->consensus = getConsensus(inFile, seqAn, prefix, wrkDir);
   icm_mesg->length = strlen(icm_mesg->consensus);

   // update quality
   icm_mesg->quality = (char *) safe_malloc(sizeof(char) * icm_mesg->length+1);
   memset(icm_mesg->quality, '1', icm_mesg->length);
   icm_mesg->quality[icm_mesg->length] = '\0';

   //update read data
   int currRead = 0;
   char line[AS_SEQAN_MAX_RESULT_LENGTH];
   char resultFile[AS_SEQAN_MAX_BUFFER_LENGTH];
   getFileName(prefix, wrkDir, AS_SEQAN_RESULT, resultFile);

   FILE *tempOut;

   tempOut = fopen(resultFile,"r");
   // skip header of output
   while (!feof(tempOut)) {
      fgets(line, AS_SEQAN_MAX_RESULT_LENGTH, tempOut);
      if (strncasecmp(line, "typ", 3) == 0) {
         break;
      }
   }

   // now read alignments of each read   
   for (currRead = 0; currRead < icm_mesg->num_pieces; currRead++) {
      if (currRead > 0) {
         // read the typ: line
         fgets(line, AS_SEQAN_MAX_RESULT_LENGTH, tempOut);
      }
      // read the seq: line         
      fgets(line, AS_SEQAN_MAX_RESULT_LENGTH, tempOut);
      // read the Pos line
      fgets(line, AS_SEQAN_MAX_RESULT_LENGTH, tempOut);
      chomp(line);
      CDS_COORD_t begin, end;
      sscanf(line,"Pos:"F_COORD","F_COORD,&begin,&end);
      icm_mesg->pieces[currRead].position.bgn = begin;
      icm_mesg->pieces[currRead].position.end = end;

      // read the dln line
      fgets(line, AS_SEQAN_MAX_RESULT_LENGTH, tempOut);
      chomp(line);
      sscanf(line,"dln:"F_S32, &icm_mesg->pieces[currRead].delta_length);

      // read the del line
      fgets(line, AS_SEQAN_MAX_RESULT_LENGTH, tempOut);
      chomp(line);

      if (icm_mesg->pieces[currRead].delta_length > 0) {
         char *dlnStr = line+AS_SEQAN_MAX_HEADER_LENGTH;

         icm_mesg->pieces[currRead].delta = (int32 *)safe_malloc(sizeof(int32) * icm_mesg->pieces[currRead].delta_length);
         int i = 0;
         while (i < icm_mesg->pieces[currRead].delta_length) {
            icm_mesg->pieces[currRead].delta[i] = (int32) strtol(dlnStr,&dlnStr,10);
            i++;
         }
      }
      // read blank line
      fgets(line, AS_SEQAN_MAX_RESULT_LENGTH, tempOut);
   }

   // now read the alignments of each unitig
   for (currRead = 0; currRead < icm_mesg->num_unitigs; currRead++) {
      // read the seq: line
      while (strncmp(line, "Pos:", 4) != 0) {
         fgets(line, AS_SEQAN_MAX_RESULT_LENGTH, tempOut);
      }
      // we read the Pos line above, process it now
      chomp(line);
      CDS_COORD_t begin, end;
      sscanf(line,"Pos:"F_COORD","F_COORD,&begin,&end);

      // read the dln line
      fgets(line, AS_SEQAN_MAX_RESULT_LENGTH, tempOut);
      chomp(line);
      sscanf(line,"dln:"F_S32, &icm_mesg->unitigs[currRead].delta_length);

      // read the del line
      char * del = readMultiLine(tempOut);

      if (currRead >= 0) {
         icm_mesg->unitigs[currRead].position.bgn = begin;
         icm_mesg->unitigs[currRead].position.end = end;   
         if (icm_mesg->unitigs[currRead].delta_length > 0) {
            char *dlnStr = del+AS_SEQAN_MAX_HEADER_LENGTH;
   
            icm_mesg->unitigs[currRead].delta = (int32 *)safe_malloc(sizeof(int32) * icm_mesg->unitigs[currRead].delta_length);
            int i = 0;
            while (i < icm_mesg->unitigs[currRead].delta_length) {
               icm_mesg->unitigs[currRead].delta[i] = (int32) strtol(dlnStr,&dlnStr,10);
               i++;
            }
         }
      }
      safe_free(del);
   }

   fclose(tempOut);
}

int main(int argc, char **argv) {
   int arg = 1;
   int err = 0;
   int hlp = 0;

   char * gkpStoreName  = NULL;
   int    gkpStorePart  = 0;
   char * msgFile       = NULL;
   char * outputFileName= NULL;
   char * seqAn         = NULL;
   char * wrkDir        = NULL;
   char * seqStoreName  = NULL;
   int    seqStoreVer   = 0;
   int    seqStorePart  = 0;     

   while (arg < argc) {
      if (strcmp(argv[arg], "-c") == 0) {
         msgFile = argv[++arg];
      } else if (strcmp(argv[arg], "-G") == 0) {
         gkpStoreName = argv[++arg];
      } else if (strcmp(argv[arg], "-S") == 0) {
         gkpStorePart = atoi(argv[++arg]);
      } else if (strcmp(argv[arg], "-o") == 0) {
         outputFileName = argv[++arg];
      } else if (strcmp(argv[arg], "-s") == 0) {
         seqAn = argv[++arg];
      } else if (strcmp(argv[arg], "-w") == 0) {
         wrkDir = argv[++arg];
      } else if (strcmp(argv[arg], "-u") == 0) {
         seqStoreName = argv[++arg];
      } else if (strcmp(argv[arg], "-V") == 0) {
         seqStoreVer = atoi(argv[++arg]);
      } else if (strcmp(argv[arg], "-p") == 0) {
         seqStorePart = atoi(argv[++arg]);
      } else {
         err++;
      }
      arg++;
   }

   if ((err) || (gkpStoreName == NULL) || (msgFile == NULL) || (outputFileName == NULL) || seqAn == NULL) {
      fprintf(stderr, "USAGE: SeqAn_CNS -G <gkpStore> -c <input.cgb> -o <output.cgi> -s <seqan_executable> [-u seqstore, required for contig consensus] [-w working directory]\n");      
      exit(1);
   }

   GateKeeperStore          *gkpStore = openGateKeeperStore(gkpStoreName, FALSE);
   loadGateKeeperPartition(gkpStore, gkpStorePart);
   
   fragRecord                fr;
   GenericMesg              *pmesg;
   tSequenceDB              *sequenceDB = NULL;   
   
   FILE *infp = fopen(msgFile,"r");
   FILE *tempReads;
   FILE *outfp = fopen(outputFileName, "w");
   char fileName[AS_SEQAN_MAX_BUFFER_LENGTH];
   char *prefix = outputFileName;
   getFileName(prefix, wrkDir, AS_SEQAN_INPUT_NAME, fileName);

   int i = 0;
   
   while ((EOF != ReadProtoMesg_AS(infp, &pmesg))) {
      int freeMem = 0;
     
      if (pmesg->t == MESG_IUM) {
         IntUnitigMesg *ium_mesg = (IntUnitigMesg *)pmesg->m;         
         
         if (strlen(ium_mesg->consensus) == 0) {
            tempReads = fopen(fileName,"w");

            for (i =0; i < ium_mesg->num_frags; i++) {
               // get the fragment sequence
               getFrag(gkpStore, ium_mesg->f_list[i].ident, &fr, FRAG_S_QLT);
               unsigned int   clrBeg = getFragRecordClearRegionBegin(&fr, AS_READ_CLEAR_OBT);
               unsigned int   clrEnd = getFragRecordClearRegionEnd  (&fr, AS_READ_CLEAR_OBT);
               char          *seqStart = getFragRecordSequence(&fr);
               char          *seq      = seqStart+clrBeg;

               seq[clrEnd] = 0;
               AS_UTL_writeFastA(tempReads,
                  seq, clrEnd-clrBeg,
                   ">"F_IID","F_IID"\n", ium_mesg->f_list[i].position.bgn, ium_mesg->f_list[i].position.end);
            }
            fclose(tempReads);
            updateRecord(ium_mesg, fileName, seqAn, prefix, wrkDir);
            freeMem = 1;
         }
         WriteProtoMesg_AS(outfp, pmesg);

         if (freeMem) {
            safe_free(ium_mesg->consensus);
            safe_free(ium_mesg->quality);
         }
      }
      else if (pmesg->t == MESG_ICM) {         
         IntConConMesg *icm_mesg = (IntConConMesg *)pmesg->m;

         if (seqStoreName == NULL) {
            fprintf(stderr, "USAGE: The -u option is required for contig consensus\n");
            exit(1);
         }
         if (sequenceDB == NULL) {
            sequenceDB = openSequenceDB(seqStoreName, FALSE, seqStoreVer);
            openSequenceDBPartition(sequenceDB, seqStorePart);
         }
         
         if (strlen(icm_mesg->consensus) == 0) {
            tempReads = fopen(fileName,"w");

            for (i =0; i < icm_mesg->num_pieces; i++) {
               // get the fragment sequence
               getFrag(gkpStore, icm_mesg->pieces[i].ident, &fr, FRAG_S_QLT);
               unsigned int   clrBeg   = getFragRecordClearRegionBegin(&fr, AS_READ_CLEAR_LATEST);
               unsigned int   clrEnd   = getFragRecordClearRegionEnd  (&fr, AS_READ_CLEAR_LATEST);
               char          *seqStart = getFragRecordSequence(&fr);
               char          *seq      = seqStart+clrBeg;

               seq[clrEnd] = 0;
               AS_UTL_writeFastA(tempReads,
                  seq, clrEnd-clrBeg,
                   ">"F_IID","F_IID"\n", icm_mesg->pieces[i].position.bgn, icm_mesg->pieces[i].position.end);
            }
            
            // now handle the unitig messages
            for (i =0; i < icm_mesg->num_unitigs; i++) {
               VA_TYPE(char) *ungappedSequence = CreateVA_char(0);
               VA_TYPE(char) *ungappedQuality  = CreateVA_char(0);
               MultiAlignT *uma = loadMultiAlignTFromSequenceDB(sequenceDB, icm_mesg->unitigs[i].ident, 1);
               assert(uma != NULL);
               
               GetMultiAlignUngappedConsensus(uma, ungappedSequence, ungappedQuality);
               char * seq = Getchar(ungappedSequence,0);

               AS_UTL_writeFastA(tempReads,
                  seq, strlen(seq),
                   ">"F_IID","F_IID"\n", icm_mesg->unitigs[i].position.bgn, icm_mesg->unitigs[i].position.end);
            }
            fclose(tempReads);

            updateICMRecord(icm_mesg, fileName, seqAn, prefix, wrkDir);
            freeMem = 1;
         }
         WriteProtoMesg_AS(outfp, pmesg);

         if (freeMem) {
            safe_free(icm_mesg->consensus);
            safe_free(icm_mesg->quality);
         }
      }
   }
      
   fclose(infp);
   fclose(outfp);

   return 0;
}
