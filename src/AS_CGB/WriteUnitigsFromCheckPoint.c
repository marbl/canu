
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
#include "AS_CGB_all.h"
#include "AS_CGB_unitigger_globals.h"

/****************************************************************************/
#define ANALYSIS_EXTENSION     ".cga"
#define TIMINGS 1

static void save_the_chunk_graph
(
 MesgWriter WriteMesg_AS,
 const char * const Graph_Store_File_Prefix,
 const char * const Graph_Store_File_Name,
 const int analysis_level,
 const char * bubble_boundaries_filename,
 const int output_fom_messages,
 const float cgb_unique_cutoff,
 const int iterator,
 const int fragment_count_target,
 int argc,
 char * argv [],
 TStateGlobals * gstate,
 THeapGlobals  * heapva
 )
{
  FILE *fcga = NULL;
  FILE *fcam = NULL;
  FILE *fwrn = stderr;
  FILE *fp_unitig_statistics = NULL;

  char strtmp2[CMD_BUFFER_SIZE-1];
  int ierr;
  time_t tp1, tp2;

  if(analysis_level > 0) {
    if(NULL != Graph_Store_File_Prefix) {
      sprintf(strtmp2,"%s%s.%d",Graph_Store_File_Prefix,
	      ANALYSIS_EXTENSION,iterator);
      SAFE_FOPEN( fcga, strtmp2, "w");
      sprintf(strtmp2,"%s%s.%d",Graph_Store_File_Prefix,".cam",iterator);
      SAFE_FOPEN( fcam, strtmp2, "w");
      sprintf(strtmp2,"%s%s.%d",Graph_Store_File_Prefix,".cus",iterator);
      SAFE_FOPEN( fp_unitig_statistics, strtmp2, "w");
    }

    if(TIMINGS) {
      tp1 = time(NULL); fprintf(stderr,"Begin chunk graph analysis\n");
      system_date();
    }
    chunk_graph_analysis
      (
       /* Input Only */
       analysis_level,
       gstate->max_frag_iid,
       heapva->frags,  /* The internal representation of
                           the fragment reads. */
       heapva->edges,  /* The internal representation of the  overlaps. */
       heapva->frag_annotations,
       gstate->nbase_in_genome,
       cgb_unique_cutoff,
       gstate->global_fragment_arrival_rate,
       bubble_boundaries_filename,
       heapva->chunkfrags,
               
       /* Modify the chunk annotation, by setting the "isrc" index into 
          chunksrc array. */
       heapva->thechunks,
       heapva->chunksrcs, 
       // The character handle used to store the annotations.
       fcga, fcam, fp_unitig_statistics, fwrn );

    SAFE_FCLOSE(fcga);
    SAFE_FCLOSE(fcam);
    SAFE_FCLOSE(fp_unitig_statistics);
    
    if(TIMINGS) {
      tp2 = time(NULL); 
      fprintf(stderr,"%10" F_TIME_TP " sec: Finished chunk graph analysis\n",
              (tp2-tp1));
      system_date();
    }
  }

  /* unitigger output */ {  

    if(TIMINGS) { 
      tp1 = time(NULL); fprintf(stderr,"Begin chunk graph message I/O\n");
      system_date();
    }

    {
      FILE *fcgb = NULL;
      char strtmp2[CMD_BUFFER_SIZE-1];
      
      fprintf(stderr,"Openning output file to append the "
              "IUM+UOM messages.\n");
      strcpy(strtmp2,Graph_Store_File_Prefix);

      if(fragment_count_target == 0) {
          sprintf(strtmp2,"%s.cgb_tmp",Graph_Store_File_Prefix);
      } else {
        int file_count = 0;
          sprintf(strtmp2,"%s_%03d.cgb_tmp",
		  Graph_Store_File_Prefix,file_count);
      }

      SAFE_FOPEN( fcgb, strtmp2, "w");

      {
        // Output the snapshot ADT message.
        GenericMesg mesg;
        AuditMesg  adt_mesg;
        
        /*  Audit Message record */
        mesg.t = MESG_ADT;
        mesg.m = &adt_mesg;
        mesg.s = sizeof(AuditMesg);
        adt_mesg.list = NULL;
        VersionStampADT(&adt_mesg, argc, argv);
        WriteMesg_AS(fcgb,&mesg);
      }
      
#ifdef USE_IBA_FILE
      {
        FILE *fiba = NULL;
        fprintf(stderr,"Openning intermediate file to read the "
                "IBA+ADT(+ADL)+IDT messages.\n");
        sprintf(strtmp2,"%s/fgb.iba",Graph_Store_File_Name);
        if(NULL == (fiba = fopen(strtmp2,"r"))){
          fprintf(stderr,"*WARNING: Can not open fgb.iba file %s\n",strtmp2);
          // exit(1);
        }
        
        if(NULL != fiba) {
          // Output the ADT messages from all the incremental processing.
          MesgReader ReadMesg_AS = (MesgReader)InputFileType_AS(fiba);
          GenericMesg *pmesg = NULL;
          while(EOF != ReadMesg_AS(fiba,&pmesg)) {
            WriteMesg_AS(fcgb,pmesg);
          }
        }
        if(NULL != fiba) { ierr = fclose(fiba); assert(ierr == 0); fiba = NULL;}
      }
#endif        
      if(NULL != fcgb) { ierr = fclose(fcgb); assert(ierr == 0); fcgb = NULL;}
    }
    
    output_the_chunks
      (/*Input Only*/
       WriteMesg_AS,
       heapva->frags,
       heapva->edges,
       heapva->frag_annotations,
       heapva->chunkfrags,
       heapva->thechunks,
       heapva->chunkseqs,
       heapva->chunkquas,
       heapva->chunksrcs,
       analysis_level,
       output_fom_messages,
       gstate->global_fragment_arrival_rate,
       fragment_count_target,
       Graph_Store_File_Prefix
       // heapva->the_imp_messages,
       // heapva->the_ium_messages
       );


    if(TIMINGS) {
      tp2 = time(NULL); 
      fprintf(stderr,
              "%10" F_TIME_TP " sec: Finished chunk graph I/O\n",(tp2-tp1));
      system_date();
    }
  }
}

int main(int argc, char * argv[]) {
  
  TStateGlobals    * gstate = calloc(1,sizeof(TStateGlobals));
  THeapGlobals     * heapva = calloc(1,sizeof(THeapGlobals));
  const size_t new_additional_number_of_frags = 0;
  const size_t new_additional_number_of_edges = 0;
  const size_t new_additional_amount_of_text  = 0;
 
  UnitiggerGlobals * rg = calloc(1,sizeof(UnitiggerGlobals));

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && 
	   ((ch = getopt(argc, argv, 
			 "A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:PQ:R:S:V:U:W:XY:Z:"
			 "ab:cd:e:fg:hi:j:l:m:n:o:p:qr:st:u:v:w:x:y:z:"
                         "156:789"
                         )) != EOF)) {
      switch(ch) {
        /* The required command line options: */
      case 'A':
        // -A <int> : Specify an analysis level.
        rg->analysis_level = atoi(optarg);
        //analysis_flag = (analysis_level > 0);
        rg->analysis_flag = rg->analysis_level;
        break;
      case 'B':
        // -B <int> : Specifies a target number of fragments in IMP
        // records per cgb file.
        rg->fragment_count_target = atoi(optarg);
        break;
      case 'h':
        // help
        // goto UsageStatement;
        // Usage( argv[0], NULL );
        break;
      case 'i':
        // -i <file> : The input fragment overlap graph store.
#ifndef USE_STATIC_CMD_STRINGS
        rg->Input_Graph_Store = optarg;
#else // USE_STATIC_CMD_STRINGS
        strcpy(rg->Input_Graph_Store,optarg);
#endif // USE_STATIC_CMD_STRINGS
        fprintf(stderr,"* The input fragment graph store is <%s>.\n",
                rg->Input_Graph_Store);
	assert(NULL != rg->Input_Graph_Store);
        rg->input_store_flag = (rg->Input_Graph_Store[0] != '\0');
        break;
      case 'j':
        // -j <float> : Astat cut-off threshold
	rg->cgb_unique_cutoff = atof(optarg);
	fprintf(stderr,"* cgb_unique_cutoff set to %f\n", rg->cgb_unique_cutoff);
	break; 
      case 'o':
        // -o <file> : The output fragment overlap graph store.
#ifndef USE_STATIC_CMD_STRINGS
        rg->Output_Graph_Store = optarg;
#else // USE_STATIC_CMD_STRINGS
        strcpy(rg->Output_Graph_Store,optarg);
#endif // USE_STATIC_CMD_STRINGS
        fprintf(stderr,"* The output fragment graph store is <%s>.\n",
                rg->Output_Graph_Store);
	assert(NULL != rg->Output_Graph_Store);
        rg->output_store_flag = (rg->Output_Graph_Store[0] != '\0');
        break;
      case 'O':
        // -O <file> : The output file.
#ifndef USE_STATIC_CMD_STRINGS
        rg->Output_Graph_Store_Prefix = optarg;
#else // USE_STATIC_CMD_STRINGS
        strcpy(rg->Output_Graph_Store_Prefix,optarg);
#endif // USE_STATIC_CMD_STRINGS
        fprintf(stderr,"* The output fragment graph store is <%s>.\n",
                rg->Output_Graph_Store_Prefix);
	assert(NULL != rg->Output_Graph_Store_Prefix);
        rg->output_store_flag = (rg->Output_Graph_Store_Prefix[0] != '\0');
        break;
      case 'P':
        // -P : Any "protoIO" messages will be in ASCII rather than
        // binary.
	fprintf(stderr,"** ASCII mode\n");
        rg->as_proto_output = TRUE;
        break;
     default:
       fprintf(stderr,"Unsupported option <%c>\n", ch);
        break;
      }
    }
//  UsageStatement:
  }

  { 
    //    rg->analysis_level            = 0; // atoi(argv[3]);
    rg->bubble_boundaries_filename = NULL;
    rg->output_fom_messages       = FALSE;
    rg->num_cgb_passes            = 2;
  }
  
  fprintf(stderr,"Made it to here 2 \n");
  open_fgb_store( gstate, heapva);
  fprintf(stderr,"Before reading checkpoint.\n");
  read_fgb_store( rg->Input_Graph_Store, gstate, heapva,
                  new_additional_number_of_frags,
                  new_additional_number_of_edges,
                  new_additional_amount_of_text);

  fprintf(stderr,"After reading checkpoint.\n");
  ReportHeapUsage_CGB( gstate, heapva, stderr); 

  {
    MesgWriter WriteMesg_AS = NULL;
    MesgWriter ErrorWriter_AS = NULL;

    // VersionStamp(argc,argv);
    WriteMesg_AS = (MesgWriter)OutputFileType_AS((rg->as_proto_output
                                      ? AS_PROTO_OUTPUT : AS_BINARY_OUTPUT));
    ErrorWriter_AS = (MesgWriter)OutputFileType_AS(AS_PROTO_OUTPUT);
    
    if(NULL != rg->Output_Graph_Store_Prefix) {
      save_the_chunk_graph
        (
         WriteMesg_AS,
         rg->Output_Graph_Store_Prefix,
         rg->Output_Graph_Store,
         rg->analysis_level,
         rg->bubble_boundaries_filename,
         rg->output_fom_messages,
         rg->cgb_unique_cutoff,
         rg->num_cgb_passes,
         rg->fragment_count_target,
         argc,
         argv,
         gstate,
         heapva
         );
    }
  }

  close_fgb_store( gstate, heapva);
  return 0;
}
