
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
#include <stdio.h>
#include <assert.h>
#include <unistd.h> /* man 3 getopt */
#include <errno.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_UTL_Var.h"
VA_DEF(OFGMesg)
VA_DEF(OverlapMesg)
VA_DEF(IntMultiPos)
VA_DEF(IntUnitigMesg)

#include "AS_CGB_miniunitigger.h"

#define FILENAME_MAX_LEN 500
#define CMD_BUFFER_SIZE  500

static int accept_tmp_as_final_file
( const char *thePath1,
  const char *thePath2
  )
{
  char a_message[CMD_BUFFER_SIZE];
  int ierr;
  ierr = rename(thePath1, thePath2);
  if( ierr != 0 ) {
    if( errno == ENOENT ) {
      perror(thePath1);
      ierr = 0; // Accept the error ENOENT code as a NFS workaround.
    }
    sprintf(a_message,"Failed to rename %s to %s.\n",
            thePath1, thePath2);
    perror(a_message);
    // </usr/include/errno.h>
  }
  return ierr;
}

static void input_mesgs_to_VA
(
 FILE                 * fovl,
 VA_TYPE(OFGMesg)     * the_ofg_messages,
 VA_TYPE(OverlapMesg) * the_ovl_messages,
 VA_TYPE(char)        * the_ofg_source,
 VA_TYPE(char)        * the_ovl_source
)
{
  assert(NULL != fovl);
  assert(NULL != the_ofg_messages);
  assert(NULL != the_ovl_messages);
  assert(NULL != the_ofg_source);
  assert(NULL != the_ovl_source);
  {
    int nadt=0,nidt=0,nilk=0,nofg=0,novl=0,nirp=0,nibc=0,niba=0;
    GenericMesg *pmesg = NULL;
    MesgReader ReadMesg_AS = InputFileType_AS(fovl);
    while( EOF != ReadMesg_AS(fovl, &pmesg)) {
      const MessageType imesgtype = pmesg->t;
      
      //printf("input_mesgs_to_VA: pmesg->t = %d\n", imesgtype);
      switch(imesgtype) {
      case MESG_ADT: 
        nadt++;
        break;
      case MESG_IDT: 
        nidt++;
        break;
      case MESG_ILK: 
        nilk++;
        break;
      case MESG_OFG:
        nofg++;
        // We need to store the "source" field in the_ofg_source here.
        AppendVA_OFGMesg( the_ofg_messages, (OFGMesg *) pmesg->m );
        break;
      case MESG_OVL:
        novl++;
        AppendVA_OverlapMesg( the_ovl_messages, (OverlapMesg *) pmesg->m );
        break;
      case MESG_IRP:
        nirp++;
        break;
      case MESG_IBC:
        nibc++;
        break;
      case MESG_IBA:
        niba++;
        break;
      default:
        {
          fprintf(stderr,"Unexpected message type %d\n",imesgtype);
          assert(FALSE);
        }
        break;
      }
    }
    
    fprintf(stderr,
            "input_mesgs_to_VA: nadt = %d\n"
            "input_mesgs_to_VA: nidt = %d\n"
            "input_mesgs_to_VA: nilk = %d\n"
            "input_mesgs_to_VA: nofg = %d\n"
            "input_mesgs_to_VA: novl = %d\n"
            "input_mesgs_to_VA: nirp = %d\n"
            "input_mesgs_to_VA: nibc = %d\n"
            "input_mesgs_to_VA: niba = %d\n"
            , nadt, nidt, nilk, nofg, novl, nirp, nibc, niba);
    
  }
}


static void output_the_IUM_to_file
(/* Input Only*/
 MesgWriter                  WriteMesg_AS,
 const VA_TYPE(char)    *    the_imp_source,
 const VA_TYPE(char)    *    the_ium_source,
 VA_TYPE(IntMultiPos)   *    the_imp_messages,
 VA_TYPE(IntUnitigMesg) *    the_ium_messages,
 const int                   fragment_count_target, // zero means ignore
 const char * const          Graph_Store_File_Prefix
)
{
  const IntChunk_ID nchunks
    = (IntChunk_ID)GetNumVA_IntUnitigMesg(the_ium_messages);
  char          filename[FILENAME_MAX_LEN] = {0};
  FILE *        fcgb = NULL;
  int           fragment_count = 0;
  int           file_count = 0;
  IntChunk_ID   chunk_index = 0;

  for(chunk_index=0;chunk_index < nchunks; chunk_index++) {
    IntUnitigMesg * mychunk
      = GetVA_IntUnitigMesg(the_ium_messages,chunk_index);

    fragment_count
      += GetVA_IntUnitigMesg(the_ium_messages,chunk_index)->num_frags;

    if(fragment_count_target > 0) {
      if((fragment_count >= fragment_count_target)) {
        if(NULL != fcgb) { fclose(fcgb); fcgb = NULL;}
        fragment_count = 0;
      }
      if(NULL == fcgb) {
        file_count ++;
        sprintf(filename,"%s_%03d.cgb_tmp",Graph_Store_File_Prefix,file_count);
        fcgb = fopen(filename,"w");
        assert(NULL != fcgb);
      }
    } else {
      if(NULL == fcgb) {
        sprintf(filename,"%s.cgb_tmp",Graph_Store_File_Prefix);
        fcgb = fopen(filename,"a");
	// We are depending on the ADT messages having just been
	// written to a new file with a file name of filename.
        assert(NULL != fcgb);
      }
    }

    {
      GenericMesg   pmesg;
      pmesg.t = MESG_IUM;
      pmesg.m = mychunk;
      WriteMesg_AS(fcgb,&pmesg);
    }
  }

  if(fragment_count_target > 0) {
    if(NULL != fcgb) { fclose(fcgb); fcgb = NULL;}
    fragment_count = 0;
    if(NULL == fcgb) {
      file_count ++;
      sprintf(filename,"%s_%03d.cgb_tmp",Graph_Store_File_Prefix,file_count);
      fcgb = fopen(filename,"w");
      assert(NULL != fcgb);
    }
  }

  fclose(fcgb);

  {
    int iii = 0;
    int ierr;

    if(fragment_count_target > 0) {
      for( iii=0; iii <= file_count ; iii++) {
        char thePath1[CMD_BUFFER_SIZE-1]={0};
	sprintf(filename,"%s_%03d.cgb",Graph_Store_File_Prefix,iii);
        sprintf(thePath1,"%s_tmp",filename);
        ierr = accept_tmp_as_final_file( thePath1, filename);
        assert(ierr == 0);
        // The temporary batch info file could not be moved into its final
        // position.
      }
    } else {
      char thePath1[CMD_BUFFER_SIZE-1]={0};
      sprintf(filename,"%s.cgb", Graph_Store_File_Prefix);
      sprintf(thePath1,"%s_tmp",filename);
      ierr = accept_tmp_as_final_file( thePath1, filename);
      assert(ierr == 0);
      // The temporary batch info file could not be moved into its final
      // position.
    }
  }

}


int main(int argc, char * argv [])
{
  MiniUnitiggerObject * muo = createMiniUnitigger(0,0,0);
  int status = 0;

  VA_TYPE(OFGMesg)     * the_ofg_messages =  CreateVA_OFGMesg(0);
  VA_TYPE(OverlapMesg) * the_ovl_messages = CreateVA_OverlapMesg(0);
  VA_TYPE(char)        * the_ofg_source = CreateVA_char(0);
  VA_TYPE(char)        * the_ovl_source = CreateVA_char(0);
  VA_TYPE(IntMultiPos) * the_imp_messages = CreateVA_IntMultiPos(0);
  VA_TYPE(IntUnitigMesg) * the_ium_messages = CreateVA_IntUnitigMesg(0);
  VA_TYPE(char)        * the_imp_source = CreateVA_char(0);
  VA_TYPE(char)        * the_ium_source = CreateVA_char(0);

  int fragment_count_target = 0;
  
  assert(NULL != the_ofg_messages);
  assert(NULL != the_ovl_messages);
  assert(NULL != the_ofg_source);
  assert(NULL != the_ovl_source);
  assert(NULL != the_imp_messages);
  assert(NULL != the_ium_messages);
  assert(NULL != the_imp_source);
  assert(NULL != the_ium_source);
  
  for( ; optind < argc; optind++) {
    char * filename = argv[optind];
    fprintf(stderr,"INPUT FILE for VA: <%s>\n", filename);
    {
      FILE * fovl = fopen( filename,"r");
      input_mesgs_to_VA( fovl, the_ofg_messages, the_ovl_messages,
			 the_ofg_source, the_ovl_source);
      fclose(fovl);
    }
  }
  
  set_cgb_unique_cutoff_MiniUnitigger( muo, 5.0);
  set_overlap_error_threshold_MiniUnitigger( muo, 0.06);
  set_as_cgb_max_frag_iid_MiniUnitigger( muo, 100000);

  {
    RunMiniUnitiggerParams params =
    { the_ofg_messages,
      the_ovl_messages,
      the_ofg_source,
      the_ovl_source };
    RunMiniUnitiggerResults results =
    { the_imp_messages,
      the_ium_messages,
      the_imp_source,
      the_ium_source };

    run_MiniUnitigger( muo, &params, &results);
  }

  {

    MesgWriter WriteMesg_AS = OutputFileType_AS(AS_PROTO_OUTPUT);
    
    output_the_IUM_to_file
      (/* Input Only*/
       WriteMesg_AS,
       the_imp_source,
       the_ium_source,
       the_imp_messages,
       the_ium_messages,
       fragment_count_target,
       "deleteme"
       );
  }

  destroyMiniUnitigger( muo);
  exit(status);
}

