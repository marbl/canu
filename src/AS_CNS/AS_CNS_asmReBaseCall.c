
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

const char *mainid = "$Id: AS_CNS_asmReBaseCall.c,v 1.31 2008-12-31 02:56:29 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <math.h>
#include <time.h>

// Celera Assembler includes:
#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_genericStore.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_Hash.h"
#include "AS_MSG_pmesg.h"
#include "AS_GKP_include.h"
#include "MultiAlignStore_CNS.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"

static const char *rcsid = "$Id: AS_CNS_asmReBaseCall.c,v 1.31 2008-12-31 02:56:29 brianwalenz Exp $";

static HashTable_AS *utgUID2IID;


#define DEBUG 0





static IntUnitigMesg* convert_UTG_to_IUM(SnapUnitigMesg* utgMesg)
     /*
       converts a  SnapUnitigMessage to an IntUnitigMesg; mostly a matter of UID->IID mapping
     */
{
  int i;
  AS_IID    iid;
  int32  iidStatus;
  AS_IID    *di;
  IntUnitigMesg *iumMesg =
    (IntUnitigMesg*) safe_malloc(sizeof(IntUnitigMesg));

#if DEBUG > 1
  fprintf(stderr,"UTG external acc %s\n", AS_UID_toString(utgMesg->eaccession));
#endif

  { // hash the uid to the iid
    if(ExistsInHashTable_AS(utgUID2IID,AS_UID_toInteger(utgMesg->eaccession), 0)){
      fprintf(stderr,"Encountered UTG UID %s more than once! DIE!\n",
	      AS_UID_toString(utgMesg->eaccession));
      exit(1);
    }
    InsertInHashTable_AS(utgUID2IID, AS_UID_toInteger(utgMesg->eaccession), 0, utgMesg->iaccession, 0);
  }

  /* Set all toplevel fields */

  iumMesg->iaccession      = utgMesg->iaccession;
  iumMesg->coverage_stat  = utgMesg->coverage_stat;
  iumMesg->microhet_prob  = utgMesg->microhet_prob;
  iumMesg->status         = utgMesg->status;
  iumMesg->length         = utgMesg->length;
  iumMesg->consensus      = strdup(utgMesg->consensus);
  iumMesg->quality        = strdup(utgMesg->quality);
  iumMesg->forced         = utgMesg->forced;
  iumMesg->num_frags      = utgMesg->num_frags;

  if( iumMesg->num_frags > 0 ){
    iumMesg->f_list = (IntMultiPos*) safe_malloc(iumMesg->num_frags*sizeof(IntMultiPos));

    for(i=0; i<iumMesg->num_frags; i++){
      iumMesg->f_list[i].type = utgMesg->f_list[i].type;
      iumMesg->f_list[i].sourceInt = 0;

      iid = getGatekeeperUIDtoIID(gkpStore, utgMesg->f_list[i].eident, NULL);

      if( iid == 0 ){
	fprintf(stderr,"Error: Unknown uid fragment ID %s at %s:%d\n",
		AS_UID_toString(utgMesg->f_list[i].eident),__FILE__,__LINE__);
	exit(1);
      }

      iumMesg->f_list[i].ident       = iid;
      iumMesg->f_list[i].delta_length = utgMesg->f_list[i].delta_length;
      iumMesg->f_list[i].position     = utgMesg->f_list[i].position;

      if( iumMesg->f_list[i].delta_length > 0 ){
	iumMesg->f_list[i].delta = utgMesg->f_list[i].delta;
      }
      else
	iumMesg->f_list[i].delta = NULL;
    }
  }

  return iumMesg;
}



static IntConConMesg* convert_CCO_to_ICM(SnapConConMesg* ccoMesg)
     /* converts a SnapConConMesg to an IntConConMessage.
	What happens ?
	- mostly UID to IID mapping
     */
{
  int i;
  AS_IID    iid;

  IntConConMesg *icmMesg = (IntConConMesg*) safe_malloc(sizeof(IntConConMesg));

  /* we assume that the numbers are in ascending order */

  icmMesg->iaccession = ccoMesg->iaccession;
  icmMesg->placed     = ccoMesg->placed;
  icmMesg->length     = ccoMesg->length;
  icmMesg->consensus  = strdup(ccoMesg->consensus);
  icmMesg->quality    = strdup(ccoMesg->quality);
  icmMesg->forced     = ccoMesg->forced;
  icmMesg->num_pieces = ccoMesg->num_pieces;
  icmMesg->num_unitigs= ccoMesg->num_unitigs;
  icmMesg->num_vars   = ccoMesg->num_vars;  // affects .asm/CCO

  if (icmMesg->num_vars > 0) {
     icmMesg->v_list = (IntMultiVar*) safe_malloc(ccoMesg->num_vars*sizeof(IntMultiVar));
     for(i=0; i<icmMesg->num_vars; i++) // i loop
     {
        icmMesg->v_list[i].position       = ccoMesg->vars[i].position;
        icmMesg->v_list[i].num_reads      = ccoMesg->vars[i].num_reads;
        icmMesg->v_list[i].min_anchor_size    = ccoMesg->vars[i].min_anchor_size;
        icmMesg->v_list[i].var_length     = ccoMesg->vars[i].var_length ;
        icmMesg->v_list[i].nr_conf_alleles=strdup(ccoMesg->vars[i].nr_conf_alleles);
        icmMesg->v_list[i].weights        = strdup(ccoMesg->vars[i].weights);
        icmMesg->v_list[i].var_seq        = strdup(ccoMesg->vars[i].var_seq);
        icmMesg->v_list[i].conf_read_iids = strdup(ccoMesg->vars[i].conf_read_iids);
     }

  } else {
    icmMesg->v_list=NULL;
  }

  if( icmMesg->num_pieces > 0 ){
    icmMesg->pieces = (IntMultiPos*) safe_malloc(icmMesg->num_pieces*sizeof(IntMultiPos));
    for(i=0; i<icmMesg->num_pieces; i++){// i loop
      icmMesg->pieces[i].type = ccoMesg->pieces[i].type;
      icmMesg->pieces[i].sourceInt = 0;

      iid = getGatekeeperUIDtoIID(gkpStore, ccoMesg->pieces[i].eident, NULL);
      if( iid == 0 ){
	fprintf(stderr,"Error: Unknown uid fragment ID %s at %s:%d\n",
		AS_UID_toString(ccoMesg->pieces[i].eident),__FILE__,__LINE__);
	exit(1);
      }
      icmMesg->pieces[i].ident       = iid;
      icmMesg->pieces[i].delta_length = ccoMesg->pieces[i].delta_length;
      icmMesg->pieces[i].position     = ccoMesg->pieces[i].position;

      if( icmMesg->pieces[i].delta_length > 0 ){
	icmMesg->pieces[i].delta = ccoMesg->pieces[i].delta; /*** COPY BY REFERENCE ***/
      }
      else
	icmMesg->pieces[i].delta = NULL;

    }
  }

  if( icmMesg->num_unitigs > 0 ){
    icmMesg->unitigs = (IntUnitigPos*) safe_malloc(icmMesg->num_unitigs*sizeof(IntUnitigPos));
    for(i=0; i<icmMesg->num_unitigs; i++){
      int32 *iid;
      icmMesg->unitigs[i].type  = ccoMesg->unitigs[i].type;
      if (!ExistsInHashTable_AS(utgUID2IID, AS_UID_toInteger(ccoMesg->unitigs[i].eident), 0)) {
	fprintf(stderr,"Error: Reference before definition for unitig UID %s at %s:%d\n",
		AS_UID_toString(ccoMesg->pieces[i].eident),__FILE__,__LINE__);
	exit(1);
      }
      icmMesg->unitigs[i].ident        = LookupValueInHashTable_AS(utgUID2IID, AS_UID_toInteger(ccoMesg->unitigs[i].eident), 0);
      icmMesg->unitigs[i].position     = ccoMesg->unitigs[i].position;
      icmMesg->unitigs[i].delta        = ccoMesg->unitigs[i].delta; /*** COPY BY REFERENCE ***/
      icmMesg->unitigs[i].delta_length = ccoMesg->unitigs[i].delta_length;
    }
  }
  return icmMesg;
}



static void
help_message(int argc, char *argv[])
{
    fprintf(stderr,"  Usage:\n\n"
    "  %s [options] -f FragStoreDir < [ASM file] > [new ASM file]\n"
    "\n Standard option flags:\n"
    "    -K           don't split alleles when calling consensus\n"
    "    -N           don't output variation record to .cns file\n"
    "    -w win_size  specify the size of the 'smoothing window' that will be used in consensus calling\n"
    "                 If two SNPs are located win_size or less bases apart one from another,\n"
    "                 then they will be treated as one block\n"
    "    -m           Load fragStorePartition into memory (default reads from disk)\n"
    "    -d int       Depth of Celera coverage below which to include external data in basecalling\n"
    "                    0 (default) indicates that external data should always be used\n"
    "                    1 yields the traditional behavior, which uses external only in absence of Celera\n"
    "                  > 1 will include publice data is the Celera depth falls below the given value\n"
    "    -X           Allow 'expert' options (following)\n"
    "\n Expert option flags:\n"
    "    -q string    Override default quality call parameters\n"
    "                    string is colon separated list of the form '%%f:%%d:%%f'\n"
    "                    where first field is estimated sequencing error rate (default: .015)\n"
    "                         second field is number of sequenced haplotypes (default: 1)\n"
    "                          third field is estimated SNP rate (default: 1/1000)\n"
    "\n Input: asm file on stdin\n"
    "\n Output: to stdout (and stderr)\n",
    argv[0]);
    exit(1);
}

int main (int argc, char *argv[]) {
    char *frgStorePath=NULL;
    char *gkpStorePath=NULL;


    /**************** Process Command Line Arguments *********************/
    /* Parse the argument list using "man 3 getopt". */
    int expert=0;
    int in_memory=0;
    CNS_Options options = { CNS_OPTIONS_SPLIT_ALLELES_DEFAULT,
                            CNS_OPTIONS_MIN_ANCHOR_DEFAULT };

    int ch,errflg=0,illegal_use=0,help_flag=0,iflags=0;

    optarg = NULL;

    argc = AS_configure(argc, argv);

    while ( !errflg &&
           ( (ch = getopt(argc, argv,
                 "f:g:hKM:mq:w:X")) != EOF))
    {
        switch(ch) {
	case 'f':
	  frgStorePath=optarg;
          iflags++;
          iflags++;
	  break;
	case 'g':
	  gkpStorePath=optarg;
          iflags++;
          iflags++;
	  break;
        case '?':
        case 'h':
          help_flag = 1;
          break;
        case 'K':
          options.split_alleles = 0;
          iflags++;
          break;
        case 'm':
          in_memory = 1;
          iflags++;
          break;
        case 'w':
          options.smooth_win = atoi(optarg);
          iflags++;
          iflags++;
          break;
        case 'X':
          expert = 1;
          iflags++;
          break;
        default :
          {
          help_flag = 1;
          fprintf(stderr,"Unrecognized option -%c",optopt);
          }
        }
    }
    if ( (argc - iflags) != 1)  help_flag = 1;

    if (help_flag)
        help_message(argc, argv);

    if ( illegal_use ) {
        fprintf(stderr,"\n %s -h provides usage information.\n",argv[0]);
        exit(1);
    }

    /****************          Open Fragment Store             ***********/

    gkpStore = openGateKeeperStore(argv[optind++], FALSE);

    if (in_memory)
      loadGateKeeperStorePartial(gkpStore, 0, 0, FRAG_S_QLT);

    /* initialize a unitig UID-to-IID hash table */
    utgUID2IID = CreateScalarHashTable_AS();

    /****************      Initialize reusable stores          ***********/
    sequenceStore = NULL;
    qualityStore = NULL;
    beadStore = NULL;
    columnStore = NULL;
    manodeStore = NULL;


    /**************** Prepare Unitig Store ****************************/
    unitigStore = CreateMultiAlignStoreT();

    /**************** Loop on Input Messages **************************/
    {
      int contig_count=0,unitig_count=0;
      VA_TYPE(char) *recalled_sequence=CreateVA_char(200000);
      VA_TYPE(char) *recalled_quality=CreateVA_char(200000);
      GenericMesg *pmesg;
      GenericMesg tmesg;
      MultiAlignT *ma;
      time_t t;
      t = time(0);
      fprintf(stderr,"# asmReBaseCall $Revision: 1.31 $ processing. Started %s\n",
	      ctime(&t));
      InitializeAlphTable();

      while ( (ReadProtoMesg_AS(stdin,&pmesg) != EOF)){

        switch(pmesg->t)
        {
          case MESG_UTG:
          {
	    SnapUnitigMesg *eunitig;
	    IntUnitigMesg *iunitig;
            eunitig = (SnapUnitigMesg *)(pmesg->m);
	    iunitig = convert_UTG_to_IUM(eunitig);
	    ma = CreateMultiAlignTFromIUM(iunitig, iunitig->iaccession, 0);
	    safe_free(iunitig->f_list);
	    safe_free(iunitig->consensus);
	    safe_free(iunitig->quality);
	    safe_free(iunitig);
	    SetMultiAlignInStore(unitigStore,ma->maID,ma);
	    WriteProtoMesg_AS(stdout,pmesg); // write out the unitig message
	    break;
          }
          case MESG_CCO:
          {
	    SnapConConMesg *econtig;
	    IntConConMesg *icontig;
            econtig = (SnapConConMesg *)(pmesg->m);
	    icontig = convert_CCO_to_ICM(econtig);
	    ma = CreateMultiAlignTFromICM(icontig, icontig->iaccession, 0);
	    MultiAlignContig_ReBasecall(ma,recalled_sequence,recalled_quality,&options);
	    {
	      char *ptr;

	      //ptr = Getchar(ma->consensus,0);
	      ptr = Getchar(recalled_sequence,0);
	      assert(strlen(ptr)==strlen(econtig->consensus));
	      strcpy(econtig->consensus,ptr);

	      //ptr = Getchar(ma->quality,0);
	      ptr = Getchar(recalled_quality,0);
	      assert(strlen(ptr)==strlen(econtig->quality));
	      strcpy(econtig->quality,ptr);

	    }
	    {
	      int i;
	      if(icontig->v_list!=NULL){
		for(i=0;i<icontig->num_vars;i++){
                  safe_free(icontig->v_list[i].nr_conf_alleles);
                  safe_free(icontig->v_list[i].weights);
		  safe_free(icontig->v_list[i].var_seq);
                  safe_free(icontig->v_list[i].conf_read_iids);
		}
		safe_free(icontig->v_list);
	      }
	    }
	    safe_free(icontig->consensus);
	    safe_free(icontig->quality);
	    safe_free(icontig->pieces);
	    safe_free(icontig->unitigs);
	    safe_free(icontig);
	    WriteProtoMesg_AS(stdout,pmesg); // write out the modified contig message
	    break;
	  }
          default:
            WriteProtoMesg_AS(stdout,pmesg);
        }
      }
    }

    DeleteMultiAlignStoreT(unitigStore);
    closeGateKeeperStore(gkpStore);
    return 0;
}
