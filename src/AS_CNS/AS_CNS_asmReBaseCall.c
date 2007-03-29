#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

#include <math.h>
#include <time.h>

#ifdef X86_GCC_LINUX
#include <fpu_control.h>
#endif

// Celera Assembler includes:
#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_genericStore.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_UTL_ID_store.h"
#include "PrimitiveVA_MSG.h"
#include "AS_UTL_version.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_Hash.h"
#include "AS_UTL_PHash.h"
#include "AS_UTL_version.h"
#include "AS_MSG_pmesg.h"
#include "AS_GKP_include.h"
#include "MultiAlignStore_CNS.h"
#include "MultiAlignment_CNS.h"
#include "Globals_CNS.h"
#include "PublicAPI_CNS.h"

static const char CM_ID[] = "$Id: AS_CNS_asmReBaseCall.c,v 1.12 2007-03-29 08:33:46 brianwalenz Exp $";

static UIDHashTable_AS *utgUID2IID;


#define DEBUG 0

float CNS_SNP_RATE   = 0.0003; // Used to calculate BIAS
int   CNS_HAPLOTYPES = 1;   // Used to calculate BIAS
int   CNS_USE_PUBLIC = 0;   // Used to direct basecalling to include public data
int   CNS_CALL_PUBLIC = 0;   // Used to direct basecalling to favor public data


/***********************/
/* conversion routines */
/***********************/


static int fraguid2iid(uint64 uid){
  PHashValue_AS value;
  if(HASH_FAILURE == getGatekeeperUIDtoIID(gkpStore, uid, &value)){
    fprintf(stderr,"Tried to look up iid of unknown uid: " F_UID " -- DIE!\n",uid);
    exit (-1);
  }
  return (value.IID);
}



static IntUnitigMesg* convert_UTG_to_IUM(SnapUnitigMesg* utgMesg)
     /*
       converts a  SnapUnitigMessage to an IntUnitigMesg; mostly a matter of UID->IID mapping
     */
{
  int i;
  CDS_IID_t iid;
  int32  iidStatus;
  CDS_IID_t *di;
  IntUnitigMesg *iumMesg = 
    (IntUnitigMesg*) safe_malloc(sizeof(IntUnitigMesg));

#if DEBUG > 1
  fprintf(stderr,"UTG external acc " F_UID "\n",utgMesg->eaccession);
#endif

  { // hash the uid to the iid
    if(LookupInUID2IIDHashTable_AS(utgUID2IID,utgMesg->eaccession)!=NULL){
      fprintf(stderr,"Encountered UTG UID " F_UID " more than once! DIE!\n",
	      utgMesg->eaccession);
      exit(-1);
    }
    InsertInUID2IIDHashTable_AS(utgUID2IID,utgMesg->eaccession,utgMesg->iaccession);
  }

  /* Set all toplevel fields */
  
  iumMesg->iaccession      = utgMesg->iaccession;

#ifdef AS_ENABLE_SOURCE
  iumMesg->source         = strdup(utgMesg->source);
#endif
  iumMesg->coverage_stat  = utgMesg->coverage_stat;
  iumMesg->status         = utgMesg->status;
  iumMesg->a_branch_point = utgMesg->a_branch_point;
  iumMesg->b_branch_point = utgMesg->b_branch_point;
  iumMesg->length         = utgMesg->length;
  iumMesg->consensus      = strdup(utgMesg->consensus);
  iumMesg->quality        = strdup(utgMesg->quality);
  iumMesg->forced         = utgMesg->forced;
  iumMesg->num_frags      = utgMesg->num_frags;
  iumMesg->num_vars       = utgMesg->num_vars; 

  assert(iumMesg->num_vars==0);

  if( iumMesg->num_frags > 0 ){
    iumMesg->f_list = (IntMultiPos*) safe_malloc(iumMesg->num_frags*sizeof(IntMultiPos));

    for(i=0; i<iumMesg->num_frags; i++){
      iumMesg->f_list[i].type = utgMesg->f_list[i].type;
#ifdef AS_ENABLE_SOURCE
      iumMesg->f_list[i].sourceInt = atoi(utgMesg->f_list[i].source);
#endif
      iid = fraguid2iid(utgMesg->f_list[i].eident);
      
      if( iid == 0 ){
	char dummy[40];
	fprintf(stderr,"Error: Unknown uid fragment ID " F_UID " at %s:%d\n",
		utgMesg->f_list[i].eident,__FILE__,__LINE__);
	exit(-1);
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
  CDS_IID_t iid;

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
        icmMesg->v_list[i].anchor_size    = ccoMesg->vars[i].anchor_size;
        icmMesg->v_list[i].var_length     = ccoMesg->vars[i].var_length ;
        icmMesg->v_list[i].nr_conf_alleles=strdup(ccoMesg->vars[i].nr_conf_alleles);
        icmMesg->v_list[i].weights        = strdup(ccoMesg->vars[i].weights);
        icmMesg->v_list[i].var_seq        = strdup(ccoMesg->vars[i].var_seq);  
     }
      
  } else {
    icmMesg->v_list=NULL;
  }

  if( icmMesg->num_pieces > 0 ){ 
    icmMesg->pieces = (IntMultiPos*) safe_malloc(icmMesg->num_pieces*sizeof(IntMultiPos));
    for(i=0; i<icmMesg->num_pieces; i++){// i loop
      icmMesg->pieces[i].type = ccoMesg->pieces[i].type;
#ifdef AS_ENABLE_SOURCE
      icmMesg->pieces[i].sourceInt = atoi(ccoMesg->pieces[i].source);
#endif

      iid = fraguid2iid(ccoMesg->pieces[i].eident);
      if( iid == 0 ){
	char dummy[40];
	fprintf(stderr,"Error: Unknown uid fragment ID " F_UID " at %s:%d\n",
		ccoMesg->pieces[i].eident,__FILE__,__LINE__);
	exit(-1);
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
      iid = LookupInUID2IIDHashTable_AS(utgUID2IID,ccoMesg->unitigs[i].eident);
      if(iid==NULL){
	char dummy[40];
	fprintf(stderr,"Error: Reference before definition for unitig UID " F_UID " at %s:%d\n",
		ccoMesg->pieces[i].eident,__FILE__,__LINE__);
	exit(-1);
      }
      icmMesg->unitigs[i].ident = *iid;
      icmMesg->unitigs[i].position = ccoMesg->unitigs[i].position;
      icmMesg->unitigs[i].delta = ccoMesg->unitigs[i].delta; /*** COPY BY REFERENCE ***/
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
                            CNS_OPTIONS_SMOOTH_WIN_DEFAULT,
                            CNS_OPTIONS_MAX_NUM_ALLELES };

#ifdef X86_GCC_LINUX
   /*
  ** Set the x86 FPU control word to force double
  ** precision rounding rather than `extended'
  ** precision rounding. This causes base
  ** calls and quality values on x86 GCC-Linux
  ** (tested on RedHat Linux) machines to be
  ** identical to those on IEEE conforming UNIX
  ** machines.
  */
  fpu_control_t fpu_cw;

  fpu_cw = ( _FPU_DEFAULT & ~_FPU_EXTENDED ) | _FPU_DOUBLE;

  _FPU_SETCW( fpu_cw );
#endif

    int ch,errflg=0,illegal_use=0,help_flag=0,iflags=0;

    fprintf(stderr,"Version: %s\n",CM_ID);

    optarg = NULL;
    ALIGNMENT_CONTEXT=AS_CONSENSUS;
    
    while ( !errflg && 
           ( (ch = getopt(argc, argv, 
                 "d:f:g:G?hKM:mq:w:X")) != EOF))
    {
        switch(ch) {
        case 'd':
          {
            CNS_USE_PUBLIC = atoi(optarg);
          }
          iflags++;
          iflags++;
          break;
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
        case 'G':
          if ( ! expert ) {
             fprintf(stderr,"Command line switch %c requires -X; try adding -X...\n",
                  ch); 
             illegal_use = 1;
          } else {
            CNS_CALL_PUBLIC = 1;
          }
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
        case 'M':
          options.max_num_alleles = atoi(optarg);
          iflags++;
          iflags++;
          break;
        case 'm':
          in_memory = 1;
          iflags++;
          break;
        case 'q':
          if ( ! expert ) {
             fprintf(stderr,"Command line switch %c requires -X; try adding -X...\n",
                  ch); 
             illegal_use = 1;
          } else {
            sscanf(optarg,"%f:%d:%f",&CNS_SEQUENCING_ERROR_EST,&CNS_HAPLOTYPES,
                &CNS_SNP_RATE);
            if (!(CNS_SEQUENCING_ERROR_EST > 0) || 
                CNS_SEQUENCING_ERROR_EST > .10 ) 
            {
              fprintf(stderr,"ERROR: Sequencing error estimate (-q flag) should be "
                  "within (0,.10) (%4f was specified\n",
                  CNS_SEQUENCING_ERROR_EST);
              illegal_use = 1;
            }
            if (CNS_HAPLOTYPES < 1) {
              fprintf(stderr,"ERROR: Haplotypes sampled (-h flag) must be > 0 "
                             "(%d was specified\n",CNS_HAPLOTYPES);
              illegal_use = 1;
            }
            if ((CNS_SNP_RATE < 0) || CNS_SNP_RATE > .10 ) {
              fprintf(stderr,
                  "ERROR: SNP rate estimate (-s flag) should be within [0,.10) "
                  "(%4f was specified\n",CNS_SNP_RATE);
              illegal_use = 1;
            }
          }
          iflags++;
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
    utgUID2IID = CreateUIDHashTable_AS(5000000);

    /****************      Initialize reusable stores          ***********/
    sequenceStore = NULL;
    qualityStore = NULL;
    beadStore = NULL;
    columnStore = NULL;
    manodeStore = NULL;

    cnslog = stderr;
    cnsout = stderr;


    /**************** Prepare Unitig Store ****************************/
    unitigStore = CreateMultiAlignStoreT(0);

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
      fprintf(stderr,"# asmReBaseCall $Revision: 1.12 $ processing. Started %s\n",
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
	    ma = CreateMultiAlignTFromIUM(iunitig,-1,0);
#ifdef AS_ENABLE_SOURCE
	    safe_free(iunitig->source);
#endif
	    safe_free(iunitig->f_list);
	    safe_free(iunitig->consensus);
	    safe_free(iunitig->quality);
	    safe_free(iunitig);
	    SetMultiAlignInStore(unitigStore,ma->id,ma);
	    WriteProtoMesg_AS(stdout,pmesg); // write out the unitig message
	    break;
          }
          case MESG_CCO:
          {
	    SnapConConMesg *econtig;
	    IntConConMesg *icontig;
            econtig = (SnapConConMesg *)(pmesg->m);
	    icontig = convert_CCO_to_ICM(econtig);
	    ma = CreateMultiAlignTFromICM(icontig,-1,0);
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
          case MESG_ADT:
          {
	    AuditMesg *adt_mesg;
	    adt_mesg = (AuditMesg *)(pmesg->m);
	    VersionStampADT(adt_mesg,argc,argv);
	    WriteProtoMesg_AS(stdout,pmesg);
          }
          break;
          default:
            WriteProtoMesg_AS(stdout,pmesg);
        }
        fflush(cnsout);
        fflush(cnslog);
      }
    }

    DeleteMultiAlignStoreT(unitigStore);
    closeGateKeeperStore(gkpStore);
    return 0;
}
