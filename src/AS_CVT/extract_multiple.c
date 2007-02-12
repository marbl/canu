
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
static char CM_ID[] = "$Id: extract_multiple.c,v 1.7 2007-02-12 22:16:57 brianwalenz Exp $";

/*
  IMPORTANT NOTE:

  For reasons as yet unknown, this program core dumps on exit. Please
  set "ulimit -c 0" before running.
*/


/*********************************************************************/
// headers
/*********************************************************************/
// Standard header files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <math.h>
#include <unistd.h>

// Project header files
#include "AS_global.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_PHash.h"
#include "AS_UTL_Var.h"
#include "AS_ALN_aligners.h"
#include "AS_CVT_hashtable.h"
#include "AS_CVT_uid_interface.h"
#include "AS_GKP_include.h"

#define STRING_LENGTH 1024

#define MAX_FILE_DESCRIPTORS  4096

#define MAX_COMPONENTS (MAX_FILE_DESCRIPTORS/4)
//#define MAX_COMPONENTS 1

#define NUM_RECRUITERS     1000000
#define NUM_REC_FRAGS   1000000000

#define MODERATE_NUMBER      10000
#define BIG_NUMBER         1000000

#define DEBUG      // minimal debugging
#define DEBUG1     // more detailed debugging
//#define DEBUG2     // voluminous detailed debugging

#define HASH_FILE  "hash.lib"
#define FRAGS_FILE  "frags.lib"
#define MATES_FILE "mates.lib"


typedef struct
{
  CDS_UID_t uid;
  char       fn[STRING_LENGTH];
  FILE     * fp;
} ComponentFile;

// component structure
typedef struct
{
  char          input[STRING_LENGTH];
  char          directory[STRING_LENGTH];
  ComponentFile internal;
  ComponentFile external;
} Component;

VA_DEF(Component)


typedef struct
{
  CDS_IID_t iid1;
  CDS_IID_t iid2;
} IIDPair;


typedef struct
{
  CDS_IID_t mate;
  CDS_IID_t dist;
  OrientType orient;
} IIDOrient;

VA_DEF(IIDPair);


typedef struct
{
  CDS_UID_t uid;
  CDS_IID_t iid;
  CDS_IID_t start;
  cds_uint32 num;
} Recruiter;
  
// structure to hold store handles
typedef struct
{
  GateKeeperStore gkp_store;
  FragStoreHandle frg_store;
} StoreSet;
typedef StoreSet * StoreSetp;


// structure to hold iids/uids of fragments & mates
typedef struct
{
  IntFragment_ID  iid;
  Fragment_ID     uid;
  FragType        type;
  IntFragment_ID  mate_iid;
  Fragment_ID     mate_uid;
} IDSet;
typedef IDSet * IDSetp;


// structure to hold command-line & other 'global' items
typedef struct
{
  char           * program_name;        // name of this program
  char           * gkp_store_name;      // grande gkpstore name
  char           * frg_store_name;      // grande fragstore name
  char           * scaffold_filename;   // scaffold UID - frag UID pairs file
  char           * input_store;         // input store
  char           * output_store;        // output store
  char           * directory;           // directory for output
  char           * input_filename;      // name of input file
  UIDInteractorp   uid_fetcher;         // object that gets real UIDs
  StoreSetp        s_set;               // set of stores
  ReadStructp      rs;                  // reusable readStruct
  VA_TYPE(Component)    * components;   // component structures
  // the following are component definition dependent
  VA_TYPE(IIDPair)      * dists;        // master list of dists-components
  VA_TYPE(IIDPair)      * frags;        // master list of frags-components
  // The following may be persistent across runs
  IIDOrient             * mates_dists;  // list of frag-mate iids
  cds_uint32              num_mates;
  HashTable             * uid_ht;       // hashtable of Recruiter items
  cds_uint32              num_ht;
  VA_TYPE(CDS_IID_t)   * rec_frags;    // pointed to by recruiter hashtable
} Globals;
typedef Globals * Globalsp;

/*
  rc = CreateVA_cds_uint16(10000);
  ResetVA_cds_uint16(rc);
  EnableRangeVA_cds_uint16(rc, GetMultiAlignLength(ma));
  
  a = GetNumVA_cds_uint16(rc);
  b = *(GetVA_cds_uint16(rc,0));
  
  cds_uint16 * temp = GetVA_cds_uint16(rc, i);
  (*temp)++;
  
  AppendVA_cds_uint16(rc, &c);

  DeleteVA_cds_uint16(rc);
*/

/****************************************************************************
 * Simple file line counting & reading
 ****************************************************************************/
int GetNextLine( FILE * fp, char * string, int length )
{
  string[0] = '\0';
  while( fgets( string, length, fp ) )
  {
    if( string[0] != '#' )
      return 0;
  }
  return 1;
}


int CountLines( FILE * fp )
{
  int num_lines = 0;
  char string[STRING_LENGTH];

  while( !GetNextLine( fp, string, STRING_LENGTH ) )
    num_lines++;

  rewind( fp );
  return num_lines;
}


int ReadNextLine( FILE * fp, char * dest )
{
  char string[STRING_LENGTH];
  
  dest[0] = '\0';
  if( ! GetNextLine( fp, string, STRING_LENGTH ) )
  {
    if( string[strlen( string ) - 1] == '\n' )
      string[strlen( string ) - 1] = '\0';
    strcpy( dest, string );
    return 0;
  }
  return 1;
}


/****************************************************************************
 * Populate a FragMesg from a fragStore given an iid
 ****************************************************************************/
int PopulateFragment( FragMesg * f,
                      FragStoreHandle fs,
                      CDS_IID_t iid,
                      ReadStructp rs )
{
  cds_uint32 bgn, end;

  if( getFragStore( fs, iid, FRAG_S_ALL, rs ) )
  {
    fprintf( stderr, "Failed to get fragment " F_IID " data\n", iid );
    return 1;
  }

  f->action = AS_ADD;
  getAccID_ReadStruct( rs, &(f->eaccession) );
  getReadType_ReadStruct( rs, &(f->type) );
  getEntryTime_ReadStruct (rs, &(f->entry_time) );
  getClearRegion_ReadStruct( rs, &bgn, &end, READSTRUCT_LATEST );
  f->clear_rng.bgn = bgn;
  f->clear_rng.end = end;
  getSource_ReadStruct( rs, f->source, STRING_LENGTH );
  getSequence_ReadStruct( rs, f->sequence, f->quality, AS_READ_MAX_LEN );
  getReadIndex_ReadStruct( rs, &(f->iaccession) );
  
  return 0;
}

/****************************************************************************
 * Code for using several stores
 ****************************************************************************/
void CloseStores( StoreSetp s_set )
{
  if( s_set )
  {
    closeGateKeeperStore( &(s_set->gkp_store) );
    free( s_set );
  }
}


StoreSetp OpenStores( char * gkp_store_name,
                      char * frg_store_name )

{
  StoreSetp s_set;

  // allocate the containing structure
  if( (s_set = (StoreSetp) calloc( 1, sizeof( StoreSet ) )) == NULL )
  {
    fprintf( stderr, "Failed to allocate StoreSet\n" );
    return NULL;
  }

  // Open grande stores
  {
    // open the regional gatekeeper store
    initGateKeeperStore( &(s_set->gkp_store), gkp_store_name );
    openReadOnlyGateKeeperStore( &(s_set->gkp_store) );
    
    // open the regional fragment store
    s_set->frg_store = openFragStore( frg_store_name, "r" );
    if( s_set->frg_store == NULLSTOREHANDLE )
    {
      fprintf( stderr, "Failed to open grande frag store %s\n",
               frg_store_name );
      CloseStores( s_set );
      return NULL;
    }
  }
  
  return s_set;
}




int WriteBATADT( Globalsp globals,
                 ComponentFile * cf,
                 char * input_filename)
{
  GenericMesg  gen;
  BatchMesg    bat;
  AuditLine    adl;
  AuditMesg    adt;
  char         comment[STRING_LENGTH];

  // first message must be BAT with true UID;
  bat.name = globals->program_name;
  bat.created = time(0);
  bat.eaccession = cf->uid;
  bat.comment = "Grande dataset from Recruiter UIDs";
  gen.t = MESG_BAT;
  gen.m = &bat;
  WriteProtoMesg_AS( cf->fp, &gen );
  
  // next message should/must be ADT
  adl.next = NULL;
  adl.name = globals->program_name;
  adl.complete = time(0);
  adl.version = "Beta";

  sprintf(comment, "Based on component definition:\n");
  strcat(comment, input_filename);
  strcat(comment, "\n");

  adl.comment = comment;
  adt.list = &adl;
  gen.t = MESG_ADT;
  gen.m = &adt;
  WriteProtoMesg_AS( cf->fp, &gen );

  return 0;
}



Fragment_ID GetFragmentUID( GateKeeperStore gkp_store,
                            IntFragment_ID iid )
{
  GateKeeperFragmentRecord frag;

  getGateKeeperFragmentStore( gkp_store.frgStore, iid, &frag );

  return frag.readUID;
}


FragType GetFragmentType( GateKeeperStore gkp_store,
                          IntFragment_ID iid )
{
  GateKeeperFragmentRecord frag;

  getGateKeeperFragmentStore( gkp_store.frgStore, iid, &frag );

  return frag.type;
}


Distance_ID GetDistanceUID( GateKeeperStore gkp_store,
                            IntDistance_ID iid )
{
  GateKeeperDistanceRecord dist;

  getIndexStore( gkp_store.dstStore, iid, &dist );
  
  return dist.UID;
}


int PopulateDistanceMesg( GateKeeperStore gkp_store,
                          IDSetp idset,
                          DistanceMesg * dst )
{
  GateKeeperDistanceRecord dist;

  getIndexStore( gkp_store.dstStore, idset->iid, &dist );
  dst->action     = AS_ADD;
  dst->eaccession = idset->uid = dist.UID;
  dst->mean       = dist.mean;
  dst->stddev     = dist.stddev;

  return 0;
}


IntFragment_ID FindMate( GateKeeperStore gkp_store,
                         IntFragment_ID iid,
                         GateKeeperLinkRecord * link )
{
  GateKeeperFragmentRecord     frag;
  GateKeeperLinkRecordIterator iterator;

  // get the fragment
  getGateKeeperFragmentStore( gkp_store.frgStore, iid, &frag );

  if( frag.linkHead != 0 )
  {
    CreateGateKeeperLinkRecordIterator( gkp_store.lnkStore,
                                        frag.linkHead,
                                        iid,
                                        &iterator );

    while( NextGateKeeperLinkRecordIterator( &iterator, link ) )
    {
      if(link->type == AS_MATE)
        return ((link->frag1 == iid) ? link->frag2 : link->frag1);
    }
  }
  
  return 0;
}

static int LookupGateKeeperRecords(StoreSetp s_set,
                                   GateKeeperBatchRecord *gkpb,
                                   GateKeeperBatchRecord *gkpb_next,
                                   GateKeeperDistanceRecord *gkpd,
                                   GateKeeperSequenceRecord *gkps)
{
  if( getGateKeeperSequenceStore( (s_set->gkp_store).seqStore, gkpl->sequenceID, gkps) )
  {
     fprintf( stderr,
              "Failed to get sequence " F_IID " data from gatekeeperstore\n",
              gkpl->sequenceID );
     return 1;
  }
  if( getGateKeeperBatchStore( (s_set->gkp_store).batStore, gkpl->birthBatch, gkpb) )
  {
     fprintf( stderr, "Failed to get batch %u data from gatekeeperstore\n",
              gkpl->birthBatch );
     return 1;
  }
  if( getGateKeeperBatchStore( (s_set->gkp_store).batStore, gkpl->birthBatch+1, gkpb_next) )
  {
     fprintf( stderr, "Failed to get batch %u data from gatekeeperstore\n",
              gkpl->birthBatch+1 );
     return 1;
  }
  if( getGateKeeperDistanceStore( (s_set->gkp_store).dstStore, gkpl->lengthID, gkpd) )
  {
     fprintf( stderr,
              "Failed to get distance " F_IID " data from gatekeeperstore\n",
              gkpl->lengthID );
     return 1;
  }
  return 0;
}


// Usage message followed by abortive exit
void Usage( char * program_name, char * message )
{
  if( message != NULL )
    fprintf( stderr, "%s: %s\n", program_name, message );
  
  fprintf( stderr, "Usage: %s\n", program_name );
  fprintf( stderr, "\t[-g grande gatekeeper store]\n" );
  fprintf( stderr, "\t[-f grande fragment store]\n" );
  fprintf( stderr, "\t[-s scaffold UID - fragment UID file]\n" );
  fprintf( stderr, "\t[-i UID-UID pair store to use (instead of -s)]\n" );
  fprintf( stderr, "\t[-o UID-UID pair store to create]\n" );
  fprintf( stderr, "\t[-d output directory]\n" );
  fprintf( stderr, "\t[-c file of component definition filenames]\n" );

  fprintf(stderr,
          "\nNew directories will be created under output directory.\n"
          "using input component names\n");
  fprintf( stderr,
           "\nThese directories will contain uniquely named .frg files\n" );

  fprintf( stderr, "\nIMPORTANT NOTE:\n"
           "For reasons as yet unknown, this program core dumps on exit.\n"
           "Please set 'ulimit -c 0' before running.\n\n" );
  exit( 1 );
}


void ProcessCommandLine( Globals * globals, int argc, char ** argv )
{
  int ch, errflg = 0;
  
  // set default values
  memset( globals, 0, sizeof( Globals ) );
  globals->program_name = argv[0];

  // do the standard loop to get parameters
  optarg = NULL;
  while( !errflg && ((ch = getopt( argc, argv, "g:f:s:b:i:o:d:c:" )) != EOF) )
  {
    switch( ch )
    {
      case 'g':
        globals->gkp_store_name = optarg;
        fprintf( stdout, "* Grande gatekeeper store: %s\n",
                 globals->gkp_store_name );
        break;
      case 'f':
        globals->frg_store_name = optarg;
        fprintf( stdout, "* Grande fragment store: %s\n",
                 globals->frg_store_name );
        break;
      case 's':
        globals->scaffold_filename = optarg;
        fprintf( stdout, "* Scaffold UID - frag UID filename: %s\n",
                 globals->scaffold_filename );
        break;
      case 'i':
        globals->input_store = optarg;
        fprintf( stdout, "* Input recruiter UID - frag UID store: %s\n",
                 globals->input_store );
        break;
      case 'o':
        globals->output_store = optarg;
        fprintf( stdout, "* Output recruiter UID - frag UID store: %s\n",
                 globals->output_store );
        break;
      case 'd':
        globals->directory = optarg;
        fprintf( stdout, "* Output directory: %s\n",
                 globals->directory );
        break;
      case 'c':
        globals->input_filename = optarg;
        fprintf( stdout, "* Component definitions input filename: %s\n",
                 globals->input_filename );
        break;
      default:
        Usage( globals->program_name, "Unknown parameter or flag." );
        break;
    }
  }

  // check that all parameters were supplied
  if( ! globals->gkp_store_name )
    Usage( globals->program_name, "Specify a grande gatekeeper store." );
  if( ! globals->frg_store_name )
    Usage( globals->program_name, "Specify a grande fragment store." );
  if( ! globals->directory )
    Usage( globals->program_name, "Specify an output directory." );
  if( ! globals->input_filename )
    Usage( globals->program_name,
           "Specify a component definitions filename." );
  if( ! globals->input_store )
  {
    if( ! globals->scaffold_filename )
      Usage( globals->program_name,
             "Specify recruiter UID - frag UID filenames or an input store." );
  }
  else
  {
    if( globals->scaffold_filename )
      Usage( globals->program_name,
             "-s is incompatible with -i." );
  }
}


Globalsp CreateGlobals(void)
{
  Globalsp globals;
  globals = (Globalsp) calloc(1, sizeof(Globals));
  return globals;
}


void FreeGlobals(Globalsp globals)
{
  if(globals)
  {
    if(globals->uid_fetcher)
      DestroyUIDInteractor( globals->uid_fetcher );
    if(globals->rs)
      delete_ReadStruct( globals->rs );
    if(globals->s_set)
      CloseStores( globals->s_set );
    if(globals->components)
      DeleteVA_Component(globals->components);
    if(globals->dists)
      DeleteVA_IIDPair(globals->dists);
    if(globals->frags)
      DeleteVA_IIDPair(globals->frags);
    if(globals->mates_dists)
      free(globals->mates_dists);
    free(globals);
  }
}


int CreateComponentFile(Globalsp globals,
                        char * subdir,
                        ComponentFile * cf,
                        char * pref)
{
  time_t seconds;
  struct tm timer;
  
  seconds = time( 0 );
  if( localtime_r( &seconds, &timer ) == NULL )
  {
    fprintf( stderr, "Failed to get time!\n" );
    return 1;
  }
  
  if(GetUID(globals->uid_fetcher, &(cf->uid)))
  {
    fprintf(stderr, "Failed to get a UID for an output file\n");
    return 1;
  }
  
  // create an output filename
  sprintf( cf->fn, "%s/%s/%s_%04d%02d%02d_%020" F_UIDP "_HUMAN_FRAGS.frg",
           globals->directory,
           subdir,
           pref,
           timer.tm_year + 1900,
           timer.tm_mon + 1,
           timer.tm_mday,
           cf->uid );

  if((cf->fp = fopen(cf->fn, "w"))==NULL)
  {
    fprintf(stderr, "Failed to open file %s for writing.\n", cf->fn);
    return 1;
  }

  // write a BAT & ADT message
  if(WriteBATADT( globals, cf, subdir ))
  {
    fprintf(stderr, "Failed to write BAT & ADT message\n" );
    return 1;
  }
  return 0;
}


int compareIIDs(const CDS_IID_t * a, const CDS_IID_t * b)
{
  if(*a > *b) return 1;
  if(*a < *b) return -1;
  return 0;
}


int AppendToMaster(VA_TYPE(IIDPair) * master,
                   VA_TYPE(CDS_IID_t) * iids,
                   CDS_IID_t component_index)
{
  IIDPair iid_pair;
  CDS_IID_t iid;
  cds_uint32 i;

  if(GetNumVA_CDS_IID_t(iids) > 1)
  {
    qsort(iids->Elements,
          GetNumVA_CDS_IID_t(iids),
          sizeof(CDS_IID_t),
          (int (*) (const void *, const void *)) compareIIDs);
  }

  // append to master iid-component lists
  iid_pair.iid2 = component_index;
  for(iid = 0, i = 0; i < GetNumVA_CDS_IID_t(iids); i++)
  {
    if(iid != *GetVA_CDS_IID_t(iids, i))
    {
      iid_pair.iid1 = *GetVA_CDS_IID_t(iids, i);
      AppendVA_IIDPair(master, &iid_pair);
      iid = iid_pair.iid1;
    }
  }
  
  return 0;
}


int AddComponent(Globalsp globals, char * component_filename)
{
  cds_uint32 c_i;
  char string[STRING_LENGTH];
  static VA_TYPE(CDS_UID_t) * uids;
  static VA_TYPE(CDS_IID_t) * frag_iids;
  static VA_TYPE(CDS_IID_t) * dist_iids;
  FILE * fp;
  Component comp;
  char * temp_char;
  CDS_IID_t * iid;

#ifdef DEBUG
  fprintf(stderr, "Adding component from definition file %s\n",
          component_filename);
#endif
  
  if((temp_char = strrchr(component_filename, (int) '/')) == NULL)
    strcpy(comp.input, component_filename);
  else
    strcpy(comp.input, &(temp_char[1]));

#ifdef DEBUG1
  fprintf(stderr, "Assembly directory is %s/%s\n",
          globals->directory, comp.input);
#endif
  
  // create a component directory
  sprintf(string,"mkdir %s/%s", globals->directory, comp.input);
  if(system(string))
  {
    fprintf(stderr, "Failed to make directory via:\n%s\n", string);
    fprintf(stderr, "Assuming directory already exists....\n");
  }
  
  // create the components variable array if it doesn't exist
  if(globals->components == NULL)
  {
#ifdef DEBUG1
    fprintf(stderr, "Creating array of %d components\n", MAX_COMPONENTS);
#endif
    globals->components = CreateVA_Component(MAX_COMPONENTS);
    if(globals->components == NULL)
    {
      fprintf(stderr, "Failed to allocate array of %d components\n",
              MAX_COMPONENTS);
      return 1;
    }
    ResetVA_Component(globals->components);
    // EnableRangeVA_Component(globals->components, MAX_COMPONENTS);
  }

  // add the component to the variable array of components
  {
    if(CreateComponentFile(globals,
                           comp.input,
                           &(comp.internal),
                           "INT"))
    {
      fprintf(stderr,
              "Failed to create internal assembly input file for %s\n",
              comp.input);
      return 1;
    }
    if(CreateComponentFile(globals,
                           comp.input,
                           &(comp.external),
                           "EXT"))
    {
      fprintf(stderr,
              "Failed to create external assembly input file for %s\n",
              comp.input);
      return 1;
    }
    // capture the index of the component in the array
    c_i = GetNumVA_Component(globals->components);
    AppendVA_Component(globals->components, &comp);
  }

  // open the file, read UIDs into a var array, sort it, close the input file
  if(uids == NULL)
  {
#ifdef DEBUG1
    fprintf(stderr, "Creating array of %d recruiter uids\n", MODERATE_NUMBER);
#endif
    uids = CreateVA_CDS_UID_t(MODERATE_NUMBER);
    if(uids == NULL)
    {
      fprintf(stderr, "Failed to allocate array of %d uids\n",
              MODERATE_NUMBER);
      return 1;
    }
  }
  ResetVA_CDS_UID_t(uids);
  // EnableRangeVA_CDS_UID_t(uids, MODERATE_NUMBER);

  // set up or reinitialize reusable frag/mate iid array
  if(frag_iids == NULL)
  {
#ifdef DEBUG1
    fprintf(stderr, "Creating array of %d fragment iids\n", BIG_NUMBER);
#endif
    frag_iids = CreateVA_CDS_IID_t(BIG_NUMBER);
    if(frag_iids == NULL)
    {
      fprintf(stderr, "Failed to allocate array of %d iids\n",
              BIG_NUMBER);
      return 1;
    }
  }
  ResetVA_CDS_IID_t(frag_iids);
  // EnableRangeVA_CDS_IID_t(frag_iids, BIG_NUMBER);

  // set up/reinitialize reusable distance iid array
  if(dist_iids == NULL)
  {
#ifdef DEBUG1
    fprintf(stderr, "Creating array of %d distance iids\n", MODERATE_NUMBER);
#endif
    dist_iids = CreateVA_CDS_IID_t(MODERATE_NUMBER);
    if(dist_iids == NULL)
    {
      fprintf(stderr, "Failed to allocate array of %d iids\n",
              MODERATE_NUMBER);
      return 1;
    }
  }
  ResetVA_CDS_IID_t(dist_iids);
  // EnableRangeVA_CDS_IID_t(dist_iids, MODERATE_NUMBER);

  // open the file
  if((fp = fopen(component_filename, "r")) == NULL)
  {
    fprintf(stderr, "Failed to open component definition file %s\n",
            component_filename);
    return 1;
  }

  while(fgets( string, STRING_LENGTH, fp))
  {
    CDS_UID_t rec_uid = STR_TO_UID(string, NULL, 10);
    Recruiter * rec;
    CDS_IID_t i;

    if((rec = (Recruiter *) LookupInHashTable(globals->uid_ht,
                                              rec_uid))==NULL)
    {
      fprintf(stderr,
              "Failed to look up recruiting entity " F_UID " in hashtable\n",
              rec_uid);
      fprintf(stderr, "Continuing...\n");
      continue;
    }

    // add frag iids, mates & distances to temp arrays
    for(i = rec->start; i < rec->start + rec->num; i++)
    {
      iid = GetVA_CDS_IID_t(globals->rec_frags, i);
      AppendVA_CDS_IID_t(frag_iids, iid);
      if(globals->mates_dists[*iid].mate != 0)
      {
        AppendVA_CDS_IID_t(dist_iids,
                            &(globals->mates_dists[*iid].dist));
      }
    }
  }
  fclose(fp);

  // append unique iids to the master lists
  AppendToMaster(globals->frags, frag_iids, c_i);
  AppendToMaster(globals->dists, dist_iids, c_i);

  return 0;
}


int compare_iid_pair(const IIDPair * a, const IIDPair * b)
{
  if(a->iid1 > b->iid1) return 1;
  if(a->iid1 < b->iid1) return -1;
  return 0;
}





int WriteComponents(Globalsp globals)
{
  cds_uint32 i;
  CDS_IID_t last_iid;
  IIDPair * iid_pair;
  GenericMesg  gen;
  DistanceMesg dst;
  FragMesg     frg;
  LinkMesg     lkg;
  char     src[STRING_LENGTH];
  char     seq[AS_READ_MAX_LEN];
  char     qvs[AS_READ_MAX_LEN];
  Component  * comp;
  int write_lkg;
  
#ifdef DEBUG1
  fprintf(stderr, "Sorting master IID-component arrays by IID\n");
#endif
  
  // sort all the master arrays by IID
  qsort(globals->frags->Elements,
        GetNumVA_IIDPair(globals->frags),
        sizeof(IIDPair),
        (int (*) (const void *, const void *)) compare_iid_pair);
  qsort(globals->dists->Elements,
        GetNumVA_IIDPair(globals->dists),
        sizeof(IIDPair),
        (int (*) (const void *, const void *)) compare_iid_pair);

#ifdef DEBUG1
  fprintf(stderr, "Finished sorting master IID-component arrays by IID\n");
#endif

  // DISTANCES
#ifdef DEBUG1
  fprintf(stderr, "Looping over distances to write\n");
#endif
  gen.t = MESG_DST;
  gen.m = &dst;
  dst.action = AS_ADD;
  for(last_iid = 0, i = 0; i < GetNumVA_IIDPair(globals->dists); i++)
  {
    iid_pair = GetVA_IIDPair(globals->dists, i);
    if(iid_pair->iid1 != last_iid)
    {
      GateKeeperDistanceRecord gkpdr;
      getIndexStore( globals->s_set->gkp_store.dstStore,
                     iid_pair->iid1,
                     &gkpdr );
      dst.eaccession = gkpdr.UID;
      dst.mean       = gkpdr.mean;
      dst.stddev     = gkpdr.stddev;
      last_iid = iid_pair->iid1;
    }
    // write to the iid2 component
    comp = GetVA_Component(globals->components, iid_pair->iid2);
    WriteProtoMesg_AS(comp->internal.fp, &gen);
  }
#ifdef DEBUG1
  fprintf(stderr, "Finished writing distances\n");
#endif


  // READS
#ifdef DEBUG1
  fprintf(stderr, "Looping over reads to write\n");
#endif
  sprintf( src, " " );
  frg.source = src;
  frg.sequence = seq;
  frg.quality = qvs;
  
  dst.action = AS_ADD;
  write_lkg = 0;
  for(last_iid = 0, i = 0; i < GetNumVA_IIDPair(globals->frags); i++)
  {
    iid_pair = GetVA_IIDPair(globals->frags, i);
    if(iid_pair->iid1 != last_iid)
    {
      if(PopulateFragment(&frg,
                          globals->s_set->frg_store,
                          iid_pair->iid1,
                          globals->rs))
      {
        fprintf(stderr, "Failed to populate fragment " F_IID "\n",
                iid_pair->iid1);
        return 1;
      }

      if(globals->mates_dists[iid_pair->iid1].mate != 0 &&
         globals->mates_dists[iid_pair->iid1].mate < iid_pair->iid1)
      {
        lkg.action = AS_ADD;
        lkg.type = AS_MATE;
        lkg.entry_time = frg.entry_time;
        lkg.link_orient = globals->mates_dists[iid_pair->iid1].orient;
        lkg.frag1 = frg.eaccession;
        lkg.frag2 =
          GetFragmentUID(globals->s_set->gkp_store,
                         globals->mates_dists[iid_pair->iid1].mate);
        lkg.distance =
          GetDistanceUID(globals->s_set->gkp_store,
                         globals->mates_dists[iid_pair->iid1].dist);
        write_lkg = 1;
      }
      else
        write_lkg = 0;
      last_iid = iid_pair->iid1;
    }
    
    // write to the iid2 component
    comp = GetVA_Component(globals->components, iid_pair->iid2);
    gen.t = MESG_FRG;
    gen.m = &frg;
    WriteProtoMesg_AS(comp->internal.fp, &gen);
    if(write_lkg)
    {
      gen.t = MESG_LKG;
      gen.m = &lkg;
      WriteProtoMesg_AS(comp->internal.fp, &gen);
    }
  }
#ifdef DEBUG1
  fprintf(stderr, "Finished writing reads\n");
#endif

  // close the files
#ifdef DEBUG1
  fprintf(stderr, "Looping over component files to close\n");
#endif
  for(i = 0; i < GetNumVA_Component(globals->components); i++)
  {
    Component * comp = GetVA_Component(globals->components, i);
    if(comp->internal.fp)
      fclose(comp->internal.fp);
    if(comp->external.fp)
      fclose(comp->external.fp);
  }
#ifdef DEBUG1
  fprintf(stderr, "Done closing component files\n");
#endif

#ifdef DEBUG
  fprintf(stderr, "Finished writing set of components\n");
#endif
  
  return 0;
}


void ResetComponents(Globalsp globals)
{
  ResetVA_Component(globals->components);
  ResetVA_IIDPair(globals->dists);
  ResetVA_IIDPair(globals->frags);
}


int GenerateComponentInputs(Globalsp globals)
{
  FILE * fp_in;
  char   component_filename[STRING_LENGTH];
  int    i;

  // assumption that the output directory exists...
#ifdef DEBUG
  fprintf(stderr,
          "Generating component input files from definitions list %s\n",
          globals->input_filename);
#endif

  // create variable arrays
  globals->dists = CreateVA_IIDPair(0);
  globals->frags = CreateVA_IIDPair(0);
  if(globals->dists == NULL ||
     globals->frags == NULL)
  {
    fprintf( stderr, "Failed to create empty variable arrays.\n" );
    return 1;
  }
  
  if( (fp_in = fopen( globals->input_filename, "r" )) == NULL )
  {
    fprintf( stderr, "Failed to open input filename %s\n",
             globals->input_filename );
    return 1;
  }
  
  /* loop over components in blocks limited by number of available
     file descriptors
  */
  for( i = 0; fgets( component_filename, STRING_LENGTH, fp_in ); i++ )
  {
    if(component_filename[strlen(component_filename)-1] == '\n')
      component_filename[strlen(component_filename)-1] = '\0';
    if(i > 0 && i % MAX_COMPONENTS == 0)
    {
      /* done with last set of MAX_COMPONENTS components
         sort, read once, write many & reset data structures
      */
#ifdef DEBUG
      fprintf(stderr, "Read %d component definitions. Writing components\n",
              i);
#endif
      if(WriteComponents(globals))
      {
        fprintf(stderr, "Failed to write set of component files\n");
        return 1;
      }
      ResetComponents(globals);
    }
    if(AddComponent(globals, component_filename))
    {
      fprintf(stderr, "Failed to add component to components data set\n");
      return 1;
    }
  }
#ifdef DEBUG
  fprintf(stderr, "Read %d component definitions. Writing final components\n",
          i);
#endif
  if(WriteComponents(globals))
  {
    fprintf(stderr, "Failed to write set of component files\n");
    return 1;
  }

  fclose( fp_in );
  return 0;
}


int SetUpMatesDistsArray(Globalsp globals)
{
  StoreStat stats;
  GateKeeperLinkRecord gkpl;
  cds_int64 i;
  CDS_IID_t iid;

#ifdef DEBUG
  fprintf(stderr, "Setting up mates & distances array\n");
#endif
  
  statsStore(globals->s_set->gkp_store.frgStore, &stats);
  globals->mates_dists = calloc(stats.lastElem + 1, sizeof(IIDOrient));
  if(globals->mates_dists == NULL)
  {
    fprintf(stderr, "Failed to allocate array of " F_S64 " iid pairs\n",
            stats.lastElem + 1);
    return 1;
  }

  globals->num_mates = stats.lastElem + 1;
  for(i = 1; i <= stats.lastElem; i++)
  {
    // if not set yet
    if(globals->mates_dists[i].mate == 0)
    {
      // if has a mate, set both & the distance
      if((iid = FindMate(globals->s_set->gkp_store, i, &gkpl)) != 0)
      {
        globals->mates_dists[i].mate = iid;
        globals->mates_dists[i].dist = gkpl.distance;
        globals->mates_dists[i].orient = getLinkOrientation(&gkpl);

        globals->mates_dists[iid].mate = i;
        globals->mates_dists[iid].dist = gkpl.distance;
        globals->mates_dists[iid].orient = globals->mates_dists[i].orient;
      }
    }
  }

#ifdef DEBUG
  fprintf(stderr, "Finished setting up mates & distances array\n");
#endif
  
  return 0;
}


int AddIIDsToRecFrags(Globalsp globals,
                      VA_TYPE(CDS_IID_t) * iids,
                      CDS_UID_t last_uid)
{
  Recruiter     rec;
  int           i;
  CDS_IID_t    iid;
  PHashValue_AS value;

  // sort the iids
  if(GetNumVA_CDS_IID_t(iids) > 1)
    qsort(iids->Elements, GetNumVA_CDS_IID_t(iids), sizeof(CDS_IID_t),
          (int (*) (const void *, const void *)) compareIIDs);
  
  // iterate through them & add to the big set
  iid = 0;
  rec.uid = last_uid;
  rec.start = GetNumVA_CDS_IID_t(globals->rec_frags);
  rec.num = 0;
  for(i = 0; i < GetNumVA_CDS_IID_t(iids); i++)
  {
    if(iid != *GetVA_CDS_IID_t(iids, i))
    {
      iid = *GetVA_CDS_IID_t(iids, i);
      AppendVA_CDS_IID_t(globals->rec_frags, &iid);
      rec.num++;
    }
  }
  
  rec.iid = 0;

  // insert the scaffold/bland UID into the hashtable
  if(InsertInHashTable(globals->uid_ht, rec.uid, (char *) &rec))
  {
    fprintf( stderr, "Failed to insert UID into hashtable\n");
    return 1;
  }
  globals->num_ht++;

  return 0;
}


int ProcessUIDPairsFile(Globalsp globals, char * filename)
{
  FILE * fp;
  char   line[STRING_LENGTH];
  CDS_UID_t last_uid = 0;
  CDS_UID_t new_uid;
  CDS_UID_t frag_uid;
  VA_TYPE(CDS_IID_t) * iids;
  PHashValue_AS value;
  cds_uint32 num_recruiters = 0;

  // create a temporary iid array
#ifdef DEBUG1
  fprintf(stderr, "Allocating array of %d iids\n", MODERATE_NUMBER);
#endif
  iids = CreateVA_CDS_IID_t(MODERATE_NUMBER);
  if(iids == NULL)
  {
    fprintf(stderr, "Failed to allocate variable array of %d iids\n",
            MODERATE_NUMBER);
    return 1;
  }
  ResetVA_Component(iids);
  // EnableRangeVA_cds_uint16(iids, MODERATE_NUMBER);

#ifdef DEBUG
  fprintf(stderr, "Reading recruiter UID - fragment UID file %s\n",
          filename);
#endif
  
  // open scaffold file
  if((fp = fopen(filename, "r")) == NULL)
  {
    fprintf(stderr, "Failed to open scaffold UID - frag UID file: %s\n",
            filename);
    return 1;
  }

  // zip through pairs
  while(fgets(line, STRING_LENGTH, fp))
  {
    // if new scaffold, append last to master array
    // reset temporary var array & start adding to it
    sscanf(line, F_UID " " F_UID, &new_uid, &frag_uid);
    if(new_uid != last_uid && last_uid != 0)
    {
#ifdef DEBUG2
      fprintf(stderr, "Adding recruiter " F_UID " frags to master set\n",
              last_uid);
#endif
      
      if(AddIIDsToRecFrags(globals, iids, last_uid))
      {
        fprintf(stderr, "Failed to add IIDs to recruited frags set\n");
        return 1;
      }
      num_recruiters++;

      // prepare for next scaffold
      ResetVA_Component(iids);
      
      if(LookupInHashTable(globals->uid_ht, new_uid))
      {
        fprintf(stderr, "Error! Scaffold UIDs in file are not sorted!\n");
        return 1;
      }
    }

    last_uid = new_uid;
    if( HASH_SUCCESS !=
        LookupTypeInPHashTable_AS( globals->s_set->gkp_store.hashTable, 
                                   UID_NAMESPACE_AS,
                                   frag_uid,
                                   AS_IID_FRG, 
                                   FALSE,
                                   stderr,
                                   &value ) )
    {
      fprintf(stderr, "Fragment UID " F_UID " of recruiter " F_UID " is not in store!\n",
              frag_uid, last_uid);
      continue;
    }
    else
    {
      AppendVA_CDS_IID_t(iids, &(value.IID));
    }
  }

#ifdef DEBUG2
  fprintf(stderr, "Adding final recruiter " F_UID " frags to master set\n",
          last_uid);
#endif
  if(AddIIDsToRecFrags(globals, iids, last_uid))
  {
    fprintf(stderr, "Failed to add IIDs to recruited frags set\n");
    return 1;
  }
  num_recruiters++;
  
  fclose(fp);
  DeleteVA_CDS_IID_t(iids);

#ifdef DEBUG
  fprintf(stderr, "Finished processing UID pairs file %s\n", filename);
  fprintf(stderr, "%u recruiters processed\n", num_recruiters);
#endif
          
  return 0;
}


int ReadUIDPairsFiles(Globalsp globals)
{
#ifdef DEBUG
  fprintf(stderr, "Reading UID pairs files\n");
#endif
#ifdef DEBUG1
  fprintf(stderr, "Creating recruiter uid hashtable and frag iid array\n");
#endif
  
  // create the hashtable
  if((globals->uid_ht = CreateHashTable(NUM_RECRUITERS,
                                        sizeof(Recruiter))) == NULL)
  {
    fprintf( stderr, "Failed to create hashtable for %d recruiters\n",
             NUM_RECRUITERS );
    return 1;
  }

  // create the huge fragment iid array pointed to by recruiter objects
  globals->rec_frags = CreateVA_CDS_IID_t(NUM_REC_FRAGS);
  if(globals->rec_frags == NULL)
  {
    fprintf(stderr, "Failed to allocate variable array of %d iids\n",
            NUM_REC_FRAGS);
    return 1;
  }
  ResetVA_Component(globals->rec_frags);
  // EnableRangeVA_cds_uint16(globals->rec_frags, NUM_REC_FRAGS);

  if(globals->scaffold_filename != NULL &&
     ProcessUIDPairsFile(globals, globals->scaffold_filename, 0))
  {
    fprintf(stderr, "Failed to read & process scaffold UID file: %s\n",
            globals->scaffold_filename);
    return 1;
  }
    
  return 0;
}


int ReadInputStore(Globalsp globals)
{
  FILE * fp;
  char filename[STRING_LENGTH];
  Recruiter * recruiters;
  cds_uint32 i;

#ifdef DEBUG
  fprintf(stderr, "Reading input store %s\n", globals->input_store);
#endif

  // read hash table elements - stupid, but quick implementation
  sprintf(filename, "%s/%s", globals->input_store, HASH_FILE);
#ifdef DEBUG1
  fprintf(stderr, "Reading hash data from %s\n", filename);
#endif
  if((fp = fopen(filename, "r")) == NULL)
  {
    fprintf(stderr, "Failed to open hash file %s for reading\n", filename);
    return 1;
  }
  if(fread(&(globals->num_ht), sizeof(globals->num_ht), 1, fp) != 1)
  {
    fprintf(stderr, "Failed to read %d items\n", 1);
    fclose(fp);
    return 1;
  }
  recruiters = (Recruiter *) malloc( sizeof(Recruiter) * globals->num_ht);
  if(recruiters == NULL)
  {
    fprintf(stderr, "Failed to allocate temporary array of %u recruiters.\n",
            globals->num_ht);
    fclose(fp);
    return 1;
  }
  if((globals->uid_ht = CreateHashTable(globals->num_ht,
                                        sizeof(Recruiter))) == NULL)
  {
    fprintf( stderr, "Failed to create hashtable for %u recruiters\n",
             globals->num_ht);
    fclose(fp);
    return 1;
  }
  if(fread(recruiters, sizeof(Recruiter), globals->num_ht, fp) !=
     globals->num_ht)
  {
    fprintf(stderr, "Failed to read %u recruiters from file %s\n",
            globals->num_ht, filename);
    fclose(fp);
    return 1;
  }
  for(i = 0; i < globals->num_ht; i++)
  {
    if(InsertInHashTable(globals->uid_ht,
                         recruiters[i].uid,
                         (char *) &(recruiters[i])))
    {
      fprintf( stderr, "Failed to insert UID into hashtable\n");
      fclose(fp);
      return 1;
    }
    
  }
  free(recruiters);
  fclose(fp);
  
  // read var array of fragments
  sprintf(filename, "%s/%s", globals->input_store, FRAGS_FILE);
#ifdef DEBUG1
  fprintf(stderr, "Reading fragment data from %s\n", filename);
#endif
  if((fp = fopen(filename, "r")) == NULL)
  {
    fprintf(stderr, "Failed to open frags file %s for reading\n", filename);
    fclose(fp);
    return 1;
  }
  globals->rec_frags = CreateVA_CDS_IID_t(1);
  LoadFromFile_VA(fp, globals->rec_frags, "CDS_IID_t", 0);
  fclose(fp);

  // read array of mates
  sprintf(filename, "%s/%s", globals->input_store, MATES_FILE);
#ifdef DEBUG1
  fprintf(stderr, "Reading mates-dists data from %s\n", filename);
#endif
  if((fp = fopen(filename, "r")) == NULL)
  {
    fprintf(stderr, "Failed to open mates file %s for reading\n", filename);
    return 1;
  }
  if(fread(&(globals->num_mates), sizeof(globals->num_mates), 1, fp) != 1)
  {
    fprintf(stderr, "Failed to read %d items\n", 1);
    fclose(fp);
    return 1;
  }
  globals->mates_dists =
    (IIDOrient *) malloc(globals->num_mates * sizeof(IIDOrient));
  if(globals->mates_dists == NULL)
  {
    fprintf(stderr, "Failed to allocate array of %u mates/dists\n",
            globals->num_mates);
    fclose(fp);
    return 1;
  }
  if(fread(globals->mates_dists, sizeof(IIDOrient), globals->num_mates, fp) !=
     globals->num_mates)
  {
    fprintf(stderr, "Failed to read %u items\n", globals->num_mates);
    fclose(fp);
    return 1;
  }
  fclose(fp);

#ifdef DEBUG
  fprintf(stderr, "Finished reading input store\n");
#endif

  return 0;
}


int WriteOutputStore(Globalsp globals)
{
  FILE * fp;
  char filename[STRING_LENGTH];
  Recruiter * recruiter;
  
#ifdef DEBUG
  fprintf(stderr, "Writing output store %s\n", globals->output_store);
#endif

#ifdef DEBUG1
  fprintf(stderr, "Creating store %s\n", globals->output_store);
#endif
  sprintf(filename, "mkdir %s", globals->output_store);
  if(system(filename))
  {
    fprintf(stderr, "Failed to create output store %s\n",
            globals->output_store);
    fprintf(stderr, "Assuming directory already exists....\n");
  }
  
  // write hash table elements - stupid, but quick implementation
  sprintf(filename, "%s/%s", globals->output_store, HASH_FILE);
#ifdef DEBUG1
  fprintf(stderr, "Writing hash data to %s\n", filename);
#endif
  if((fp = fopen(filename, "w")) == NULL)
  {
    fprintf(stderr, "Failed to open hash file %s for writing\n", filename);
    return 1;
  }
  // space filler for now
  if(fwrite(&(globals->num_ht), sizeof(globals->num_ht), 1, fp) != 1)
  {
    fprintf(stderr, "Failed to write %d items\n", 1);
    return 1;
  }
  RewindHashTable(globals->uid_ht);
  while((recruiter = (Recruiter *) GetNextHashTableRet(globals->uid_ht)) !=
         NULL)
  {
    if(fwrite(recruiter, sizeof(Recruiter), 1, fp) != 1)
    {
      fprintf(stderr, "Failed to write %d items\n", 1);
      return 1;
    }
  }
  fclose(fp);
  
  // write var array of fragments
  sprintf(filename, "%s/%s", globals->output_store, FRAGS_FILE);
#ifdef DEBUG1
  fprintf(stderr, "Writing fragment data to %s\n", filename);
#endif
  if((fp = fopen(filename, "w")) == NULL)
  {
    fprintf(stderr, "Failed to open frags file %s for writing\n", filename);
    return 1;
  }
  if(!CopyToFile_VA(globals->rec_frags, fp))
  {
    fprintf(stderr, "Failed to write variable array of " F_SIZE_T " items to file %s\n",
            GetNumVA_CDS_IID_t(globals->rec_frags), filename);
    return 1;
  }
  fclose(fp);

  // write array of mates
  sprintf(filename, "%s/%s", globals->output_store, MATES_FILE);
#ifdef DEBUG1
  fprintf(stderr, "Writing mates-dists data to %s\n", filename);
#endif
  if((fp = fopen(filename, "w")) == NULL)
  {
    fprintf(stderr, "Failed to open mates file %s for writing\n", filename);
    return 1;
  }
  if(fwrite(&(globals->num_mates), sizeof(globals->num_mates), 1, fp) !=
     1)
  {
    fprintf(stderr, "Failed to write %d items\n", 1);
    return 1;
  }
  if(fwrite(globals->mates_dists, sizeof(IIDOrient), globals->num_mates, fp) !=
     globals->num_mates)
  {
    fprintf(stderr, "Failed to write %u items\n", globals->num_mates);
    return 1;
  }
  fclose(fp);

#ifdef DEBUG
  fprintf(stderr, "Finished writing output store\n");
#endif
  
  return 0;
}


int main( int argc, char ** argv )
{
  Globals * globals = CreateGlobals();

  // parse parameters & exit on any problem
  ProcessCommandLine( globals, argc, argv );
  
  // open the stores
  if( (globals->s_set = OpenStores( globals->gkp_store_name,
                                    globals->frg_store_name )) == NULL )
  {
    fprintf( stderr, "Failed to open one or more stores.\n" );
    return 1;
  }
  
  // allocate a reusable readstruct;
  if( (globals->rs = new_ReadStruct()) == NULL )
  {
    fprintf( stderr, "Failed to allocate readStruct.\n" );
    return 1;
  }

  if( (globals->uid_fetcher = CreateUIDInteractor( 1 )) == NULL )
  {
    fprintf( stderr, "Failed to connect to UID generator.\n" );
    return 1;
  }

  // read in the uid pairs
  if(globals->scaffold_filename)
  {
    if(ReadUIDPairsFiles(globals))
    {
      fprintf( stderr, "Failed to read in UID pairs files\n" );
      return 1;
    }
    if(SetUpMatesDistsArray(globals))
    {
      fprintf(stderr, "Failed to set up fragment mates array\n" );
      return 1;
    }
  }
  else
  {
    if(ReadInputStore(globals))
    {
      fprintf( stderr, "Failed to read input store %s\n",
               globals->input_store);
      return 1;
    }
  }

  // write output store
  if(globals->output_store && WriteOutputStore(globals))
  {
    fprintf( stderr, "Failed to write output store %s\n",
             globals->output_store );
    return 1;
  }

  if(GenerateComponentInputs(globals))
  {
    fprintf( stderr, "Failed to generate component input files\n");
    return 1;
  }

  FreeGlobals(globals);

  return 0;
}
