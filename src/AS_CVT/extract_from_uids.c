
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
static char CM_ID[] = "$Id: extract_from_uids.c,v 1.2 2004-09-23 20:25:21 mcschatz Exp $";

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
#include <unistd.h> /* man 3 getopt */
#include <memory.h>
#include <assert.h>
#include <math.h>

// Project header files
#include "AS_global.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_PHash.h"
#include "AS_UTL_Var.h"
#include "AS_ALN_aligners.h"
#include "AS_CVT_hashtable.h"
#include "AS_CVT_uid_interface.h"
#include "AS_GKP_include.h"


#define STRING_LENGTH 1024


// structure to hold an array of strings - used for filenames
typedef struct
{
  int     num_strings;
  char ** strings;
} StringSet;
typedef StringSet * StringSetp;


// structure to hold store handles
typedef struct
{
//  GateKeeperStore rb_gkp_store;
//  FragStoreHandle r_bac_store;
//  GateKeeperStore rf_gkp_store;
//  FragStoreHandle r_frg_store;
  GateKeeperStore g_gkp_store;
  FragStoreHandle g_frg_store;
} StoreSet;
typedef StoreSet * StoreSetp;


// structure to hold command-line & other 'global' items
typedef struct
{
  char           * program_name;        // name of this program
  cds_uint64       batch_int_uid;       // UID for the batch generated
  cds_uint64       batch_ext_uid;       // UID for the batch generated
  char             out_file_name[1024]; // output filename
  FILE           * fp_int_out;          // pointer to internal data output file
  FILE           * fp_ext_out;          // pointer to internal data output file
/*
  char           * list_file_name;      // file listing .cns filenames
  char           * r_frg_store_name;    // regional fragstore name
  char           * r_bac_store_name;    // regional bactigstore name
*/
  char           * g_gkp_store_name;    // grande gkpstore name
  char           * g_frg_store_name;    // grande fragstore name
  int              include_mates;       // 1 if including unrecruited mates
  UIDInteractorp   uid_fetcher;         // object that gets real UIDs
  StringSetp       f_set;               // list of files
  StoreSetp        s_set;               // set of stores
  ReadStructp      rs;                  // reusable readStruct
} Globals;
typedef Globals * Globalsp;


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


/****************************************************************************
 * StringSet code
 ****************************************************************************/
// free a set of strings - fully or partially allocated
void FreeStringSet( StringSetp string_set )
{
  int i;
  
  if( string_set )
  {
    if( string_set->strings )
    {
      for( i = 0; i < string_set->num_strings; i++ )
      {
        if( string_set->strings[i] )
          free( string_set->strings[i] );
      }
      free( string_set->strings );
    }
    free( string_set );
  }
}

// allocate a set of strings
StringSetp AllocateStringSet( int num_strings )
{
  StringSetp string_set = NULL;

  if( (string_set = (StringSetp) calloc( 1, sizeof( StringSet ) )) != NULL )
  {
    if( (string_set->strings =
         (char **) calloc( num_strings, sizeof( char * ) )) != NULL )
    {
      for( string_set->num_strings = 0;
           string_set->num_strings < num_strings;
           string_set->num_strings++ )
      {
        if( (string_set->strings[string_set->num_strings] =
             calloc( STRING_LENGTH, sizeof( char ) )) == NULL )
        {
          FreeStringSet( string_set );
          return NULL;
        }
      }
    }
    else
    {
      FreeStringSet( string_set );
      string_set = NULL;
    }
  }
  return string_set;
}


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
 * Read a list of files (strings) from a file
 ****************************************************************************/
StringSetp ReadListOfFiles( int argc, char ** argv, int first )
{
  StringSetp f_set = NULL;

  if( argc - first )
  {
    int i;

    f_set = AllocateStringSet( argc - first );
    if( f_set == NULL )
    {
      fprintf( stderr, "Failed to allocate set of %d files\n", argc - first );
      return NULL;
    }

    for( i = 0; i < f_set->num_strings; i++ )
    {
      strcpy( f_set->strings[i], argv[i + first] );
    }
  }
  return f_set;
}



/****************************************************************************
 * Populate a FragMesg from a fragStore given an iid
 ****************************************************************************/
int PopulateFragment( FragMesg * f, FragStoreHandle fs, CDS_CID_t iid,
                      ReadStructp rs )
{
  cds_uint32 bgn, end;
  Locale_ID eilocale;

  if( getFragStore( fs, iid, FRAG_S_ALL, rs ) )
  {
    fprintf( stderr, "Failed to get fragment " F_IID " data\n", iid );
    return 1;
  }

  // populate the FRG message fields
  f->action = AS_ADD;
  getAccID_ReadStruct( rs, &(f->eaccession) );
  getReadType_ReadStruct( rs, &(f->type) );
  // NOTE: The locale IID is stored in the fragStore, NOT the UID
  getLocID_ReadStruct( rs, &eilocale );
  f->ilocale = (IntLocale_ID) eilocale;
  // elocale is unavailable
  // eseq_id is unavailable
  // ebactig_id is unavailable
  getLocalePos_ReadStruct( rs, &bgn, &end );
  f->locale_pos.bgn = bgn;
  f->locale_pos.end = end;
  getEntryTime_ReadStruct (rs, &(f->entry_time) );
  getClearRegion_ReadStruct( rs, &bgn, &end, READSTRUCT_LATEST );
  f->clear_rng.bgn = bgn;
  f->clear_rng.end = end;
  getSource_ReadStruct( rs, f->source, STRING_LENGTH );
  getSequence_ReadStruct( rs, f->sequence, f->quality, AS_READ_MAX_LEN );
  getReadIndex_ReadStruct( rs, &(f->iaccession) );
  // f->screened is moot
  // ilocale is unavailable
  // iseq_id is unavailable
  // ibactig_id is unavailable
  
  return 0;
}

void FreeFragment( FragMesg f) {
  if ( f.source != NULL ) free(f.source);
  if ( f.sequence != NULL ) free(f.sequence);
  if ( f.quality != NULL ) free(f.quality);
}

/****************************************************************************
 * Code for using several stores
 ****************************************************************************/
void CloseStores( StoreSetp s_set )
{
  if( s_set )
  {
    /*
    CloseGateKeeperStore( &(s_set->rb_gkp_store) );
    CloseGateKeeperStore( &(s_set->rf_gkp_store) );
    if( s_set->r_frg_store )
      closeFragStore( s_set->r_frg_store );
    if( s_set->r_bac_store )
      closeFragStore( s_set->r_bac_store );
    */
    CloseGateKeeperStore( &(s_set->g_gkp_store) );
    if( s_set->g_frg_store )
      closeFragStore( s_set->g_frg_store );
    
    free( s_set );
  }
}


StoreSetp OpenStores( char * g_gkp_store_name,
                      char * g_frg_store_name )

{
  StoreSetp s_set;

  // allocate the containing structure
  if( (s_set = (StoreSetp) calloc( 1, sizeof( StoreSet ) )) == NULL )
  {
    fprintf( stderr, "Failed to allocate StoreSet\n" );
    return NULL;
  }

  /*
  // Open regional stores
  {
    // open the regional gatekeeper stores
    InitGateKeeperStore( &(s_set->rb_gkp_store), r_bac_store_name );
    OpenReadOnlyGateKeeperStore( &(s_set->rb_gkp_store) );
    InitGateKeeperStore( &(s_set->rf_gkp_store), r_frg_store_name );
    OpenReadOnlyGateKeeperStore( &(s_set->rf_gkp_store) );
    
    // open the regional fragment store
    s_set->r_frg_store = openFragStore( r_frg_store_name, "r" );
    if( s_set->r_frg_store == NULLSTOREHANDLE )
    {
      fprintf( stderr, "Failed to open regional frag store %s\n",
               r_frg_store_name );
      CloseStores( s_set );
      return NULL;
    }
    
    // open the regional fragment store
    s_set->r_bac_store = openFragStore( r_bac_store_name, "r" );
    if( s_set->r_bac_store == NULLSTOREHANDLE )
    {
      fprintf( stderr, "Failed to open regional bactig store %s\n",
               r_bac_store_name );
      CloseStores( s_set );
      return NULL;
    }
  }
  */

  // Open grande stores
  {
    // open the regional gatekeeper store
    InitGateKeeperStore( &(s_set->g_gkp_store), g_gkp_store_name );
    OpenReadOnlyGateKeeperStore( &(s_set->g_gkp_store) );
    
    // open the regional fragment store
    s_set->g_frg_store = openFragStore( g_frg_store_name, "r" );
    if( s_set->g_frg_store == NULLSTOREHANDLE )
    {
      fprintf( stderr, "Failed to open grande frag store %s\n",
               g_frg_store_name );
      CloseStores( s_set );
      return NULL;
    }
  }
  
  return s_set;
}


int ConstructAndWriteUBACMesg( FILE * fp,
                               StoreSetp ss,
                               Locale_ID distuid,
                               Locale_ID esequid,
                               GateKeeperLocaleRecord gkpl)
{
  GenericMesg   gen;
  BacMesg       bac;
  int           i;
  
  bac.action = AS_ADD;
  bac.ebac_id = gkpl.UID;
  bac.eseq_id = esequid;
  bac.entry_time = time(0);
  bac.source = " ";
  bac.elength = distuid;

  if( gkpl.type == AS_FBAC )
  {
    bac.type = AS_FBAC;
    bac.num_bactigs = 0;
    bac.bactig_list = NULL;
  }
  else
  {
    bac.type = AS_UBAC;
    bac.num_bactigs = gkpl.numBactigs;
    bac.bactig_list = (BactigMesg *) calloc( gkpl.numBactigs,
                                             sizeof( BactigMesg ) );
    if( bac.bactig_list == NULL )
    {
      fprintf( stderr, "Error: Failed to allocate list of %d bactigs\n",
               gkpl.numBactigs );
      return 1;
    }
    i = 0;
    while( i < bac.num_bactigs )
    {
      GateKeeperBactigRecord gkpb;
      
      if( getGateKeeperBactigStore(ss->g_gkp_store.btgStore,
                                   gkpl.firstBactig+i,
                                   &gkpb ) )
      {
        fprintf( stderr,
                 "Error: Failed to access bactig %d in bactigstore!\n",
                 gkpl.firstBactig+i );
        return 1;
      }

      bac.bactig_list[i].eaccession = gkpb.UID;
      bac.bactig_list[i].length = gkpb.length;
      i++;
    }
  }

  // write BAC
  gen.t = MESG_BAC;
  gen.m = &bac;
  WriteProtoMesg_AS( fp, &gen );

  if( bac.bactig_list )
    free(bac.bactig_list);
  return 0;
}


int ConstructAndWriteEBACMesg( FILE * fp,
                               FragMesg * frag,
                               cds_uint64 elength )
{
  GenericMesg   gen;
  BacMesg       bac;
  
  bac.action = AS_ADD;
  bac.ebac_id = frag->elocale;
  bac.type = AS_EBAC;
  bac.entry_time = frag->entry_time;
  bac.source = " ";
  bac.elength = elength;
  bac.num_bactigs = 0;
  bac.bactig_list = NULL;
  
  // write BAC
  gen.t = MESG_BAC;
  gen.m = &bac;
  WriteProtoMesg_AS( fp, &gen );
  return 0;
}


int WriteBATADT( Globalsp globals, int internal )
{
  GenericMesg  gen;
  BatchMesg    bat;
  AuditLine    adl;
  AuditMesg    adt;
  char       * comment;
  int          i;
  FILE *fp = (internal)?globals->fp_int_out:globals->fp_ext_out;

  // first message must be BAT with true UID;
  bat.name = globals->program_name;
  bat.created = time(0);
  bat.eaccession = (internal)?globals->batch_int_uid: globals->batch_ext_uid;
  bat.comment = "Grande dataset from fragment and BAC UIDs";
  gen.t = MESG_BAT;
  gen.m = &bat;
  WriteProtoMesg_AS( fp, &gen );
  
  // next message should/must be ADT
  adl.next = NULL;
  adl.name = globals->program_name;
  adl.complete = time(0);
  adl.version = "Beta";

  comment = calloc( globals->f_set->num_strings * 1024, sizeof( char ) );
  if( comment == NULL )
  {
    fprintf( stderr, "Failed to allocate comment line\n" );
    return 1;
  }
  sprintf( comment, "Based on UID list files:\n" );
  for(i = 0; i < globals->f_set->num_strings; i++ )
  {
    strcat( comment, globals->f_set->strings[i] );
    strcat( comment, "\n" );
  }
  adl.comment = comment;
  adt.list = &adl;
  gen.t = MESG_ADT;
  gen.m = &adt;
  WriteProtoMesg_AS( fp, &gen );

  free( comment );
  return 0;
}


int InitializeOutput( Globalsp globals )
{
  time_t seconds;
  struct tm timer;

  // get the time
  seconds = time( 0 );
  if( localtime_r( &seconds, &timer ) == NULL )
  {
    fprintf( stderr, "Failed to get time!\n" );
    return 1;
  }
  
  // get a Internal batch UID
  if( GetUID( globals->uid_fetcher, &(globals->batch_int_uid) ) )
  {
    fprintf( stderr, "Failed to get Internal batch UID.\n" );
    return 1;
  }
  // get a External batch UID
  if( GetUID( globals->uid_fetcher, &(globals->batch_ext_uid) ) )
  {
    fprintf( stderr, "Failed to get External batch UID.\n" );
    return 1;
  }
  
  // create an output filename
  sprintf( globals->out_file_name,
           "INT_%04d%02d%02d_%020" F_UIDP "_HUMAN_FRAGS.frg",
           timer.tm_year + 1900,
           timer.tm_mon + 1,
           timer.tm_mday,
           globals->batch_int_uid );

  // open the output file for writing
  globals->fp_int_out = fopen( globals->out_file_name, "w" );
  if( globals->fp_int_out == NULL )
  {
    fprintf( stderr, "Failed to open %s for writing.\n",
             globals->out_file_name );
    return 1;
  }
  fprintf( stdout, "Output filename is %s\n", globals->out_file_name );
  sprintf( globals->out_file_name,
           "EXT_%04d%02d%02d_%020" F_UIDP "_HUMAN_BACS.frg",
           timer.tm_year + 1900,
           timer.tm_mon + 1,
           timer.tm_mday,
           globals->batch_ext_uid );

  // open the output file for writing
  globals->fp_ext_out = fopen( globals->out_file_name, "w" );
  if( globals->fp_ext_out == NULL )
  {
    fprintf( stderr, "Failed to open %s for writing.\n",
             globals->out_file_name );
    return 1;
  }
  fprintf( stdout, "Output filename is %s\n", globals->out_file_name );
  if ( WriteBATADT( globals, 1) ) return 0;
  return( WriteBATADT( globals, 0) );
}


Locale_ID GetBACUID( GateKeeperStore gkp_store, IntLocale_ID iid )
{
  GateKeeperBactigRecord gkpf;
  GateKeeperLocaleRecord gkpl;
  CDS_CID_t bacID;
  
  if( getGateKeeperBactigStore( gkp_store.btgStore, iid, &gkpf) )
  {
    fprintf( stderr, "Failed to get fragment " F_IID " data from gatekeeperstore\n",
             iid );
    return 0;
  }

  bacID = gkpf.bacID;  
  if( getGateKeeperLocaleStore( gkp_store.locStore, bacID, &gkpl) )
  {
    fprintf( stderr, "Failed to get locale " F_IID " data from gatekeeperstore\n",
             bacID );
    return 0;
  }
  
  return gkpl.UID;
}

Locale_ID GetLocaleUID( GateKeeperStore gkp_store, IntLocale_ID iid )
{
  GateKeeperLocaleRecord gkpl;
  
  if( getGateKeeperLocaleStore( gkp_store.locStore, iid, &gkpl) )
  {
    fprintf( stderr, "Failed to get locale " F_IID " data from gatekeeperstore\n",
             iid );
    return 0;
  }
  
  return gkpl.UID;
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


int GetInputUIDs( Globalsp globals, HashTablep * rht, HashTablep * bht )
{
  int    fs_i;
  int    read_count = 0;
  int    bac_end_count = 0;
  int    bac_count = 0;
  int    missing_count = 0;
  IDSet  idset;

  // create read hashtable
  if( (*rht = CreateHashTable( 500000, sizeof( IDSet ) )) == NULL )
  {
    fprintf( stderr, "Failed to create read hashtable\n" );
    return 1;
  }

  // create a bac hashtable
  if( (*bht = CreateHashTable( 10000, sizeof( IDSet ) )) == NULL )
  {
    fprintf( stderr, "Failed to create BAC hashtable\n" );
    FreeHashTable( *rht );
    return 1;
  }

  // loop over input files
  read_count = bac_end_count = bac_count = 0;
  memset( &idset, 0, sizeof( IDSet ) );
  for( fs_i = 0; fs_i < globals->f_set->num_strings; fs_i++ )
  {
    FILE          * fp_in;
    char            line[STRING_LENGTH];

    fprintf( stdout, "Reading data from %s\n",
             globals->f_set->strings[fs_i] );
    
    // open the input file
    if( (fp_in = fopen( globals->f_set->strings[fs_i], "r" )) == NULL )
    {
      fprintf( stderr, "Failed to open file %s for reading\n",
               globals->f_set->strings[fs_i] );
      FreeHashTable( *rht );
      FreeHashTable( *bht );
      return 1;
    }
    
    // iterate through input messages
    while( fgets( line, STRING_LENGTH, fp_in ) )
    {
      PHashValue_AS value;

      // get the UID
      idset.uid = STR_TO_UID( line, NULL, 10 );

      // see if it's a BAC or fragment
      if( HASH_SUCCESS !=
          LookupTypeInPHashTable_AS( globals->s_set->g_gkp_store.hashTable, 
                                     UID_NAMESPACE_AS,
                                     idset.uid,
                                     AS_IID_FRAG, 
                                     FALSE,
                                     stderr,
                                     &value ) )
      {
        // see if it's a BAC
        if( HASH_SUCCESS !=
            LookupTypeInPHashTable_AS( globals->s_set->g_gkp_store.hashTable, 
                                       UID_NAMESPACE_AS,
                                       idset.uid,
                                       AS_IID_LOC,
                                       FALSE,
                                       stderr,
                                       &value ) )
        {
          fprintf( stderr,
                   F_UID " is not a fragment or BAC UID "
                   "in grande gatekeeper store\n",
                   idset.uid );
          missing_count++;
          continue;
        }
        else
        {
          // it's a BAC UID
          // add it if it's not already there
          if( ! LookupInHashTable( *bht, idset.uid ) )
          {
            if( InsertInHashTable( *bht, idset.uid, (char *) &idset ) )
            {
              fprintf( stderr,
                       "Failed to insert UID into locale hashtable\n" );
              fclose( fp_in );
              FreeHashTable( *rht );
              FreeHashTable( *bht );
              return 1;
            }
            bac_count ++;
          }
        }
      }
      else
      {
        // it's a fragment
        // get it's type
        // populate idset.iid, idset.type
        idset.iid = value.IID;
        idset.type = GetFragmentType( globals->s_set->g_gkp_store, idset.iid );
        
        // add it if it's not already present
        if( ! LookupInHashTable( *rht, idset.uid ) )
        {
          if( InsertInHashTable( *rht, idset.uid, (char *) &idset ) )
          {
            fprintf( stderr,
                     "Failed to insert UID into read hasthable\n" );
            fclose( fp_in );
            FreeHashTable( *rht );
            FreeHashTable( *bht );
            return 1;
          }
          read_count += (idset.type == AS_READ ||
                         idset.type == AS_EXTR ||
                         idset.type == AS_TRNR) ? 1 : 0;
          bac_end_count += (idset.type == AS_EBAC) ? 1 : 0;
        }
      }
    }
    fclose( fp_in );
  }
  
  fprintf( stdout, "%10d     reads referenced in input files\n",
           read_count );
  fprintf( stdout, "%10d  BAC ends referenced in input files\n",
           bac_end_count );
  fprintf( stdout, "%10d      BACs referenced in input files\n",
           bac_count );
  fprintf( stdout, "%10d     UIDs not found in gkpStore\n",
           missing_count );
  
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

  if( frag.linkHead == 0 )
    return 0;
  
  CreateGateKeeperLinkRecordIterator( gkp_store.lnkStore,
                                      frag.linkHead,
                                      iid,
                                      &iterator );

  while( NextGateKeeperLinkRecordIterator( &iterator, link ) )
  {
    return ((link->frag1 == iid) ? link->frag2 : link->frag1);
  }

  return 0;
}


static int LookupGrandeGateKeeperRecords(StoreSetp s_set, Locale_ID bacuid, 
                                         CDS_CID_t *baciid,
                                         GateKeeperBatchRecord *gkpb,
                                         GateKeeperBatchRecord *gkpb_next,
                                         GateKeeperLocaleRecord *gkpl,
                                         GateKeeperDistanceRecord *gkpd,
                                         GateKeeperSequenceRecord *gkps)
{
  // Lookup this UID in the Grande gkpstore of locales
  PHashValue_AS value;
  if( HASH_SUCCESS !=
      LookupTypeInPHashTable_AS(s_set->g_gkp_store.hashTable, 
                                UID_NAMESPACE_AS,
                                bacuid,
                                AS_IID_LOC,
                                FALSE,
                                stderr,
                                &value ) )
  {
    fprintf( stderr,
             F_UID " is not a locale UID in grande gatekeeper store\n",
             bacuid );
    return 1;
  }
  else
  {
    fprintf( stderr,
             F_UID " is a locale UID with IID " F_IID " in grande gatekeeper store\n",
             bacuid, value.IID );
  }

  *baciid = value.IID;
  // now, given that locale information

  if( getGateKeeperLocaleStore( (s_set->g_gkp_store).locStore, *baciid, gkpl) )
  {
    fprintf( stderr, "Failed to get locale " F_IID " data from gatekeeperstore\n",
             value.IID );
    return 1;
  }

//  fprintf( stderr," Found locale " F_IID " in Grande GateKeeperStore\n",*baciid);

  if( getGateKeeperSequenceStore( (s_set->g_gkp_store).seqStore, gkpl->sequenceID, gkps) )
  {
     fprintf( stderr, "Failed to get sequence " F_IID " data from gatekeeperstore\n",
              gkpl->sequenceID );
     return 1;
  }
//  fprintf( stderr," Found sequence " F_IID " in Grande GateKeeperStore\n",gkpl->sequenceID);

  if( getGateKeeperBatchStore( (s_set->g_gkp_store).batStore, gkpl->birthBatch, gkpb) )
  {
     fprintf( stderr, "Failed to get batch %u data from gatekeeperstore\n",
              gkpl->birthBatch );
     return 1;
  }

//  fprintf( stderr," Found batch %u in Grande GateKeeperStore\n",gkpl->birthBatch);

  if( getGateKeeperBatchStore( (s_set->g_gkp_store).batStore, gkpl->birthBatch+1, gkpb_next) )
  {
     fprintf( stderr, "Failed to get batch %u data from gatekeeperstore\n",
              gkpl->birthBatch+1 );
     return 1;
  }

//  fprintf( stderr," Found batch %u in Grande GateKeeperStore\n",gkpl->birthBatch);

  if( getGateKeeperDistanceStore( (s_set->g_gkp_store).dstStore, gkpl->lengthID, gkpd) )
  {
     fprintf( stderr, "Failed to get distance " F_IID " data from gatekeeperstore\n",
              gkpl->lengthID );
     return 1;
  }

//  fprintf( stderr," Found distance " F_IID " in Grande GateKeeperStore\n",gkpl->lengthID);

  return 0;
}


int WriteGrandeReadsAndEBACs( Globalsp globals, HashTablep r_rht )
{
  IDSet         idset2;
  IDSet         ddset1;
  IDSetp        ddsetp1;
  IDSetp        idsetp1;
  IDSetp        idsetp2;
  PHashValue_AS value;
  HashTablep    dht;
  GenericMesg   gen;
  FragMesg      frag;
  LinkMesg      lkg;
  DistanceMesg  dist;
  char          src[STRING_LENGTH];
  char          seq[AS_READ_MAX_LEN];
  char          qvs[AS_READ_MAX_LEN];
  GateKeeperLinkRecord link;
  int           write_link;
  int           missing_frags = 0;
  int           missing_ebac_mates = 0;
  int           recruited_mates = 0;
  int           recruited_ebac_mates = 0;
  GateKeeperStore this_gkp_store;
  FragStoreHandle this_frg_store;
  FILE          * this_fp;
  
  frag.source = src;
  frag.sequence = seq;
  frag.quality = qvs;
  
  // allocate a temporary distance IID/UID hashtable
  // to avoid writing out the same distance message more than once
  if( (dht = CreateHashTable( 1000, sizeof( IDSet ) )) == NULL )
  {
    fprintf( stderr, "Failed to allocate hashtable\n" );
    return 1;
  }
  
  // loop through the read UIDs to get gatekeeper IIDs
  RewindHashTable( r_rht );
  while( (idsetp1 = (IDSetp) GetNextHashTableRet( r_rht )) != NULL )
  {
    write_link = 0;
    this_gkp_store = globals->s_set->g_gkp_store;
    this_frg_store = globals->s_set->g_frg_store;
    this_fp = globals->fp_int_out;

    // try to get the grande fragment IID in grande
    // NOTE: one could look up all frags in regional gkpstore, but
    //       regional labels EBACs as reads. Grande cannot.
    //       if the fragment is not in grande, use the regional one
    if( LookupTypeInPHashTable_AS( this_gkp_store.hashTable,
                                   UID_NAMESPACE_AS,
                                   idsetp1->uid,
                                   AS_IID_FRAG, 
                                   FALSE,
                                   stderr,
                                   &value ) != HASH_SUCCESS )
    {
      fprintf( stderr,
               "Fragment " F_UID " is not in gatekeeper store.\n",
               idsetp1->uid );
      continue;
    }
    
    idsetp1->iid = value.IID;

    // write the fragment out
    if( PopulateFragment( &frag,
                          this_frg_store,
                          idsetp1->iid,
                          globals->rs ) )
    {
      fprintf( stderr, "Failed to populate fragment " F_IID "\n",
               idsetp1->iid );
      FreeHashTable( dht );
      return 1;
    }

    // if it's an EBAC, need a DST & BAC message before FRG
    if( frag.type == AS_EBAC )
    {
      this_fp = globals->fp_ext_out;

      // get correct fragment elocale
      if( (frag.elocale = GetLocaleUID( this_gkp_store,
                                        frag.ilocale )) == 0 )
      {
        fprintf( stderr,
           "Failed to get BAC UID from IID " F_IID " in grande gatekeeper store\n",
                 frag.ilocale );
        FreeHashTable( dht );
        return 1;
      }

      // look up mate to find distance iid
      if( (idset2.iid = FindMate( this_gkp_store,
                                  idsetp1->iid,
                                  &link )) == 0 )
      {
        fprintf( stderr,
                 "Failed to find EBAC " F_UID " mate in grande gatekeeper store.\n",
                 idsetp1->uid );
        missing_ebac_mates++;
        continue;
      }

      if( (idset2.uid = GetFragmentUID( this_gkp_store, idset2.iid )) == 0 )
      {
        fprintf( stderr,
                 "Failed to get fragment IID " F_IID " from grande fragstore.\n",
                 idset2.iid );
        FreeHashTable( dht );
        return 1;
      }
      
      // if mate is not in r_ht, or it is but it's mate_uid == 0, then
      // need to write distance & bac message upstream
      // of this first fragment in pair
      if( (idsetp2 = (IDSetp) LookupInHashTable( r_rht,
                                                 idset2.uid )) == NULL ||
          idsetp2->mate_uid == 0 )
      {
        // put together distance data
        ddset1.iid = link.distance;
        if( (ddset1.uid = GetDistanceUID( this_gkp_store, ddset1.iid )) == 0 )
        {
          fprintf( stderr, "Failed to get distance UID from IID " F_IID "\n",
                   ddset1.iid );
          FreeHashTable( dht );
          return 1;
        }
        if( (ddsetp1 = (IDSetp) LookupInHashTable( dht,
                                                   ddset1.uid )) == NULL )
        {
          ddsetp1 = &ddset1;
          // add the distance UID/IID to the g_dht hashtable
          // first, get the UID
          if( PopulateDistanceMesg( this_gkp_store, &ddset1, &dist ) )
          {
            fprintf( stderr, "Failed to populate dist from IID " F_IID "\n",
                     ddset1.iid );
            FreeHashTable( dht );
            return 1;
          }
          
          // insert the distance iid/uid in hashtable
          if( InsertInHashTable( dht, ddset1.uid, (char *) &ddset1 ) )
          {
            fprintf( stderr, "Failed to insert distance iid in hashtable\n" );
            return 1;
          }
          
          // write the distance message
          gen.t = MESG_DST;
          gen.m = &dist;
          WriteProtoMesg_AS( globals->fp_int_out, &gen );
        
        }
        // write ebac message
        ConstructAndWriteEBACMesg( globals->fp_ext_out, &frag, ddsetp1->uid );
      }
      gen.t = MESG_FRG;
      gen.m = &frag;
      WriteProtoMesg_AS( globals->fp_ext_out, &gen );
    }
    else
    {
      // otherwise, just write out the read
      gen.t = MESG_FRG;
      gen.m = &frag;
      WriteProtoMesg_AS( this_fp, &gen );
    }
      
    // see if the read has a mate & get the link & distance data
    if( (idset2.iid = FindMate( this_gkp_store,
                                idsetp1->iid,
                                &link )) != 0 )
    {
      if( (idset2.uid = GetFragmentUID( this_gkp_store, idset2.iid )) == 0 )
      {
        fprintf( stderr,
                 "Failed to get fragment IID " F_IID " from grande fragstore.\n",
                 idset2.iid );
        FreeHashTable( dht );
        return 1;
      }
      // if it's already in r_rht, the mate will get get written later
      if( (idsetp2 = (IDSetp) LookupInHashTable( r_rht,
                                                 idset2.uid )) == NULL )
      {
        if( globals->include_mates )
        {
          if( PopulateFragment( &frag,
                                this_frg_store,
                                idset2.iid,
                                globals->rs ) )
          {
            fprintf( stderr, "Failed to populate grande fragment " F_IID "\n",
                     idset2.iid );
            return 1;
          }
          gen.t = MESG_FRG;
          gen.m = &frag;
          if( this_fp == globals->fp_int_out )
          {
            recruited_mates++;
          }
          else
          {
            recruited_ebac_mates++;
          }
          WriteProtoMesg_AS( this_fp, &gen );
          
          idsetp2 = &idset2;
          write_link = 1;
        }
      }
      else
      {
        // populate the mate_uid, so the link will be written when the
        // mate_uid fragment is encountered
        idsetp1->mate_uid = idsetp2->uid;
        if( idsetp2->mate_uid == idsetp1->uid )
          write_link = 1;
      }

      // write a link if command-line directs all mates to be included or
      // both fragments are in component & both have been written
      if( write_link == 1 )
      {
        // lookup the distance iid in dht to see if we need to write it
        ddset1.iid = link.distance;
        if( (ddset1.uid = GetDistanceUID( this_gkp_store, ddset1.iid )) == 0 )
        {
          fprintf( stderr, "Failed to get distance UID from IID " F_IID "\n",
                   ddset1.iid );
          FreeHashTable( dht );
          return 1;
        }
        if( (ddsetp1 = (IDSetp) LookupInHashTable( dht,
                                                   ddset1.uid )) == NULL )
        {
          ddsetp1 = &ddset1;
          // add the distance UID/IID to the g_dht hashtable
          // first, get the UID
          if( PopulateDistanceMesg( this_gkp_store, &ddset1, &dist ) )
          {
            fprintf( stderr, "Failed to populate dist from IID " F_IID "\n",
                     ddset1.iid );
            FreeHashTable( dht );
            return 1;
          }
          
          // insert the distance iid/uid in hashtable
          if( InsertInHashTable( dht, ddset1.uid, (char *) &ddset1 ) )
          {
            fprintf( stderr, "Failed to insert distance iid in hashtable\n" );
            return 1;
          }

          // write the distance message to the internal file
          gen.t = MESG_DST;
          gen.m = &dist;
          WriteProtoMesg_AS( globals->fp_int_out, &gen );
        }

        // write the link
        lkg.action = AS_ADD;
        if( this_fp == globals->fp_ext_out )
          lkg.type = AS_BAC_GUIDE;
        else
          lkg.type = AS_MATE;
        if( frag.type == AS_TRNR )
          lkg.link_orient = AS_OUTTIE;
        else
          lkg.link_orient = AS_INNIE;
        lkg.distance = ddsetp1->uid;
        lkg.entry_time = frag.entry_time;
        lkg.frag1 = idsetp1->uid;
        lkg.frag2 = idsetp2->uid;
        lkg.distance = ddsetp1->uid;
        gen.t = MESG_LKG;
        gen.m = &lkg;
        WriteProtoMesg_AS( this_fp, &gen );
      }
    }
  }

  FreeHashTable( dht );
  fprintf( stdout, "%10d reads missing from grande gatekeeper store.\n",
           missing_frags );
  fprintf( stdout, "%10d EBAC mates missing from grande gatekeeper store.\n",
           missing_ebac_mates );
  fprintf( stdout, "%10d EBAC mates recruited.\n", recruited_ebac_mates );
  if( globals->include_mates )
    fprintf( stdout, "%10d read mates recruited.\n", recruited_mates );
  return 0;
}


int WriteGrandeShreddedData( Globalsp globals, HashTablep r_bht )
{
  HashTablep bht;
  IDSet      idset1;
  IDSetp     idsetp1;
  IDSetp     idsetp2;
  int        bad_bactig_ids = 0;
  int        missing_ubacs = 0;

  // create a grande BAC hashtable
  if( (bht = CreateHashTable( r_bht->num_nodes,
                              sizeof( IDSet ) )) == NULL )
  {
    fprintf( stderr, "Failed to create BAC/sequence hashtable\n" );
    return 1;
  }

  // loop through the bactig? UIDs and write distance, BAC, & all fragments
  RewindHashTable( r_bht );
  while( (idsetp1 = (IDSetp) GetNextHashTableRet( r_bht )) != NULL )
  {
    // look up this UID in grande gatekeeper store to get BAC UID
    idset1.uid =  idsetp1->uid;

    // look up BAC UID in bht to see if we've already written it
    if( (idsetp2 = (IDSetp) LookupInHashTable( bht, idset1.uid )) == NULL )
    {
      // add it
      GateKeeperBatchRecord gkpb,gkpb_next;
      GateKeeperLocaleRecord gkpl;
      GateKeeperDistanceRecord gkpd;
      GateKeeperSequenceRecord gkps;
      DistanceMesg dist;
      GenericMesg p;
      CDS_CID_t baciid;
      int i;
      if( InsertInHashTable( bht, idset1.uid, (char *) &idset1 ) )
      {
        fprintf( stderr, "Failed to insert BAC UID into hashtable\n" );
        return 1;
      }
      if ( LookupGrandeGateKeeperRecords( globals->s_set, idset1.uid, &baciid, 
                                          &gkpb, &gkpb_next, &gkpl, &gkpd, &gkps ) ){
        fprintf( stderr,
               "Failed to find UBAC info for " F_UID " in Grande GateKeeperStore\n",
                 idset1.uid );
        missing_ubacs++;
        continue;
      }

      dist.action = AS_ADD;
      dist.eaccession = gkpd.UID;
      dist.mean = gkpd.mean;
      dist.stddev = gkpd.stddev;
      p.t = MESG_DST;
      p.m = &dist; 
      WriteProtoMesg_AS(globals->fp_ext_out, &p);

      ConstructAndWriteUBACMesg( globals->fp_ext_out, globals->s_set,
                                 gkpd.UID, gkps.UID, gkpl);

      {
        int32 begin_frag;
        int num_frags=0;
        int locale_found=0;
        FragMesg frag;
        char     src[STRING_LENGTH];
        char     seq[AS_READ_MAX_LEN];
        char     qvs[AS_READ_MAX_LEN];
        GateKeeperFragmentRecord gkpf;
        GateKeeperBactigRecord   gkpbt;

        sprintf( src, " " );
        frag.source = src;
        frag.sequence = seq;
        frag.quality = qvs;
        
        begin_frag = gkpb.numFragments;
        for (i=begin_frag;
             !getGateKeeperFragmentStore( globals->s_set->g_gkp_store.frgStore,
                                          i,
                                          &gkpf );
             i++)
        {
          // Get the fragment from the gatekeeper store
          // check whether frag is in the right locale, and if so, fill in
          // missing info and write it out to the EXT stream
          // if not, we're done with this locale, and need to break
          if( gkpf.localeID == baciid )
          {
            // here's a keeper 
            locale_found = 1;

            // populate the fragment from the fragment store
            if( PopulateFragment( &frag, globals->s_set->g_frg_store, i,
                                globals->rs ) )
            {
              fprintf( stderr, "Failed to populate grande fragment %u\n",
                       i );
              return 1;
            }

            if( gkpf.bactigID == 0 )
            {
              if( frag.type != AS_FBAC )
              {
                fprintf( stderr,
                         "Found fragment " F_UID " with 0 bactig ID.\n",
                         gkpf.readUID );
                bad_bactig_ids++;
                continue;
              }
              else
              {
                frag.ebactig_id = 0;
              }
            }
            else
            {
              // need to determine which bactig this is
              getGateKeeperBactigStore( globals->s_set->g_gkp_store.btgStore,
                                        gkpf.bactigID,
                                        &gkpbt );
              frag.ebactig_id = gkpbt.UID;
            }
            
            frag.elocale = gkpl.UID;
            frag.eseq_id = gkps.UID;
            p.t = MESG_FRG;
            p.m = &frag;
            WriteProtoMesg_AS(globals->fp_ext_out, &p);
            num_frags++;
          }
          else
          {
            if (locale_found) break;
          }
        } 
      }
    }
  }

  FreeHashTable( bht );
  fprintf( stdout, "%10d reads with bad bactig IIDs.\n",
           bad_bactig_ids );
  fprintf( stdout, "%10d UBACs missing from grande gatekeeper store.\n",
           missing_ubacs );
  
  return 0;
}


int GenerateFRGFile( Globalsp globals )
{
  HashTablep   r_rht = NULL;  // regional fragment UIDs
  HashTablep   r_bht = NULL;  // regional locale/BAC UIDs

  // get IIDs of reads & bactigs from input files & convert to UIDs
  // by looking up in regional fragstore & bactigstore
  if( GetInputUIDs( globals, &r_rht, &r_bht ) )
  {
    fprintf( stderr, "Failed to read input files.\n" );
    return 1;
  }

  // open & start writing output file
  if( InitializeOutput( globals ) )
  {
    fprintf( stderr, "Failed to initialize output file\n" );
    return 1;
  }

  // write Celera data - DSTs, FRGs, LKGs
  // this includes BAC ends...
  if( WriteGrandeReadsAndEBACs( globals, r_rht ) )
  {
    fprintf( stderr, "Failed to write Celera internal data\n" );
    return 1;
  }

  // write External data - DSTs, BACs, FRGs
  if( WriteGrandeShreddedData( globals, r_bht ) )
  {
    fprintf( stderr, "Failed to write external data\n" );
    return 1;
  }
  
  FreeHashTable( r_rht );
  FreeHashTable( r_bht );
  fclose( globals->fp_int_out );
  fclose( globals->fp_ext_out );
  return 0;
}


// Usage message followed by abortive exit
void Usage( char * program_name, char * message )
{
  if( message != NULL )
    fprintf( stderr, "%s: %s\n", program_name, message );
  
  fprintf( stderr, "Usage: %s\n", program_name );
  /*
  fprintf( stderr, "\t[-l file listing .cns filenames]\n" );
  fprintf( stderr, "\t[-f regional fragment/gatekeeper store]\n" );
  fprintf( stderr, "\t[-b regional bactig/gatekeeper store]\n" );
  */
  fprintf( stderr, "\t[-G grande gatekeeper store]\n" );
  fprintf( stderr, "\t[-F grande fragment store]\n" );
  fprintf( stderr, "\t[-m include unrecruited mates]\n" );
  fprintf( stderr, "\t<UID list files>\n" );
  
  fprintf( stderr,
           "\nOutput will be written to uniquely named .frg files\n" );

  fprintf( stderr,
           "\nUIDs may be mixture of fragments and BACs\n" );

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
  while( !errflg && ((ch = getopt( argc, argv, "G:F:m" )) != EOF) )
  {
    switch( ch )
    {
      /*
      case 'l':
        globals->list_file_name = optarg;
        fprintf( stdout, "* File listing filenames: %s\n",
                 globals->list_file_name );
        break;
      case 'f':
        globals->r_frg_store_name = optarg;
        fprintf( stdout, "* Regional fragment store: %s\n",
                 globals->r_frg_store_name );
        break;
      case 'b':
        globals->r_bac_store_name = optarg;
        fprintf( stdout, "* Regional bactig store: %s\n",
                 globals->r_bac_store_name );
        break;
        */
      case 'G':
        globals->g_gkp_store_name = optarg;
        fprintf( stdout, "* Grande gatekeeper store: %s\n",
                 globals->g_gkp_store_name );
        break;
      case 'F':
        globals->g_frg_store_name = optarg;
        fprintf( stdout, "* Grande fragment store: %s\n",
                 globals->g_frg_store_name );
        break;
      case 'm':
        globals->include_mates = 1;
        fprintf( stdout, "* Including unrecruited mates\n" );
        break;
      default:
        Usage( globals->program_name, "Unknown parameter or flag." );
        break;
    }
  }

  // check that all parameters were supplied
  /*
  if( ! globals->list_file_name )
    Usage( globals->program_name, "Specify a file listing .cns filename." );
  if( ! globals->r_frg_store_name )
    Usage( globals->program_name, "Specify an regional fragment store." );
  if( ! globals->r_bac_store_name )
    Usage( globals->program_name, "Specify an regional bactig store." );
  */
  if( ! globals->g_gkp_store_name )
    Usage( globals->program_name, "Specify a grande gatekeeper store." );
  if( ! globals->g_frg_store_name )
    Usage( globals->program_name, "Specify a grande fragment store." );
  
  // read the list of UBACs
  if( (globals->f_set = ReadListOfFiles( argc, argv, optind )) == NULL )
  {
    fprintf( stderr, "Failed to open UID list files\n" );
    exit( 1 );
  }
}


int main( int argc, char ** argv )
{
  Globals   globals;
  int       i;

  // parse parameters & exit on any problem
  ProcessCommandLine( &globals, argc, argv );
  
  fprintf( stdout, "Getting data from %d uid list files:\n",
           globals.f_set->num_strings );
  for( i = 0; i < globals.f_set->num_strings; i++ )
    fprintf( stdout, "%s\n", globals.f_set->strings[i] );
  
  // open the stores
  if( (globals.s_set = OpenStores( globals.g_gkp_store_name,
                                   globals.g_frg_store_name )) == NULL )
  {
    fprintf( stderr, "Failed to open one or more stores.\n" );
    return 1;
  }
  
  // allocate a reusable readstruct;
  if( (globals.rs = new_ReadStruct()) == NULL )
  {
    fprintf( stderr, "Failed to allocate readStruct.\n" );
    return 1;
  }

  if( (globals.uid_fetcher = CreateUIDInteractor( 1 )) == NULL )
  {
    fprintf( stderr, "Failed to connect to UID generator.\n" );
    return 1;
  }
  
  // read & write
  if( GenerateFRGFile( &globals ) )
  {
    fprintf( stderr, "Failed to create FRG file from inputs.\n" );
    return 1;
  }

  DestroyUIDInteractor( globals.uid_fetcher );
  delete_ReadStruct( globals.rs );
  CloseStores( globals.s_set );
  FreeStringSet( globals.f_set );

  return 0;
}
