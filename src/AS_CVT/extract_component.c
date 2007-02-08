
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
static char CM_ID[] = "$Id: extract_component.c,v 1.8 2007-02-08 06:48:52 brianwalenz Exp $";

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
#include <unistd.h> /* man 3 getopt */
#include <string.h>
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
  GateKeeperStore rb_gkp_store;
  GateKeeperStore rf_gkp_store;
  FragStoreHandle r_frg_store;
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
  char           * list_file_name;      // file listing .cns filenames
  char           * r_frg_store_name;    // regional fragstore name
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
StringSetp ReadListOfFiles( char * filename )
{
  FILE * fp;
  StringSetp f_set = NULL;
  int num_files;

  if( (fp = fopen( filename, "r" )) == NULL )
  {
    fprintf( stderr, "Failed to open file %s for reading.\n", filename );
    return NULL;
  }

  num_files = CountLines( fp );
  if( num_files )
  {
    int i;
    f_set = AllocateStringSet( num_files );
    if( f_set == NULL )
    {
      fprintf( stderr, "Failed to allocate set of %d files\n",
               num_files );
      fclose( fp );
      return NULL;
    }

    for( i = 0; i < f_set->num_strings; i++ )
    {
      if( ReadNextLine( fp, f_set->strings[i] ) )
      {
        fprintf( stderr, "Failed to read file %s\n", f_set->strings[i] );
        FreeStringSet( f_set );
        fclose( fp );
        return NULL;
      }
    }
  }

  fclose( fp );
  return f_set;
}



/****************************************************************************
 * Populate a FragMesg from a fragStore given an iid
 ****************************************************************************/
int PopulateFragment( FragMesg * f, FragStoreHandle fs, cds_uint32 iid,
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
    CloseGateKeeperStore( &(s_set->rb_gkp_store) );
    CloseGateKeeperStore( &(s_set->rf_gkp_store) );
    if( s_set->r_frg_store )
      closeFragStore( s_set->r_frg_store );
    
    CloseGateKeeperStore( &(s_set->g_gkp_store) );
    if( s_set->g_frg_store )
      closeFragStore( s_set->g_frg_store );
    
    free( s_set );
  }
}


StoreSetp OpenStores( char * r_frg_store_name,
                      char * g_gkp_store_name,
                      char * g_frg_store_name )

{
  StoreSetp s_set;

  // allocate the containing structure
  if( (s_set = (StoreSetp) calloc( 1, sizeof( StoreSet ) )) == NULL )
  {
    fprintf( stderr, "Failed to allocate StoreSet\n" );
    return NULL;
  }

  // Open regional stores
  {
    // open the regional gatekeeper stores
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
    
  }

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



int WriteBATADT( Globalsp globals)
{
  GenericMesg  gen;
  BatchMesg    bat;
  AuditLine    adl;
  AuditMesg    adt;
  char       * comment;
  int          i;
  FILE *fp = globals->fp_int_out;

  // first message must be BAT with true UID;
  bat.name = globals->program_name;
  bat.created = time(0);
  bat.eaccession = globals->batch_int_uid;
  bat.comment = "Grande dataset from regional component(s)";
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
  sprintf( comment, "Based on regional assembly files:\n" );
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

  if ( WriteBATADT( globals) ) return 0;
}



Fragment_ID GetFragmentUID( GateKeeperStore gkp_store,
                            IntFragment_ID iid )
{
  GateKeeperFragmentRecord frag;

  getGateKeeperFragmentStore( gkp_store.frgStore, iid, &frag );

  return frag.readUID;
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


int GetInputUIDs( Globalsp globals, HashTablep * rht )
{
  int    fs_i;
  int    read_count;
  IDSet  idset;

  // create read hashtable
  if( (*rht = CreateHashTable( 400000, sizeof( IDSet ) )) == NULL )
  {
    fprintf( stderr, "Failed to create read hashtable\n" );
    return 1;
  }

  // loop over input files
  read_count = 0;
  memset( &idset, 0, sizeof( IDSet ) );
  for( fs_i = 0; fs_i < globals->f_set->num_strings; fs_i++ )
  {
    FILE          * fp_in;
    GenericMesg   * genp;
    IntUnitigMesg * ium;
    int             i;

    fprintf( stdout, "Reading data from %s\n",
             globals->f_set->strings[fs_i] );
    
    // open the input file
    if( (fp_in = fopen( globals->f_set->strings[fs_i], "r" )) == NULL )
    {
      fprintf( stderr, "Failed to open file %s for reading\n",
               globals->f_set->strings[fs_i] );
      FreeHashTable( *rht );
      return 1;
    }
    
    // iterate through input messages
    while( ReadProtoMesg_AS( fp_in, &genp ) != EOF )
    {
      switch( genp->t )
      {
        case MESG_IUM:
          ium = genp->m;
          for( i = 0; i < ium->num_frags; i++ )
          {
            switch( ium->f_list[i].type )
            {
              case AS_READ:
              case AS_EXTR:
              case AS_TRNR:
                // read the fragment store's fixed data for the fragment
                idset.iid = ium->f_list[i].ident;
                if( (idset.uid =
                     GetFragmentUID( globals->s_set->rf_gkp_store,
                                     ium->f_list[i].ident )) == 0 )
                {
                  fprintf( stderr,
                           "Failed to get fragment " F_IID " data from"
                           " regional fragstore\n",
                           ium->f_list[i].ident );
                  fclose( fp_in );
                  FreeHashTable( *rht );
                  return 1;
                }
                idset.type = ium->f_list[i].type;

                // add it if it's not already present
                if( ! LookupInHashTable( *rht, idset.uid ) )
                {
                  if( InsertInHashTable( *rht, idset.uid, (char *) &idset ) )
                  {
                    fprintf( stderr,
                             "Failed to insert UID into read hasthable\n" );
                    fclose( fp_in );
                    FreeHashTable( *rht );
                    return 1;
                  }
                  read_count += (ium->f_list[i].type == AS_READ ||
                                 ium->f_list[i].type == AS_EXTR ||
                                 ium->f_list[i].type == AS_TRNR) ? 1 : 0;
                }
                break;
              default:
                fprintf( stderr,
                         "Unknown IMP type %u\n",
                         ium->f_list[i].type );
                break;
            }
          }
        default:
          break;
      }
    }
    fclose( fp_in );
  }
  
  fprintf( stdout, "%10d     reads referenced in input files\n",
           read_count );

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



int WriteGrandeReads( Globalsp globals, HashTablep r_rht )
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
  int           recruited_mates = 0;
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
    if( LookupTypeInPHashTable_AS( this_gkp_store.hashTable,
                                   UID_NAMESPACE_AS,
                                   idsetp1->uid,
                                   AS_IID_FRG, 
                                   FALSE,
                                   stderr,
                                   &value ) != HASH_SUCCESS )
    {
      fprintf( stdout,
               "Retrieving fragment (" F_UID "," F_IID ") from regional stores.\n",
               idsetp1->uid, idsetp1->iid );
      
      this_gkp_store = globals->s_set->rf_gkp_store;
      this_frg_store = globals->s_set->r_frg_store;

      // NOTE: regional iid is already in idsetp1
    }
    else
    {
      idsetp1->iid = value.IID;
    }

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

    {
      // just write out the read
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
            recruited_mates++;
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
  if( globals->include_mates )
    fprintf( stdout, "%10d read mates recruited.\n", recruited_mates );
  return 0;
}


int GenerateFRGFile( Globalsp globals )
{
  HashTablep   r_rht = NULL;  // regional fragment UIDs

  // get IIDs of reads from input files & convert to UIDs
  // by looking up in regional fragstore
  if( GetInputUIDs( globals, &r_rht ) )
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
  if( WriteGrandeReads( globals, r_rht ) )
  {
    fprintf( stderr, "Failed to write Celera internal data\n" );
    return 1;
  }

  FreeHashTable( r_rht );
  fclose( globals->fp_int_out );
  return 0;
}


// Usage message followed by abortive exit
void Usage( char * program_name, char * message )
{
  if( message != NULL )
    fprintf( stderr, "%s: %s\n", program_name, message );
  
  fprintf( stderr, "Usage: %s\n", program_name );
  fprintf( stderr, "\t[-l file listing .cns filenames]\n" );
  fprintf( stderr, "\t[-f regional fragment/gatekeeper store]\n" );
  fprintf( stderr, "\t[-G grande gatekeeper store]\n" );
  fprintf( stderr, "\t[-F grande fragment store]\n" );
  fprintf( stderr, "\t[-m include unrecruited mates]\n" );
  fprintf( stderr,
           "\nOutput will be written to a uniquely named .frg file\n" );

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
  while( !errflg && ((ch = getopt( argc, argv, "l:f:b:G:F:m" )) != EOF) )
  {
    switch( ch )
    {
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
  if( ! globals->list_file_name )
    Usage( globals->program_name, "Specify a file listing .cns filename." );
  if( ! globals->r_frg_store_name )
    Usage( globals->program_name, "Specify an regional fragment store." );
  if( ! globals->g_gkp_store_name )
    Usage( globals->program_name, "Specify a grande gatekeeper store." );
  if( ! globals->g_frg_store_name )
    Usage( globals->program_name, "Specify a grande fragment store." );

}


int main( int argc, char ** argv )
{
  Globals   globals;
  int       i;

  // parse parameters & exit on any problem
  ProcessCommandLine( &globals, argc, argv );
  
  // read the list of UBACs
  if( (globals.f_set = ReadListOfFiles( globals.list_file_name )) == NULL )
  {
    fprintf( stderr, "Failed to open ubac_list_file %s\n",
             globals.list_file_name );
    return 1;
  }
  fprintf( stdout, "Getting data from %d consensus files:\n",
           globals.f_set->num_strings );
  for( i = 0; i < globals.f_set->num_strings; i++ )
    fprintf( stdout, "%s\n", globals.f_set->strings[i] );
  
  // open the stores
  if( (globals.s_set = OpenStores( globals.r_frg_store_name,
                                   globals.g_gkp_store_name,
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
