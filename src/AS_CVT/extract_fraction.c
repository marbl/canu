
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
// Standard header files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <assert.h>

// Project header files
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_ReadStruct.h"
#include "AS_UTL_PHash.h"

#define STRING_LENGTH 1024

cds_uint64 global_uid=100000000;



typedef struct
{
  cds_uint32  iid;
  cds_uint64  uid;
  FragType    type;
  cds_float32 mean;
  cds_float32 stddev;
  cds_uint32  num_mates;
  cds_float32 num_used;
} LibraryStats;

VA_DEF(LibraryStats)

typedef struct
{
  VA_TYPE(LibraryStats) * lib_stats;
  cds_uint32 num_bac_ends;
  cds_uint32 num_reads;
  cds_uint32 num_extrs;
  cds_uint32 num_trnrs;
  cds_uint32 num_frags;
  cds_uint32 num_unmated;
  cds_float32 num_unmated_used;
} GateKeeperStats;


void FreeGateKeeperStats(GateKeeperStats * stats)
{
  if(stats)
  {
    if(stats->lib_stats)
      DeleteVA_LibraryStats(stats->lib_stats);
    free(stats);
  }
}


GateKeeperStats * CollectGateKeeperStats(GateKeeperStore * gkp_store)
{
  GateKeeperStats * stats;
  StoreStat store_stat;
  cds_uint32 iid;
  GateKeeperFragmentRecord gkfr;
  GateKeeperDistanceRecord gkdr;
  GateKeeperLinkRecord gklr;

  stats = calloc(1, sizeof(GateKeeperStats));
  if(stats == NULL)
  {
    fprintf(stderr, "Failed to allocate a stats structure\n");
    return NULL;
  }

  // loop over fragments to count reads vs bac ends
  statsStore( gkp_store->frgStore, &store_stat );
  for( iid = 1; iid <= store_stat.lastElem; iid++ )
  {
    if( getGateKeeperFragmentStore( gkp_store->frgStore, iid, &gkfr) )
    {
      fprintf( stderr,
               "Failed to get fragment IID " F_IID " in gatekeeper store.\n",
               iid );
      continue;
    }
    if(!gkfr.deleted)
    {
      stats->num_frags++;
      switch(gkfr.type)
      {
        case AS_READ:
          stats->num_reads++;
          break;
        case AS_EXTR:
          stats->num_extrs++;
          break;
        case AS_TRNR:
          stats->num_trnrs++;
          break;
        default:
          break;
      }
    }
  }

  // now populate libraries
  statsStore( gkp_store->dstStore, &store_stat );
  stats->lib_stats = CreateVA_LibraryStats(store_stat.lastElem + 1);
  ResetVA_LibraryStats(stats->lib_stats);
  EnableRangeVA_LibraryStats(stats->lib_stats, store_stat.lastElem + 1);
  if(stats->lib_stats == NULL)
  {
    fprintf(stderr,
            "Failed to allocate variable array of " F_S64 " distances\n",
            store_stat.lastElem);
    FreeGateKeeperStats(stats);
    return NULL;
  }

  for( iid = 1; iid <= store_stat.lastElem; iid++ )
  {
    LibraryStats * ls;
    if( getGateKeeperDistanceStore( gkp_store->dstStore, iid, &gkdr) )
    {
      fprintf( stderr,
               "Failed to get distance IID " F_IID " in gatekeeper store.\n",
               iid );
      continue;
    }
    if(!gkdr.deleted)
    {
      ls = (GetVA_LibraryStats(stats->lib_stats,iid));
      if(ls == NULL)
      {
        fprintf(stderr,
               "Something went wrong accessing distance iid " F_IID " in var array\n",
                iid);
        FreeGateKeeperStats(stats);
        return NULL;
      }
      ls->iid = iid;
      ls->uid = gkdr.UID;
      ls->mean = gkdr.mean;
      ls->stddev = gkdr.stddev;
    }
  }

  // now populate mate numbers
  stats->num_unmated = stats->num_frags;
  statsStore( gkp_store->lnkStore, &store_stat );
  for(iid = 1; iid <= store_stat.lastElem; iid++)
  {
    if( getGateKeeperLinkStore( gkp_store->lnkStore, iid, &gklr) )
    {
      fprintf( stderr,
               "Failed to get link " F_IID " in gatekeeper store.\n",
               iid );
      continue;
    }

    if(!gklr.deleted && (gklr.type == AS_MATE || gklr.type == AS_BAC_GUIDE))
    {
      LibraryStats * ls;
      ls = (GetVA_LibraryStats(stats->lib_stats,gklr.distance));
      if(ls != NULL)
      {
        ls->num_mates += 2;
        stats->num_unmated -= 2;
      }
    }
  }

  return(stats);
}


int PopulateFragment( FragMesg * f,
                      FragStoreHandle fs,
                      cds_uint32 iid,
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
  getEntryTime_ReadStruct (rs, &(f->entry_time) );
  getClearRegion_ReadStruct( rs, &bgn, &end, READSTRUCT_LATEST );
  f->clear_rng.bgn = bgn;
  f->clear_rng.end = end;
  getSource_ReadStruct( rs, f->source, STRING_LENGTH );
  getSequence_ReadStruct( rs, f->sequence, f->quality, AS_READ_MAX_LEN );
  getReadIndex_ReadStruct( rs, &(f->iaccession) );
  
  return 0;
}

int PopulateFile( char * prog_name,
                  char * out_name,
                  GateKeeperStore * gkp_store,
                  FragStoreHandle frg_store,
                  float fraction,
                  GateKeeperStats * stats)
{
  FILE * fp;
  GenericMesg gen;

  fp = fopen(out_name, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "Failed to open file %s for writing.\n", out_name);
    return 1;
  }

  fprintf(stdout, "Output file %s being written.\n", out_name);

  // write BAT/ADT messages
  {
    BatchMesg bat;
    AuditMesg adt;
    AuditLine adl;
    char comment[STRING_LENGTH];

    bat.name = prog_name;
    bat.created = time(0);
    bat.eaccession = global_uid++;
    sprintf(comment, "%f fraction of %s\n", fraction, gkp_store->storePath);
    bat.comment = comment;

    gen.m = &bat;
    gen.t = MESG_BAT;
    WriteProtoMesg_AS(fp, &gen);

    adt.list = &adl;
    adl.next = NULL;
    adl.name = prog_name;
    adl.complete = time(0);
    adl.version = "$Revision: 1.4 $";
    adl.comment = "Contact Ian @x3036 to report problems";

    gen.m = &adt;
    gen.t = MESG_ADT;
    WriteProtoMesg_AS(fp, &gen);
  }

  // write distance messages
  {
    cds_uint32 iid;
    DistanceMesg dst;
    LibraryStats * ls;

    for(iid = 1; iid < GetNumVA_LibraryStats(stats->lib_stats); iid++)
    {
      ls = (GetVA_LibraryStats(stats->lib_stats,iid));
      if(ls == NULL)
      {
        fprintf(stderr,
                "Failed to get library " F_IID " from variable array\n",
                iid);
        return 1;
      }
    
      dst.action = AS_ADD;
      dst.eaccession = ls->uid;
      dst.mean = ls->mean;
      dst.stddev = ls->stddev;

      gen.m = &dst;
      gen.t = MESG_DST;
      WriteProtoMesg_AS(fp, &gen);
    }
  }

  // loop over links to write output
  {
    cds_uint32 iid;
    StoreStat store_stat;
    GateKeeperLinkRecord gklr;
    FragMesg frg;
    LinkMesg lkg;
    ReadStructp rs;
    char src[STRING_LENGTH];
    char seq[AS_READ_MAX_LEN];
    char qvs[AS_READ_MAX_LEN];
    VA_TYPE(char) * used;
    char used_val=(char) 1;

    // set up a variable array to mark which fragments have been used
    // needed to determine which to use as unmated fragments at the end
    statsStore( gkp_store->frgStore, &store_stat );
    used = CreateVA_char(store_stat.lastElem + 1);
    ResetVA_char(used);
    EnableRangeVA_char(used, store_stat.lastElem + 1);
    
    sprintf(src, " ");
    frg.source = src;
    frg.sequence = seq;
    frg.quality = qvs;

    if((rs = new_ReadStruct()) == NULL)
    {
      fprintf( stderr, "Failed to allocate readStruct.\n" );
      return 1;
    }
    
    statsStore( gkp_store->lnkStore, &store_stat );
    for(iid = 1; iid <= store_stat.lastElem; iid++)
    {
      if( getGateKeeperLinkStore( gkp_store->lnkStore, iid, &gklr) )
      {
        fprintf( stderr,
                 "Failed to get link " F_IID " in gatekeeper store.\n",
                 iid );
        continue;
      }

      if(!gklr.deleted && (gklr.type == AS_MATE || gklr.type == AS_BAC_GUIDE))
      {
        LibraryStats * ls;
        ls = (GetVA_LibraryStats(stats->lib_stats,gklr.distance));
        if(ls != NULL)
        {
          if(ls->num_used / ls->num_mates < fraction)
          {
            // write fragment 1
            PopulateFragment(&frg, frg_store, gklr.frag1, rs);

            lkg.frag1 = frg.eaccession;
            gen.m = &frg;
            gen.t = MESG_FRG;
            WriteProtoMesg_AS(fp, &gen);
            
            // write fragment 2
            PopulateFragment(&frg, frg_store, gklr.frag2, rs);
            lkg.frag2 = frg.eaccession;
            gen.m = &frg;
            gen.t = MESG_FRG;
            WriteProtoMesg_AS(fp, &gen);

            lkg.action = AS_ADD;
            lkg.type = gklr.type;
            lkg.entry_time = frg.entry_time;
            lkg.link_orient = getLinkOrientation(&gklr);
            lkg.distance = ls->uid;
            gen.m = &lkg;
            gen.t = MESG_LKG;
            WriteProtoMesg_AS(fp, &gen);

            SetVA_char(used, gklr.frag1, &used_val);
            SetVA_char(used, gklr.frag2, &used_val);
            
            ls->num_used += 2;
          }
        }
      }
    }
    for(iid = 1;
        iid < GetNumVA_char(used) &&
          (stats->num_unmated_used / stats->num_unmated < fraction);
        iid++)
    {
      if(*GetVA_char(used, iid) != used_val)
      {
        // use it as a pair of unmated fragments
        // write fragment 1
        PopulateFragment(&frg, frg_store, iid, rs);
        gen.m = &frg;
        gen.t = MESG_FRG;
        WriteProtoMesg_AS(fp, &gen);
        stats->num_unmated_used++;
      }
    }
    delete_ReadStruct(rs);
    DeleteVA_char(used);
  }

  fclose(fp);
  
  return 0;
}

void WriteStats( GateKeeperStats * stats, char * gkp_name)
{
  cds_uint32 iid;
  
  fprintf(stdout, "Gatekeeper store: %s\n", gkp_name );

  fprintf(stdout, "\nFragment types:\n");
  fprintf(stdout, "\tCelera reads:     %10u\n", stats->num_reads);
  fprintf(stdout, "\tExternal reads:   %10u\n", stats->num_extrs);
  fprintf(stdout, "\tTransposon reads: %10u\n", stats->num_trnrs);
  fprintf(stdout, "\tBAC ends:         %10u\n", stats->num_bac_ends);
  fprintf(stdout, "\tTotal fragments===%10u\n", stats->num_frags);

  for(iid = 1; iid < GetNumVA_LibraryStats(stats->lib_stats); iid++)
  {
    LibraryStats * ls;
    ls = (GetVA_LibraryStats(stats->lib_stats,iid));

    fprintf(stdout, "\nLibrary %20" F_UIDP ":\n", ls->uid);
    fprintf(stdout, "\ttype: %c\n", ls->type);
    fprintf(stdout, "\tmean:   %9.2f\n", ls->mean);
    fprintf(stdout, "\tstddev: %9.2f\n", ls->stddev);
    fprintf(stdout, "\tmate pairs: %9u\n", ls->num_mates / 2);
  }

  fprintf(stdout, "\nUnmated fragments: %10u\n", stats->num_unmated);
}


int main( int argc, char ** argv )
{
  cds_uint32 i;
  cds_uint32 first_var=1;
  cds_uint32 num_stores=0;
  GateKeeperStore gkp_store;
  GateKeeperStats * stats;

  fprintf(stderr, "Lightly shotgunned, unfinished, &"
          " finished BAC data is ignored.\n");
  if(argc > 2 && strcmp(argv[1], "-i") == 0)
  {
    fprintf(stdout, "Just getting store stats\n");
    first_var=2;
  }
  else if((argc - 1) % 3 != 0)
  {
    fprintf( stderr,
             "Usage: %s <-i> (<gkp_store> <frg_store> <fraction>)*\n",
             argv[0] );
    fprintf( stderr, "\t-i just dumps store stats\n");
    exit( 1 );
  }

  num_stores = (argc - first_var) / 3;
  if(num_stores == 0 && first_var == 2)
    num_stores = 1;
  
  // loop over stores
  for(i = 0; i < num_stores; i++)
  {
    char * gkp_name = argv[i * 3 + first_var];
    
    // open the gatekeeper store
    InitGateKeeperStore( &gkp_store, gkp_name );
    OpenReadOnlyGateKeeperStore( &gkp_store );

    // collect stats on which to base proportional extraction
    stats = CollectGateKeeperStats(&gkp_store);

    if(first_var==1)
    {
      char * frg_name = argv[i * 3 + first_var + 1];
      float fraction = atof(argv[i * 3 + first_var + 2]);
      FragStoreHandle frg_store;
      char out_name[STRING_LENGTH];
      char * temp_char;
      
      // set up an output .frg file for the store
      temp_char = rindex(frg_name, (int) '/');
      if(temp_char)
        strcpy(out_name, &(temp_char[1]));
      else
        strcpy(out_name, frg_name);
      temp_char = rindex(out_name, (int) '.');
      if(temp_char)
        temp_char[0] = '\0';
      strcat(out_name, ".frg");
      
      // open the fragment store
      frg_store = openFragStore( frg_name, "r" );
      if( frg_store == NULLSTOREHANDLE )
      {
        fprintf( stderr,
                 "Failed to open regional frag store %s\n",
                 frg_name );
        exit( 1 );
      }

      // generate output
      if(PopulateFile(argv[0], out_name,
                      &gkp_store, frg_store, fraction,
                      stats))
      {
        fprintf(stderr, "Failed to write output file %s\n", out_name);
        return 1;
      }

      // close fragstore
      closeFragStore(frg_store);
    }
    else
    {
      // dump stats for library
      WriteStats(stats, gkp_name);
    }

    // free stats
    FreeGateKeeperStats(stats);
    
    // close gatekeeper store
    CloseGateKeeperStore(&gkp_store);
  }
  
  return 0;
}
