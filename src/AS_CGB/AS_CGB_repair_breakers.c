
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
static char CM_ID[] = "$Id: AS_CGB_repair_breakers.c,v 1.1.1.1 2004-04-14 13:50:02 catmandew Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_ALN_aligners.h"
#include "AS_CGB_all.h"
#include "AS_CGB_breakers.h"
#include "AS_CGB_util.h"
#include "AS_CGB_store.h"
#include "AS_CGB_unitigger_globals.h"
#include "AS_CGB_repair_breakers.h"
#include "AS_FGB_buildFragmentHash.h"

/****************************************************************************/
/* Globals */

static int TIMINGS = TRUE;
static MesgWriter WriteMesg_AS = NULL, ErrorWriter_AS = NULL;

#define INITIAL_SIZE 100000

/****************************************************************************/

/*
  Repair spurs outline
  0. parse command line to identify the following:
     a. spurs and/or chimeras files to process
     b. flag specifying repair vs. mark/remove spur/chimeric fragments/edges
        and hard removal (high coverage) or temporary removal (low coverage)
     c. fgbStore - to be modified
     d. fgb parameters - needed for relabeling edges
        NOTE: intent is to append to existing store (-a -i fgbStore)
     e. verbose level, -h option
  1. read .spr (spur) & .chi (chimera) files into arrays of data structures
  2. read .cgb file & identify relevant fragments involved
  3. if repairing, then compute overlaps to add to fgbStore
  4. read the fgbStore, allocating room for new overlaps
  5. sort the fragments & overlaps to be added & marked/removed - last to first
     in terms of order in fgbStore
  6. loop over fragments in spur/chimera list from last to first
     if repairing (i.e., adding edges), shift fgbStore edges down as we go
       & insert new edges in sorted order.
     for all fragments to be marked/removed, label as we go
  7. check_store
  8. save checkpoint
  9. relabel edges in fgbStore
 10. check_store
 11. save & exit

 Verbose options:
 0. always regurgitate command-line options
 2. possibly list relevant fragments & whether spur or chimera
 3. possibly list overlaps found & not found
 General: provide feedback on where program is in process
*/



  /* repair_cgb_spurs parameters
     FGB parameters:
           "[-p filename] parameters filename.\n"
           "-i <GraphStore> The input GraphStore.\n"
           "-o <GraphStore> The input GraphStore.\n"
           "-a append to fgbStore.\n"
           "[-P] Specify ASCII output for protoIO messages.\n"
           "[-X] Turn on developer mode.\n"
           "[-D <int>] Specify a debugging level.\n"
           "[-v <int>] Specify a verbosity level.\n"
?           "[-t <int>] Specify the number of threads.\n"
           "[-A <int>] Specify an analysis level.\n"
?           "[-E <int>] Set the allowed number of input errors.\n"
?           "[-C <int>] Set the checkpoint level.\n"
?           "[-R <string>] The restart instructions.\n"
?           "[-T] Output transitively removeble overlaps.\n"
           "[-N <max IID of the fragments> ]\n"
?           "[-w <int>]The work limit per candidate edge for de-chording.\n"
?           "[-x <int>] Threshold adjacency list position for double thresholding.\n"
?           "[-y <int>] Threshold adjacency list position for single thresholding.\n"
?           "[-z <int>] Adjacency list maximum.\n"

     Additional repair_cgb_spurs parameters
           [-H chimeras]
           [-S spurs]
           [-F fragStore]
           [-O ovl file]
           [-r iid file]
           [-g <int>] repair approach
             1: mark all
             2: remove all
             3: add overlaps where found & mark the rest
             4: add overlaps where found & remove the rest

     Other:
       need [-h] for help - to display usage
  */

/*
typedef struct
{
  cds_uint32    num_ovls;
  cds_uint32    num_allocated;
  OverlapMesg * ovls;
} OverlapSet;
typedef OverlapSet * OverlapSetp;
*/

typedef struct
{
  IntFragment_ID iid;
  int            suffix;
} FragmentEnd;
VA_DEF(FragmentEnd)

// structure to hold an array of IntFragment_ID's
/*
typedef struct
{
  cds_uint32       num_iids;
  cds_uint32       num_allocated;
  FragmentEnd    * iids;
} IIDSet;
typedef IIDSet * IIDSetp;
*/

// types to indicate if overlap(s) between fragment & chunk worked out
typedef enum
{
  RS_PROBLEM = 0,
  RS_FINAL_OVERLAP,
  RS_INTERMEDIATE_OVERLAP,
  RS_NO_OVERLAP
} RSType;


#ifdef STANDALONE
// Globals structure for this application
typedef struct
{
  char         * program_name;
  char         * params_file;
  char         * Input_Graph_Store;
  char         * Output_Graph_Store;
  int            input_store_flag;
  int            output_store_flag;
  int            append_store_flag;
  int            markup_the_graph;
  int            check_point_level;
  int            compress_the_graph;
  int            walk_depth;
  IntFragment_ID iv_start;
  int            developer_mode_flag;
  int            debug_level;
  int            verbosity_level;
  int            analysis_level;
  int            num_threads;
  int            analysis_flag;
  IntFragment_ID as_cgb_max_frag_iid;
  IntFragment_ID maxfrags;
  IntEdge_ID     maxedges;
  size_t         maxtext;
  int            work_limit_per_candidate_edge;
  int            remove_blizzard_overlaps;
  int            use_all_overlaps_in_reaper_pass;
  int            dvt_double_sided_threshold_fragment_end_degree;
  int            con_double_sided_threshold_fragment_end_degree;
  int            cutoff_fragment_end_degree;
  char         * chimeras_file;
  char         * spurs_file;
  char         * frag_store;
  char        ** cgb_files;
  int            num_cgb_files;
  char         * ovl_file;
  char         * iid_file;
  int            breaker_fix;
  BreakerSet   * chims;
  BreakerSet   * spurs;
} RepairGlobals;
typedef RepairGlobals * RepairGlobalsp;
#endif // STANDALONE


#ifdef STANDALONE
// set default or uninitialized values for important variables
static void InitializeGlobals
(
  RepairGlobalsp rg,
  char * program_name
  )
{
  memset( rg, 0, sizeof( RepairGlobals ) );
  rg->program_name = program_name;
  rg->num_threads = 4;
  rg->work_limit_per_candidate_edge = 1000;
  rg->remove_blizzard_overlaps = FALSE;
  rg->use_all_overlaps_in_reaper_pass = TRUE;
  rg->as_cgb_max_frag_iid = 100000 * CGB_MULTIPLIER;
  rg->dvt_double_sided_threshold_fragment_end_degree = 100;
  rg->con_double_sided_threshold_fragment_end_degree = 100;
  rg->cutoff_fragment_end_degree = 1000000;
  rg->breaker_fix = 3;
  rg->markup_the_graph = TRUE;
  rg->check_point_level = 0;
  rg->compress_the_graph = TRUE;
  rg->walk_depth = 100;
  rg->iv_start = 0;
}
#endif // STANDALONE

#ifdef STANDALONE
static void Usage
(
  char * program_name,
  char * message
  )
{
  if( message )
    fprintf( stderr, "%s - %s\n", program_name, message );

  fprintf( stderr, "Usage: %s [options] cgb_files\n"
         "FGB parameters:\n"
           "\t[-a]            Append to fgbStore.\n"
           "\t[-i <path>]     Input fgbStore.\n"
           "\t[-o <path>]     Output fgbStore.\n"
           "\t[-p <filename>] Parameters filename.\n"
           "\t[-P]            ASCII output for protoIO messages.\n"
           "\t[-X]            Turn on developer mode.\n"
           "\t[-D <int>]      Debugging level.\n"
           "\t[-v <int>]      Verbosity level.\n"
           "\t[-t <int>]      Number of threads - ignored.\n"
           "\t[-A <int>]      Analysis level.\n"
           "\t[-T]            Output transitively removeable overlaps.\n"
           "\t[-N <int>]      Max IID of the fragments.\n"
           "\t[-w <int>]      Work limit per candidate edge for de-chording.\n"
           "\t[-x <int>]      Adjacency list position for double thresholding.\n"
           "\t[-y <int>]      Adjacency list position for single thresholding.\n"
           "\t[-z <int>]      Adjacency list maximum.\n"
         "Additional repair_cgb_spurs parameters:\n"
           "\t[-H <filename>] Chimeras file.\n"
           "\t[-S <filename>] Spurs file.\n"
           "\t[-F <path>]     FragStore.\n"
           "\t[-O <filename>] Overlap file of found overlaps to create.\n"
           "\t[-r <filename>] IID file of 'bad' fragments to create.\n"
           "\t[-g <int>]      Type of breaker repairs to make:\n"
           "\t                  1: mark all (default),\n"
           "\t                  2: remove all,\n"
           "\t                  3: add overlaps where found & mark the rest,\n"
           "\t                  4: add overlaps where found & remove the rest.\n"
         "General:\n"
           "\t[-h]            Display this help message.\n",
           program_name );

  exit( 1 );
}
#endif // STANDALONE


#ifdef STANDALONE
static int ParseCommandLine
(
  RepairGlobalsp rg,
  int argc,
  char ** argv
  )
{
  int ierr;
  
  InitializeGlobals( rg, argv[0] );
  
  /**************** Process Command Line Arguments *********************/
  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && 
           ((ch = getopt(argc, argv, 
                         "p:i:o:aPXD:v:t:"
                         "A:g:"
                         "N:w:x:z:"
                         "H:S:F:O:r:"
                         "Z:"
                         "h"
                         )) != EOF))
    {
      switch(ch)
      {
        /* The required command line options: */
        case 'p':
          // -p <file> : The parameters file.
          rg->params_file = optarg;
          fprintf( stderr, "* Parameters file is <%s>.\n", rg->params_file );
          loadParams( rg->params_file );
          break;
        case 'i':
          // -i <file> : The input fragment overlap graph store.
          rg->Input_Graph_Store = optarg;
          fprintf( stderr, "* Input fragment graph store is <%s>.\n",
                   rg->Input_Graph_Store );
          rg->input_store_flag = (rg->Input_Graph_Store[0] != '\0');
          break;
        case 'o':
          // -o <file> : The output fragment overlap graph store.
          rg->Output_Graph_Store = optarg;
          fprintf( stderr, "* Output fragment graph store is <%s>.\n",
                   rg->Output_Graph_Store );
          rg->output_store_flag = (rg->Output_Graph_Store[0] != '\0');
          break;
        case 'a':
          // -a : Append to the input store.
          rg->append_store_flag = TRUE;
          fprintf( stderr, "* Appending the input store.\n" );
          break;
        case 'P':
          // -P : Any "protoIO" messages will be in ASCII rather than
          // binary.
          WriteMesg_AS = OutputFileType_AS( AS_PROTO_OUTPUT );
          break;
        case 'X':
          // -X : Turn on developer mode.
          rg->developer_mode_flag = TRUE;
          fprintf( stderr, "* -X: Turn on developer mode.\n" );
          break;
        case 'D':
          // -D <int> : Specify a debugging level.
          rg->debug_level = atoi( optarg );
          break;
        case 'v':
          // -v <int> : Specify a verbosity level.
          rg->verbosity_level = atoi( optarg );
          break;
        case 't':
          // -t <int> : Specify the number of threads.
          rg->num_threads = atoi( optarg );
#ifdef _OPENMP
          omp_set_num_threads( rg->num_threads );
#endif
          break;
          
          /* Recommened command line options for developer mode: */
        case 'A':
          // -A <int> : Specify an analysis level.
          rg->analysis_level = atoi( optarg );
          rg->analysis_flag = (rg->analysis_level > 0);
          break;

        case 'T':
          // -T : unused
          break;
        case 'N':
          // -N <int> : max IID of the fragments.
          rg->as_cgb_max_frag_iid = atol( optarg );
          fprintf(stderr,"** as_cgb_max_frag_iid = " F_IID "\n",
                  rg->as_cgb_max_frag_iid);
          break;
        case 'w':
          // -w <int> : The work limit per candidate edge for
          // de-chording.
          rg->work_limit_per_candidate_edge = atoi( optarg );
          fprintf( stderr,
                   "* work_limit_per_candidate_edge set to %d\n",
                   rg->work_limit_per_candidate_edge );
          break;
      case 'x':
        // -x <int> : The threshold adjaceny list position for "Double
        // sided thresholding".
	rg->dvt_double_sided_threshold_fragment_end_degree = atoi(optarg);
	fprintf(stderr,
                "* dvt_double_sided_threshold_fragment_end_degree set to %d\n",
		rg->dvt_double_sided_threshold_fragment_end_degree);
	break;
        // case 'y':
        // -y <int> : unused
	// break;
      case 'z':
        // -z <int> : The threshold adjaceny list position for "Containment degree
        // trimming".
	rg->con_double_sided_threshold_fragment_end_degree = atoi(optarg);
	fprintf(stderr,
                "* con_double_sided_threshold_fragment_end_degree set to %d\n",
		rg->con_double_sided_threshold_fragment_end_degree);
	break;
      case 'Q':
        // -Q <int>: Specify Reaper pass.
        {
          int reaper_pass = atoi(optarg);
          rg->use_all_overlaps_in_reaper_pass = (reaper_pass == 0);
        }
        break;
      case 'Z':
        // -Z <boolean> : Remove the blizzard overlaps
        rg->remove_blizzard_overlaps = atoi( optarg );
        break;
        
      case 'H':
        // -H <filename> : identify chimeras file
        rg->chimeras_file = optarg;
        fprintf( stderr, " * chimeras file is %s\n", rg->chimeras_file );
        break;
      case 'S':
        // -S <filename> : identify spurs file
        rg->spurs_file = optarg;
        fprintf( stderr, " * spurs file is %s\n", rg->spurs_file );
        break;
      case 'F':
        // -F <path> : identify fragment store
        rg->frag_store = optarg;
        fprintf( stderr, " * frag store is %s\n", rg->frag_store );
        break;
      case 'O':
        // -O <filename> : identify overlap file
        rg->ovl_file = optarg;
        fprintf( stderr, " * overlap file is %s\n", rg->ovl_file );
        break;
      case 'r':
        // -r <filename> : identify iid file
        rg->iid_file = optarg;
        fprintf( stderr, " * iid file is %s\n", rg->iid_file );
        break;
      case 'g':
        // -g <int> : specify modifications to be made
        rg->breaker_fix = atoi( optarg );
        switch( rg->breaker_fix )
          {
          case 1:
            fprintf( stderr,
                     " * modifications: marking all\n" );
            break;
          case 2:
            fprintf( stderr,
                     " * modifications: removing all\n" );
            break;
          case 3:
            fprintf( stderr,
                     " * modifications: adding overlaps & marking\n" );
            break;
          case 4:
            fprintf( stderr,
                     " * modifications: adding overlaps & removing\n" );
            break;
          default:
            Usage( argv[0], "ERROR: Invalid -B repair mode\n" );
            break;
          }
        break;
      case 'h':
        // help
        Usage( argv[0], NULL );
        break;
      default :
        fprintf( stderr, "Unrecognized option -%c\n", optopt );
        Usage( argv[0], NULL );
        break;
      }
    }
    
    // get the list of cgb files
    rg->num_cgb_files = argc - optind;
    if( rg->num_cgb_files <= 0 )
      Usage( argv[0], "ERROR: Need one or more .cgb files\n" );
    rg->cgb_files = &(argv[optind]);
                      
    // Need input fgb store
    if( !rg->input_store_flag && !rg->ovl_file )
      Usage( argv[0], "ERROR: Either -i or -O option is required.\n" );

    // Need either -o or -a or ovl file
    if( (rg->output_store_flag + rg->append_store_flag != 1) &&
        !rg->ovl_file )
      Usage( argv[0],
             "ERROR: Either (-a or -o) or -O options is required.\n" );

    // need one or both chimera/spur file
    if( rg->chimeras_file == NULL && rg->spurs_file == NULL )
      Usage( argv[0], "ERROR: Need one or both of chimeras & spurs file\n" );

    // need frag store
    if( rg->frag_store == NULL )
      Usage( argv[0], "ERROR: Need to specify frag store\n" );

    // Set up append/output fgb store
    if( rg->Input_Graph_Store )
    {
      char buffer[CMD_BUFFER_SIZE - 1];
      if( rg->append_store_flag )
      {
        sprintf( buffer, "chmod -R u+w %s", rg->Input_Graph_Store );
        ierr = system( buffer );
        if( ierr )
        {
          fprintf( stderr,
                   "ERROR: Trouble getting write access to the input graph "
                   "store directory.\n" );
          exit( 1 );
        }
      }
      sprintf( buffer, "(cd %s)", rg->Input_Graph_Store );
      ierr = system( buffer );
      if( ierr )
      {
        fprintf( stderr,
                 "ERROR: Trouble validating "
                 "the input graph store directory.\n" );
        exit( 1 );
      }
      if( rg->append_store_flag )
      {
        rg->Output_Graph_Store = rg->Input_Graph_Store;
      }
      
      if( rg->output_store_flag )
      {
        sprintf( buffer, "mkdir %s", rg->Output_Graph_Store );
        ierr = system( buffer );
        if( ierr )
        {
          fprintf( stderr,
                   "ERROR: Trouble creating the output graph store directory.\n"
                   "       Did the directory already exist?\n" );
          exit( 1 );
        }
      }
    }
  }
  return 0;
}
#endif // STANDALONE

/*
static void FreeOverlapSet
(
  VA_TYPE(OverlapMesg) * os
  )
{
  if( os )
    DeleteVA_OverlapMesg(os);
}

static OverlapSetp AllocateOverlapSet
(
  cds_uint32 num_ovls
  )
{
  OverlapSetp os;

  if( (os = (OverlapSetp) calloc( 1, sizeof( OverlapSet ) )) == NULL )
    return NULL;

  os->ovls = (OverlapMesg *) calloc( num_ovls, sizeof( OverlapMesg ) );
  if( os->ovls == NULL )
  {
    FreeOverlapSet( os );
    return NULL;
  }
  os->num_allocated = num_ovls;
  
  return os;
}
*/

/*
static void FreeIIDSet
(
  IIDSetp is
  )
{
  if( is )
  {
    if( is->iids )
      free( is->iids );
    free( is );
  }
}

static IIDSetp AllocateIIDSet
(
  cds_uint32 num_iids
  )
{
  IIDSetp is;

  if( (is = (IIDSetp) calloc( 1, sizeof( IIDSet ) )) == NULL )
    return NULL;

  is->iids = (FragmentEnd *) calloc( num_iids, sizeof( FragmentEnd ) );
  if( is->iids == NULL )
  {
    FreeIIDSet( is );
    return NULL;
  }
  is->num_allocated = num_iids;
  
  return is;
}
*/

static int AddIIDToIIDSet( IntFragment_ID iid, int suffix, VA_TYPE(FragmentEnd) * is )
{
  FragmentEnd fe;
  fe.iid = iid;
  fe.suffix = suffix;

  AppendVA_FragmentEnd(is, &fe);
  /*
  if( is->num_allocated <= is->num_iids )
  {
    is->iids =
      (FragmentEnd *) realloc( is->iids, (++(is->num_allocated)) *
                               sizeof( FragmentEnd ) );
    if( is->iids == NULL )
      return 1;
  }
  is->iids[is->num_iids].suffix = suffix;
  is->iids[is->num_iids].iid = iid;
  is->num_iids++;
  */
  return 0;
}


static void AddOvlToOvls
(
  VA_TYPE(OverlapMesg) * ovls,
  OverlapMesg * ovl
  )
{
  AppendVA_OverlapMesg(ovls, ovl);
/*  
  // allocate or reallocate
  if( ovls_in == NULL )
  {
    ovls_out = (OverlapMesg *) calloc( 1, sizeof( OverlapMesg ) );
    if( ovls_out == NULL )
      return NULL;
  }
  else
  {
    ovls_out = (OverlapMesg *) realloc( ovls_in, (*num_ovls + 1) *
                                        sizeof( OverlapMesg ) );
    if( ovls_out == NULL )
      return NULL;
  }
  
  // copy the overlap's memory into the last element of the ovls_out array
  memcpy( &(ovls_out[*num_ovls]), ovl, sizeof( OverlapMesg ) );
  (*num_ovls)++;
  return ovls_out;
*/
}


static int AddOvlsToOvls
(
  VA_TYPE(OverlapMesg) * ovls_out,
  VA_TYPE(OverlapMesg) * ovls_in
  )
{
  int i;
  for(i = 0; i < GetNumVA_OverlapMesg(ovls_in); i++)
    AppendVA_OverlapMesg(ovls_out, GetVA_OverlapMesg(ovls_in, i));
  /*
  // allocate or reallocate
  if( *ovls_out == NULL )
  {
    *ovls_out = (OverlapMesg *) calloc( num_ovls_in, sizeof( OverlapMesg ) );
    if( *ovls_out == NULL )
      return 1;
  }
  else
  {
    *ovls_out = (OverlapMesg *) realloc( *ovls_out,
                                         (*num_ovls_out + num_ovls_in) *
                                         sizeof( OverlapMesg ) );
    if( *ovls_out == NULL )
      return 1;
  }
  
  // copy the overlap array's memory into the last element of the ovls_out array
  memcpy( &((*ovls_out)[*num_ovls_out]),
          ovls_in,
          num_ovls_in * sizeof( OverlapMesg ) );
  (*num_ovls_out) += num_ovls_in;
  */
  return 0;
}


static int AddOvlsToOvlSet
(
  VA_TYPE(OverlapMesg) * ovls_out,
  VA_TYPE(OverlapMesg) * ovls_in
  )
{
  int i;
  for(i = 0; i < GetNumVA_OverlapMesg(ovls_in); i++)
    AppendVA_OverlapMesg(ovls_out, GetVA_OverlapMesg(ovls_in, i));
  /*
  // see if OverlapSet's memory needs to grow
  if( os->num_allocated - os->num_ovls < num_ovls )
  {
    os->ovls =
      (OverlapMesg *) realloc( os->ovls,
                               (os->num_ovls + num_ovls) *
                               sizeof( OverlapMesg ) );
    if( os->ovls == NULL )
      return 1;
  }
  memcpy( &(os->ovls[os->num_ovls]), ovls, num_ovls * sizeof( OverlapMesg ) );
  os->num_ovls += num_ovls;
  */
  return 0;
}



static int PopulateFragment
(
  IntFragment_ID iid,
  InternalFragMesg * ifg,
  FragStoreHandle fstore
  )
{
  static ReadStructp rs = NULL;
  char temp_seq[AS_READ_MAX_LEN];
  char temp_qvs[AS_READ_MAX_LEN];
  cds_uint32 bgn, end;

  rs = (rs == NULL) ? new_ReadStruct() : rs;
  if( rs == NULL )
  {
    fprintf( stderr, "Failed to get new ReadStruct\n" );
    return 1;
  }
    
  if( getFragStore( fstore, iid, FRAG_S_ALL, rs ) )
  {
    fprintf( stderr, "Failed to get fragment " F_IID " sequence\n", iid );
    return 1;
  }
  
  getSequence_ReadStruct( rs, temp_seq, temp_qvs, AS_READ_MAX_LEN );
  getClearRegion_ReadStruct( rs, &bgn, &end, READSTRUCT_LATEST );
  ifg->clear_rng.bgn = 0;
  ifg->clear_rng.end = end - bgn;

  // allocate/copy sequence
  if( (ifg->sequence = (char *) malloc( ifg->clear_rng.end + 1 *
                                        sizeof( char ) )) == NULL )
  {
    fprintf( stderr, "Failed to allocate frag " F_IID " sequence of " F_COORD " bases\n",
             iid, ifg->clear_rng.end );
    return 1;
  }
  strncpy( (char *) ifg->sequence, &(temp_seq[bgn]), ifg->clear_rng.end );
  ifg->sequence[ifg->clear_rng.end] = '\0';

  // allocate/copy quality
  if( (ifg->quality = (char *) malloc( ifg->clear_rng.end + 1 *
                                        sizeof( char ) )) == NULL )
  {
    fprintf( stderr, "Failed to allocate frag " F_IID " quality of " F_COORD " bases\n",
             iid, ifg->clear_rng.end );
    return 1;
  }
  strncpy( (char *) ifg->quality, &(temp_qvs[bgn]), ifg->clear_rng.end );
  ifg->quality[ifg->clear_rng.end] = '\0';

  ifg->source = NULL;
  ifg->iaccession = iid;
  
  return 0;
}

  
static RSType OverlapFragments
(
  InternalFragMesg * ifg1,
  int s1,
  IntMultiPos      * ium2,
  int s2,
  FragStoreHandle fstore,
  VA_TYPE(OverlapMesg) * ovls
  )
{
  OverlapMesg * ovl = NULL;
  InternalFragMesg ifg2;
  int where;

  // get the fragment sequence & clear range (modified for DP_Compare)
  if( PopulateFragment( ium2->ident, &ifg2, fstore ) )
  {
    fprintf( stderr, "ERROR: Failed to populate fragment " F_IID "\n", ium2->ident );
    return RS_PROBLEM;
  }
  
  // get an overlap, if there is one
  ovl = DP_Compare_AS( ifg1, &ifg2,
                       -ifg2.clear_rng.end, ifg1->clear_rng.end,
                       (s1 == s2 ),
                       STANDARD_ERATE, STANDARD_THRESH, STANDARD_RANGE,
                       AS_FIND_ALIGN, &where );
  free( ifg2.sequence );
  free( ifg2.quality );
  
  // if there isn't one, add the fragment to the fragment_list
  // otherwise, add it to the overlap_list
  // and stop if ifg1 is contained or the overlap is a dovetail
  if( !ovl )
  {
    return RS_NO_OVERLAP;
  }
  else
  {
    // set some missing overlap message fields
    ovl->min_offset = ovl->ahg;
    ovl->max_offset = ovl->ahg;
    ovl->polymorph_ct = 0;
    ovl->delta = (signed char *) "\0";
    
    if( ovl->aifrag == ifg1->iaccession )
      ovl->quality = ((cds_float32) getDeltaLength( ovl->delta )) /
        ((ifg1->clear_rng.end - ifg1->clear_rng.bgn) - ovl->ahg );
    else
      ovl->quality = ((cds_float32) getDeltaLength( ovl->delta )) /
        ((ifg1->clear_rng.end - ifg1->clear_rng.bgn) - ovl->bhg );

    // add the overlap to a set
    AddOvlToOvls( ovls, ovl );

    // if ifg1 is contained by ifg2, it's a final overlap
    // if ifg1 contains ifg2, it's an intermediate overlap
    if( ovl->overlap_type == AS_CONTAINMENT )
    {
      if( ovl->bifrag == ifg1->iaccession )
        return RS_FINAL_OVERLAP;
      else
        return RS_INTERMEDIATE_OVERLAP;
    }

    /* Check the dovetail orientation
       This is necessary because you might expect:
            ifg1: ----------->
            ifg2:       ----------->
       but get:
            ifg1:       ----------->
            ifg2: ----------->
    */
    if( s1 )
    {
      if( s2 )
      {
        // expecting AS_INNIE
        if( ovl->orientation == AS_INNIE )
          return RS_FINAL_OVERLAP;
        else
          return RS_INTERMEDIATE_OVERLAP;
      }
      else
      {
        // expecting AS_NORMAL
        if( (ovl->orientation == AS_NORMAL &&
             ovl->aifrag == ifg1->iaccession) ||
            (ovl->orientation == AS_ANTI &&
             ovl->bifrag == ifg1->iaccession) )
          return RS_FINAL_OVERLAP;
        else
          return RS_INTERMEDIATE_OVERLAP;
      }
    }
    else
    {
      if( s2 )
      {
        // expecting AS_NORMAL
        if( (ovl->orientation == AS_NORMAL &&
             ovl->bifrag == ifg1->iaccession) ||
            (ovl->orientation == AS_ANTI &&
             ovl->aifrag == ifg1->iaccession) )
          return RS_FINAL_OVERLAP;
        else
          return RS_INTERMEDIATE_OVERLAP;
      }
      else
      {
        // expecting AS_OUTTIE
        if( ovl->orientation == AS_OUTTIE )
          return RS_FINAL_OVERLAP;
        else
          return RS_INTERMEDIATE_OVERLAP;
      }
    }
  }
}


static int GetFragmentChunkOverlaps
(
  InternalFragMesg * ifg1,
  int s1,
  IntUnitigMesg * chunk2,
  int s2,
  FragStoreHandle fstore,
  VA_TYPE(OverlapMesg) * ovls
)
{
  int           i;
  VA_TYPE(OverlapMesg) *   local_ovls = CreateVA_OverlapMesg(10000);
  
  // loop over fragments in chunk2
  // if suffix of chunk2, then loop from last to first
  if( s2 )
  {
    for( i = chunk2->num_frags - 1; i >= 0; i-- )
    {
      // only overlap with non-contained fragments
      if( chunk2->f_list[i].contained == 0 )
      {
        switch( OverlapFragments( ifg1,
                                  s1,
                                  &(chunk2->f_list[i]),
                                  (s2 == (chunk2->f_list[i].position.end >
                                          chunk2->f_list[i].position.bgn)),
                                  fstore,
                                  local_ovls ) )
        {
          case RS_PROBLEM:
            // something bad happened
            fprintf( stderr, "ERROR: Problem overlapping fragments " F_IID " & " F_IID "\n",
                     ifg1->iaccession,chunk2->f_list[i].ident );
            if( local_ovls != NULL )
              DeleteVA_OverlapMesg( local_ovls );
            return 1;
            break;
          case RS_NO_OVERLAP:
            // no overlap this time means no overlaps at all
            if( local_ovls != NULL )
              DeleteVA_OverlapMesg( local_ovls );
            return 0;
            break;
          case RS_FINAL_OVERLAP:
            // either dovetail or ifg1 is contained - done
            if( AddOvlsToOvls( ovls, local_ovls ) )
            {
              fprintf( stderr, "ERROR: Failed to combine overlap arrays\n" );
              return 1;
            }
            if( local_ovls != NULL )
              DeleteVA_OverlapMesg( local_ovls );
            return 0;
            break;
          case RS_INTERMEDIATE_OVERLAP:
          default:
            // ifg1 contains chunk2 fragment, keep looking
            break;
        }
      }
    }
  }
  else
  {
    // !s2, so loop through chunk2 frags from first to last
    for( i = 0; i < chunk2->num_frags; i++ )
    {
      // only overlap with non-contained fragments
      if( chunk2->f_list[i].contained == 0 )
      {
        switch( OverlapFragments( ifg1,
                                  s1,
                                  &(chunk2->f_list[i]),
                                  (s2 == (chunk2->f_list[i].position.end >
                                          chunk2->f_list[i].position.bgn)),
                                  fstore,
                                  local_ovls ) )
        {
          case RS_PROBLEM:
            // something bad happened
            fprintf( stderr, "ERROR: Problem overlapping fragments " F_IID " & " F_IID "\n",
                     ifg1->iaccession,chunk2->f_list[i].ident );
            if( local_ovls != NULL )
              DeleteVA_OverlapMesg( local_ovls );
            return 1;
            break;
          case RS_NO_OVERLAP:
            // no overlap this time means no overlaps at all
            if( local_ovls != NULL )
              DeleteVA_OverlapMesg( local_ovls );
            return 0;
            break;
          case RS_FINAL_OVERLAP:
            // either dovetail or ifg1 is contained - done
            if( AddOvlsToOvls( ovls, local_ovls ) )
            {
              fprintf( stderr, "ERROR: Failed to combine overlap arrays\n" );
              return 1;
            }
            if( local_ovls != NULL )
              DeleteVA_OverlapMesg( local_ovls );
            return 0;
            break;
          case RS_INTERMEDIATE_OVERLAP:
          default:
            // ifg1 contains chunk2 fragment, keep looking
            break;
        }
      }
    }
  }

  // if here, ran out of fragments to loop through
  if( local_ovls != NULL )
    DeleteVA_OverlapMesg( local_ovls );
  return 0;
}


static int FindSpurAndChimeraOverlaps
(
  BreakerSetp chims,
  BreakerSetp spurs,
  char * frag_store,
  VA_TYPE(OverlapMesg) * osp,
  VA_TYPE(FragmentEnd) * isp
  )
{
  int              i;
  FragStoreHandle  fstore;
  ReadStructp      rs;

  // allocate an initial overlap set
  /*
  *osp =
    AllocateOverlapSet( ((chims != NULL) ? chims->num_breakers : 0) *
                        CHIMERA_OVERLAPS +
                        ((spurs != NULL) ? spurs->num_breakers : 0) *
                        SPUR_OVERLAPS );
  if( *osp == NULL )
  {
    fprintf( stderr, "ERROR: Failed to allocate overlap set\n" );
    return 1;
  }
  */

  // allocate an initial iid set
  /*
  *isp = AllocateIIDSet( ((chims != NULL) ? chims->num_breakers : 0) +
                         ((spurs != NULL) ? spurs->num_breakers : 0) );
  if( *isp == NULL )
  {
    fprintf( stderr, "ERROR: Failed to allocate IID set\n" );
    DeleteVA_OverlapMesg( osp );
    return 1;
  }
  */

  // open the fragment store to get ifg data
  fstore = openFragStore( frag_store, "r" );
  if( fstore == NULLSTOREHANDLE )
  {
    fprintf( stderr, "Failed to open frag store %s\n", frag_store );
    DeleteVA_OverlapMesg( osp );
    DeleteVA_FragmentEnd( isp );
    return 1;
  }
  if( (rs = new_ReadStruct()) == NULL )
  {
    fprintf( stderr, "Failed to get new ReadStruct\n" );
    closeFragStore( fstore );
    DeleteVA_OverlapMesg( osp );
    DeleteVA_FragmentEnd( isp );
    return 1;
  }
    
  // loop over chimeras
  for( i = 0; chims != NULL && i < chims->num_breakers; i++ )
  {
    VA_TYPE(OverlapMesg) * sd_ovls = CreateVA_OverlapMesg(1000);
    VA_TYPE(OverlapMesg) * sa_ovls = CreateVA_OverlapMesg(1000);
    InternalFragMesg ifg;

    // get s fragment
    if( chims->breakers[i].chunks[Chunk_s].num_frags == 0 )
      continue;
    
    if( PopulateFragment( chims->breakers[i].chunks[Chunk_s].f_list[0].ident,
                          &ifg, fstore ) )
    {
      fprintf( stderr,
               "ERROR: Failed to populate chunk 's' fragment (IID=" F_IID ")\n",
               chims->breakers[i].chunks[Chunk_s].f_list[0].ident );
      delete_ReadStruct( rs );
      return 1;
    }
    
    // try to overlap s-d
    // loop over containents (keeping them) until a dovetail is found
    if( GetFragmentChunkOverlaps( &ifg,
//                                  chims->breakers[i].suffixes[Chunk_s],
                                  1,
                                  &(chims->breakers[i].chunks[Chunk_d]),
                                  chims->breakers[i].suffixes[Chunk_d],
                                  fstore,
                                  sd_ovls ) )
    {
      fprintf( stderr, "ERROR: Failed to get fragment/chunk overlaps\n" );
      delete_ReadStruct( rs );
      return 1;
    }
    
    if( !sd_ovls )
    {
      // no overlaps, so add fragment to iid_list to mark/remove
      // & don't bother trying s-a overlap
      if( AddIIDToIIDSet( ifg.iaccession, 1, isp ) )
      {
        fprintf( stderr, "ERROR: Failed to add iid to iid set\n" );
        delete_ReadStruct( rs );
        return 1;
      }
      continue;
    }
        
    // try to overlap s-a
    // loop over containents (keeping them) until a dovetail is found
    if( GetFragmentChunkOverlaps( &ifg,
//                                  chims->breakers[i].suffixes[Chunk_s],
                                  0,
                                  &(chims->breakers[i].chunks[Chunk_a]),
                                  chims->breakers[i].suffixes[Chunk_a],
                                  fstore,
                                  sa_ovls ) )
    {
      fprintf( stderr, "ERROR: Failed to get fragment/chunk overlaps\n" );
      delete_ReadStruct( rs );
      return 1;
    }
    free( ifg.sequence );
    free( ifg.quality );

    // if sa_ovls (implies sd_ovls), add them to ovl set
    if( sa_ovls )
    {
      if( AddOvlsToOvlSet( osp, sd_ovls ) )
      {
        fprintf( stderr, "ERROR: Failed to add overlap to overlap set\n" );
        delete_ReadStruct( rs );
        return 1;
      }
      DeleteVA_OverlapMesg( sd_ovls );
      if( AddOvlsToOvlSet( osp, sa_ovls ) )
      {
        fprintf( stderr, "ERROR: Failed to add overlap to overlap set\n" );
        delete_ReadStruct( rs );
        return 1;
      }
      DeleteVA_OverlapMesg( sa_ovls );
    }
    else
    {
      // no overlaps, so add fragment to iid_list to mark/remove
      // & don't bother trying s-a overlap
      if( AddIIDToIIDSet( ifg.iaccession, 0, isp ) )
      {
        fprintf( stderr, "ERROR: Failed to add iid to iid set\n" );
        delete_ReadStruct( rs );
        return 1;
      }
      DeleteVA_OverlapMesg( sd_ovls );
      continue;
    }
  }

  // loop over spurs
  for( i = 0; spurs && i < spurs->num_breakers; i++ )
  {
    VA_TYPE(OverlapMesg) * sd_ovls = CreateVA_OverlapMesg(1000);
    InternalFragMesg ifg;

    // get s fragment
    if( spurs->breakers[i].chunks[Chunk_s].num_frags == 0 )
      continue;
    if( PopulateFragment( spurs->breakers[i].chunks[Chunk_s].f_list[0].ident,
                          &ifg, fstore ) )
    {
      fprintf( stderr,
               "ERROR: Failed to populate chunk 's' fragment (IID=" F_IID ")\n",
               spurs->breakers[i].chunks[Chunk_s].f_list[0].ident );
      delete_ReadStruct( rs );
      return 1;
    }
    
    // try to overlap s-d
    // loop over containents (keeping them) until a dovetail is found
    // negate the suffix of the s chunk, since wrt D, the other end is involved
    if( GetFragmentChunkOverlaps( &ifg,
                                  1 - spurs->breakers[i].suffixes[Chunk_s],
                                  &(spurs->breakers[i].chunks[Chunk_d]),
                                  spurs->breakers[i].suffixes[Chunk_d],
                                  fstore,
                                  sd_ovls ) )
    {
      fprintf( stderr, "ERROR: Failed to get fragment/chunk overlaps\n" );
      delete_ReadStruct( rs );
      return 1;
    }
    free( ifg.sequence );
    free( ifg.quality );
    
    // if sd_ovls, add them to ovl set
    if( sd_ovls )
    {
      if( AddOvlsToOvlSet( osp, sd_ovls ) )
      {
        fprintf( stderr, "ERROR: Failed to add overlap to overlap set\n" );
        delete_ReadStruct( rs );
        return 1;
      }
      DeleteVA_OverlapMesg( sd_ovls );
    }
    else
    {
      // no overlaps, so add fragment to iid_list to mark/remove
      // & don't bother trying s-a overlap
      if( AddIIDToIIDSet( ifg.iaccession,
                          1 - spurs->breakers[i].suffixes[Chunk_s],
                          isp ) )
      {
        fprintf( stderr, "ERROR: Failed to add iid to iid set\n" );
        delete_ReadStruct( rs );
        return 1;
      }
      continue;
    }
  }
  delete_ReadStruct( rs );
  return 0;
}


static int CreateFragmentSet(
  BreakerSetp chims,
  BreakerSetp spurs,
  VA_TYPE(FragmentEnd) * isp
  )
{
  int i;
  
  // allocate an iid set
  /*
  *isp = AllocateIIDSet( ((chims != NULL) ? chims->num_breakers : 0) +
                         ((spurs != NULL) ? spurs->num_breakers : 0) );
  if( *isp == NULL )
  {
    fprintf( stderr, "ERROR: Failed to allocate IID set\n" );
    return 1;
  }
  */

  // loop over chimeras & add fragments
  for( i = 0; chims && i < chims->num_breakers; i++)
  {
    FragmentEnd fe;
    fe.iid = chims->breakers[i].chunks[Chunk_s].f_list[0].ident;
    fe.suffix = chims->breakers[i].suffixes[Chunk_s];
    AppendVA_FragmentEnd(isp, &fe);
  }

  // loop over spurs & add fragments
  for( i = 0; spurs && i < spurs->num_breakers; i++)
  {
    FragmentEnd fe;
    fe.iid = spurs->breakers[i].chunks[Chunk_s].f_list[0].ident;
    fe.suffix = 1 - spurs->breakers[i].suffixes[Chunk_s];
    AppendVA_FragmentEnd(isp, &fe);
  }

  return 0;
}

static int LabelFragments
(
  THeapGlobals * heapva,
  VA_TYPE(FragmentEnd) * is,
  FragmentHashObject * afr_to_avx,
  int breaker_fix
)
{
  Tlab frag_label;
  Tnes edge_label;
  IntFragment_ID p_brk;

  if( is == NULL || GetNumVA_FragmentEnd(is) == 0 )
    return 0;
  
  // set the labels
  switch( breaker_fix )
  {
    case 1:
    case 3:
      frag_label = AS_CGB_MARKED_BREAKER_FRAG;
      edge_label = AS_CGB_MARKED_BY_BREAKER;
      break;
    case 2:
    case 4:
      frag_label = AS_CGB_REMOVED_BREAKER_FRAG;
      edge_label = AS_CGB_REMOVED_BY_BREAKER;
      break;
    default:
      fprintf( stderr, "ERROR: Invalid modification rule value\n" );
      return 1;
      break;
  }

  // label the fragments appropriately
  // p_* for proximal end, d_* for distal end
  for( p_brk = 0; p_brk < GetNumVA_FragmentEnd(is); p_brk++ )
  {
    FragmentEnd * fe = GetVA_FragmentEnd(is, p_brk);
    const IntFragment_ID p_id = get_vid_FragmentHash(afr_to_avx,fe->iid);
    int is0;
    
    // label the fragment
    set_lab_fragment( heapva->frags, p_id, frag_label );
    
    // label the fragment's edges & label the other fragment's edge
    for( is0 = 0; is0 < 2; is0++ )
    {
      const IntEdge_ID p_start = get_segstart_vertex( heapva->frags,
                                                      p_id, is0 );
      const IntEdge_ID p_finish = get_seglen_vertex( heapva->frags,
                                                     p_id, is0 );
      IntEdge_ID p_ei;
      
      // loop over all edges of the fragment & label them
      for( p_ei = p_start; p_ei < p_finish; p_ei++ )
      {
        const IntFragment_ID a_id = get_avx_edge( heapva->edges, p_ei );
        const IntFragment_ID b_id = get_bvx_edge( heapva->edges, p_ei );
        const IntFragment_ID d_id = (a_id == p_id) ? b_id : a_id;
        const int js0 = (a_id == p_id) ?
          get_bsx_edge( heapva->frags, p_ei ) :
          get_asx_edge( heapva->frags, p_ei );
        const IntEdge_ID d_start = get_segstart_vertex( heapva->frags,
                                                        d_id, js0 );
        const IntEdge_ID d_finish = get_seglen_vertex( heapva->frags,
                                                       d_id, js0 );
        IntEdge_ID d_ei;
        
        set_nes_edge( heapva->edges, p_ei, edge_label );
#if 0
        fix_overlap_edge_mate(frags,edges,p_ei);
#else        
        // now label the reverse direction edge
        // find the edge connecting back to the proximal fragment
        for( d_ei = d_start; d_ei < d_finish; d_ei++ )
        {
          if( get_avx_edge( heapva->edges, d_ei ) == p_id ||
              get_bvx_edge( heapva->edges, d_ei ) == p_id )
          {
            set_nes_edge( heapva->edges, d_ei, edge_label );
            break;
          }
        }
#endif        
      }
    }
  }
  return 0;
}


static int write_overlap_file
( UnitiggerGlobals * rg,
  VA_TYPE(OverlapMesg) * os )
{
  FILE * fp;
  GenericMesg gen;
  size_t i;
  
  if( (fp = fopen( rg->ovl_file, "w" )) == NULL )
  {
    fprintf( stderr, "Failed to open %s for writing\n", rg->ovl_file );
    return 1;
  }

  // write an audit message
  {
    AuditMesg   adt;
    AuditLine   adl;

    adl.next = NULL;
    adl.name = rg->program_name;
    adl.complete = time(0);
    adl.version = "$Revision: 1.1.1.1 $";
    adl.comment = "";
    adt.list = &adl;
    gen.m = &adt;
    gen.t = MESG_ADT;
    if( WriteProtoMesg_AS( fp, &gen ) )
    {
      fprintf( stderr, "ERROR: Call to WriteProtoMesg_AS failed.\n" );
      return 1;
    }
  }

  // loop to write out overlap messages
  gen.t = MESG_OVL;
  for( i = 0; i < GetNumVA_OverlapMesg(os); i++ )
  {
    OverlapMesg * ovl = GetVA_OverlapMesg(os, i);
    gen.m = ovl;
    // check to prevent writing out garbage overlaps
    // NOTE: needs further debugging
    if((ovl->overlap_type == AS_DOVETAIL ||
        ovl->overlap_type == AS_CONTAINMENT) &&
       (ovl->orientation == AS_NORMAL ||
        ovl->orientation == AS_INNIE ||
        ovl->orientation == AS_OUTTIE ||
        ovl->orientation == AS_ANTI) &&
       ((ovl->quality >= 0.0f) &&
        (ovl->quality <  1.0f) )
       ) {
      if( WriteProtoMesg_AS( fp, &gen ) ) 
        {
          fprintf( stderr, "ERROR: Call to WriteProtoMesg_AS failed.\n" );
          return 1;
        }
    }
    else
    {
      fprintf(stderr, "Warning: Encountered problem overlap:\n");
      WriteProtoMesg_AS(stderr, &gen);
    }
  }
    
  fclose( fp );
  return 0;
}


static int write_iid_file( UnitiggerGlobals * rg, VA_TYPE(FragmentEnd) * is )
{
  FILE * fp;
  size_t i;
  
  if( (fp = fopen( rg->iid_file, "w" )) == NULL )
  {
    fprintf( stderr, "Failed to open file %s for writing\n",
             rg->iid_file );
    return 1;
  }

  for( i = 0; i < GetNumVA_FragmentEnd(is); i++ )
  {
    FragmentEnd * fe = GetVA_FragmentEnd(is, i);
    fprintf( fp, "%10 " F_IIDP " %d\n", fe->iid, fe->suffix );
  }

  fclose( fp );
  
  return 0;
}


int main_repair_globals
(
 int argc,
 char * argv[],
 TStateGlobals * gstate,
 THeapGlobals  * heapva,
 UnitiggerGlobals * rg
 )
{
  VA_TYPE(FragmentEnd) * iid_list = NULL;
  VA_TYPE(OverlapMesg) * overlap_list = CreateVA_OverlapMesg(INITIAL_SIZE);
  BreakerSet   * the_chims = NULL;
  BreakerSet   * the_spurs = NULL;

#if 0
  ParseCommandLine( rg, argc, argv );
#endif  

  if( rg->chimeras_file != NULL )
  {
    fprintf( stderr, "Reading chimeras file %s\n", rg->chimeras_file );
    assert(NULL == the_chims);
    the_chims = ReadBreakersFile( rg->chimeras_file, Chimera );
    if( the_chims == NULL )
    {
      fprintf( stderr, "ERROR: Failed to get chimera data from %s\n",
               rg->chimeras_file );
#if 0
      return 1;
#endif
    }
  }

  if( rg->spurs_file != NULL )
  {
    fprintf( stderr, "Reading spurs file %s\n", rg->spurs_file );
    assert(NULL == the_spurs);
    the_spurs = ReadBreakersFile( rg->spurs_file, Spur );
    if( the_spurs == NULL )
    {
      fprintf( stderr, "ERROR: Failed to get spur data from %s\n",
               rg->spurs_file );
      return 1;
    }
  }

  // read in the chunks & fragment IDs from the .cgb file
  fprintf( stderr, "Reading unitig data from cgb files\n" );
  {
    int ii;
    for(ii=0; ii < rg->num_cgb_files; ii++) {
      fprintf( stderr, ">> %d %s\n", ii, (rg->the_cgb_files)[ii]);
    }
  }
  if( GetUnitigData( the_chims, the_spurs,
                     rg->the_cgb_files, rg->num_cgb_files ) )
  {
    fprintf( stderr, "Failed to get chimera & spur unitig data\n" );
    return 1;
  }

  iid_list =
    CreateVA_FragmentEnd
    (((the_chims != NULL) ? the_chims->num_breakers : 0) * CHIMERA_OVERLAPS +
     ((the_spurs != NULL) ? the_spurs->num_breakers : 0) * SPUR_OVERLAPS );
  
  // if instructed to add overlaps where possible, find them
  if( rg->breaker_fix > 2 )
  {
    fprintf( stderr, "Looking for chimera & spur overlaps\n" );
    if( FindSpurAndChimeraOverlaps
        ( the_chims, the_spurs, rg->frag_store,
          overlap_list, iid_list ) )
    {
      fprintf( stderr, "ERROR: Failed to create overlap & fragment lists\n" );
      return 1;
    }
    
    // increase the number of maxedges for opening the frag store later
    rg->maxedges += GetNumVA_OverlapMesg(overlap_list);
  }
  else
  {
    fprintf( stderr, "Creating list of bad chimera & spur fragments\n" );
    // otherwise just create list of fragments to mark/remove
    if( CreateFragmentSet( the_chims, the_spurs, iid_list ) )
    {
      fprintf( stderr, "ERROR: Failed to creage fragment list\n" );
      return 1;
    }
  }

  // free the breaker sets
  FreeBreakerSet( the_chims ); the_chims = NULL;
  FreeBreakerSet( the_spurs ); the_spurs = NULL;

  // if user specified modification of fgb store
  if( rg->Input_Graph_Store )
  {
    FragmentHashObject * afr_to_avx = NULL;
    
    // read the fgb store
    fprintf( stderr, "Reading fgb store %s\n", rg->Input_Graph_Store );
    read_fgb_store( rg->Input_Graph_Store,
                    gstate, heapva,
                    rg->maxfrags, rg->maxedges, rg->maxtext);
    
    // modify the fgb data
    afr_to_avx = build_FragmentHash( heapva->frags, 0);
    if( afr_to_avx == NULL )
    {
      fprintf( stderr, "ERROR: failed to map iids to indices\n" );
      return 1;
    }
      
    // make the initial changes to the fgb store
    if( rg->breaker_fix > 2 )
    {
      InsertOverlapsIntoGraph
        ( heapva, overlap_list, NULL,
          afr_to_avx,
          rg->use_all_overlaps_in_reaper_pass,
          rg->dvt_double_sided_threshold_fragment_end_degree,
          rg->con_double_sided_threshold_fragment_end_degree,
          rg->intrude_with_non_blessed_overlaps_flag,
          rg->overlap_error_threshold
          );
    }
    
    LabelFragments( heapva, iid_list, afr_to_avx, rg->breaker_fix );
    
    destroy_FragmentHash(afr_to_avx);
    SAFE_FREE(afr_to_avx);

    if( rg->Output_Graph_Store )
    {
      char buffer[CMD_BUFFER_SIZE-1] = {0};
      
      fprintf( stderr, "Creating file to append the new "
               "IBA+ADT(+ADL) messages.\n" );
      sprintf( buffer, "cp %s/fgb.iba %s/fgb.iba_tmp" ,
               rg->Input_Graph_Store, rg->Output_Graph_Store );
      if( system( buffer ) )
      {
        fprintf( stderr, "ERROR: Trouble copying the IBA+ADT file." );
        return 1;
      }
    }

    // process the fgb
    fprintf( stderr, "Rebuilding the fragment overlap graph\n" );
    processing_phase( gstate, heapva,
                      rg->check_point_level,
                      rg->Output_Graph_Store,
                      rg->analysis_flag,
                      rg->remove_blizzard_overlaps,
                      rg->compress_the_graph,
                      rg->walk_depth,
                      rg->cutoff_fragment_end_degree,
                      rg->work_limit_per_candidate_edge,
                      rg->iv_start,
                      rg->markup_the_graph );

    // write the fgb store
    {
      char the_path1[CMD_BUFFER_SIZE-1] = {0};
      char the_path2[CMD_BUFFER_SIZE-1] = {0};
      char buffer[CMD_BUFFER_SIZE-1] = {0};
      
      sprintf(the_path1,"%s/%s",rg->Output_Graph_Store,"fgb.ckp_tmp");
      sprintf(the_path2,"%s/%s",rg->Output_Graph_Store,"fgb.ckp");
      
      write_fgb_store( the_path1 , gstate, heapva );
      sprintf( buffer, "mv -f %s %s", the_path1, the_path2 );
      fprintf( stderr, "%s\n", buffer );
      if( system( buffer ) )
      {
        fprintf( stderr,
                 "ERROR: The temporary checkpoint file could not be moved "
                 "into its final position.\n" );
        exit( 1 );
      }
      system_date();
    }
  }  // if user specified modification of fgb store

  // if user specified write ovl file
  if( rg->ovl_file )
  {
    if( write_overlap_file( rg, overlap_list ) )
    {
      fprintf( stderr, "ERROR: Failed to write overlap file.\n" );
      return 1;
    }
  }
  
  if( rg->iid_file )
  {
    if( write_iid_file( rg, iid_list ) )
    {
      fprintf( stderr, "ERROR: Failed to write iid file.\n" );
      return 1;
    }
  }
  
  DeleteVA_OverlapMesg( overlap_list ); overlap_list = NULL;
  DeleteVA_FragmentEnd( iid_list ); iid_list = NULL;
  fprintf( stderr, "%s done.\n", rg->program_name );
  return 0;
}

#ifdef STANDALONE
int main ( int argc, char ** argv )
{
  RepairGlobals    rg     = {0};
  TStateGlobals    gstate = {0};
  THeapGlobals     heapva = {0};

  ParseCommandLine( &rg, argc, argv );

  main_repair_globals
    (
     argc,
     argv,
     &gstate,
     &heapva,
     &rg
     );
}
#endif
