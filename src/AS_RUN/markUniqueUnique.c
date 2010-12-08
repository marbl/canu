const char *mainid = "$Id: markUniqueUnique.c,v 1.9 2010-12-08 12:43:28 skoren Exp $";

//  Assembly terminator module. It is the backend of the assembly
//  pipeline and replaces internal accession numbers by external
//  accession numbers.

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <unistd.h>
#include  <assert.h>

#include "AS_global.h"
#include "AS_UTL_UID.h"
#include "AS_UTL_Hash.h"
#include "AS_UTL_Var.h"
#include "AS_MSG_pmesg.h"
#include "MultiAlignStore.h"


#define DEFAULT_UNITIG_LENGTH           2000
#define DEFAULT_NUM_INSTANCES              1
#define DEFAULT_DISTANCE_TO_ENDS        1000
#define NUM_INSTANCES_AT_SCAFFOLD_ENDS     2
#define TIG_STORE_UTGS                     2
VA_DEF(uint32)

int main (int argc, char *argv[]) {
   int      firstFileArg   = 0;
   char    *asmFileName    = NULL;
   int      minLength      = DEFAULT_UNITIG_LENGTH;
   int      numInstances   = DEFAULT_NUM_INSTANCES;
   int      distanceToEnds = DEFAULT_DISTANCE_TO_ENDS;
   uint32   numToggled     = 0;
   uint32   i              = 0;

   argc = AS_configure(argc, argv);
  
   int arg=1;
   int err=0;
   while (arg < argc) {
      if (strcmp(argv[arg], "-a") == 0) {
         asmFileName = argv[++arg];
      } else if (strcmp(argv[arg], "-l") == 0) {
         minLength = atoi(argv[++arg]);
         if (minLength <= 0) err++;
      } else if (strcmp(argv[arg], "-n") == 0) {
         numInstances = atoi(argv[++arg]);
         if (numInstances < 0) err++;
      } else if (strcmp(argv[arg], "-d") == 0) {
         distanceToEnds = atoi(argv[++arg]);
         if (distanceToEnds <= 0) err++;
      } else if ((argv[arg][0] != '-') && (firstFileArg == 0)) {
         firstFileArg = arg;
         arg = argc;
      } else if (strcmp(argv[arg], "-h") == 0) {
         err++;
      } else {
         fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
         err++;
      }
      arg++;
   }

   if ((asmFileName == NULL) || firstFileArg == 0 || (err)) {
      fprintf(stderr, "usage: %s -a asmFile [-l minLength] [-n numInstances] [-d distanceToEnd] CGI file list\n", argv[0]);
      fprintf(stderr, "  -a asmFile       mandatory path to the assembly asm file\n");
      fprintf(stderr, "  CGI fileList     mandatory list of CGI files containing unitigs to be toggled\n");
      fprintf(stderr, "  -l minLength     minimum size of a unitig to be toggled, default=%d)\n", DEFAULT_UNITIG_LENGTH);
      fprintf(stderr, "  -n numInstances  number of instances of a surrogate that is toggled, default = %d\n", DEFAULT_NUM_INSTANCES);
      fprintf(stderr, "  -d distanceToEnd max number of bases the surrogate can be from the end of a scaffold for toggling, default = %d\n", DEFAULT_DISTANCE_TO_ENDS);
      fprintf(stderr, "\n");
      fprintf(stderr, "  Reads assembly,\n");
      fprintf(stderr, "  finds surrogates matching specified parameters.\n");
      fprintf(stderr, "  There are two special cases in addition to the standard length and number of instances parameters.\n");
      fprintf(stderr, "  1. If numInstances is set to 0, all surrogate unitigs with more than one read will be toggled to unique.\n");
      fprintf(stderr, "  2. If a surrogate appears twice, both times at the ends of a scaffold (within distanceToEnd bases) then it is toggled\n");
      fprintf(stderr, "  Scans the input tigStore for unitigs representing the identified surrogates and changes them to be unique in the store.\n");
      exit(1);
   }
  
   HashTable_AS      *UIDtoIID         = CreateScalarHashTable_AS();
   HashTable_AS      *CTGtoFirstUTG    = CreateScalarHashTable_AS();
   HashTable_AS      *CTGtoLastUTG     = CreateScalarHashTable_AS();
   VA_TYPE(int32)    *unitigLength	   = CreateVA_int32(8192);
   VA_TYPE(uint32)   *surrogateCount   = CreateVA_uint32(8192);
   VA_TYPE(uint32)   *surrogateAtScaffoldEnds   = CreateVA_uint32(8192);
   
   GenericMesg    *pmesg;
   FILE           *infp;
   infp = fopen(asmFileName, "r");   
   while ((EOF != ReadProtoMesg_AS(infp, &pmesg))) {
      SnapUnitigMesg    *utg     = NULL;
      SnapConConMesg    *ctg     = NULL;
      SnapScaffoldMesg  *scf     = NULL;
      uint32             count   = 0;
      uint32             forward = TRUE;
      uint32             lastCtg = 0;

      switch(pmesg->t)
      {
         case MESG_UTG:
            utg = (SnapUnitigMesg*)(pmesg->m);
            Setint32(unitigLength, utg->iaccession, &utg->length);

            if (utg->length >= minLength && (utg->status == AS_NOTREZ || utg->status == AS_SEP)) {               
               // store the mapping for this unitig's UID to IID and initialize it's instance counter at 0
               count = 0;
               InsertInHashTable_AS(UIDtoIID, AS_UID_toInteger(utg->eaccession), 0, (uint64)utg->iaccession, 0);               
               Setuint32(surrogateCount, utg->iaccession, &count);
            }
            break;    
         case MESG_CCO:
            ctg = (SnapConConMesg *)(pmesg->m);
            
            for (i = 0; i < (uint32) ctg->num_unitigs; i++) {
               // increment the surrogate unitigs instance counter
               if (ExistsInHashTable_AS(UIDtoIID, AS_UID_toInteger(ctg->unitigs[i].eident), 0)) {
                  uint32 *ret = Getuint32(surrogateCount, (uint32) LookupValueInHashTable_AS(UIDtoIID, AS_UID_toInteger(ctg->unitigs[i].eident), 0));
                  assert(ret != NULL);
                  (*ret)++;

                  // store first surrogate in a contig
                  if (!ExistsInHashTable_AS(CTGtoFirstUTG, AS_UID_toInteger(ctg->eaccession), 0) && 
                        MIN(ctg->unitigs[i].position.bgn, ctg->unitigs[i].position.end) < distanceToEnds) {
                     InsertInHashTable_AS(CTGtoFirstUTG, AS_UID_toInteger(ctg->eaccession), 0, LookupValueInHashTable_AS(UIDtoIID, AS_UID_toInteger(ctg->unitigs[i].eident), 0), 0); 
                  }

                  // also store the last
                  if ((ctg->length - MAX(ctg->unitigs[i].position.bgn, ctg->unitigs[i].position.end)) < distanceToEnds) {
                     ReplaceInHashTable_AS(CTGtoLastUTG, AS_UID_toInteger(ctg->eaccession), 0, LookupValueInHashTable_AS(UIDtoIID, AS_UID_toInteger(ctg->unitigs[i].eident), 0), 0);
                  }
               }
            }
            break;
         case MESG_SCF:
            scf = (SnapScaffoldMesg *)(pmesg->m);
            
            count = scf->iaccession;
            if (scf->contig_pairs[0].orient.isAnti() || scf->contig_pairs[0].orient.isOuttie()) {
               forward = FALSE;
            }
            lastCtg = MAX(scf->num_contig_pairs - 1, 0);
            
            // All four cases below follow the same pattern
            // The first time a surrogate is found at the end of a scaffold, we record the scaffold ID
            // When the surrogate is seen at the end of a second scaffold, we record that it has been found at the ends of two scaffolds (UINT32_MAX)
            // If the surrogate is seen more than once in a single scaffold, it is eliminated (it can't connect two scaffolds)
            // If the surrogate is only seen once at the end of a scaffold (and again in the middle), it is eliminated
            // 1. Contig is first in scaffold and is forward, take the surrogate from the beginning of contig, if it exists                        
            if (ExistsInHashTable_AS(CTGtoFirstUTG, AS_UID_toInteger(scf->contig_pairs[0].econtig1), 0) && forward) {
               uint32 *myval = Getuint32(surrogateAtScaffoldEnds, (uint32) LookupValueInHashTable_AS(CTGtoFirstUTG, AS_UID_toInteger(scf->contig_pairs[0].econtig1), 0));
               if (myval != NULL && (*myval) == scf->iaccession) {
                  count = 0;
               } else if (myval != NULL && (*myval) != 0 && (*myval) != scf->iaccession) {
                  count = UINT32_MAX;
               }
               Setuint32(surrogateAtScaffoldEnds, (uint32) LookupValueInHashTable_AS(CTGtoFirstUTG, AS_UID_toInteger(scf->contig_pairs[0].econtig1), 0), &count);
               count = scf->iaccession;
            }
            // 2. Contig is last in scaffold and is reversed, take the surrogate from the beginning of the contig, if it exists
            if (ExistsInHashTable_AS(CTGtoFirstUTG, AS_UID_toInteger(scf->contig_pairs[lastCtg].econtig2), 0) && !forward) {
               uint32 *myval = Getuint32(surrogateAtScaffoldEnds, (uint32) LookupValueInHashTable_AS(CTGtoFirstUTG, AS_UID_toInteger(scf->contig_pairs[lastCtg].econtig2), 0));
               if (myval != NULL && (*myval) == scf->iaccession) {
                  count = 0;
               } else if (myval != NULL && (*myval) != 0 && (*myval) != scf->iaccession) {
                  count = UINT32_MAX;
               }
               Setuint32(surrogateAtScaffoldEnds, (uint32) LookupValueInHashTable_AS(CTGtoFirstUTG, AS_UID_toInteger(scf->contig_pairs[lastCtg].econtig2), 0), &count);
               count = scf->iaccession;
            }
            // 3. Contig is first in scaffold and is reversed, take the surrogate from the end of the contig, if it exists            
            if (ExistsInHashTable_AS(CTGtoLastUTG, AS_UID_toInteger(scf->contig_pairs[0].econtig1), 0) && !forward) {
               uint32 *myval = Getuint32(surrogateAtScaffoldEnds, (uint32) LookupValueInHashTable_AS(CTGtoLastUTG, AS_UID_toInteger(scf->contig_pairs[0].econtig1), 0));
               if (myval != NULL && (*myval) == scf->iaccession) {
                  count = 0;
               } else if (myval != NULL && (*myval) != 0 && (*myval) != scf->iaccession) {
                  count = UINT32_MAX;
               }
               Setuint32(surrogateAtScaffoldEnds, (uint32) LookupValueInHashTable_AS(CTGtoLastUTG, AS_UID_toInteger(scf->contig_pairs[0].econtig1), 0), &count);
               count = scf->iaccession;
            }
            // 4. Contig is last in scaffold and is forward, take the surrogate from the end of the contig, if it exists
            if (ExistsInHashTable_AS(CTGtoLastUTG, AS_UID_toInteger(scf->contig_pairs[lastCtg].econtig2), 0) && forward) {
               uint32 *myval = Getuint32(surrogateAtScaffoldEnds, (uint32) LookupValueInHashTable_AS(CTGtoLastUTG, AS_UID_toInteger(scf->contig_pairs[lastCtg].econtig2), 0));
               if (myval != NULL && (*myval) == scf->iaccession) {
                  count = 0;
               } else if (myval != NULL && (*myval) != 0 && (*myval) != scf->iaccession) {
                  count = UINT32_MAX;
               }
               Setuint32(surrogateAtScaffoldEnds, (uint32) LookupValueInHashTable_AS(CTGtoLastUTG, AS_UID_toInteger(scf->contig_pairs[lastCtg].econtig2), 0), &count);
               count = scf->iaccession;
            }
            break;
         default:
            break;
      }
   }
   fclose(infp);
  
   int     toggled   = FALSE;
   uint32 *ret       = NULL;
   uint32 *atScfEnd  = NULL;
   // open the tig store for in-place writing (we don't increment the version since CGW always reads a fixed version initially)
   // this also removes any partitioning
   MultiAlignStore *tigStore = new MultiAlignStore(argv[firstFileArg], TIG_STORE_UTGS, 0, 0, TRUE, TRUE);

   for (i = 0; i < tigStore->numUnitigs(); i++) {
      //MultiAlignT  *ma = tigStore->loadMultiAlign(i, TRUE);
      ret = Getuint32(surrogateCount, i);
      atScfEnd = Getuint32(surrogateAtScaffoldEnds, i);
      toggled = FALSE;
                     
      if (ret != NULL && (*ret) == (uint32)numInstances && numInstances != 0) {
         toggled = TRUE;
      } 
      // if we find a surrogate that has two instances and it is at scaffold ends mark toggle it as well
      else if (ret != NULL && (*ret) == NUM_INSTANCES_AT_SCAFFOLD_ENDS && atScfEnd != NULL && (*atScfEnd) == UINT32_MAX) {
         toggled = TRUE;
      }   
      // special case, mark non-singleton unitigs as unique if we are given no instances
      else if (numInstances == 0 && (*Getint32(unitigLength, i)) >= minLength && tigStore->getNumFrags(i, TRUE) > 1) {
         toggled = TRUE;
      }
      
      if (toggled) {
         tigStore->setUnitigFUR(i, AS_FORCED_UNIQUE);
         numToggled++;
      }      
   }
   
   DeleteHashTable_AS(UIDtoIID);
   DeleteHashTable_AS(CTGtoFirstUTG);
   DeleteHashTable_AS(CTGtoLastUTG);
   delete tigStore;
   
   fprintf(stderr, "Toggled %d\n", numToggled);
   
   return 0;
}
