const char *mainid = "$Id: markUniqueUnique.c,v 1.2 2009-06-10 18:05:14 brianwalenz Exp $";

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

#define DEFAULT_UNITIG_LENGTH    2000
#define DEFAULT_NUM_INSTANCES       1

VA_DEF(uint32);

int main (int argc, char *argv[]) {
   int      firstFileArg   = 0;
   char    *asmFileName    = NULL;
   int32    minLength      = DEFAULT_UNITIG_LENGTH;
   int32    numInstances   = DEFAULT_NUM_INSTANCES;
   int32    numToggled     = 0;

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
         if (numInstances <= 0) err++;
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
      fprintf(stderr, "usage: %s -a asmFile [-l minLength] [-n numInstances] CGI file list\n", argv[0]);
      fprintf(stderr, "  -a asmFile       mandatory path to the assembly asm file\n");
      fprintf(stderr, "  CGI fileList     mandatory list of CGI files containing unitigs to be toggled\n");
      fprintf(stderr, "  -l minLength       minimum size of a unitig to be toggled, default=%d)\n", DEFAULT_UNITIG_LENGTH);
      fprintf(stderr, "  -n numInstances  number of instances of a surrogate that is toggled, default = %d\n", DEFAULT_NUM_INSTANCES);
      fprintf(stderr, "\n");
      fprintf(stderr, "  Reads assembly,\n");
      fprintf(stderr, "  finds surrogates matches specified parameters.\n");
      fprintf(stderr, "  Scans the input CGI files for unitigs representing the identified surrogates and changes them to be unique.\n");
      exit(1);
   }
  
   HashTable_AS      *UIDtoIID         = CreateScalarHashTable_AS();
   VA_TYPE(uint32)   *surrogateCount   = CreateVA_uint32(8192);
   int                i = 0;
   
   GenericMesg    *pmesg;
   FILE           *infp;
   infp = fopen(asmFileName, "r");   
   while ((EOF != ReadProtoMesg_AS(infp, &pmesg))) {
      SnapUnitigMesg *utg = NULL;
      SnapConConMesg *ctg = NULL;

      switch(pmesg->t)
      {
         case MESG_UTG:
            utg = (SnapUnitigMesg*)(pmesg->m);
            if (utg->length >= minLength && (utg->status == AS_NOTREZ || utg->status == AS_SEP)) {
               uint32 count = 0;

               // store the mapping for this unitig's UID to IID and initialize it's instance counter at 0
               InsertInHashTable_AS(UIDtoIID, AS_UID_toInteger(utg->eaccession), 0, (uint64)utg->iaccession, 0);               
               Setuint32(surrogateCount, utg->iaccession, &count);
            }
            break;    
         case MESG_CCO:
            ctg = (SnapConConMesg *)(pmesg->m);
            for (i = 0; i < ctg->num_unitigs; i++) {
               // increment the surrogate unitigs instance counter
               if (ExistsInHashTable_AS(UIDtoIID, AS_UID_toInteger(ctg->unitigs[i].eident), 0)) {
                  uint32 *ret = Getuint32(surrogateCount, (uint32) LookupValueInHashTable_AS(UIDtoIID, AS_UID_toInteger(ctg->unitigs[i].eident), 0));
                  assert(ret != NULL);
                  (*ret)++;
               }
            }
            break;
      }
   }
   fclose(infp);
  
   for(i = firstFileArg; i < argc; i++){
      infp = fopen(argv[i], "r");

      while ((EOF != ReadProtoMesg_AS(infp, &pmesg))) {
         IntUnitigMesg *utg = NULL;
         switch(pmesg->t)
         {
            case MESG_IUM:
               utg = (IntUnitigMesg *)(pmesg->m);
               uint32 *ret = Getuint32(surrogateCount, utg->iaccession);
               if (ret != NULL && (*ret) == numInstances) {
                  utg->unique_rept = AS_FORCED_UNIQUE;
                  numToggled++;
               }
               break;
         }
         WriteProtoMesg_AS(stdout, pmesg);
      }
      
      fclose(infp);
   }
   
   DeleteHashTable_AS(UIDtoIID);
   fprintf(stderr, "Toggled %d\n", numToggled);
   
   return 0;
}
