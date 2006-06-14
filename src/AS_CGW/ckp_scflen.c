static char CM_ID[] = "$Id:";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "cds.h"
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_timer.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "DiagnosticsCGW.h"
#include "ScaffoldGraph_CGW.h"
#include "AS_UTL_Var.h"

int main (int argc , char * argv[] ) {

  Global_CGW *data;
  char *prefix;
  int setFragStore = FALSE;
  int setGatekeeperStore = FALSE;
  int setPrefixName = FALSE;
  int ckptNum = NULLINDEX;
  int i, index;
  char subset_map[1000];
  char full_map[1000];
  char setSubsetMap=0;
  char setFullMap=0;
  char full_ovlPath[1000];
  int setFullOvl=0;

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->timefp = stderr;

  setbuf(stdout,NULL);

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,"c:f:g:n:1:2:o:")) != EOF)){
      switch(ch) {
      case 'c':
	strcpy( data->File_Name_Prefix, argv[optind - 1]);
	setPrefixName = TRUE;		  
	break;
      case 'f':
	strcpy( data->Frag_Store_Name, argv[optind - 1]);
	setFragStore = TRUE;
	break;
      case 'g':
	strcpy( data->Gatekeeper_Store_Name, argv[optind - 1]);
	setGatekeeperStore = TRUE;
	break;	  
      case 'n':
	ckptNum = atoi(argv[optind - 1]);
	break;
      case '?':
	fprintf(stderr,"Unrecognized option -%c",optopt);
      default :
	errflg++;
      }
    }

    if((setPrefixName == FALSE) || (setFragStore == 0) || (setGatekeeperStore == 0))
      {
	fprintf(stderr,"* argc = %d optind = %d setFragStore = %d setGatekeeperStore = %d\n",
		argc, optind, setFragStore,setGatekeeperStore);

	exit (-1);
      }
  }

  ScaffoldGraph = 
    LoadScaffoldGraphFromCheckpoint( data->File_Name_Prefix, ckptNum, FALSE);

  // over all scfs in graph
  { int sid;
  int numScaf = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);
  for (sid = 0; sid < numScaf ; sid++){
    NodeCGW_T *scf = GetGraphNode(ScaffoldGraph->ScaffoldGraph,sid);
    if(scf->flags.bits.isDead)continue;
    printf("Scaffold %d length = %f\n",sid, scf->bpLength.mean);
  }}
}

