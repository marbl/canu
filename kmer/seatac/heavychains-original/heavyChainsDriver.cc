// This file is part of A2Amapper.
// Copyright (c) 2004 Applera Corporation
// Author: Clark Mobarry
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


// The input file must come in "sequence axis pair segments".  
// Assuming Brian's EST mapper file format:

#ifndef __CYGWIN__
#include <libgen.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
using namespace std;

#include "heavyChainsMod.h"
#define BUFFERSIZE 1024

int get_args(
             int argc,
             char *argv[],
             map<string,string> &globals) {
  char * whoami = NULL;
  for(int ii=0;ii<argc;ii++){
    fprintf(stderr,"argv[%d]=<%s>\n", ii, argv[ii]);
  }
  int errflg = 0;
#ifndef __CYGWIN__
  whoami = strdup(basename(*argv));
#endif
  argv++;
#define OPT_ARG ( (argv[0][1])?(++(*argv)):(*(++argv)) )
  while(!errflg && *argv != NULL ) {
    if (**argv == '-') {
      (*argv)++;
      while (!errflg && **argv) {
        switch (**argv) {
#if 0
        case 'b':
          fprintf(stderr,"parse brians format\n");
          BRIANS_FORMAT = atoi(OPT_ARG);
          goto loopin;
#endif
        case 'v':
          fprintf(stderr,"parse verbose\n");
          //VERBOSE++;
          goto loopin;
          //case 'o':
          //outfile = strdup(OPT_ARG);
          //goto loopout;
        case 'g':
          fprintf(stderr,"parse globals format\n");
          {
            char * buffer = OPT_ARG;
            char * key1 = strtok(buffer,"/=");
            //if(key1 != NULL) { fprintf(stderr,"key1=<%s> %d\n", key1, strlen(key1));}
            char * key2 = strtok(NULL,"/=");
            //if(key2 != NULL) { fprintf(stderr,"key2=<%s> %d\n", key2, strlen(key2));}
            if(key1 != NULL && key2 != NULL) {
              globals[string(key1)] = string(key2);
            } else {
              fprintf(stderr,"Failed to parse ATAC global: %s\n", OPT_ARG);
              assert(0);
            }
          }
          goto loopout;
          
        case '1':
          fprintf(stderr,"parse -1\n");
          globals[string("assemblyId1")] = string(OPT_ARG);
          // fprintf(stderr,"assemblyId1=<%s>\n",
          goto loopout;
        case '2':
          fprintf(stderr,"parse -2\n");
          globals[string("assemblyId2")] = string(OPT_ARG);
          goto loopout;
        case 'S':
          fprintf(stderr,"parse -S\n");
          globals[string("heavyMinFill")] = string(OPT_ARG);
          goto loopout;
        case 'J':
          fprintf(stderr,"parse -J\n");
          globals[string("heavyMaxJump")] = string(OPT_ARG);
          goto loopout;
        case 'h':
#if 0
          usage(stdout);
#endif
          exit(0);
          goto loopin;
        default:
          // std::cerr << whoami << ": unknown flag '-" << *argv <<"'"<< std::endl;
          fprintf(stderr,"%s : unknown flag '-%s'\n", whoami, *argv);
          errflg++;
        }
      loopin:
        (*argv)++;
	fprintf(stderr,"loopin\n");
      }
    loopout:
      argv++;
	fprintf(stderr,"loopout\n");
    }  else {
      {globals[string("inpname")] = string(*argv); ++argv;}
      {globals[string("outname")] = string(*argv); ++argv;}
    }
  }
  fprintf(stderr,"return from get_args\n");
  return errflg;
} /* end of  get_args */


int main(int argc, char *argv[]) {
  map<string,string> globals;
  int VERBOSE = 3;
  char *inpfile=NULL, *outfile=NULL;
  int errflg = get_args( argc, argv, globals);
  TheStats   * ts = new TheStats( globals );
  FILE *inpF = fopen(globals[string("inpname")].c_str(),"r");
  FILE *outF= fopen(globals[string("outname")].c_str(),"w");

  if (inpfile) {
    if (!(inpF = fopen(inpfile,"r"))) {
      fprintf(stderr,"unable to open %s\n",inpfile);
      exit(2);
    }
  }
  if (outfile) {
    if (!(outF = fopen(outfile,"w"))) {
      fprintf(stderr,"unable to open %s\n",outfile);
      exit(3);
    }
  }

#if 1
  if (VERBOSE>1) {
    map<string,string>::const_iterator mi;
    fprintf(stderr,"Globals:\n");
    for(mi=globals.begin(); mi!=globals.end(); mi++){
      fprintf(stderr," globals[%s]=%s\n", (mi->first).c_str(), (mi->second).c_str());
    }
  }

  if (VERBOSE>1) {
    fprintf(stderr,"The input strings are <%s> <%s> <%s> <%s>\n",
            globals["assemblyId1"].c_str(),  globals["assemblyId2"].c_str(),
            globals["heavyMinFill"].c_str(), globals["heavyMaxJump"].c_str());
    // Note that even this print statement set the key value to non-null!
  }
#endif
  
  char *assemblyId1 = NULL;
  char *assemblyId2 = NULL;
  int verbose = 1;
  int maxjump = 100000;  // Default maximum intra-run jump allowed in a good run.
  double minscore = 100.; // Default minimum of bp filled in a good run.

  if(globals.find("assemblyId1") != globals.end()){
    fprintf(stderr,"setting assemblyId1\n");
    assemblyId1 = strdup((globals["assemblyId1"]).c_str());
  } else {
    fprintf(stderr,"error: assemblyId1 is not set\n");
    exit(3);
  }
    
  if(globals.find("assemblyId2") != globals.end()){
    fprintf(stderr,"setting assemblyId2\n");
    assemblyId2 = strdup((globals["assemblyId2"]).c_str());
  } else {
    fprintf(stderr,"error: assemblyId2 is not set\n");
    exit(3);
  }
  if(globals.find("heavyMinFill") != globals.end()){
    fprintf(stderr,"setting minscore\n");
    minscore = atof((globals["heavyMinFill"]).c_str());
  }
  if(globals.find("heavyMaxJump") != globals.end()){
    fprintf(stderr,"setting maxjump\n");
    maxjump  = atoi((globals["heavyMaxJump"]).c_str());
  }

  if (VERBOSE>1) {
    fprintf(stderr,"The input parameters are assemblyId1=%s,assemblyId2=%s,MINSCORE=%f,MAXJUMP=%d\n",
            assemblyId1, assemblyId2, minscore, maxjump);
  }

  int  old_stra1=-1, old_stra2=-1; // True strand ordinals are non-negative.
  int endOfInput = 0;
  char linebuffer[BUFFERSIZE] = {0};
  long matchid = 0;
  StrandPair * sp = new StrandPair(verbose,maxjump,minscore);

  while( !endOfInput ) {
    endOfInput = 1;

    int  new_stra1=-1, new_stra2=-1;
    int xln, yln;
    int tmp_xlo, tmp_ylo, tmp_filled;
    char tmp_ori;

    if(NULL != fgets(linebuffer,BUFFERSIZE,inpF)) {
      //if(PANICDEBUG){ printf("linebuffer=%s",linebuffer);}

      if(linebuffer[0] == '-' || linebuffer[0] == 'M'){
        if(linebuffer[0] != 'M') {
          // Now input using Brian's EST mapper format ....
          int scanc = sscanf(linebuffer,
                             "-%c -e %d %d %d -D %d %d %d -F %d\n",
                             &tmp_ori,
                             &new_stra1,&tmp_xlo,&xln,&new_stra2,&tmp_ylo,&yln,&tmp_filled);
          assert(tmp_ori == 'f' || tmp_ori == 'r');
          assert(scanc >= 7);
          endOfInput = 0;
	} else {
	  // Now input using ATAC format ....
	  // M r r6379 . B33:9 134205509 488863 1 R27M4:9 128704627 484912 -1
	  char classCode, subtype;
	  char selfId[100], parentId[100];
          char new_ass1[100], new_ass2[100];
	  int xfl, yfl;
          int scanc = sscanf(linebuffer,
                             "%c %c %s %s "
                             "%s %d %d %d %s %d %d %d\n",
                             &classCode, &subtype, selfId, parentId,
                             new_ass1,&tmp_xlo,&xln,&xfl,
                             new_ass2,&tmp_ylo,&yln,&yfl);

          // Check if new_ass1 equals assemblyId1 etc.
          assert(scanc >= 12);
#if 0
          printf("classCode=%c\n", classCode);
          printf("subtype  =%c\n", subtype);
          printf("selfId   =%s\n", selfId);
          printf("parentId =%s\n", parentId);
          printf("new_ass1 =%s\n", new_ass1);
          printf("xfl      =%d\n", xfl);
          printf("new_ass2 =%s\n", new_ass2);
          printf("yfl      =%d\n", yfl);
#endif
          assert(xfl == 1 || xfl == -1);
          assert(yfl == 1 || yfl == -1);
          tmp_ori = (xfl == yfl ? 'f' : 'r');
          {
            char *p;
            for(p = new_ass1; *p != 0; p++){ if(*p == ':'){ p++; break;}}
            sscanf(p,"%d", &new_stra1);
          }
          {
            char *p;
            for(p = new_ass2; *p != 0; p++){ if(*p == ':'){ p++; break;}}
            sscanf(p,"%d", &new_stra2);
          }
          endOfInput=0;
        }
#if 0
        if(!endOfInput){
          // Assert that the data came in with the axis numerically or lexigraphically sorted:
          assert((old_stra1 <= new_stra1) ||
                 (old_stra1 == new_stra1 && old_stra2 <= new_stra2));
        }
#endif
        
      } else if(linebuffer[0] == '#') {
        fprintf(outF,"%s",linebuffer);
        endOfInput = 0;
      } else if(linebuffer[0] == '/') {
        //fprintf(outF,"%s",linebuffer);
        char * key1 = strtok(linebuffer,"/=");
        //if(key1 != NULL) { fprintf(stderr,"key1=<%s> %d\n", key1, strlen(key1));}
        char * key2 = strtok(NULL,"\n");
        //if(key2 != NULL) { fprintf(stderr,"key2=<%s> %d\n", key2, strlen(key2));}
        if(key1 != NULL && key2 != NULL) {
          globals[string(key1)] = string(key2);
        }
        endOfInput = 0;
      } else if(linebuffer[0] == '!') {
        fprintf(outF,"%s",linebuffer);
        endOfInput = 0;
      } else {
        fprintf(stderr,"Unrecognized input\n");
        fprintf(outF,"%s",linebuffer);
        endOfInput = 0;
      }
    }

    if ( (new_stra1 != old_stra1 || new_stra2 != old_stra2 || endOfInput)
	 && sp->size() >0 ) {
      sp->process();
      matchid = sp->print(outF, 0, assemblyId1, assemblyId2, matchid );
      ts->add(sp);
      delete sp;
      sp = new StrandPair(verbose,maxjump,minscore);
    }

    if(linebuffer[0] == '-' || linebuffer[0] == 'M') {
      // If we just processed a new point.
      sp->addHit(
		 tmp_ori,
                 new_stra1, tmp_xlo, xln,
		 new_stra2, tmp_ylo, yln,
		 tmp_filled);

      old_stra1 = new_stra1;
      old_stra2 = new_stra2;
    }

  }

  ts->add(sp);
  delete sp;

  ts->print(outF);
  delete ts;

  if (inpfile) fclose(inpF);
  if (outfile) fclose(outF);
}
