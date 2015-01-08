
#ifndef MULTIALIGNUNITIG_H
#define MULTIALIGNUNITIG_H

static const char *rcsid_MULTIALIGNUNITIG_H = "$Id:  $";

#include "AS_global.H"

class multiAlignUnitigOptions {
  multiAlignUnitigOptions() {
    verbose      = 0;      //  extern uint32 VERBOSE_MULTIALIGN_OUTPUT;
    splitAlleles = true;   //  #define CNS_OPTIONS_SPLIT_ALLELES_DEFAULT  1
    smoothWindow = 11;     //  #define CNS_OPTIONS_MIN_ANCHOR_DEFAULT    11
    doPhasing    = true;   //  #define CNS_OPTIONS_DO_PHASING_DEFAULT     1
  };
  ~multiAlignUnitigOptions() {
  };

  uint32  verbose;
  bool    splitAlleles;
  uint32  smoothWindow;
  bool    doPhasing;
};



//  uint32  printWidth;
//  uint32  printSpacing;


//extern uint32 MULTIALIGN_PRINT_WIDTH;       //  In MultiAlignPrint.c
//extern uint32 MULTIALIGN_PRINT_SPACING;


//typedef enum {
//  CNS_SMOOTH = 1, // only eliminate pairwise construction artifacts
//  CNS_POLYX  = 2, // align poly-X regions
//  CNS_INDEL  = 4  // push apart mushed block indels
//}  CNS_RefineLevel;



bool
generateMultiAlignment(tgTig                     *tig,
                       gkStore                   *gkpStore,
                       multiAlignUnitigOptions&   options,
                       int32                     *failed);

//  Options to things in MultiAligment_CNS.c
//extern int32 allow_neg_hang;

#endif
