#ifndef GLIMMER_SPLICE_H
#define GLIMMER_SPLICE_H

struct  ICM_Score_Node_t
  {
   short int  mut_info_pos;
   float  mut_info;
   float  *prob; // was prob[ALPHABETSIZE];
  };

struct  ICM_t
{
   int  empty;
   int  model_len;     
   int  model_depth;  
   int  periodicity;
   int  num_nodes;
   struct ICM_Score_Node_t  **score;
};

struct  Fixed_Length_ICM_t
{
   int  length;
   int  max_depth;
   int  special_position;
   int  model_type;
   int  * permutation;
   struct ICM_t *sub_model;
};

extern void    readModel(struct Fixed_Length_ICM_t *fixed, const char *path);
extern double  Score_Window (struct Fixed_Length_ICM_t fixed, char * w, int left);
extern int     getModelLength(struct Fixed_Length_ICM_t fixed);

#endif    /* GLIMMER_SPLICE_H */

