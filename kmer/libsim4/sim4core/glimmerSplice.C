#include "sim4.H"


const char      ALPHA_STRING [] = "acgt";
const int       DEFAULT_PERIODICITY = 3;
const int       DEFAULT_MODEL_DEPTH = 7;

const int       DEFAULT_MODEL_LEN = 12;
const int       ALPHABETSIZE = 4;
const int       MAX_ERROR_MSG_LEN = 1000;
const int       ICM_VERSION_ID = 200;
const unsigned int NUM_FIXED_LENGTH_PARAMS = 6;

const int       ID_STRING_LEN = 400;
//const unsigned  NUM_FIXED_LENGTH_PARAMS = 6;

#define  PARENT(x) ((int) ((x) - 1) / ALPHABETSIZE)

int  Filter(char Ch);
int  Subscript(char ch);
void Permute_String(char * s, int * perm, int n);
void *Safe_malloc(size_t len, const char * src_fname, size_t line_num);
void *Safe_realloc(void * q, size_t len, const char * src_fname, size_t line_num);
void *Safe_calloc(size_t n, size_t len, const char * src_fname, size_t line_num);

int  Int_Power(int a, int b);
void Input(struct ICM_t *p, FILE *fp, int model_len,int model_depth, int periodicity);
int  Get_Model_Depth (struct ICM_t p) { return  p.model_depth; }
int  Get_Model_Len (struct ICM_t p)   { return  p.model_len; }
int  Get_Periodicity (struct ICM_t p) { return  p.periodicity; }
double  Full_Window_Prob(struct ICM_t p, const char * string, int frame);

int     Get_Length(struct Fixed_Length_ICM_t fixed){  return  fixed.length; }
double  Score_Window (struct Fixed_Length_ICM_t fixed, char * w, int left);
int  getModelLength(struct Fixed_Length_ICM_t fixed) { return  fixed.length;}
int  getModelType(struct Fixed_Length_ICM_t fixed)   { return  fixed.model_type;}
int  getSpecialPosition(struct Fixed_Length_ICM_t fixed){ return  fixed.special_position; }

void readModel(struct Fixed_Length_ICM_t *fixed, const char *path);

void readModel(struct Fixed_Length_ICM_t *fixed, const char *path) 
  {
   FILE  * fp;
   char  line [ID_STRING_LEN];
   int  param [NUM_FIXED_LENGTH_PARAMS];
   int  i;

   if ((fp = fopen (path, "r"))==NULL) {
      fprintf(stderr, "Error: Could not open Glimmer model file for reading (%s).\n", path); 
      exit(1);
   }

   fread (line, sizeof (char), ID_STRING_LEN, fp); // skip the text header line

   if  (fread (param, sizeof (int), NUM_FIXED_LENGTH_PARAMS, fp)
          != NUM_FIXED_LENGTH_PARAMS)
       {
        fprintf (stderr, "ERROR reading file \"%s\"\n", path);
        exit (-1);
       }

   if  (ICM_VERSION_ID != param [0])
       {
        fprintf (stderr, "Bad ICM version = %d  should be %d\n",
                 param [0], ICM_VERSION_ID);
        exit (-1);
       }
   if  (ID_STRING_LEN != param [1])
       {
        fprintf (stderr, "Bad ID_STRING_LEN = %d  should be %d\n",
                 param [1], ID_STRING_LEN);
        exit (-1);
       }
   (*fixed).length = param [2];
   (*fixed).max_depth = param [3];
   (*fixed).special_position = param [4];
   (*fixed).model_type = param [5];

   (*fixed).permutation = (int *) Safe_malloc((*fixed).length*sizeof(int), __FILE__,__LINE__);

   for(i=0;i<(*fixed).length;i++) {
     (*fixed).permutation[i] = 0;
   }
   fread ((*fixed).permutation, sizeof (int), (*fixed).length, fp);

   (*fixed).sub_model = (struct ICM_t *) Safe_malloc
             ((*fixed).length * sizeof (struct ICM_t ), __FILE__, __LINE__);

   for  (i = 0;  i < (*fixed).length;  i ++) 
   {
      (*fixed).sub_model[i].score = (struct ICM_Score_Node_t * *)
          Safe_calloc (1, sizeof (struct ICM_Score_Node_t *), __FILE__, __LINE__);
//      for  (j = 0;  j < 1;  j ++) {
//        (*fixed).sub_model[i].score[j] = (struct ICM_Score_Node_t *)
//                   Safe_calloc (12, sizeof (struct ICM_Score_Node_t),
//                   __FILE__, __LINE__);
//        for(k=0;k<4;k++) {
//          (*fixed).sub_model[i].score[j][k].prob = (float *)
//             Safe_calloc (4, sizeof(float), __FILE__, __LINE__);
//        }
//      } 
//
   }

   for  (i = 0;  i < (*fixed).length;  i ++)
   {
      Input(&((*fixed).sub_model[i]), fp,1,0,1);
   }
}

//  Input the contents of this model from  fp , which has already been opened.

void  Input(struct ICM_t *p, FILE *fp, int model_len,int model_depth, int periodicity)
{
 char  line [ID_STRING_LEN];
 int  param [NUM_FIXED_LENGTH_PARAMS];
 int  node_id;
 int  prev_node;
 int  period;
 int  i,j;
 (*p).model_len = model_len;
 (*p).model_depth = model_depth;
 (*p).periodicity = periodicity;

 (*p).empty = 1;

 // skip the text header line
 if  (fread (line, sizeof (char), ID_STRING_LEN, fp) != (unsigned) (ID_STRING_LEN))
 {
   fprintf (stderr, "ERROR reading ICM header\n");
   exit (-1);
 }

 if  (fread (param, sizeof (int), NUM_FIXED_LENGTH_PARAMS, fp) != NUM_FIXED_LENGTH_PARAMS)
 {
    fprintf (stderr, "ERROR reading parameters\n");
    exit (-1);
 }

 if  (ICM_VERSION_ID != param [0])
 {
   fprintf (stderr, "Bad ICM version = %d  should be %d\n", param [0], ICM_VERSION_ID);
   exit (-1);
 }
 if  (ID_STRING_LEN != param [1])
 {
    fprintf (stderr, "Bad ID_STRING_LEN = %d  should be %d\n",
             param [1], ID_STRING_LEN);
    exit (-1);
 }

 (*p).model_len   = param [2];
 (*p).model_depth = param [3];
 (*p).periodicity = param [4];
 (*p).num_nodes   = param [5];


 (*p).score = (struct ICM_Score_Node_t **) Safe_malloc
             ((*p).periodicity * sizeof (struct ICM_Score_Node_t *), __FILE__, __LINE__);
 for  (i = 0;  i < (*p).periodicity;  i ++) {
   (*p).score [i] = (struct ICM_Score_Node_t *) Safe_calloc
                 ((*p).num_nodes, sizeof (struct ICM_Score_Node_t), __FILE__, __LINE__);
   for(j=0;j<(*p).num_nodes;j++) {
     (*p).score[i][j].prob = (float *) Safe_malloc(ALPHABETSIZE*sizeof(float), __FILE__, __LINE__);
   }
 }


 period = -1;
 prev_node = 0;
 while  (fread (& node_id, sizeof (int), 1, fp) != 0)
 {
   if  (node_id < 0) break;
   if  (node_id == 0)  period++;

   // read in the probabilities
   if  (fread ((*p).score [period] [node_id] . prob, sizeof (float), ALPHABETSIZE, fp) != 
             (unsigned) (ALPHABETSIZE))
   {
      fprintf (stderr, "ERROR reading icm node = %d  period = %d\n", node_id, period);
      exit (-1);
   }


   // read in the max mutual information position
   if  (fread (& ((*p).score [period] [node_id] . mut_info_pos), sizeof (short int), 1, fp) != 1)
   {
      fprintf (stderr, "ERROR reading mut_info_pos for node = %d  period = %d\n", node_id, period);
      exit (-1);
   }

   // check for cut nodes
   if  (node_id != 0 && prev_node != node_id - 1)
      for  (i = prev_node + 1;  i < node_id;  i ++)
         (*p).score [period] [i] . mut_info_pos = -2;

   if  (node_id == 0 && period > 0)
      for  (i = prev_node + 1;  i < (*p).num_nodes;  i ++)
         (*p).score [period - 1] [i] . mut_info_pos = -2;

   prev_node = node_id;
 }

 if  (period != periodicity - 1)
 {
   fprintf (stderr, "ERROR:  Too few nodes for periodicity = %d\n", periodicity);
   exit (-1);
 }

 // check for cut nodes in last period
 if  (prev_node != (*p).num_nodes - 1)
     for  (i = prev_node + 1;  i < (*p).num_nodes;  i ++)
        (*p).score [period] [i] . mut_info_pos = -2;

 (*p).empty = 0;
}

//  Rearrange the characters in  s  according
//  to the permutation in  perm .

void  Permute_String(char * s, int * perm, int n)
  {
   static char  * buff = NULL;
   static int  buff_len = 0;
   int  i;

   if  (n > buff_len)
       {
        buff = (char *) Safe_realloc (buff, n, __FILE__, __LINE__);
        buff_len = n;
       }

   for  (i = 0;  i < n;  i ++)
     buff [i] = s [perm [i]];
   strncpy (s, buff, n);

   return;
  }

//  Return a single  a, c, g or t  for  Ch .

int  Filter(char Ch)
  {
   switch  (tolower (Ch))
     {
      case  'a' :
      case  'c' :
      case  'g' :
      case  't' :
        return  Ch;
      case  'r' :     // a or g
        return  'g';
      case  'y' :     // c or t
        return  'c';
      case  's' :     // c or g
        return  'c';
      case  'w' :     // a or t
        return  't';
      case  'm' :     // a or c
        return  'c';
      case  'k' :     // g or t
        return  't';
      case  'b' :     // c, g or t
        return  'c';
      case  'd' :     // a, g or t
        return  'g';
      case  'h' :     // a, c or t
        return  'c';
      case  'v' :     // a, c or g
        return  'c';
      default :       // anything
        return  'c';
    }
  }

//  Return the subscript equivalent (used in offsets of the
//  model) for character  ch .

int  Subscript(char ch)
  {
   char  * p;

   p = strchr ((char *)ALPHA_STRING, tolower (Filter (ch)));
   if  (p == NULL)
       {
        fprintf (stderr, "ERROR:  Bad character %c in subscript conversion",
                 ch);
        exit (-1);
       }

   return  (int) (p - ALPHA_STRING);
  }


//  Return the log-probability of the last character in the first
//  model_len  bases of  string  conditioned on the preceding characters
//  using the entries in  score [frame] .

double  Full_Window_Prob (struct ICM_t icm, const char * string, int frame) 
  {
   double  prob;
   int  num_node, i, pos, sub;

   num_node = 0;

   for  (i = 0;  i < icm.model_depth;  i ++)
     {
      pos = icm.score [frame] [num_node] . mut_info_pos;

      if  (pos == -1)
          break;

      if  (pos < -1)  // No information here or below in tree, go back up
                      // Shouldn't happen
          {
           num_node = PARENT (num_node);
           pos = icm.score [frame] [num_node] . mut_info_pos;
           break;
          }

      sub = Subscript (string [pos]);

      num_node = (num_node * ALPHABETSIZE) + sub + 1;
     }

   pos = icm.score [frame] [num_node] . mut_info_pos;
   if  (pos < -1)
       {
        num_node = PARENT (num_node);
        pos = icm.score [frame] [num_node] . mut_info_pos;
       }

   sub = Subscript (string [icm.model_len - 1]);

   prob = (double) icm.score [frame] [num_node] . prob [sub];

   if  (pos < -1)
       {
        fprintf (stderr, "WARNING:  prob = %.4f  pos = %d in  Full_Window_Prob\n",
                 prob, pos);
        fprintf (stderr, "num_node = %d\n",
                 num_node);
       }

   return  prob;
  }


//  Return the score of this model on string  w
double  Score_Window (struct Fixed_Length_ICM_t fixed, char * w, int left)
{
   static char  * buff = NULL;
   static int  buff_len = 0;
   double  score = 0.0;
   int  i;

   if  (fixed.length > buff_len)
   {
     buff = (char *) Safe_realloc (buff, fixed.length+1, __FILE__, __LINE__);
     buff_len = fixed.length;
   }

   strncpy (buff, w, fixed.length);
// strncpy (buff, w, left);
// strncpy (buff+left, w+left+2, fixed.length-left);

   if  (fixed.permutation != NULL)
       Permute_String (buff, fixed.permutation, fixed.length);

   for  (i = 0;  i < fixed.length;  i ++)
     {
      if  (buff [i] == '\0')
      {
        fprintf (stderr, "ERROR:  String \"%s\" too short in Score_Window\n",
                  buff);
        exit (-1);
      }
      score += Full_Window_Prob (fixed.sub_model[i], buff, 0);
     
     }

   return  score;
}

void  Clean_Exit
    (const char * msg, const char * src_fname, size_t line_num)

//  Write string  msg  to  stderr  and also a line indicating
//  the error happen in source file  src_fname  at line  line_num
//  if they are not  NULL  and  0  respectively.
//  Then exit with an error condition.

  {
   fprintf (stderr, "%s\n", msg);
   if  (src_fname != NULL)
       fprintf (stderr, "  in file  %s", src_fname);
   if  (line_num != 0)
       fprintf (stderr, "  at line  %lu", (long unsigned) (line_num));
   fprintf (stderr, "  errno = %d\n", errno);

   exit (-1);
  }

void *  Safe_calloc
    (size_t n, size_t len, const char * src_fname, size_t line_num)

//  Allocate and return a pointer to enough memory to hold an
//  array with  n  entries of  len  bytes each.  All memory is
//  cleared to 0.  If fail, print a message and exit, assuming the
//  call came from  source file  src_fname  at line  line_num .

  {
   void  * p;
   char  Clean_Exit_Msg_Line [MAX_ERROR_MSG_LEN];

   p = calloc (n, len);
   if  (p == NULL)
       {
        sprintf (Clean_Exit_Msg_Line,
                 "ERROR:  calloc failed  %lu x %lu",
                 (long unsigned) (n), (long unsigned) (len));
        Clean_Exit (Clean_Exit_Msg_Line, src_fname, line_num);
       }

   return  p;
  }


void *  Safe_malloc
    (size_t len, const char * src_fname, size_t line_num)

//  Allocate and return a pointer to  len  bytes of memory.
//  If fail, print a message and exit, assuming the call came from
//  source file  src_fname  at line  line_num .

  {
   void  * p;
   char  Clean_Exit_Msg_Line [MAX_ERROR_MSG_LEN];

   p = malloc (len);
   if  (p == NULL)
       {
        sprintf (Clean_Exit_Msg_Line,
                 "ERROR:  malloc failed  %lu  bytes",
                 (long unsigned) (len));
        Clean_Exit (Clean_Exit_Msg_Line, src_fname, line_num);
       }

   return  p;
  }

void *  Safe_realloc
    (void * q, size_t len, const char * src_fname, size_t line_num)

//  Reallocate memory for  q  to  len  bytes and return a
//  pointer to the new memory.  If fail, print a message and exit,
//  assuming the call came from source file  src_fname  at line  line_num .

  {
   char  Clean_Exit_Msg_Line [MAX_ERROR_MSG_LEN];

   void  * p;
   p = realloc (q, len);
   if  (p == NULL)
       {
        sprintf (Clean_Exit_Msg_Line,
                 "ERROR:  realloc failed  %lu  bytes",
                 (long unsigned) (len));
        Clean_Exit (Clean_Exit_Msg_Line, src_fname, line_num);
       }

   return  p;
  }


int  Int_Power(int a, int b)
{
   int  result = 1;
   int  p = a;

   while  (b > 0)
     {
      if  (b & 1)
          result *= p;
      p = p * p;
      b >>= 1;
     }

   return  result;
}

