
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
#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>


#define  HASH_KMER_SKIP    0
  // Number of bytes skipped between kmers put in the hash table.


static int  Binomial_Hit_Limit
    (double p, double n, double prob_bound);



int  main
    ()

  {
   double  genome_len, copies, prob;
   double  total_len, num_frags;
   double  avg_frag_kmers, genome_kmers;
   double  p, s;
   int  kmer_len, hit_limit;
   
   printf ("Enter genome len, copy number and prob\n");
   scanf ("%lf %lf %lf", & genome_len, & copies, & prob);

   printf ("Enter total len, num frags and kmer len\n");
   scanf ("%lf %lf %d", & total_len, & num_frags, & kmer_len);

   avg_frag_kmers
       = (total_len / num_frags - kmer_len)
            / (1 + HASH_KMER_SKIP);
   genome_kmers = 2 * (genome_len - (kmer_len - 1));

   p = 2 * copies * avg_frag_kmers/genome_kmers;
   s = sqrt (num_frags * p * (1.0 - p));
   hit_limit = 1 + Binomial_Hit_Limit (p, num_frags - 1, prob);
   printf ("Hit limit = %d for k-mer/reverse complement counted together\n",
           hit_limit);
   printf ("   p = %e   s = %e\n", p, s);
           

   p = copies * avg_frag_kmers/genome_kmers;
   s = sqrt (num_frags * p * (1.0 - p));
   hit_limit = 1 + Binomial_Hit_Limit (p, num_frags - 1, prob);
   printf ("Hit limit = %d for k-mers counted separately\n",
           hit_limit);
   printf ("   p = %e   s = %e\n", p, s);

   

   return  0;
  }



static int  Binomial_Hit_Limit
    (double p, double n, double prob_bound)

//  Return  k  such that the probability of  >= k  successes in
//   n  trials is  < prob_bound , where  p  is the probability
//  of success of each trial.

  {
   double  lambda, target, sum, term, q;
   int  i;

   lambda = n * p;
   q = 1.0 - p;

   if  (lambda <= 5.0 && n >= 50.0)
       {  // use Poisson approximation
        target = (1.0 - prob_bound) * exp (lambda);
        sum = term = 1.0;
        for  (i = 1;  sum <= target && i < 50;  i ++)
          {
           term = term * lambda / (double) i;
           sum += term;
          }
        if  (sum > target)
            return  i;
       }
   if  (n >= 30.0)
       {  // use Normal approximation
        double  t, z;
        double  c [3] = {2.515517, 0.802853, 0.010328};
        double  d [4] = {1.0, 1.432788, 0.189269, 0.001308};

        if  (prob_bound <= 0.5)
            t = sqrt (-2.0 * log (prob_bound));
          else
            t = sqrt (-2.0 * log (1.0 - prob_bound));

        z = t - ((c [2] * t + c [1]) * t + c [0])
                  / (((d [3] * t + d [2]) * t + d [1]) * t + d [0]);

        if  (prob_bound <= 0.5)
            target = z;
          else
            target = -z;

        return  (int) ceil (lambda + sqrt (lambda * q) * target);
       }
     else
       {  // brute force
        target = 1.0 - prob_bound;
        sum = term = pow (q, n);
        for  (i = 1;  sum <= target && i < n;  i ++)
          {
           term *= (n + 1 - i) / i;
           term *= p / q;
           sum += term;
          }
        return  i;
       }
  }
