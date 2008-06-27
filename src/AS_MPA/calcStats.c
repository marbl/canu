/**************************************************************************
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
/* $Id: calcStats.c,v 1.6 2008-06-27 06:29:16 brianwalenz Exp $ */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <float.h>
#include <math.h>

/*
  flags:
    -c means two columns per line
      col1 value
      col2 number of instances of the value
    else
      one column per line of values
 */

void usage(char * progname)
{
  fprintf(stderr, "Usage: %s [-h] [-c] [-f file] [-l m] [-g n]\n"
          "  -h        print help\n"
          "  -c        two-column data. 1st col is value, 2nd is count\n"
          "              default is one-column data - of values\n"
          "  -f file   name of file containing data\n"
          "              default is stdin\n"
          "  -l m      include only numbers less than m\n"
          "  -g n      include only numbers greater than n\n",
          progname);
  exit(1);
}


int main(int argc, char ** argv)
{
  int type = 1;
  char * fn = NULL;
  FILE * fp;
  char line[2048];
  unsigned long n;
  double sumX;
  double sumX2;
  double val;
  double mean;
  double sigma;
  double greaterThan = DBL_MIN;
  double lessThan = DBL_MAX;
  {
    int ch;
    while((ch = getopt(argc, argv, "chf:g:l:")) != EOF)
    {
      switch(ch)
      {
        case 'l':
          lessThan = atof(optarg);
          break;
        case 'g':
          greaterThan = atof(optarg);
          break;
        case 'c':
          type = 2;
          break;
        case 'h':
          usage(argv[0]);
          break;
        case 'f':
          fn = optarg;
          break;
        default:
          usage(argv[0]);
          break;
      }
    }
  }
  if(fn == NULL)
  {
    fp = stdin;
  }
  else
  {
    fp = fopen(fn, "r");
    if(fp == NULL)
      usage(argv[0]);
  }

  sumX = sumX2 = n = 0;
  if(type == 1)
  {
    while(fgets(line, 2047, fp) != NULL)
    {
      sscanf(line, "%lf", &val);

      if(val <= greaterThan || val >= lessThan) continue;

      n++;
      sumX += val;
      sumX2 += val * val;
    }
  }
  else if(type == 2)
  {
    unsigned long instances = 0;

    while(fgets(line, 2047, fp) != NULL)
    {
      sscanf(line, "%lf\t%lu", &val, &instances);

      if(val <= greaterThan || val >= lessThan) continue;

      n += instances;
      sumX += val * instances;
      sumX2 += val * val * instances;
    }

  }
  fclose(fp);

  mean = sumX / n;
  sigma = sqrt((sumX2 - n * mean * mean) / n);
  fprintf(stdout, "%f = mean\n%f = sigma\n%lu = n\n", mean, sigma, n);
  // fprintf(stderr, "sumX2 = %e\nsumX = %e\n", sumX2, sumX);
  return 0;
}

