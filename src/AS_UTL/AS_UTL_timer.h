
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
/* 	$Id: AS_UTL_timer.h,v 1.2 2004-09-23 20:25:29 mcschatz Exp $	 */
/***************************************************************************
 *  AS_UTL_timer
 *  
 *  Saul A. Kravitz 7/99
 *
 *  Simple timer object using ISO-C's clock() API.
 *  Resolution of clock is machine dependent (CLOCKS_PER_SEC)
 *  
 **************************************************************************/
#ifndef AS_UTL_TIMER_H
#define AS_UTL_TIMER_H

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#define GRANULARITY_CUTOFF 2147 // 2147 seconds is where the clock_t wraps around

typedef struct{
  uint64 cumTime;   // cumulative times for all start/stop sequences since last reset

  // Unfortunately, clock_t is only 32 bits
  // Each tick is scaled by CLOCKS_PER_SEC
  clock_t startTime; // Beginning of current start/stop sequence
  clock_t stopTime;  // Time of last stop. stop-start was then added to cumTime;
  time_t startTimeInSeconds;
  time_t stopTimeInSeconds;	  
  int running;
  int64 cycles;
}TimerT;

struct rusage timeUsage;

/* Initialize a Timer struction */
static void InitTimerT(TimerT *timer){
  timer->cycles = 0;
  timer->cumTime = 0;
  timer->startTime = 0;
  timer->stopTime = 0;
  timer->startTimeInSeconds = 0;
  timer->stopTimeInSeconds = 0;
  timer->running = 0;
}

#define ResetTimerT InitTimerT

/* Start the Timer */
static void StartTimerT(TimerT *timer){
  timer->running = 1;
  timer->startTime = clock();
  timer->stopTime = 0;

  getrusage(RUSAGE_SELF, &timeUsage);
  timer->startTimeInSeconds = timeUsage.ru_utime.tv_sec + timeUsage.ru_stime.tv_sec;;
  timer->stopTimeInSeconds = 0;
}

/* Stop the Timer
   Only has impact if timer is running (after a StartTimerT).
   Increments cycles count
   Saves stop time in stopTime
   Increments cumulative timer
   Returns lap time as a double
*/
static double StopTimerT(TimerT *timer){
  clock_t delta;
  long stopTime, startTime;
  time_t deltaInSeconds;
  long stopTimeInSeconds, startTimeInSeconds;
  
  if(timer->running == 0)
    return (double)0.0;
  timer->running = 0;
  timer->cycles++;
  timer->stopTime = clock();
  getrusage(RUSAGE_SELF, &timeUsage);
  timer->stopTimeInSeconds = timeUsage.ru_utime.tv_sec + timeUsage.ru_stime.tv_sec;

  startTime = (timer->startTime);
  stopTime = (timer->stopTime);
  delta = (stopTime - startTime);

  startTimeInSeconds = (timer->startTimeInSeconds);
  stopTimeInSeconds = (timer->stopTimeInSeconds);
  deltaInSeconds = (stopTimeInSeconds - startTimeInSeconds);

  if (deltaInSeconds < GRANULARITY_CUTOFF) 
  {
	timer->cumTime += (uint64)delta;
	return (double)delta/ (double)CLOCKS_PER_SEC;
  }
  // else
  {
	timer->cumTime += (uint64)deltaInSeconds * CLOCKS_PER_SEC;
	return deltaInSeconds;
  }
}


/* Cumulative timer time in seconds */
static double TotalTimerT(TimerT *timer, long *cycles){

  if(cycles)
    *cycles = timer->cycles;
  return (double)timer->cumTime/ (double)CLOCKS_PER_SEC ;

}

/* Last cycles' time time in seconds */
static double LapTimerT(TimerT *timer){
  clock_t delta = (timer->stopTime - timer->startTime);
  time_t deltaInSeconds = (timer->stopTimeInSeconds - timer->startTimeInSeconds);

  if(timer->running)
    return 0.0;

  if (deltaInSeconds < GRANULARITY_CUTOFF) {
    return (double)delta/ (double)CLOCKS_PER_SEC ;
  }
  // else 
  return (double) deltaInSeconds;
}



#endif






