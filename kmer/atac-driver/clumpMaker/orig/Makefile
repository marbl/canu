#######################################################################
#
#           Copyright (c) 2002 TIGR Center for the Advancement of Genomics
#                          All rights reserved.
#
#######################################################################

######################## Variable Definitions #########################
LOCAL_WORK = $(shell cd ../; pwd)

CPPSOURCES = tClumpHits.cc clumpHits.cc

SOURCES = $(CPPSOURCES)

OBJECTS = $(CPPSOURCES:%.cc=%.o)

PROGS = clumpHits

LIBRARIES =

AR=ar
ARFLAGS=-rvs



########################### Include variables etc.# #################

LIBDIRS =  -L$(LEDAROOT) -L/usr/X11R6/lib -L/usr/local/lib -L/usr/lib

LIBS =  -lW -lP -lG -lL -lD3 -lX11 -lm -lstdc++


INCLUDES = -I- -I. -I$(LOCAL_WORK)/src -I$(LEDAROOT)/incl 

CPPFLAGS = -Wall -D__unix__ -O3
#CPPFLAGS = -Wall -D__unix__ -g

CPP = g++ $(CPPFLAGS) $(INCLUDES) $(LIBDIRS)

CC = g++ $(CPPFLAGS) $(INCLUDES) $(LIBDIRS)

################### Project Rule Makefile Inclusion ###################
# include Makeconf

######################### Target Definitions ##########################
all: $(OBJECTS) $(CPPPROGS)

del:		
		rm -f $(OBJECTS) $(CPPPROGS) *~

clean: del
		rm -f *~ *core *.core *.o


$(OBJECTS):	%.o:%.cc

%.o:%.cc
		$(CPP) -c $<

$(PROGS):	
		$(CPP) -o $@ $+  $(LIBDIRS) $(LIBS)


lib: $(OBJECTS) Makefile
	$(AR) $(ARFLAGS) $(OBJECTS)

###############Specific Program Dependencies############################

clumpHits: clumpHits.o tClumpHits.o

tClumpHits.o: tClumpHits.cc tClumpHits.h