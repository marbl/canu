#
###########################################################################
#
# This file is part of Celera Assembler, a software program that 
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received (LICENSE.txt) a copy of the GNU General Public 
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
###########################################################################
#

## Set this to JCVI to use CURL to access JCVI's UID server.
## Set this to TIGR to use SOAP to access TIGR's UID server.
## Add your own favorite server....

include $(LOCAL_WORK)/src/site_name.as

UID_SERVER = $(SITE_NAME)

ifeq ($(UID_SERVER), TIGR)
  USE_SOAP_UID = 1
endif
