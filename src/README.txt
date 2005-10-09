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

The makefiles in this and sub-directories were created to work with
gmake. The c_make.as file sets several variables and options based on
the OSTYPE environment variable. This version of the makefiles check
for linux, aix, and osf3.2. If you wish to build the assembler, you
may need to either edit the OSTYPE environment variable or the
c_make.as file.


Before building the assembler, you'll need to set the SITE_NAME make variable.
This is used for controlling site specific configurations, especially for
the terminiator and EUIDs. To set it, either create a file called site_name.as
with a single line:

SITE_NAME=Your_site_name

or run make as 'make SITE_NAME=Your_site_name'


Scripts (perl, csh, bash, etc) in the scripts and example directories
as well as in subdirectories of src specify the path of the
shell/interpreter for the script. You may need to edit these to match
the locations of interpreters and shells installed on your UNIX
system.
