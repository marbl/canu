=============================================================================

This file is part of Celera Assembler Software, a program that 
assembles whole-genome shotgun reads into contigs and scaffolds.
Copyright 1999-2004, Applera Corporation. All rights reserved.

The Celera Assembler Software (Copyright 1999-2004, Applera Corporation. All 
rights reserved.) (the "Software") is covered by one or more U.S. patents and is 
being made available free of charge by Applera Corporation subject to the terms 
and conditions of the GNU General Public License, version 2, as published by the 
Free Software Foundation (the "GNU General Public License"). The Software is 
free software; you can redistribute it and/or modify it under the terms of the 
GNU General Public License. The Software is distributed in the hope that it will 
be useful, but WITHOUT ANY WARRANTIES, EXPRESS OR IMPLIED (INCLUDING, WITHOUT 
LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR 
PURPOSE). You should have received (LICENSE.txt) a copy of the GNU General Public 
License along with the Software; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

=============================================================================

The source code for Celera Assembler can be obtained from Applera Corporation.
The prefered method to request source is on the Web at 
  http://myscience.appliedbiosystems.com/publications/compass/index.jsp
Another way to request source is to send e-mail to
  genomesoftware@celera.com
If electronic methods fail, you can request source by mail sent to 
  Informatics Research, Celera, 45 West Gude Drive, Rockville, MD 20850

=============================================================================

The Celera Assembler is provided "as-is".
Compiling it, parameterizing it, and running it may be non-trivial tasks.
Requests for consulting services will be considered if addressed to
  genomesoftware@celera.com
or
  Informatics Research, Celera, 45 West Gude Drive, Rockville, MD 20850

=============================================================================

When referencing this work, please cite the following papers.

Istrail et al. (2004) 
Whole-Genome Shotgun Assembly and Comparison of Human Genome Assemblies. 
PNAS 101 1916-1921.

Myers et al. (2000) 
A Whole-Genome Assembly of Drosophila. 
Science 287 2196-2204.

Venter et al. (2001) 
The Sequence of the Human Genome. 
Science 291 1304-1351.

=============================================================================

The Celera Assembler consists of 12 main components.
Each component was designed to run in a unix environment such as
Compaq(R) Tru64(R) or IBM(R) AIX(R) or Red Hat Linux.
Some components require large amounts of RAM for large datasets.
Each component can be compiled from the C source using
a compiler such as GNU gcc.
Each component can be built using the makefiles provided, by running
a make utility such as GNU gmake in the src subdirectory.
The makefiles assume an environment variable $OSTYPE is defined
as either "osf3.2", "aix", or "linux".
Some documentation is provided in the doc subdirectory.
Some useful scripts are provided in the scripts subdirectory.
Users should start with the wga.pl Perl script in the examples subdirectory.
This script executes the Celera Assembler pipeline on a simulated dataset.
The dataset can be generated with the component called celsim.
The dataset is also provided with the Celera Assembler distribution.
Note the script does not execute every component of the Celera Assembler,
and in particular it does not demonstrate multi-file operation.

=============================================================================

