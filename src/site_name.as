#
# $Id: site_name.as,v 1.8 2008-03-11 03:34:43 brianwalenz Exp $
#
# Set the name of the site/location where the assembler is buing built.
#
# This mostly controls the behavior of the UID client/server, but
# may affect other components.  When adding a new site, modifications
# are needed in:
#   AS_RUN -- Makefile and runCA/runCAutil/site*
#   AS_ARD -- Makefile (experimental, only enabled for JCVI)
#   AS_UID -- Makefile
#
# Site-specific compilation settings can be made in c_make.as.
#
# The default site name is invalid to force the user to pick the
# correct site.  A reasonable default would be the LOCAL site, but
# this would generate numerous errors at JCVI when those users forget
# to set the proper site name.

ifeq ($(SITE_NAME),)

  #  The default site.
  #
  SITE_NAME=CELERA

  #  JCVI uses a web server via curl to get UIDs.
  #
  #SITE_NAME=JCVI

  #  The "LOCAL" site uses a text file in the current directory
  #  to remember what the last UID used is.
  #
  #SITE_NAME=LOCAL

  #  "BRI" is like "LOCAL", but adds some customization to runCA.
  #
  #SITE_NAME=LOCAL
endif

ifeq (,$(findstring $(SITE_NAME),"JCVI LOCAL BRI"))
  $(error Invalid SITE_NAME $(SITE_NAME).  Must be one of JCVI, LOCAL, BRI.  Please edit this file)
endif
