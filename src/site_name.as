#
# $Id: site_name.as,v 1.6 2006-11-15 13:27:18 eliv Exp $
#
# Set the name of the site/location where the assembler is buing built.
#
# This mostly controls the behavior of the UID client/server, but
# may affect other components.
#
# Other site-specific settings can be made in c_make.as.
#
# The default site name is invalid to force the user to pick the correct
# site.  A reasonable default would be the LOCAL site, but this would
# generate numerous errors at TIGR and JCVI when those users forget to
# set the proper site name.

ifeq ($(SITE_NAME),)

  #  The default site.
  #
  SITE_NAME=CELERA

  #  TIGR uses a SOAP server to get UIDs.
  #
  #SITE_NAME=TIGR

  #  JCVI uses a web server via curl to get UIDs.
  #
  #SITE_NAME=JCVI

  #  The "LOCAL" site uses a text file in the current directory
  #  to remember what the last UID used is.
  #
  #SITE_NAME=LOCAL
endif

ifeq (,$(findstring $(SITE_NAME),"TIGR JCVI LOCAL"))
  $(error Invalid SITE_NAME $(SITE_NAME).  Must be one of TIGR, JCVI, LOCAL.  Please edit this file)
endif
