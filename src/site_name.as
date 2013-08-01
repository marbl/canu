#
# $Id$
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

ifeq ($(SITE_NAME),)
  SITE_NAME=LOCAL
endif

ifneq ($(SITE_NAME),LOCAL)
 ifneq ($(SITE_NAME),JCVI)
  ifneq ($(SITE_NAME),BRI)
   $(error Invalid SITE_NAME $(SITE_NAME).  If set, must be JCVI)
  endif
 endif
endif

