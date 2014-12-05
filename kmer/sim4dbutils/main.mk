
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif


SUBMAKEFILES := cleanPolishes.mk \
                comparePolishes.mk \
                fixPolishesIID.mk \
                convertToAtac.mk \
                convertToExtent.mk \
                convertPolishes.mk \
                detectChimera.mk \
                depthOfPolishes.mk \
                filterPolishes.mk \
                headPolishes.mk \
                mappedCoverage.mk \
                mergePolishes.mk \
                parseSNP.mk \
                pickBestPolish.mk \
                pickBestPair.mk \
                pickUniquePolish.mk \
                plotCoverageVsIdentity.mk \
                removeDuplicate.mk \
                sortPolishes.mk \
                summarizePolishes.mk \
                uniqPolishes.mk \
                vennPolishes.mk \
                realignPolishes.mk \
                removeRedundant.mk \
                reportAlignmentDifferences.mk
