
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := $(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := $(OSTYPE)-$(MACHINETYPE)/bin
endif

SUBMAKEFILES := libutil/main.mk \
                libbio/main.mk \
                libseq/main.mk \
                libmeryl/main.mk \
                libkmer/main.mk \
                libsim4/main.mk \
                meryl/main.mk \
                leaff/main.mk \
                sim4dbutils/main.mk
