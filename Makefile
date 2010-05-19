ifeq (${QUIP_ROOT},)
  QUIP_ROOT = ${PWD}/../..
  export QUIP_ROOT
endif

all:
	$(MAKE) -f Makefile.atomeye -I${QUIP_ROOT} -I${QUIP_ROOT}/Makefiles -I${QUIP_ROOT}/build.${QUIP_ARCH} all

.DEFAULT:
	$(MAKE) -f Makefile.atomeye -I${QUIP_ROOT} -I${QUIP_ROOT}/Makefiles -I${QUIP_ROOT}/build.${QUIP_ARCH} $@
