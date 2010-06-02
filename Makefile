ifeq (${QUIP_ROOT},)
  QUIP_ROOT = ${PWD}/../..
  export QUIP_ROOT
endif

ifneq (${QUIP_ARCH},)
  HAVE_QUIP = 1
  export HAVE_QUIP
endif

all:
ifeq (${HAVE_QUIP},1)
	$(MAKE) -f Makefile.atomeye -I${QUIP_ROOT} -I${QUIP_ROOT}/Makefiles -I${QUIP_ROOT}/build.${QUIP_ARCH}
else
	$(MAKE) -f Makefile.atomeye
endif

.DEFAULT: 
ifeq (${HAVE_QUIP},1)
	$(MAKE) -f Makefile.atomeye -I${QUIP_ROOT} -I${QUIP_ROOT}/Makefiles -I${QUIP_ROOT}/build.${QUIP_ARCH} $@
else
	$(MAKE) -f Makefile.atomeye $@
endif

