ifeq (${QUIP_ROOT},)
  QUIP_ROOT = ${PWD}/../..
  export QUIP_ROOT
endif

ifneq (${QUIP_ARCH},)
  HAVE_QUIP = 1
  export HAVE_QUIP
  FOX = FoX-4.0.3
  export FOX_LIBDIR=${QUIP_ROOT}/src/FoX-4.0.3/objs.${QUIP_ARCH}/lib
  export FOX_INCDIR=${QUIP_ROOT}/src/FoX-4.0.3/objs.${QUIP_ARCH}/finclude
endif

all:
ifeq (${HAVE_QUIP},1)
	$(MAKE) -f Makefile.atomeye -I${QUIP_ROOT} -I${QUIP_ROOT}/arch -I${QUIP_ROOT}/build/${QUIP_ARCH}
else
	$(MAKE) -f Makefile.atomeye
endif

.DEFAULT: 
ifeq (${HAVE_QUIP},1)
	$(MAKE) -f Makefile.atomeye -I${QUIP_ROOT} -I${QUIP_ROOT}/arch -I${QUIP_ROOT}/build/${QUIP_ARCH} $@
else
	$(MAKE) -f Makefile.atomeye $@
endif

