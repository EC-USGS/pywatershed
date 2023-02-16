#
# Makefile --
#
# Top-level makefile for the PRMS
#

# Arguments
#   FC: name of fortran compiler
#   CC: name of C compiler
#   DBL_PREC: [True|False] if doudble precision build is requested
#
# Example:
#   make FC=ifort CC=icc DBL_PREC=true

# Choose if all floats should be promoted to doubles
ifeq (, $(DBL_PREC))
     include ./makelist
else ifeq ($(DBL_PREC), false)
     include ./makelist
else ifeq ($(DBL_PREC), true)
     include ./makelist_double_precision
else
     $(error      Argument DBL_PREC must be "true" or "false")
endif

# Standard Targets for Users

all: directories standard

directories:
	@mkdir -p $(BINDIR)

standard:
	cd $(MMFDIR); $(MAKE);
	cd $(PRMSDIR); $(MAKE);

clean:
	cd $(MMFDIR); $(MAKE) clean;
	cd $(PRMSDIR); $(MAKE) clean;

cleanall:
	cd $(MMFDIR); $(MAKE) clean;
	cd $(PRMSDIR); $(MAKE) clean;
	cd $(PRMSDIR); $(MAKE) cleantarget;
