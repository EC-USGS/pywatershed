#
# Makefile --
#
# Top-level makefile for the PRMS
#

include ./makelist

#
# Standard Targets for Users
#

all: standard

standard:
	cd $(MMFDIR); $(MAKE);
	cd $(PRMSDIR); $(MAKE);

clean:
	cd $(MMFDIR); $(MAKE) clean;
	cd $(PRMSDIR); $(MAKE) clean;

