#
# Makefile --
#
# Top-level makefile for the PRMS
#

include ./makelist

#
# Standard Targets for Users
#

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
