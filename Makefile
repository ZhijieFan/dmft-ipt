HERE  =`pwd`
EXE=hmmpt_af_matsubara
DIR=./drivers/hm
DIREXE= $(HOME)/.bin

#########################################################################
include $(SFDIR)/etc/lib.mk
include $(SFDIR)/etc/libdmft.mk
#########################################################################


all: 	version
	@echo " ........... compile: optimized ........... "
	@echo $(VER)
	$(FC) $(STD) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE) $(LIBDMFT) $(SFMODS) $(SFLIBS) 
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)

opt: 	version
	@echo " ........... compile: optimized ........... "
	$(FC) $(OPT) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE) $(LIBDMFT) $(SFMODS) $(SFLIBS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)


debug: 	version
	@echo " ........... compile : debug   ........... "
	$(FC) $(DEB) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE) $(LIBDMFT_DEB) $(SFMODS_DEB) $(SFLIBS_DEB) 
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)


clean: 
	@echo 'removing *.mod *.o *~'
	@rm -vf *.mod *.o *~ revision.inc

#########################################################################
include $(SFDIR)/etc/version.mk
#########################################################################
