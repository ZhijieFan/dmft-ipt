HERE  =`pwd`
EXE=ahmmpt_matsubara
DIR=./drivers/ahm
DIREXE= $(HOME)/.bin

#########################################################################
include $(HOME)/lib/lib.mk
include $(HOME)/lib/libdmft.mk
#########################################################################


all: 	version
	@echo " ........... compile: optimized ........... "
	@echo $(VER)
	$(FC) $(STD) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE) $(LIBDMFT) $(MODS) $(LIBS) 
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)

opt: 	version
	@echo " ........... compile: optimized ........... "
	$(FC) $(OPT) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE) $(LIBDMFT) $(MODS) $(LIBS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)


debug: 	version
	@echo " ........... compile : debug   ........... "
	$(FC) $(DEB) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE) $(LIBDMFT_DEB) $(MODS_DEB) $(LIBS_DEB) 
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)


clean: 
	@echo 'removing *.mod *.o *~'
	@rm -vf *.mod *.o *~ revision.inc

#########################################################################
include $(HOME)/lib/version.mk
#########################################################################
