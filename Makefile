#########################################################################
include sfmake.inc
#########################################################################
#EXE=ahmbcs_keldysh
#EXE=ahmipt_keldysh
#EXE=ahmipt_real
#EXE=pammpt_real_fixdens
EXE=ahmipt_matsubara
#EXE=ahmmpt_matsubara_2dsquare
#EXE=ahmmpt_real_2dsquare
#EXE=hmipt_keldysh_2dsquare
#EXE=hmipt_real
#EXE=hmipt_matsubara
#EXE=hmipt_matsubara_hypercubic
#EXE=hmipt_matsubara_dmft_loop_g0
#EXE=hmipt_matsubara_2dsquare
#EXE=hmmpt_matsubara

DIR=./drivers
DIREXE= $(HOME)/.bin

.SUFFIXES: .f90 
OBJS=IPT_VARS_GLOBAL.o \
IPT_MATS.o    \
IPT_KELDYSH.o \
IPT_SOPT.o    \
IPT_SC_SOPT.o \
IPT_AF_MATS.o \
DMFT_IPT.o

ARGS= $(SFLIBS)
ARGS_DEB=$(SFLIBS_DEB)
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)

#=================STANDARD COMPILATION====================================
all: FLAG=$(STD)
all: version $(OBJS)
	@echo " ........... compile: optimized ........... "
	@echo $(VER)
	$(FC) $(FLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)

#================OPTIMIZED COMPILATION====================================
opt: FLAG=$(OPT)
opt: 	version $(OBJS)
	@echo " ........... compile: optimized ........... "
	@echo $(VER)
	$(FC) $(FLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)


#================DEBUGGIN COMPILATION=====================================
debug: FLAG=$(DEB)
debug: 	version $(OBJS)
	@echo " ........... compile : debug   ........... "
	$(FC) $(FLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS_DEB)
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)


.f90.o:	
	$(FC) $(FLAG) -c $< $(SFINCLUDE) 

clean: 
	@echo 'removing *.mod *.o *~'
	@rm -vf *.mod *.o *~ revision.inc

version:
	@echo $(VER)
