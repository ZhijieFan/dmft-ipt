##HUBBARD MODELS 
#EXE=ipt_hm_matsubara
EXE=ipt_hm2d_matsubara
#EXE=mpt_hm_matsubara
#EXE=ipt_hm_real
#EXE=ipt_hm_keldysh
##ATTRACTIVE HUBBARD
#EXE=ipt_ahm_matsubara
#EXE=ipt_ahm_real
#EXE=ipt_ahm_keldysh
#EXE=ipt_ahm_keldysh_bias

DIR=./drivers
DIREXE= $(HOME)/.bin

.SUFFIXES: .f90

#REVISION SOFTWARE GIT:
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

OBJS=IPT_GF.o IPT_VARS_GLOBAL.o IPT_MATSUBARA.o IPT_REAL.o IPT_KELDYSH.o DMFT_IPT.o #IPT_AF_MATS.o

#MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

#INCARGS=-I/opt/scifor/gnu/include -I/opt/dmft_tools/gnu/include
#FFLAG +=-ffree-line-length-none $(INCARGS)

#ARGS=-ldmftt -lscifor $(MKLARGS) -lminpack -larpack -lparpack 
#ARGS= -ldmftt -lscifor -lfftpack -llapack -lblas -lminpack
ARGS= -ldmftt -lscifor -lfftpack -llapack -lblas -lminpack

all:compile


compile: version $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FFLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)

.f90.o:	
	$(FC) $(FFLAG) -c $< 

completion:
	src_completion.sh $(DIR)/$(EXE).f90
	@echo "run: . .bash_completion.d/$(EXE) to add completion for $(EXE) in this shell"

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

version:
	@echo $(VER)
