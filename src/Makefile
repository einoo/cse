## no implicit rules
.SUFFIXES : .f90

## definitions
FC=gfortran

NR = nr
OBS = obs
DA = enkf
CSE0 = cse_ms
CSE1 = fcst
CSE2 = csec
CSE3 = cse_org
CSET = cse_tt
CSER0 = cse_rev0

NAT_SRC = nature_run.f90 rk4.f90
NAT_OBJ = nature_run.o rk4.o

nature:  $(NR) clean

$(NR) : $(NAT_OBJ)
	$(FC) -o $(NR) $(NAT_OBJ)

$(NAT_OBJ):
	$(FC) -c $(NAT_SRC)

OBS_SRC = make_obs.f90 normal.f90
OBS_OBJ = make_obs.o normal.o

observation: $(OBS) clean

$(OBS) : $(OBS_OBJ)
	$(FC) -o $(OBS) $(OBS_OBJ) rk4.o

$(OBS_OBJ):
	$(FC) -c $(OBS_SRC) rk4.f90

DA_SRC = EnKF.f90 da_enkf.f90
DA_OBJ = EnKF.o da_enkf.o

dataassimilation : $(DA) clean

$(DA) : $(DA_OBJ)
	$(FC) -o $(DA) $(DA_OBJ) rk4.o normal.o -llapack

$(DA_OBJ):
	$(FC) -c $(DA_SRC) rk4.f90 normal.f90
	
CSE0_SRC = ETKF.f90 da_etkf.f90
CSE0_OBJ = ETKF.o da_etkf.o

original_cse : $(CSE0) clean

$(CSE0) : $(CSE0_OBJ)
	$(FC) -o $(CSE0) $(CSE0_OBJ) rk4.o -llapack

$(CSE0_OBJ):
	$(FC) -c $(CSE0_SRC) rk4.f90

CSE1_SRC = cse_fcst.f90
CSE1_OBJ = cse_fcst.o

free_forecast :$(CSE1) clean

$(CSE1) : $(CSE1_OBJ)
	$(FC) -o $(CSE1) $(CSE1_OBJ) rk4.o da_etkf.o -llapack

$(CSE1_OBJ):
	$(FC) -c $(CSE1_SRC) rk4.f90 da_etkf.f90

CSE2_SRC = cse_demo.f90
CSE2_OBJ = cse_demo.o

cse_case :$(CSE2) clean

$(CSE2) : $(CSE2_OBJ)
	$(FC) -o $(CSE2) $(CSE2_OBJ) rk4.o da_etkf.o -llapack

$(CSE2_OBJ):
	$(FC) -c $(CSE2_SRC) rk4.f90 da_etkf.f90

CSE3_SRC = cse_org.f90 cse_control.f90
CSE3_OBJ = cse_org.o cse_control.o

cse_run :$(CSE3) clean

$(CSE3) : $(CSE3_OBJ)
	$(FC) -o $(CSE3) $(CSE3_OBJ) rk4.o da_etkf.o -llapack

$(CSE3_OBJ):
	$(FC) -c $(CSE3_SRC) rk4.f90 da_etkf.f90

CSET_SRC = cse_tes.f90
CSET_OBJ = cse_tes.o

cse_tesv :$(CSET) clean

$(CSET) : $(CSET_OBJ)
	$(FC) -o $(CSET) $(CSET_OBJ) rk4.o -llapack

$(CSET_OBJ):
	$(FC) -c $(CSET_SRC) rk4.f90 

cse_rev0 : rev0_obj
	$(FC) -o $(CSER0) cse_rev0.o rk4.o cse_control.o da_etkf.o -llapack

rev0_obj:
	$(FC) -c cse_rev0.f90 rk4.f90 cse_control.f90 da_etkf.f90

clean:
	rm -f *.o *~ *.mod core

