! Cell cycle from Basse2002

module cycle_mod
use real_kind_mod
use global

implicit none

integer, parameter :: G1_phase      = 1
integer, parameter :: G1_checkpoint = 2
integer, parameter :: S_phase       = 3
integer, parameter :: S_checkpoint  = 4
integer, parameter :: G2_phase      = 5
integer, parameter :: G2_checkpoint = 6
integer, parameter :: M_phase       = 7
integer, parameter :: dividing      = 8

logical :: use_volume_based_transition = .false.
real(REAL_KIND) :: starvation_arrest_threshold = 5
real(REAL_KIND) :: max_arrest_time = 6*3600
logical :: inhibit_misrepair = .false.
logical :: use_rad_state = .true.   ! irradiation before G2 causes G2 checkpoint delay = dose*1 hours

contains

!--------------------------------------------------------------------------
! Phase transitions are now based on cell volume, cp%V, to allow for delay
! when growth is slowed by starvation of oxygen and/or glucose.
! Note that the volumes required for the transitions (cp%G1_V,..)  never change.
! Now treat G1, S, G2 in the same way regarding checkpoints - no fixed delay
!--------------------------------------------------------------------------
subroutine log_timestep(cp, ccp, dt, dies)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: dt
logical :: dies
integer :: phase, ityp, nPL, kpar=0
real(REAL_KIND) :: cf_O2, cf_glucose, pcp_O2, pcp_glucose, pcp_starvation, R, fV
logical :: switch, S_switch, M_switch

S_switch = .false.
M_switch = .false.
if (cp%dVdt == 0) then
	write(nflog,*) 'dVdt=0, kcell: ',kcell_now,cp%phase
	stop
endif
nPL = cp%N_PL
!if (kcell_now == 2 .and. nPL > 0) write(*,*) 'kcell,phase,nPL: ',kcell_now,cp%phase,nPL
ityp = cp%celltype
10 continue
phase = cp%phase
!write(*,'(a,2i6,f8.0)') 'kcell, phase, phase: ',kcell_now, phase,tnow 
if (phase == G1_phase) then
    switch = (tnow > cp%G1_time)
    if (switch) then
        cp%phase = G1_checkpoint
        cp%G1_flag = .true.     ! no fixed checkpoint delay
        cp%G1S_time = tnow + f_TCP(ccp,nPL)
        goto 10
    endif
elseif (phase == G1_checkpoint) then  ! this checkpoint combines the release from G1 delay and the G1S repair check
!    if (.not.cp%G1_flag) then
!        R = par_uni(kpar)
!        cp%G1_flag = (R < ccp%Pk_G1*dt)
!    endif
    cp%G1S_flag = (nPL == 0 .or. tnow > cp%G1S_time)
    if (use_metabolism) then
		cp%G1S_flag = cp%G1S_flag .and. (cp%metab%A_rate > r_Ag)
	endif
    if (cp%G1_flag .and. cp%G1S_flag) then  ! switch to S-phase
        S_switch = .true.
        cp%phase = S_phase
! Note: now %I_rate has been converted into equivalent %dVdt, to simplify code 
	    cp%S_duration = (max_growthrate(ityp)/cp%dVdt)*cp%fg*ccp%T_S
	    cp%S_time = 0   ! this is now the amount of time progress through S phase: 0 -> %S_duration
        goto 10
    endif
elseif (phase == S_phase) then
    cp%arrested = (cp%dVdt/max_growthrate(ityp) < ccp%arrest_threshold)
    if (.not.cp%arrested) then
        cp%S_time = cp%S_time + dt
    endif
    switch = (cp%S_time >= cp%S_duration)
    if (switch) then
!        if (kcell_now == 20) write(nflog,*) 'switch S - chkpt: tnow: ',tnow/3600
        cp%phase = S_checkpoint
        cp%S_flag = .true. ! no fixed checkpoint delay
        cp%SG2_time = tnow + f_TCP(ccp,nPL)
        goto 10
    endif
elseif (phase == S_checkpoint) then
    cp%SG2_flag = (nPL == 0 .or. tnow > cp%SG2_time)
    if (use_metabolism) then
		cp%SG2_flag = cp%SG2_flag .and. (cp%metab%A_rate > r_Ag)
	endif
    if (cp%S_flag .and. cp%SG2_flag) then
!        if (kcell_now == 20) write(nflog,*) 'switch chkpt - G2: tnow: ',tnow/3600
        cp%phase = G2_phase
! Note: now %I_rate has been converted into equivalent %dVdt, to simplify code 
		cp%G2_time = tnow + (max_growthrate(ityp)/cp%dVdt)*cp%fg*ccp%T_G2
        goto 10
    endif
elseif (phase == G2_phase) then
!    if (tnow > cp%G2_time .and. cp%V < cp%divide_volume) then
!        fV = (cp%divide_volume-cp%V)/cp%divide_volume
!        if (fV > 0.012) write(nflog,'(a,i4,4e12.3)') 'G2 exit: V, divide_volume: ',kcell_now,tnow,cp%V,cp%divide_volume,fV
!    endif
	switch = (tnow > cp%G2_time .and. cp%V > cp%divide_volume) ! to prevent volumes decreasing 
    if (switch) then
!        if (kcell_now == 20) write(nflog,*) 'switch G2 - chkpt: tnow: ',tnow/3600
!        fV = (cp%V - cp%divide_volume)/cp%divide_volume
        cp%V = cp%divide_volume     ! correct for slight volume discrepancy here, to maintain correct cell volume
        cp%phase = G2_checkpoint
        cp%G2_flag = .false.
        cp%G2M_time = tnow + f_TCP(ccp,nPL)		!ccp%Tcp(nPL)
        if (use_rad_state .and. cp%rad_state > 0 .and. cp%rad_state < 3) then
            cp%G2M_time = tnow + rad_dose*ccp%G2_delay_factor*3600    ! default 1h/Gy
        endif
        goto 10
    endif
elseif (phase == G2_checkpoint) then ! this checkpoint combines the release from G2 delay and the G2M repair check
!    if (.not.cp%G2_flag) then
!        R = par_uni(kpar)
!        cp%G2_flag = (R < ccp%Pk_G2*dt)
!    endif
    cp%G2_flag = .true.     ! no checkpoint delay
    cp%G2M_flag = (nPL == 0 .or. tnow > cp%G2M_time)
    if (use_metabolism) then
		cp%G2M_flag = cp%G2M_flag .and. (cp%metab%A_rate > r_Ag)
	endif
    if (cp%G2_flag .and. cp%G2M_flag) then  ! switch to M-phase
        M_switch = .true.
        cp%phase = M_phase
        cp%M_time = tnow + ccp%T_M   
        if (cp%rad_state > 0) then     ! count cells, count DSBs, count Ch2 lesions
            rad_count(1) = rad_count(1) + 1
            rad_count(2) = rad_count(2) + cp%N_PL
            rad_count(3) = rad_count(3) + cp%N_Ch2
        endif
        goto 10
    endif
elseif (phase == M_phase) then
!    if (tnow > cp%M_time) then
!        cp%phase = dividing
!        cp%doubling_time = tnow
!    endif
    cp%phase = M_phase      ! signals start of mitosis
endif
dies = .false.
if (S_switch .and. (cp%N_PL > 0 .or. cp%N_IRL > 0)) then
    dies = mortality(cp,ccp,'S')
endif
if (M_switch .and. (cp%N_PL > 0 .or. cp%N_IRL > 0 .or. cp%N_Ch1 > 0 .or. cp%N_Ch2 > 0)) then
    dies = mortality(cp,ccp,'M')
    if (dies .and. cp%rad_state > 0) rad_count(6) = rad_count(6) + 1
endif
if (dies) return
!if (nPL > 0 .and. .not.cp%irrepairable) then
if (nPL > 0) then
    call radiation_repair(cp, ccp, dt)
endif
end subroutine

!--------------------------------------------------------------------------
! This uses exponentially distributed checkpoint times for G1, S, G2,
! fixed (with growth scaling) G1, S, G2 times.
!--------------------------------------------------------------------------
subroutine exp_timestep(cp, ccp, dt, dies)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: dt
logical :: dies
integer :: phase, ityp, nPL, kpar=0
real(REAL_KIND) :: delay, duration
logical :: switch, S_switch, M_switch

phase = cp%phase
S_switch = .false.
M_switch = .false.
!if (colony_simulation) write(*,*) 'colony_simulation: ID,phase,N_PL: ',cp%ID,phase,cp%N_PL,cp%N_Ch1,cp%N_Ch2
if (cp%dVdt == 0) then
!	if (phase == G1_phase .or. phase == S_phase .or. phase == G2_phase) then
		write(nflog,*) 'dVdt=0, kcell, phase: ',kcell_now,phase
		stop
!	endif
endif
nPL = cp%N_PL
ityp = cp%celltype

if (phase == G1_phase) then
    if (tnow > cp%G1_time) then
        cp%phase = G1_checkpoint
        cp%G1ch_entry_time = tnow
        cp%G1ch_max_delay = max(cp%G1ch_time, f_TCP(ccp,nPL))
    endif
elseif (phase == G1_checkpoint) then  ! this checkpoint combines the release from G1 delay and the G1S repair check
    if (nPL == 0) then
        delay = cp%G1ch_time
    else
        delay = cp%G1ch_max_delay
    endif
    if (tnow > cp%G1ch_entry_time + delay) then
        S_switch = .true.
        cp%phase = S_phase
	    duration = (max_growthrate(ityp)/cp%dVdt)*ccp%T_S
        cp%S_time = tnow + duration
    endif
elseif (phase == S_phase) then
    if (tnow > cp%S_time) then
        cp%phase = S_checkpoint
        cp%Sch_entry_time = tnow
        cp%Sch_max_delay = max(cp%Sch_time, f_TCP(ccp,nPL))
    endif
elseif (phase == S_checkpoint) then
    if (nPL == 0) then
        delay = cp%Sch_time
    else
        delay = cp%Sch_max_delay
    endif
    if (tnow > cp%Sch_entry_time + delay) then
        cp%phase = G2_phase
	    duration = (max_growthrate(ityp)/cp%dVdt)*ccp%T_G2
        cp%G2_time = tnow + duration
    endif
elseif (phase == G2_phase) then
    if (tnow > cp%G2_time .and. cp%V < cp%divide_volume) then
        write(*,'(a,3e12.3)') 'G2 exit: V, divide_volume: ',cp%V,cp%divide_volume,(cp%divide_volume-cp%V)/cp%divide_volume
    endif
	switch = (tnow > cp%G2_time .and. cp%V > cp%divide_volume) ! try this to prevent volumes decreasing 
    if (switch) then
        cp%phase = G2_checkpoint
        cp%G2ch_entry_time = tnow
        cp%G2ch_max_delay = max(cp%G2ch_time, f_TCP(ccp,nPL))
    endif
elseif (phase == G2_checkpoint) then ! this checkpoint combines the release from G2 delay and the G2M repair check
    if (nPL == 0) then
        delay = cp%G2ch_time
    else
        delay = cp%G2ch_max_delay
    endif
    if (tnow > cp%G2ch_entry_time + delay) then
        M_switch = .true.
        cp%phase = M_phase
        cp%M_time = tnow + ccp%T_M   
    endif
elseif (phase == M_phase) then
    if (tnow > cp%M_time) then
        cp%phase = dividing
    endif
endif
!if (cp%ID == 1 .and. (nPL > 0 .or. cp%N_Ch1 > 0)) then
!    write(nflog,*) 'ID=1: phase,nPL: ',cp%phase,nPL,cp%N_Ch1,cp%N_Ch2,M_switch
!endif
dies = .false.  
if (S_switch .and. (cp%N_PL > 0 .or. cp%N_IRL > 0)) then
    dies = mortality(cp,ccp,'S')
endif
if (M_switch .and. (cp%N_PL > 0 .or. cp%N_IRL > 0 .or. cp%N_Ch1 > 0 .or. cp%N_Ch2 > 0)) then
    dies = mortality(cp,ccp,'M')
    if (cp%ID == 1) write(nflog,*) 'ID=1: M_switch: dies: ',dies
endif
if (dies) return
!if (nPL > 0 .and. .not.cp%irrepairable) then
if (nPL > 0) then
    call radiation_repair(cp, ccp, dt)
endif
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function mortality(cp,ccp,pchar) result(dies)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
character :: pchar
integer :: N, ityp, kpar=0
real(REAL_KIND) :: R, psurvive, psurvive_IRL, pdeath    ! , psurvive_PL
logical :: dies

N = cp%N_PL + cp%N_IRL
!psurvive_PL = 0.999  !ccp%psurvive_Ch1
psurvive_IRL = ccp%psurvive_PL
psurvive = (ccp%psurvive_PL**cp%N_PL)*(psurvive_IRL**cp%N_IRL)
if (pchar == 'M') then
!    psurvive = (psurvive/10)*(ccp%psurvive_Ch1**cp%N_Ch1)*(ccp%psurvive_Ch2**cp%N_Ch2) ! why /10?????
    psurvive = (psurvive)*(ccp%psurvive_Ch1**cp%N_Ch1)*(ccp%psurvive_Ch2**cp%N_Ch2)
endif
pdeath = 1 - psurvive
R = par_uni(kpar)
dies = (R < pdeath)
!write(*,'(a,a,2i4,f6.3,L2)') 'mortality: pchar,psurvive: ',pchar,cp%N_PL,cp%N_IRL,psurvive,dies
end function

!--------------------------------------------------------------------------
! We want:
!    Prel = 1 for Fc >= Fmax
!    Prel = 0 for Fc <= 1
!    Prel to vary between 0 and 1 as Fc varies between 1 and Fmax
! For intermediate Fc values, make the mean duration of arrest an
! exponentially-distributed random variable with mean Tarr, then in a time
! step of length dt the probability of release is Prel = dt/Tarr.
! For Prel to vary between 0 and 1 we would need Tarr to vary between 
! infinity and dt.  Since infinity is not an option, instead we can specify
! a maximum arrest time, e.g. Tarr = Tmax (= 10h, say).  This means that as 
! Fc approaches 1, the minimum value of Prel = dt/Tmax, which for dt = 10 min 
! and Tmax = 10h gives Prel = 1/60.
!--------------------------------------------------------------------------
function getPcp_release(fc,dt) result(pcp)
real(REAL_KIND) :: fc, dt, pcp
real(REAL_KIND) :: mu

if (fc <= 1) then
    pcp = 0
elseif (fc >= starvation_arrest_threshold) then
    pcp = 1
else
    mu = dt + (max_arrest_time - dt)*(starvation_arrest_threshold - fc)/(starvation_arrest_threshold - 1)
    pcp = dt/mu
endif
end function


!--------------------------------------------------------------------------
! Damage from radiation dose following Curtis1986
! dose is the Gy received, tmin is the duration (min) over which it is delivered.
! Krepair, Kmisrepair are Curtis's epsilon_PL, epsilon_2PL
! As modified according to Bill Wilson:
! There are three classes of lesions:
!	%NL1 = number of potentially lethal (i.e. repairable) lesions = Curtis's n_PL
!	%NL2(1) = number of lethal lesions of type b
!	%NL2(2) = number of lethal lesions of type c
! During the radiation exposure, damage of different kinds occurs at a rate
! that is given by the product dose intensity (Gy/h), eta, and SER_OER.
! eta_PL and eta_L(:) are the rate constants for PL and L lesions
! SER_OER is the product of the SER and OER, where
! SER = Sensitivity Enhancement Ratio which is drug dependent
! OER = Oxygen Enhancement Ratio which is a function of intracellular O2 conc.
! (Note: OER for PL uses OER_alpha, for L uses OER_beta)
! Over the duration of the exposure repair is also occurring (this will be
! insignificant if the duration is very short).
! Note that the dose (Gy) occurs over the time tmin (minutes).  The implicit 
! assumption is that tmin < DELTA_T.
! Since repair is simulated in the following call to timestep() we can
! ignore repair in the (assumed short) dose period.
! Now eta_PL doubles during S-phase
!--------------------------------------------------------------------------
subroutine radiation_damage(cp, ccp, dose0, SER_OER, tmin)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: dose0, SER_OER, tmin
real(REAL_KIND) :: dose, dDdt, dtmin, dthour, R
real(REAL_KIND) :: p_PL, p_IRL, Krepair, Kmisrepair, misrepair_factor, fraction, inhibition, eta_PL
integer :: nt, it, nPL, nPL0, nIRL, ityp, kpar=0
logical :: do_repair = .false.

if (cp%phase < S_phase) then
    eta_PL = ccp%eta_PL
elseif (cp%phase > S_phase) then
    eta_PL = 2*ccp%eta_PL
else
    fraction = cp%S_time/cp%S_duration
	fraction = max(0.0, fraction)
	fraction = min(1.0, fraction)
    eta_PL = (1 + fraction)*ccp%eta_PL
endif

dose = dose0*SER_OER
nt = dose*eta_PL/0.01
dtmin = tmin/nt
!write(*,*) 'from dose: nt, dtmin: ',dose,nt,dtmin 
!dtmin = 0.0001
!nt = max(tmin/dtmin, 1.0)
dthour = dtmin/60
!nt = max(tmin/dtmin, 1.0)
dDdt = dose/(nt*dthour)     ! Gy/h
nPL0 = cp%N_PL
if (cp%irrepairable) then
	write(nflog,*) 'radiation_damage: irrepairable cell!'
	return	! this should not happen
endif
nIRL = 0
if (do_repair) then
    call getRepairParameters(cp,ccp,Krepair,Kmisrepair)
endif

if (.not.do_repair) then
	nPL = dose*eta_PL
	p_PL =  dose*eta_PL - nPL
	R = par_uni(kpar)
	if (R < p_PL) then
		nPL = nPL + 1
	endif
	nIRL = dose*ccp%eta_IRL
	p_IRL = dose*ccp%eta_IRL - nIRL
	R = par_uni(kpar)
	if (R < p_IRL) then
		nIRL = nIRL + 1
	endif
	
	cp%N_PL = nPL0 + nPL
	cp%N_IRL = nIRL
	cp%irrepairable = (nIRL > 0)
		
	return
endif

inhibition = 1
if (use_inhibiter) then
    C_inhibiter = cp%Cin(drug_A)
    inhibition = repairInhibition(C_inhibiter)
endif
Krepair = (1 - inhibition)*Krepair
if (inhibit_misrepair) then
    Kmisrepair = (1 - inhibition)*Kmisrepair
endif

p_PL = eta_PL*dose/nt
p_IRL = ccp%eta_IRL*dose/nt
do it = 1,nt
	R = par_uni(kpar)
	if (R < p_PL + p_IRL) then
		if (R < p_PL) then
			nPL = nPL + 1
		else
			nIRL = nIRL + 1
		endif
	endif
	if (do_repair) then
		! repair/misrepair
		R = par_uni(kpar)
		if (R < nPL*Krepair*dthour) then
			nPL = nPL - 1
		endif
		R = par_uni(kpar)
		if (R < (nPL**2)*Kmisrepair*dthour) then	! -> scalar
			nPL = nPL - 1
			R = par_uni(kpar)
			if (R < ccp%fraction_Ch1) then
				cp%N_Ch1 = cp%N_Ch1 + 1
			else
				cp%N_Ch2 = cp%N_Ch2 + 1
			endif
		endif
	endif
enddo
cp%N_PL = nPL
cp%N_IRL = nIRL
cp%irrepairable = (nIRL > 0)
end subroutine

!--------------------------------------------------------------------------
! Time unit = hour
! This may need to be changed, because it implicitly assumes that no more
! than one repair and one misrepair of each type can occur within a time step.
! The fix would be to subdivide the time step, as in the damage subroutine.
! The alternate computation of the number of misrepairs is based on integrating
! the rate of misrepair, which is dN/dt = -N^2.Kmisrepair.dt
! dN/N^2 = - K.dt
! -1/N = -K.t + c
! -1/N0 = K.0 + c ==> c = -1/N0
! N(t) = N(0)/(Kmisrepair*N(0)*t + 1)
!--------------------------------------------------------------------------
subroutine radiation_repair(cp, ccp, dt)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: dt
integer :: i, k, ityp, nPL, nPL0, nmis, kpar=0
real(REAL_KIND) :: dthour, fraction, Krepair, Kmisrepair, misrepair_factor, R, p_rep, p_mis
real(REAL_KIND) :: rnPL, rnPL0, dPL, inhibition
integer :: nt, it
logical :: use_prob = .false.

nPL = cp%N_PL
if (nPL == 0) return
if (use_prob) then
	nt = 100
else
	nt = 1
endif
dthour = dt/(nt*3600)    ! sec -> hour
call getRepairParameters(cp,ccp,Krepair,Kmisrepair)

! First allow true repair to occur
rnPL = nPL*exp(-Krepair*nt*dthour)
nPL = rnPL
R = par_uni(kpar)
if (R < (rnPL - nPL)) nPL = nPL + 1
!write(nflog,*) 'true repair: ', nPL0,nPL,nPL0-nPL

! Then misrepair occurs on remaining nPL
if (use_prob) then
	do it = 1,nt
		if (nPL == 0) exit
		p_mis = nPL**2*Kmisrepair*dthour
		R = par_uni(kpar)
		if (R < p_mis) then
			nPL = nPL - 1
			R = par_uni(kpar)
			if (R < ccp%fraction_Ch1) then
				cp%N_Ch1 = cp%N_Ch1 + 1
			else
				cp%N_Ch2 = cp%N_Ch2 + 1
			endif
		endif
	enddo
else
	if (nPL > 0) then
		rnPL0 = nPL
		rnPL = rnPL0/(rnPL0*Kmisrepair*nt*dthour + 1)
		dPL = rnPL0 - rnPL		! this is the number of misrepairs - convert to integer
		nmis = dPL
!		write(*,*) 'misrepair: ',dPL,nmis
		R = par_uni(kpar)
		if (R < (dPL - nmis)) nmis = nmis + 1
		nPL = nPL - nmis
		do k = 1,nmis
			R = par_uni(kpar)
			if (R < ccp%fraction_Ch1) then
				cp%N_Ch1 = cp%N_Ch1 + 1
			else
				cp%N_Ch2 = cp%N_Ch2 + 1
			endif
		enddo
	endif
endif
cp%N_PL = nPL
!write(*,*) 'stopping in radiation_repair'
!stop
!if (cp%ID == 1) then
!    write(nflog,'(a,4i6)') 'repair: cell #, N_PL, N_Ch2, N_Ch2: ',cp%ID,cp%N_PL,cp%N_CH1,cp%N_Ch2
!    write(nflog,'(a,4e12.3)') 'C, Krepair, inhibition, Kmisrepair: ',C_inhibiter,Krepair, inhibition, Kmisrepair
!endif
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine getRepairParameters(cp,ccp,Krepair,Kmisrepair)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: Krepair, Kmisrepair
real(REAL_KIND) :: Krepair_HRR, Krepair_NHEJ
real(REAL_KIND) :: Kmisrepair_NHEJ, Kmisrepair_DIM
real(REAL_KIND) :: fraction, inhibition, C_inhibiter
integer :: ityp

ityp = cp%celltype
if (cp%phase < S_phase) then
    Krepair_HRR = ccp%HRR_repair_base
elseif (cp%phase > S_phase) then
    Krepair_HRR = ccp%HRR_repair_max
else
    if (use_volume_based_transition) then   ! fraction = fraction of passage through S phase
        fraction = (cp%V - cp%G1_V)/(cp%S_V - cp%G1_V)
    else
!        fraction = 1 - (cp%S_time - tnow)/(cp%S_time - cp%S_start_time)
        fraction = cp%S_time/cp%S_duration
		fraction = max(0.0, fraction)
		fraction = min(1.0, fraction)
    endif
    Krepair_HRR = ccp%HRR_repair_base + fraction*(ccp%HRR_repair_max - ccp%HRR_repair_base)
endif
Krepair_NHEJ = ccp%NHEJ_repair
Kmisrepair_NHEJ = ccp%NHEJ_misrepair
Kmisrepair_DIM = ccp%DIM_misrepair

inhibition = 0
if (use_inhibiter) then
    C_inhibiter = cp%Cin(drug_A)
    inhibition = repairInhibition(C_inhibiter)
endif
Krepair_NHEJ = (1 - inhibition)*Krepair_NHEJ
if (inhibit_misrepair) then
    Kmisrepair_NHEJ = (1 - inhibition)*Kmisrepair_NHEJ
endif
Krepair = Krepair_HRR + Krepair_NHEJ
Kmisrepair = Kmisrepair_NHEJ + Kmisrepair_DIM

end subroutine

!--------------------------------------------------------------------------
! Hill function for max checkpoint delay TCP
!--------------------------------------------------------------------------
function f_TCP(ccp,n) result(TCP)
integer :: n
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: TCP

TCP = ccp%bTCP*n/(ccp%aTCP + n)
!if (n > 0 .and. kcell_now < 100) write(*,'(a,2i8,f8.3)') 'f_TCP: ',kcell_now,n,TCP
TCP = 3600*TCP	! hours -> seconds
end function

!--------------------------------------------------------------------------
! C is the concentration of the repair inhibiting drug.
! More inhibition ==> less repair, because Krepair ==> (1 - inhibition)*Krepair
!--------------------------------------------------------------------------
function repairInhibition(C) result(inhibition)
real(REAL_KIND) :: C, inhibition

inhibition = a_inhibit*C/(b_inhibit + C)
end function

end module

