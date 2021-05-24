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

contains

!--------------------------------------------------------------------------
! Phase transitions are now based on cell volume, cp%V, to allow for delay
! when growth is slowed by starvation of oxygen and/or glucose.
! Note that the volumes required for the transitions (cp%G1_V,..)  never change.
! Now treat G1, S, G2 in the same way regarding checkpoints
!--------------------------------------------------------------------------
subroutine log_timestep(cp, ccp, dt)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: dt
integer :: phase, ityp, nPL, kpar=0
real(REAL_KIND) :: cf_O2, cf_glucose, pcp_O2, pcp_glucose, pcp_starvation, R
logical :: switch

phase = cp%phase
if (cp%dVdt == 0) then
	write(nflog,*) 'dVdt=0, kcell, phase: ',kcell_now,phase
	stop
endif
nPL = cp%N_PL
ityp = cp%celltype
if (phase == G1_phase) then
    if (use_volume_based_transition) then
        switch = (cp%V > cp%G1_V)
    else
        switch = (tnow > cp%G1_time)
    endif
    if (switch) then
        cp%phase = G1_checkpoint
        cp%G1_flag = .false.
        cp%G1S_time = tnow + f_TCP(ccp,nPL)		!ccp%Tcp(nPL)
    endif
elseif (phase == G1_checkpoint) then  ! this checkpoint combines the release from G1 delay and the G1S repair check
    if (.not.cp%G1_flag) then
        R = par_uni(kpar)
        cp%G1_flag = (R < ccp%Pk_G1*dt)
    endif
    cp%G1S_flag = (nPL == 0 .or. tnow > cp%G1S_time)
    if (use_metabolism) then
		cp%G1S_flag = cp%G1S_flag .and. (cp%metab%A_rate > r_Ag)
	endif
    if (cp%G1_flag .and. cp%G1S_flag) then
        cp%phase = S_phase
! Note: now %I_rate has been converted into equivalent %dVdt, to simplify code 
	    cp%S_duration = (max_growthrate(ityp)/cp%dVdt)*ccp%T_S
	    cp%S_time = 0   ! this is now the amount of time progress through S phase: 0 -> %S_duration
    endif
elseif (phase == S_phase) then
    cp%arrested = (cp%dVdt/max_growthrate(ityp) < ccp%arrest_threshold)
    if (.not.cp%arrested) then
        cp%S_time = cp%S_time + dt
    endif
    switch = (cp%S_time >= cp%S_duration)
!   switch = (tnow > cp%S_time)
    if (switch) then
!        cp%phase = G2_phase
        cp%phase = S_checkpoint
        cp%S_flag = .false.
        cp%SG2_time = tnow + f_TCP(ccp,nPL)
    endif
elseif (phase == S_checkpoint) then
    if (.not.cp%S_flag) then
        R = par_uni(kpar)
        cp%S_flag = (R < ccp%Pk_S*dt)
    endif
    cp%SG2_flag = (nPL == 0 .or. tnow > cp%SG2_time)
    if (use_metabolism) then
		cp%SG2_flag = cp%SG2_flag .and. (cp%metab%A_rate > r_Ag)
	endif
    if (cp%S_flag .and. cp%SG2_flag) then
        cp%phase = G2_phase
! Note: now %I_rate has been converted into equivalent %dVdt, to simplify code 
!	    cp%S_duration = (max_growthrate(ityp)/cp%dVdt)*ccp%T_S
!	    cp%S_time = 0   ! this is now the amount of time progress through S phase: 0 -> %S_duration
		cp%G2_time = tnow + (max_growthrate(ityp)/cp%dVdt)*ccp%T_G2
    endif

elseif (phase == G2_phase) then
!    if (use_volume_based_transition) then
!        switch = (cp%V > cp%G2_V)
!    else
		switch = (tnow > cp%G2_time .and. cp%V > cp%divide_volume) ! try this to prevent volumes decreasing 
!    endif
    if (switch) then
        cp%phase = G2_checkpoint
        cp%G2_flag = .false.
        cp%G2M_time = tnow + f_TCP(ccp,nPL)		!ccp%Tcp(nPL)
    endif
elseif (phase == G2_checkpoint) then ! this checkpoint combines the release from G2 delay and the G2M repair check
    if (.not.cp%G2_flag) then
        R = par_uni(kpar)
        cp%G2_flag = (R < ccp%Pk_G2*dt)
    endif
    cp%G2M_flag = (nPL == 0 .or. tnow > cp%G2M_time)
    if (use_metabolism) then
		cp%G2M_flag = cp%G2M_flag .and. (cp%metab%A_rate > r_Ag)
	endif
    if (cp%G2_flag .and. cp%G2M_flag) then
        cp%phase = M_phase
        cp%M_time = tnow + ccp%T_M   
    endif
elseif (phase == M_phase) then
    if (tnow > cp%M_time) then
        cp%phase = dividing
!        cp%doubling_time = tnow
    endif
endif    
if (nPL > 0 .and. .not.cp%irrepairable) then
    call radiation_repair(cp, ccp, dt)
endif
end subroutine

!--------------------------------------------------------------------------
! This uses exponentially distributed checkpoint times for G1, S, G2,
! fixed (with growth scaling) G1, S, G2 times.
!--------------------------------------------------------------------------
subroutine exp_timestep(cp, ccp, dt)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: dt
integer :: phase, ityp, nPL, kpar=0
real(REAL_KIND) :: delay, duration
logical :: switch

phase = cp%phase
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
        cp%phase = M_phase
        cp%M_time = tnow + ccp%T_M   
    endif
elseif (phase == M_phase) then
    if (tnow > cp%M_time) then
        cp%phase = dividing
    endif
endif    
if (nPL > 0 .and. .not.cp%irrepairable) then
    call radiation_repair(cp, ccp, dt)
endif
end subroutine


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
!--------------------------------------------------------------------------
subroutine radiation_damage(cp, ccp, dose0, SER_OER, tmin)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: dose0, SER_OER, tmin
real(REAL_KIND) :: dose, dDdt, dtmin, dthour, R
real(REAL_KIND) :: p_PL, p_IRL, Krepair, Kmisrepair, misrepair_factor, fraction, inhibition
integer :: nt, it, nPL, nPL0, nIRL, ityp, kpar=0
logical :: do_repair = .false.

dose = dose0*SER_OER
nt = dose*ccp%eta_PL/0.01
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
	! For repair/misrepair
	ityp = cp%celltype
	if (cp%phase < S_phase) then
		Krepair = ccp%Krepair_base
	elseif (cp%phase > S_phase) then
		Krepair = ccp%Krepair_max
	else
		if (use_volume_based_transition) then   ! fraction = fraction of passage through S phase
			fraction = (cp%V - cp%G1_V)/(cp%S_V - cp%G1_V)
		else
!			fraction = 1 - (cp%S_time - tnow)/(cp%S_time - cp%S_start_time)
            fraction = cp%S_time/cp%S_duration
			fraction = max(0.0, fraction)
			fraction = min(1.0, fraction)
		endif
		Krepair = ccp%Krepair_base + fraction*(ccp%Krepair_max - ccp%Krepair_base)
	endif
	if (cp%phase == M_phase) then
		misrepair_factor = ccp%mitosis_factor
	else
		misrepair_factor = 1
	endif
	Kmisrepair = misrepair_factor*ccp%Kmisrepair
endif

if (.not.do_repair) then
	nPL = dose*ccp%eta_PL
	p_PL =  dose*ccp%eta_PL - nPL
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
	
    if (cp%ID == 1) then
        write(nflog,*) 'damage: cell #, nPL0, N_PL, N_IRL: ',cp%ID,nPL0,cp%N_PL,cp%N_IRL
        write(nflog,*) 'eta_PL, dose0, SER_OER, C: ',ccp%eta_PL,dose0,SER_OER,C_inhibiter
    endif
	
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

p_PL = ccp%eta_PL*dose/nt
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
integer :: i, k, ityp, nPL, nmis, kpar=0
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
ityp = cp%celltype
if (cp%phase < S_phase) then
    Krepair = ccp%Krepair_base
elseif (cp%phase > S_phase) then
    Krepair = ccp%Krepair_max
else
    if (use_volume_based_transition) then   ! fraction = fraction of passage through S phase
        fraction = (cp%V - cp%G1_V)/(cp%S_V - cp%G1_V)
    else
!        fraction = 1 - (cp%S_time - tnow)/(cp%S_time - cp%S_start_time)
        fraction = cp%S_time/cp%S_duration
		fraction = max(0.0, fraction)
		fraction = min(1.0, fraction)
    endif
    Krepair = ccp%Krepair_base + fraction*(ccp%Krepair_max - ccp%Krepair_base)
endif

inhibition = 1
if (use_inhibiter) then
    C_inhibiter = cp%Cin(drug_A)
    inhibition = repairInhibition(C_inhibiter)
endif
if (cp%ID == 1) write(nflog,'(a,f8.4)') 'inhibition: ',inhibition
Krepair = (1 - inhibition)*Krepair
if (cp%phase == M_phase) then
	misrepair_factor = ccp%mitosis_factor
else
	misrepair_factor = 1
endif
Kmisrepair = misrepair_factor*ccp%Kmisrepair	! -> scalar
if (inhibit_misrepair) then
    Kmisrepair = (1 - inhibition)*Kmisrepair
endif

!do it = 1,nt
!	R = par_uni(kpar)
!	if (R < nPL*Krepair*dthour) then
!		nPL = nPL - 1
!		if (nPL == 0) exit
!	endif
!	R = par_uni(kpar)
!	if (R < nPL**2*Kmisrepair*dthour) then
!		nPL = nPL - 1
!		R = par_uni(kpar)
!		if (R < ccp%fraction_Ch1) then
!			cp%N_Ch1 = cp%N_Ch1 + 1
!		else
!			cp%N_Ch2 = cp%N_Ch2 + 1
!		endif
!		if (nPL == 0) exit
!	endif
!enddo

! First allow true repair to occur
rnPL = nPL*exp(-Krepair*nt*dthour)
nPL = rnPL
R = par_uni(kpar)
if (R < (rnPL - nPL)) nPL = nPL + 1

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
!if (cp%ID == 1) then
!    write(nflog,'(a,4i6)') 'repair: cell #, N_PL, N_Ch2, N_Ch2: ',cp%ID,cp%N_PL,cp%N_CH1,cp%N_Ch2
!    write(nflog,'(a,4e12.3)') 'C, Krepair, inhibition, Kmisrepair: ',C_inhibiter,Krepair, inhibition, Kmisrepair
!endif
end subroutine

!--------------------------------------------------------------------------
! Hill function for max checkpoint delay TCP
!--------------------------------------------------------------------------
function f_TCP(ccp,n) result(TCP)
integer :: n
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: TCP

TCP = ccp%bTCP*n/(ccp%aTCP + n)
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

#if 0

!--------------------------------------------------------------------------
subroutine old_timestep(cp, ccp, dt)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: dt
integer :: phase, ityp, nPL, kpar=0
real(REAL_KIND) :: cf_O2, cf_glucose, pcp_O2, pcp_glucose, pcp_starvation, R
logical :: switch

phase = cp%phase
if (cp%dVdt == 0) then
!	if (phase == G1_phase .or. phase == S_phase .or. phase == G2_phase) then
		write(nflog,*) 'dVdt=0, kcell, phase: ',kcell_now,phase
		stop
!	endif
endif
!nPL = cp%NL1
nPL = cp%N_PL
ityp = cp%celltype
!if (.not.use_metabolism) then
!	if (.not.colony_simulation .and. (phase == Checkpoint1)) then    ! check for starvation arrest
!		cf_O2 = cp%Cin(OXYGEN)/anoxia_threshold
!		pcp_O2 = getPcp_release(cf_O2,dt)
!		cf_glucose = cp%Cin(GLUCOSE)/aglucosia_threshold
!		pcp_glucose = getPcp_release(cf_glucose,dt)
!		pcp_starvation = pcp_O2*pcp_glucose
!		if (pcp_starvation == 0) then
!			return
!		elseif (pcp_starvation < 1) then
!			R = par_uni(kpar)
!			if (R < pcp_starvation) return
!		endif
!	endif
!endif
if (phase == G1_phase) then
    if (use_volume_based_transition) then
        switch = (cp%V > cp%G1_V)
    else
        switch = (tnow > cp%G1_time)
    endif
    if (switch) then
        cp%phase = G1_checkpoint
        cp%G1_flag = .false.
        cp%G1S_time = tnow + f_TCP(ccp,nPL)		!ccp%Tcp(nPL)
    endif
elseif (phase == G1_checkpoint) then  ! this checkpoint combines the release from G1 delay and the G1S repair check
    if (.not.cp%G1_flag) then
        R = par_uni(kpar)
        cp%G1_flag = (R < ccp%Pk_G1*dt)
    endif
    cp%G1S_flag = (nPL == 0 .or. tnow > cp%G1S_time)
    if (use_metabolism) then
		cp%G1S_flag = cp%G1S_flag .and. (cp%metab%A_rate > ATPg)
	endif
    if (cp%G1_flag .and. cp%G1S_flag) then
        cp%phase = S_phase
! Note: now %I_rate has been converted into equivalent %dVdt, to simplify code 
!	    cp%S_time = tnow + (max_growthrate(ityp)/cp%dVdt)*cp%fg*ccp%T_S
!	    cp%S_time = tnow + (max_growthrate(ityp)/cp%dVdt)*ccp%T_S
	    cp%S_duration = (max_growthrate(ityp)/cp%dVdt)*ccp%T_S
	    cp%S_time = 0   ! this is now the amount of time progress through S phase: 0 -> %S_duration
!	    endif
    endif
elseif (phase == S_phase) then
    cp%arrested = (cp%dVdt/max_growthrate(ityp) < ccp%arrest_threshold)
!    if (kcell_now == 1) then
!        write(*,'(a,2e12.3,2x,L)') 'S_phase: dVdt, fraction, arrested: ',cp%dVdt,cp%dVdt/max_growthrate(ityp),cp%arrested
!    endif
    if (.not.cp%arrested) then
        cp%S_time = cp%S_time + dt
    endif
    switch = (cp%S_time >= cp%S_duration)
!   switch = (tnow > cp%S_time)
    if (switch) then
        cp%phase = G2_phase
! Note: now %I_rate has been converted into equivalent %dVdt, to simplify code
		cp%G2_time = tnow + (max_growthrate(ityp)/cp%dVdt)*ccp%T_G2
    endif
elseif (phase == G2_phase) then
    if (use_volume_based_transition) then
        switch = (cp%V > cp%G2_V)
    else
! Note: now %I_rate has been converted into equivalent %dVdt, to simplify code
!		if (use_metabolism) then
!			switch = (tnow > cp%G2_time .and. cp%metab%Itotal > cp%metab%I2divide) ! try this to prevent volumes decreasing 
!		else
			switch = (tnow > cp%G2_time .and. cp%V > cp%divide_volume) ! try this to prevent volumes decreasing 
!		endif
    endif
    if (switch) then
        cp%phase = G2_checkpoint
        cp%G2_flag = .false.
        cp%G2M_time = tnow + f_TCP(ccp,nPL)		!ccp%Tcp(nPL)
    endif
elseif (phase == G2_checkpoint) then ! this checkpoint combines the release from G2 delay and the G2M repair check
    if (.not.cp%G2_flag) then
        R = par_uni(kpar)
        cp%G2_flag = (R < ccp%Pk_G2*dt)
    endif
    cp%G2M_flag = (nPL == 0 .or. tnow > cp%G2M_time)
    if (use_metabolism) then
		cp%G2M_flag = cp%G2M_flag .and. (cp%metab%A_rate > ATPg)
	endif
    if (cp%G2_flag .and. cp%G2M_flag) then
        cp%phase = M_phase
        cp%M_time = tnow + ccp%T_M   
    endif
elseif (phase == M_phase) then
    if (tnow > cp%M_time) then
        cp%phase = dividing
!        cp%doubling_time = tnow
    endif
endif    
if (nPL > 0 .and. .not.cp%irrepairable) then
    call radiation_repair(cp, ccp, dt)
endif
end subroutine

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
!--------------------------------------------------------------------------
subroutine old_radiation_damage(cp, ccp, dose, SER_OER, tmin)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: dose, SER_OER(2), tmin
real(REAL_KIND) :: dDdt, rnL(2), rnPL, drnPLdt, drnLdt(2), dtmin, dthour, fraction, Krepair, Kmissum, dr, R
integer :: nt, it, n, dn, k, ityp, kpar=0

dtmin = 0.01
dthour = dtmin/60
nt = max(tmin/dtmin, 1.0)
dDdt = dose/(nt*dthour)     ! Gy/h
rnPL = cp%NL1
rnL = cp%NL2
Kmissum = ccp%Kmisrepair
ityp = cp%celltype
if (cp%phase < S_phase) then
    Krepair = ccp%Krepair_base
elseif (cp%phase > S_phase) then
    Krepair = ccp%Krepair_max
else
    if (use_volume_based_transition) then   ! fraction = fraction of passage through S phase
        fraction = (cp%V - cp%G1_V)/(cp%S_V - cp%G1_V)
    else
!        fraction = 1 - (cp%S_time - tnow)/ccp%T_S
        fraction = 1 + min(cp%S_time/cp%S_duration, 1.0)
    endif
    Krepair = ccp%Krepair_base + fraction*(ccp%Krepair_max - ccp%Krepair_base)
endif
do it = 1,nt
    drnPLdt = ccp%eta_PL*SER_OER(1)*dDdt - Krepair*rnPL - Kmissum*rnPL**2
    rnPL = rnPL + drnPLdt*dthour
    drnLdt = ccp%eta_L*SER_OER(2)*dDdt + ccp%Kmisrepair*rnPL**2
    rnL = rnL + drnLdt*dthour
enddo
n = rnPL
dr = rnPL - n
R = par_uni(kpar)
dn = 0
if (R < dr) dn = 1
cp%NL1 = cp%NL1 + n + dn
do k = 1,2
    n = rnL(k)
    dr = rnL(k) - n
    R = par_uni(kpar)
    dn = 0
    if (R < dr) dn = 1
    cp%NL2(k) = cp%NL2(k) + n + dn
enddo
end subroutine

!--------------------------------------------------------------------------
! Time unit = hour
! This may need to be changed, because it implicitly assumes that no more
! than one repair and one misrepair of each type can occur within a time step.
! The fix would be to subdivide the time step, as in the damage subroutine.
!--------------------------------------------------------------------------
subroutine old_radiation_repair(cp, ccp, dt)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: dt
integer :: i, ityp, kpar=0
real(REAL_KIND) :: dthour, fraction, Krepair, R

dthour = dt/3600    ! sec -> hour
ityp = cp%celltype
if (cp%phase < S_phase) then
    Krepair = ccp%Krepair_base
elseif (cp%phase > S_phase) then
    Krepair = ccp%Krepair_max
else
    if (use_volume_based_transition) then   ! fraction = fraction of passage through S phase
        fraction = (cp%V - cp%G1_V)/(cp%S_V - cp%G1_V)
    else
!        fraction = 1 - (cp%S_time - tnow)/ccp%T_S
        fraction = 1 + min(cp%S_time/cp%S_duration, 1.0)
    endif
    Krepair = ccp%Krepair_base + fraction*(ccp%Krepair_max - ccp%Krepair_base)
endif
R = par_uni(kpar)
if (R < cp%NL1*Krepair*dthour) then
    cp%NL1 = cp%NL1 - 1
    if (cp%NL1 == 0) return
endif
do i = 1,2
    R = par_uni(kpar)
    if (R < cp%NL1**2*ccp%Kmisrepair(i)*dthour) then
        cp%NL1 = cp%NL1 - 1
        cp%NL2(i) = cp%NL2(i) + 1
        if (cp%NL1 == 0) return
    endif
enddo
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
program main
use cycle_mod
use Tcp_mod
implicit none
integer :: istep, kcell, nd
real(REAL_KIND) :: tnow, Td
real(REAL_KIND) :: DELTA_T, hours, dose, tmin
integer :: Ncells, Nsteps, i
type(cell_type), allocatable, target :: cell_list(:)
type(cell_type), pointer :: cp
type(cycle_parameters_type), target :: cc_parameters
type(cycle_parameters_type), pointer :: ccp
logical :: all_divided

ccp => cc_parameters(1)
ccp%eta_L = 0.14                    ! /Gy
ccp%eta_PL = 0.6                    ! /Gy
ccp%T_G1 = 12                       ! h
ccp%T_S = 9                         ! h
ccp%T_G2 = 1                        ! h
ccp%T_M = 0.5                       ! h
ccp%G1_mean_delay = 3               ! h
ccp%G2_mean_delay = 2               ! h
ccp%Pk_G1 = 1./ccp%G1_mean_delay    ! /h
ccp%Pk_G2 = 1./ccp%G2_mean_delay    ! /h
ccp%Krepair = 0.1
ccp%Kmisrepair(1:3) = 0.01852 = (1./3)*(0.5/9)
ccp%Kcp = 0.13

DELTA_T = 0.05   ! hours
Ncells = 100
allocate(cell_list(Ncells))

! test radiation dose
cp => cell_list(1)
dose = 5
tmin = 5
do i = 1,10
    cp%NL1 = 0
    cp%NL2 = 0
    call radiation_damage(cp, ccp, dose, tmin)
    write(*,*) 'radiation damage: ',cp%NL1,cp%NL2
enddo
stop

!write(*,*) 'makeTCP'
!call makeTCP(ccp%tcp,NTCP,ccp%epsilon_PL,ccp%epsilon_2PL,ccp%Kcp) ! set checkpoint repair time limits 
!write(*,*) 'did makeTCP'

hours = 10*24
Nsteps = hours/DELTA_T
tnow = 0
do icell = 1,Ncells
    cp => cell_list(icell)
    cp%phase = G1_phase
    cp%G1_flag = .false.
    cp%G1_time = tnow + ccp%T_G1
    cp%NL1 = 100
    cp%NL2 = 0
enddo
do istep = 1,Nsteps
    all_divided = .true.
    do kcell = 1,Ncells
        cp => cell_list(kcell)
        if (cp%phase == divided) cycle
        all_divided = .false.
        tnow = istep*DELTA_T
        if (use_exponential_cycletime) then
            call exp_timestep(cp, ccp, tnow, DELTA_T)
        else
            call log_timestep(cp, ccp, tnow, DELTA_T)
        endif
    enddo
    if (all_divided) then
        write(*,*) 'All cells have divided'
        exit
    endif
enddo
do kcell = 1,Ncells
    cp => cell_list(kcell)
    write(*,'(a,i6,f6.2,5i6)') 'Cell: ',kcell,cp%doubling_time, cp%phase, cp%NL1, cp%NL2
enddo
stop

Td = 0
nd = 0
do kcell = 1,Ncells
    if (cell_list(kcell)%phase == divided) then
        nd = nd+1
        Td = Td + cell_list(kcell)%doubling_time
    endif
enddo
write(*,'(a,i6,f6.2)') 'Average doubling time: ',nd,Td/nd
end program

#endif