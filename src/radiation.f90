! Radiation damage/repair
module radiation
use global
!use cycle_mod
implicit none

! Radiation damage/repair
!real(REAL_KIND) :: eta_PL, eta_L(2), Kcp
!real(REAL_KIND) :: Krepair_base, Krepair_max, Kmisrepair(2)
!real(REAL_KIND) :: Tcp(0:NTCP)

! Revised:
!drop eta_L(2)
!add eta_IRL	(irrepairable lesion)
!Kmisrepair(2) -> Kmisrepair, fraction_Ch1	! Kmisrepair = total misrepair rate (NHEJ + MMEJ), fraction_Ch1 = fraction of Ch1 generated
!mitosis_factor = increase of Kmisrepair during mitosis

! Cell properties
!integer :: NL1, NL2(2)

! Now the cell needs:
! integers:
!N_PL, N_Ch1, N_Ch2
! logical:
!irrepairable

! Notes: 
! High fidelity repair can occur only in S, G2

contains

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
subroutine new_radiation_damage(cp, ccp, dose, SER_OER, tmin)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: dose, SER_OER, tmin
real(REAL_KIND) :: dDdt, dtmin, dthour, R
real(REAL_KIND) :: p_PL, p_IRL, Krepair, Kmisrepair, misrepair_factor, fraction
integer :: nt, it, nPL, nIRL, ityp, kpar=0
logical :: do_repair = .true.

dtmin = 0.01
dthour = dtmin/60
nt = max(tmin/dtmin, 1.0)
dDdt = dose/(nt*dthour)     ! Gy/h
nPL = cp%N_PL
if (cp%irrepairable) then
	write(*,*) 'radiation_damage: irrepairable cell!'
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
			fraction = 1 - (cp%S_time - tnow)/ccp%T_S(ityp)
		endif
		Krepair = ccp%Krepair_base + fraction*(ccp%Krepair_max - ccp%Krepair_base)
	endif
	if (cp%phase == M_phase) then
		misrepair_factor = ccp%mitosis_factor
	else
		misrepair_factor = 1
	endif
	Kmisrepair = misrepair_factor*ccp%Kmisrepair(1)	! -> scalar
endif

p_PL = ccp%eta_PL*SER_OER*dDdt*dthour
p_IRL = ccp%eta_IRL*SER_OER*dDdt*dthour
do it = 1,nt
	R = par_uni(kpar)
	if (R < p_PL) nPL = nPL + 1
	R = par_uni(kpar)
	if (R < p_IRL) nIRL = nIRL + 1
	if (do_repair) then
		! repair/misrepair
		R = par_uni(kpar)
		if (R < nPL*Krepair*dthour) then
			nPL = nPL - 1
		endif
		R = par_uni(kpar)
		if (R < nPL**2*Kmisrepair*dthour) then	! -> scalar
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
cp%irrepairable = (nIRL > 0)
end subroutine

!--------------------------------------------------------------------------
! Time unit = hour
! This may need to be changed, because it implicitly assumes that no more
! than one repair and one misrepair of each type can occur within a time step.
! The fix would be to subdivide the time step, as in the damage subroutine.
!--------------------------------------------------------------------------
subroutine new_radiation_repair(cp, ccp, dt)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: dt
integer :: i, ityp, nPL, kpar=0
real(REAL_KIND) :: dthour, fraction, Krepair, Kmisrepair, misrepair_factor, R
integer :: nt, it

nt = 10
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
        fraction = 1 - (cp%S_time - tnow)/ccp%T_S(ityp)
    endif
    Krepair = ccp%Krepair_base + fraction*(ccp%Krepair_max - ccp%Krepair_base)
endif
if (cp%phase == M_phase) then
	misrepair_factor = ccp%mitosis_factor
else
	misrepair_factor = 1
endif
Kmisrepair = misrepair_factor*ccp%Kmisrepair(1)	! -> scalar
nPL = cp%N_PL
do it = 1,nt
	R = par_uni(kpar)
	if (R < nPL*Krepair*dthour) then
		nPL = nPL - 1
		if (nPL == 0) exit
	endif
	R = par_uni(kpar)
	if (R < nPL**2*Kmisrepair*dthour) then
		nPL = nPL - 1
		R = par_uni(kpar)
		if (R < ccp%fraction_Ch1) then
			cp%N_Ch1 = cp%N_Ch1 + 1
		else
			cp%N_Ch2 = cp%N_Ch2 + 1
		endif
		if (nPL == 0) exit
	endif
enddo
cp%N_PL = nPL
end subroutine

end module