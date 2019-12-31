! To test simple cell metabolism model
! Concentration of ATP varies by < 10%  https://en.wikipedia.org/wiki/Glycolysis#Intermediates_for_other_pathways

! Units:
!     time				s = seconds
!     distance			cm
!     volume			cm^3
!     mass				micromole = 10^-6 mol = mumol
!     flux				mumol/s
!     concentration		mumol/cm^3 = mM
!
! Need to comment out 'use chemokine' when used in the test program metab.exe
!
! Question: How do the results of this model translate into cell rate of volume growth? 
!
! This version includes lumped amino acids - AA
!----------------------------------------------------------------------------------------------------------------
module FCN_mod
use real_kind_mod
implicit none

real(REAL_KIND) :: FCN_a1, FCN_b1, FCN_c1
real(REAL_KIND) :: FCN_r_Pm_base, FCN_fPn, FCN_Km_P, FCN_fGn, FCN_r_G, FCN_N_GA, FCN_V, FCN_K1, FCN_K2, FCN_C_L

end module

!----------------------------------------------------------------------------------------------------------------
module metabolism
use real_kind_mod
use global
#if .not. EXCEL
use chemokine
#endif
implicit none

integer, parameter :: FGP_SOLVER_MAXATP_TANDEM = 1
integer, parameter :: FGP_SOLVER_MAXATP_STAGED = 2
integer, parameter :: FGP_SOLVER_SURVIVAL_STAGED = 3

logical :: use_glutamine, from_excel, solved 

real(REAL_KIND), parameter :: f_I_threshold = 0.5	! NOT USED
real(REAL_KIND) :: I_threshold
! From spheroid-abm, the max rates of consumption of oxygen and glucose are: 
!   oxygen:  6.25e-17 moles/cell/s
!   glucose: 6.80e-17 moles/cell/s
! We work with mumol/cell/sec, and convert these values by scaling by 1.0e6, to give
!   oxygen:  6.25e-11 mumol/cell/s
!   glucose: 6.80e-11 mumol/cell/s

real(REAL_KIND) :: Hill_Km_O2
real(REAL_KIND) :: Hill_N_O2
real(REAL_KIND) :: Hill_Km_G	! Hill Km for dependence of glycolysis rate on glucose
real(REAL_KIND) :: Hill_N_G		! Hill N for dependence of glycolysis rate on glucose 
real(REAL_KIND) :: Hill_Km_Gln	! Hill Km for dependence of glutamine metabolism rate on glutamine
real(REAL_KIND) :: Hill_N_Gln	! Hill N for dependence of glutamine metabolism rate on glutamine 
real(REAL_KIND) :: Hill_Km_P    ! Hill Km for dependence of pyruvate oxidation rate on pyruvate
real(REAL_KIND) :: Hill_N_P		! Hill N for dependence of pyruvate oxidation rate on pyruvate
real(REAL_KIND) :: Hill_Km_C	! Hill N for dependence of glycolysis on normalised oxygen consumption rate (Pasteur effect)
real(REAL_KIND) :: K_H1			! HIF-1 k1
real(REAL_KIND) :: K_H2			! HIF-1 k2
real(REAL_KIND) :: K_Hb			! HIF-1 kb
real(REAL_KIND) :: K_PDK		! K_PDK
real(REAL_KIND) :: PDKmin		! PDKmin
real(REAL_KIND) :: C_O2_norm    ! Note: these _norm values correspond to unconstrained growth
real(REAL_KIND) :: C_G_norm
real(REAL_KIND) :: C_Gln_norm
real(REAL_KIND) :: C_L_norm
real(REAL_KIND) :: C_A_norm
real(REAL_KIND) :: O2_baserate
real(REAL_KIND) :: G_baserate
real(REAL_KIND) :: Gln_baserate
real(REAL_KIND) :: K_PL			! P -> L
real(REAL_KIND) :: K_LP			! L -> P
real(REAL_KIND) :: r_Pu, r_Gu, r_Glnu, r_Au, r_Iu, r_Ou, r_Lu, C_Pu   ! unconstrained conditions
real(REAL_KIND) :: G_maxrate, O2_maxrate, Gln_maxrate
real(REAL_KIND) :: f_N
integer :: fgp_solver

type param_set_type
real(REAL_KIND) :: a0, b0, c0, d0, a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3, p, q, h
end type
!type param_set_type
!real(REAL_KIND) :: a0, b0, c0, a1, b1, a2, b2, a3, b3, c3, d3, a4, b4, c4, h
!end type
logical :: use_nitrogen = .true.

type(param_set_type) :: ps

real(REAL_KIND), parameter :: average_volume = 1.2
real(REAL_KIND), parameter :: r_H = 1.21e-4
real(REAL_KIND), parameter :: d_H = 2.3e-3
real(REAL_KIND), parameter :: Km_H = 20*0.18/160	! Kelly had 2 mmHg, Bill suggests 20

real(REAL_KIND), parameter :: C_Gln_Km_factor = 0.2
#if EXCEL
real(REAL_KIND), parameter :: PI = 4*atan(1.0)
#endif
logical :: use_Kelly = .false.

contains

!--------------------------------------------------------------------------
! This is NOT the normalised rate
! H = HIF1
! C_G = glucose concentration
! r_O = current O2 consumption rate, i.e. from previous time step
! Note: need to set to parent cell value (or r_O_norm) for a newly created cell
! The idea is that when the rate of pyruvate and glutamine oxidation goes low,
! (which is indicated by r_O going low), the citrate and ATP production rates go low,
! and this enhances glycolysis.
! Pasteur Effect.
!--------------------------------------------------------------------------
function get_glycosis_rate(H, C_G, C_Gln, r_O) result(rate)
integer :: ityp
real(REAL_KIND) :: H, C_G, C_Gln, r_O, rate
real(REAL_KIND) :: metab, r_O_fract, cfactor, MM_Gln
integer :: N_Gln
!real(REAL_KIND) :: Km = 0.02	! guess

metab = glucose_metab(C_G)
!O2factor = LowO2Factor(C_O2)
r_O_fract = 1   !min(r_O/r_Ou,1.0)  !TODO !!!!!!!!!!!!!!!!!!!!!
cfactor = 1     !r_O_fract/(Hill_Km_C + r_O_fract)
!write(nflog,'(a,2e12.3,2f6.3)') 'r_O: ',r_O,r_O_norm,r_O_fract,O2factor
!metab = cfactor*metab
!N_Gln = int(Hill_N_Gln)
!MM_Gln = f_MM(C_Gln,Hill_Km_Gln,N_Gln)      ! need a switch for applying this factor
rate = G_maxrate*(1 + K_Hb*H)*metab !*MM_Gln
end function

!--------------------------------------------------------------------------
! MM_O2 is the Michaelis-Menten factor for O2
!--------------------------------------------------------------------------
function get_glutamine_rate(C_Gln, fPDK, MM_O2, C_G) result(r_Gln)
real(REAL_KIND) :: C_Gln, fPDK, MM_O2, C_G
real(REAL_KIND) :: r_Gln, MM_G
integer :: N_Gln, N_G

N_Gln = int(Hill_N_Gln)
N_G = int(Hill_N_G)
MM_G = f_MM(C_G,Hill_Km_G,N_G)
r_Gln = fPDK*MM_O2*Gln_maxrate*f_MM(C_Gln,Hill_Km_Gln,N_Gln)    !*MM_G
end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function LowO2Factor(C_O2) result(factor)
real(REAL_KIND) :: C_O2, factor
real(REAL_KIND) :: Km = 0.001 ! mM
factor = C_O2/(Km + C_O2)
end function

!--------------------------------------------------------------------------
! Currently this sets up parameters for type 1 cells only.
! For tandem case we need to incorporate glutamine.
!--------------------------------------------------------------------------
subroutine SetupMetabolism(mp,ok)
type(metabolism_type), pointer :: mp
logical :: ok
integer :: ityp, it, N_P, N_O2
real(REAL_KIND) :: MM_O2, MM_P, V, K1, K2, Km_P, C_P, C_O2, C_G, C_L, C_Gln !, C_Gu, r_Gu, r_Pu, r_Lu, r_An, r_In
real(REAL_KIND) :: Km_O2_factor = 1
real(REAL_KIND) :: average_volume = 1.2
real(REAL_KIND) :: H, fPDK, rP1, rP2, rP3

if (noSS) then 
    write(nflog,*) 'SetupMetabolism: not using SS solver'
else
    write(nflog,*) 'SetupMetabolism: using SS solver'
endif
ok = .true.
V = Vcell_cm3*average_volume
!mp => phase_metabolic(1)

Hill_Km_O2 = chemo(OXYGEN)%MM_C0
Hill_N_O2 = chemo(OXYGEN)%Hill_N
Hill_Km_G = chemo(GLUCOSE)%MM_C0
Hill_N_G = chemo(GLUCOSE)%Hill_N
Hill_Km_Gln = chemo(GLUTAMINE)%MM_C0
Hill_N_Gln = chemo(GLUTAMINE)%Hill_N
Hill_N_P = 1
Hill_Km_O2 = Km_O2_factor*Hill_Km_O2
N_O2 = chemo(OXYGEN)%Hill_N
O2_maxrate = chemo(OXYGEN)%max_cell_rate
write(nflog,*) 'O2_maxrate: ',O2_maxrate
G_maxrate = chemo(GLUCOSE)%max_cell_rate
Gln_maxrate = chemo(GLUTAMINE)%max_cell_rate
O2_baserate = 0		! are these needed at all??
G_baserate = 0
H = get_HIF1steadystate(C_O2_norm)
call analyticSetPDK1(H, fPDK, 1.0d10)
mp%HIF1 = H
mp%PDK1 = fPDK

if (use_nitrogen) then
    r_Glnu = f_N*Gln_maxrate
    call get_unconstrained_rates
!    write(nflog,'(a,6e12.3)') 'got unconstrained rates: G,P,I,Gln,A,O: ',r_Gu,r_Pu,r_Iu,r_Glnu,r_Au,r_Ou
    call set_param_set(r_Gu,C_L_norm)
    mp%f_G = f_Gu
    mp%f_P = f_Pu
    mp%f_Gln = f_Glnu
    
else
K1 = K_PL
K2 = K_LP
N_P = 1
Km_P = Hill_Km_P
C_O2 = C_O2_norm
C_G = C_G_norm
C_L = C_L_norm
C_Gln = C_Gln_norm
MM_O2 = f_MM(C_O2,Hill_Km_O2,N_O2)

r_Ou = O2_maxrate	! initial guess, for f_metab

mp%f_G = f_Gu
mp%f_P = f_Pu
mp%f_Gln = f_Glnu
mp%O_rate = r_Ou
call f_metab(mp, C_O2, C_G, C_L, C_Gln)
!write(nflog,*) 'f_Gln: ',mp%f_Gln
!write(nflog,'(a,6e11.3)') 'G,P,L,A,I,Gn rates: ',mp%G_rate,mp%P_rate,mp%L_rate,mp%A_rate,mp%I_rate,mp%Gln_rate

r_Gu = mp%G_rate
r_Pu = mp%P_rate
r_Au = mp%A_rate
r_Iu = mp%I_rate
r_Glnu = mp%Gln_rate
r_Ou = mp%O_rate
C_Pu = mp%C_P
!C_A_norm = 20		!#################  From Bill #################  now an input parameter
!mp%f_G = f_Gln
!mp%f_P = f_Pn
!mp%C_P = C_Pu
!mp%C_A = C_A_norm
!mp%A_rate = r_Au
!mp%I_rate = r_Iu
!mp%G_rate = r_Gu
!mp%I_rate_max = r_Iu
write(nflog,'(a,4e12.3)') 'f_Gu,f_Pu,r_Ag,C_Pu: ',f_Gu,f_Pu,r_Ag,C_Pu
endif
r_Ag = f_ATPg*r_Au
r_As = f_ATPs*r_Au
mp%recalcable = -1
end subroutine

!--------------------------------------------------------------------------
! Use K_H1 for the exponent, K_H2 for the rate coefficient
!--------------------------------------------------------------------------
function get_HIF1steadystate(C_O) result(H)
real(REAL_KIND) :: C_O, H
real(REAL_KIND) :: x
real(REAL_KIND) :: C_O_max = 0.18

x = min(C_O/C_O_max,1.0)
H = (1-x)**K_H1
!write(*,*) 'get_HIF1steadystate: K_H1, K_H2, K_HB: ', K_H1, K_H2, K_Hb, K_PDK, PDKmin
!write(*,*) 'C_O, x, H: ',C_O, x, H
end function

!--------------------------------------------------------------------------
!With the Kelly2008 model:
!  a = rH
!  b = 1 + dH*C_O/(rH*(Km + C_O))
!--------------------------------------------------------------------------
subroutine analyticSetHIF1z(C_O, H, dt)
real(REAL_KIND) :: C_O, H, dt
real(REAL_KIND) :: a, b, c, ee, H0

if (use_Kelly) then
	a = r_H
	b = 1 + d_H*C_O/(r_H*(Km_H + C_O))
else
	ee = K_H1*C_O
	if (ee > 100) then
		H = 0
		return
	endif
	a = K_H2
	b = exp(ee)
endif
H0 = H
c = 1 - b*H0
H = (1 - c*exp(-a*b*dt))/b
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine analyticSetHIF1(C_O, H, dt)
integer :: ityp
real(REAL_KIND) :: C_O, H, dt
real(REAL_KIND) :: x, H0, Heq
real(REAL_KIND) :: C_O_max = 0.18

x = min(C_O/C_O_max,1.0)
Heq = (1-x)**K_H1
H0 = H
H = Heq + (H0 - Heq)*exp(-K_H2*dt)
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine analyticSetPDK1(H, P, dt)
integer :: ityp
real(REAL_KIND) :: P, H, dt
real(REAL_KIND) :: a, b, c, d, P0

a = K_PDK
d = 1 - PDKmin
b = (1 - d*H)
P0 = P
c = P0 - b
P = b + c*exp(-a*dt)
end subroutine

!--------------------------------------------------------------------------
! The fraction of pyruvate that is converted to acetyl-CoA depends on the
! rate of glycolysis.  The normalised rate has a maximum value of 1 under
! normoxic conditions. This corresponds to the minimum pyruvate oxidation fraction
! (equivalently, the maximum lactate fraction).  The minimum pyruvate oxidation
! fraction is the fraction of pyruvate that is directed to the Krebs cycle
! when both glucose and oxygen are in ample supply.  
! In fact this is the nominal minimum, and the fraction can be even less 
! if hypoxia elevates glycosis rate.
! The basic idea is that over a range of glycolysis rate, the total rate of
! production of ATP is a constant.
! The total is the sum of ATP produced by glycolysis and by pyruvate oxidation.
! If the glycolysis rate is reduced, the intermediate production rate can be maintained 
! by increasing the pyruvate oxidation fraction (reducing lactate production), to the limit
! of fraction = 1.  Further reductions in glycolysis will then reduce ATP production.
! REVISED
! Glycolysis:
! A fraction N_GI goes to make intermediates, I_rate = N_GI*G_rate
! the remainder makes pyruvate, PP_rate = N_GA*(1 - N_GI)*G_rate
! and A_rate = PP_rate
! Pyruvate oxidation:
! A fraction N_PI goes to make intermediates via the TCA cycle, I_rate = N_PI*P_rate
! the remainder (1 - N_PI) is fully oxidised, producing N_PA ATP/mole, (N_PA = 18)
! and A_rate from pyruvate oxidation = N_PA*(1 - N_PI)*P_rate
! Bill:
! The intermediates used for anabolism are downstream of acetyl-CoA (including acetyl-CoA itself).
! Yes, PDK1 reduces pyruvate utilisation (r) but this is upstream of acetyl-CoA. So the formalism 
! needs to have the intermediates coming from acetyl-CoA rather than pyruvate itself.
!--------------------------------------------------------------------------
subroutine get_metab_rates(mp, Cin)
integer :: ityp
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: Cin(:)

!if (noSS) then
!    call f_metab_noSS(mp,Cin(OXYGEN),Cin(GLUCOSE),Cin(LACTATE),Cin(4))
!else
    call f_metab(mp,Cin(OXYGEN),Cin(GLUCOSE),Cin(LACTATE), Cin(GLUTAMINE))
!endif
end subroutine

!--------------------------------------------------------------------------
! Use:
!	r_G for dG/dt, the rate of glycolysis = rate of glucose consumption
!	f_G for N_GI, the fraction of dG/dt that goes to intermediates
!	r_P for dP/dt, rate of utilisation of pyruvate
!	f_P for N_PI, the fraction of dP/dt that goes to intermediates
!	r_A for dA/dt, rate of ATP production
!	r_I for dI/dt, rate of intermediates production
!
!	r_Gu for dG/dt under normal conditions, no nutrient constraints, H = 1
!	f_Gu for the f_G under normal conditions (upper bound of f_G)
!	r_Pu for dP/dt under normal conditions
!	f_Pu for f_P under normal conditions (upper bound of f_P)
!	r_Au for r_A under normal conditions
!	r_Iu for r_I under normal conditions
!	alpha = r_P as a fraction of dP/dt under normal conditions
!
!	C_G for IC glucose concentration
!	C_O2 for IC oxygen concentration
!	C_P for IC pyruvate concentration
!	C_L for IC lactate concentration
!	
! r_P = N_GA*(1 - f_G)*r_G + V*(K2*C_L - K1*C_P - dC_P/dt) = fPDK*r_P_max*MM(O2)*MM(C_P)
! with the constraint that C_P >= 0
! Steady-state approach may not be feasible, because if there is 
! plenty of glucose but O2 is very low, r_G will be high but r_P will tend
! towards 0.  This must lead to an increase in C_P.
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Use f_Gln for glutamine factor, to avoid confusion with f_Gln
! (Note that this assumes that C_Gn is given by C_Gn_norm - no longer)
!--------------------------------------------------------------------------
subroutine f_metab(mp, C_O2_, C_G_, C_L_, C_Gln_)
integer :: ityp
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C_O2_, C_G_, C_L_, C_Gln_
real(REAL_KIND) :: C_O2, C_G, C_L, C_Gln
real(REAL_KIND) :: r_G, fPDK
real(REAL_KIND) :: f_G, f_P, f_Gln, r_P, r_A, r_I, r_L, r_Gln
real(REAL_KIND) :: K1, K2, C_P
real(REAL_KIND) :: r_GP, r_GA, r_PA, r_Pm, V, Km_O2, Km_P, Km_G, Km_Gln, a, b, c, d, e, MM_P, MM_O2, MM_G, MM_Gln, Km_GO
real(REAL_KIND) :: r_GI, r_PI, r_O2
real(REAL_KIND) :: F, p, q, scale, dr_I, dr_P, dC_P, w, dw, r_I_prev
real(REAL_KIND) :: FF, z
!real(REAL_KIND) :: w0, w1, F0, F1, dFdw, rtol, w0prev, z, H, q1, q2, q3, q4, q5
integer :: N_O2, N_P, it, k
logical :: dbug = .false.
logical, save :: first = .true.
logical :: bills_idea = .true.

z = 1   ! Try this first
!write(nflog,*)
C_O2 = max(0.0,C_O2_)
C_G = max(0.0,C_G_)
C_L = max(0.0,C_L_)
C_Gln = max(0.0,C_Gln_)

if (use_nitrogen) then
!    if (mp%recalcable > 0) then
!        call f_metab_recalc(mp, C_O2, C_G, C_L, C_Gln)
!    else
        call f_metab_nitrogen(mp, C_O2, C_G, C_L, C_Gln)
!    endif
    return
endif

N_O2 = Hill_N_O2
Km_O2 = Hill_Km_O2
N_P = 1
Km_P = Hill_Km_P
Km_Gln = Hill_Km_Gln
if (dbug) write(nflog,*) 'Hill_Km_Gln,Hill_N_Gln: ',Hill_Km_Gln,Hill_N_Gln
V = Vcell_cm3*average_volume		! should be actual cell volume cp%V
K1 = K_PL
K2 = K_LP
if (use_glutamine) then
    MM_Gln = C_Gln/(Km_Gln + C_Gln)
else
    MM_Gln = 1
endif
f_G = mp%f_G!*MM_Gln
f_P = mp%f_P!*MM_Gln
f_Gln = mp%f_Gln!*MM_Gln
C_P = mp%C_P
!write(nflog,'(a,L1,6e12.3)') 'metab: C_Gln,Km_Gln,MM_Gln,f_G,f_P,f_Gln: ',use_glutamine,C_Gln,Km_Gln,MM_Gln,f_G,f_P,f_Gln
if (dbug) write(nflog,'(a,5e12.3)') 'metab: C_O2,C_G,C_L,C_Gln,C_P: ',C_O2,C_G,C_L,C_Gln,C_P
r_O2 = mp%O_rate
mp%G_rate = get_glycosis_rate(mp%HIF1,C_G,C_Gln,r_O2)	! Note: this is the previous O_rate
r_G = mp%G_rate
!write(nflog,*) 'f_metab r_G: ',r_G
fPDK = mp%PDK1
MM_O2 = f_MM(C_O2,Km_O2,N_O2)
!r_Pm_base = fPDK*MM_O2*O2_maxrate/N_PO	! note that MM_P is not here, since it varies it is added as needed
!	write(nflog,'(a,2e12.3)') 'r_Pm_base: ',r_Pm_base
! r_O2 = fact_O2*C_P/(Km_P + C+P)

if (use_glutamine) then
    MM_Gln = C_Gln/(Km_Gln + C_Gln)
	r_Gln = get_glutamine_rate(C_Gln, fPDK, MM_O2, C_G)
    if (bills_idea) then
        FF = O2_maxrate*MM_O2/(1+z)
    else
        F = fPDK*O2_maxrate*MM_O2
    endif
else
    r_Gln = 0
    MM_Gln = 0
    if (bills_idea) then
        FF = O2_maxrate*MM_O2
    else
        F = fPDK*O2_maxrate*MM_O2
    endif
endif

!write(nflog,'(a,4e12.3)') 'fPDK, MM_O2, MM_Gln, F: ',fPDK,MM_O2,MM_Gln,F

p = (r_G*(1 - f_G)*N_GA + V*K2*C_L)*(1 - f_P)*N_PO + r_Gln*(1 - f_Gln)*N_GlnO
if (bills_idea) p = p - FF*z*MM_Gln
q = V*K1*(1 - f_P)*N_PO
if (dbug) write(nflog,'(a,5e12.3)') 'r_O2,r_G,r_Gln, p,q: ',r_O2,r_G,r_Gln, p,q
a = q
if (bills_idea) then
    b = FF*fPDK + q*Km_P - p
else
    b = F + q*Km_P - p
endif
c = -p*Km_P
d = sqrt(b*b - 4*a*c) 
if (dbug) write(nflog,'(a,4e12.3)') 'a,b,c,d: ',a,b,c,d
C_P = (-b + d)/(2*a)
!write(nflog,*) 'f_metab: C_P: ',C_P

r_P = r_G*(1 - f_G)*N_GA - V*(K1*C_P - K2*C_L)
if (r_P < 0) then
!    write(nflog,'(a,9e12.3)') 'C_P, C_L, C_P-C_L, r_G, f_G, V, K1, K2, r_P: ',C_P, C_L, C_P-C_L, r_G, f_G, V, K1, K2, r_P
    r_P = 0
    C_P = (r_G*(1 - f_G)*N_GA + V*K2*C_L)/(V*K1)
!    C_P = C_P_min  ! doesn't work
!    r_P = r_G*(1 - f_G)*N_GA - V*(K1*C_P - K2*C_L)
!    write(nflog,'(a,2e12.3)') '=========== r_P < 0, adjust C_P (f_G): ',C_P,f_G
    if (C_P < 0) write(nflog,*) 'C_P < 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
endif
r_O2 = (1 - f_P)*r_P*N_PO
r_L = N_GA*(1 - f_G)*r_G - r_P
r_A = (1-f_G)*r_G*N_GA + (1-f_P)*r_P*N_PA	! production
r_I = f_G*r_G*N_GI + f_P*r_P*N_PI			! production
if (use_glutamine) then
	r_O2 = r_O2 + (1 - f_Gln)*r_Gln*N_GlnO
	r_A = r_A + (1 - f_Gln)*r_Gln*N_GlnA
	r_I = r_I + f_Gln*N_GlnI*r_Gln
	r_I = r_I*MM_Gln
!	if (r_P == 0) write(nflog,'(a,5e12.3)') 'r_P = 0: F, r_G,r_Gln,r_A,C_P: ',F,r_G,r_Gln,r_A,C_P
endif

if (.not.first .and. r_I > r_Iu) then
    write(nflog,'(a,2e12.3)') 'r_I > r_Iu: ',r_I,r_Iu
    ! reduce w until r_I = r_Iu  (crude but foolproof - could be converted to Newton's method)
    dw = 0.001
    do k = 1,100
        r_I_prev = r_I
        w = 1.0 - (k-1)*dw
        f_G = w*f_Gu
        f_P = w*f_Pu
        f_Gln = w*f_Glnu
        
        p = (r_G*(1 - f_G)*N_GA + V*K2*C_L)*(1 - f_P)*N_PO + r_Gln*(1 - f_Gln)*N_GlnO
        if (bills_idea) p = p - FF*z*MM_Gln
        q = V*K1*(1 - f_P)*N_PO
        a = q
        if (bills_idea) then
            b = FF*fPDK + q*Km_P - p
        else
            b = F + q*Km_P - p
        endif
        c = -p*Km_P
        d = sqrt(b*b - 4*a*c) 
        C_P = (-b + d)/(2*a)
        
        r_P = r_G*(1 - f_G)*N_GA - V*(K1*C_P - K2*C_L)
        r_I = f_G*r_G*N_GI + f_P*r_P*N_PI			! production
        if (use_glutamine) then
	        r_I = r_I + f_Gln*N_GlnI*r_Gln
	    endif
!	    write(nflog,'(a,i4,2e12.3)') 'k,r_I,r_Iu: ',k,r_I,r_Iu
	    if (r_I < r_Iu) then
	        w = w - dw*(r_Iu - r_I)/(r_I_prev - r_I)
	        exit
	    endif
	enddo
    f_G = w*f_Gu
    f_P = w*f_Pu
    f_Gln = w*f_Glnu
    
    r_P = r_G*(1 - f_G)*N_GA - V*(K1*C_P - K2*C_L)
    r_I = f_G*r_G*N_GI + f_P*r_P*N_PI			! production
    r_O2 = (1 - f_P)*r_P*N_PO
    r_L = N_GA*(1 - f_G)*r_G - r_P
    r_A = (1-f_G)*r_G*N_GA + (1-f_P)*r_P*N_PA	! production
    if (use_glutamine) then
	    r_O2 = r_O2 + (1 - f_Gln)*r_Gln*N_GlnO
	    r_A = r_A + (1 - f_Gln)*r_Gln*N_GlnA
	    r_I = r_I + f_Gln*N_GlnI*r_Gln
	endif
!    write(nflog,'(a,3e15.6)') 'w,r_I,r_Iu: ',w,r_I,r_Iu
endif

if (dbug) write(nflog,'(a,5e12.3)') 'rates: O2,G,P,A,I: ',r_O2,r_G,r_P,r_A,r_I
if (dbug) write(nflog,*) 'other r_O2: ',F*C_P/(Km_P + C_P)

mp%f_G = f_G
mp%f_P = f_P
mp%f_Gln = f_Gln
mp%A_rate = r_A									! production
mp%I_rate = r_I									! production
mp%P_rate = r_P									! utilisation
mp%O_rate = r_O2								! consumption
mp%Gln_rate = r_Gln								! consumption
mp%L_rate = r_L									! production
mp%C_P = C_P
solved = (from_excel .and. mp%A_rate >= r_Ag)
! Add base rate correction
mp%O_rate = mp%O_rate + O2_baserate
mp%G_rate = mp%G_rate + G_baserate
!write(nflog,'(a,i8,2e12.3)') 'G_rate: ',istep,mp%HIF1,mp%G_rate
if (dbug) write(nflog,*)
first = .false.
end subroutine



!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function solve_C_P(a,b,c) result(x)
real(REAL_KIND) :: a, b, c, x
real(REAL_KIND) :: d

d = b*b - 4*a*c
if (d < 0) then
	write(*,*) 'Error: solve_C_P: a,b,c,d: ',a,b,c,d
!	write(*,'(a,3e12.3)') 'a,b,c: ',a,b,c
!	write(*,'(a,e12.3)') '-b/2a: ',-b/(2*a)
	x = 0
	return
else
	d = sqrt(d)
endif
x = (-b + d)/(2*a)
if (x < 0) then
	write(*,*) 'solve_C_P: x < 0: ',x
	stop
endif
end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function f_MM(C,Km,N) result(v)
real(REAL_KIND) :: C, Km, v
integer :: N

v = C**N/(Km**N + C**N)
end function

!----------------------------------------------------------------------------------
subroutine Set_f_GP(mp,C)
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C(:)

if (fgp_solver == FGP_SOLVER_MAXATP_TANDEM) then
	call Set_f_GP_tandem(mp,C)
!elseif (fgp_solver == FGP_SOLVER_MAXATP_STAGED) then
!	call Set_f_GP_maxATP(mp,C)
!elseif (fgp_solver == FGP_SOLVER_SURVIVAL_STAGED) then
!	call Set_f_GP_survival(mp,C)
endif
end subroutine


!----------------------------------------------------------------------------------
subroutine Set_f_GP_tandem(mp,C)
integer :: ityp
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C(:)
real(REAL_KIND) :: C_O2, C_G, C_L, C_Gln, C_P, r_Atarget
real(REAL_KIND) :: Km_P, fPDK, MM_O2, MM_Gln, r_G, r_Gln, r_P, r_A, r_I, r_O2, f_G, f_P, f_Gln
real(REAL_KIND) :: F, V, K1, K2, a, b, cc, a1, b1, c1, q1, q2, q3, q4, q5
real(REAL_KIND) :: rtol, x0, x1, dx, F0, F1, dFdx, w, z
real(REAL_KIND) :: average_volume = 1.2
integer :: N_O2, N_Gln, k

!write(nflog,'(a,4e12.3)') 'Set_f_GP_tandem: C: ',C(1:4)
if (N1D == 0) then
	! Cex for 3D case - vspheroid
	C_O2 = C(1)
	C_G = C(2)
	C_L = C(3)
	C_Gln = C(4)
else
	! Cex for 1D case - monolayer
	C_O2 = C(1)
	C_G = C(N1D+2)
	C_L = C(2*N1D+3)
	C_Gln = C(3*N1D+4)
endif
solved = .false.

! TESTING=================================================================================
!C_O2 = 0.1
!C_G = 5.0
!C_L = 1.0
!C_Gln = 1.0
!=========================================================================================

V = Vcell_cm3*average_volume		! should be actual cell volume cp%V 
K1 = K_PL
K2 = K_LP
Km_P = Hill_Km_P
fPDK = mp%PDK1
N_O2 = Hill_N_O2
N_Gln = Hill_N_Gln

f_G = f_Gu
f_P = f_Pu
f_Gln = f_Glnu
mp%f_G = f_Gu
mp%f_P = f_Pu
mp%f_Gln = f_Glnu
C_P = mp%C_P
r_Atarget = r_Ag
!r_Atarget = r_Au
!write(nflog,'(a,6e12.3)') 'C_O2,C_G,C_L,C_Gln,C_P,Km_P: ',C_O2,C_G,C_L,C_Gln,C_P,Km_P
call f_metab(mp,C_O2,C_G,C_L,C_Gln)
if (solved) return

!write(nflog,'(a,2e15.6)') 'Set_f_GP_tandem: did f_metab: A_rate,r_Au: ',mp%A_rate,r_Au
if (mp%A_rate < r_Atarget) then      ! need to reduce f_G, f_P, f_Gln
    solved = from_excel
    MM_O2 = f_MM(C_O2,Hill_Km_O2,N_O2)
!    write(nflog,*) 'Hill_Km_Gln,N_Gln: ',Hill_Km_Gln,N_Gln
    if (use_glutamine) then
        MM_Gln = 1  !f_MM(C_Gln,Hill_Km_Gln,N_Gln)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        F = fPDK*O2_maxrate*MM_O2*MM_Gln
    	r_Gln = get_glutamine_rate(C_Gln, fPDK, MM_O2, C_G)
    else
        F = fPDK*O2_maxrate*MM_O2
        r_Gln = 0
    endif
!    write(nflog,'(a,5e12.3)') 'fPDK,O2_maxrate,MM_O2,MM_Gln,f: ',fPDK,O2_maxrate,MM_O2,MM_Gln,f
	r_G = get_glycosis_rate(mp%HIF1,C_G,C_Gln,mp%O_rate)  ! Note: this is the previous O_rate
!	write(nflog,'(a,2e12.3)') 'r_G,r_Gln: ',r_G,r_Gln
    a = r_Atarget - r_G*N_GA - r_Gln*N_GlnA + r_Gln*N_GlnO*N_PA/N_PO
    b = f_Gu*r_G*N_GA + f_Glnu*r_Gln*N_GlnA - f_Glnu*r_Gln*N_GlnO*N_PA/N_PO
    cc = N_PA*f/N_PO
!    write(nflog,'(a,3e12.3)') 'a,b,cc: ',a,b,cc
    a1 = -a*Km_P/b
    b1 = (cc-a)/b
    c1 = Km_P
    w = (a1 + b1*C_P)/(c1 + C_P)
!    write(nflog,'(a,4e12.3)') 'a1,b1,c1,w: ',a1,b1,c1,w
    
    q1 = r_Gln*f_Glnu*N_GlnO + r_G*(f_Pu + f_Gu)*N_PO*N_GA + f_Pu*N_PO*V*K2*C_L
    q2 = -f_Pu*N_PO*V*K1
    q3 = N_PO*V*K1
    q4 = -r_G*f_Pu*f_Gu*N_PO*N_GA
    q5 = -(r_Gln*N_GlnO + r_G*N_PO*N_GA + N_PO*V*K2*C_L)
    
    ! Newtons method
	rtol = 1.0d-6
	x0 = mp%C_P
	dx = max(0.00001,x0/1000)
	do k = 1,100
		w = (a1 + b1*x0)/(c1 + x0)
		w = min(w,1.0)
		F0 = F*x0/(Km_P + x0) + q1*w + q2*w*x0 + q3*x0 + q4*w*w + q5
		x1 = x0 + dx
!		write(nflog,'(a,i4,3e12.3)') 'k,w,x0,x1: ',k,w,x0,x1
		w = (a1 + b1*x1)/(c1 + x1)
		w = min(w,1.0)
		F1 = F*x1/(Km_P + x1) + q1*w + q2*w*x1 + q3*x1 + q4*w*w + q5
		dFdx = (F1 - F0)/dx
		x0 = x0 - F0/dFdx
		if (abs(F0/dFdx) < rtol) then
		    exit
		endif
	enddo
	z = x0
    if (z < 0) then
        write(nflog,*) 'Error in Set_f_GP_tandem: z<0: ',z
        stop
    endif
	w = (a1 + b1*z)/(c1 + z)
!	write(nflog,'(a,2e12.3)') 'Set_f_GP_tandem: z,w: ',z,w
	w = max(w,0.0)
	f_G = w*f_Gu
	f_P = w*f_Pu
	f_Gln = w*f_Glnu
	C_P = z
	r_P = (1-f_G)*r_G*N_GA - V*(K1*C_P - K2*C_L)
	if (r_P < 0) then
	    solved = .false.
!	    write(nflog,'(a,3e12.3)') 'Set_f_GP_tandem: z,w,r_P: ',z,w,r_P
	    ! set r_P = 0, resolve for z,w
	    r_P = 0
	    w = (r_G*N_GA + r_Gln*N_GlnA - r_Atarget)/(f_Gu*r_G*N_GA + f_Glnu*r_Gln*N_GlnA)
	    w = max(w,0.0)
	    z = ((1 - w*f_Gu)*r_G*N_GA + V*K2*C_L)/(V*K1)
    	f_G = w*f_Gu
    	f_P = w*f_Pu
    	f_Gln = w*f_Glnu
!	    write(nflog,'(a,3e12.3)') 'resolved with r_P = 0: z,w: ',z,w
	endif
endif
mp%f_G = f_G
mp%f_P = f_P
mp%f_Gln = f_Gln
if (solved) then
    r_P = r_G*(1 - f_G)*N_GA - V*(K1*C_P - K2*C_L)
    r_I = f_G*r_G*N_GI + f_P*r_P*N_PI			! production
    r_O2 = (1 - f_P)*r_P*N_PO
    r_A = (1-f_G)*r_G*N_GA + (1-f_P)*r_P*N_PA	! production
    if (use_glutamine) then
        r_O2 = r_O2 + (1 - f_Gln)*r_Gln*N_GlnO
        r_A = r_A + (1 - f_Gln)*r_Gln*N_GlnA
        r_I = r_I + f_Gln*N_GlnI*r_Gln
    endif
    mp%P_rate = r_P
    mp%I_rate = r_I
    mp%O_rate = r_O2
    mp%A_rate = r_A
    mp%L_rate = N_GA*(1 - f_G)*r_G - r_P
    mp%C_P = C_P
endif
end subroutine



!----------------------------------------------------------------------------------
! Computes metabolism rate as a fraction of the maximum cell rate
! Use the "soft landing" option for Hill_N = 1 if MM_threshold = 0
!----------------------------------------------------------------------------------
function O2_metab(C) result(metab)
integer :: ichemo
real(REAL_KIND) :: C
real(REAL_KIND) :: metab

ichemo = OXYGEN
if (ichemo == OXYGEN) then
	if (chemo(ichemo)%Hill_N == 2) then
		if (C > 0) then
			metab = C*C/(chemo(ichemo)%MM_C0*chemo(ichemo)%MM_C0 + C*C)
		else
			metab = 0
		endif
	else
!		if (MM_THRESHOLD > 0) then
!			if (C > ODEdiff%C1_soft) then
!				metab = (C-ODEdiff%deltaC_soft)/(chemo(ichemo)%MM_C0 + C - ODEdiff%deltaC_soft)
!			elseif (C > 0) then
!				metab = ODEdiff%k_soft*C*C
!			else
!				metab = 0
!			endif
!		else
			if (C > 0) then
				metab = C/(chemo(ichemo)%MM_C0 + C)
			else
				metab = 0
			endif
!		endif
	endif
endif
end function

!----------------------------------------------------------------------------------
! Computes metabolism rate as a fraction of the maximum cell rate
!----------------------------------------------------------------------------------
function glucose_metab(C) result(metab)
real(REAL_KIND) :: C, metab
real(REAL_KIND) :: Kmin, Kmax, Km1
real(REAL_KIND) :: Vmax1, Vmax2, Km2, n1, n2
real(REAL_KIND) :: fV = 0.6
real(REAL_KIND) :: fK = 0.08
real(REAL_KIND) :: fboost = 2
real(REAL_KIND) :: Cboost = 0.1
logical :: variable_Km = .false.
logical :: double_Km = .false.
logical :: use_boost = .false.

if (C == 0) then
	metab = 0
	return
endif
if (use_boost) then

elseif (double_Km) then
	Km1 = Hill_Km_G
	Km2 = fK*Km1
	n1 = Hill_N_G
	n2 = 1
	metab = fV*C**n1/(Km1**n1 + C**n1) + (1 - fV)*C**n2/(Km2**n2 + C**n2)
elseif (variable_Km) then
	Kmax = Hill_Km_G	! These are completely arbitrary values
	Kmin = Kmax/15
	Km1 = 1*Kmin
	metab = C*(Km1 + C)/(Kmin*Km1 + Kmax*C + C*(Km1 + C))
else
	metab = C**Hill_N_G /(C**Hill_N_G + Hill_Km_G**Hill_N_G)
endif
end function

#if EXCEL
!--------------------------------------------------------------------------------------
! Determine real roots r(:) of the cubic equation:
! x^3 + a.x^2 + b.x + c = 0
! If there is one real root, n=1 and the root is r(1)
! If there are three distinct real roots, n=3 and the roots are r(1), r(2), r(3)
! If there is a repeated root, n=2 and the single root is r(1), the repeated root is r(2)
!--------------------------------------------------------------------------------------
subroutine cubic_roots(a, b, c, r, n)
real(REAL_KIND) :: a, b, c, r(3)
integer :: n
real(REAL_KIND) :: QQ, RR, theta, R2, Q3, AA, BB

QQ = (a*a - 3*b)/9
RR = (2*a*a*a - 9*a*b + 27*c)/54
Q3 = QQ*QQ*QQ
R2 = RR*RR
if (R2 < Q3) then
	n = 3
	theta = acos(RR/sqrt(Q3))
	r(1) = -2*sqrt(QQ)*cos(theta/3) - a/3
	r(2) = -2*sqrt(QQ)*cos((theta+2*PI)/3) - a/3
	r(3) = -2*sqrt(QQ)*cos((theta-2*PI)/3) - a/3
else
	n = 1
	AA = -sign(1.d0,RR)*(abs(RR) + sqrt(R2 - Q3))**(1.d0/3.d0)
	if (AA == 0) then
		BB = 0
	else
		BB = QQ/AA
	endif
	r(1) = AA + BB - a/3
endif
end subroutine
#endif
!=====================================================================================================

!#include 'nitrogen.f90'

!=====================================================================================================
!=====================================================================================================
subroutine get_unconstrained_rates
real(REAL_KIND) :: V
type(metabolism_type), target :: metab
type(metabolism_type), pointer :: mp

V = Vcell_cm3*average_volume		! should be actual cell volume cp%V
 
r_Gu = G_maxrate*glucose_metab(C_G_norm)
r_Pu = (O2_maxrate - r_Glnu*(1 - f_Glnu)*N_GlnO)/((1 - f_Pu)*N_PO)
!a0 = r_G*N_GA + ps%p*(1 - f_Glnu)*N_GlnA - ps%h*ps%p*f_Glnu*N_GlnI - r_Ag
!b0 = -r_G*f_Gu*(N_GA + ps%h*N_GI)
!c0 = ps%q*(1 - f_Glnu)*N_GlnA - ps%h*ps%q*f_Glnu*N_GlnI
!d0 = f_Pu*(N_PA + ps%h*N_PI)

ps%a1 = ps%p*f_Glnu*N_GlnI
ps%b1 = r_Gu*f_Gu*N_GI
ps%c1 = -ps%q*f_Glnu*N_GlnI
ps%d1 = f_Pu*N_PI

r_Iu = r_Gu*f_Gu*N_GI + r_Pu*f_Pu*N_PI + r_Glnu*f_Glnu*N_GlnI
r_Au = r_Gu*(1 - f_Gu)*N_GA + r_Pu*(1 - f_Pu)*N_PA + r_Glnu*(1 - f_Glnu)*N_GlnA
r_Ou = r_Pu*(1 - f_Pu)*N_PO + r_Glnu*(1 - f_Glnu)*N_GlnO
C_Pu = ((1 - f_Gu)*r_Gu*N_GP + V*K_LP*C_L_norm - r_Pu)/(V*K_PL)
r_Lu = (1 - f_Gu)*r_Gu*N_GP - r_Pu
write(nflog,'(a,7e12.3)') 'unconstrained rates: A,I,G,Gln,P,L,O: ',r_Au,r_Iu,r_Gu,r_Glnu,r_Pu,r_Lu,r_Ou
write(nflog,'(a,e12.3)') 'O2_maxrate: ',O2_maxrate
r_Ag = f_ATPg*r_Au
r_As = f_ATPs*r_Au
write(nflog,'(a,e12.3)') 'r_Ag: ',r_Ag
return

! Now check by solving
write(nflog,*) 'Checking unconstrained rates'
mp => metab
mp%HIF1 = get_HIF1steadystate(C_O2_norm)
call analyticSetPDK1(mp%HIF1, mp%PDK1, 1.0d10)
mp%O_rate = O2_maxrate
call f_metab_nitrogen(mp, C_O2_norm, C_G_norm, C_L_norm, C_Gln_norm)
write(nflog,*) 'check ended'
close(nflog)

end subroutine

!==================================================================================================
!==================================================================================================
subroutine set_param_set(r_G,C_Gln)
real(REAL_KIND) :: r_G, C_Gln
real(REAL_KIND) :: V, MM_Gln, Km_Gln, r1, r2
integer :: N_Gln
logical :: use_wxfGlnu = .true.

V = Vcell_cm3*average_volume
N_Gln = 1
Km_Gln = Hill_Km_Gln
MM_Gln = f_MM(C_Gln,Km_Gln,N_Gln)

ps%h = (r_Au - r_Ag)/r_Iu
r1 = min(r_Glnu,Gln_maxrate*MM_Gln)
r2 = max(r1,Gln_maxrate*MM_Gln)
ps%p = r2
ps%q = (r2 - r1)/r_Pu

ps%a0 = r_G*N_GA + ps%p*(1 - f_Glnu)*N_GlnA - ps%h*ps%p*f_Glnu*N_GlnI - r_Ag
ps%b0 = -r_G*f_Gu*(N_GA + ps%h*N_GI)
ps%c0 = ps%q*(1 - f_Glnu)*N_GlnA - ps%h*ps%q*f_Glnu*N_GlnI - N_PA
ps%d0 = f_Pu*(N_PA + ps%h*N_PI)

if (use_wxfGlnu) then
    ps%a1 = 0
    ps%b1 = r_G*f_Gu*N_GI + ps%p*f_Glnu*N_GlnI
    ps%c1 = 0
    ps%d1 = f_Pu*N_PI - ps%q*f_Glnu*N_GlnI

    ps%a2 = r_G*N_GA + ps%p*N_GlnA
    ps%b2 = -(r_G*f_Gu*N_GA + ps%p*f_Glnu*N_GlnA)
    ps%c2 = N_PA - ps%q*N_GlnA
    ps%d2 = ps%q*f_Glnu*N_GlnA - f_Pu*N_PA
    
    ps%a3 = ps%p*N_GlnO
    ps%b3 = -ps%p*f_Glnu*N_GlnO
    ps%c3 = N_PO - ps%q*N_GlnO
    ps%d3 = ps%q*f_Glnu*N_GlnO - f_Pu*N_PO
else
    ps%a1 = ps%p*f_Glnu*N_GlnI
    ps%b1 = r_G*f_Gu*N_GI
    ps%c1 = -ps%q*f_Glnu*N_GlnI
    ps%d1 = f_Pu*N_PI

    ps%a2 = r_G*N_GA + ps%p*(1 - f_Glnu)*N_GlnA
    ps%b2 = -r_G*f_Gu*N_GA
    ps%c2 = -ps%q*(1 - f_Glnu)*N_GlnA + N_PA
    ps%d2 = -f_Pu*N_PA
endif

end subroutine

!==================================================================================================
!==================================================================================================
subroutine f_metab_nitrogen(mp, C_O2, C_G, C_L, C_Gln)
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C_O2, C_G, C_L, C_Gln
real(REAL_KIND) :: w, a, b, c, d, aa, bb, cc, dd, w1, w2, f_NG, q, V, tol, r_P0, r_P1, r_PA, a1, b1, r_I1, r_Gln1
real(REAL_KIND) :: Km_O2, Km_Gln, MM_O2, MM_Gln, L1, L2, C_P, r_G, r_P, r_O2, r_Gln, r_A, r_I, r_L
real(REAL_KIND) :: clean, check
integer :: N_O2, N_Gln, iw, Nw
real(REAL_KIND) :: dw, r_Imax, w_max, r1, r2, a0, b0, c0, d0, r_O2limit, alfa
logical :: checking = .false.
logical :: done, dbug
clean = -1

dbug = .false.
!dbug = (C_G < 0.05)
r_G = get_glycosis_rate(mp%HIF1,C_G,C_Gln,mp%O_rate)  ! Note: this is the previous O_rate
call set_param_set(r_G,C_Gln)
if (dbug) write(nflog,'(a,3e12.3)') 'ps%h, ps%p, ps%q: ',ps%h, ps%p, ps%q

N_O2 = Hill_N_O2
Km_O2 = Hill_Km_O2
N_Gln = 1
Km_Gln = Hill_Km_Gln
MM_O2 = f_MM(C_O2,Km_O2,N_O2)
L1 = mp%PDK1*O2_maxrate*MM_O2
MM_Gln = f_MM(C_Gln,Km_Gln,N_Gln)
L2 = Gln_maxrate*MM_Gln
r1 = min(r_Glnu,Gln_maxrate*MM_Gln)
r2 = max(r1,Gln_maxrate*MM_Gln)
if (dbug) then
!    write(nflog,'(a,5e12.3)') 'r_Glnu,MM_Gln,r1,r2: ',r_Glnu,MM_Gln,r1,r2,Gln_maxrate*MM_Gln
    a0 = ps%h*ps%a1 - ps%a2 + r_Ag
    b0 = ps%h*ps%b1 - ps%b2
    c0 = ps%c2 - ps%h*ps%c1
    d0 = ps%d2 - ps%h*ps%d1
!    write(nflog,'(a,4e12.3)') 'ps% a0,b0,c0,d0: ',ps%a0,ps%b0,ps%c0,ps%d0
!    write(nflog,'(a,4e12.3)') '    a0,b0,c0,d0: ',a0,b0,c0,d0
    ! For w = 1
    r_P = (ps%h*(ps%a1 + ps%b1) + r_Ag - (ps%a2 + ps%b2))/(ps%c2 + ps%d2 - ps%h*(ps%c1 + ps%d1))
!    write(nflog,'(a,2e12.3)') 'r_P for w=1, r_G: ',r_P,r_G
    ! Get w from r_Pu
    w = (ps%a0 - ps%c0*r_Pu)/(ps%d0*r_Pu - ps%b0)
!    write(nflog,'(a,e12.3)') 'w from r_Pu: ',w
    ! unconstrained rates
    r_A = r_Gu*(1-f_Gu)*N_GA + r_Pu*(1-f_Pu)*N_PA + r_Glnu*(1-f_Glnu)*N_GlnA    ! good
    r_I = r_Gu*f_Gu*N_GI + r_Pu*f_Pu*N_PI + r_Glnu*f_Glnu*N_GlnI
!    write(nflog,'(a,4e12.3)') 'from unconstrained: r_A,r_Ag,r_I,h: ',r_A,r_Ag,r_I,(r_A-r_Ag)/r_I
!    write(nflog,'(a,2e12.3)') 'p-q.r_P: ',ps%p - ps%q*r_Pu,r_Glnu
    r_I = ps%a1 + ps%b1 + r_Pu*(ps%c1 + ps%d1)
    r_A = ps%a2 + ps%b2 + r_Pu*(ps%c2 + ps%d2)  ! bad
!    write(nflog,'(a,2e12.3)') 'from w=1, r_Pu: r_I,r_A: ',r_I,r_A
!    write(nflog,'(a,2e12.3)') 'a2: ',ps%a2,r_G*N_GA + ps%p*(1-f_Glnu)*N_GlnA
!    write(nflog,'(a,2e12.3)') 'b2: ',ps%b2,-r_G*f_Gu*N_GA
!    write(nflog,'(a,2e12.3)') 'c2: ',ps%c2,-ps%q*(1-f_Glnu)*N_GlnA + N_PA
!    write(nflog,'(a,2e12.3)') 'd2: ',ps%d2,-f_Pu*N_PA
endif
V = Vcell_cm3*average_volume
r_P = 0
r_Gln = 0
r_I = 0
r_A = 0
r_O2 = 0
r_L = 0
C_P = 0

r_Imax = 0
w_max = -1
Nw = 10
dw = 1.0/Nw
do iw = 1,Nw+1
    w = (iw-1)*dw
    r_P = (ps%a0 + ps%b0*w)/(ps%c0 + ps%d0*w)
    r_I = ps%a1 + ps%b1*w + ps%c1*r_P + ps%d1*w*r_P
    if (dbug) write(nflog,'(a,f6.3,2e12.3)') 'w,r_P,r_I: ',w,r_P,r_I
    if (r_I > r_Imax) then
        r_Imax = r_I
        w_max = w
    endif
enddo
w = w_max
r_P = (ps%a0 + ps%b0*w)/(ps%c0 + ps%d0*w)
C_P = ((1 - w*f_Gu)*r_G*N_GP + V*K_LP*C_L - r_P)/(V*K_PL)
!write(nflog,'(a,f6.3,2e12.3)') 'w,r_P,C_P: ',w,r_P,C_P
if (C_P < 0.0) then
    C_P = 0
    a = ps%d0*r_G*f_Gu*N_GP
    b = ps%b0 + ps%c0*r_G*f_Gu*N_GP - ps%d0*(r_G*N_GP + V*K_LP*C_L)
    c = ps%a0 - ps%c0*(r_G*N_GP + V*K_LP*C_L)
    d = sqrt(b*b - 4*a*c)
    w1 = (-b + d)/(2*a)
    w2 = (-b - d)/(2*a)
    w = max(w2,0.0)
    r_P = r_G*(1 - w*f_Gu)*N_GP + V*K_LP*C_L
!    write(nflog,'(a,3f8.4,2e12.3)') 'w1,w2,w,r_P,C_P: ',w1,w2,w,r_P,C_P
endif    
r_Gln = ps%p - ps%q*r_P
r_I = ps%a1 + ps%b1*w + ps%c1*r_P + ps%d1*w*r_P
r_A = ps%a2 + ps%b2*w + ps%c2*r_P + ps%d2*w*r_P
r_O2 = r_P*(1 - w*f_Pu)*N_PO + r_Gln*(1 - w*f_Glnu)*N_GlnO
r_O2limit = mp%PDK1*O2_maxrate*MM_O2
if (r_O2 > r_O2limit) then
#if 0
    ! solve for w, r_P with r_O2 = r_O2limit
    a = ps%b3*ps%d0 + ps%d3*ps%b0
    b = ps%b3*ps%c0 + ps%a3*ps%d0 + ps%c3*ps%b0 + ps%d3*ps%a0 - r_o2limit*ps%d0
    c = ps%a3*ps%c0 + ps%c3*ps%a0 - r_O2limit*ps%c0
    d = sqrt(b*b - 4*a*c)
    w1 = (-b + d)/(2*a)
    w2 = (-b - d)/(2*a)
!    write(nflog,'(a,3f8.4,2e12.3)') 'w1,w2: ',w1,w2
    w = w2
    if (w > 0) then
        r_P = (ps%a0 + ps%b0*w)/(ps%c0 + ps%d0*w)
        r_O2 = ps%a3 + ps%b3*w + ps%c3*r_P + ps%d3*w*r_P
        write(nflog,'(a,f8.4,3e12.3)') 'w,r_P,r_O2,a0/c0: ',w,r_P,r_O2,ps%a0/ps%c0
    else
        write(nflog,'(a,f8.4)') 'for r_O2 bad w: ',w
        w = 0
    endif
#endif
    alfa = r_O2limit/r_O2
    r_O2 = r_O2limit
    r_P = alfa*r_P
    r_Gln = alfa*r_Gln
!    w = alfa*w
    w = (r_G*N_GA + r_P*N_PA + r_Gln*N_GlnA - r_Ag)/ &
        (ps%h*(r_G*f_Gu*N_GI + r_P*f_Pu*N_PI + r_Gln*f_Glnu*N_GlnI) + (r_G*f_Gu*N_GA + r_P*f_Pu*N_PA + r_Gln*f_Glnu*N_GlnA))
    r_A = r_G*(1 - w*f_Gu)*N_GA + r_P*(1 - w*f_Pu)*N_PA + r_Gln*(1 - w*f_Glnu)*N_GlnA
    r_I = r_G*w*f_Gu*N_GI + r_P*w*f_Pu*N_PI + r_Gln*w*f_Glnu*N_GlnI
    C_P = ((1 - w*f_Gu)*r_G*N_GP + V*K_LP*C_L - r_P)/(V*K_PL)
    write(nflog,'(a,f8.4,3e12.3)') 'for r_O2: w,r_P,r_Gln,C_P: ',w,r_P,r_Gln,C_P
    if (C_P < 0) then
        write(nflog,*) 'C_P < 0 !!!!!!!!!!!!!!!!!!!!!!!!!'
    endif
    C_P = max(C_P,0.0)
endif
r_L = r_G*(1 - w*f_Gu)*N_GP - r_P
!write(nflog,'(a,8e12.3)') 'ps% a0,b0,c0,d0,p,q,r1,r2: ',ps%a0,ps%b0,ps%c0,ps%d0,ps%p,ps%q,r1,r2
done = .true.
if (done) then
    mp%P_rate = r_P
    mp%G_rate = r_G
    mp%Gln_rate = r_Gln
    mp%I_rate = r_I
    mp%A_rate = r_A
    mp%O_rate = r_O2
    mp%L_rate = r_L
    mp%C_P = C_P
    mp%f_G = w*f_Gu
    mp%f_P = w*f_Pu
    mp%recalcable = clean
endif
end subroutine


#if 0
!==================================================================================================
!==================================================================================================
subroutine get_unconstrained_rates
real(REAL_KIND) :: f_NG, V, qu, MM_Gln, Km_Gln, L1
integer :: N_Gln
type(metabolism_type), target :: metab
type(metabolism_type), pointer :: mp

V = Vcell_cm3*average_volume		! should be actual cell volume cp%V
f_NG = f_N/f_Glnu
N_Gln = 1
Km_Gln = Hill_Km_Gln
MM_Gln = f_MM(C_Gln_norm,Km_Gln,N_Gln)
!L2 = Gln_maxrate*MM_Gln
 
write(nflog,'(a,3f6.3)') 'f_N, f_Glnu, f_NG: ',f_N, f_Glnu, f_NG
r_Gu = G_maxrate*glucose_metab(C_G_norm)
qu = (f_NG*(1 - f_Glnu)*N_GlnO)/(1 - f_N*N_GlnI)
r_Pu = (O2_maxrate - qu*r_Gu*f_Gu*N_GI)/((1 - f_Pu)*N_PO + qu*f_Pu*N_PI)
r_Iu = (r_Gu*f_Gu*N_GI + r_Pu*f_Pu*N_PI)/(1 - f_N*N_GlnI) 
r_Glnu = f_NG*r_Iu
L1 = O2_maxrate

!if (r_Glnu > L2) then
!    write(nflog,*) 'ERROR: unconstrained Gln rate exceeds constraint L2'    ! Can this be fixed?
!    stop
!endif
! Instead, adjust Gln_maxrate to match r_Glnu
Gln_maxrate = r_Glnu/MM_Gln
write(nflog,'(a,2e12.3)') 'r_Glnu, L1: ',r_Glnu,L1

r_Au = r_Gu*(1 - f_Gu)*N_GA + r_Pu*(1 - f_Pu)*N_PA + r_Glnu*(1 - f_Glnu)*N_GlnA
r_Ou = r_Pu*(1 - f_Pu)*N_PO + r_Glnu*(1 - f_Glnu)*N_GlnO
C_Pu = ((1 - f_Gu)*r_Gu*N_GP + V*K_LP*C_L_norm - r_Pu)/(V*K_PL)
r_Lu = (1 - f_Gu)*r_Gu*N_GP - r_Pu
write(nflog,'(a,7e12.3)') 'unconstrained rates: A,I,G,Gln,P,L,O: ',r_Au,r_Iu,r_Gu,r_Glnu,r_Pu,r_Lu,r_Ou
r_Ag = f_ATPg*r_Au
r_As = f_ATPs*r_Au
return

! Now check by solving
write(nflog,*) 'Checking unconstrained rates'
mp => metab
mp%HIF1 = get_HIF1steadystate(C_O2_norm)
call analyticSetPDK1(mp%HIF1, mp%PDK1, 1.0d10)
mp%O_rate = O2_maxrate
call f_metab_nitrogen(mp, C_O2_norm, C_G_norm, C_L_norm, C_Gln_norm)
write(nflog,*) 'check ended'
close(nflog)
end subroutine

!==================================================================================================
!==================================================================================================
subroutine set_param_set(r_G,C_L)
real(REAL_KIND) :: r_G, C_L
real(REAL_KIND) :: V, f_NG, q, fGlnA
!real(REAL_KIND) :: r_P1, r_P2, C_P, w

V = Vcell_cm3*average_volume
f_NG = f_N/f_Glnu
!q = 1/(1 - f_NG*f_Glnu*N_GlnI)
fGlnA = 1 - f_Glnu

ps%h = (r_Au - r_Ag)/r_Iu
ps%a0 = r_G*N_GP + V*K_LP*C_L
ps%b0 = -r_G*f_Gu*N_GP
ps%c0 = -V*K_PL
!a1 = (f_Pu*N_PI)/(1 - f_N*N_GlnI)
ps%a1 = f_Pu*N_PI/(1 - f_N*N_GlnI)
ps%b1 = r_G*f_Gu*N_GI/(1 - f_N*N_GlnI)
ps%a2 = f_NG*ps%a1
ps%b2 = f_NG*ps%b1
ps%a3 = N_PA
ps%b3 = ps%a2*fGlnA*N_GlnA - f_Pu*N_PA
ps%c3 = r_G*N_GA
ps%d3 = ps%b2*fGlnA*N_GlnA - r_G*f_Gu*N_GA
ps%a4 = N_PO
ps%b4 = ps%a2*fGlnA*N_GlnO - f_Pu*N_PO
ps%c4 = ps%b2*fGlnA*N_GlnO

!C_Pu = (r_Pu - ps%a0 - ps%b0)/ps%c0
!C_P = ((1 - f_Gu)*r_G*N_GP + V*K_LP*C_L - r_Pu)/(V*K_PL)
!write(nflog,'(a,4e12.3)') 'set_param_set: V, C_Pu, C_P: ',V,C_Pu,C_P
!r_P1 = (1 - f_Gu)*r_Gu*N_GP - V*(K_PL*C_Pu - K_LP*C_L)
!r_P2 = (1 - f_Gu)*r_Gu*N_GP - V*(K_PL*C_P - K_LP*C_L)
!write(nflog,'(a,2e12.3)') 'computed r_P1, r_P2: ',r_P1,r_P2

! check value of g(w,x)
!w = 1
!r_P1 = (ps%c3 - f_ATPg*r_Au + w*(ps%d3 - ps%h*ps%b1))/(w*(ps%h*ps%a1 - ps%b3) - ps%a3)
!write(nflog,'(a,2e12.3)') 'computed r_P1: ',r_P1
!write(nflog,'(a,2e12.3)') 'Numerator is a + bw, where a, b: ',ps%c3 - f_ATPg*r_Au,ps%d3 - ps%h*ps%b1
!write(nflog,'(a,2e12.3)') 'Denominator is cw - d, where c, d: ',ps%h*ps%a1 - ps%b3,ps%a3

end subroutine

!==================================================================================================
! Let r_P = f(w) = (a + bw)/(c + dw)
!==================================================================================================
subroutine f_metab_nitrogen(mp, C_O2, C_G, C_L, C_Gln)
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C_O2, C_G, C_L, C_Gln
real(REAL_KIND) :: w, a, b, c, d, aa, bb, cc, dd, w1, w2, f_NG, q, V, tol, r_P0, r_P1, r_PA, a1, b1, r_I1, r_Gln1
real(REAL_KIND) :: Km_O2, Km_Gln, MM_O2, MM_Gln, L1, L2, C_P, r_G, r_P, r_O2, r_Gln, r_A, r_I, r_L
real(REAL_KIND) :: clean, check
integer :: N_O2, N_Gln, iw
logical :: checking = .false.

logical :: done
clean = 1

r_G = get_glycosis_rate(mp%HIF1,C_G,C_Gln,mp%O_rate)  ! Note: this is the previous O_rate
call set_param_set(r_G,C_L)

N_O2 = Hill_N_O2
Km_O2 = Hill_Km_O2
N_Gln = 1
Km_Gln = Hill_Km_Gln
MM_O2 = f_MM(C_O2,Km_O2,N_O2)
L1 = mp%PDK1*O2_maxrate*MM_O2
MM_Gln = f_MM(C_Gln,Km_Gln,N_Gln)
L2 = Gln_maxrate*MM_Gln

V = Vcell_cm3*average_volume
f_NG = f_N/f_Glnu

if (checking) then
! checking code
f_NG = f_N/f_Glnu
q = f_NG*(1 - f_Glnu)*N_GlnO/(1 - f_N*N_GlnI)
w = 1
r_P = (r_Ou - q*r_G*w*f_Gu*N_GI)/((1 - w*f_Pu)*N_PO + q*w*f_Pu*N_PI)
write(nflog,'(a,e12.3)') 'r_P with w=1: ',r_P
C_P = (r_P - ps%a0 - ps%b0*w)/ps%c0
write(nflog,'(a,2e12.3)') 'C_P,C_Pu: ',C_P,C_Pu
r_Gln = w*(r_P*ps%a2 + ps%b2)
r_I = r_P*ps%a1*w + ps%b1*w
!r_I = w*(r_G*f_Gu*N_GI + r_P*f_Pu*N_PI)/(1 - f_N*N_GlnI) 
a1 = (f_Pu*N_PI)/(1 - f_N*N_GlnI)
b1 = (r_G*f_Gu*N_GI)/(1 - f_N*N_GlnI)
r_I1 = w*(r_P*a1 + b1)
r_Gln1 = f_NG*r_I1
r_P1 = r_P
write(nflog,'(a,6e12.3)') 'a1,ps%a1,b1,ps%b1,r_I1,r_Gln1: ',a1,ps%a1,b1,ps%b1,r_I1,r_Gln1
r_O2 = r_P*(ps%a4 + ps%b4*w) + ps%c4*w 
r_A = r_P*(ps%a3 + ps%b3*w) + ps%c3 + ps%d3*w
r_L = r_G*(1 - w*f_Gu)*N_GP - r_P
write(nflog,'(a,5e12.3)') 'r_Gln,r_I,r_O2,r_A,r_L: ',r_Gln,r_I,r_O2,r_A,r_L
write(nflog,'(a,23e12.3)') 'h,(r_A-r_Ag)/r_I,(r_Au-r_Ag)/r_Iu: ',ps%h,(r_A-r_Ag)/r_I,(r_Au-r_Ag)/r_Iu
write(nflog,'(a,2e12.3)') '(r_A-r_Ag)/(r_Au-r_Ag),r_I/r_Iu: ',(r_A-r_Ag)/(r_Au-r_Ag),r_I/r_Iu
endif
done = .false.
a = ps%c3 - r_Ag
b = ps%d3 - ps%h*ps%b1
c = - ps%a3
d = ps%h*ps%a1 - ps%b3

! Check for solution to eqtn (5): r_O2 = L1
! aa.w^2 + bb.w + cc = 0
aa = b*ps%b4 + d*ps%c4
bb = b*ps%a4 + a*ps%b4 + c*ps%c4 - d*L1
cc = a*ps%a4 - c*L1
dd = sqrt(bb*bb - 4*aa*cc)
w1 = min(1.0,(-bb + dd)/(2*aa))
w2 = min(1.0,(-bb - dd)/(2*aa))
!write(nflog,'(a,2e12.3)') 'L1 solution for w: ',(-bb + dd)/(2*aa),(-bb - dd)/(2*aa)
w = w2
! Check eqtn (6): is r_Gln < L2?
r_P = (a + b*w)/(c + d*w)

!w = 1 - 0.0001       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!r_P = r_P1

if (r_P < 0) then
!    write(nflog,'(a,e12.3)') 'with w=1, r_P: ',r_P
    clean = -1
    r_P = 0
    q = r_G*f_Gu*N_GI/(1 - f_N*N_GlnI)
    w = min(L1/(q*f_NG*(1 - f_Glnu)*N_GlnO),L2/(q*f_NG))
!    write(nflog,'(a,f6.3)') 'In L1: set r_P=0 ==================> w: ',w
endif

C_P = (r_P - ps%a0 - ps%b0*w)/ps%c0
!r_P0 = ps%a0 + ps%b0 + ps%c0*C_Pu
!r_P1 = (a + b)/(c + d)  ! with w = 1
!r_PA = (r_Au - ps%c3 - ps%d3)/(ps%a3 + ps%b3)   ! with r_A = r_Au, w = 1
!write(nflog,'(a,7e12.3)') 'r_P,r_P0,r_P1,r_PA,C_L,C_P,C_Pu: ',r_P,r_P0,r_P1,r_PA,C_L,C_P,C_Pu
if (C_P < 0) then
    clean = -1
    w = zero_C_P(a,b,c,d)
    C_P = 0
    r_P = r_G*(1 - w*f_Gu)*N_GP + V*K_LP*C_L
!    write(nflog,'(a,f8.3,e12.3)') 'In L1, C_P < 0, set w: ',w,r_P
endif
r_Gln = w*(r_P*ps%a2 + ps%b2)
if (w < 0) then
    clean = -1
    w = 0
    r_I = 0
    r_Gln = 0
    r_O2 = L1
    r_P = r_O2/N_PO
    C_P = (r_G*N_GP + V*K_LP*C_L - r_P)/(V*K_PL)
    r_L = r_G*N_GP - r_P
    r_A = r_G*N_GA + r_P*N_PA
    done = .true.
!    write(nflog,'(a,2e12.3)') 'w < 0: r_P, C_P: ',r_P,C_P
else
    if (r_Gln <= L2) then
!        write(nflog,'(a,f8.3,2x,2e12.3)') 'L1 solution for w: ',w,r_Gln,L2
        done = .true.
        r_O2 = r_P*(ps%a4 + ps%b4*w) + ps%c4*w
        r_I = r_P*ps%a1*w + ps%b1*w
        r_A = r_P*(ps%a3 + ps%b3*w) + ps%c3 + ps%d3*w
        r_L = r_G*(1 - w*f_Gu)*N_GP - r_P
    else
!        write(nflog,'(a,2f8.3,2x,2e15.6)') 'L2 violated: ',w1,w2,r_Gln,L2
    endif
endif
if (.not.done) then
    ! Check for solution to eqtn (6): r_Gln = L2
    aa = b*ps%a2 + d*ps%b2
    bb = a*ps%a2 + c*ps%b2 - d*L2
    cc = -c*L2
    dd = sqrt(bb*bb - 4*aa*cc)
    w1 = (-bb + dd)/(2*aa)
    w2 = (-bb - dd)/(2*aa)
!    write(nflog,'(a,2f8.3)') 'w1,w2: ',w1,w2
    w = min(1.0,w2)
!    w = max(0.0,w)
!    if (w < 0) write(nflog,*) 'In L2: w < 0'
    ! Check eqtn (5): is r_O2 < L1?
    r_P = (a + b*w)/(c + d*w)
    if (r_P < 0) then
!        write(nflog,'(a,f6.3)') 'In L2: set r_P=0, ====================================> w: ',-a/b
!        do iw = 1,11
!            w = (iw-1)*0.1
!            r_P = (a + b*w)/(c + d*w)
!            C_P = (r_P - ps%a0 - ps%b0*w)/ps%c0
!            r_O2 = r_P*(ps%a4 + ps%b4*w) + ps%c4*w
!            r_Gln = w*(r_P*ps%a2 + ps%b2)
!            write(nflog,'(i2,f6.1,f6.3,e12.3,2x,2e12.3,2x,2e12.3)') iw,w,C_P,r_P,r_O2,L1,r_Gln,L2
!        enddo
        clean = -1
        r_P = 0      
        q = r_G*f_Gu*N_GI/(1 - f_N*N_GlnI)
        w = min(L1/(q*f_NG*(1 - f_Glnu)*N_GlnO),L2/(q*f_NG))
        r_I = q*w
        r_Gln = f_NG*r_I
 !       write(nflog,'(a,3e12.3)') 'r_I parts: ',r_G*(1 - w*f_Gu)*N_GI, r_Gln*(1 - f_Glnu)*N_GlnI
    else
        r_Gln = w*(r_P*ps%a2 + ps%b2)
        r_I = r_P*ps%a1*w + ps%b1*w
    endif
    C_P = (r_P - ps%a0 - ps%b0*w)/ps%c0
    if (C_P < 0) then
!        write(nflog,*) 'In L2, C_P < 0, set w = 1'
        clean = -1
        w = 1
        C_P = 0
        r_P = r_G*(1 - w*f_Gu)*N_GP + V*K_LP*C_L
        r_Gln = w*(r_P*ps%a2 + ps%b2)
        r_I = r_P*ps%a1*w + ps%b1*w
    endif
    r_O2 = r_P*(ps%a4 + ps%b4*w) + ps%c4*w 
    r_A = r_P*(ps%a3 + ps%b3*w) + ps%c3 + ps%d3*w
    r_L = r_G*(1 - w*f_Gu)*N_GP - r_P
    if (r_O2 <= L1) then
!        write(nflog,'(a,f8.3,2x,2e12.3)') 'L1 solution for w: ',w,r_O2,L1
        done = .true.
    else
!        write(nflog,'(a,2f8.3,2x,2e12.3)') 'L1 violated: ',w,r_O2,L1
    endif
endif
if (done) then
    mp%P_rate = r_P
    mp%G_rate = r_G
    mp%Gln_rate = r_Gln
    mp%I_rate = r_I
    mp%A_rate = r_A
    mp%O_rate = r_O2
    mp%L_rate = r_L
    mp%C_P = C_P
    mp%f_G = w*f_Gu
    mp%f_P = w*f_Pu
    mp%recalcable = clean
!    return
    ! Check that clean really is recalcable
    tol = 1.0e-4
    check = 1
    w = mp%f_G/f_Gu
    r_P = (a + b*w)/(c + d*w)
    if (abs((r_P-mp%P_rate)/r_P) > tol) check = -1
    r_O2 = r_P*(ps%a4 + ps%b4*w) + ps%c4*w
    if (abs((r_O2-mp%O_rate)/r_O2) > tol) check = -1
    r_I = r_P*ps%a1*w + ps%b1*w
    if (abs((r_I-mp%I_rate)/r_I) > tol) check = -1
    r_A = r_P*(ps%a3 + ps%b3*w) + ps%c3 + ps%d3*w
    if (abs((r_A-mp%A_rate)/r_A) > tol) check = -1
    r_L = r_G*(1 - w*f_Gu)*N_GP - r_P
    if (abs((r_L-mp%L_rate)/r_L) > tol) check = -1
    if (check < 0) then
        write(nflog,'(a,f3.0,5e15.6)') 'Bad check: ',check,r_P,r_O2,r_I,r_A,r_L
        write(nflog,'(a,f3.0,5e15.6)') '           ',clean,mp%P_rate,mp%O_rate,mp%I_rate,mp%A_rate,mp%L_rate
    endif
    return
endif
write(nflog,*) 'What now???????????????????????????????????????????????????????????????'
end subroutine


!==================================================================================================
!==================================================================================================
subroutine f_metab_recalc(mp, C_O2, C_G, C_L, C_Gln)
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C_O2, C_G, C_L, C_Gln
real(REAL_KIND) :: w, a, b, c, d, aa, bb, cc, dd, f_NG, q, V, w1, w2
real(REAL_KIND) :: Km_O2, Km_Gln, MM_O2, MM_Gln, L1, L2, C_P, r_G, r_P, r_O2, r_Gln, r_A, r_I, r_L
integer :: N_O2, N_Gln, iw
real(REAL_KIND) :: clean
logical :: done

done = .false.
clean = 1
a = ps%c3 - r_Ag
b = ps%d3 - ps%h*ps%b1
c = - ps%a3
d = ps%h*ps%a1 - ps%b3
N_O2 = Hill_N_O2
Km_O2 = Hill_Km_O2
N_Gln = 1
Km_Gln = Hill_Km_Gln
MM_O2 = f_MM(C_O2,Km_O2,N_O2)
L1 = mp%PDK1*O2_maxrate*MM_O2
MM_Gln = f_MM(C_Gln,Km_Gln,N_Gln)
L2 = Gln_maxrate*MM_Gln

V = Vcell_cm3*average_volume
f_NG = f_N/f_Glnu

r_G = get_glycosis_rate(mp%HIF1,C_G,C_Gln,mp%O_rate)  ! Note: this is the previous O_rate 
call set_param_set(r_G,C_L)
w = mp%f_G/f_Gu
r_P = (a + b*w)/(c + d*w)
if (r_P < 0) then
    clean = -1
    r_P = 0
    q = r_G*f_Gu*N_GI/(1 - f_N*N_GlnI)
    w = min(L1/(q*f_NG*(1 - f_Glnu)*N_GlnO),L2/(q*f_NG))
!    write(nflog,'(a,f6.3)') 'In L1: set r_P=0 ==================> w: ',w
endif
C_P = (r_P - ps%a0 - ps%b0*w)/ps%c0
if (C_P < 0) then
    clean = -1
    w = zero_C_P(a,b,c,d)
    C_P = 0
    r_P = r_G*(1 - w*f_Gu)*N_GP + V*K_LP*C_L
!    write(nflog,'(a,f8.3,e12.3)') 'In L1, C_P < 0, set w: ',w,r_P
endif
r_Gln = w*(r_P*ps%a2 + ps%b2)
if (r_Gln < L2) then
    done = .true.
    r_O2 = r_P*(ps%a4 + ps%b4*w) + ps%c4*w
    r_I = r_P*ps%a1*w + ps%b1*w
    r_A = r_P*(ps%a3 + ps%b3*w) + ps%c3 + ps%d3*w
    r_L = r_G*(1 - w*f_Gu)*N_GP - r_P
else
    aa = b*ps%a2 + d*ps%b2
    bb = a*ps%a2 + c*ps%b2 - d*L2
    cc = -c*L2
    dd = sqrt(bb*bb - 4*aa*cc)
    w1 = (-bb + dd)/(2*aa)
    w2 = (-bb - dd)/(2*aa)
!    write(nflog,'(a,2f8.3)') 'w1,w2: ',w1,w2
    w = min(1.0,w2)
!    w = max(0.0,w)
    if (w < 0) write(nflog,*) 'In L2: w < 0'
    ! Check eqtn (5): is r_O2 < L1?
    r_P = (a + b*w)/(c + d*w)
    if (r_P < 0) then
        clean = -1
        r_P = 0      
        q = r_G*f_Gu*N_GI/(1 - f_N*N_GlnI)
        w = min(L1/(q*f_NG*(1 - f_Glnu)*N_GlnO),L2/(q*f_NG))
        r_I = q*w
        r_Gln = f_NG*r_I
 !       write(nflog,'(a,3e12.3)') 'r_I parts: ',r_G*(1 - w*f_Gu)*N_GI, r_Gln*(1 - f_Glnu)*N_GlnI 
    else
        r_Gln = w*(r_P*ps%a2 + ps%b2)
        r_I = r_P*ps%a1*w + ps%b1*w
    endif
    C_P = (r_P - ps%a0 - ps%b0*w)/ps%c0
    if (C_P < 0) then
        write(nflog,*) 'In L2, C_P < 0, set w = 1'
        clean = -1
        w = 1
        C_P = 0
        r_P = r_G*(1 - w*f_Gu)*N_GP + V*K_LP*C_L
        r_Gln = w*(r_P*ps%a2 + ps%b2)
        r_I = r_P*ps%a1*w + ps%b1*w
    endif
    r_O2 = r_P*(ps%a4 + ps%b4*w) + ps%c4*w 
    r_A = r_P*(ps%a3 + ps%b3*w) + ps%c3 + ps%d3*w
    r_L = r_G*(1 - w*f_Gu)*N_GP - r_P
endif   
if (done) then
    mp%P_rate = r_P
    mp%G_rate = r_G
    mp%Gln_rate = r_Gln
    mp%I_rate = r_I
    mp%A_rate = r_A
    mp%O_rate = r_O2
    mp%L_rate = r_L
    mp%C_P = C_P
    mp%f_G = w*f_Gu
    mp%f_P = w*f_Pu
    mp%recalcable = clean
else
    write(nflog,*) 'What now???????????????????????????????????????????????????????????????'
endif
end subroutine

!==================================================================================================
!==================================================================================================
function zero_C_P(a,b,c,d) result(w)
real(REAL_KIND) :: a, b, c, d, w
real(REAL_KIND) :: aa, bb, cc, dd, w1, w2

aa = ps%b0*d
bb = ps%b0*c + ps%a0*d - b
cc = ps%a0*c - a
dd = sqrt(bb*bb - 4*aa*cc)
w1 = min(1.0,(-bb + dd)/(2*aa))
w2 = min(1.0,(-bb - dd)/(2*aa))
!write(nflog,'(a,2f8.3)') 'zero_C_P: w1,w2: ',w1,w2
w = w1
if (w < 0) then
    w = 0
endif
end function

#endif
end module


