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
real(REAL_KIND) :: Hill_Km_ON	! Hill Km for dependence of ON metabolism rate on ON 
real(REAL_KIND) :: Hill_N_ON	! Hill N for dependence of ON metabolism rate on ON 
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
real(REAL_KIND) :: C_ON_norm
real(REAL_KIND) :: O2_baserate
real(REAL_KIND) :: G_baserate
real(REAL_KIND) :: Gln_baserate
real(REAL_KIND) :: K_PL			! P -> L
real(REAL_KIND) :: K_LP			! L -> P
real(REAL_KIND) :: r_Pu, r_Gu, r_Glnu, r_Au, r_Iu, r_Ou, r_Lu, r_Onu, C_Pu   ! unconstrained conditions
real(REAL_KIND) :: G_maxrate, O2_maxrate, Gln_maxrate, ON_maxrate
real(REAL_KIND) :: f_N, f_NG, r_Abase, r_Ibase, C_Gln_cut, Km_rGln_factor
integer :: fgp_solver

!logical :: first_metab

type param_set_type
real(REAL_KIND) :: a0, b0, c0, d0, a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3, p, q, h
end type
!type param_set_type
!real(REAL_KIND) :: a0, b0, c0, a1, b1, a2, b2, a3, b3, c3, d3, a4, b4, c4, h
!end type

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
logical :: use_wxfGlnu = .false.
logical :: use_nitrogen = .true.
logical :: use_ON = .true.
integer :: knt
logical :: hyper_simple = .true.
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
!write(nflog,'(a,6e12.3)') 'get_glycosis_rate: H,C_G,(1 + K_Hb*H),metab,G_maxrate,rate: ',H,C_G,(1 + K_Hb*H),metab,G_maxrate,rate
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
integer :: ityp, it, N_P, N_O2, res
real(REAL_KIND) :: MM_O2, MM_P, V, K1, K2, Km_P, C_P, C_O2, C_G, C_L, C_Gln, C_ON, Cin(5) !, C_Gu, r_Gu, r_Pu, r_Lu, r_An, r_In
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
f_NG = f_N/(f_Glnu*N_GlnI)   ! this assumes that f_N is defined correctly: r_GlnI = f_N*r_I     ! not used now
!mp => phase_metabolic(1)

Hill_Km_O2 = chemo(OXYGEN)%MM_C0
Hill_N_O2 = chemo(OXYGEN)%Hill_N
Hill_Km_G = chemo(GLUCOSE)%MM_C0
Hill_N_G = chemo(GLUCOSE)%Hill_N
Hill_Km_Gln = chemo(GLUTAMINE)%MM_C0
Hill_N_Gln = chemo(GLUTAMINE)%Hill_N
Hill_Km_ON = chemo(OTHERNUTRIENT)%MM_C0
Hill_N_ON = chemo(OTHERNUTRIENT)%Hill_N
Hill_N_P = 1
Hill_Km_O2 = Km_O2_factor*Hill_Km_O2
N_O2 = chemo(OXYGEN)%Hill_N
O2_maxrate = chemo(OXYGEN)%max_cell_rate
write(nflog,*) 'O2_maxrate: ',O2_maxrate
G_maxrate = chemo(GLUCOSE)%max_cell_rate
Gln_maxrate = chemo(GLUTAMINE)%max_cell_rate
ON_maxrate = chemo(OTHERNUTRIENT)%max_cell_rate

O2_baserate = 0		! are these needed at all??
G_baserate = 0
H = get_HIF1steadystate(C_O2_norm)
call analyticSetPDK1(H, fPDK, 1.0d10)
mp%HIF1 = H
mp%PDK1 = fPDK

if (hyper_simple) then
    call get_unconstrained_rates_simple(res)
    mp%f_G = f_Gu
    mp%f_P = f_Pu
    mp%f_Gln = f_Glnu
    mp%G_rate = r_Gu
    mp%P_rate = r_Pu
    mp%A_rate = r_Au
    mp%I_rate = r_Iu
    mp%Gln_rate = r_Glnu
    mp%O_rate = r_Ou
    mp%ON_rate = r_ONu
!    mp%tagged = .false.
elseif (use_ON) then
    r_Glnu = Gln_maxrate
!    call get_unconstrained_rates_ON(res)
!    call set_param_set_ON(r_Gu,C_L_norm)
    mp%f_G = f_Gu
    mp%f_P = f_Pu
    mp%f_Gln = f_Glnu
    if (res /= 0) then
        write(nflog,*) 'get_unconstrained_rates_ON returned: ',res
        ok = .false.
        return
    endif
elseif (use_nitrogen) then
    r_Glnu = f_N*Gln_maxrate
!    call get_unconstrained_rates2
!    write(nflog,'(a,6e12.3)') 'got unconstrained rates: G,P,I,Gln,A,O: ',r_Gu,r_Pu,r_Iu,r_Glnu,r_Au,r_Ou
!    call set_param_set2(r_Gu,C_L_norm)
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
    C_ON = C_ON_norm
    MM_O2 = f_MM(C_O2,Hill_Km_O2,N_O2)
    r_Ou = O2_maxrate	! initial guess, for f_metab
    mp%f_G = f_Gu
    mp%f_P = f_Pu
    mp%f_Gln = f_Glnu
    mp%O_rate = r_Ou
    Cin = [C_O2, C_G, C_L, C_Gln, C_ON]
!    call f_metab(mp, Cin, C_Gln_norm, res)
    r_Gu = mp%G_rate
    r_Pu = mp%P_rate
    r_Au = mp%A_rate
    r_Iu = mp%I_rate
    r_Glnu = mp%Gln_rate
    r_Ou = mp%O_rate
    r_ONu = mp%ON_rate
    C_Pu = mp%C_P
    write(nflog,'(a,4e12.3)') 'f_Gu,f_Pu,r_Ag,C_Pu: ',f_Gu,f_Pu,r_Ag,C_Pu
endif
r_Ag = f_ATPg*r_Au
r_As = f_ATPs*r_Au
write(nflog,'(a,e12.3)') 'r_Ag: ',r_Ag
rIA = r_Iu/r_Au
write(nflog,'(a,3e12.3)') 'Unconstrained: r_G, r_Gln, r_ON: ',r_Gu,r_Glnu,r_ONu
write(nflog,'(a,2e12.3,f8.4)') 'r_I, r_A, r_I/r_A: ',r_Iu,r_Au,rIA
write(nflog,'(a,f8.4)') 'f_Gln with only Gln to preserve rIA: ',N_GlnA*rIA/(N_GlnI + N_GlnA*rIA)
mp%recalcable = -1
knt = 0
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
subroutine get_metab_rates(cp, Cin, C_GlnEx, res)
integer :: res
type(cell_type), pointer :: cp
real(REAL_KIND) :: Cin(:), C_GlnEx

!if (noSS) then
!    call f_metab_noSS(mp,Cin(OXYGEN),Cin(GLUCOSE),Cin(LACTATE),Cin(4))
!else
    call f_metab(cp,Cin(:), C_GlnEx, res)
!endif
end subroutine

!--------------------------------------------------------------------------
! Only for vmonolayer
!--------------------------------------------------------------------------
subroutine get_unconstrained_rates_simple(res)
integer :: res
real(REAL_KIND) :: r_GP, r_GIu, r_PIu, r_GPIu, r_GlnIu, r_ONIu, f_PP
type(metabolism_type), target :: metab
type(metabolism_type), pointer :: mp

write(nflog,*) 'get_unconstrained_rates_simple'
f_PP = f_PPu    ! was 5./85.
mp => metab
mp%HIF1 = get_HIF1steadystate(C_O2_norm)
call analyticSetPDK1(mp%HIF1, mp%PDK1, 1.0d10)
r_Ou = O2_maxrate
r_Gu = get_glycosis_rate(mp%HIF1,C_G_norm,C_Gln_norm,O2_maxrate)
r_Glnu = Gln_maxrate
r_ONu = f_rON_base*ON_maxrate
r_GP = (1 - f_Gu)*r_Gu*N_GP
r_Pu = f_PP*r_GP
r_Lu = (1 - f_PP)*r_GP
r_GIu = r_Gu*f_Gu*N_GI
r_PIu = r_Pu*f_Pu*N_PI
r_GPIu = r_GIu + r_PIu
r_GlnIu = r_Glnu*f_Glnu*N_GlnI
r_ONIu = r_ONu*f_ONu*N_ONI
r_Iu = r_GPIu + r_GlnIu + r_ONIu
r_Au = r_Gu*(1 - f_Gu)*N_GA + r_Pu*(1 - f_Pu)*N_PA + r_Glnu*(1 - f_Glnu)*N_GlnA + r_ONu*(1 - f_ONu)*N_ONA
write(nflog,'(a,4e12.3)') 'r_Gu, r_Pu, r_GIu, r_PIu: ',r_Gu, r_Pu, r_GIu, r_PIu
write(nflog,'(a,4e12.3)') 'r_GPIu, r_GlnIu, r_ONIu, r_Iu: ',r_GPIu, r_GlnIu, r_ONIu, r_Iu
write(nflog,'(a)') '---------------------------------------------------------------'
write(nflog,'(a,3f6.3)') 'fractions of ATP from: G, P, Gln: ', &
            r_Gu*(1 - f_Gu)*N_GA/r_Au, r_Pu*(1 - f_Pu)*N_PA/r_Au, r_Glnu*(1 - f_Glnu)*N_GlnA/r_Au
write(nflog,'(a)') '---------------------------------------------------------------'
res = 0
end subroutine

!--------------------------------------------------------------------------
! Determine q to make r_NI a fraction f of r_GlnI in unconstrained situation
! r_NIu = q*(r_GIu + r_PIu + r_GlnIu) = f*r_GlnIu
! q = f*r_GlnIu/(r_GPIu + r_GlnIu)
!f, q:  0.100 0.046
!f, q:  0.300 0.139
!f, q:  0.400 0.185
!f, q:  0.500 0.231
!f, q:  0.600 0.277
!f, q:  0.700 0.323
!f, q:  0.800 0.369
! NOT USED
!--------------------------------------------------------------------------
subroutine getq(f, q)
real(REAL_KIND) :: f, q
real(REAL_KIND) :: r_GlnIu, r_GPIu, f_PP

f_PP = f_PPu    ! was 5./85.
r_Glnu = Gln_maxrate
r_GPIu = (f_Gu*N_GI + f_Pu*N_PI*f_PP*(1 - f_Gu)*N_GP)*r_Gu
r_GlnIu = f_Glnu*N_GlnI*r_Glnu
q = f*r_GlnIu/(r_GPIu + r_GlnIu)
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function get_f_Gln(C_Gln) result(f)
real(REAL_KIND) :: C_Gln
real(REAL_KIND) :: f
real(REAL_KIND) :: Km_Gln, fcorrect
integer :: N_Gln

!C_Gln_lo = 0.25
!C_Gln_hi = 0.3
!f_rGln_lo = 0.2
if (C_Gln > C_Gln_hi) then
    f = 1.0
elseif (C_Gln > C_Gln_lo) then
    f = f_rGln_lo + (1.0 - f_rGln_lo)*(C_Gln - C_Gln_lo)/(C_Gln_hi - C_Gln_lo)
else
    N_Gln = chemo(GLUTAMINE)%Hill_N
    Km_Gln = chemo(GLUTAMINE)%MM_C0     ! Michaelis-Menten Km
    fcorrect =  f_rGln_lo/f_MM(C_Gln_lo,Km_Gln,N_Gln)
    f = fcorrect*f_MM(C_Gln,Km_Gln,N_Gln)
endif
end function


!--------------------------------------------------------------------------
subroutine f_metab(cp, Cin, C_GlnEx, res)
integer :: res
real(REAL_KIND) :: Cin(:), C_GlnEx
type(cell_type), pointer :: cp
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C_O2, C_G, C_L, C_Gln, C_ON
real(REAL_KIND) :: r_G, fPDK, w, f, f_PP, v, z, zmin, wlim
real(REAL_KIND) :: f_G, f_P, f_Gln, f_ON, r_P, r_A, r_I, r_L, r_Gln
real(REAL_KIND) :: r_GP, r_GA, r_PA, r_GlnA, Km_O2, MM_O2, Km_ON
real(REAL_KIND) :: r_GI, r_PI, r_GlnI, r_NI, r_GPI, r_GlnONI, r_ONI, r_ONIu, r_ON, r_ONA, r_GPA, r_O2, r_N, r_Nu
real(REAL_KIND) :: a, b, cc, d, e, dw, r_Atest, r_Atestq, w1, w2
integer :: N_O2, N_Gln, N_ON, k, Nw, iw
real(REAL_KIND) :: C, C0, C_Gln_min, f_Gln_C0, r_Gln_max, r_GlnI_max, r_GlnIu, f_Gln_max, r_ON_max, r_ONI_max
real(REAL_KIND) :: C_GlnEx_prev
logical :: use_ON = .true.

!cp => cell_list(kcell_now)     ! BAD can't use a global variable that changes and use OMP 
mp => cp%metab
if (cp%ATP_tag .or. cp%GLN_tag) then    ! the cell has been tagged to die
    res = 0
    return
endif

f_ON = f_ONu
f_PP = f_PPu    ! was 5./85.
f_Gln = f_Glnu
!C_Gln_min = C_Gln_cut    ! 0.02  ! growth suppressed below this extra-cellular conc  NOT USED
C0 = chemo(GLUTAMINE)%MM_C0
f_Gln_max = Km_rGln_factor  !2.0
!write(nflog,'(a,2f8.4)') 'DEBUG: C_GlnEx, C_Gln_min: ',C_GlnEx,C_Gln_min
fPDK = mp%PDK1
!r_Gln_max = fPDK*f_Gln_max*r_Glnu
r_Gln_max = fPDK*r_Glnu
r_GlnI_max = r_Gln_max*f_Gln*N_GlnI
r_GlnIu = r_Glnu*f_Gln*N_GlnI   ! no fPDK!
res = 0
C_O2 = max(0.0,Cin(OXYGEN))
C_G = max(0.0,Cin(GLUCOSE))
C_L = max(0.0,Cin(LACTATE))
C_Gln = max(0.0,Cin(GLUTAMINE))
C_ON = max(0.0,Cin(OTHERNUTRIENT))

N_O2 = Hill_N_O2
Km_O2 = Hill_Km_O2
r_O2 = mp%O_rate
r_G = get_glycosis_rate(mp%HIF1,C_G,C_Gln,r_O2)	! Note: r_O2 is the previous O_rate - not used
                                                ! dependence on C_Gln not wanted now - not used
v = min(1.0,r_G/r_Gu)
MM_O2 = f_MM(C_O2,Km_O2,N_O2)

! This is probably valid only for vmonolayer, not when different cells see different C_GlnEx

if (mp%C_GlnEx_prev == 0) then
    C = C_GlnEx
else
    C = (C_GlnEx + mp%C_GlnEx_prev)/2
endif
mp%C_GlnEx_prev = C_GlnEx

f_G = f_Gu
f_P = f_Pu

w = get_f_Gln(C)    ! this is the fraction of r_Glnu 
f_Gln = f_Glnu
r_Gln = w*r_Glnu

r_GlnI = r_Gln*f_Gln*N_GlnI

r_GI = f_G*r_G*N_GI
r_GA = (1 - f_G)*r_G*N_GA
r_GP = (1 - f_G)*r_G*N_GP
r_P = f_PP*r_GP
r_L = (1 - f_PP)*r_GP
r_PI = f_P*r_P*N_PI
r_PA = (1 - f_P)*r_P*N_PA
r_GPI = r_GI + r_PI

Km_ON = chemo(OTHERNUTRIENT)%MM_C0
N_ON = chemo(OTHERNUTRIENT)%Hill_N
r_ON_max = ON_maxrate*f_MM(C_ON,Km_ON,N_ON) 

r_ONI_max = r_ON_max*f_ON*N_ONI
r_ONI = min(r_Iu - r_GPI - r_GlnI, r_ONI_max) 
r_ONI = max(r_ONI,0.0)
! Try turning this off
!if (w < f_rGln_lo) then
!    r_ONI = r_ONI*w/(f_rGln_lo)   ! to reduce r_ON when r_Gln goes low
!endif

r_ON = r_ONI/(f_ON*N_ONI)

r_I = r_GPI + r_GlnI + r_ONI
! Making ON also a Nitrogen contributor
r_N = f_IN*(Gln_Nshare*r_GlnI + (1-Gln_Nshare)*r_ONI)
!r_Nu = (GLN_Nshare*r_Glnu + (1-Gln_Nshare)*r_ONu)
!write(nflog,'(a,4e12.3)') 'r_Gln, r_GlnI, f_IN, r_N: ',r_Gln, r_GlnI, f_IN, r_N
!write(nflog,'(a,3e12.3)') 'r_N, f_rGln_threshold*r_Iu: ',r_N,f_rGln_threshold*r_Iu

!if (r_N < f_rGln_threshold*r_Nu) then    ! death
if (r_N < f_rGln_threshold*r_Iu) then    ! death
!    write(nflog,'(a,4e12.3)') 'tagged for death from low r_Gln: ',r_N,f_rGln_threshold,r_Iu,f_rGln_threshold*r_Iu
    cp%GLN_tag = .true.
    mp%f_G = f_G
    mp%f_P = f_P
    mp%f_Gln = f_Gln
    mp%G_rate = 0
    mp%A_rate = 0
    mp%I_rate = 0
    mp%P_rate = 0
    mp%O_rate = 0
    mp%Gln_rate = 0
    mp%ON_rate = 0
    mp%L_rate = 0
    res = 0
!    write(nflog,'(a,i6,5f8.3,2e12.3)') 'GLN_tag: ',cp%ID,C_GlnEx,C,w,cp%Cex(GLUTAMINE),cp%Cin(GLUTAMINE),r_GlnI,r_N
    return
endif
!write(nflog,'(a,f6.3,5e12.3)') 'w,C_GlnEx,C,r_Gln,r_ON: ',w,C_GlnEx,C,r_Gln,r_ON
!write(nflog,'(a,4e12.3)') 'C_ON, Km_ON, r_ON, r_ON_max: ',C_ON, Km_ON, r_ON, r_ON_max
r_ONA = (1 - f_ON)*r_ON*N_ONA

!r_Gln = r_GlnI/(f_Gln*N_GlnI)
r_GlnA = (1 - f_Gln)*r_Gln*N_GlnA

r_A = r_GA + r_PA + r_GlnA + r_ONA
if (r_A < r_Ag) then    ! solve for w s.t. with w*f_G, w*f_P, w*f_Gln, r_A = r_Ag
!   r_A = r_GlnA + ((1-w*f_G)*N_GA + (1-w*f_P)*N_PA*f_PP*(1-w*f_G)*N_GP)*r_G
!   now have added in (1 - w*f_ON)*N_ONA*r_ON
! => quadratic in w
!    write(nflog,'(a,3e12.3)') 'r_GA, r_PA, r_GlnA: ',r_GA, r_PA, r_GlnA
!    write(nflog,'(a,f8.3,2e12.3)') 'w, r_A, r_Ag: ',w, r_A, r_Ag
!    if (w > 0) then
        e = N_PA*f_PP*N_GP*r_G
        a = e*f_P*f_G
!        b = -(f_Gln*N_GlnA*r_Gln + f_ON*N_ONA*r_ON + f_G*N_GA*r_G + e*(f_G+f_P))
!        cc = N_GlnA*r_Gln + N_ONA*r_ON + N_GA*r_G + e - r_Ag
        b = -(f_ON*N_ONA*r_ON + f_G*N_GA*r_G + e*(f_G+f_P))
        cc = r_GlnA + N_ONA*r_ON + N_GA*r_G + e - r_Ag
        d = sqrt(b*b - 4*a*cc) 
        w1 = (-b + d)/(2*a)
        w2 = (-b - d)/(2*a)
!        write(nflog,'(a,6e11.3)') 'a,b,cc,d,w1,w2: ',a,b,cc,d,w1,w2
        if (w2 < 0) then
            w = 0
        elseif (w2 > 1) then
            write(nflog,*) 'ERROR: w to adjust r_A > 1: ',w2
!            res = 1
!            return
            w = 1
        else
            w = w2
        endif
!    endif
    
    write(nflog,'(a,f8.3)') 'Adjusting r_A: w: ',w
    f_G = w*f_G
    f_P = w*f_P
!    f_Gln = w*f_Gln
    f_ON = w*f_ON   ! added
    
    r_GI = f_G*r_G*N_GI
    r_GA = (1 - f_G)*r_G*N_GA
    r_GP = (1 - f_G)*r_G*N_GP
    r_P = f_PP*r_GP
    r_L = (1 - f_PP)*r_GP
    r_PI = f_P*r_P*N_PI
    r_PA = (1 - f_P)*r_P*N_PA
!    r_GlnI = r_Gln*f_Gln*N_GlnI
!    r_GlnA = (1 - f_Gln)*r_Gln*N_GlnA
    r_ONI = r_ON*f_ON*N_ONI         ! added
    r_ONA = (1 - f_ON)*r_ON*N_ONA   ! added
    r_GPI = r_GI + r_PI
endif

r_A = r_GA + r_PA + r_GlnA + r_ONA
if (r_A < r_As) then
    write(nflog,*) 'death from r_A'
    cp%ATP_tag = .true.
    mp%f_G = f_G
    mp%f_P = f_P
    mp%f_Gln = f_Gln
    mp%G_rate = 0
    mp%A_rate = 0
    mp%I_rate = 0
    mp%P_rate = 0
    mp%O_rate = 0
    mp%Gln_rate = 0
    mp%ON_rate = 0
    mp%L_rate = 0
    res = 0
    return
endif

r_I = r_GPI + r_GlnI + r_ONI
!write(nflog,'(a,f7.3)') 'Fraction of intermediates from glucose: ',r_GPI/r_I
r_O2 = (1 - f_P)*r_P*N_PO + (1 - f_Gln)*r_Gln*N_GlnO

mp%f_G = f_G
mp%f_P = f_P
mp%f_Gln = f_Gln
mp%G_rate = r_G
mp%A_rate = r_A									! production
mp%I_rate = r_I									! production
mp%P_rate = r_P									! utilisation
mp%O_rate = r_O2								! consumption
mp%Gln_rate = r_Gln								! consumption
mp%ON_rate = r_ON								! consumption
mp%L_rate = r_L									! production
!if (r_N < f_rGln_threshold*r_Iu) then    ! death
!    write(nflog,*) 'death'
!    mp%tagged = .true.
!endif
!write(nflog,'(a,5e12.3)') 'r_G, r_Gln, r_Glnu, r_I, r_A: ',r_G,r_Gln,r_Glnu,r_I,r_A 
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
! Test new simple model 22/09/2020
! Try adding in N-intermediates from glutamine, fraction q of total intermediates
! Need to make r_I -> 0 as r_A -> r_Ag
! Crude method: if r_A < r_Ag, adjust all f_G, f_P, f_Gln 
! to bring r_A up to r_Ag, reducing r_GI, r_PI, r_GlnI accordingly.
!--------------------------------------------------------------------------
subroutine f_metab8(mp, Cin, C_GlnEx, res)
integer :: res
real(REAL_KIND) :: Cin(:), C_GlnEx
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C_O2, C_G, C_L, C_Gln
real(REAL_KIND) :: r_G, fPDK, w, q, f, f_PP, v, z, zmin, wlim
real(REAL_KIND) :: f_G, f_P, f_Gln, f_ON, r_P, r_A, r_I, r_L, r_Gln
real(REAL_KIND) :: r_GP, r_GA, r_PA, r_GlnA, Km_O2, MM_O2
real(REAL_KIND) :: r_GI, r_PI, r_GlnI, r_NI, r_GPI, r_GlnONI, r_ONI, r_ONIu, r_ON, r_ONA, r_GPA, r_O2
real(REAL_KIND) :: a, b, cc, d, e, dw, r_Atest, r_Atestq, w1, w2
integer :: N_O2, N_Gln, k, Nw, iw
real(REAL_KIND) :: C, C0, C_Gln_min, f_Gln_C0, r_Gln_max, r_GlnI_max, r_GlnIu, f_Gln_max, r_ON_max, r_ONI_max
logical :: use_ON = .true.

f_ON = f_ONu
f_PP = f_PPu    ! was 5./85.
q = f_IN
f_Gln = f_Glnu
C_Gln_min = C_Gln_cut    ! 0.02  ! growth suppressed below this extra-cellular conc  NOT USED
C0 = chemo(GLUTAMINE)%MM_C0
f_Gln_max = Km_rGln_factor  !2.0
!write(nflog,'(a,2f8.4)') 'DEBUG: C_GlnEx, C_Gln_min: ',C_GlnEx,C_Gln_min
fPDK = mp%PDK1
!r_Gln_max = fPDK*f_Gln_max*r_Glnu
r_Gln_max = fPDK*r_Glnu
r_GlnI_max = r_Gln_max*f_Gln*N_GlnI
r_GlnIu = r_Glnu*f_Gln*N_GlnI   ! no fPDK!
!r_ON_max = r_ONu
r_ON_max = ON_maxrate       ! f_Gln_max*r_ONu      ! use same max rate factor as for Gln
r_ONI_max = r_ON_max*f_ON*N_ONI
write(nflog,'(a,4e12.3)') 'r_Glnu,r_ONu: ',r_Glnu,r_ONu
res = 0
C_O2 = max(0.0,Cin(OXYGEN))
C_G = max(0.0,Cin(GLUCOSE))
C_L = max(0.0,Cin(LACTATE))
C_Gln = max(0.0,Cin(GLUTAMINE))

N_O2 = Hill_N_O2
Km_O2 = Hill_Km_O2
r_O2 = mp%O_rate
r_G = get_glycosis_rate(mp%HIF1,C_G,C_Gln,r_O2)	! Note: r_O2 is the previous O_rate - not used
                                                        ! dependence on C_Gln not wanted now - not used
v = min(1.0,r_G/r_Gu)
MM_O2 = f_MM(C_O2,Km_O2,N_O2)

if (mp%C_GlnEx_prev == 0) then
    C = C_GlnEx
else
    C = (C_GlnEx + mp%C_GlnEx_prev)/2
endif
mp%C_GlnEx_prev = C_GlnEx

N_Gln = chemo(GLUTAMINE)%Hill_N
if (C < 0) then
    w = 0
else
    w = C**N_Gln/(C0**N_Gln + C**N_Gln)
endif

! To test slowly reducing r_ON_max
!if (w < 0.01) then
!    r_ON_max = r_ON_max*(1 - (1-w)/1000)
!    r_ONI_max = r_ON_max*f_ON*N_ONI
!endif
!write(nflog,'(a,f8.3,4e12.3)') 'w,r_ON_max,r_ONI_max: ',w,r_ON_max,r_ONI_max,f_ON,N_ONI

!r_G = ((1+w)/2)*r_G     !!!!!!!!!!!!!!!!!!!!!!!!!! test !!!!!!!!!!!!!!!!!!!!!!!!!!!!

f_G = w*f_Gu
f_P = w*f_Pu
f_Gln = f_Glnu

r_GI = f_G*r_G*N_GI
r_GA = (1 - f_G)*r_G*N_GA
r_GP = (1 - f_G)*r_G*N_GP
f_PP = w*f_PP      !!!!!!!!!!!!!!!!!!!!!!!!!! test !!!!!!!!!!!!!!!!!!!!!!!!
r_P = f_PP*r_GP
r_L = (1 - f_PP)*r_GP
r_PI = f_P*r_P*N_PI
r_PA = (1 - f_P)*r_P*N_PA
r_GPI = r_GI + r_PI

! Note: 
!   r_GI = f_G*N_GI*r_G
!   r_PI = f_P*N_PI*r_P = f_P*N_PI*f_PP*r_GP = f_P*N_PI*f_PP*(1 - f_G)*N_GP*r_G
! therefore:
!   r_GPI = (f_G*N_GI + f_P*N_PI*f_PP*(1 - f_G)*N_GP)*r_G
!
!   r_GA = (1-f_G)*N_GA*r_G
!   r_PA = (1-f_P)*N_PA*r_P = (1-f_P)*N_PA*f_PP*r_GP = (1-f_P)*N_PA*f_PP*(1-f_G)*N_GP*r_G
! therefore:
!   r_GPA = ((1-f_G)*N_GA + (1-f_P)*N_PA*f_PP*(1-f_G)*N_GP)*r_G

!r_Gln = Min(r_Gln_max, (r_Iu - r_GI - r_PI)/(f_Gln*N_GlnI)) 

if (use_ON) then
    r_GlnI = min(r_Iu - r_GPI, w*r_GlnIu)
    r_Gln = r_GlnI/(f_Gln*N_GlnI)
    write(nflog,'(a,2f8.4,e12.3)') 'w, C_GlnEx, r_GlnI: ',w,C_GlnEx,r_GlnI
    z = 1
    wlim = 0.02
    zmin = 0.3
    if (w < wlim) then      ! explain
        z = zmin + (1.0 - zmin)*w/wlim
    endif
    r_ONI = min(r_Iu - r_GPI - r_GlnI, r_ONI_max) 
    r_ONI = z*r_ONI
!    r_ONI = ((w+1)/2)*r_ONI     !!!!!!!!!!!!!!!!!!!!!!!! testing !!!!!!!!!!!!!!!!!!!!!!
    r_ON = r_ONI/(f_ON*N_ONI)
    write(nflog,'(a,2e12.3)') 'r_Gln, r_ON: ',r_Gln, r_ON
    r_ONA = (1 - f_ON)*r_ON*N_ONA
    r_I = r_GPI + r_GlnI + r_ONI
    write(nflog,'(a,2f6.3,5e12.3)') 'w,z,r_GPI,r_GlnI,r_ONI,r_I,r_Iu: ',w,z,r_GPI,r_GlnI,r_ONI,r_I,r_Iu
else
    r_GlnI = Min(r_GlnI_max, r_GlnONI)
    r_ONI = 0
    r_ON = 0
    r_ONA = 0
endif

!write(nflog,'(a,4e12.3,f8.3)') 'r_Iu, r_GPI, r_GlnI, r_I, w: ',r_Iu, r_GPI, r_GlnI, r_GPI + r_GlnI, w
!if (r_GlnI >= r_GlnI_max) write(nflog,'(a,4e12.3)') 'DEBUG: r_G,r_GPI,r_GlnI,r_GlnI_max: ',r_G,r_GPI,r_GlnI,r_GlnI_max
!write(nflog,'(a,5e12.3)') 'r_Gln_max,r_Iu,r_GI,r_PI,f_Gln*N_GlnI: ',r_Gln_max,r_Iu,r_GI,r_PI,f_Gln*N_GlnI
!write(nflog,'(a,e12.3)') 'r_Gln: ',r_Gln
!r_Gln = w*fPDK*r_Gln'
!r_GlnI = f_Gln*r_Gln*N_GlnI

! Now account for N-intermediates r_NI fraction q of total r_I
! NOT USED NOW
if (.not.use_ON .and. r_GlnI < r_GPI*q/(1-q)) then
    if (r_GlnI_max > r_GPI*q/(1-q)) then
        r_GlnI = r_GPI*q/(1-q)
        write(nflog,'(a,e12.3)') 'DEBUG: r_GlnI: ',r_GlnI
    else    ! reduce r_G
        r_GlnI = r_GlnI_max
        r_GPI = r_GlnI*(1-q)/q
        r_G = r_GPI/(f_G*N_GI + f_P*N_PI*f_PP*(1 - f_G)*N_GP)
        r_GI = f_G*r_G*N_GI
        r_GA = (1 - f_G)*r_G*N_GA
        r_GP = (1 - f_G)*r_G*N_GP
        r_P = f_PP*r_GP
        r_L = (1 - f_PP)*r_GP
        r_PI = f_P*r_P*N_PI
        r_PA = (1 - f_P)*r_P*N_PA
        write(nflog,'(a,3e12.3)') 'DEBUG: reduce r_G: r_GlnI,r_G,r_GPI: ',r_GlnI,r_G,r_GPI
    endif
endif

r_Gln = r_GlnI/(f_Gln*N_GlnI)
r_GlnA = (1 - f_Gln)*r_Gln*N_GlnA

r_A = r_GA + r_PA + r_GlnA + r_ONA
if (r_A < r_Ag) then    ! solve for w s.t. with w*f_G, w*f_P, w*f_Gln, r_A = r_Ag
!   r_A = (1 - w*f_Gln)*N_GlnA*r_Gln + ((1-w*f_G)*N_GA + (1-w*f_P)*N_PA*f_PP*(1-w*f_G)*N_GP)*r_G
!   now have added in (1 - w*f_ON)*N_ONA*r_ON
! => quadratic in w
    write(nflog,'(a,3e12.3)') 'r_GA, r_PA, r_GlnA: ',r_GA, r_PA, r_GlnA
    write(nflog,'(a,f8.3,2e12.3)') 'w, r_A, r_Ag: ',w, r_A, r_Ag
    if (w > 0) then
        e = N_PA*f_PP*N_GP*r_G
        a = e*f_P*f_G
        b = -(f_Gln*N_GlnA*r_Gln + f_ON*N_ONA*r_ON + f_G*N_GA*r_G + e*(f_G+f_P))
        cc = N_GlnA*r_Gln + N_ONA*r_ON + N_GA*r_G + e - r_Ag
        if (.true.) then
        d = sqrt(b*b - 4*a*cc) 
        w1 = (-b + d)/(2*a)
        w2 = (-b - d)/(2*a)
        write(nflog,'(a,6e11.3)') 'a,b,cc,d,w1,w2: ',a,b,cc,d,w1,w2
        if (w2 < 0) then
            w = 0
        elseif (w2 > 1) then
            write(nflog,*) 'ERROR: w to adjust r_A > 1: ',w2
!            res = 1
!            return
            w = 1
        else
            w = w2
        endif
    endif
    endif
    
!    Nw = 10 
!    dw = 1.0/Nw
!    do iw = Nw,0,-1
!        w = iw*dw
!        r_Atest = (1 - w*f_Gln)*N_GlnA*r_Gln + ((1-w*f_G)*N_GA + (1-w*f_P)*N_PA*f_PP*(1-w*f_G)*N_GP)*r_G
!        r_Atestq = a*w*w + b*w + cc + r_Ag
!        write(nflog,'(a,i4,f8.3,3e12.3)') 'iw,w,r_Atestq,r_Atest,r_Ag: ',iw,w,r_Atestq,r_Atest,r_Ag
!        if (r_Atest > r_Ag) exit
!    enddo

    write(nflog,'(a,f8.3)') 'Adjusting r_A: w: ',w
    f_G = w*f_G
    f_P = w*f_P
    f_Gln = w*f_Gln
    f_ON = w*f_ON   ! added
    
    r_GI = f_G*r_G*N_GI
    r_GA = (1 - f_G)*r_G*N_GA
    r_GP = (1 - f_G)*r_G*N_GP
    r_P = f_PP*r_GP
    r_L = (1 - f_PP)*r_GP
    r_PI = f_P*r_P*N_PI
    r_PA = (1 - f_P)*r_P*N_PA
    r_GlnI = r_Gln*f_Gln*N_GlnI
    r_GlnA = (1 - f_Gln)*r_Gln*N_GlnA
    r_ONI = r_ON*f_ON*N_ONI         ! added
    r_ONA = (1 - f_ON)*r_ON*N_ONA   ! added
    r_GPI = r_GI + r_PI
endif

r_A = r_GA + r_PA + r_GlnA + r_ONA
!write(nflog,'(a,2e12.3)') 'Checking r_A, r_Ag: ',r_A,r_Ag
r_I = r_GPI + r_GlnI + r_ONI
r_NI = q*r_I
r_O2 = (1 - f_P)*r_P*N_PO + (1 - f_Gln)*r_Gln*N_GlnO

!if (r_I == 0) then
!    write(nflog,*) 'r_I = 0: STOPPING'
!    res = 1
!    return
!endif
!write(nflog,'(a,4e12.3)') 'r_GI, r_PI, r_GlnI, r_I: ',r_GI, r_PI, r_GlnI, r_I
!write(nflog,'(a)') '----------------------------------------------------------------------------'

mp%f_G = f_G
mp%f_P = f_P
mp%f_Gln = f_Gln
mp%G_rate = r_G
mp%A_rate = r_A									! production
mp%I_rate = r_I									! production
mp%P_rate = r_P									! utilisation
mp%O_rate = r_O2								! consumption
mp%Gln_rate = r_Gln								! consumption
mp%ON_rate = r_ON								! consumption
mp%L_rate = r_L									! production
!write(nflog,'(a,5e12.3)') 'r_G, r_Gln, r_Glnu, r_I, r_A: ',r_G,r_Gln,r_Glnu,r_I,r_A 
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

return
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
real(REAL_KIND) :: C_O2, C_G, C_L, C_Gln, C_ON, C_P, r_Atarget, Cin(5)
real(REAL_KIND) :: Km_P, fPDK, MM_O2, MM_Gln, r_G, r_Gln, r_P, r_A, r_I, r_O2, f_G, f_P, f_Gln
real(REAL_KIND) :: F, V, K1, K2, a, b, cc, a1, b1, c1, q1, q2, q3, q4, q5
real(REAL_KIND) :: rtol, x0, x1, dx, F0, F1, dFdx, w, z
real(REAL_KIND) :: average_volume = 1.2
integer :: N_O2, N_Gln, k, res

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
Cin = [C_O2,C_G,C_L,C_Gln,C_ON]
!call f_metab(mp,Cin,1.0d0,res)
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
! Not used
subroutine get_unconstrained_rates_ON(res)
integer :: res
type(metabolism_type), target :: metab
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C(NUTS), C_GlnEx, r_P, r_G, MM_Gln, MM_ON, f_Gln, r_Gln, r_ON, r_GlnI, r_ONI, r_I, r_Pc
real(REAL_KIND) :: w, h, hp, f1, f0, f_cutoff, MM_rGln, r_A, ratio
integer :: k
logical :: iter = .true.

mp => metab
mp%HIF1 = get_HIF1steadystate(C_O2_norm)
call analyticSetPDK1(mp%HIF1, mp%PDK1, 1.0d10)
mp%O_rate = O2_maxrate
C = [C_O2_norm, C_G_norm, C_L_norm, C_Gln_norm, C_ON_norm]
C_GlnEx = C_Gln_norm
write(nflog,'(a,5f8.3)') 'C: ', C
! Try to set Gln_maxrate from O2_maxrate, assuming that maximum r_P occurs when w = 1, and r_L = 0
!r_G = get_glycosis_rate(mp%HIF1,C_G_norm,C_Gln_norm,O2_maxrate)
!r_P = r_G*(1 - f_Gu)*N_GP
!Gln_maxrate = (O2_maxrate - r_P*N_PO)/N_GlnO
!write(nflog,'(a,e12.3)') 'Gln_maxrate computed from r_O: ',Gln_maxrate

r_G = get_glycosis_rate(mp%HIF1,C_G_norm,C_Gln_norm,O2_maxrate)
f_Gln = f_Glnu
if (hyper_simple) then
    r_Gu = r_G
    r_Glnu = Gln_maxrate
    r_ONu = ON_maxrate
    r_P = r_G*(1 - f_Gu)*N_GA
    r_GlnI = r_Glnu*f_Glnu*N_GlnI
    r_Iu = r_G*f_Gu*N_GI + r_P*f_Pu*N_PI + r_GlnI
    r_Au = r_G*(1 - f_Gu)*N_GA + r_P*(1 - f_Pu)*N_PA + r_Gln*(1 - f_Glnu)*N_GlnA
    res = 0
    return
endif
MM_Gln = f_MM(C_Gln_norm,Hill_Km_Gln,int(Hill_N_Gln))
MM_ON = f_MM(C_ON_norm,Hill_Km_ON,int(Hill_N_ON))
ratio = (C_Gln_norm/C_ON_norm)    ! ratio of r_Gln/r_ON
f_cutoff = 1
MM_rGln = 1
f0 = f_cutoff*MM_Gln
r_Gln = f0*Gln_maxrate
r_GlnI = r_Gln*f_Gln*N_GlnI
!r_ON = MM_rGln*MM_ON*ON_maxrate      ! Note: = 0 if ON_maxrate = 0.
r_ON = r_Gln/ratio
r_ONI = r_ON*N_ONI

if (iter) then
! Initial guess: use w = 1, r_L = 0, C_P = C_L_norm 
w = 0.5
r_P = r_G*(1 - w*f_Gu)*N_GA
r_GlnI = r_Gln*w*f_Gln*N_GlnI
r_Iu = r_G*w*f_Gu*N_GI + r_P*w*f_Pu*N_PI + r_GlnI + r_ONI
r_Au = r_G*(1 - w*f_Gu)*N_GA + r_P*(1 - w*f_Pu)*N_PA + r_Gln*(1 - w*f_Gln)*N_GlnA
r_Ag = f_ATPg*r_Au
h = (r_Au - r_Ag)/r_Iu
hp = h
write(nflog,'(a,5e12.3)') 'initial guess: r_G,r_P,r_Iu,r_Au,h: ',r_G,r_P,r_Iu,r_Au,h
do k = 1,10
!    call f_metab(mp, C, C_GlnEx, res)
    if (res /= 0) then
        write(nflog,*) 'f_metab failed'
        return
    endif
! r_Ou = r_Pu*N_PO + r_Glnu*N_GlnO  (assuming that ON utilisation does not consume O2) 
! Set Gln_maxrate = r_Glnu
!    Gln_maxrate = (mp%O_rate - mp%P_rate*N_PO)/N_GlnO
! oops! drives Gln_maxrate -> 0 
! Let N_GlnO = N_PO, solve for this
    N_PO = O2_maxrate/(mp%P_rate + mp%Gln_rate)
    N_GlnO = N_PO
    r_Au = mp%A_rate
    r_Iu = mp%I_rate
    r_Ag = f_ATPg*r_Au
    h = (r_Au - r_Ag)/r_Iu
    h = (h + hp)/2
    hp = h
    write(nflog,'(a,6e12.3)') 'Unconstrained: r_Au, r_Ag, h, N_PO, mp%O_rate,O2_maxrate: ',r_Au, r_Ag, h, N_PO, mp%O_rate, O2_maxrate
    write(nflog,*)
enddo
r_Gu = mp%G_rate
r_Glnu = mp%Gln_rate
r_ONu = mp%ON_rate
r_Lu = mp%L_rate
r_Pu = mp%P_rate
r_Iu = mp%I_rate
r_Au = mp%A_rate
r_Ou = mp%O_rate
else    ! Try direct solution for r_Au, r_Iu, h 
    f1 = 0.0
    w = 0.0
    r_P = r_G*((1 - w*f_Gu)*N_GP - f_GL)
    r_A = r_G*(1 - w*f_Gu)* N_GA + r_P*(1 - w*f_Pu)*N_PA + r_Gln*(1 - f1)*N_GlnA
    r_I = r_G*w*f_Gu*N_GI + r_P*w*f_Pu*N_PI + r_Gln*f1*N_GlnI + r_ON*N_ONI
    r_Au = r_A
    r_Iu = r_I
    r_Ag = f_ATPg*r_Au
    h = (r_Au - r_Ag)/r_Iu
    write(nflog,'(a,3e12.3)') 'r_Au,r_Iu,h: ',r_Au,r_Iu,h
    r_Gu = r_G
    r_Glnu = r_Gln
    r_ONu = r_ON
    r_Lu = f_GL*r_Gu
    r_Pu = r_P
    r_Iu = r_I
    r_Au = r_A
    r_Ou = 0
    res = 0
endif
write(nflog,'(a10,8a12)') 'rates: ','r_G', 'r_Gln', 'r_ON', 'r_L', 'r_P', 'r_I', 'r_A', 'r_O2'
write(nflog,'(a10,8e12.3)') 'unconstr: ',r_Gu,r_Glnu,r_Onu,r_Lu,r_Pu,r_Iu,r_Au,r_Ou
write(nflog,'(a,f8.3)') 'C_P: ',mp%C_P
end subroutine


end module


