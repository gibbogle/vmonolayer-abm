module optimise
use metabolism
implicit none

integer, parameter :: N_r_G = 20
integer, parameter :: N_C_L = 20
integer, parameter :: N_r_Pm = 20
real(REAL_KIND) :: r_G_array(N_r_G)
real(REAL_KIND) :: C_L_array(N_C_L)
real(REAL_KIND) :: r_Pm_array(N_r_Pm)

contains


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine optimiser(ityp, r_G, C_L, r_Pm_base, mp, x, y, C_P)
integer :: ityp
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: r_G, r_Pm_base, C_L, x, y, C_P
real(REAL_KIND) :: f_G, f_P, r_P, r_A, r_I, r_L, r_O2, f_PO, f_PA
real(REAL_KIND) :: f_G_prev, f_P_prev
real(REAL_KIND) :: K1, K2
real(REAL_KIND) :: r_Pt, r_GP, r_GA, r_PA, r_Pm, V, Km_O2, Km_P, a, b, c, d, e, MM_P, r_A_target
real(REAL_KIND) :: r_PA_max, MM_P_max, c1, c2, C_P_temp, C_O2, MM_O2, dA, f_G_new
integer :: N_O2, N_P, it
real(REAL_KIND) :: C_P_max, C_P_limit = 10

N_P = 1
Km_P = Hill_Km_P
N_O2 = chemo(OXYGEN)%Hill_N
Km_O2 = chemo(OXYGEN)%MM_C0
V = Vcell_cm3

f_PO = N_PO
f_PA = N_PA
K1 = K_PL
K2 = K_LP
f_G = f_G_norm
f_P = f_P_norm
r_A_target = r_A_norm

if (r_G < 0.01) then
	r_G = max(0.0d0,r_G)
	f_G = 0
	r_GA = 2*r_G
!	r_P = r_GA - V*(K1*C_P - K2*C_L) = MM(C_P)*r_Pm_base/(1-f_P) = r_Pm_base*C_P/((1-f_P)*(Km_P + C_P))
	e = r_GA + V*K2*C_L
	a = V*K1
	b = r_Pm_base/(1-f_P) - e + V*K1*Km_P
	c = -e*Km_P
	d = sqrt(b*b - 4*a*c)
	C_P = (-b + d)/(2*a)
	r_P = r_GA - V*(K1*C_P - K2*C_L)	
	r_PA = f_PA*(1-f_P)*r_P
else

it = 0
do
	it = it + 1
	write(*,'(a,i2,e12.3)') 'r_A_target: ',it,r_A_target
	c1 = r_Pm_base*(f_PA + 1/(1-f_P))
	c2 = r_A_target + V*K2*C_L
	a = V*K1
	b = V*K1*Km_P + c1 - c2
	c = -Km_P*c2
!	C_P = solve_C_P(a,b,c)	! solves a*x^2 + b*x + c = 0
	d = sqrt(b*b - 4*a*c)
	C_P = (-b + d)/(2*a)
	r_GA = r_A_target - f_PA*r_Pm_base*C_P/(Km_P + C_P)
	r_P = r_GA - V*(K1*C_P - K2*C_L)
	C_P_max = (r_GA + V*K2*C_L)/(V*K1)
	if (C_P > C_P_max) then
		write(*,*) 'C_P_max is exceeded!!!!'
	endif
	write(*,'(a,4f8.4)') 'f_P, f_G, C_P, C_P_max: ',f_P, f_G, C_P, C_P_max
	write(*,'(a,2e12.3)') 'r_P, r_GA: ',r_P,r_GA
	dA = r_GA - 2*(1-f_G)*r_G
	if (dA < 0) then
		r_PA = f_PA*(1-f_P)*r_P
		r_A = r_GA + r_PA
		write(*,'(a,3e12.3)') 'r_GA,r_PA,r_A: ',r_GA,r_PA,r_A/r_A_norm
		write(*,*)
		exit
	else
		write(*,'(a,e12.3)') 'dA: ',dA
		f_G_new = max(0.0d0,f_G - dA/(2*r_G))
		dA = dA - 2*(f_G - f_G_new)*r_G
		if (f_G_new == 0) dA = 1.01*dA
		write(*,'(a,2e12.3)') 'f_G_new, dA: ',f_G_new, dA
		write(*,*)
		f_G = f_G_new
		r_A_target = r_A_target - dA
		cycle
	endif
enddo
endif

mp%A_rate = r_GA + r_PA				! production
mp%I_rate = f_G*r_G + f_P*r_P		! production
mp%P_rate = r_P						! utilisation
mp%O_rate = f_PO*r_P*(1-f_P)		! consumption
mp%L_rate = V*(K1*C_P - K2*C_L)		! production

x = f_G
y = f_P
write(*,'(a,3f10.6)') 'f_G,f_P,C_P: ',f_G,f_P,C_P
write(*,'(a,6e10.3)') 'r_G,A,P,I,L,O2: ',mp%G_rate,mp%A_rate,mp%P_rate,mp%I_rate,mp%L_rate,mp%O_rate
write(*,'(a,6f8.3)') 'normalised P,A,I: ',mp%P_rate/r_P_norm,mp%A_rate/r_A_norm,mp%I_rate/r_I_norm
end subroutine

!--------------------------------------------------------------------------
! HIF1, C_G      -> r_G (analyticSetHIF1)
! PDK1, C_O2     -> r_Pm (analyticSetPDK1)
! r_Pm, C_L, r_G -> f_G, f_P, C_P (table lookup)
! f_G, f_P, C_P  -> r_P, r_A, r_I, r_L, r_O
! Time trials comparing optimiser with interpolator show that using
! lookup tables gives a speed improvement factor of only 9/8 - not worth it.
!--------------------------------------------------------------------------
subroutine run_optimiser
real(REAL_KIND) :: HIF1, fPDK, f_PO, f_PA, V, K1, K2, Km_O2, C_O2, C_G, r_G_max, C_L_max, r_Pm_base, r_Pm_max
real(REAL_KIND) :: MM_O2
real(REAL_KIND) :: r_G, r_Pm, C_L, alfa
real(REAL_KIND) :: x, y, C_P, ysum
real(REAL_KIND) :: r_P, r_A, r_I, r_O, r_L
real(REAL_KIND) :: dr_G, dC_L, dr_Pm
integer :: ityp, N_O2, N, i
integer :: i_r_G, i_C_L, i_r_Pm
type(metabolism_type), target :: met
type(metabolism_type), pointer :: mp

!if (allocated(f_G_lookup)) deallocate(f_G_lookup)
!if (allocated(f_P_lookup)) deallocate(f_P_lookup)
!if (allocated(C_P_lookup)) deallocate(C_P_lookup)
!allocate(f_G_lookup(N_r_G,N_C_L,N_r_Pm))
!allocate(f_P_lookup(N_r_G,N_C_L,N_r_Pm))
!allocate(C_P_lookup(N_r_G,N_C_L,N_r_Pm))

mp => met
V = Vcell_cm3
ityp = 1
met = metabolic
HIF1 = 0.5
C_O2 = 0.1
C_G = 0.0045
fPDK = 0.5
N_O2 = chemo(OXYGEN)%Hill_N
Km_O2 = chemo(OXYGEN)%MM_C0
K1 = K_PL
K2 = K_LP
f_PO = N_PO
f_PA = N_PA
MM_O2 = f_MM(C_O2,Km_O2,N_O2)
r_G_max = get_glycosis_rate(HIF1,C_G)
!r_Pm_max = fPDK*MM_O2*O2_maxrate/(f_PO*(1-f_P_norm))	! note that MM_P is not here, since it varies it is added in optimiser()
! This is the rate of oxidation of pyruvate, i.e. (1-f_P)*r_P, excluding the effect of MM(C_P)
!r_Pm_base = fPDK*MM_O2*(O2_maxrate-base_O_rate)/f_PO	! note that MM_P is not here, since it varies it is added in optimiser()
r_Pm_base = fPDK*MM_O2*O2_maxrate/f_PO	! note that MM_P is not here, since it varies it is added in optimiser()
C_L_max = 3.0

r_G = r_G_max
C_L = 3*C_L_max
r_Pm = r_Pm_base
C_P = C_P_norm		! initial guess
write(*,'(a,3e12.3)') 'r_G, C_L, r_Pm: ',r_G, C_L, r_Pm
write(*,*)
call optimiser(ityp, r_G, C_L, r_Pm, mp, x, y, C_P)
!mp%O_rate = mp%O_rate + MM_O2*base_O_rate
write(*,'(a,4f8.4)') 'optimiser: x, y, C_P: ',x,y,C_P
stop

r_Pm_max = r_Pm_base
dr_G = r_G_max/(N_r_G - 1)
dC_L = C_L_max/(N_C_L - 1)
dr_Pm = r_Pm_max/(N_r_Pm - 1)
do i_r_G = 1,N_r_G
	r_G_array(i_r_G) = (i_r_G - 1)*dr_G
enddo
do i_C_L = 1,N_C_L
	C_L_array(i_C_L) = (i_C_L - 1)*dC_L
enddo
do i_r_Pm = 1,N_r_Pm
	r_Pm_array(i_r_Pm) = (i_r_Pm - 1)*dr_Pm
enddo
C_L = C_L_max
r_G = r_G_max
N = 21
!do i = 1,N
!	alfa = (i-1)*1.0/(N-1)
!	C_L = alfa*C_L_max
do i_r_G = 1,N_r_G
do i_C_L = 1,N_C_L
do i_r_Pm = 1,N_r_Pm
	r_G = r_G_array(i_r_G)
	C_L = C_L_array(i_C_L)
	r_Pm = r_Pm_array(i_r_Pm)
	call optimiser(ityp, r_G, C_L, r_Pm, mp, x, y, C_P)
!	f_G_lookup(i_r_G,i_C_L,i_r_Pm) = x
!	f_P_lookup(i_r_G,i_C_L,i_r_Pm) = y
!	C_P_lookup(i_r_G,i_C_L,i_r_Pm) = C_P
!	write(*,'(a,4f8.4,e12.3)') 'x, y, C_L, C_P: ',x,y,C_L,C_P
!	write(*,'(a,4e12.3)') 'A, I, O, L: ',mp%A_rate, mp%I_rate, mp%O_rate, mp%L_rate
	r_P = 2*(1-x)*r_G - V*(K1*C_P - K2*C_L)
	r_A = 2*(1-x)*r_G + f_PA*(1-y)*r_P
	r_I = x*r_G + y*r_P
	r_O = f_PO*r_P
	r_L = V*(K1*C_P - K2*C_L)
!	write(*,'(a,4e12.3)') 'r_A, r_I: ',r_A,r_I,r_O,r_L
enddo
enddo
enddo

r_G = 0.9*r_G_max
C_L = 0.7*C_L_max
r_Pm = r_Pm_max
write(*,*) 'start optimiser'
ysum = 0
do i = 1,50000000
	call optimiser(ityp, r_G, C_L, r_Pm, mp, x, y, C_P)
	ysum = ysum + y
enddo
!write(*,'(a,4f8.4)') 'optimiser: x, y, C_P: ',x,y,C_P
write(*,*) 'done: ',ysum
write(*,*) 'start interpolator:'
ysum = 0
do i = 1,50000000
!	call interpolator(r_G, C_L, r_Pm, x, y, C_P)
	ysum = ysum + y
enddo
!write(*,'(a,4f8.4)') 'interpolator: x, y, C_P: ',x,y,C_P
write(*,*) 'done: ',ysum
end subroutine

#IF 0
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine interpolator(r_G, C_L, r_Pm, f_G, f_P, C_P)
real(REAL_KIND) :: r_G, C_L, r_Pm, f_G, f_P, C_P
integer :: i1, i2, i_r_G1, i_r_G2, i_C_L1, i_C_L2, i_r_Pm1, i_r_Pm2
real(REAL_KIND) :: a_r_G, a_C_L, a_r_Pm, fac(8), r_G_max, C_L_max, r_Pm_max
real(REAL_KIND) :: dr_G, dC_L, dr_Pm

dr_G = r_G_array(N_r_G)/(N_r_G - 1)
dC_L = C_L_array(N_C_L)/(N_C_L - 1)
dr_Pm = r_Pm_array(N_r_Pm)/(N_r_Pm - 1)

!i1 = 0
!do i2 = 1,N_r_G
!	if (r_G_array(i2) > r_G) then
!		i1 = i2 - 1
!		exit
!	endif
!enddo
!if (i1 > 0) then
i1 = r_G/dr_G
if (i1 < N_r_G) then
	i_r_G1 = i1
	i_r_G2 = i1+1
	a_r_G = (r_G - r_G_array(i1))/(r_G_array(i1+1) - r_G_array(i1))
else
	i_r_G1 = N_r_G
	i_r_G2 = N_r_G
	a_r_G = 1
endif

!write(*,*) 'r_G: ',i_r_G1,i_r_G2,a_r_G

!i1 = 0
!do i2 = 1,N_C_L
!	if (C_L_array(i2) > C_L) then
!		i1 = i2 - 1
!		exit
!	endif
!enddo
!if (i1 > 0) then
i1 = C_L/dC_L
if (i1 < N_r_G) then
	i_C_L1 = i1
	i_C_L2 = i1+1
	a_C_L = (C_L - C_L_array(i1))/(C_L_array(i1+1) - C_L_array(i1))
else
	i_C_L1 = N_C_L
	i_C_L2 = N_C_L
	a_C_L = 1
endif
!write(*,*) 'C_L: ',i_C_L1,i_C_L2,a_C_L

!i1 = 0
!do i2 = 1,N_r_Pm
!	if (r_Pm_array(i2) > r_Pm) then
!		i1 = i2 - 1
!		exit
!	endif
!enddo
!if (i1 > 0) then
i1 = r_Pm/dr_Pm
if (i1 < N_r_Pm) then
	i_r_Pm1 = i1
	i_r_Pm2 = i1+1
	a_r_Pm = (r_Pm - r_Pm_array(i1))/(r_Pm_array(i1+1) - r_Pm_array(i1))
else
	i_r_Pm1 = N_r_Pm
	i_r_Pm2 = N_r_Pm
	a_r_Pm = 1
endif
!write(*,*) 'r_Pm: ',i_r_Pm1,i_r_Pm2,a_r_Pm

!f_G = (1 - a_r_G)*(1 - a_C_L)*(1 - a_r_Pm)*f_G_lookup(i_r_G1,i_C_L1,i_r_Pm1) &
!    + (1 - a_r_G)*     a_C_L *(1 - a_r_Pm)*f_G_lookup(i_r_G1,i_C_L2,i_r_Pm1) &
!    + (1 - a_r_G)*(1 - a_C_L)*     a_r_Pm *f_G_lookup(i_r_G1,i_C_L1,i_r_Pm2) &
!    + (1 - a_r_G)*     a_C_L *     a_r_Pm *f_G_lookup(i_r_G1,i_C_L2,i_r_Pm2) &
!         + a_r_G *(1 - a_C_L)*(1 - a_r_Pm)*f_G_lookup(i_r_G2,i_C_L1,i_r_Pm1) &
!         + a_r_G *     a_C_L *(1 - a_r_Pm)*f_G_lookup(i_r_G2,i_C_L2,i_r_Pm1) &
!         + a_r_G *(1 - a_C_L)*     a_r_Pm *f_G_lookup(i_r_G2,i_C_L1,i_r_Pm2) &
!         + a_r_G *     a_C_L *     a_r_Pm *f_G_lookup(i_r_G2,i_C_L2,i_r_Pm2) 
fac(1) = (1 - a_r_G)*(1 - a_C_L)*(1 - a_r_Pm)
fac(2) = (1 - a_r_G)*     a_C_L *(1 - a_r_Pm)
fac(3) = (1 - a_r_G)*(1 - a_C_L)*     a_r_Pm 
fac(4) = (1 - a_r_G)*     a_C_L *     a_r_Pm 
fac(5) =      a_r_G *(1 - a_C_L)*(1 - a_r_Pm)
fac(6) =      a_r_G *     a_C_L *(1 - a_r_Pm)
fac(7) =      a_r_G *(1 - a_C_L)*     a_r_Pm 
fac(8) =      a_r_G *     a_C_L *     a_r_Pm 
f_G = fac(1)*f_G_lookup(i_r_G1,i_C_L1,i_r_Pm1) &
    + fac(2)*f_G_lookup(i_r_G1,i_C_L2,i_r_Pm1) &
    + fac(3)*f_G_lookup(i_r_G1,i_C_L1,i_r_Pm2) &
    + fac(4)*f_G_lookup(i_r_G1,i_C_L2,i_r_Pm2) &
    + fac(5)*f_G_lookup(i_r_G2,i_C_L1,i_r_Pm1) &
    + fac(6)*f_G_lookup(i_r_G2,i_C_L2,i_r_Pm1) &
    + fac(7)*f_G_lookup(i_r_G2,i_C_L1,i_r_Pm2) &
    + fac(8)*f_G_lookup(i_r_G2,i_C_L2,i_r_Pm2) 
f_P = fac(1)*f_P_lookup(i_r_G1,i_C_L1,i_r_Pm1) &
    + fac(2)*f_P_lookup(i_r_G1,i_C_L2,i_r_Pm1) &
    + fac(3)*f_P_lookup(i_r_G1,i_C_L1,i_r_Pm2) &
    + fac(4)*f_P_lookup(i_r_G1,i_C_L2,i_r_Pm2) &
    + fac(5)*f_P_lookup(i_r_G2,i_C_L1,i_r_Pm1) &
    + fac(6)*f_P_lookup(i_r_G2,i_C_L2,i_r_Pm1) &
    + fac(7)*f_P_lookup(i_r_G2,i_C_L1,i_r_Pm2) &
    + fac(8)*f_P_lookup(i_r_G2,i_C_L2,i_r_Pm2) 
C_P = fac(1)*C_P_lookup(i_r_G1,i_C_L1,i_r_Pm1) &
    + fac(2)*C_P_lookup(i_r_G1,i_C_L2,i_r_Pm1) &
    + fac(3)*C_P_lookup(i_r_G1,i_C_L1,i_r_Pm2) &
    + fac(4)*C_P_lookup(i_r_G1,i_C_L2,i_r_Pm2) &
    + fac(5)*C_P_lookup(i_r_G2,i_C_L1,i_r_Pm1) &
    + fac(6)*C_P_lookup(i_r_G2,i_C_L2,i_r_Pm1) &
    + fac(7)*C_P_lookup(i_r_G2,i_C_L1,i_r_Pm2) &
    + fac(8)*C_P_lookup(i_r_G2,i_C_L2,i_r_Pm2) 
end subroutine
#ENDIF

end module