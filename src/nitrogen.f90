! Model with r_Gln = (f_N/f_Gln)r_I

!==================================================================================================
!==================================================================================================
subroutine get_unconstrained_rates
real(REAL_KIND) :: f_NG, V, qu, MM_Gln, Km_Gln, L2
integer :: N_Gln

V = Vcell_cm3*average_volume		! should be actual cell volume cp%V
f_NG = f_N/f_Glnu
N_Gln = 1
Km_Gln = Hill_Km_Gln
MM_Gln = f_MM(C_Gln_norm,Km_Gln,N_Gln)
L2 = Gln_maxrate*MM_Gln
write(nflog,'(a,3f6.3)') 'f_N, f_Glnu, f_NG: ',f_N, f_Glnu, f_NG
r_Gu = G_maxrate*glucose_metab(C_G_norm)
qu = (f_NG*(1 - f_Glnu)*N_GlnO)/(1 - f_N*N_GlnI)
r_Pu = (O2_maxrate - qu*r_Gu*f_Gu*N_GI)/((1 - f_Pu)*N_PO + qu*f_Pu*N_PI)
r_Iu = (r_Gu*f_Gu*N_GI + r_Pu*f_Pu*N_PI)/(1 - f_N*N_GlnI) 
r_Glnu = f_NG*r_Iu
write(nflog,'(a,2e12.3)') 'r_Glnu, L2: ',r_Glnu,L2
if (r_Glnu > L2) then
    write(nflog,*) 'ERROR: unconstrained Gln rate exceeds constraint L2'    ! Can this be fixed?
    stop
endif
r_Au = r_Gu*(1 - f_Gu)*N_GA + r_Pu*(1 - f_Pu)*N_PA + r_Glnu*(1 - f_Glnu)*N_GlnA
r_Ou = r_Pu*(1 - f_Pu)*N_PO + r_Glnu*(1 - f_Glnu)*N_GlnO
C_Pu = ((1 - f_Gu)*r_Gu*N_GP + V*K_LP*C_L_norm - r_Pu)/(V*K_PL)
r_Lu = (1 - f_Gu)*r_Gu*N_GP - r_Pu
write(nflog,'(a,7e12.3)') 'unconstrained rates: A,I,G,Gln,P,L,O: ',r_Au,r_Iu,r_Gu,r_Glnu,r_Pu,r_Lu,r_Ou
end subroutine

!==================================================================================================
!==================================================================================================
subroutine set_param_set(r_G,C_L)
real(REAL_KIND) :: r_G, C_L
real(REAL_KIND) :: V, f_NG, q, fGlnA
!real(REAL_KIND) :: r_P1, r_P2, C_P, w

V = Vcell_cm3*average_volume
f_NG = f_N/f_Glnu
q = 1/(1 - f_NG*f_Glnu*N_GlnI)
fGlnA = 1 - f_Glnu

ps%h = (1 - f_ATPg)*r_Au/r_Iu
ps%a0 = r_G*N_GP + V*K_LP*C_L
ps%b0 = -r_G*f_Gu*N_GP
ps%c0 = -V*K_PL
ps%a1 = q*f_Pu*N_PI
ps%b1 = q*r_G*f_Gu*N_GI
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
real(REAL_KIND) :: w, a, b, c, d, aa, bb, cc, dd, w1, w2, f_NG, q, V
real(REAL_KIND) :: Km_O2, Km_Gln, MM_O2, MM_Gln, L1, L2, C_P, r_G, r_P, r_O2, r_Gln, r_A, r_I, r_L
integer :: N_O2, N_Gln, iw
logical :: done, clean, check

clean = .true.

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
if (r_P < 0) then
    clean = .false.
    r_P = 0
    q = r_G*f_Gu*N_GI/(1 - f_N*N_GlnI)
    w = min(L1/(q*f_NG*(1 - f_Glnu)*N_GlnO),L2/(q*f_NG))
!    write(nflog,'(a,f6.3)') 'In L1: set r_P=0 ==================> w: ',w
endif
C_P = (r_P - ps%a0 - ps%b0*w)/ps%c0
if (C_P < 0) then
    clean = .false.
    w = zero_C_P(a,b,c,d)
    C_P = 0
    r_P = r_G*(1 - w*f_Gu)*N_GP + V*K_LP*C_L
!    write(nflog,'(a,f8.3,e12.3)') 'In L1, C_P < 0, set w: ',w,r_P
endif
r_Gln = w*(r_P*ps%a2 + ps%b2)
if (w < 0) then
    clean = .false.
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
    if (r_Gln < L2) then
!        write(nflog,'(a,f8.3,2x,2e12.3)') 'L1 solution for w: ',w,r_Gln,L2
        done = .true.
        r_O2 = r_P*(ps%a4 + ps%b4*w) + ps%c4*w
        r_I = r_P*ps%a1*w + ps%b1*w
        r_A = r_P*(ps%a3 + ps%b3*w) + ps%c3 + ps%d3*w
        r_L = r_G*(1 - w*f_Gu)*N_GP - r_P
    else
!        write(nflog,'(a,2f8.3,2x,2e12.3)') 'L2 violated: ',w1,w2,r_Gln,L2
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
    if (w < 0) write(nflog,*) 'In L2: w < 0'
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
        clean = .false.
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
        clean = .false.
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
        write(nflog,'(a,2f8.3,2x,2e12.3)') 'L1 violated: ',w,r_O2,L1
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
    
    ! Check that clean really is recalcable
    check = .true.
    w = mp%f_G/f_Gu
    r_P = (a + b*w)/(c + d*w)
    if (r_P /= mp%P_rate) check = .false.
    r_O2 = r_P*(ps%a4 + ps%b4*w) + ps%c4*w
    if (r_O2 /= mp%O_rate) check = .false.  
    r_I = r_P*ps%a1*w + ps%b1*w
    if (r_I /= mp%I_rate) check = .false.
    r_A = r_P*(ps%a3 + ps%b3*w) + ps%c3 + ps%d3*w
    if (r_A /= mp%A_rate) check = .false.
    r_L = r_G*(1 - w*f_Gu)*N_GP - r_P
    if (r_L /= mp%L_rate) check = .false.
    if (check /= clean) then
        write(nflog,'(a,L4,5e12.3)') 'Bad check: ',check,r_P,r_O2,r_I,r_A,r_L
        write(nflog,'(a,L4,5e12.3)') '           ',clean,mp%P_rate,mp%O_rate,mp%I_rate,mp%A_rate,mp%L_rate
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
real(REAL_KIND) :: w, a, b, c, d, aa, bb, cc, dd, w1, w2, f_NG, q, V
real(REAL_KIND) :: Km_O2, Km_Gln, MM_O2, MM_Gln, L1, L2, C_P, r_G, r_P, r_O2, r_Gln, r_A, r_I, r_L
integer :: N_O2, N_Gln, iw
logical :: done, clean, check

clean = .true.

r_G = get_glycosis_rate(mp%HIF1,C_G,C_Gln,mp%O_rate)  ! Note: this is the previous O_rate
call set_param_set(r_G,C_L)

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




