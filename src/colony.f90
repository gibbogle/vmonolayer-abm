! To determine distribution of colony size

module colony

use global
use cellstate
implicit none

!integer, parameter :: n_colony_days=10
integer, parameter :: max_trials = 1000
integer :: nmax
!integer, allocatable :: perm_index(:)
!logical :: use_permute

contains

!---------------------------------------------------------------------------------------------------
! Simulate fate of cells grown with no nutrient constraints.
! The only determinants of the colony size for a cell are (considering radiation only):
! volume
! divide_time_mean(ityp)
! radiation_tag
! p_rad_death
! growth_delay
! G2_M
! dt_delay
! t_growth_delay_end
! N_delayed_cycles_left
! The new method simply continues the simulation from where it ended, for 10 days. 
!---------------------------------------------------------------------------------------------------
subroutine make_colony_distribution(n_colony_days,dist,ddist,ndist) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: make_colony_distribution
use, intrinsic :: iso_c_binding
real(c_double) :: n_colony_days, dist(*), ddist
integer(c_int) :: ndist
integer, parameter :: ddist50 = 50
integer, parameter :: ndist50 = 1000/ddist50
real(REAL_KIND) :: V0, dVdt, dt, t, tend, sum1, sum2, SD, SE, ave, dist50(ndist50), dmin
real(REAL_KIND) :: tnow_save
integer :: k, kk, kcell, ityp, n, idist, ncycmax, ntot, nlist_save, ntrials, ndays, nt, idist50, kmin
type (cell_type), pointer :: cp
logical :: ok
integer :: dist_cutoff = 200

colony_simulation = .true.
ndays = (Nsteps*DELTA_T)/(24.*60*60)
nlist_save = nlist
tnow_save = tnow
ncycmax = 24*3600*n_colony_days/divide_time_mean(1) + 1
nmax = 2**ncycmax
allocate(ccell_list(nmax))
if (allocated(perm_index)) deallocate(perm_index)
allocate(perm_index(nlist_save))
if (Ncells > max_trials) then
    use_permute = .true.
    ntrials = max_trials
else
    use_permute = .false.
    ntrials = Ncells
endif
call make_perm_index(ok)
if (.not.ok) then
    call logger('Error: make_perm_index')
    dist(1:ndist) = 0
    return
endif
ddist = nmax/ndist
dmin = 1.0e10
do k = 1,ndist
    if (abs(k*100 - ddist) < dmin) then
        dmin = abs(k*100 - ddist)
        kmin = k
    endif
enddo
ddist = kmin*100
write(nfout,*)
write(logmsg,'(a,f8.1,2i5)') 'make_colony_distribution: n_colony_days: ',n_colony_days
call logger(logmsg)
write(nfout,'(a,f8.1,2i5)') 'make_colony_distribution: n_colony_days: ',n_colony_days
write(nfout,'(a,i5,f8.1,i5)') 'ndist,ddist,ntrials: ',ndist,ddist,ntrials
dist(1:ndist) = 0
dist50 = 0
ntot = 0
sum1 = 0
sum2 = 0
tend = tnow + n_colony_days*24*3600    ! plate for 10 days
kk = 0
k = 0
nt = 0  ! count of runs giving colony size n > dist_threshold
do while(k < ntrials)
    kk = kk+1
    kcell = perm_index(kk)
	cp => cell_list(kcell)
	if (cp%state == DEAD) then
	    write(nflog,*) 'colony: cell dead: ',kcell
	    cycle
	else
	    k = k+1
	endif
	ityp = cp%celltype
	! Now simulate colony growth from a single cell
	tnow = tnow_save
	call make_colony(kcell,tend,n)
	ntot = ntot + n
	if (n > dist_cutoff) then
	    nt = nt + 1
	    sum1 = sum1 + n
    	sum2 = sum2 + n*n
    endif
	idist = n/ddist + 1
	dist(idist) = dist(idist) + 1
	if (mod(k,100) == 0) then
	    write(logmsg,'(a,3i8)') 'cell: n, idist: ',k,n,idist
	    call logger(logmsg)
	endif
	idist50 = n/ddist50 + 1
	if (idist50 <= ndist50) then
	    dist50(idist50) = dist50(idist50) + 1
	endif
enddo

ave = sum1/nt
SD = sqrt((sum2 - nt*ave**2)/(nt-1))
SE = SD/sqrt(real(nt))
write(nfout,*)
write(nfout,'(a,i4,a,i5)') 'With cutoff colony size: ',dist_cutoff, ' number of colonies: ',nt
write(nfout,'(a,3f8.1)') 'average size, SD, SE: ',ave,SD,SE
write(nfout,*)
dist50(1:ndist50) = dist50(1:ndist50)/ntrials
write(nfout,'(a,2i8,f8.1)') 'Colony size distribution < 1000:'
do idist50 = 1,ndist50
	write(nfout,'(i6,a,i6,f7.4)') int((idist50-1)*ddist50),'-',int(idist50*ddist50),dist50(idist50)
enddo
write(nfout,*)

dist(1:ndist) = dist(1:ndist)/ntrials
write(logmsg,'(a,2i8,f8.1)') 'Colony size distribution: ', nlist_save,ntot,real(ntot)/nlist_save
call logger(logmsg)
write(nfout,'(a,2i8,f8.1)') 'Colony size distribution:'
do idist = 1,ndist
	write(logmsg,'(i6,a,i6,f7.4)') int((idist-1)*ddist),'-',int(idist*ddist),dist(idist)
    call logger(logmsg)
	write(nfout,'(i6,a,i6,f7.4)') int((idist-1)*ddist),'-',int(idist*ddist),dist(idist)
enddo
deallocate(ccell_list)
deallocate(perm_index)

colony_simulation = .false.
nlist = nlist_save
tnow = tnow_save

end subroutine

!---------------------------------------------------------------------------------------------------
! The cell is at the point of division - possibly G2_M (arrested at G2/M checkpoint)
! For now only radiation tagging is handled
! Growth rate dVdt (mean) is used only to estimate the time of next division
!---------------------------------------------------------------------------------------------------
subroutine make_colony(kcell,tend,n)
integer :: kcell, n
real(REAL_KIND) :: tend, dt 
integer :: icell, ityp, nlist0, kpar=0
real(REAL_KIND) :: V0, Tdiv0, r_mean, c_rate, dVdt, Tmean, R
logical :: changed, ok
type (cell_type), pointer :: cp

!write(*,'(a,i6,2f12.0)') 'make_colony: ',kcell,tnow,tend
ccell_list(1) = cell_list(kcell)
cp => ccell_list(1)
!ccell_list(1)%anoxia_tag = .false.
!ccell_list(1)%aglucosia_tag = .false.
!ccell_list(1)%ATP_tag = .false.
!ccell_list(1)%radiation_tag = .false.
!ccell_list(1)%drug_tag = .false.
ityp = cp%celltype
if (cp%state /= DYING) cp%dVdt = max_growthrate(ityp)
dt = DELTA_T
!write(*,*) ccell_list(1)%dVdt,max_growthrate(1),ccell_list(1)%G2_time
nlist = 1
ncells = 1
ngaps = 0
do while (tnow < tend)
	tnow = tnow + dt
    call new_grower(dt,changed,ok)
enddo
!write(*,*) nlist,ccell_list(nlist)%dVdt,max_growthrate(1),ccell_list(nlist)%V
n = 0
do icell = 1,nlist
	if (ccell_list(icell)%state /= DEAD) n = n+1
enddo
end subroutine

end module
