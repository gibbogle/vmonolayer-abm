! To determine distribution of colony size

module colony

use global
use cellstate
implicit none

!integer, parameter :: n_colony_days=10
integer, parameter :: max_trials = 1000
integer :: nmax
integer, allocatable :: perm_index(:)
logical :: use_permute

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
real(REAL_KIND) :: V0, dVdt, dt, t, tend
real(REAL_KIND) :: tnow_save
integer :: k, kcell, ityp, n, idist, ncycmax, ntot, nlist_save, ntrials, ndays
type (cell_type), pointer :: cp
logical :: ok

write(logmsg,'(a,f8.1)') 'make_colony_distribution: n_colony_days: ',n_colony_days
call logger(logmsg)
colony_simulation = .true.
ndays = (Nsteps*DELTA_T)/(24.*60*60)
nlist_save = nlist
tnow_save = tnow
ncycmax = 24*3600*n_colony_days/divide_time_mean(1) + 3
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
write(logmsg,*) 'Number of trials: ',ntrials
call logger(logmsg)
if (n_colony_days <= 4) then
    ddist = 1
elseif (n_colony_days <= 5) then
    ddist = 2
elseif (n_colony_days <= 6) then
    ddist = 5
elseif (n_colony_days <= 7) then
    ddist = 10
elseif (n_colony_days <= 8) then
    ddist = 20
else
    ddist = 50
endif
dist(1:ndist) = 0
ntot = 0
tend = tnow + n_colony_days*24*3600    ! plate for 10 days
!do kcell = 1, nlist_save
do k = 1, ntrials
    kcell = perm_index(k)
	cp => cell_list(kcell)
!	if (cp%state == DEAD) cycle
	ityp = cp%celltype
	
	! Now simulate colony growth from a single cell
	tnow = tnow_save
	call make_colony(kcell,tend,n)
	if (mod(k,100) == 0) then
	    write(logmsg,*) 'cell: n: ',k,n
	    call logger(logmsg)
	endif
!	if (n < 50) then
!	    write(*,*) 'small colony: ',kcell,n
!	    stop
!	endif
	ntot = ntot + n
	idist = n/ddist + 1
	dist(idist) = dist(idist) + 1
enddo 
dist(1:ndist) = dist(1:ndist)/sum(dist(1:ndist))
write(logmsg,'(a,2i8,f8.1)') 'Colony size distribution: ', nlist_save,ntot,real(ntot)/nlist_save
call logger(logmsg)
write(nfout,'(a,2i8,f8.1)') 'Colony size distribution: ', nlist_save,ntot,real(ntot)/nlist_save
do idist = 1,ndist
	write(logmsg,'(i4,a,i4,f6.3)') int((idist-1)*ddist),'-',int(idist*ddist),dist(idist)
    call logger(logmsg)
	write(nfout,'(i4,a,i4,f6.3)') int((idist-1)*ddist),'-',int(idist*ddist),dist(idist)
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
ccell_list(1)%anoxia_tag = .false.
ccell_list(1)%aglucosia_tag = .false.
ccell_list(1)%drug_tag = .false.
dt = DELTA_T

!write(*,*) ccell_list(1)%dVdt,max_growthrate(1),ccell_list(1)%G2_time
nlist = 1
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

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine make_perm_index(ok)
logical :: ok
integer :: np, kcell, kpar=0

np = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	np = np + 1
	perm_index(np) = kcell
enddo
if (np /= ncells) then
	write(logmsg,*) 'Error: make_perm_index: np /= Ncells: ',np,ncells,nlist
	call logger(logmsg)
	ok = .false.
	return
endif
if (use_permute) then
	call permute(perm_index,np,kpar)
endif
ok = .true.
end subroutine

end module
