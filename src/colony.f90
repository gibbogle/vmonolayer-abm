! To determine distribution of colony size

module colony

use global
use cellstate
implicit none

!integer, parameter :: n_colony_days=10
integer, parameter :: max_trials = 20000     ! 20000
!integer, allocatable :: perm_index(:) 
!logical :: use_permute
real(REAL_KIND), parameter :: min_log10PE = -4.5    ! -4.0
! For daily colony size distributions
real(REAL_KIND) :: bin_size
integer, allocatable :: bin_count(:,:)
integer :: nbins
integer :: kcell0

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
!
! When the number of surviving cells at the start of the colony simulation is very small, the results
! of the simulation very greatly with the random seed.  Need to perform multiple simulations with
! these cells to reduce the variability.
! First create a list of the surviving cells.
!---------------------------------------------------------------------------------------------------
subroutine make_colony_distribution(n_colony_days,dist,ddist,ndist,PE) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: make_colony_distribution
use, intrinsic :: iso_c_binding
real(c_double) :: n_colony_days, dist(*), ddist, PE
integer(c_int) :: ndist
integer, parameter :: ddist50 = 50
integer, parameter :: ndist50 = 1000/ddist50
real(REAL_KIND) :: V0, dVdt, dt, t, tend, sum1, sum2, SD, SE, ave, dist50(ndist50), dmin, log10PE
real(REAL_KIND) :: tnow_save
integer :: k, kk, kcell, ityp, n, idist, ncycmax, ntot, nlist_save, ntrials, ndays, nt, idist50, kmin, nrepeat, krep, kp, kpmax
type (cell_type), pointer :: cp
logical :: ok
integer :: dist_cutoff = 50
integer, allocatable :: ncolony(:),ntcolony(:)
integer, allocatable :: survivor(:)
integer :: iday

nbins = 20
bin_size = 10
if (allocated(bin_count)) deallocate(bin_count)
if (allocated(ccell_list)) deallocate(ccell_list)
if (allocated(survivor)) deallocate(survivor)
if (allocated(perm_index)) deallocate(perm_index)
if (allocated(ncolony)) deallocate(ncolony)
if (allocated(ntcolony)) deallocate(ntcolony)

allocate(bin_count(10,0:nbins+1))
bin_count = 0
simulate_colony = .true.
colony_simulation = .true.
ndays = int(n_colony_days)
nlist_save = nlist
tnow_save = tnow
ncycmax = 24*3600*ndays/divide_time_mean(1) + 1
nColonyMax = 2**(ncycmax+1)
write(nflog,'(a,e12.3,2i8)') 'divide_time_mean(1),ncycmax,nColonyMax: ',divide_time_mean(1),ncycmax,nColonyMax
allocate(ccell_list(nColonyMax))
allocate(survivor(Ncells))
allocate(perm_index(nlist_save))
!if (Ncells > max_trials) then
!    use_permute = .true.
!    ntrials = max_trials
!else
!    use_permute = .false.
!    ntrials = Ncells
!endif
!call make_perm_index(ok)
!if (.not.ok) then
!    call logger('Error: make_perm_index')
!    dist(1:ndist) = 0
!    return
!endif
k = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) then
	    cycle
	else
	    k = k+1
	    if (k > Ncells) then
	        write(*,*) 'Error: dimension of survivor(:) exceeded'
	        stop
	    endif
	    survivor(k) = kcell
	endif
enddo
if (k /= Ncells) then
    write(*,*) 'Error: inconsistent survivor numbers: ',Ncells, k
    stop
endif
! At this point survivor(:) holds the cell IDs of all surviving cells
write(nflog,*) 'Created survivor list'
if (Ncells > max_trials) then
    nrepeat = 1
    ntrials = max_trials
    use_permute = .true.
    call make_perm_index(ok)
    if (.not.ok) then
        call logger('Error: make_perm_index')
        dist(1:ndist) = 0
        return
    endif
    kpmax = 0
    do kk = 1,Ncells
        kpmax = max(kpmax,perm_index(kk))
    enddo
else
    nrepeat = real(max_trials)/Ncells
    ! This is the number of repeat simulations per survivor
    ntrials = nrepeat*Ncells
    use_permute = .false.
endif
write(nflog,*) 'Ncells,nrepeat,ntrials: ',Ncells,nrepeat,ntrials
ddist = (nColonyMax)/ndist
dmin = 1.0e10
do k = 1,ndist
    if (abs(k*100 - ddist) < dmin) then
        dmin = abs(k*100 - ddist)
        kmin = k
    endif
enddo
ddist = kmin*100
allocate(ncolony(0:ndays))
allocate(ntcolony(0:ndays))
write(logmsg,'(a,i4,i6)') 'make_colony_distribution: ndays, ntrials: ',ndays,ntrials
call logger(logmsg)
if (.not.use_PEST) then
    write(nfout,*)
    write(nfout,'(a,i4)') 'make_colony_distribution: ndays: ',ndays
    write(nfout,'(a,i5,f8.1,i6)') 'ndist,ddist,ntrials: ',ndist,ddist,ntrials
endif
dist(1:ndist) = 0
dist50 = 0
ntot = 0
sum1 = 0
sum2 = 0
tend = tnow + ndays*24*3600    ! plate for 10 days
kk = 0
k = 0
nt = 0  ! count of runs giving colony size n > dist_threshold
ntcolony = 0
!do while(k < ntrials)
!    kk = kk+1
!    kcell = perm_index(kk)
!	cp => cell_list(kcell)
!	if (cp%state == DEAD) then
!	    cycle
!	else
!	    k = k+1
!	endif
!	write(nflog,*) 'colony: ',k,kk,kcell
do kk = 1,min(max_trials,Ncells)
    if (use_permute) then
!        kp = perm_index(kk)
!        kcell = survivor(kp)
        kcell = perm_index(kk)
    else
        kcell = survivor(kk)
    endif
    cp => cell_list(kcell)
	ityp = cp%celltype
do krep = 1,nrepeat
    k = k+1
	! Now simulate colony growth from a single cell 
	tnow = tnow_save
	ncolony = 0
	kcell0 = kcell
	call make_colony(kcell,tend,ncolony)
!	if (ncolony(4) > 10 .and. ncolony(4) < 20) write(nflog,'(a,i8,10i6)') 'kcell0,ncolony: ',kcell0,ncolony(:)
	n = ncolony(ndays)
!	write(*,'(a,i6,i3)') 'trial, n: ',krep,n
!	if (n == 0) stop
!	if (ncolony(5) > 0 .and. ncolony(5) < 50) write(*,*) 'ncolony(5): ',ncolony(5)	
	ntcolony = ntcolony + ncolony
	ntot = ntot + n
	if (n > dist_cutoff) then
	    nt = nt + 1
	    sum1 = sum1 + n
    	sum2 = sum2 + n*n
    endif
	idist = n/ddist + 1
	dist(idist) = dist(idist) + 1
	if (mod(k,ntrials/10) == 0) then
	    write(logmsg,'(a,3i8)') 'cell: n, idist: ',k,n,idist
	    call logger(logmsg)
	endif
	idist50 = n/ddist50 + 1
	if (idist50 <= ndist50) then
	    dist50(idist50) = dist50(idist50) + 1
	endif
enddo
enddo
!	write(*,*) 'stopping'
!	stop

if (nt > 0) then
    ave = sum1/nt
    SD = sqrt((sum2 - nt*ave**2)/(nt-1))
    SE = SD/sqrt(real(nt))
else
    ave = 0
    SD = 0
    SE = 0
endif
PE = real(nt)/ntrials
dist(1:ndist) = dist(1:ndist)/ntrials
dist50(1:ndist50) = dist50(1:ndist50)/ntrials

if (.not.use_PEST) then
    write(nfout,*)
    write(nfout,'(a,i4,a,i5)') 'With cutoff colony size: ',dist_cutoff, ' number of colonies: ',nt
    write(nfout,'(a,3f8.1)') 'average size, SD, SE: ',ave,SD,SE
    write(nfout,*)
    write(nfout,'(a,2i8,f8.1)') 'Colony size distribution < 1000:'
    do idist50 = 1,ndist50
	    write(nfout,'(i6,a,i6,f7.4)') int((idist50-1)*ddist50),'-',int(idist50*ddist50),dist50(idist50)
    enddo
    write(nfout,*)
    write(nfout,'(a)') 'Colony size distribution:'
    do idist = 1,ndist
	    write(nfout,'(i6,a,i6,f7.4)') int((idist-1)*ddist),'-',int(idist*ddist),dist(idist)
    enddo
    write(nfout,*)
    write(nfout,*) 'Total colony population multiplication factor by day:'
    write(nfout,'(10f9.2)') ntcolony(1:ndays)/real(ntrials)
    write(nfout,*) 'Daily colony size distributions'
    do iday = 1,ndays
        write(nfout,'(a,i3)') 'Day: ',iday
        write(nfout,'(22i6)') bin_count(iday,:)
    enddo
endif
!if (PE > 0.000001) then
!    log10PE = log10(PE)
!else
!    log10PE = -4
!endif
log10PE = max(log10(PE),min_log10PE)
write(nfout,'(a,f8.4)') 'log10PE: ',log10PE


write(logmsg,'(a,2i8,f8.1)') 'Colony size distribution < 1000:' 
call logger(logmsg)
do idist50 = 1,ndist50
    write(logmsg,'(i6,a,i6,f7.4)') int((idist50-1)*ddist50),'-',int(idist50*ddist50),dist50(idist50)
    call logger(logmsg)
enddo
write(logmsg,'(a,2i12,f8.1)') 'Colony size distribution: ', nlist_save,ntot,real(ntot)/nlist_save
call logger(logmsg)
do idist = 1,ndist
	write(logmsg,'(i6,a,i6,f7.4)') int((idist-1)*ddist),'-',int(idist*ddist),dist(idist)
    call logger(logmsg)
enddo
write(logmsg,'(a)') 'Total colony population multiplication factor by day:'
call logger(logmsg)
write(logmsg,'(10f9.2)') ntcolony(1:ndays)/real(ntrials)
call logger(logmsg)
write(logmsg,'(a,2f8.4)') 'Plating efficiency PE, log10PE: ',PE,log10PE
call logger(logmsg)
call logger("------------------------------------")

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
subroutine make_colony(kcell,tend,ncolony)
integer :: kcell, ncolony(0:)
real(REAL_KIND) :: t0, tend, dt 
integer :: icell, ityp, nlist0, iday, n, k, kpar=0
real(REAL_KIND) :: V0, Tdiv0, r_mean, c_rate, dVdt, Tmean, R, totgen
logical :: changed, ok
type (cell_type), pointer :: cp
integer :: cnt(3)

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
!if (cp%ID == 1) write(*,*) 'ID=1: N_PL,N_IRL,N_Ch1,N_Ch2: ',cp%N_PL,cp%N_IRL,cp%N_Ch1,cp%N_Ch2
nlist = 1
ncells = 1
ngaps = 0
iday = 1
t0 = tnow
do while (tnow < tend)
	tnow = tnow + dt
    call new_grower(dt,changed,ok)
    call CellDeath(dt,ok)
    if (tnow >= t0 + iday*24*3600) then
!        if (cp%ID == 1) write(*,*) 'ID=1: nlist: ',nlist
        n = countColony()   ! counts cells in the colony
!        write(*,*) 'iday, n: ',iday,n,tnow,tend
!        write(*,*) 'ncells, nlist, n: ',ncells, nlist, n
        ncolony(iday) = n
        call binner(iday,n)
        iday = iday+1
    endif
enddo
ncolony(iday) = countColony()
!write(*,*) 'ending day, n, nlist: ',iday, ncolony(iday), nlist
return

!write(*,*) nlist,ccell_list(nlist)%dVdt,max_growthrate(1),ccell_list(nlist)%V 
cnt = 0
totgen = 0
n = 0
do icell = 1,nlist
    k = ccell_list(icell)%state
    cnt(k) = cnt(k) + 1
	if (ccell_list(icell)%state /= DEAD) then
	    totgen = totgen + ccell_list(icell)%generation
	endif
enddo
n = cnt(1) + cnt(2)
!write(*,*) 'Cell counts, average generation: ',cnt,totgen/n
end subroutine

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
function countColony() result(n)
integer :: n
integer :: icell, k, cnt(3)

cnt = 0
n = 0
do icell = 1,nlist
    k = ccell_list(icell)%state
	if (k /= DEAD) then
	    n = n+1
!	    write(*,*) 'cell state: ',ccell_list(icell)%state
	endif
	cnt(k) = cnt(k) + 1
enddo
!write(*,*) 'cnt: ',cnt
end function

!---------------------------------------------------------------------------------------------------
! Increment the bin count for day iday, with colony count n
!---------------------------------------------------------------------------------------------------
subroutine binner(iday,n)
integer :: iday, n
integer :: ibin

if (n == 0) then
    bin_count(iday,0) = bin_count(iday,0) + 1
elseif (n > bin_size*nbins) then
    bin_count(iday,nbins+1) = bin_count(iday,nbins+1) + 1
else
    ibin = int(n/bin_size) + 1
    bin_count(iday,ibin) = bin_count(iday,ibin) + 1
endif
end subroutine

end module
