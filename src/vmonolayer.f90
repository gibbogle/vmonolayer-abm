module monolayer_mod
use global
use chemokine
use ode_diffuse
use cellstate
use winsock  
use colony
use transfer
use metabolism
!use Tcp_mod

IMPLICIT NONE

contains 

!-----------------------------------------------------------------------------------------
! This subroutine is called to initialize a simulation run. 
! ncpu = the number of processors to use
! infile = file with the input data
! outfile = file to hold the output 
!-----------------------------------------------------------------------------------------
subroutine Setup(ncpu,infile,outfile,ok)
integer :: ncpu
character*(*) :: infile, outfile
logical :: ok
character*(64) :: msg
integer :: ichemo, error, kcell, idrug, ityp
real(REAL_KIND) :: tgrowth(MAX_CELLTYPES)
type(cycle_parameters_type),pointer :: ccp
type(metabolism_type), pointer :: mp

ok = .true.
initialized = .false.
par_zig_init = .false.
colony_simulation = .false.

inputfile = infile
outputfile = outfile
call logger("ReadCellParams new")
call ReadCellParams(ok)
if (.not.ok) return
call logger("did ReadCellParams")

start_wtime = wtime()

if (ncpu == 0) then
	ncpu = ncpu_input
endif
Mnodes = ncpu
write(logmsg,*) 'ncpu: ',ncpu 
call logger(logmsg)

#if defined(OPENMP) || defined(_OPENMP)
    call logger("OPENMP defined")
    call omp_initialisation(ok)
    if (.not.ok) return
#else
    call logger("OPENMP NOT defined")
    if (Mnodes > 1) then
        write(logmsg,'(a)') 'No OpenMP, using one thread only'
        call logger(logmsg)
        Mnodes = 1
    endif
#endif

! Set up grid alignment
!NY = NX
!NZ = NX
!x0 = (NX + 1.0)/2.
!y0 = (NY + 1.0)/2.
!z0 = (NZ + 1.0)/2.
!blob_centre = [x0,y0,z0]   ! (units = grids)
!
!DXB = 1.0e-4*DXB	! um -> cm
!ixb0 = (1 + NXB)/2
!iyb0 = (1 + NYB)/2
!izb0 = 6
!xb0 = (ixb0-1)*DXB
!yb0 = (iyb0-1)*DXB 
!zb0 = (izb0-1)*DXB
!centre_b = [xb0, yb0, zb0]
!dxb3 = dxb*dxb*dxb
!
!write(nflog,'(a,3e12.3)') 'x0,y0,z0: ',x0,y0,z0
!write(nflog,'(a,3e12.3)') 'xb0,yb0,zb0: ',xb0,yb0,zb0
!grid_offset(1) = (ixb0-1)*DXB - ((NX+1)/2)*DELTA_X
!grid_offset(2) = (iyb0-1)*DXB - ((NY+1)/2)*DELTA_X
!grid_offset(3) = (izb0-1)*DXB - ((NZ+1)/2)*DELTA_X
!write(nflog,'(a,3e12.3)') 'grid_offset: ',grid_offset
!write(nflog,'(a,3e12.3)') 'blob_centre: ',blob_centre*DELTA_X + grid_offset

call ArrayInitialisation(ok)
if (.not.ok) return
call logger('did ArrayInitialisation')

if (use_metabolism) then
	use_cell_cycle = .true.
	chemo(OXYGEN)%controls_growth = .false.
	chemo(OXYGEN)%controls_death = .false.
	chemo(GLUCOSE)%controls_growth = .false.
	chemo(GLUCOSE)%controls_death = .false.
endif

call SetupChemo

! New cell cycle formulation - need a value for max (unconstrained) growth rate
use_volume_method = .not.use_cell_cycle
!if (use_cell_cycle .and. .not.use_volume_based_transition) then
!    use_constant_growthrate = .true.
!endif
! Growth occurs during G1, S and G2, not in checkpoints
do ityp = 1,2
	ccp => cc_parameters(ityp)
	tgrowth = ccp%T_G1 + ccp%T_S + ccp%T_G2
	max_growthrate = Vdivide0/(2*tgrowth)
	write(nflog,*) 'ityp, Vdivide0, max_growthrate: ',ityp,Vdivide0, max_growthrate
enddo

!is_dropped = .false.
!adrop = 1
!bdrop = 1
!cdrop = 0
!zmin = 1
mp => phase_metabolic(1)
call SetupMetabolism(mp,ok)
call PlaceCells(ok)
call setTestCell(kcell_test)
!call show_volume_data
!call SetRadius(Nsites)
!call getVolume(blob_volume,blob_area)
!blob_radius = sqrt(blob_area/PI)
!blob_centre = getCentre()
write(logmsg,*) 'did PlaceCells: Ncells: ',Ncells
call logger(logmsg)
if (.not.ok) return

istep = 0
do ichemo = 1,TRACER
	if (chemo(ichemo)%used) then
		call InitConcs(ichemo)
		call SetupMedium(ichemo)
	endif
enddo
call UpdateChemomap
!call AdjustMM
call SetInitialGrowthRate
NATP_tag = 0
Nradiation_tag = 0
Ndrug_tag = 0
Ndrug_tag = 0
Nradiation_dead = 0
Ndrug_dead = 0
NATP_dead = 0
ndivided = 0
Ndying = 0
Ndead = 0

ndoublings = 0
doubling_time_sum = 0
ncells_mphase = 0

!radiation_dosed = .false.
t_simulation = 0
total_dMdt = 0
chemo(:)%total_flux_prev = 0
t_lastmediumchange = 0
medium_change_step = .false.
limit_stop = .false.
!Vex_min = 1.0
!Vex_max = 0
kcell_dbug = 0
write(logmsg,'(a,i6)') 'Startup procedures have been executed: initial T cell count: ',Ncells0
call logger(logmsg)
call averages
!if (is_radiation) then
!	ccp => cc_parameters(selected_celltype)
!	call logger('makeTCPradiation')
!	call makeTCPradiation(selected_celltype,NTCP) ! set checkpoint repair time limits 
!	call logger('did makeTCPradiation')
!else
!	ccp%tcp = 0
!endif
end subroutine

!----------------------------------------------------------------------------------------- 
!----------------------------------------------------------------------------------------- 
subroutine show_volume_data
integer :: kcell
real(REAL_KIND) :: Vsum, Vdivsum

write(nfout,*) 'Volume data:'
write(nfout,'(a,L)') 'use_divide_time_distribution: ',use_divide_time_distribution
write(nfout,'(a,L)') 'use_V_dependence: ',use_V_dependence
Vsum = 0
Vdivsum = 0
do kcell = 1,nlist
	write(nfout,'(i6,2f6.2)') kcell,cell_list(kcell)%V, cell_list(kcell)%divide_volume
	Vsum = Vsum + cell_list(kcell)%V
	Vdivsum = Vdivsum + cell_list(kcell)%divide_volume
enddo
write(nfout,*)
write(nfout,'(a,f6.2)') 'Average initial volume: ',Vsum/nlist
write(nfout,'(a,f6.2)') 'Average divide volume: ',Vdivsum/nlist
end subroutine

!----------------------------------------------------------------------------------------- 
! Initialise medium concentration
!-----------------------------------------------------------------------------------------
subroutine SetupMedium(ichemo)
integer :: ichemo, im

if (chemo(ichemo)%present) then
    Caverage(MAX_CHEMO+ichemo) = chemo(ichemo)%bdry_conc
    Cmediumave(ichemo) = chemo(ichemo)%bdry_conc
else
    Caverage(MAX_CHEMO+ichemo) = 0
    Cmediumave(ichemo) = 0
endif
if (ichemo <= 3) then
    C_OGL(ichemo,:) = chemo(ichemo)%bdry_conc
endif
if (ichemo == GLUCOSE) then
!	Cglucose = chemo(ichemo)%bdry_conc
	chemo(ichemo)%Cmedium = chemo(ichemo)%bdry_conc
endif
if (ichemo >= DRUG_A) then
!	im = ichemo - DRUG_A
!	Cdrug(im,:) = 0
	chemo(ichemo)%Cmedium = 0
endif
end subroutine

!----------------------------------------------------------------------------------------- 
!-----------------------------------------------------------------------------------------
subroutine omp_initialisation(ok)
logical :: ok
integer :: npr, nth

ok = .true.
!if (Mnodes == 1) return
#if defined(OPENMP) || defined(_OPENMP)
write(logmsg,'(a,i2)') 'Requested Mnodes: ',Mnodes
call logger(logmsg)
npr = omp_get_num_procs()
write(logmsg,'(a,i2)') 'Machine processors: ',npr
call logger(logmsg)

nth = omp_get_max_threads()
write(logmsg,'(a,i2)') 'Max threads available: ',nth
call logger(logmsg)
if (nth < Mnodes) then
    Mnodes = nth
    write(logmsg,'(a,i2)') 'Setting Mnodes = max thread count: ',nth
	call logger(logmsg)
endif

call omp_set_num_threads(Mnodes)
!$omp parallel
nth = omp_get_num_threads()
write(logmsg,*) 'Threads, max: ',nth,omp_get_max_threads()
call logger(logmsg)
!$omp end parallel
#endif

call logger('did omp_initialisation')
!call test_omp1

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_omp
integer, parameter :: n = 10
integer :: i

integer :: sum1, sum2
integer, allocatable :: y1(:)
integer :: y2(n)

allocate(y1(n))
y1 = 1
y2 = 1

sum1 = 0
sum2 = 0
!$omp parallel do
do i = 1,n
	sum1 = sum1 + y1(i)
	sum2 = sum2 + y2(i)
enddo
!$omp end parallel do
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ArrayInitialisation(ok)
logical :: ok
integer :: x,y,z,k, ichemo
integer :: MAXX, z1, z2, nc0, inflow
integer :: cog_size
real(REAL_KIND) :: d, rr(3)

ok = .false.
call RngInitialisation

! These are deallocated here instead of in subroutine wrapup so that when a simulation run ends
! it will still be possible to view the cell distributions and chemokine concentration fields.
!if (allocated(occupancy)) deallocate(occupancy)
if (allocated(cell_list)) deallocate(cell_list)
!if (allocated(allstate)) deallocate(allstate)
!if (allocated(ODEdiff%ivar)) deallocate(ODEdiff%ivar)
if (allocated(gaplist)) deallocate(gaplist)
!if (allocated(Cslice)) deallocate(Cslice)
call logger('did deallocation')

!nsteps_per_min = 1.0/DELTA_T
ngaps = 0
nlist = 0

write(logmsg,*) 'Initial count, max_nlist: ',initial_count, max_nlist
call logger(logmsg)

! How big is cell_list?
!write(*,*) 'size of cell_list: ',sizeof(cell_list(1)),sizeof(cell_list(1))*max_nlist
allocate(cell_list(max_nlist))
allocate(gaplist(max_ngaps))

ok = .true.

end subroutine

!----------------------------------------------------------------------------------------- 
!-----------------------------------------------------------------------------------------
subroutine RngInitialisation
integer, allocatable :: zig_seed(:)
integer :: i
integer :: npar, grainsize = 32

npar = Mnodes
write(logmsg,*) 'npar = ',npar,seed
call logger(logmsg)
allocate(zig_seed(0:npar-1))
do i = 0,npar-1
    zig_seed(i) = seed(1)*seed(2)*(i+1)
enddo
call par_zigset(npar,zig_seed,grainsize)
par_zig_init = .true.
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine ReadCellParams(ok)
logical :: ok
integer :: i, idrug, imetab, nmetab, im, itestcase, Nmm3, ichemo, itreatment, iuse_extra, iuse_relax, iuse_par_relax, iuse_FD
integer :: iuse_oxygen, iuse_glucose, iuse_lactate, iuse_tracer, iuse_drug, iuse_metab, iV_depend, iV_random, iuse_gd_all, iuse_divide_dist
!integer ::  idrug_decay, imetab_decay
integer :: ictype, idisplay, isconstant, ioxygengrowth, iglucosegrowth, ilactategrowth, ioxygendeath, iglucosedeath
integer :: iuse_drop, iconstant, isaveprofiledata, isaveslicedata, iusecellcycle, iusemetabolism, ifullymixed, isynchronise
logical :: use_metabolites
integer :: isaveFACSdata
real(REAL_KIND) :: days, bdry_conc, percent, d_n_limit
real(REAL_KIND) :: sigma(2), DXmm, anoxia_tag_hours, anoxia_death_hours, aglucosia_tag_hours, aglucosia_death_hours
character*(12) :: drug_name
character*(1) :: numstr
type(cycle_parameters_type),pointer :: ccp

ok = .true.
chemo(:)%used = .false.
Vsite_cm3 = 0.2000E-08  ! from spheroid - to retain same scaling 

open(nfcell,file=inputfile,status='old')
read(nfcell,'(a)') header
if (header(1:3) == 'GUI') then
	gui_run_version = header
	header = 'DD/MM/YYYY header_string'
else
	read(nfcell,*) gui_run_version				! program run version number
endif
read(nfcell,*) dll_run_version				! DLL run version number 
!read(nfcell,*) NX							! size of grid
read(nfcell,*) initial_count				! initial number of tumour cells
read(nfcell,*) iuse_divide_dist
read(nfcell,*) divide_time_median(1)
read(nfcell,*) divide_time_shape(1)
read(nfcell,*) divide_time_median(2)
read(nfcell,*) divide_time_shape(2)
read(nfcell,*) iV_depend
read(nfcell,*) iV_random
read(nfcell,*) days							! number of days to simulate
!read(nfcell,*) d_n_limit					! possible limit on diameter or number of cells
!diam_count_limit = d_n_limit
read(nfcell,*) DELTA_T						! time step size (sec)
!read(nfcell,*) NXB							! size of coarse grid = NXB = NYB
!read(nfcell,*) NZB							! size of coarse grid = NZB
!read(nfcell,*) DXF							! fine grid spacing, off-lattice model (um)
!read(nfcell,*) a_separation
!read(nfcell,*) a_force
!read(nfcell,*) c_force
!read(nfcell,*) x0_force
!read(nfcell,*) x1_force
!read(nfcell,*) kdrag
!read(nfcell,*) frandom
read(nfcell,*) NT_CONC						! number of subdivisions of DELTA_T for diffusion computation
!read(nfcell,*) Nmm3							! number of cells/mm^3
!DXmm = 1.0/(Nmm3**(1./3))
!DELTA_X = DXmm/10							! mm -> cm
!Vsite_cm3 = DELTA_X*DELTA_X*DELTA_X			! total site volume (cm^3)
!read(nfcell,*) fluid_fraction				! fraction of the (non-necrotic) tumour that is fluid
read(nfcell,*) Vcell_pL                     ! nominal cell volume in pL
read(nfcell,*) well_area                    ! well bottom area (cm^2)
read(nfcell,*) medium_volume0				! initial total volume (cm^3)
read(nfcell,*) ifullymixed					! medium is fully mixed
fully_mixed = (ifullymixed == 1)
read(nfcell,*) Vdivide0						! nominal cell volume multiple for division
read(nfcell,*) dVdivide						! variation about nominal divide volume
read(nfcell,*) MM_THRESHOLD					! O2 concentration threshold Michaelis-Menten "soft-landing" (uM)
read(nfcell,*) anoxia_threshold			    ! O2 threshold for anoxia (uM)
read(nfcell,*) anoxia_tag_hours				! hypoxic time leading to tagging to die by anoxia (h)
read(nfcell,*) anoxia_death_hours			! time after tagging to death by anoxia (h)
read(nfcell,*) aglucosia_threshold			! O2 threshold for aglucosia (uM)
read(nfcell,*) aglucosia_tag_hours			! hypoxic time leading to tagging to die by aglucosia (h)
read(nfcell,*) aglucosia_death_hours		! time after tagging to death by aglucosia (h)
read(nfcell,*) itestcase                    ! test case to simulate
read(nfcell,*) seed(1)						! seed vector(1) for the RNGs
read(nfcell,*) seed(2)						! seed vector(2) for the RNGs
read(nfcell,*) ncpu_input					! for GUI just a placeholder for ncpu, used only when execute parameter ncpu = 0
read(nfcell,*) Ncelltypes					! maximum number of cell types in the spheroid
do ictype = 1,Ncelltypes
	read(nfcell,*) percent
	celltype_fraction(ictype) = percent/100
!	read(nfcell,*) idisplay
!	celltype_display(ictype) = (idisplay == 1)
enddo
read(nfcell,*) NT_GUI_OUT					! interval between GUI outputs (timesteps)
read(nfcell,*) show_progeny                 ! if != 0, the number of the cell to show descendents of
read(nfcell,*) iuse_oxygen		! chemo(OXYGEN)%used
read(nfcell,*) ioxygengrowth
chemo(OXYGEN)%controls_growth = (ioxygengrowth == 1)
read(nfcell,*) ioxygendeath
chemo(OXYGEN)%controls_death = (ioxygendeath == 1)
read(nfcell,*) chemo(OXYGEN)%diff_coef
read(nfcell,*) chemo(OXYGEN)%medium_diff_coef
read(nfcell,*) chemo(OXYGEN)%membrane_diff_in
!chemo(OXYGEN)%membrane_diff_out = chemo(OXYGEN)%membrane_diff_in
read(nfcell,*) chemo(OXYGEN)%membrane_diff_out
read(nfcell,*) chemo(OXYGEN)%bdry_conc
read(nfcell,*) iconstant
chemo(OXYGEN)%constant = (iconstant == 1)
read(nfcell,*) chemo(OXYGEN)%max_cell_rate
read(nfcell,*) chemo(OXYGEN)%MM_C0
read(nfcell,*) chemo(OXYGEN)%Hill_N
read(nfcell,*) iuse_glucose		!chemo(GLUCOSE)%used
read(nfcell,*) iglucosegrowth
chemo(GLUCOSE)%controls_growth = (iglucosegrowth == 1)
read(nfcell,*) iglucosedeath
chemo(GLUCOSE)%controls_death = (iglucosedeath == 1)
read(nfcell,*) chemo(GLUCOSE)%diff_coef
read(nfcell,*) chemo(GLUCOSE)%medium_diff_coef
read(nfcell,*) chemo(GLUCOSE)%membrane_diff_in
read(nfcell,*) chemo(GLUCOSE)%membrane_diff_out
read(nfcell,*) chemo(GLUCOSE)%bdry_conc
chemo(GLUCOSE)%dose_conc = chemo(GLUCOSE)%bdry_conc
read(nfcell,*) iconstant
chemo(GLUCOSE)%constant = (iconstant == 1)
read(nfcell,*) chemo(GLUCOSE)%max_cell_rate
read(nfcell,*) chemo(GLUCOSE)%MM_C0
read(nfcell,*) chemo(GLUCOSE)%Hill_N

read(nfcell,*) iuse_lactate		
!read(nfcell,*) ilactategrowth
!chemo(LACTATE)%controls_growth = (ilactategrowth == 1)
read(nfcell,*) chemo(LACTATE)%diff_coef
read(nfcell,*) chemo(LACTATE)%medium_diff_coef
read(nfcell,*) chemo(LACTATE)%membrane_diff_in
read(nfcell,*) chemo(LACTATE)%membrane_diff_out
read(nfcell,*) chemo(LACTATE)%bdry_conc
chemo(LACTATE)%bdry_conc = max(0.001,chemo(LACTATE)%bdry_conc)
chemo(LACTATE)%dose_conc = chemo(LACTATE)%bdry_conc
read(nfcell,*) chemo(LACTATE)%max_cell_rate
read(nfcell,*) chemo(LACTATE)%MM_C0
read(nfcell,*) chemo(LACTATE)%Hill_N

read(nfcell,*) iuse_tracer		!chemo(TRACER)%used
read(nfcell,*) chemo(TRACER)%diff_coef
read(nfcell,*) chemo(TRACER)%medium_diff_coef
read(nfcell,*) chemo(TRACER)%membrane_diff_in
read(nfcell,*) chemo(TRACER)%membrane_diff_out
read(nfcell,*) chemo(TRACER)%bdry_conc
read(nfcell,*) iconstant
chemo(TRACER)%constant = (iconstant == 1)
read(nfcell,*) chemo(TRACER)%max_cell_rate
read(nfcell,*) chemo(TRACER)%MM_C0
read(nfcell,*) chemo(TRACER)%Hill_N
! removed old read of TPZ and DNB drug data

read(nfcell,*) LQ(1)%alpha_H
read(nfcell,*) LQ(1)%beta_H
read(nfcell,*) LQ(1)%OER_am
read(nfcell,*) LQ(1)%OER_bm
read(nfcell,*) LQ(1)%K_ms
read(nfcell,*) LQ(1)%death_prob
read(nfcell,*) LQ(1)%growth_delay_factor
read(nfcell,*) LQ(1)%growth_delay_N
read(nfcell,*) LQ(2)%alpha_H
read(nfcell,*) LQ(2)%beta_H
read(nfcell,*) LQ(2)%OER_am
read(nfcell,*) LQ(2)%OER_bm
read(nfcell,*) LQ(2)%K_ms
read(nfcell,*) LQ(2)%death_prob
read(nfcell,*) LQ(2)%growth_delay_factor
read(nfcell,*) LQ(2)%growth_delay_N
read(nfcell,*) iuse_gd_all
use_radiation_growth_delay_all = (iuse_gd_all == 1)
read(nfcell,*) iusecellcycle
use_cell_cycle = (iusecellcycle == 1)
read(nfcell,*) isynchronise
synchronise = (isynchronise == 1)
call ReadCellCycleParameters(nfcell)
read(nfcell,*) iusemetabolism
!use_metabolism = (iusemetabolism == 1)
use_metabolism = .true.
call ReadMetabolismParameters(nfcell)
read(nfcell,*) O2cutoff(1)
read(nfcell,*) O2cutoff(2)
read(nfcell,*) O2cutoff(3)
read(nfcell,*) hypoxia_threshold
read(nfcell,*) growthcutoff(1)
read(nfcell,*) growthcutoff(2)
read(nfcell,*) growthcutoff(3)
read(nfcell,*) Cthreshold
read(nfcell,*) Clabel_threshold
read(nfcell,*) spcrad_value

!read(nfcell,*) iuse_extra
!read(nfcell,*) iuse_relax
!read(nfcell,*) iuse_par_relax
!read(nfcell,*) iuse_FD
!read(nfcell,*) iuse_drop
!read(nfcell,*) Ndrop
!read(nfcell,*) alpha_shape
!read(nfcell,*) beta_shape
!read(nfcell,*) isaveprofiledata
!read(nfcell,*) saveprofile%filebase
!read(nfcell,*) saveprofile%dt
!read(nfcell,*) saveprofile%nt
!read(nfcell,*) isaveslicedata
!read(nfcell,*) saveslice%filebase
!read(nfcell,*) saveslice%dt
!read(nfcell,*) saveslice%nt 

! For the simplified treatment of solving for intra- and extracellular concentration (Caverage) in the monolayer need only a single cell type
!Ncelltypes = 1

read(nfcell,*) Ndrugs_used
if (Ndrugs_used > 0) then
    call ReadDrugData(nfcell)
endif

read(nfcell,*) isaveFACSdata
read(nfcell,*) saveFACS%filebase
read(nfcell,*) saveFACS%tstart
read(nfcell,*) saveFACS%dt
read(nfcell,*) saveFACS%nt

is_radiation = .false.
if (use_events) then
	call ReadProtocol(nfcell)
	use_treatment = .false.
endif
close(nfcell)
call logger('Finished reading input data')

!if (use_PEST) then
!	call ReadPESTParameters
!endif

! Rescale
chemo(OXYGEN)%membrane_diff_in = chemo(OXYGEN)%membrane_diff_in*Vsite_cm3/60		! /min -> /sec
chemo(OXYGEN)%membrane_diff_out = chemo(OXYGEN)%membrane_diff_out*Vsite_cm3/60		! /min -> /sec
chemo(OXYGEN)%max_cell_rate = chemo(OXYGEN)%max_cell_rate*1.0e6						! mol/cell/s -> mumol/cell/s
chemo(GLUCOSE)%membrane_diff_in = chemo(GLUCOSE)%membrane_diff_in*Vsite_cm3/60		! /min -> /sec
chemo(GLUCOSE)%membrane_diff_out = chemo(GLUCOSE)%membrane_diff_out*Vsite_cm3/60	! /min -> /sec
chemo(GLUCOSE)%max_cell_rate = chemo(GLUCOSE)%max_cell_rate*1.0e6					! mol/cell/s -> mumol/cell/s
chemo(LACTATE)%membrane_diff_in = chemo(LACTATE)%membrane_diff_in*Vsite_cm3/60		! /min -> /sec
chemo(LACTATE)%membrane_diff_out = chemo(LACTATE)%membrane_diff_out*Vsite_cm3/60	! /min -> /sec
chemo(LACTATE)%max_cell_rate = chemo(LACTATE)%max_cell_rate*1.0e6					! mol/cell/s -> mumol/cell/s
chemo(TRACER)%membrane_diff_in = chemo(TRACER)%membrane_diff_in*Vsite_cm3/60		! /min -> /sec
chemo(TRACER)%membrane_diff_out = chemo(TRACER)%membrane_diff_out*Vsite_cm3/60		! /min -> /sec


if (celltype_fraction(1) == 1.0) then
	write(nflog,*) 'Type 1 cells'
	selected_celltype = 1
elseif (celltype_fraction(2) == 1.0) then
	write(nflog,*) 'Type 2 cells'
	selected_celltype = 2
else
	write(logmsg,*) 'Error: cells must all be of the same type'
	call logger(logmsg)
	ok = .false.
	return
endif

if (chemo(OXYGEN)%Hill_N /= 1 .and. chemo(OXYGEN)%Hill_N /= 2) then
	call logger('Error: OXYGEN_HILL_N must be 1 or 2')
	ok = .false.
	return
endif
!if (chemo(GLUCOSE)%Hill_N /= 1 .and. chemo(GLUCOSE)%Hill_N /= 2) then
!	call logger('Error: GLUCOSE_HILL_N must be 1 or 2')
!	ok = .false.
!	return
!endif
if (chemo(LACTATE)%Hill_N /= 1 .and. chemo(LACTATE)%Hill_N /= 2) then
	call logger('Error: LACTATE_HILL_N must be 1 or 2')
	ok = .false.
	return
endif
!DXB = 4*DXF
MM_THRESHOLD = MM_THRESHOLD/1000					! uM -> mM
anoxia_threshold = anoxia_threshold/1000			! uM -> mM
aglucosia_threshold = aglucosia_threshold/1000		! uM -> mM
O2cutoff = O2cutoff/1000							! uM -> mM
hypoxia_threshold = hypoxia_threshold/1000			! uM -> mM
!relax = (iuse_relax == 1)
!use_parallel = (iuse_par_relax == 1)
!use_FD = (iuse_FD == 1)
chemo(OXYGEN)%used = (iuse_oxygen == 1)
chemo(GLUCOSE)%used = (iuse_glucose == 1)
chemo(LACTATE)%used = use_metabolism .and. (iuse_lactate == 1)
chemo(TRACER)%used = (iuse_tracer == 1)
chemo(OXYGEN)%MM_C0 = chemo(OXYGEN)%MM_C0/1000		! uM -> mM
chemo(GLUCOSE)%MM_C0 = chemo(GLUCOSE)%MM_C0/1000	! uM -> mM
chemo(LACTATE)%MM_C0 = chemo(LACTATE)%MM_C0/1000	! uM -> mM
if (.not.chemo(OXYGEN)%used) then
    chemo(OXYGEN)%controls_growth = .false.
    chemo(OXYGEN)%controls_death = .false.
endif
if (.not.chemo(GLUCOSE)%used) then
    chemo(GLUCOSE)%controls_growth = .false.
    chemo(GLUCOSE)%controls_death = .false.
endif

!mitosis_duration = ccp%T_M  ! seconds 

LQ(:)%growth_delay_factor = 60*60*LQ(:)%growth_delay_factor	! hours -> seconds
divide_dist(1:2)%class = LOGNORMAL_DIST
divide_time_median(1:2) = 60*60*divide_time_median(1:2)			! hours -> seconds
sigma(1:2) = log(divide_time_shape(1:2))
!divide_dist%p1 = log(divide_time_mean/exp(sigma*sigma/2))	
divide_dist(1:2)%p1 = log(divide_time_median(1:2))	
divide_dist(1:2)%p2 = sigma
divide_time_mean(1:2) = exp(divide_dist(1:2)%p1 + 0.5*divide_dist(1:2)%p2**2)	! mean = median.exp(sigma^2/2)
write(logmsg,'(a,24e12.4)') 'shape, sigma: ',divide_time_shape(1:2),sigma(1:2)
call logger(logmsg)
write(logmsg,'(a,4e12.4)') 'Median, mean divide time: ',divide_time_median(1:2)/3600,divide_time_mean(1:2)/3600
call logger(logmsg)
call AdjustCycleTimes

use_divide_time_distribution = (iuse_divide_dist == 1)
use_V_dependence = (iV_depend == 1)
randomise_initial_volume = (iV_random == 1)
use_constant_divide_volume = (dVdivide == 0)
!use_extracellular_O2 = (iuse_extra == 1)
t_anoxia_limit = 60*60*anoxia_tag_hours				! hours -> seconds
anoxia_death_delay = 60*60*anoxia_death_hours		! hours -> seconds
t_aglucosia_limit = 60*60*aglucosia_tag_hours		! hours -> seconds
aglucosia_death_delay = 60*60*aglucosia_death_hours	! hours -> seconds
Vcell_cm3 = 1.0e-9*Vcell_pL							! nominal cell volume in cm3
Vdivide0 = Vdivide0*Vcell_cm3
total_volume = medium_volume0

write(nflog,*) 'Vdivide0: ',Vdivide0

!write(logmsg,'(a,3e12.4)') 'DELTA_X, cell_radius: ',DELTA_X,cell_radius
!call logger(logmsg)
!write(logmsg,'(a,4e12.4)') 'Volumes: site, extra, cell (average, base): ',Vsite_cm3, Vextra_cm3, Vsite_cm3-Vextra_cm3, Vcell_cm3
!call logger(logmsg)

!saveprofile%active = (isaveprofiledata == 1)
!saveprofile%it = 1
!saveprofile%dt = 60*saveprofile%dt		! mins -> seconds
!saveslice%active = (isaveslicedata == 1)
!saveslice%it = 1
!saveslice%dt = 60*saveslice%dt			! mins -> seconds
saveFACS%active = (isaveFACSdata == 1)
saveFACS%it = 1
saveFACS%dt = 60*saveFACS%dt			! mins -> seconds
saveFACS%tstart = 60*saveFACS%tstart	! mins -> seconds

!use_dropper = (iuse_drop == 1)

! Setup test_case
test_case = .false.
if (itestcase /= 0) then
    test_case(itestcase) = .true.
endif

!if (mod(NX,2) /= 0) NX = NX+1					! ensure that NX is even
!NYB = NXB

open(nfout,file=outputfile,status='replace')
write(nfout,'(a,a)') 'GUI version: ',gui_run_version
write(nfout,'(a,a)') 'DLL version: ',dll_run_version
write(nfout,*)

write(nflog,*)
write(nflog,'(a,a)') 'GUI version: ',gui_run_version
write(nflog,'(a,a)') 'DLL version: ',dll_run_version
write(nflog,*)

if (.not.use_PEST) then
	open(nfres,file='vmonolayer_ts.out',status='replace')
else
	open(nfres,file=PEST_outputfile,status='replace')
endif
!write(nfres,'(a,a)') 'GUI version: ',gui_run_version
!write(nfres,'(a,a)') 'DLL version: ',dll_run_version
!write(nfres,*)
write(nfres,'(a)') 'date info GUI_version DLL_version &
istep hour Ncells(1) Ncells(2) Nviable Nnonviable &
NATP_dead(1) NATP_dead(2) NdrugA_dead(1) NdrugA_dead(2) NdrugB_dead(1) NdrugB_dead(2) &
Nradiation_dead(1) Nradiation_dead(2) Ntotal_dead &
Ntagged_ATP(1) Ntagged_ATP(2) Ntagged_drugA(1) Ntagged_drugA(2) Ntagged_drugB(1) Ntagged_drugB(2) &
Ntagged_radiation(1) Ntagged_radiation(2) &
f_viable f_hypoxic(1) f_hypoxic(2) f_hypoxic(3) f_clonohypoxic(1) f_clonohypoxic(2) f_clonohypoxic(3) f_growth(1) f_growth(2) f_growth(3) &
f_nogrow f_clonogenic plating_efficiency(1) plating_efficiency(2) &
EC_oxygen EC_glucose EC_lactate EC_drugA EC_drugA_metab1 EC_drugA_metab2 EC_drugB EC_drugB_metab1 EC_drugB_metab2 &
IC_oxygen IC_glucose IC_lactate IC_pyruvate IC_drugA IC_drugA_metab1 IC_drugA_metab2 IC_drugB IC_drugB_metab1 IC_drugB_metab2 &
medium_oxygen medium_glucose medium_lactate medium_drugA medium_drugA_metab1 medium_drugA_metab2 medium_drugB medium_drugB_metab1 medium_drugB_metab2 &
G1_phase G1_checkpoint S_phase G2_phase G2_checkpoint M_phase Nmutations &
doubling_time glycolysis_rate pyruvate_oxidation_rate ATP_rate intermediates_rate Ndivided pyruvate_oxidised_fraction'
write(logmsg,*) 'Opened nfout: ',trim(outputfile)
! Note order change
call logger(logmsg)

!if (use_PEST) then
!	open(nfPESTout,file=PEST_outputfile,status='replace')
!	write(logmsg,*) 'Opened PEST outputfile: ',trim(PEST_outputfile)
!	call logger(logmsg)
!endif

Nsteps = days*24*60*60/DELTA_T		! DELTA_T in seconds
NT_DISPLAY = 2						! This is the updating interval (calls to get_summary) in the GUI version.  Not used by command-line version.
DT_DISPLAY = NT_DISPLAY*DELTA_T
write(logmsg,'(a,2i6,f6.0)') 'nsteps, NT_CONC, DELTA_T: ',nsteps,NT_CONC,DELTA_T
call logger(logmsg)
write(logmsg,'(a,i6,f6.0)') 'NT_DISPLAY, DT_DISPLAY: ',NT_DISPLAY, DT_DISPLAY
call logger(logmsg)

if (.not.use_new_drugdata) then
	call DetermineKd	! Kd is now set or computed in the GUI 
endif
ok = .true.

end subroutine

!-----------------------------------------------------------------------------------------
! The cell cycle parameters include the parameters for radiation damage and repair, 
! and for the associated checkpoint duration limits Tcp(:).
! Time unit = hour
!-----------------------------------------------------------------------------------------
subroutine ReadCellCycleParameters(nf)
integer :: nf
type(cycle_parameters_type),pointer :: ccp
integer :: ityp
real(REAL_KIND) :: total

do ityp = 1,2
ccp => cc_parameters(ityp)

read(nf,*) ccp%T_G1
read(nf,*) ccp%T_S
read(nf,*) ccp%T_G2
read(nf,*) ccp%T_M
read(nf,*) ccp%G1_mean_delay
read(nf,*) ccp%G2_mean_delay
read(nf,*) ccp%Apoptosis_rate
read(nf,*) ccp%eta_PL
read(nf,*) ccp%eta_IRL
read(nf,*) ccp%Krepair_base
read(nf,*) ccp%Krepair_max
read(nf,*) ccp%Kmisrepair
read(nf,*) ccp%mitosis_factor
read(nf,*) ccp%fraction_Ch1
read(nf,*) ccp%psurvive_Ch1
read(nf,*) ccp%psurvive_Ch2
read(nf,*) ccp%aTCP
read(nf,*) ccp%bTCP
!read(nf,*) ccp%Kcp

total = ccp%T_G1 + ccp%T_S + ccp%T_G2 + ccp%T_M + ccp%G1_mean_delay + ccp%G2_mean_delay
write(nflog,'(a,7f8.2)') 'T_G1,T_S,T_G2,T_M,G1_delay,G2_delay, total: ',ccp%T_G1,ccp%T_S,ccp%T_G2,ccp%T_M, &
						ccp%G1_mean_delay,ccp%G2_mean_delay,total
						
ccp%T_G1 = 3600*ccp%T_G1                    ! hours -> seconds
ccp%T_S = 3600*ccp%T_S
ccp%T_G2 = 3600*ccp%T_G2
ccp%T_M = 3600*ccp%T_M
ccp%G1_mean_delay = 3600*ccp%G1_mean_delay
ccp%G2_mean_delay = 3600*ccp%G2_mean_delay
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine AdjustCycleTimes
integer :: ityp
real(REAL_KIND) :: tmean, tsum, tfactor
type(cycle_parameters_type),pointer :: ccp

do ityp = 1,2
	ccp => cc_parameters(ityp)
	tmean = divide_time_mean(ityp)
	tsum = ccp%T_G1 + ccp%T_S + ccp%T_G2 + ccp%T_M + ccp%G1_mean_delay + ccp%G2_mean_delay
	tfactor = tmean/tsum
	ccp%T_G1 = tfactor*ccp%T_G1
	ccp%T_S = tfactor*ccp%T_S
	ccp%T_G2 = tfactor*ccp%T_G2
	ccp%T_M = tfactor*ccp%T_M
	ccp%G1_mean_delay = tfactor*ccp%G1_mean_delay
	ccp%G2_mean_delay= tfactor*ccp%G2_mean_delay
	ccp%Pk_G1 = 1./ccp%G1_mean_delay    ! /sec
	ccp%Pk_G2 = 1./ccp%G2_mean_delay    ! /sec
	write(nflog,'(a,4e12.3)') 'Pk_G1, Pk_G2: ',ccp%Pk_G1,ccp%Pk_G2
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ReadMetabolismParameters(nf)
integer :: nf
integer :: ityp

read(nf,*) f_G_norm
read(nf,*) f_P_norm
read(nf,*) N_GA
read(nf,*) N_PA
read(nf,*) N_GI
read(nf,*) N_PI
read(nf,*) N_PO
read(nf,*) K_H1
read(nf,*) K_H2
read(nf,*) K_HB
read(nf,*) K_PDK
read(nf,*) PDKmin
read(nf,*) C_O2_norm
read(nf,*) C_G_norm
read(nf,*) C_L_norm
read(nf,*) f_ATPs
read(nf,*) f_ATPg
read(nf,*) f_ATPramp
read(nf,*) K_PL
read(nf,*) K_LP
read(nf,*) Hill_Km_P
Hill_N_P = 1
Hill_Km_P = Hill_Km_P/1000		! uM -> mM
!ATP_Km = ATP_Km/1000			! uM -> mM
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ReadDrugData(nf)
integer :: nf
integer :: idrug, im, ictyp, ival
character*(16) :: drugname

write(logmsg,*) 'ReadDrugData'
call logger(logmsg)
if (allocated(drug)) then
	deallocate(drug)
endif
allocate(drug(Ndrugs_used))
do idrug = 1,Ndrugs_used
	read(nf,'(a)') drug(idrug)%classname
	if (drug(idrug)%classname == 'TPZ') then
		drug(idrug)%drugclass = TPZ_CLASS
	elseif (drug(idrug)%classname == 'DNB') then
		drug(idrug)%drugclass = DNB_CLASS
	endif
	drug(idrug)%nmetabolites = 2			! currently all drugs have 2 metabolites
	drug(idrug)%use_metabolites = .true.	! currently simulate metabolites
	drug(idrug)%phase_dependent = .false.
	drug(idrug)%active_phase = .false.
    do im = 0,2			! 0 = parent, 1 = metab_1, 2 = metab_2
		read(nf,'(a)') drugname
		if (im == 0) then
			drug(idrug)%name = drugname
		endif
		read(nf,*) drug(idrug)%diff_coef(im)
		read(nf,*) drug(idrug)%medium_diff_coef(im)
		read(nf,*) drug(idrug)%membrane_diff_in(im)
		read(nf,*) drug(idrug)%membrane_diff_out(im)
		read(nf,*) drug(idrug)%halflife(im)
		drug(idrug)%membrane_diff_in(im) = drug(idrug)%membrane_diff_in(im)*Vsite_cm3/60	! /min -> /sec
		drug(idrug)%membrane_diff_out(im) = drug(idrug)%membrane_diff_out(im)*Vsite_cm3/60	! /min -> /sec
		do ictyp = 1,ncelltypes
            read(nf,*) drug(idrug)%Kmet0(ictyp,im)
            read(nf,*) drug(idrug)%C2(ictyp,im)
            read(nf,*) drug(idrug)%KO2(ictyp,im)
            read(nf,*) drug(idrug)%Vmax(ictyp,im)
            read(nf,*) drug(idrug)%Km(ictyp,im)
            read(nf,*) drug(idrug)%Klesion(ictyp,im)
            read(nf,*) drug(idrug)%kill_O2(ictyp,im)
            read(nf,*) drug(idrug)%kill_drug(ictyp,im)
            read(nf,*) drug(idrug)%kill_duration(ictyp,im)
            read(nf,*) drug(idrug)%kill_fraction(ictyp,im)
            read(nf,*) drug(idrug)%SER_max(ictyp,im)
            read(nf,*) drug(idrug)%SER_Km(ictyp,im)
            read(nf,*) drug(idrug)%SER_KO2(ictyp,im)
            read(nf,*) drug(idrug)%n_O2(ictyp,im)
            read(nf,*) drug(idrug)%death_prob(ictyp,im)
            if (use_new_drugdata) then
	            read(nf,*) drug(idrug)%Kd(ictyp,im)
	        endif
            read(nf,*) ival
            drug(idrug)%kills(ictyp,im) = (ival == 1)
            read(nf,*) ival
            drug(idrug)%kill_model(ictyp,im) = ival
            read(nf,*) ival
            drug(idrug)%sensitises(ictyp,im) = (ival == 1)
            drug(idrug)%Vmax(ictyp,im) = drug(idrug)%Vmax(ictyp,im)/60						! /min -> /sec
            drug(idrug)%Kmet0(ictyp,im) = drug(idrug)%Kmet0(ictyp,im)/60					! /min -> /sec
            drug(idrug)%KO2(ictyp,im) = 1.0e-3*drug(idrug)%KO2(ictyp,im)					! um -> mM
            drug(idrug)%kill_duration(ictyp,im) = 60*drug(idrug)%kill_duration(ictyp,im)	! min -> sec
		enddo
	    if (drug(idrug)%name == 'EDU') then
			drug(idrug)%nmetabolites = 1
			drug(idrug)%phase_dependent = .true.
			drug(idrug)%active_phase(S_phase) = .true.
		endif
	    if (drug(idrug)%name == 'PI') then
			drug(idrug)%nmetabolites = 1
			drug(idrug)%phase_dependent = .true.
			drug(idrug)%active_phase(1:6) = .true.
		endif
    enddo
    write(nflog,*) 'drug: ',idrug,drug(idrug)%classname,'  ',drug(idrug)%name
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Skip lines until the 'PROTOCOL' line
!-----------------------------------------------------------------------------------------
subroutine ReadProtocol(nf)
integer :: nf
integer :: itime, ntimes, kevent, ichemo, idrug, im
character*(64) :: line
character*(16) :: drugname
character*(1)  :: numstr
character*(1) :: fullstr
real(REAL_KIND) :: t, dt, vol, conc, O2conc, O2flush, dose, O2medium, glumedium
type(event_type) :: E

write(logmsg,*) 'ReadProtocol'
call logger(logmsg)
chemo(TRACER+1:)%used = .false.
do
	read(nf,'(a)') line
	if (trim(line) == 'PROTOCOL') exit
enddo
read(nf,*) ntimes
if (ntimes == 0) then
	call logger('no events')
	Nevents = 0
	return
endif
Nevents = ntimes
if (allocated(event)) deallocate(event)
allocate(event(2*ntimes))
kevent = 0
do itime = 1,ntimes
	read(nf,'(a)') line
	write(nflog,'(a)') line
	if (trim(line) == 'DRUG') then
		kevent = kevent + 1
		event(kevent)%etype = DRUG_EVENT
		read(nf,'(a)') line
		write(nflog,'(a)') line
		drugname = trim(line)
		write(nflog,*) 'ndrugs_used: ',ndrugs_used
		do idrug = 1,ndrugs_used
			if (drugname == drug(idrug)%name) then
				ichemo = TRACER + 1 + 3*(idrug-1)
				exit
			endif
		enddo
		write(nflog,*) 'ichemo: ',ichemo
		! Need to copy drug(idrug) parameters to chemo(ichemo) 
		call CopyDrugParameters(idrug,ichemo)
		read(nf,*) t
		read(nf,*) dt
		read(nf,*) vol
		read(nf,*) O2conc
		read(nf,*) O2flush
		read(nf,*) conc
		event(kevent)%time = t
		event(kevent)%ichemo = ichemo
		event(kevent)%idrug = idrug
		event(kevent)%volume = vol
		event(kevent)%conc = conc
		event(kevent)%O2conc = O2conc
		event(kevent)%dose = 0
		event(kevent)%full = .false.	
		chemo(ichemo)%used = .true.
		write(nflog,'(a,i3,2f8.3)') 'define DRUG_EVENT: volume, O2conc: ',kevent,event(kevent)%volume,event(kevent)%O2conc
		if (drug(idrug)%use_metabolites) then
			do im = 1,drug(idrug)%nmetabolites
				chemo(ichemo+im)%used = .true.
			enddo
		endif

		kevent = kevent + 1
		event(kevent)%etype = MEDIUM_EVENT
		event(kevent)%time = t + dt
		event(kevent)%ichemo = 0
!		event(kevent)%volume = medium_volume0
		event(kevent)%volume = total_volume
		event(kevent)%conc = 0
		event(kevent)%O2medium = O2flush		
		event(kevent)%glumedium = chemo(GLUCOSE)%dose_conc		
		event(kevent)%lacmedium = chemo(LACTATE)%dose_conc	
		event(kevent)%full = .false.	
		event(kevent)%dose = 0
		write(nflog,'(a,i3,2f8.3)') 'define MEDIUM_EVENT: volume: ',kevent,event(kevent)%volume,event(kevent)%O2medium
	elseif (trim(line) == 'MEDIUM') then
		kevent = kevent + 1
		event(kevent)%etype = MEDIUM_EVENT
		read(nf,*) t
		read(nf,*) vol
		read(nf,*) fullstr
		read(nf,*) O2medium
		read(nf,*) glumedium
		event(kevent)%time = t
		event(kevent)%volume = vol	
		event(kevent)%ichemo = 0
		event(kevent)%O2medium = O2medium
		event(kevent)%glumedium = glumedium
		event(kevent)%lacmedium = chemo(LACTATE)%dose_conc
		event(kevent)%full = (trim(fullstr) == 'Y' .or. trim(fullstr) == 'y')
		event(kevent)%dose = 0
		write(nflog,'(a,i3,2f8.3)') 'define MEDIUM_EVENT: volume: ',kevent,event(kevent)%volume,event(kevent)%O2medium
	elseif (trim(line) == 'RADIATION') then
        is_radiation = .true.
		kevent = kevent + 1
		event(kevent)%etype = RADIATION_EVENT
		read(nf,*) t
		read(nf,*) dose
		event(kevent)%time = t
		event(kevent)%dose = dose	
		event(kevent)%ichemo = 0
		event(kevent)%volume = 0
		event(kevent)%conc = 0
	endif
enddo
Nevents = kevent
! Set events not done
! convert time from hours to seconds
! convert volume from uL to cm^3  NO LONGER - now using cm^3 everywhere
write(logmsg,*) 'nevents: ',nevents
call logger(logmsg)
do kevent = 1,Nevents
	event(kevent)%done = .false.
	event(kevent)%time = event(kevent)%time*60*60
!	event(kevent)%volume = event(kevent)%volume*1.0e-3
	E = event(kevent)
!	write(*,'(a,i3,f8.0,2i3,3f8.4)') 'event: ',kevent,E%time,E%etype,E%ichemo,E%volume,E%conc,E%dose
enddo
! Check that events are sequential
do kevent = 1,Nevents-1
	if (event(kevent)%time >= event(kevent+1)%time) then
		write(logmsg,*) 'Error: non-sequential event: ',kevent,event(kevent)%time
		call logger(logmsg)
		stop
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine CopyDrugParameters(idrug,ichemo)
integer :: idrug,ichemo
integer :: im, im1, im2
character*(1) :: numstr

im1 = 0
chemo(ichemo)%name = drug(idrug)%name
if (drug(idrug)%use_metabolites) then
	do im = 1,drug(idrug)%nmetabolites
		chemo(ichemo+im)%used = .true.
		chemo(ichemo+im)%name = trim(chemo(ichemo)%name) // '_metab'
		write(numstr,'(i1)') im
		chemo(ichemo+im)%name = trim(chemo(ichemo+im)%name) // numstr
	enddo
	im2 = 2
else
	im2 = 0
endif
do im = im1, im2
	chemo(ichemo+im)%diff_coef = drug(idrug)%diff_coef(im)
	chemo(ichemo+im)%medium_diff_coef = drug(idrug)%medium_diff_coef(im)
	chemo(ichemo+im)%membrane_diff_in = drug(idrug)%membrane_diff_in(im)
	chemo(ichemo+im)%membrane_diff_out = drug(idrug)%membrane_diff_out(im)
	chemo(ichemo+im)%halflife = drug(idrug)%halflife(im)
!	chemo(ichemo+im)%medium_dlayer = d_layer
	chemo(ichemo+im)%decay = (chemo(ichemo+im)%halflife > 0)
	if (chemo(ichemo+im)%decay) then
		chemo(ichemo+im)%decay_rate = DecayRate(chemo(ichemo+im)%halflife)
	else
		chemo(ichemo+im)%decay_rate = 0
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! SN30000 CT model
! ----------------
! The rate of cell killing, which becomes a cell death probability rate, is inferred from
! the cell kill experiment.
! The basic assumption is that the rate of killing depends on the drug metabolism rate.
! There are five models:
! kill_model = 1:
!   killing rate = c = Kd.dM/dt
! kill_model = 2:
!   killing rate = c = Kd.Ci.dM/dt
! kill_model = 3:
!   killing rate = c = Kd.(dM/dt)^2
! kill_model = 4:
!   killing rate = c = Kd.Ci
! kill_model = 5:
!   killing rate = c = Kd.Ci^2
! where dM/dt = F(O2).kmet0.Ci
! In the kill experiment both O2 and Ci are held constant:
! O2 = CkillO2, Ci = Ckill
! In this case c is constant and the cell population N(t) is given by:
! N(t) = N(0).exp(-ct), i.e.
! c = -log(N(T)/N(0))/T where T is the duration of the experiment
! N(T)/N(0) = 1 - f, where f = kill fraction
! kill_model = 1:
!   c = Kd.F(CkillO2).kmet0.Ckill => Kd = -log(1-f)/(T.F(CkillO2).kmet0.Ckill)
! kill_model = 2:
!   c = Kd.F(CkillO2).kmet0.Ckill^2 => Kd = -log(1-f)/(T.F(CkillO2).kmet0.Ckill^2)
! kill_model = 3:
!   c = Kd.(F(CkillO2).kmet0.Ckill)^2 => Kd = -log(1-f)/(T.(F(CkillO2).kmet0.Ckill)^2)
! kill_model = 4:
!   c = Kd.Ckill => Kd = -log(1-f)/(T.Ckill)
! kill_model = 5:
!   c = Kd.Ckill^2 => Kd = -log(1-f)/(T.Ckill^2)
!-----------------------------------------------------------------------------------------
subroutine DetermineKd
real(REAL_KIND) :: C2, KO2, n_O2, Kmet0, kmet 
real(REAL_KIND) :: f, T, Ckill, Ckill_O2, Kd
integer :: idrug, ictyp, im, kill_model

do idrug = 1,ndrugs_used
!	if (idrug == 1 .and. .not.chemo(TPZ_DRUG)%used) cycle
!	if (idrug == 2 .and. .not.chemo(DNB_DRUG)%used) cycle
	do ictyp = 1,Ncelltypes
		do im = 0,2
			if (drug(idrug)%kills(ictyp,im)) then
				C2 = drug(idrug)%C2(ictyp,im)
				KO2 = drug(idrug)%KO2(ictyp,im)
				n_O2 = drug(idrug)%n_O2(ictyp,im)
				Kmet0 = drug(idrug)%Kmet0(ictyp,im)
				kill_model = drug(idrug)%kill_model(ictyp,im)
				Ckill_O2 = drug(idrug)%kill_O2(ictyp,im)
				f = drug(idrug)%kill_fraction(ictyp,im)
				T = drug(idrug)%kill_duration(ictyp,im)
				Ckill = drug(idrug)%kill_drug(ictyp,im)
				kmet = (1 - C2 + C2*(KO2**n_O2)/(KO2**n_O2 + Ckill_O2**n_O2))*Kmet0
				if (kill_model == 1) then
					Kd = -log(1-f)/(T*kmet*Ckill)
				elseif (kill_model == 2) then
					Kd = -log(1-f)/(T*kmet*Ckill**2)
				elseif (kill_model == 3) then
					Kd = -log(1-f)/(T*(kmet*Ckill)**2)
				elseif (kill_model == 4) then
					Kd = -log(1-f)/(T*Ckill)
				elseif (kill_model == 5) then
					Kd = -log(1-f)/(T*Ckill**2)
				endif
				drug(idrug)%Kd(ictyp,im) = Kd
			endif
!			if (idrug == 1) then
!				TPZ%Kd(i) = Kd
!			elseif (idrug == 2) then
!				DNB%Kd(i,im) = Kd
!			endif
		enddo
	enddo
enddo

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine PlaceCells(ok)
logical :: ok
integer :: kcell, k, ichemo, ityp, site(3)
real(REAL_KIND) :: rsite(3)

lastID = 0
kcell = 0
Ncells_type = 0
!if (kcell > initial_count) then
!	write(logmsg,*) 'Cell count already exceeds specified number: ',kcell,initial_count
!	call logger(logmsg)
!	ok = .false.
!	return
!endif
!! Now add cells to make the count up to the specified initial_count
!if (kcell < initial_count) then
!	call AddBdryCells(kcell)
!	kcell = initial_count
!endif

rsite = [0.,0.,0.]
do kcell = 1,initial_count
	call AddCell(kcell,rsite)
enddo
nlist = kcell-1
Ncells = nlist
Ncells0 = Ncells
Nviable = Ncells_type
!Nreuse = 0	
ok = .true.
!write(logmsg,*) 'idbug: ',idbug
!call logger(logmsg)
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine AddCell(kcell, rsite)
integer :: kcell
real(REAL_KIND) :: rsite(3)
integer :: ityp, k, kpar = 0
real(REAL_KIND) :: v(3), c(3), R1, R2, V0, Tdiv, Vdiv, p(3), R, gfactor
type(cell_type), pointer :: cp
type(cycle_parameters_type),pointer :: ccp
type(metabolism_type), pointer :: metabolic
	
metabolic => phase_metabolic(1)
cp => cell_list(kcell)
cp%ID = kcell
cp%state = ALIVE
cp%generation = 1
cp%celltype = random_choice(celltype_fraction,Ncelltypes,kpar)
ityp = cp%celltype
ccp => cc_parameters(ityp)
Ncells_type(ityp) = Ncells_type(ityp) + 1
cp%Iphase = .true.
!cp%nspheres = 1

V0 = Vdivide0/2
cp%divide_volume = get_divide_volume(ityp, V0, Tdiv, gfactor)
cp%divide_time = Tdiv
cp%fg = gfactor
cp%dVdt = max_growthrate(ityp)
cp%metab%I_rate = metabolic%I_rate_max	! this is just to ensure that initial growth rate is not 0
if (use_volume_method) then
    !cp%divide_volume = Vdivide0
    if (initial_count == 1) then
	    cp%V = 0.9*cp%divide_volume
    else
    !	cp%V = (0.5 + 0.49*par_uni(kpar))*cp%divide_volume
        if (randomise_initial_volume) then
	        cp%V = cp%divide_volume*0.5*(1 + par_uni(kpar))
        else
	        cp%V = cp%divide_volume/1.6
        endif
    endif
    cp%t_divide_last = 0    ! not correct
else	! use cell cycle
    cp%NL1 = 0
    cp%NL2 = 0
    cp%N_PL = 0
    cp%N_Ch1 = 0
    cp%N_Ch2 = 0
    cp%irrepairable = .false.
    ! Need to assign phase, volume to complete phase, current volume
    call SetInitialCellCycleStatus(cp)
    cp%starved = .false.
endif
if (use_metabolism) then	! Fraction of I needed to divide = fraction of volume needed to divide
!	cp%metab%I2Divide = get_I2Divide(cp)
!	cp%metab%Itotal = cp%metab%I2Divide*(cp%V - V0)/(cp%divide_volume - V0)
	cp%dVdt = max_growthrate(ityp)
endif
!cp%radius(1) = (3*cp%V/(4*PI))**(1./3.)
!cp%centre(:,1) = rsite 
!cp%site = rsite/DELTA_X + 1
!cp%d = 0
!cp%birthtime = 0
!cp%growthrate = test_growthrate
!cp2%divide_volume = get_divide_volume()
!cp%d_divide = (3*cp%divide_volume/PI)**(1./3.)
!cp%mitosis = 0

cp%ATP_tag = .false.
cp%drug_tag = .false.
cp%radiation_tag = .false.
!cp%anoxia_tag = .false.
!cp%aglucosia_tag = .false.
cp%growth_delay = .false.
cp%G2_M = .false.
cp%p_rad_death = 0

!cp%t_anoxia = 0
!cp%t_aglucosia = 0

call get_random_vector3(v)	! set initial axis direction
!cp%d = 0.1*small_d
!c = cp%centre(:,1)
!cp%centre(:,1) = c + (cp%d/2)*v
!cp%centre(:,2) = c - (cp%d/2)*v
!cp%nbrs = 0
!cp%Cex = Caverage(1,1,1,:)	! initially the concentrations are the same everywhere
cp%Cin(OXYGEN) = chemo(OXYGEN)%bdry_conc
cp%Cin(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
cp%Cin(LACTATE) = chemo(LACTATE)%bdry_conc
cp%CFSE = generate_CFSE(1.d0)

!cp%growth_rate_factor = get_growth_rate_factor()
cp%ATP_rate_factor = get_ATP_rate_factor()
!cp%ndt = ndt
end subroutine

!--------------------------------------------------------------------------------
! %divide_time and %fg have been generated
! Assuming growth rate is max_growthrate
! Assuming NOT volume_based_transition, need to determine: 
! for G1:
!	%G1_time
! for Checkpoint1:
!	%G1_flag
!	%G1S_time
! for S:
!	%S_time
! for G2:
!	%G2_time
! for Checkpoint2:
!	%G2_flag
!	%G2M_time
! for M:
!	%M_time
! for all phases:
!	%V
!	%t_divide_last
! Note: no cells start in mitosis - set all phase=6 cells at the end of G2 checkpoint 
!--------------------------------------------------------------------------------
subroutine SetInitialCellCycleStatus(cp)
type(cell_type), pointer :: cp
type(cycle_parameters_type),pointer :: ccp
integer :: ityp, iphase
integer :: kpar = 0
real(REAL_KIND) :: Tdiv, Tmean, V0, fg, rVmax, fsum, R, x, y, z, phase_fraction(6), phase_time(6)

ityp = cp%celltype
ccp => cc_parameters(ityp)
Tdiv = cp%divide_time
V0 = Vdivide0/2
Tmean = divide_time_mean(ityp)
rVmax = max_growthrate(ityp)
fg = cp%fg
R = par_uni(kpar)
x = (4 - sqrt(16 - 12*R))/2	
! This is the level of progress through the cell cycle,  0 -> 1, for prob. density f(x) = 1 - 2(x-0.5)/3
! Cumulative prob. function F(x) = -x^2/3 + 4x/3, then from R U(0,1): (-x^2 + 4x)/3 = R, x^2 - 4x + 3R = 0
phase_time(1) = fg*ccp%T_G1
phase_time(2) = ccp%G1_mean_delay
phase_time(3) = fg*ccp%T_S
phase_time(4) = fg*ccp%T_G2
phase_time(5) = ccp%G2_mean_delay
phase_time(6) = ccp%T_M
phase_fraction = phase_time/Tdiv
! These fractions must sum to 1 because of get_divide_volume (check)

if (synchronise) then	! all cells are starting M phase
	cp%V = V0 + (phase_time(1) + phase_time(3) + phase_time(4))*rVmax 
	cp%phase = Checkpoint2
	cp%G2M_time = 0
	cp%G2_flag = .true.
	cp%t_divide_last = -(phase_time(1) + phase_time(2) + phase_time(3) + phase_time(4) + phase_time(5))
	return
endif
fsum = 0
do iphase = 1,6
	if (fsum + phase_fraction(iphase) > x) then	! this is the phase
		y = (x - fsum)/phase_fraction(iphase)	! this is fractional progress through the phase
		z = 1 - y								! this is fraction phase left to complete, => time until phase transition
		if (iphase == G1_phase) then
			cp%phase = G1_phase
			cp%G1_time = z*phase_time(1)
			cp%V = V0 + y*phase_time(1)*rVmax
			cp%t_divide_last = -y*phase_time(1)
		elseif (iphase == Checkpoint1) then
			cp%phase = Checkpoint1
			cp%G1S_time = z*phase_time(2)
			cp%G1_flag = .false.
			cp%V = V0 + phase_time(1)*rVmax
			cp%t_divide_last = -(phase_time(1) + y*phase_time(2))
		elseif (iphase == S_phase) then
			cp%phase = S_phase
			cp%S_start_time = -y*phase_time(3)
			cp%S_time = z*phase_time(3)
			cp%V = V0 + (phase_time(1) + y*phase_time(3))*rVmax 
			cp%t_divide_last = -(phase_time(1) + phase_time(2) + y*phase_time(3))
		elseif (iphase == G2_phase) then
			cp%phase = G2_phase
			cp%G2_time = z*phase_time(4)
			cp%V = V0 + (phase_time(1) + phase_time(3) + y*phase_time(4))*rVmax 
			cp%t_divide_last = -(phase_time(1) + phase_time(2) + phase_time(3) + y*phase_time(4))
		elseif (iphase == Checkpoint2) then
			cp%phase = Checkpoint2
			cp%G2M_time = z*phase_time(5)
			cp%G2_flag = .false.
			cp%V = V0 + (phase_time(1) + phase_time(3) + phase_time(4))*rVmax 
			cp%t_divide_last = -(phase_time(1) + phase_time(2) + phase_time(3) + phase_time(4) + y*phase_time(5))
		elseif (iphase == M_phase) then
!			cp%phase = M_phase
!			cp%M_time = z*phase_time(6)
			cp%V = V0 + (phase_time(1) + phase_time(3) + phase_time(4))*rVmax 
!			cp%t_divide_last = -(phase_time(1) + phase_time(2) + phase_time(3) + phase_time(4) + phase_time(5) + y*phase_time(6))
			cp%phase = Checkpoint2
			cp%G2M_time = 0
			cp%G2_flag = .false.
			cp%t_divide_last = -(phase_time(1) + phase_time(2) + phase_time(3) + phase_time(4) + phase_time(5))
		else
			write(nfout,*) 'Error in SetInitialCellCycleStatus' 
			stop
		endif
		exit
	endif
	fsum = fsum + phase_fraction(iphase)
enddo
!write(*,*)
!write(*,'(a,3f8.3)') 'Tdiv, Tmean, fg: ',Tdiv/3600,Tmean/3600,fg
!write(*,'(a,6f8.3)') 'phase_time: ',phase_time/3600
!write(*,'(a,6f8.3)') 'phase_fraction: ',phase_fraction
!write(*,'(a,i2,5f8.3,e12.3)') 'iphase, R,x,y,z,tlast,V: ',iphase,R,x,y,z,cp%t_divide_last/3600,cp%V
!if (iphase >= 5) write(*,'(a,2e12.3)') 'Vdiv, V: ',cp%divide_volume, cp%V
!write(*,'(a,2f8.3)') 'Tdiv,sum of phases: ',Tdiv/3600,(phase_time(1) + phase_time(2) + phase_time(3) + phase_time(4) + phase_time(5) + phase_time(6))/3600
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine setTestCell(kcell)
integer :: kcell
integer :: ityp
real(REAL_KIND) :: V0, Tdiv, Tgrowth, Tgrowth0, Tfixed, rVmax, phase_time1
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp

ityp = 1
ccp => cc_parameters(ityp)
cp => cell_list(kcell)
V0 = Vdivide0/2
Tdiv = divide_time_mean(ityp)
cp%divide_time = Tdiv
rVmax = max_growthrate(ityp)
write(nflog,'(a,3f8.0)') 'ccp%T_G1, ccp%T_S, ccp%T_G2: ',ccp%T_G1, ccp%T_S, ccp%T_G2
Tgrowth0 = ccp%T_G1 + ccp%T_S + ccp%T_G2
Tfixed = ccp%T_M + ccp%G1_mean_delay + ccp%G2_mean_delay
write(nflog,'(a,3f8.0)') 'ccp%T_M, ccp%G1_mean_delay, ccp%G2_mean_delay: ',ccp%T_M, ccp%G1_mean_delay, ccp%G2_mean_delay
Tgrowth = Tdiv - Tfixed
write(nflog,'(a,4f8.0)') 'Tdiv, Tfixed, Tgrowth0, Tgrowth: ',Tdiv, Tfixed, Tgrowth0, Tgrowth
cp%fg = Tgrowth/Tgrowth0
cp%divide_volume = V0 + Tgrowth*rVmax
phase_time1 = cp%fg*ccp%T_G1
cp%phase = G1_phase
cp%G1_time = phase_time1
cp%V = V0
cp%t_divide_last = 0
write(nflog,*) 'setTestCell'
write(nflog,*) 'phase: ',cp%phase
write(nflog,*) 'divide_time:', cp%divide_time
write(nflog,*) 'V: ',cp%V
write(nflog,*) 'divide_volume: ',cp%divide_volume
write(nflog,*) 'fg: ',cp%fg
write(nflog,*) 'G1_time: ',cp%G1_time
end subroutine

!-------------------------------------------------------------------------------- 
!--------------------------------------------------------------------------------
subroutine oldAddCell(k,site)
integer :: k, site(3)
integer :: ityp, kpar = 0
real(REAL_KIND) :: V0, Tdiv, R, gfactor

lastID = lastID + 1
cell_list(k)%ID = lastID
cell_list(k)%celltype = random_choice(celltype_fraction,Ncelltypes,kpar)
ityp = cell_list(k)%celltype
Ncells_type(ityp) = Ncells_type(ityp) + 1
!cell_list(k)%site = site
cell_list(k)%state = 1
cell_list(k)%generation = 1
!cell_list(k)%drugA_tag = .false.
!cell_list(k)%drugB_tag = .false.
cell_list(k)%drug_tag = .false.
cell_list(k)%radiation_tag = .false.
cell_list(k)%ATP_tag = .false.
!cell_list(k)%anoxia_tag = .false.
!cell_list(k)%aglucosia_tag = .false.
!cell_list(k)%exists = .true.
cell_list(k)%active = .true.
cell_list(k)%growth_delay = .false.
cell_list(k)%G2_M = .false.
cell_list(k)%p_rad_death = 0
!R = par_uni(kpar)
!cell_list(k)%divide_volume = Vdivide0 + dVdivide*(2*R-1)
V0 = Vdivide0/2
cell_list(k)%divide_volume = get_divide_volume(ityp,V0, Tdiv, gfactor)
cell_list(k)%divide_time = Tdiv
cell_list(k)%fg = gfactor
R = par_uni(kpar)
if (randomise_initial_volume) then
	cell_list(k)%V = cell_list(k)%divide_volume*0.5*(1 + R)
else
	cell_list(k)%V = 1.0*Vcell_cm3
endif
!write(nflog,'(a,i6,f6.2)') 'volume: ',k,cell_list(k)%V
cell_list(k)%t_divide_last = 0		! used in colony growth
!cell_list(k)%t_anoxia = 0
cell_list(k)%Cin = 0
cell_list(k)%Cin(OXYGEN) = chemo(OXYGEN)%bdry_conc
cell_list(k)%Cin(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
cell_list(k)%Cin(TRACER) = chemo(TRACER)%bdry_conc
cell_list(k)%CFSE = generate_CFSE(1.d0)
cell_list(k)%M = 0
!occupancy(site(1),site(2),site(3))%indx(1) = k
end subroutine

!--------------------------------------------------------------------------------
! Add cells at the boundary to bring the total count from k up to initial_count
! (1) Make a list of all boundary sites (sites in contact with an OUTSIDE site)
! (2) Iteratively traverse the list to select the adjacent OUTSIDE site closest 
! to the centre.
!--------------------------------------------------------------------------------
!subroutine AddBdryCells(klast)
!integer :: klast
!integer :: kcell, i, kb, site(3), nbsite(3), nbt, kbmin, imin
!integer, allocatable :: sitelist(:,:)
!real(REAL_KIND) :: r2, r2min
!
!nbt = 0
!do kcell = 1,klast
!	site = cell_list(kcell)%site
!	do i = 1,27
!		if (i == 14) cycle
!		nbsite = site + jumpvec(:,i)
!		if (occupancy(nbsite(1),nbsite(2),nbsite(3))%indx(1) == OUTSIDE_TAG) then
!			nbt = nbt+1
!			exit
!		endif
!	enddo
!enddo
!
!allocate(sitelist(3,nbt))
!
!nbt = 0
!do kcell = 1,klast
!	site = cell_list(kcell)%site
!	do i = 1,27
!		if (i == 14) cycle
!		nbsite = site + jumpvec(:,i)
!		if (occupancy(nbsite(1),nbsite(2),nbsite(3))%indx(1) == OUTSIDE_TAG) then
!			nbt = nbt+1
!			sitelist(:,nbt) = site
!			exit
!		endif
!	enddo
!enddo
!	
!
!do kcell = klast+1,initial_count
!	r2min = 1.0e10
!	do kb = 1,nbt
!		site = sitelist(:,kb)
!		do i = 1,27
!			if (i == 14) cycle
!			nbsite = site + jumpvec(:,i)
!			if (occupancy(nbsite(1),nbsite(2),nbsite(3))%indx(1) == OUTSIDE_TAG) then
!				r2 = (nbsite(1) - blob_centre(1))**2 + (nbsite(2) - blob_centre(2))**2 + (nbsite(3) - blob_centre(3))**2
!				if (r2 < r2min) then
!					kbmin = kb
!					imin = i
!					r2min = r2
!				endif
!			endif
!		enddo
!	enddo
!	site = sitelist(:,kbmin) + jumpvec(:,imin)
!	call AddCell(kcell,site)
!enddo
!
!deallocate(sitelist)
!		
!end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine InitConcs(ichemo)
integer :: ichemo
!integer :: i, kcell, site(3), ntvars
real(REAL_KIND) :: c0

if (istep == 0) then
	if (ichemo == OXYGEN .or. ichemo == GLUCOSE .or. ichemo == LACTATE .or. ichemo == TRACER) then
		c0 = chemo(ichemo)%bdry_conc
	else	! drug or metabolite
		c0 = 0
	endif
	Caverage(ichemo) = c0
endif
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine ProcessEvent(radiation_dose)
real(REAL_KIND) :: radiation_dose
integer :: kevent, ichemo, idrug, im, nmetab
real(REAL_KIND) :: V, C(MAX_CHEMO)
type(event_type) :: E
logical :: full

!write(logmsg,*) 'ProcessEvent'
!call logger(logmsg)
radiation_dose = 0
do kevent = 1,Nevents
	E = event(kevent)
	if (t_simulation >= E%time .and. .not.E%done) then
		write(nflog,'(a,i3,2f8.0,i3,2f10.4)') 'Event: ',E%etype,t_simulation,E%time,E%ichemo,E%volume,E%conc
		if (E%etype == RADIATION_EVENT) then
			radiation_dose = E%dose
			write(logmsg,'(a,f8.0,f8.3)') 'RADIATION_EVENT: time, dose: ',t_simulation,E%dose
			call logger(logmsg)
		elseif (E%etype == MEDIUM_EVENT) then
			write(logmsg,'(a,f8.0,f8.3,2f8.4)') 'MEDIUM_EVENT: time, volume, O2medium: ',t_simulation,E%volume,E%O2medium
			call logger(logmsg)
			C = 0
			C(OXYGEN) = E%O2medium
			C(GLUCOSE) = E%glumedium
			C(LACTATE) = E%lacmedium
			C(DRUG_A:DRUG_A+5) = 0
			V = E%volume
			full = E%full
			call MediumChange(V,C,full)
		elseif (E%etype == DRUG_EVENT) then
			C = 0
			C(OXYGEN) = E%O2conc
			C(GLUCOSE) = chemo(GLUCOSE)%dose_conc
			C(LACTATE) = chemo(LACTATE)%dose_conc
			ichemo = E%ichemo
			idrug = E%idrug
			C(ichemo) = E%conc
			V = E%volume
			write(logmsg,'(a,2f8.3)') 'DRUG_EVENT: volume, conc: ',E%volume,E%conc
			call logger(logmsg)
			! set %present
			chemo(ichemo)%present = .true.
			chemo(ichemo)%bdry_conc = 0
			nmetab = drug(idrug)%nmetabolites
			do im = 1,nmetab
				if (chemo(ichemo + im)%used) then
					chemo(ichemo + im)%present = .true.
					chemo(ichemo + im)%bdry_conc = 0
				endif
			enddo
			full = E%full
			call MediumChange(V,C,full)
			call UpdateChemomap
!			call UpdateCbnd(0.0d0)
		endif
		event(kevent)%done = .true.
	endif
enddo
end subroutine

!----------------------------------------------------------------------------------
! Radiation treatment is stored in protocol(0)
! NOT USED NOW
!----------------------------------------------------------------------------------
subroutine Treatment(radiation_dose)
real(REAL_KIND) :: radiation_dose
integer :: i, idrug, ichemo, nmetab, im	!, ichemo_metab

radiation_dose = 0
do i = 1,protocol(0)%n
	if (t_simulation >= protocol(0)%tstart(i) .and. .not.protocol(0)%started(i)) then
		radiation_dose = protocol(0)%dose(i)
		protocol(0)%started(i) = .true.
		protocol(0)%ended(i) = .true.
		write(nflog,*) 'Radiation started: dose: ',radiation_dose
		exit
	endif
enddo
do idrug = 1,2
	ichemo = protocol(idrug)%ichemo
	if (idrug == 1) then
		nmetab = 2
	elseif (idrug == 2) then
		nmetab = 2
	endif
	do i = 1,protocol(idrug)%n
		if (i == 1 .and. t_simulation < protocol(idrug)%tstart(i)) then
			chemo(ichemo)%bdry_conc = 0
			do im = 1,nmetab
				chemo(ichemo + im)%bdry_conc = 0
			enddo
			exit
		endif
		if (t_simulation >= protocol(idrug)%tstart(i) .and. .not.protocol(idrug)%started(i)) then
			chemo(ichemo)%bdry_conc = protocol(idrug)%conc(i)
			protocol(idrug)%started(i) = .true.
			protocol(idrug)%ended(i) = .false.
			chemo(ichemo)%present = .true.
			call InitConcs(ichemo)
			call SetupMedium(ichemo)
			do im = 1,nmetab
				if (chemo(ichemo + im)%used) then
					chemo(ichemo + im)%present = .true.
					call InitConcs(ichemo + im)
					call SetupMedium(ichemo + im)
				endif
			enddo
			write(nflog,*) 'Started DRUG: ',chemo(ichemo)%name,chemo(ichemo)%bdry_conc, i
			exit
		endif
	enddo
	do i = 1,protocol(idrug)%n
		if (t_simulation >= protocol(idrug)%tend(i) .and. .not.protocol(idrug)%ended(i)) then
			chemo(ichemo)%bdry_conc = 0
			protocol(idrug)%ended(i) = .true.
			call InitConcs(ichemo)
			call SetupMedium(ichemo)
			do im = 1,nmetab
				chemo(ichemo + im)%bdry_conc = 0
				if (chemo(ichemo + im)%used) then
					call InitConcs(ichemo + im)
					call SetupMedium(ichemo + im)
				endif
			enddo
			write(nflog,*) 'Ended DRUG: ',chemo(ichemo)%name,i
			exit
		endif
	enddo
enddo	
end subroutine

!-----------------------------------------------------------------------------------------
! If the volume removed is Vr, the fraction of constituent mass that is retained
! in the medium is (Vm - Vr)/Vm.  The calculation does not apply to oxygen.
! Usually Vr = Ve.
! Revised treatment of concentrations, to allow for setting O2 in medium change
!
! Now only medium concentrations are stored, in Caverage(MAX_CHEMO+1:2*MAX_CHEMO)
! (oxygen is a special case, Caverage is actually conc at well bottom)
!-----------------------------------------------------------------------------------------
subroutine MediumChange(Ve,Ce,full)
real(REAL_KIND) :: Ve, Ce(:)
logical :: full
real(REAL_KIND) :: R, Vm, Vr, Vcells, mass(MAX_CHEMO), C
integer :: ichemo, idrug, im, iparent

write(nflog,*) 'MediumChange:'
write(nflog,'(a,f8.4)') 'Ve: ',Ve
write(nflog,'(a,13f8.4)') 'Ce: ',Ce
write(nflog,'(a,13e12.3)')'medium_M: ',chemo(OXYGEN+1:)%medium_M
if (full) then
	total_volume = Ve
	Caverage(MAX_CHEMO+1:2*MAX_CHEMO) = Ce
else
	Vcells = Ncells*Vcell_cm3
	Vm = total_volume - Vcells
	Vr = min(Vm,Ve)
	!write(nflog,'(a,4f8.4)') 'total_volume, Vcells, Vm, Vr: ',total_volume, Vcells, Vm, Vr 
	mass = (Vm - Vr)*Caverage(MAX_CHEMO+1:2*MAX_CHEMO) + Ve*Ce(:)
	total_volume = Vm - Vr + Ve + Vcells
	Caverage(MAX_CHEMO+1:2*MAX_CHEMO) = mass/(total_volume - Vcells)
endif
chemo(OXYGEN)%bdry_conc = Ce(OXYGEN)
do ichemo = GLUCOSE,LACTATE
	C = Caverage(MAX_CHEMO+ichemo)
	Caverage(ichemo) = C
	chemo(ichemo)%Cmedium = C
	chemo(ichemo)%bdry_conc = C
enddo
do idrug = 1,2
	iparent = DRUG_A + 3*(idrug-1)
	do im = 0,2
		ichemo = iparent + im	
		chemo(ichemo)%Cmedium = Caverage(MAX_CHEMO+ichemo)
	enddo
enddo

call SetConstLevels	
do ichemo = 1,3
	C_OGL(ichemo,:) = chemo(ichemo)%Cmedium(:)
enddo
write(nflog,'(a,3e12.3)') 'Const Cmedium: ',Caverage(MAX_CHEMO+1:MAX_CHEMO+3)
write(nflog,'(a,2e12.3)') 'Drug Cmedium:  ',Caverage(MAX_CHEMO+DRUG_A:MAX_CHEMO+DRUG_A+1)
write(nflog,*) 'glucose bdry_conc: ',chemo(GLUCOSE)%bdry_conc
t_lastmediumchange = istep*DELTA_T
medium_change_step = .true.
end subroutine

!-----------------------------------------------------------------------------------------
! Total cell O2 flux determines O2 at the well bottom Cex, but flux is determined by Cex.
! Use the current values for Cin and Cex to get the mass flux, which is equated to the
! flux that corresponds to the area, Kdiff and concentration gradient (Cbnd - Cex)/depth.
! Kd.A.(Cbnd - Cex)/d = Ncells.(Kin.Cex - Kout.Cin)
! => Cex = (A.Kd.Cbnd/d + Ncells.Kout.Cin)/(A.Kd/d + Ncells.Kin) 
!-----------------------------------------------------------------------------------------
subroutine SetConstLevels
integer :: ichemo, k, kcell, idrug, iparent, im
real(REAL_KIND) :: Kin, Kout, Kd, Cex, Cin, Cbnd, A, d, flux, Cin_prev, alpha
real(REAL_KIND) :: tol = 1.0e-6

!ichemo = OXYGEN

do ichemo = 1,3
if (chemo(ichemo)%constant .or. fully_mixed) then
    Cex = chemo(ichemo)%bdry_conc
    Cin = getCin(ichemo,Cex)
else
    Kin = chemo(ichemo)%membrane_diff_in
    Kout = chemo(ichemo)%membrane_diff_out
    Cex = Caverage(MAX_CHEMO+ichemo)
    Cin = Caverage(ichemo)
    Cbnd = chemo(ichemo)%bdry_conc
    Kd = chemo(ichemo)%medium_diff_coef
    A = well_area
    d = total_volume/A
	if (ichemo == OXYGEN) then
		do k = 1,100
			Cin_prev = Cin
			Cex = (A*Kd*Cbnd/d + Ncells*Kout*Cin)/(A*Kd/d + Ncells*Kin)
			Cin = getCin(ichemo,Cex)
		!    write(*,'(a,i4,2e15.6)') 'SetMediumOxygen: ',k,Cin,Cex
			if (abs(Cin-Cin_prev)/Cin_prev < tol) exit
		enddo
    endif
endif
Caverage(ichemo) = Cin
Caverage(MAX_CHEMO+ichemo) = Cex
do k = 1,N1D
	alpha = real(k-1)/(N1D-1)
	chemo(ichemo)%Cmedium(k) = alpha*chemo(ichemo)%bdry_conc + (1-alpha)*Cex
enddo
do kcell = 1,nlist
    if (cell_list(kcell)%state == DEAD) cycle
    cell_list(kcell)%Cin(ichemo) = Cin
!	do idrug = 1,2
!		iparent = DRUG_A + 3*(idrug-1)
!		do im = 0,2
!			ichemo = iparent + im
!			if (.not.chemo(ichemo)%present) cycle
!			cell_list(kcell)%Cin(ichemo) = chemo(ichemo)%Cmedium(1)		! set IC conc to initial medium conc 
!		enddo
!	enddo
enddo
write(nflog,'(a,i4,2e12.3)') 'SetConstLevels: Cex, Cin: ',ichemo,Cex,Cin
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Advance simulation through one big time step (DELTA_T)
! The concentration fields are first solved through NT_CONC subdivisions of DELTA_T,
! then the cell states are updated, which includes cell death and division.
! On either death or division cell positions are adjusted, and site concentrations
! (if necessary), and the ODE solver mappings in ODEdiff.
!-----------------------------------------------------------------------------------------
subroutine simulate_step(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: simulate_step  
use, intrinsic :: iso_c_binding
integer(c_int) :: res
integer :: kcell, site(3), hour, nthour, kpar=0
real(REAL_KIND) :: r(3), rmax, tstart, dt, dts, radiation_dose, diam_um, framp, area, diam
!integer, parameter :: NT_CONC = 6
integer :: i, ic, ichemo, ndt, iz, idrug, ityp, idiv, ndiv, Nmetabolisingcells
integer :: nvars, ns
real(REAL_KIND) :: dxc, ex_conc(120*CYCLE_PHASE+1)		! just for testing
real(REAL_KIND) :: DELTA_T_save, t_sim_0
real(REAL_KIND) :: HIF1, PDK1
type(metabolism_type), pointer :: mp
logical :: ok = .true.
logical :: dbug
type(metabolism_type), pointer :: metabolic
	
metabolic => phase_metabolic(1)
!call testmetab2

dbug = (istep < 0)
Nmetabolisingcells = Ncells - (Ndying(1) + Ndying(2))
!if (Nmetabolisingcells == 0) then
!	call logger('# of metabolising cells = 0')
!    res = 0
!    return
!endif
if (Ncells == 0) then
	call logger('# of cells = 0')
    res = 2
    return
endif
if (limit_stop) then
	call logger('Spheroid size limit reached')
	res = 6
	return
endif

nthour = 3600/DELTA_T

if (istep == -100) then
	stop
endif

if (ngaps > 200) then
	call squeezer
endif

drug_gt_cthreshold = .false.

if (medium_change_step .or. chemo(DRUG_A)%present .or. chemo(DRUG_B)%present) then
	ndiv = 6
else
	ndiv = 1
endif
dt = DELTA_T/ndiv
dts = dt/NT_CONC
DELTA_T_save = DELTA_T
DELTA_T = dt
t_sim_0 = t_simulation
do idiv = 0,ndiv-1
	t_simulation = t_sim_0 + idiv*DELTA_T

	if (dbug) write(nflog,*) 'Solver'
	do it_solve = 1,NT_CONC
		tstart = (it_solve-1)*dts
	!	t_simulation = (istep-1)*DELTA_T + tstart
		t_simulation = t_simulation + tstart
		call Solver(it_solve,t_simulation,dts,Ncells,ok)
		if (.not.ok) then
			res = 5
			return
		endif
	enddo	! end it_solve loop
	if (dbug) write(nflog,*) 'did Solver'
	if (use_metabolism) then
		do ityp = 1,Ncelltypes
			HIF1 = metabolic%HIF1
			call analyticSetHIF1(Caverage(OXYGEN),HIF1,DELTA_T)
			metabolic%HIF1 = HIF1
			PDK1 = metabolic%PDK1
			call analyticSetPDK1(HIF1,PDK1,dt)
			metabolic%PDK1 = PDK1
		enddo
	endif
	!write(nflog,*) 'did Solver'
	call CheckDrugConcs
	call CheckDrugPresence

	if (dbug) write(nflog,*) 'GrowCells'
	call GrowCells(DELTA_T,t_simulation,ok)
	if (dbug) write(nflog,*) 'did GrowCells'
	if (.not.ok) then
		res = 3
		return
	endif
enddo	! end idiv loop

DELTA_T = DELTA_T_save
medium_change_step = .false.

!istep = istep + 1
t_simulation = (istep-1)*DELTA_T	! seconds

!!write(nflog,*) 'GrowCells'
!call GrowCells(DELTA_T,t_simulation,ok)
!!write(nflog,*) 'did GrowCells'
!if (.not.ok) then
!	res = 3
!	return
!endif

radiation_dose = 0
if (use_treatment) then     ! now we use events
	call treatment(radiation_dose)
endif
if (use_events) then
	call ProcessEvent(radiation_dose)
endif
if (radiation_dose > 0) then
	write(logmsg,'(a,f6.1)') 'Radiation dose: ',radiation_dose
	call logger(logmsg)
	call Irradiation(radiation_dose, ok)
	if (.not.ok) then
		res = 3
		return
	endif
endif

res = 0

!call test_CellDivision
!if (.not.use_TCP .and. (mod(istep,6) == 0)) then
!	call get_concdata(nvars, ns, dxc, ex_conc)
!endif

call getNviable

if (saveFACS%active) then
	if (istep*DELTA_T >= saveFACS%tstart + (saveFACS%it-1)*saveFACS%dt) then
		call WriteFACSData
		saveFACS%it = saveFACS%it + 1
		if (saveFACS%it > saveFACS%nt) then
			saveFACS%active = .false.
		endif
	endif
endif

if (dbug .or. mod(istep,nthour) == 0) then
!	mp => metabolic
	mp => phase_metabolic(1)
	write(logmsg,'(a,i6,i4,a,i8,a,i8)') 'did istep, hour: ',istep,istep/nthour,' Nlive: ',Ncells,'   Nviable: ',sum(Nviable)
	call logger(logmsg)
!	write(logmsg,'(a,4e12.3)') 'G_rate, A_rate, PO_rate, O_rate: ',mp%G_rate,mp%A_rate,mp%P_rate,mp%O_rate
!	call logger(logmsg)
!	write(logmsg,'(a,3e12.3)') 'C_O2, C_G, C_L: ',Caverage(1:3)
!	call logger(logmsg)
!	write(*,'(a,3f8.4)') 'lactate, HIF1: ',Caverage(LACTATE),Caverage(LACTATE+MAX_CHEMO),HIF1
!	call showcells 
!	if (use_metabolism) then
!		call show_metabolism(1)
!	endif
endif
! write(nflog,'(a,f8.3)') 'did simulate_step: time: ',wtime()-start_wtime

istep = istep + 1
!call averages
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine averages
real(REAL_KIND) :: ave_V, ave_dVdt, ave_fg
integer :: kcell, n
type(cell_type), pointer :: cp

ave_V = 0
ave_dVdt = 0
ave_fg = 0
n = 0
do kcell = 1,Ncells
	cp => cell_list(kcell)
	if (cp%state == DEAD) continue
	n = n+1
	ave_V = ave_V + cp%V
	ave_dVdt = ave_dVdt + cp%dVdt
	ave_fg = ave_fg + cp%fg
enddo
write(nflog,'(a,3e12.3)') 'averages: V,dVdt,fg: ',ave_V/n,ave_dVdt/n, ave_fg/n
end subroutine

!-----------------------------------------------------------------------------------------
! 
!-----------------------------------------------------------------------------------------
subroutine show_metabolism(kcell)
integer :: kcell
type(cell_type), pointer :: cp
type (metabolism_type), pointer :: mp
real(REAL_KIND) :: Cin(MAX_CHEMO)
!real(REAL_KIND) :: HIF1, G_rate, PP_rate, P_rate
!real(REAL_KIND) :: L_rate, A_rate, I_rate, O_rate

cp =>cell_list(kcell)
mp => cp%metab
Cin = Caverage(1:MAX_CHEMO)
!write(*,'(a,3f8.4)') 'O2, glucose, lactate: ',Cin(1:3) 
call get_metab_rates(mp, Cin)
return

write(*,'(a,i2,3e12.3)') 'phase, V: ',cp%phase,cp%V		!I2Divide,Itotal,mp%Itotal,mp%I2Divide
write(*,'(a,3e11.3)') 'G_rate, P_rate, O_rate: ',mp%G_rate, mp%P_rate, mp%O_rate
write(*,'(a,4e11.3)') 'L_rate, A_rate, I_rate: ',mp%L_rate, mp%A_rate, mp%I_rate
write(*,'(a,4f8.4)') 'O2, glucose, lactate, H: ',Caverage(OXYGEN),Caverage(GLUCOSE),Caverage(LACTATE),mp%HIF1
if (mp%A_rate < ATPg) then
	write(*,*) 'Not growing'
endif
write(*,*)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine showcell(kcell)
integer :: kcell
type(cell_type), pointer :: cp

cp => cell_list(kcell)
write(nflog,'(a,i6,4e12.3)') 'kcell, volume, divide_volume, dVdt, divide_time: ', &
                kcell, cp%V, cp%divide_volume, cp%dVdt, cp%divide_time
end subroutine

!-----------------------------------------------------------------------------------------
! Average volumes etc
!-----------------------------------------------------------------------------------------
subroutine showcells
integer :: kcell, n
real(REAL_KIND) :: Vsum,divVsum,dVdtsum,divtsum
!real(REAL_KIND) :: Vn   ! to normalise volumes
type(cell_type), pointer :: cp

!Vn = Vdivide0/1.6
Vsum=0
divVsum=0
dVdtsum=0
divtsum=0
n=0
do kcell = 1,nlist
    cp => cell_list(kcell)
    if (cp%state == DEAD) cycle
    n = n+1
    Vsum = Vsum + cp%V
    divVsum = divVsum + cp%divide_volume
    dVdtsum = dVdtsum + cp%dVdt
    divtsum = divtsum + cp%divide_time
enddo
write(nflog,'(a,4e12.3)') 'ave volume, divide_volume, dVdt, divide_time: ', Vsum/n,divVsum/n,dVdtsum/n,divtsum/n
!write(*,'(a,4e12.3)') 'ave volume, divide_volume, dVdt, divide_time: ', Vsum/n,divVsum/n,dVdtsum/n,divtsum/n
end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Execute(ncpu,infile_array,inbuflen,outfile_array,outbuflen,res) BIND(C) 
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
character(c_char) :: infile_array(*), outfile_array(*)
integer(c_int) :: ncpu, inbuflen, outbuflen, res
character*(2048) :: infile, outfile
logical :: ok, success, isopen
integer :: i

!use_PEST = (.not.use_TCP .and. PEST_outputfile(1:1) /= ' ')
if (use_TCP) then
	use_PEST = .false.
endif
infile = ''
do i = 1,inbuflen
	infile(i:i) = infile_array(i)
enddo
outfile = ''
do i = 1,outbuflen
	outfile(i:i) = outfile_array(i)
enddo

inquire(unit=nflog,OPENED=isopen)
if (isopen) then
	close(nflog)
endif
open(nflog,file='vmonolayer.log',status='replace')

#ifdef GFORTRAN
    write(logmsg,'(a)') 'Built with GFORTRAN'
	call logger(logmsg)
#endif

logmsg = 'OS??'
#ifdef LINUX
    write(logmsg,'(a)') 'OS is Linux'
#endif
#ifdef OSX
    write(logmsg,'(a)') 'OS is OS-X'
#endif
#ifdef _WIN32
    write(logmsg,'(a)') 'OS is Windows'
#endif
#ifdef WINDOWS
    write(logmsg,'(a)') 'OS is Windows'
#endif
call logger(logmsg)

!#ifdef OPENMP
#if defined(OPENMP) || defined(_OPENMP)
    write(logmsg,'(a)') 'Executing with OpenMP'
	call logger(logmsg)
#endif

write(logmsg,*) 'inputfile:  ', trim(infile)
call logger(logmsg)
write(logmsg,*) 'outputfile: ', trim(outfile)
call logger(logmsg)
if (use_TCP) then
	call connecter(ok)
	if (.not.ok) then
		call logger('Failed to make TCP connections')
		return
	endif
endif

DELTA_T = 600
nsteps = 100
res=0

call Setup(ncpu,infile,outfile,ok)
if (ok) then
!	clear_to_send = .true.
!	simulation_start = .true.
else
	call logger('=== Setup failed ===')
endif
if (ok) then
	res = 0
else
	res = 1
endif
execute_t1 = wtime()

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine DisableTCP
!DEC$ ATTRIBUTES DLLEXPORT :: disableTCP
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"DISABLETCP" :: disableTCP

use_TCP = .false.   ! because this is called from monolayer_main()	
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Connection(awp,port,error)
TYPE(winsockport) :: awp
integer :: port, error
integer :: address = 0
!!!character*(64) :: ip_address = "127.0.0.1"C      ! need a portable way to make a null-terminated C string
character*(64) :: host_name = "localhost"

if (.not.winsock_init(1)) then
    call logger("winsock_init failed")
    stop
endif

awp%handle = 0
awp%host_name = host_name
awp%ip_port = port
awp%protocol = IPPROTO_TCP
call Set_Winsock_Port (awp,error)

if (.not.awp%is_open) then
    write(nflog,*) 'Error: connection: awp not open: ',port
else
    write(nflog,*) 'connection: awp open: ',port, error
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Connecter(ok)
logical :: ok
integer :: error

! Main connection
ok = .true.
error = 0
call Connection(awp_0,TCP_PORT_0,error)
if (awp_0%handle < 0 .or. error /= 0) then
    write(logmsg,'(a)') 'TCP connection to TCP_PORT_0 failed'
    call logger(logmsg)
    ok = .false.
    return
endif
if (.not.awp_0%is_open) then
	write(logmsg,'(a)') 'No connection to TCP_PORT_0'
    call logger(logmsg)
    ok = .false.
    return
endif
write(logmsg,'(a)') 'Connected to TCP_PORT_0  '
call logger(logmsg)

if (use_CPORT1) then
	call connection(awp_1,TCP_PORT_1,error)
	if (awp_1%handle < 0 .or. error /= 0) then
		write(logmsg,'(a)') 'TCP connection to TCP_PORT_1 failed'
		call logger(logmsg)
		ok = .false.
		return
	endif
	if (.not.awp_1%is_open) then
		write(logmsg,'(a)') 'No connection to TCP_PORT_1'
		call logger(logmsg)
		ok = .false.
		return
	endif
	write(logmsg,'(a)') 'Connected to TCP_PORT_1  '
	call logger(logmsg)
endif
! Allow time for completion of the connection
call sleeper(2)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine terminate_run(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: terminate_run 
use, intrinsic :: iso_c_binding
integer(c_int) :: res
character*(8), parameter :: quit = '__EXIT__'
integer :: error, i

call logger('terminate_run')
call Wrapup

if (res == 0) then
	call logger(' Execution successful!')
elseif (res == -1) then
	call logger(' Execution stopped')
elseif (res == 2) then
	call logger(' No more live cells')
elseif (res == 6) then
	call logger(' Spheroid size limit reached')
elseif (res == 3) then
	call logger(' === Execution failed === ERROR in GrowCells')
elseif (res == 4) then
	call logger(' === Execution failed === ERROR in diff_solver')
elseif (res == 5) then
	call logger(' === Execution failed === ERROR in Solver')
endif
write(logmsg,'(a,f10.2)') 'Execution time (min): ',(wtime() - execute_t1)/60
call logger(logmsg)

!close(nflog)

if (use_TCP) then
	if (stopped) then
	    call winsock_close(awp_0)
	    if (use_CPORT1) call winsock_close(awp_1)
	else
	    call winsock_send(awp_0,quit,8,error)
	    call winsock_close(awp_0)
		if (use_CPORT1) then
			call winsock_send(awp_1,quit,8,error)
			call winsock_close(awp_1)
		endif
	endif
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Wrapup
integer :: ierr, ichemo, idrug
logical :: isopen

call logger('doing wrapup ...')
ierr = 0
!if (allocated(zoffset)) deallocate(zoffset)
!if (allocated(zdomain)) deallocate(zdomain)
if (allocated(gaplist)) deallocate(gaplist,stat=ierr)
!if (allocated(occupancy)) deallocate(occupancy)
!if (allocated(cell_list)) deallocate(cell_list)
!if (allocated(allstate)) deallocate(allstate)
!if (allocated(allstatep)) deallocate(allstatep)
!if (allocated(work_rkc)) deallocate(work_rkc)
!do ichemo = 1,MAX_CHEMO
!	if (allocated(chemo(ichemo)%coef)) deallocate(chemo(ichemo)%coef)
!	if (allocated(chemo(ichemo)%conc)) deallocate(chemo(ichemo)%conc)
!	if (allocated(chemo(ichemo)%grad)) deallocate(chemo(ichemo)%grad)
!enddo
!if (allocated(ODEdiff%ivar)) deallocate(ODEdiff%ivar)
!if (allocated(ODEdiff%varsite)) deallocate(ODEdiff%varsite)
!if (allocated(ODEdiff%icoef)) deallocate(ODEdiff%icoef)
if (allocated(protocol)) then
	do idrug = 0,2	!<------  change this to a variable
		if (allocated(protocol(idrug)%tstart)) deallocate(protocol(idrug)%tstart)
		if (allocated(protocol(idrug)%tend)) deallocate(protocol(idrug)%tend)
		if (allocated(protocol(idrug)%conc)) deallocate(protocol(idrug)%conc)
		if (allocated(protocol(idrug)%dose)) deallocate(protocol(idrug)%dose)
		if (allocated(protocol(idrug)%started)) deallocate(protocol(idrug)%started)
		if (allocated(protocol(idrug)%ended)) deallocate(protocol(idrug)%ended)
	enddo
	deallocate(protocol)
endif
call logger('deallocated all arrays')

! Close all open files
inquire(unit=nfout,OPENED=isopen)
if (isopen) then
	close(nfout)
	call logger('closed nfout')
endif
inquire(nfres,OPENED=isopen)
if (isopen) close(nfres)
!inquire(nfPESTout,OPENED=isopen)
!if (isopen) close(nfPESTout)
call logger('closed files')

if (par_zig_init) then
	call par_zigfree
endif
call logger('freed par_zig')
end subroutine

end module
