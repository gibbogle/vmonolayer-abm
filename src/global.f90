! Global definitions

module global

use omp_lib
use real_kind_mod
use par_zig_mod
use winsock
use, intrinsic :: ISO_C_BINDING

implicit none

integer, parameter :: TCP_PORT_0 = 5000		! main communication port (logging)
integer, parameter :: TCP_PORT_1 = 5001		! data transfer port (plotting)
integer, parameter :: NORMAL_DIST      = 1
integer, parameter :: LOGNORMAL_DIST   = 2
integer, parameter :: EXPONENTIAL_DIST = 3
integer, parameter :: CONSTANT_DIST    = 4

!integer, parameter :: DIVIDING  = 1
!integer, parameter :: QUIESCENT = 2
!integer, parameter :: DEAD      = 3
integer, parameter :: ALIVE = 1
integer, parameter :: DYING = 2
integer, parameter :: DEAD = 3

integer, parameter :: OUTSIDE_TAG  = -1
integer, parameter :: UNREACHABLE_TAG  = -2

integer, parameter :: DIVIDE_ALWAYS_PUSH  = 1
integer, parameter :: DIVIDE_USE_CLEAR_SITE  = 2
integer, parameter :: DIVIDE_USE_CLEAR_SITE_RANDOM  = 3

integer, parameter :: nfin=10, nfout=11, nflog=12, nfres=13, nfrun=14, nfcell=15, nftreatment=16, nfFACS=17, &
					  nfpar=18, nftcp=20
integer, parameter :: neumann(3,6) = reshape((/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /), (/3,6/))

integer, parameter :: CFSE = 0
integer, parameter :: OXYGEN = 1
integer, parameter :: GLUCOSE = 2
integer, parameter :: LACTATE = 3
integer, parameter :: GLUTAMINE = 4
integer, parameter :: OTHERNUTRIENT = 5
integer, parameter :: DRUG_A = 6
integer, parameter :: TPZ_DRUG = DRUG_A
integer, parameter :: TPZ_DRUG_METAB_1 = TPZ_DRUG + 1
integer, parameter :: TPZ_DRUG_METAB_2 = TPZ_DRUG + 2
integer, parameter :: DRUG_B = DRUG_A + 3
integer, parameter :: DNB_DRUG = DRUG_B
integer, parameter :: DNB_DRUG_METAB_1 = DNB_DRUG + 1
integer, parameter :: DNB_DRUG_METAB_2 = DNB_DRUG + 2
integer, parameter :: MAX_CHEMO = DRUG_B + 2
integer, parameter :: NUTS = DRUG_A - 1

integer, parameter :: GROWTH_RATE = MAX_CHEMO + 1	! (not used here, used in the GUI)
integer, parameter :: CELL_VOLUME = MAX_CHEMO + 2
integer, parameter :: O2_BY_VOL = MAX_CHEMO + 3
integer, parameter :: CYCLE_PHASE = MAX_CHEMO + 4

integer, parameter :: N_EXTRA = CYCLE_PHASE - MAX_CHEMO + 1	! = 5 = total # of variables - MAX_CHEMO
integer, parameter :: NCONST = MAX_CHEMO

integer, parameter :: TPZ_CLASS = 1
integer, parameter :: DNB_CLASS = 2
integer, parameter :: DRUG_EVENT = 1
integer, parameter :: RADIATION_EVENT = 2
integer, parameter :: MEDIUM_EVENT = 3

integer, parameter :: NTCP = 1000

integer, parameter :: DIST_NV = 20

integer, parameter :: EXTRA = 1
integer, parameter :: INTRA = 2
integer, parameter :: MAX_CELLTYPES = 2
integer, parameter :: MAX_DRUGTYPES = 2
integer, parameter :: max_nlist = 200000
integer, parameter :: NRF = 4
integer, parameter :: LIMIT_THRESHOLD = 1500

logical, parameter :: use_ODE_diffusion = .true.
logical, parameter :: compute_concentrations = .true.
logical, parameter :: use_division = .true.
logical, parameter :: use_death = .true.
logical, parameter :: use_react = .true.
logical, parameter :: use_migration = .false.	! causing an error with vacant site becoming bdry 
logical, parameter :: use_medium_flux = .true.	! flux of constituents between spheroid and medium is accounted for.
logical, parameter :: use_metabolites = .true.
logical, parameter :: use_celltype_colour = .true.

logical, parameter :: use_Cex_Cin = .true.		! assume equilibrium to derive Cin from Cex
logical, parameter :: suppress_growth = .false.

logical, parameter :: OFF_LATTICE = .false.

real(REAL_KIND), parameter :: PI = 4.0*atan(1.0)
real(REAL_KIND), parameter :: CFSE_std = 0.05
real(REAL_KIND), parameter :: small_d = 0.1e-4          ! 0.1 um -> cm

!type occupancy_type
!	integer :: indx(2)
!	real(REAL_KIND) :: C(MAX_CHEMO)
!	type(boundary_type), pointer :: bdry
!	! for FD grid weighting
!	integer :: cnr(3,8)
!	real(REAL_KIND) :: wt(8)
!end type

type metabolism_type
	real(REAL_KIND) :: HIF1
	real(REAL_KIND) :: PDK1
	real(REAL_KIND) :: I_rate_max
	real(REAL_KIND) :: G_rate
	real(REAL_KIND) :: PP_rate
	real(REAL_KIND) :: P_rate 
	real(REAL_KIND) :: L_rate
	real(REAL_KIND) :: A_rate
	real(REAL_KIND) :: I_rate
	real(REAL_KIND) :: O_rate
	real(REAL_KIND) :: Gln_rate
	real(REAL_KIND) :: ON_rate
	real(REAL_KIND) :: GA_rate
	real(REAL_KIND) :: f_G
	real(REAL_KIND) :: f_P
	real(REAL_KIND) :: f_Gln
	real(REAL_KIND) :: C_P
	real(REAL_KIND) :: C_A
	real(REAL_KIND) :: recalcable     ! if > 0, can determine rates from w = f_G/f_Gu, otherwise need full procedure
	real(REAL_KIND) :: C_GlnEx_prev
end type

type cell_type
	integer :: ID
	integer :: celltype
	integer :: site(3)
	integer :: ivin
	logical :: active
	integer :: state
	logical :: Iphase
!    integer :: nspheres             ! =1 for Iphase, =2 for Mphase
!	real(REAL_KIND) :: radius(2)	! sphere radii (um) -> cm
!	real(REAL_KIND) :: centre(3,2)  ! sphere centre positions
!	real(REAL_KIND) :: d			! centre separation distance (um) -> cm
	integer :: generation
!	real(REAL_KIND) :: conc(MAX_CHEMO)
	real(REAL_KIND) :: Cin(MAX_CHEMO)
!	real(REAL_KIND) :: Cex(MAX_CHEMO)
	real(REAL_KIND) :: dCdt(MAX_CHEMO)
	real(REAL_KIND) :: dMdt(MAX_CHEMO)      ! mumol/s
	real(REAL_KIND) :: CFSE
	real(REAL_KIND) :: dVdt             ! actual growth rate
	real(REAL_KIND) :: V			    ! actual volume cm3
	real(REAL_KIND) :: growth_rate_factor	! to introduce some random variation 
	real(REAL_KIND) :: ATP_rate_factor	! to introduce some random variation 
	real(REAL_KIND) :: divide_volume	! actual divide volume
	real(REAL_KIND) :: divide_time      ! cycle time
	real(REAL_KIND) :: fg				! to make sum(T_G1, T_S, T_G2) consistent with Tdivide
	real(REAL_KIND) :: t_divide_last	! these two values are used for colony simulation
	real(REAL_KIND) :: t_divide_next
	real(REAL_KIND) :: birthtime
!	real(REAL_KIND) :: t_anoxia
!	real(REAL_KIND) :: t_anoxia_die
!	real(REAL_KIND) :: t_aglucosia
!	real(REAL_KIND) :: t_aglucosia_die
	real(REAL_KIND) :: M
	real(REAL_KIND) :: p_rad_death
	real(REAL_KIND) :: p_drug_death(MAX_DRUGTYPES)
	real(REAL_KIND) :: t_start_mitosis
	real(REAL_KIND) :: mitosis
!	logical :: growth_delay
!	real(REAL_KIND) :: dt_delay
!	real(REAL_KIND) :: t_growth_delay_end	! this is for suppression of growth before first division
	real(REAL_KIND) :: tag_time             ! time cell is tagged to die metabolically
!	integer :: N_delayed_cycles_left		! decremented by 1 at each cell division
	logical :: radiation_tag, ATP_tag, GLN_tag
	logical :: drug_tag(MAX_DRUGTYPES)
	real(REAL_KIND) :: apoptosis_delay
	logical :: G2_M
!	logical :: exists
!	integer :: cnr(3,8)
!	real(REAL_KIND) :: wt(8)

	! Cell cycle 
    integer :: phase
    logical :: G1_flag, G1S_flag, S_flag, SG2_flag, G2_flag, G2M_flag
    real(REAL_KIND) :: G1_time, S_time, G2_time, S_duration
    real(REAL_KIND) :: G1_V, S_V, G2_V
    real(REAL_KIND) :: G1S_time, SG2_time, G2M_time, M_time
    real(REAL_KIND) :: doubling_time
    logical :: arrested ! for S-phase arrest
    real(REAL_KIND) :: S_start_time	    ! for PI labelling
    integer :: NL1, NL2(2)
    
    integer :: N_PL, N_IRL, N_Ch1, N_Ch2
    logical :: irrepairable
    
    ! exponential cycle time
    real(REAL_KIND) :: G1ch_entry_time, G1ch_time, G1ch_max_delay
    real(REAL_KIND) :: Sch_entry_time, Sch_time, Sch_max_delay
    real(REAL_KIND) :: G2ch_entry_time, G2ch_time, G2ch_max_delay
    
	type(metabolism_type) :: metab
	logical :: have_rates
    
	integer :: ndt

end type

type cycle_parameters_type
    real(REAL_KIND) :: T_G1, T_S, T_G2, T_M
    real(REAL_KIND) :: G1_mean_delay, S_mean_delay, G2_mean_delay
    real(REAL_KIND) :: Pk_G1, Pk_S, Pk_G2
    real(REAL_KIND) :: arrest_threshold
    ! Radiation damage/repair
    real(REAL_KIND) :: eta_PL, eta_IRL
    real(REAL_KIND) :: HRR_repair_base, HRR_repair_max
    real(REAL_KIND) :: NHEJ_repair, NHEJ_misrepair
    real(REAL_KIND) :: DIM_misrepair
    real(REAL_KIND) :: fraction_Ch1, psurvive_Ch1, psurvive_Ch2		! prob of surviving mitosis
    real(REAL_KIND) :: psurvive_PL
    real(REAL_KIND) :: aTCP, bTCP
    ! Apoptosis
    real(REAL_KIND) :: apoptosis_rate
    real(REAL_KIND) :: apoptosis_median
    real(REAL_KIND) :: apoptosis_shape
    real(REAL_KIND) :: f_apoptosis_rate_lo  ! multiplying factor for low apoptosis rate
    real(REAL_KIND) :: t_apoptosis_hi       ! duration of high rate of apoptosis
end type

type drug_type
	character*(3)   :: classname
	integer         :: drugclass
	character*(16)  :: name
	integer         :: nmetabolites
	logical         :: use_metabolites
	real(REAL_KIND) :: diff_coef(0:2)
	real(REAL_KIND) :: medium_diff_coef(0:2)
	real(REAL_KIND) :: membrane_diff_in(0:2)
	real(REAL_KIND) :: membrane_diff_out(0:2)
	real(REAL_KIND) :: halflife(0:2)
	logical         :: kills(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Kmet0(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: C2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: KO2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: n_O2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Vmax(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Km(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Klesion(MAX_CELLTYPES,0:2)
	integer         :: kill_model(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: death_prob(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Kd(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_O2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_drug(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_duration(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_fraction(MAX_CELLTYPES,0:2)
	logical         :: sensitises(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: SER_max(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: SER_Km(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: SER_KO2(MAX_CELLTYPES,0:2)
	
	logical :: phase_dependent
	logical :: active_phase(6)
end type

!type boundary_type
!    integer :: site(3)
!    type (boundary_type), pointer :: next
!end type

type dist_type
	integer :: class
	real(REAL_KIND) :: p1, p2, p3
end type

!type, bind(C) :: field_data
!	integer(c_int) :: site(3)
!	integer(c_int) :: state
!	real(c_double) :: volume
!	real(c_double) :: conc(0:MAX_CHEMO+N_EXTRA)	! Must agree with the definition in field.h: 0 = CFSE, MAX_CHEMO+1 = dVdt, MAX_CHEMO+2 = cellvolume
!end type

type, bind(C) :: dist_data
	logical(c_bool) :: used
	real(c_double) :: dv
	real(c_double) :: v0
	real(c_double) :: prob(DIST_NV)
end type

type treatment_type
	integer :: ichemo
	integer :: n
!	character*(16) :: name
	real(REAL_KIND), allocatable :: tstart(:)
	real(REAL_KIND), allocatable :: tend(:)
	real(REAL_KIND), allocatable :: conc(:)
	real(REAL_KIND), allocatable :: dose(:)
	logical, allocatable :: started(:)
	logical, allocatable :: ended(:)
end type

type event_type
	integer :: etype
	real(REAL_KIND) :: time
	integer :: idrug				! DRUG
	integer :: ichemo				! DRUG CHEMO INDEX
	real(REAL_KIND) :: volume		! DRUG MEDIUM
	real(REAL_KIND) :: conc			! DRUG
	real(REAL_KIND) :: O2conc		! DRUG
	real(REAL_KIND) :: O2flush		! DRUG
	real(REAL_KIND) :: dose			! RADIATION
	real(REAL_KIND) :: O2medium		! MEDIUM
	real(REAL_KIND) :: glumedium	! MEDIUM
	real(REAL_KIND) :: lacmedium	! MEDIUM
	logical :: full
	logical :: done
end type	

type LQ_type
	real(REAL_KIND) :: OER_am, OER_bm
	real(REAL_KIND) :: alpha_H, beta_H
	real(REAL_KIND) :: K_ms
	real(REAL_KIND) :: death_prob
!	real(REAL_KIND) :: growth_delay_factor
!	real(REAL_KIND) :: growth_delay_N
end type

type savedata_type
    logical :: active
    character*(128) :: filebase
    real(REAL_KIND) :: dt, tstart
    integer :: nt, it
end type

type(dist_type) :: divide_dist(MAX_CELLTYPES)
!type(occupancy_type), allocatable :: occupancy(:,:,:)
type(cell_type), allocatable, target :: cell_list(:)
type(treatment_type), allocatable :: protocol(:)
type(event_type), allocatable :: event(:)
!real(REAL_KIND), allocatable, target :: Cslice(:,:,:,:)
type(cell_type), target, allocatable :: ccell_list(:)
integer, allocatable :: perm_index(:)
type(cell_type), target :: master_cell

character*(12) :: dll_version, dll_run_version
character*(12) :: gui_version, gui_run_version
integer :: initial_count

integer :: nlist, Ncells, Ncells0, ncells_mphase, lastID, Ncelltypes
integer :: nColonyMax
integer :: Ncells_type(MAX_CELLTYPES), Ndying(MAX_CELLTYPES), Nviable(MAX_CELLTYPES), Ndead(MAX_CELLTYPES)
logical :: limit_stop
!integer :: nadd_sites, Nsites, Nreuse
integer :: Ndrugs_used
integer :: Nradiation_tag(MAX_CELLTYPES), NATP_tag(MAX_CELLTYPES), NGLN_tag(MAX_CELLTYPES)
integer :: Ndrug_tag(MAX_DRUGTYPES,MAX_CELLTYPES)
integer :: Nradiation_dead(MAX_CELLTYPES), NATP_dead(MAX_CELLTYPES), NGLN_dead(MAX_CELLTYPES)
integer :: Ndrug_dead(MAX_DRUGTYPES,MAX_CELLTYPES)
!logical :: use_radiation_growth_delay_all = .true.

integer :: ndoublings
real(REAL_KIND) :: doubling_time_sum

type(cycle_parameters_type), target :: cc_parameters(MAX_CELLTYPES)

logical :: drug_gt_cthreshold(MAX_DRUGTYPES)
real(REAL_KIND) :: Cthreshold
real(REAL_KIND) :: Clabel_threshold		! for labelling drug, e.g. EDU

type(savedata_type) :: saveFACS	!, saveprofile, saveslice

integer :: istep, nsteps, it_solve, NT_CONC, NT_GUI_OUT, show_progeny, ichemo_curr, NT_DISPLAY
integer :: Mnodes, ncpu_input
integer :: Nevents
real(REAL_KIND) :: DELTA_T, DELTA_X, fluid_fraction, Vsite_cm3, Vextra_cm3, Vcell_pL, tnow, DT_DISPLAY
!real(REAL_KIND) :: dxb, dxb3, dxf, dx3
!real(REAL_KIND) :: grid_offset(3)
real(REAL_KIND) :: Vcell_cm3, medium_volume0, total_volume, well_area, t_lastmediumchange
real(REAL_KIND) :: celltype_fraction(MAX_CELLTYPES)
integer :: selected_celltype
logical :: celltype_display(MAX_CELLTYPES)
real(REAL_KIND) :: MM_THRESHOLD, anoxia_threshold, t_anoxia_limit, anoxia_death_delay, Vdivide0, dVdivide
real(REAL_KIND) :: aglucosia_threshold, t_aglucosia_limit, aglucosia_death_delay, max_growthrate(MAX_CELLTYPES)
real(REAL_KIND) :: divide_time_median(MAX_CELLTYPES), divide_time_shape(MAX_CELLTYPES), divide_time_mean(MAX_CELLTYPES)
real(REAL_KIND) :: t_simulation, execute_t1
real(REAL_KIND) :: O2cutoff(3), hypoxia_threshold
real(REAL_KIND) :: growthcutoff(3)
real(REAL_KIND) :: spcrad_value
real(REAL_KIND) :: total_dMdt
!real(REAL_KIND) :: total_flux_prev, medium_Cbnd_prev
real(REAL_KIND) :: start_wtime

! Metabolism parameters
real(REAL_KIND) :: f_GL     ! ratio of lactate production to glucose consumption (not used in simple metab)
real(REAL_KIND) :: f_IN     ! fraction of total intermediates rate that is N-type (from glutamine only)
real(REAL_KIND) :: f_PPu    ! normal fraction of pyruvate that is processed by the TCA (rest -> lactate)
real(REAL_KIND) :: f_Gu     ! normal fraction of glycosis (r_G) going to make intermediates
real(REAL_KIND) :: f_Glnu   ! normal fraction of glutamine (r_Gn) going to make intermediates
real(REAL_KIND) :: f_ONu    ! normal fraction of ON (r_ONn) going to make intermediates
real(REAL_KIND) :: f_Pu     ! normal fraction of pyruvates (r_P) going to make intermediates
real(REAL_KIND) :: N_GA		! number of ATP molecules generated per glucose molecule in glycosis
real(REAL_KIND) :: N_GI		! number of intermediate molecules generated per glucose molecule in glycosis
real(REAL_KIND) :: N_GP		! number of pyruvate molecules generated per glucose molecule in glycosis
real(REAL_KIND) :: N_GlnA	! number of ATP molecules generated per glutamine molecule
real(REAL_KIND) :: N_GlnI	! number of intermediate molecules generated per glutamine molecule
real(REAL_KIND) :: N_PA		! number of ATP molecules generated per pyruvate molecule in pyruvate oxidation
real(REAL_KIND) :: N_PI		! number of intermediate molecules generated per pyruvate molecule in pyruvate oxidation
real(REAL_KIND) :: N_PO		! number of O2 molecules consumed per pyruvate molecule in pyruvate oxidation
real(REAL_KIND) :: N_GlnO	! number of O2 molecules consumed per glutamine molecule in glutamine oxidation
real(REAL_KIND) :: N_ONO	! number of O2 molecules consumed per ON molecule in ON oxidation
real(REAL_KIND) :: N_ONI	! number of intermediate molecules generated per ON molecule
real(REAL_KIND) :: N_ONA	! number of ATP molecules generated per ON molecule
real(REAL_KIND) :: f_ATPg	! threshold ATP production rate fractions for cell growth
real(REAL_KIND) :: f_ATPs	! threshold ATP production rate fractions for cell survival
real(REAL_KIND) :: f_ATPramp ! multiplying factor for ramp start for reducing r_G, r_P
real(REAL_KIND) :: r_Ag		! threshold ATP production rates for cell growth
real(REAL_KIND) :: r_As		! threshold ATP production rates for cell survival
real(REAL_KIND) :: CO_H		! threshold O2 for Ofactor
real(REAL_KIND) :: CG_H		! threshold glucose for Gfactor
real(REAL_KIND) :: I_rate_max   ! = r_Iu
real(REAL_KIND) :: rIA      ! = r_Iu/r_Au
real(REAL_KIND) :: C_Gln_lo, C_Gln_hi, f_rGln_lo, f_rGln_threshold, f_rON_base
real(REAL_KIND) :: Gln_Nshare
real(REAL_KIND) :: k_glutamine_decay

! By cell
type(metabolism_type), target :: phase_metabolic(3)

type(drug_type), allocatable, target :: drug(:)

integer, allocatable :: gaplist(:)
integer :: ngaps, ndivided
integer, parameter :: max_ngaps = 200000

logical :: bdry_changed
type(LQ_type) :: LQ(MAX_CELLTYPES)
character*(2048) :: inputfile
character*(2048) :: outputfile
character*(2048) :: logmsg
character*(1024) :: header
logical :: test_case(4)

TYPE(winsockport) :: awp_0, awp_1
logical :: use_TCP = .true.         ! turned off in para_main()
logical :: use_CPORT1 = .false.
logical :: stopped, clear_to_send
logical :: simulation_start, par_zig_init, initialized
logical :: use_radiation, use_treatment
!logical :: use_growth_suppression = .true.	! see usage in subroutine CellGrowth
logical :: use_extracellular_O2 = .false.
logical :: use_V_dependence
logical :: use_divide_time_distribution
logical :: use_exponential_cycletime
logical :: use_constant_divide_volume
!logical :: use_volume_method
logical :: use_cell_cycle
logical :: use_constant_growthrate = .false. 
logical :: use_new_drugdata = .true.
logical :: randomise_initial_volume = .true.
logical :: synchronise
logical :: is_radiation
!logical :: use_FD = .true.
logical :: use_permute
logical :: use_gaplist = .true.
!logical :: relax
logical :: medium_change_step
logical :: fully_mixed
logical :: use_parallel
logical :: simulate_colony      ! colony growth will be simulated
logical :: colony_simulation    ! colony growth is being simulated
logical :: use_metabolism
logical :: use_glutamine_decay
logical :: noSS = .false.   ! if true, SS solution is not used for C_P
logical :: dbug = .false.
logical :: bdry_debug

logical :: use_events = .true.

real(REAL_KIND) :: ysave(100000),dCreactsave(100000)

integer :: divide_option = DIVIDE_USE_CLEAR_SITE
!integer :: divide_option = DIVIDE_ALWAYS_PUSH
integer :: idbug = 0
integer :: Nbnd
integer :: seed(2)
integer :: kcell_dbug
integer :: kcell_now
integer :: kcell_test = 1

! PEST variables
logical :: use_PEST = .false.
character*(128) :: PEST_parfile, PEST_outputfile

! C_G base concentration (for checking expt. results)
real(REAL_KIND) :: C_G_base = 10.359

! DNA-PK parameters (temporary)
logical :: DRUG_A_inhibiter = .false.   ! set true if drug is in the input file
logical :: use_inhibiter = .false.
real(REAL_KIND) :: C_inhibiter
real(REAL_KIND) :: a_inhibit = 0.83, b_inhibit = 0.5

!integer :: icentral !extracellular variable index corresponding to a central site (NX/2,NY/2,NZ/2)

! Off-lattice parameters, in the input file but unused here
!real(REAL_KIND) :: a_separation
!real(REAL_KIND) :: a_force, b_force, c_force, x0_force, x1_force, kdrag, frandom

!real(REAL_KIND), allocatable :: omp_x(:), omp_y(:), omp_z(:)

!DEC$ ATTRIBUTES DLLEXPORT :: nsteps, DELTA_T, use_PEST, PEST_outputfile, simulate_colony

contains

!-----------------------------------------------------------------------------------------
! WTIME returns a reading of the wall clock time.
!-----------------------------------------------------------------------------------------
real(DP) function wtime()
!DEC$ ATTRIBUTES DLLEXPORT :: wtime
  integer :: clock_max, clock_rate, clock_reading

  call system_clock ( clock_reading, clock_rate, clock_max )
  wtime = real(clock_reading,kind=DP)/clock_rate
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine logger(msg)
character*(*) :: msg
integer :: error
logical :: logfile_isopen
character*(1) :: LF = char(94)

error = 0
inquire(unit=nflog,OPENED=logfile_isopen)
if (use_TCP) then
    if (awp_0%is_open) then
        call winsock_send(awp_0,trim(msg)//LF,len_trim(msg)+1,error)
    elseif (logfile_isopen) then
        write(nflog,*) trim(msg)
    else
        write(99,*) trim(msg)
    endif
else
	write(*,*) trim(msg)
endif
if (logfile_isopen) then
	write(nflog,'(a)') trim(msg)
	if (error /= 0) then
	    write(nflog,'(a,i4)') 'winsock_send error: ',error
	    close(nflog)
	endif
endif
if (error /= 0) stop
end subroutine

!---------------------------------------------------------------------
! Uniform randomly generates an integer I: n1 <= I <= n2
!---------------------------------------------------------------------
integer function random_int(n1,n2,kpar)
integer :: n1,n2,kpar
integer :: k,R

if (n1 == n2) then
    random_int = n1
elseif (n1 > n2) then
    write(logmsg,*) 'ERROR: random_int: n1 > n2: ',n1,n2
    call logger(logmsg)
    stop
endif
R = par_shr3(kpar)
k = abs(R)
random_int = n1 + mod(k,(n2-n1+1))
end function

!--------------------------------------------------------------------------------
! Make a random choice of an integer from 1 - N on the basis of probabilities in
! the array p(:) (assumed to be normalized).
!--------------------------------------------------------------------------------
integer function random_choice(p,N,kpar)
integer :: N,kpar
real(REAL_KIND) :: p(:)
integer :: k
real(REAL_KIND) :: R, psum

R = par_uni(kpar)
psum = 0
do k = 1,N
    psum = psum + p(k)
    if (R <= psum) then
        random_choice = k
        return
    endif
enddo
write(logmsg,*) 'ERROR: random_choice: ',N,p
call logger(logmsg)
stop
end function

!-----------------------------------------------------------------------------------------
! Returns a unit vector with random 3D direction
!-----------------------------------------------------------------------------------------
subroutine get_random_vector3(v)
real(REAL_KIND) :: v(3)
real(REAL_KIND) :: R1, R2, s, a
integer :: kpar=0

R1 = par_uni(kpar)
R2 = par_uni(kpar)
s = sqrt(R2*(1-R2))
a = 2*PI*R1
v(1) = 2*cos(a)*s
v(2) = 2*sin(a)*s
v(3) = 1 - 2*R2
end subroutine

!--------------------------------------------------------------------------------
! Returns a permutation of the elements of a()
!--------------------------------------------------------------------------------
subroutine permute(a,n,kpar)
integer :: a(*),n,kpar
integer :: i,k,tmp

do i = 1,n
    k = random_int(1,n,kpar)
	tmp = a(i)
	a(i) = a(k)
	a(k) = tmp
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine waste_time(n,dummy)
integer :: k, n
real(REAL_KIND) :: dummy
real(REAL_KIND) :: rsum,R
integer :: kpar=0

rsum = 0
do k = 1,n
    R = par_uni(kpar)
    rsum = rsum + R
enddo
dummy = rsum
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function norm(r)
real(REAL_KIND) :: r(3)

norm = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function norm2(r)
real(REAL_KIND) :: r(3)

norm2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine normalize(r)
real(REAL_KIND) :: r(3)

r = r/norm(r)
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine get_vnorm(v,vnorm)
real(REAL_KIND) :: v(3), vnorm(3)
real(REAL_KIND) :: d

d = dot_product(v,v)
vnorm = v/sqrt(d)
end subroutine

!-----------------------------------------------------------------------------------------
! Changed from kcell argument to cp, because it can be a colony cell, in ccell_list.
! Assumes maximum growth rate.
!-----------------------------------------------------------------------------------------
!subroutine set_divide_volume(kcell,V0)
subroutine set_divide_volume(cp,V0)
!integer :: kcell
real(REAL_KIND) :: V0
real(REAL_KIND) :: Tdiv, fg, Tfixed, Tgrowth0, Tgrowth, rVmax
integer :: ityp, kpar=0
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp

!cp => cell_list(kcell)
ityp = cp%celltype
ccp => cc_parameters(ityp)

rVmax = max_growthrate(ityp)
if (use_exponential_cycletime) then
    if (ccp%Pk_G1 > 0) then
        cp%G1ch_time = rv_exponential(ccp%Pk_G1,kpar)
    else
        cp%G1ch_time = 1.0e10
    endif
    if (ccp%Pk_S > 0) then
        cp%Sch_time = rv_exponential(ccp%Pk_S,kpar)
    else
        cp%Sch_time = 1.0e10
    endif
    if (ccp%Pk_G2 > 0) then
        cp%G2ch_time = rv_exponential(ccp%Pk_G2,kpar)
    else
        cp%G2ch_time = 1.0e10
    endif
    Tgrowth = ccp%T_G1 + ccp%T_S + ccp%T_G2
    Tfixed = cp%G1ch_time + cp%Sch_time + cp%G2ch_time + ccp%T_M
    Tdiv = Tgrowth + Tfixed
    fg = 1
else
    Tgrowth = ccp%T_G1 + ccp%T_S + ccp%T_G2
!    Tfixed = ccp%T_M + ccp%G1_mean_delay + ccp%G2_mean_delay
    Tfixed = ccp%T_M  ! no fixed checkpoint delays with lognormal time step
	Tdiv = DivideTime(ityp)     ! log-normally distributed
	! Try changing this
	!Tgrowth = Tdiv - Tfixed
	!fg = Tgrowth/Tgrowth0
	fg = (Tdiv-Tfixed)/Tgrowth
	! In get_dVdt(), dVdt is divided by cp%fg
	! In cycle, need to multiply T_G1, T_S, T_G2 by cp%fg
!	if (kcell_now == 55) write(*,'(a,i6,2f10.0,f8.4)') 'set_divide_volume: kcell, Tgrowth, Tdiv, fg: ',cp%ID,Tgrowth,Tdiv,fg
endif
cp%divide_volume = V0 + Tgrowth*rVmax   ! this uses the average growth time, to generate the divide volume for all cells
cp%divide_time = Tdiv   ! cycle time, varies with cell
cp%fg = fg
!if (kcell == 1) then
!    write(nflog,'(a,5e12.3)') 'set_divide_volume: V0,Tgrowth,Tdiv,div_vol,rVmax: ',V0,Tgrowth,Tdiv,cp%divide_volume,rVmax
!endif
end subroutine	

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real(REAL_KIND) function DivideTime(ityp)
integer :: ityp
real(REAL_KIND) :: p1, p2
integer :: kpar = 0

dividetime = 0
p1 = divide_dist(ityp)%p1
p2 = divide_dist(ityp)%p2
select case (divide_dist(ityp)%class)
case (NORMAL_DIST)
	DivideTime = rv_normal(p1,p2,kpar)
case (LOGNORMAL_DIST)
	DivideTime = rv_lognormal(p1,p2,kpar)
case (CONSTANT_DIST)
	DivideTime = p1
end select
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real function DivisionTime(ityp)
integer :: ityp
integer :: kpar = 0
real(REAL_KIND), parameter :: rndfraction = 0.2

DivisionTime = rv_lognormal(divide_dist(ityp)%p1,divide_dist(ityp)%p2,kpar)
end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function get_I2Divide(cp) result(I2div)
type(cell_type), pointer :: cp
real(REAL_KIND) :: I2div
integer :: ityp
type(cycle_parameters_type), pointer :: ccp
type(metabolism_type), pointer :: metabolic
	
metabolic => phase_metabolic(1)

ityp = cp%celltype
ccp => cc_parameters(ityp)
!I2div = cp%divide_time*metabolic(ityp)%I_rate_max
I2div = (ccp%T_G1 + ccp%T_S + ccp%T_G2)*I_rate_max
!		*metabolic%I_rate_max
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real(REAL_KIND) function rv_normal(p1,p2,kpar)
integer :: kpar
real(REAL_KIND) :: p1,p2
real(REAL_KIND) :: R

R = par_rnor(kpar)
rv_normal = p1+R*p2
end function

!--------------------------------------------------------------------------------------
! When Y is normal N(p1,p2) then X = exp(Y) is lognormal with
!   median = m = exp(p1)
!   shape  = s = exp(p2)
! Also median = m = mean/(s^2/2)
! kpar = parallel process number
!--------------------------------------------------------------------------------------
real(REAL_KIND) function rv_lognormal(p1,p2,kpar)
integer :: kpar
real(REAL_KIND) :: p1,p2
real(REAL_KIND) :: R,z

R = par_rnor(kpar)
z = p1 + R*p2
rv_lognormal = exp(z)
end function

!--------------------------------------------------------------------------------------
! For testing.
!--------------------------------------------------------------------------------------
real(REAL_KIND) function my_rnor()
real(REAL_KIND) :: sum, R
integer :: k
integer :: kpar=0

sum = 0
do k = 1,12
    R = par_uni(kpar)
    sum = sum + R
enddo
my_rnor = sum - 6.0
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real(REAL_KIND) function rv_exponential(lambda,kpar)
real(REAL_KIND) :: lambda
integer :: kpar
real(REAL_KIND) :: R

R = par_uni(kpar)
rv_exponential = -log(1-R)/lambda
end function

!--------------------------------------------------------------------------------------
! Cumulative probability distribution for a lognormal variate with median m, shape s
! Computes Pr(X < a) where X = exp(Y) and Y is N(p1,p2), p1 = log(m), p2 = log(s)
! Pr(X < a) = Pr(Y < log(a)) = Pr(p1 + p2*R < log(a)) = Pr(R < (log(a)-p1)/p2)
! where R is N(0,1)
!--------------------------------------------------------------------------------------
real(REAL_KIND) function cum_prob_lognormal(a,p1,p2)
real(REAL_KIND) :: a, p1, p2
real(REAL_KIND) :: b, prob

b = (log(a) - p1)/p2
prob = 0.5 + 0.5*erf(b/sqrt(2.0))
cum_prob_lognormal = prob
end function

!-----------------------------------------------------------------------------------------
! Generate a random value for CFSE from a distribution with mean = average
! In the simplest case we can allow a uniform distribution about the average.
! Multiplying factor in the range (1-a, 1+a)
! Better to make it a Gaussian distribution: 
!  = average*(1+s*R)
! where R = N(0,1), s = std deviation
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function generate_CFSE(average)
real(REAL_KIND) :: average, std
integer :: kpar = 0
real(REAL_KIND) :: R

! Uniform distribution
!R = par_uni(kpar)
!generate_CFSE = (1 - a + 2*a*R)*average
! Gaussian distribution
R = par_rnor(kpar)	! N(0,1)
generate_CFSE = (1 + CFSE_std*R)*average
end function

!-----------------------------------------------------------------------------------------
! ityp = cell type
! V0 = cell starting volume (after division) = %volume
! Two approaches:
! 1. Use Vdivide0 and dVdivide to generate a volume
! 2. Use the divide time log-normal distribution 
!    (use_V_dependence = false)
!-----------------------------------------------------------------------------------------
function get_divide_volume(ityp,V0,Tdiv,fg) result(Vdiv)
integer :: ityp
real(REAL_KIND) :: V0, Tdiv, fg
real(REAL_KIND) :: Vdiv, Tfixed, Tgrowth0, Tgrowth, rVmax
real(REAL_KIND) :: b, R
integer :: kpar=0
type(cycle_parameters_type), pointer :: ccp

ccp => cc_parameters(ityp)

rVmax = max_growthrate(ityp)
Tgrowth0 = ccp%T_G1 + ccp%T_S + ccp%T_G2
Tfixed = ccp%T_M + ccp%G1_mean_delay + ccp%G2_mean_delay
if (use_divide_time_distribution) then
	Tdiv = DivideTime(ityp)
	Tgrowth = Tdiv - Tfixed
	fg = Tgrowth/Tgrowth0
	Vdiv = V0 + Tgrowth*rVmax
else
	if (use_constant_divide_volume) then
		Vdiv = Vdivide0
		Tdiv = Tgrowth0 + Tfixed
		fg = 1
	else
		R = par_uni(kpar)
		Vdiv = Vdivide0 + dVdivide*(2*R-1)
		Tgrowth = (Vdiv - V0)/rVmax
		fg = Tgrowth/Tgrowth0
		Tdiv = Tgrowth + Tfixed
	endif
endif
!write(*,'(a,4e11.3)') 'get_divide_volume: rVmax,Vdiv,Tdiv, Tmean: ',rVmax,Vdiv,Tdiv/3600,divide_time_mean(ityp)/3600

end function	

!-----------------------------------------------------------------------------------------
! ityp = cell type
! V0 = cell starting volume (after division) = %volume
! Two approaches:
! 1. Use Vdivide0 and dVdivide to generate a volume
! 2. Use the divide time log-normal distribution 
!    (a) use_V_dependence = true
!    (b) use_V_dependence = false
!-----------------------------------------------------------------------------------------
function get_divide_volume1(ityp,V0,Tdiv) result(Vdiv)
integer :: ityp
real(REAL_KIND) :: V0, Tdiv
real(REAL_KIND) :: Vdiv
real(REAL_KIND) :: Tmean, b, R
integer :: kpar=0

Tmean = divide_time_mean(ityp)
if (use_divide_time_distribution) then
	Tdiv = DivideTime(ityp)
	if (use_V_dependence) then
		b = log(2.0)*(Tdiv/Tmean)
		Vdiv = V0*exp(b)
	else
		Vdiv = V0 + (Vdivide0/2)*(Tdiv/Tmean)
	endif
else
	if (use_constant_divide_volume) then
		Vdiv = Vdivide0
	else
		R = par_uni(kpar)
		Vdiv = Vdivide0 + dVdivide*(2*R-1)
	endif
	Tdiv = Tmean
endif
end function	

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

!-----------------------------------------------------------------------------------------
! Squeeze gaps out of cellist array, adjusting occupancy array.
!-----------------------------------------------------------------------------------------
subroutine squeezer()
integer :: last, kcell, site(3), indx(2), i, j, idc, n, region

!write(*,*) 'squeezer'
!call logger('squeezer')
if (ngaps == 0) return
last = nlist
kcell = 0
n = 0
do
    kcell = kcell+1
    if (cell_list(kcell)%state == DEAD) then    ! a gap
        if (kcell == last) exit
        do
            if (last == 0) then
                write(nflog,*) 'last = 0: kcell: ',kcell
                stop
            endif
            if (cell_list(last)%state == DEAD) then
                last = last-1
                n = n+1
                if (n == ngaps) exit
            else
                exit
            endif
        enddo
        if (n == ngaps) exit
        cell_list(kcell) = cell_list(last)
        last = last-1
        n = n+1
    endif
    if (n == ngaps) exit
enddo
nlist = nlist - ngaps
ngaps = 0
if (dbug) write(nflog,*) 'squeezed: ',n,nlist

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
if (np /= Ncells) then
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
