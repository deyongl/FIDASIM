!!FIDASIM4.0 
!Original version from W. W. Heidbrink 2007
!rewritten in F90 by Benedikt Geiger 2012 (IPP-Garching, ASDEX Upgrade)
!The main routine (fidasim) is at the end of the file!
module application
  implicit none
  !!                      Definition for the kind of the variables: 
  integer , parameter   :: long      = kind(int(1))
  integer , parameter   :: float     = kind(1.e0)
  integer , parameter   :: double    = kind(1.d0) 
  !!                      Indizes of the different components:
  character(100)        :: result_dir
  character(100)        :: root_dir
  integer , parameter   :: nbif_type = 1 ! full energy NBI spectra/density
  integer , parameter   :: nbih_type = 2 ! half energy NBI spectra/density
  integer , parameter   :: nbit_type = 3 ! third energy NBI spectra/density
  integer , parameter   :: halo_type = 4 ! halo spectra/density
  integer , parameter   :: afida_type = 5 ! fida spectra/density 
  integer , parameter   :: pfida_type = 6 ! fida spectra/density 
  integer , parameter   :: ntypes    = 6 ! number of different types of neutrals
  !! random number geneator save variables
  real(double)           :: ran_am
  integer                :: ran_ix=-1,ran_iy=-1,ran_k
  integer, parameter     :: ran_IA=16807,ran_IM=2147483647
  integer, parameter     :: ran_IQ=127773,ran_IR=2836
  !! eigenvalue decomposition values
  real(double),parameter:: ONE=1.d0,TWO=2.d0,ZERO=0.d0
  real(double),parameter:: XMACH_EPS=2.22d-16,MAXIT=50
  !!                      Physical units:
  real(double),parameter:: mass_u    = 1.6605402d-27  ! [kg]
  real(double),parameter:: e0        = 1.60217733d-19 ! [C]
  real(double),parameter:: pi        = 3.14159265358979323846264d0
  real(double),parameter:: c0        = 2.99792458d+08 !! [m/s]
  real(double),parameter:: h_planck  = 4.135667516d-15 !![eV/s]
  real(double),parameter:: lambda0   = 6561.d0        !!D-alpha [A]
  real(double),parameter:: v_to_E    = mass_u/(2.*e0*1.d3)*1.d-4 !!conversion cm^2/s^2 to keV
  !! ---- Stark splitting, wavelength and intenisty of all 15 lines ---- !!
  integer,parameter::n_stark = 15
  real(double),parameter,dimension(n_stark):: stark_wavel = &
       (/ -2.20200d-06,-1.65200d-06,-1.37700d-06,-1.10200d-06 &
       ,  -8.26400d-07,-5.51000d-07,-2.75600d-07, 0.00000d0 &
       ,   2.75700d-07, 5.51500d-07, 8.27400d-07, 1.10300d-06 &
       ,   1.38000d-06, 1.65600d-06, 2.20900d-06/)
  real(double),parameter,dimension(n_stark)::stark_intens= &
       (/   1.d0,   18.d0,   16.d0, 1681.d0, 2304.d0 &
       ,  729.d0, 1936.d0, 5490.d0, 1936.d0,  729.d0 &
       , 2304.d0, 1681.d0,   16.d0,   18.d0,    1.d0/)
  integer,parameter,dimension(n_stark)::stark_pi= &
       (/1,0,0,1,1,1,0,0,0,1,1,1,0,0,1/)
  integer,parameter,dimension(n_stark)::     stark_sigma=1 - stark_pi

  integer, dimension(1) :: minpos  !! dummy array to determine minloc
  !!Numerical Settings
  integer          :: nlevs             !!nr of quantum states  
  integer,parameter:: npitch_birth=100    !!nr pitches in birth profile
  real(double),parameter :: nr_halo_neutrate=200. !! to average halo neut-rate 
  real(double) :: colrad_threshold=1.d5 !! to speed up simulation!
  !! 1.d6 is a very low density rate per marker (typially ~1.d12)!
  !! it should not affect the calculation
  integer :: nbi_outside=0 !! counter for NBI markers that do not
  !! enter the simulation grid
  !! save velocity vectors of individual species
  integer,parameter :: nvelocities=100 !! n velocity vectors
  integer, dimension(nvelocities/2) :: every_second_index
  integer, dimension(nvelocities/2) :: first_half_indices 
 type prof_type
     !! kinetic profiles
     integer::nrho
     real(double) :: drho
     real(double),dimension(:)  ,allocatable :: rho  !! rho grid of profiles
     real(double),dimension(:)  ,allocatable :: te   !! Electron temperature
     real(double),dimension(:)  ,allocatable :: ti   !! Ion temperature
     real(double),dimension(:)  ,allocatable :: dene !! electron density
     real(double),dimension(:)  ,allocatable :: denp !! D-ions density
     real(double),dimension(:)  ,allocatable :: deni !! Carbon density
     real(double),dimension(:)  ,allocatable :: zeff !! zeff
     real(double),dimension(:)  ,allocatable :: vtor !! Plasma rotation [rad/s]
     real(double),dimension(:)  ,allocatable :: bckgrnd_dens !! [1/cm^3]
  end type prof_type
 type grid_type
     integer(long)  :: Nx,Ny,Nz,Nr    !! Nr. of cells in x direction
     real(double)   :: dx,dy,dz,dr
     real(double)   :: dl
     real(double)   :: xmin,xmax,ymin,ymax
     real(double)   :: zmin,zmax,rmin,rmax
     real(double)   :: dv    !! volume of cells
     integer(long)  :: ntrack!! Maximum Nr. of cells for tracking
     real(double), dimension(:), allocatable     :: xx,xxc
     real(double), dimension(:), allocatable     :: yy,yyc 
     real(double), dimension(:), allocatable     :: zz,zzc
     real(double), dimension(:), allocatable     :: rr,rrc
     real(double), dimension(:,:), allocatable   :: rho !! Normalized flux
     real(double), dimension(:,:,:), allocatable :: Brzt,E
  end type grid_type

  type nbi_type 
     integer(long)  :: isource           !! number of the NBI
     real(double)   :: dv                !! half width in y direction
     real(double)   :: dw                !! half width in z direction
     real(double)   :: focy              !! focal lenght in y direction
     real(double)   :: focz              !! focal lenght in z direction
     real(double), dimension(3)  :: divay!! divergence in y direction
     real(double), dimension(3)  :: divaz!! divergence in z direction
     real(double), dimension(3)   :: species_mix
     real(double), dimension(3)   :: xyz_pos !! position of source
     real(double)                 :: einj    !! NBI voltage  [kV]
     real(double)                 :: pinj    !! NBI power    [MW]
     real(double)                 :: vinj    !! NBI velocity [cm/s]
     real(double), dimension(3,3) :: Arot    !! Rotation matrizes of NBI
     real(double), dimension(3,3) :: Brot  
     real(double), dimension(3,3) :: Crot 
  end type nbi_type

  type tables_type
     !! energy array
     real(double)                                  :: d_eb_qp
     real(double)                                  :: d_eb_qi
     real(double)                                  :: d_eb_qe
     real(double)                                  :: d_eb_neut
     real(double)                                  :: d_eb_cx
     integer(long)                                 :: nr_eb_qp
     integer(long)                                 :: nr_eb_qi
     integer(long)                                 :: nr_eb_qe
     integer(long)                                 :: nr_eb_neut
     integer(long)                                 :: nr_eb_cx
     !! temperature array
     real(double)                                  :: d_ti_qp
     real(double)                                  :: d_ti_qi
     real(double)                                  :: d_te_qe
     real(double)                                  :: d_ti_neut
     real(double)                                  :: d_ti_cx

     integer(long)                                 :: nr_ti_qp
     integer(long)                                 :: nr_ti_qi
     integer(long)                                 :: nr_te_qe
     integer(long)                                 :: nr_ti_neut
     integer(long)                                 :: nr_ti_cx
     !! TABLES
     real(double), dimension(:,:,:,:) , allocatable :: qp 
     real(double), dimension(:,:,:,:) , allocatable :: qi
     real(double), dimension(:,:,:,:) , allocatable :: qe
     real(double), dimension(:,:,:,:) , allocatable :: cx
     real(double), dimension(:,:,:)   , allocatable :: neut
     real(double), dimension(:,:)     , allocatable :: einstein
  end type tables_type

  type fbm_type
     integer(long) :: nr
     real(double)  :: rmin
     real(double)  :: dr   
     integer(long) :: nz
     real(double)  :: zmin
     real(double)  :: dz   
     !! energy
     integer(long) :: nenergy
     real(double)  :: emax
     real(double)  :: emin
     real(double)  :: dE
     real(double)  :: eran    
     !!pitch
     integer(long) :: npitch
     real(double)  :: pmax
     real(double)  :: pmin
     real(double)  :: dP
     real(double)  :: pran
     real(double),dimension(:)      ,allocatable :: energy  !! Energy array
     real(double),dimension(:)      ,allocatable :: pitch   !! Pitch array
     real(double),dimension(:)      ,allocatable :: rgrid
     real(double),dimension(:)      ,allocatable :: zgrid   
     real(double),dimension(:,:,:,:),allocatable :: fbm !! fast-ion distribution
     real(double),dimension(:,:)    ,allocatable :: denf!! fast-ion density
  end type fbm_type

  type diag_type
     real(double),dimension(:,:),  allocatable :: xyzlos
     real(double),dimension(:,:),  allocatable :: xyzhead
     real(double),dimension(:),    allocatable :: headsize
     real(double),dimension(:),    allocatable :: opening_angle
     real(double),dimension(:),    allocatable :: sigma_pi
     real(double),dimension(:),    allocatable :: instfu
     character(10),dimension(:),   allocatable :: los_name
     integer(long) :: nchan
     integer(long) :: nlambda
     real(double)  :: dlambda
     real(double)  :: lambdamin
     real(double)  :: lambdamax
     real(double), dimension(:,:,:,:) ,allocatable :: los_wght!! Weights of LOS
  end type diag_type

  type npa_type
     real(double), dimension(:,:)     ,allocatable  :: v    !! velocity array
     real(double), dimension(:,:)     ,allocatable  :: ipos !! initial position arra
     real(double), dimension(:,:)     ,allocatable  :: fpos !! final position array
     real(double), dimension(:)       ,allocatable  :: kind !! kind
     real(double), dimension(:)       ,allocatable  :: wght !! weight
     real(double), dimension(:)       ,allocatable  :: size !! active area of detect
     integer(long)               :: counter
     real(double)                :: npa_loop
     real(double)                :: opening_angle
     real(double), dimension(3)  :: los, los_norm
     real(double)                :: dlos
  end type npa_type

  type result_type
     real(double),dimension(:,:,:,:,:),allocatable:: neut_dens! Density 
     real(double),dimension(:,:,:,:)    ,allocatable:: spectra
     real(double),dimension(:,:,:,:,:),allocatable:: birth_dens!Deposition prof
     real(double), dimension(:,:,:,:,:,:), allocatable :: velocity_vectors
     integer,     dimension(:,:,:,:)     , allocatable :: velocity_counter!!
     integer,     dimension(:,:,:,:)     , allocatable :: vel_vec_red_lev!!
  end type result_type

  type inputs_type
     integer(long) :: shot_number
     real(double)  :: time
     character(15) :: runid
     character(4)  :: diag 
     !! Monte Carlo Settings
     integer(long) :: nr_fida
     integer(long) :: nr_ndmc
     integer(long) :: nr_dcx
     integer(long) :: nr_halo 
     integer(long) :: nr_npa
     !! general settings
     integer(long) :: simfa !! simulate fast ion switch
     integer(long) :: calc_spec
     integer(long) :: load_neutrals
     integer(long) :: background_density
     integer(long) :: npa 
     integer(long) :: calc_wght
     integer(long) :: calc_birth
     !! Plasma parameters
     integer(long) :: impurity_charge
     real(double)  :: btipsign 
     real(double)  :: ai   !! atomic mass of plasma ions
     real(double)  :: ab   !! atomic mass of beam neutrals
     !! Settings for weight function calculation
     integer(long) :: wght_gyro_resolved
     integer(long) :: nr_energy_wght
     integer(long) :: nr_pitch_wght
     integer(long) :: nr_gyro_wght  
     integer(long) :: ichan_wght
     real(double)  :: emax_wght
     real(double)  :: dwav_wght
     real(double)  :: wavel_start_wght
     real(double)  :: wavel_end_wght
  end type inputs_type
  !! definition of the structures:
  type(nbi_type)     :: nbi      
  type(grid_type)    :: grid      
  type(tables_type)  :: tables    
  type(fbm_type)     :: fbm
  type(prof_type)    :: prof
  type(npa_type)     :: npa
  type(diag_type)    :: diag
  type(inputs_type):: inputs
  type(result_type)  :: result
  !! routines:
  public :: ndmc      
  public :: halo
  public :: fida
  !! structures:
  public :: nbi         
  public :: grid         
  public :: tables      
  public :: fbm     
  public :: inputs  
  public :: npa  
  !! routines to read inputs
  public :: check
  public :: read_namelist
  public :: read_diag
  public :: read_grid
  public :: read_profiles
  public :: read_tables
  public :: read_fbm  
  public :: read_neutrals
  public :: write_neutrals 
  public :: write_birth_profile
  public :: write_fida_spectra
  public :: write_nbi_halo_spectra
contains  
  !****************************************************************************
  subroutine check(stat)
    use netcdf
    integer, intent ( in) :: stat
    if(stat /= nf90_noerr) then
       print *, trim(nf90_strerror(stat))
       stop "Stopped: Failed to Write/Read netCDF file"
    end if
  end subroutine check
  !****************************************************************************
  subroutine read_namelist
    character(120)   :: filename
    integer :: i,j,k
    print*,'---- loading namelist -----' 
    filename=trim(adjustl(result_dir))//"/namelist.dat"
    open(66,file=filename)
    read(66,*) !# FIDASIM input file created...
    read(66,"(A100)") root_dir
    read(66,*) !
    read(66,*) inputs%shot_number
    read(66,*) inputs%time
    read(66,*) inputs%runid
    read(66,*) inputs%diag
    read(66,*) inputs%calc_birth
    read(66,*) inputs%simfa
    read(66,*) inputs%calc_spec
    read(66,*) !# this was nofida...
    read(66,*) inputs%npa
    read(66,*) inputs%load_neutrals
    read(66,*) inputs%background_density
    read(66,*) inputs%calc_wght
    read(66,*) !# weight function inputs
    read(66,*) inputs%wght_gyro_resolved
    read(66,*) inputs%nr_energy_wght
    read(66,*) inputs%nr_pitch_wght
    read(66,*) inputs%nr_gyro_wght 
    read(66,*) inputs%ichan_wght
    read(66,*) inputs%emax_wght
    read(66,*) inputs%dwav_wght
    read(66,*) inputs%wavel_start_wght
    read(66,*) inputs%wavel_end_wght
    read(66,*) !# Monte Carlo settings:
    read(66,*) inputs%nr_fida
    read(66,*) inputs%nr_ndmc   
    read(66,*) inputs%nr_halo   
    inputs%nr_dcx = inputs%nr_halo
    read(66,*) inputs%impurity_charge 
    read(66,*) !# Location of transp cdf file:   
    read(66,*) !cdf-file    
    read(66,*) !# discharge parameters:         
    read(66,*) inputs%btipsign
    read(66,*) inputs%ab       
    read(66,*) inputs%ai    
    read(66,*) !# wavelength grid:    
    read(66,*) diag%nlambda 
    read(66,*) diag%lambdamin
    read(66,*) diag%lambdamax
    read(66,*) !# rotate
    read(66,*) !# rotate value 
    read(66,*) !# background density parametrisation
    read(66,*) !# 1
    read(66,*) !# 2
    read(66,*) !# 3
    read(66,*) !# 4
    read(66,*) !# 5
    read(66,*) !# 6
    diag%dlambda=(diag%lambdamax-diag%lambdamin)/diag%nlambda
    close(66)
    print*,                         'Shot  :',inputs%shot_number
    print*,                         'Time:',int(inputs%time*1.d3),'ms'
  end subroutine read_namelist
  !****************************************************************************


    
  !****************************************************************************
  subroutine read_grid
    character(120)   :: filename
    integer :: i,j,k
    print*,'---- reading grid -----' 
    filename=trim(adjustl(result_dir))//"/grid.bin"
    open(66,file=filename,access='stream')
    read(66) !# simulation grid: 
    read(66) grid%Nx 
    read(66) grid%dx
    read(66) grid%xmin
    read(66) grid%xmax
    read(66) grid%Ny 
    read(66) grid%dy
    read(66) grid%ymin
    read(66) grid%ymax
    read(66) grid%Nz
    read(66) grid%dz
    read(66) grid%zmin
    read(66) grid%zmax
    read(66) grid%Nr
    read(66) grid%dr
    read(66) grid%rmin
    read(66) grid%rmax
    read(66) grid%dv
    read(66) grid%ntrack
    read(66) grid%dl
    allocate(grid%xx(grid%Nx),grid%yy(grid%Ny) &
         ,grid%zz(grid%Nz),grid%rr(grid%Nr))
    allocate(grid%xxc(grid%Nx),grid%yyc(grid%Ny) &
         ,grid%zzc(grid%Nz),grid%rrc(grid%Nr))
    read(66) grid%xx
    read(66) grid%yy
    read(66) grid%zz
    read(66) grid%rr
    read(66) grid%xxc
    read(66) grid%yyc
    read(66) grid%zzc
    read(66) grid%rrc
    close(66)
    allocate(grid%Brzt(grid%Nr,grid%Nz,3)) 
    allocate(grid%E(grid%Nr,grid%Nz,3)) 
    allocate(grid%rho(grid%Nr,grid%Nz)) 
    grid%Brzt=0.d0 ; grid%E=0.d0 ; grid%rho=0.d0
  end subroutine read_grid
  !****************************************************************************


  subroutine read_nbi
    character(120)   :: filename
    integer :: i,j,k
    print*,'---- reading NBI data -----' 
    filename=trim(adjustl(result_dir))//"/nbi.dat"
    open(66,file=filename)
    read(66,*) !# Neutral beam injection
    read(66,*) nbi%dv
    read(66,*) nbi%dw  
    read(66,*) nbi%isource
    read(66,*) nbi%divay(1) 
    read(66,*) nbi%divay(2) 
    read(66,*) nbi%divay(3) 
    read(66,*) nbi%divaz(1)  
    read(66,*) nbi%divaz(2)  
    read(66,*) nbi%divaz(3)  
    read(66,*) nbi%focy  
    read(66,*) nbi%focz   
    read(66,*) nbi%einj 
    read(66,*) nbi%pinj 
    read(66,*) !# Species-mix (Particles):
    read(66,*) nbi%species_mix(1)
    read(66,*) nbi%species_mix(2)
    read(66,*) nbi%species_mix(3)
    read(66,*) !#position of NBI source in xyz coords:
    read(66,*) nbi%xyz_pos(1)
    read(66,*) nbi%xyz_pos(2)
    read(66,*) nbi%xyz_pos(3)
    read(66,*) !# 3 rotation matrizes 3x3
    do j=1,3 
       do k=1,3
          read(66,*) nbi%Arot(j,k)
          read(66,*) nbi%Brot(j,k)
          read(66,*) nbi%Crot(j,k)
       enddo
    enddo
    close(66) 
    nbi%vinj=sqrt(2.d0*nbi%einj*1.d3 &
         *e0/(inputs%ab*mass_u))*1.d2 !! [cm/s]
    print*, 'NBI #',nbi%isource+1
    print*,'NBI power   :', real(nbi%pinj,float)
    print*,'NBI voltage :', real(nbi%einj,float)

!!$    print*,nbi%dv,nbi%dw
!!$    print*,nbi%divay,nbi%divaz
!!$    print*,nbi%focy,nbi%focz
!!$    print*,nbi%einj,nbi%pinj
!!$    print*,nbi%species_mix
!!$    print*, nbi%Crot(:,:)
!!$    print*, nbi%Brot(:,:)
!!$    print*, nbi%Arot(:,:)
!!$    print*, nbi%vinj
  end subroutine read_nbi

  subroutine read_beam
    use netcdf
    character(120) :: filename
    integer :: ncid,ab_varid,divy_varid,divz_varid,focy_varid,focz_varid,bn_varid
    integer :: einj_varid,pinj_varid,sm_varid,xyzsrc_varid,bwra_varid,bwza_varid
    integer :: arot_varid,brot_varid,crot_varid
    filename=trim(adjustl(result_dir))//"/"//trim(adjustl(inputs%runid))//"_inputs.cdf"
    print*,'---- loading beam ----'
    !!OPEN netCDF file
    call check( nf90_open(filename, nf90_nowrite, ncid) )
    !!Get the varids
    call check( nf90_inq_varid(ncid, "beam", bn_varid) )
    call check( nf90_inq_varid(ncid, "ab", ab_varid) )
    call check( nf90_inq_varid(ncid, "divy", divy_varid) )
    call check( nf90_inq_varid(ncid, "divz", divz_varid) )
    call check( nf90_inq_varid(ncid, "focy", focy_varid) )
    call check( nf90_inq_varid(ncid, "focz", focz_varid) )
    call check( nf90_inq_varid(ncid, "bmwidra", bwra_varid) )
    call check( nf90_inq_varid(ncid, "bmwidza", bwza_varid) )
    call check( nf90_inq_varid(ncid, "einj", einj_varid) )
    call check( nf90_inq_varid(ncid, "pinj", pinj_varid) )
    call check( nf90_inq_varid(ncid, "species_mix", sm_varid) )
    call check( nf90_inq_varid(ncid, "xyz_src", xyzsrc_varid) )
    call check( nf90_inq_varid(ncid, "Arot", arot_varid) )
    call check( nf90_inq_varid(ncid, "Brot", brot_varid) )
    call check( nf90_inq_varid(ncid, "Crot", crot_varid) )
    !!Read the variables
    call check( nf90_get_var(ncid, bn_varid, nbi%isource) )
    call check( nf90_get_var(ncid, bwra_varid, nbi%dv) )
    call check( nf90_get_var(ncid, bwza_varid, nbi%dw) )
    call check( nf90_get_var(ncid, divy_varid, nbi%divay(:)) )
    call check( nf90_get_var(ncid, divz_varid, nbi%divaz(:)) )
    call check( nf90_get_var(ncid, focy_varid, nbi%focy) )
    call check( nf90_get_var(ncid, focz_varid, nbi%focz) )
    call check( nf90_get_var(ncid, einj_varid, nbi%einj) )
    call check( nf90_get_var(ncid, pinj_varid, nbi%pinj) )
    call check( nf90_get_var(ncid, sm_varid, nbi%species_mix(:)) )
    call check( nf90_get_var(ncid, xyzsrc_varid, nbi%xyz_pos(:)) )
    call check( nf90_get_var(ncid, arot_varid, nbi%Arot(:,:)) )
    call check( nf90_get_var(ncid, brot_varid, nbi%Brot(:,:)) )
    call check( nf90_get_var(ncid, crot_varid, nbi%Crot(:,:)) )


    
    !!CLOSE netCDF FILE
    call check( nf90_close(ncid) )
    nbi%vinj=sqrt(2.d0*nbi%einj*1.d3 &
         *e0/(inputs%ab*mass_u))*1.d2 !! [cm/s]
    print*,'NBI #',nbi%isource+1
    print*,'NBI power :', real(nbi%pinj,float)
    print*,'NBI voltage :', real(nbi%einj,float)
  end subroutine read_beam



  subroutine read_diag
    character(120)  :: filename
    integer         :: ichan
    filename=trim(adjustl(result_dir))//"/diag.bin"
    print*,'---- loading diagnostic ----'
    open(66,file=filename,access='stream')
    read(66)diag%nchan
    allocate(diag%xyzhead(diag%nchan,3))
    allocate(diag%xyzlos(diag%nchan,3))
    allocate(diag%headsize(diag%nchan))
    allocate(diag%opening_angle(diag%nchan))
    allocate(diag%sigma_pi(diag%nchan))
    allocate(diag%instfu(diag%nchan))
    allocate(diag%los_name(diag%nchan))
    allocate(diag%los_wght(grid%Nx,grid%Ny,grid%Nz,diag%nchan))
    read(66) diag%xyzhead(:,:)
    read(66) diag%xyzlos(:,:)
    read(66) diag%headsize(:)
    read(66) diag%opening_angle(:)
    print*,'test2'
    read(66) diag%sigma_pi(:)
    read(66) diag%instfu(:)
    read(66) diag%los_name(:)
    print*,'test3'
    read(66) diag%los_wght(:,:,:,:)
    print*,'test4'
    close(66)
    if (inputs%npa.eq.0) then 
       npa%npa_loop=1.
    else
       allocate(npa%size(diag%nchan))
       npa%size=diag%headsize
       npa%npa_loop=100.
       npa%opening_angle=diag%opening_angle(1)
       npa%los=diag%xyzhead(1,:)-diag%xyzlos(1,:)
       npa%dlos=sqrt(dot_product(npa%los,npa%los))
       npa%los_norm=npa%los/npa%dlos
    endif
  end subroutine read_diag

  !***************************************************************************!
  subroutine read_profiles
    character(120):: filename
    filename=trim(adjustl(result_dir))//"/profiles.bin"
    print*,'---- reading profiles ---------------'
    open(66,file=filename,access='stream')
    read(66)prof%nrho
    allocate(prof%rho(prof%nrho), prof%te(prof%nrho), prof%ti(prof%nrho), &
         prof%dene(prof%nrho), prof%denp(prof%nrho), prof%deni(prof%nrho), &
         prof%vtor(prof%nrho), prof%zeff(prof%nrho), &
         prof%bckgrnd_dens(prof%nrho))
    read(66)prof%drho
    read(66)prof%rho
    read(66) prof%te
    read(66) prof%ti
    read(66) prof%dene
    read(66) prof%denp
    read(66) prof%deni
    read(66) prof%vtor
    read(66) prof%zeff
    read(66) prof%bckgrnd_dens
    close(66)
  end subroutine read_profiles

  subroutine read_field
    character(120):: filename
    filename=trim(adjustl(result_dir))//"/field.bin"
    print*,'---- reading magnetic field ------ '
    open(66,file=filename,access='stream')
    read(66) grid%Brzt(:,:,:)
    read(66) grid%E(:,:,:)
    read(66) grid%rho(:,:)
    close(66)
  end subroutine read_field
  !****************************************************************************
  subroutine read_tables
    character(120)  :: filename
    integer         :: n,m !! initial/final state
    integer         :: ie,iti !! energy/ti index
    integer(long)   :: nlev
   !-------------------ELECTRON EXCITATION/IONIZATION TABLE--------
    filename=trim(adjustl(root_dir))//"TABLES/qetable.bin"
    open(66,file=filename,access='stream')
    read(66) tables%nr_te_qe
    read(66) tables%d_te_qe
    read(66) tables%nr_eb_qe
    read(66) tables%d_eb_qe
    read(66) nlevs
    print*, 'levels in colrad: ', nlevs
    allocate(tables%qe(nlevs+1,nlevs,tables%nr_eb_qe,tables%nr_te_qe))
    read(66) tables%qe(:,:,:,:)
    close(66)

    !-------------------Deuterium EXCITATION/IONIZATION/CX TABLE------
    filename=trim(adjustl(root_dir))//"TABLES/qptable.bin"
    open(66,file=filename,access='stream')
    read(66) tables%nr_ti_qp
    read(66) tables%d_ti_qp
    read(66) tables%nr_eb_qp
    read(66) tables%d_eb_qp
    read(66) nlev
    if(nlev.ne.nlevs)then
       print*, tables%nr_ti_qp,tables%d_ti_qp, tables%nr_eb_qp,tables%d_eb_qp
       print*, nlev,nlevs
       stop 'stop at "read qptable"'
    endif
    allocate(tables%qp(nlevs+1,nlevs,tables%nr_eb_qp,tables%nr_ti_qp))
    read(66) tables%qp(:,:,:,:)
    close(66)


   !-------------------Deuterium CHARGE EXCHANGE BEAM-THERM-RATE ------
    filename=trim(adjustl(root_dir))//"TABLES/cxtable.bin"
    open(66,file=filename,access='stream')
    read(66) tables%nr_ti_cx
    read(66) tables%d_ti_cx
    read(66) tables%nr_eb_cx
    read(66) tables%d_eb_cx
    read(66) nlev
    if(nlev.ne.nlevs)then
       print*, tables%nr_ti_cx,tables%d_ti_cx, tables%nr_eb_cx,tables%d_eb_cx
       print*, nlev,nlevs
       stop 'stop at "read cxtable"'
    endif
    allocate(tables%cx(nlevs+1,nlevs,tables%nr_eb_cx,tables%nr_ti_cx))

    read(66) tables%cx(:,:,:,:)
    close(66)


    !------------------ m-resolved CHARGE EXCHANGE cross-sections  ---
    ! H(+) + H(n) --> H(m) + H(+)
    ! energy in keV/amu
    filename=trim(adjustl(root_dir))//"TABLES/neuttable.bin"
    open(66,file=filename,access='stream')
    read(66) tables%nr_eb_neut
    read(66) tables%d_eb_neut
    read(66) nlev
    if(nlev.ne.nlevs) stop 'stop at read tables (neut_rates)'
    allocate(tables%neut(nlevs,nlevs,tables%nr_eb_neut))
    read(66) tables%neut(:,:,:)
    close(66)

    !-------------------Impurity EXCITATION/IONIZATION/CX TABLE--------
   if(inputs%impurity_charge.lt.5.or.inputs%impurity_charge.gt.7) &
         stop 'wrong impurity charge!'
    if(inputs%impurity_charge.eq.5) &
         filename=trim(adjustl(root_dir))//"TABLES/qbtable.bin"
    if(inputs%impurity_charge.eq.6) &
         filename=trim(adjustl(root_dir))//"TABLES/qctable.bin" 
    if(inputs%impurity_charge.eq.7) &
         filename=trim(adjustl(root_dir))//"TABLES/qntable.bin"
    open(66,file=filename,access='stream')
    read(66) tables%nr_ti_qi
    read(66) tables%d_ti_qi
    read(66) tables%nr_eb_qi
    read(66) tables%d_eb_qi
    read(66) nlev
    if(nlev.ne.nlevs) stop 'stop at "read atomics qptable"'
    allocate(tables%qi(nlevs+1,nlevs,tables%nr_eb_qi,tables%nr_ti_qi))
    read(66) tables%qi(:,:,:,:)
    close(66)

    !-------------------EINSTEIN COEFFICIENTS ----------------------
    filename='TABLES/einstein.dat'
    filename=trim(adjustl(root_dir))//"TABLES/einstein.dat"
    open(66,file=filename)
    read(66,*)! 
    read(66,*)!
    read(66,*) nlev
    if(nlev.ne.nlevs) stop 'stop at "read atomics"'
    allocate(tables%einstein(nlevs,nlevs))
    do n=1,nlevs !! lower level
       do m=1,nlevs  !! upper level
          read(66,*)tables%einstein(m,n) 
       enddo
    enddo
    close(66)

  end subroutine read_tables
  !****************************************************************************
  !-------------------FASTION DISTRIBUTION FUNCTION ----------------------
  subroutine read_fbm
    character(120)  :: filename
    character(17)   :: cdf_file
    real(double)    :: cdf_time  
    logical   :: ok
    filename=trim(adjustl(result_dir))//"/fast-ion_distribution.bin"
    inquire(file=filename,exist=ok)
    if(ok)then
       print*,'---- reading fast ion distribution function ----'
       open(66,file=filename,access='stream')
       read(66) cdf_file
       read(66) cdf_time
       read(66) fbm%nr
       read(66) fbm%rmin
       read(66) fbm%dr
       allocate(fbm%rgrid(fbm%nr))
       read(66) fbm%rgrid
       read(66) fbm%nz
       read(66) fbm%zmin
       read(66) fbm%dz
       allocate(fbm%zgrid(fbm%nz))
       read(66) fbm%zgrid
       !!denf
       allocate(fbm%denf(fbm%nr,fbm%nz))
       read(66) fbm%denf
       !! energy
       read(66) fbm%nenergy
       read(66) fbm%emax
       read(66) fbm%emin
       read(66) fbm%dE
       allocate(fbm%energy(fbm%nenergy))
       read(66)fbm%energy(:)
       !! pitch
       read(66) fbm%npitch
       read(66) fbm%pmax
       read(66) fbm%pmin
       read(66) fbm%dP
       allocate(fbm%pitch(fbm%npitch))
       read(66) fbm%pitch(:)
       allocate(fbm%fbm(fbm%nenergy,fbm%npitch,fbm%nr,fbm%nz))
       read(66) fbm%fbm(:,:,:,:)
       close(66)
       fbm%eran   = fbm%emax-fbm%emin
       fbm%pran   = fbm%pmax-fbm%pmin
    else
       !!no fast-ion distribution function
       inputs%nr_fida=0 !! no fast-ions to be simulated
    endif
  end subroutine read_fbm




  subroutine write_birth_profile
    integer         :: i,j,k,p 
    character(120)  :: filename
    filename=trim(adjustl(result_dir))//"/birth.bin"    
    open (66, file =filename,access='stream')
    write(66)real(inputs%shot_number,float )
    write(66)real(inputs%time)
    write(66)real(grid%Nx,float)
    write(66)real(grid%Ny,float)
    write(66)real(grid%Nz,float) 
    write(66)real(npitch_birth,float) 
    write(66)real(result%birth_dens(:,:,:,:,:),float)
    close (66)
    print*, 'birth profile written to:      ',filename
  end subroutine write_birth_profile

  subroutine write_neutrals
    integer         :: i,j,k
    character(120)  :: filename
    real(double), dimension(3)                  :: ipos      !! start position

    filename=trim(adjustl(result_dir))//"/neutrals.bin"    
    open (66, file =filename,access='stream')
    write(66)real(inputs%shot_number,float )
    write(66)real(inputs%time)
    write(66)real(grid%Nx,float)
    write(66)real(grid%Ny,float)
    write(66)real(grid%Nz,float) 
    write(66)real(nlevs  ,float) 
    write(66)real(result%neut_dens(:,:,:,:,nbif_type),float)
    write(66)real(result%neut_dens(:,:,:,:,nbih_type),float)
    write(66)real(result%neut_dens(:,:,:,:,nbit_type),float)
    write(66)real(result%neut_dens(:,:,:,:,halo_type),float)
    close (66)
    !!print*, 'neutral density written to:      ',filename
  end subroutine write_neutrals

  subroutine write_npa
    integer         :: i,j,k
    character(120)  :: filename
    real(double)    :: bckgrnd_dens_tmp
    real(double), dimension(3)                  :: ipos      !! start position
    real(float),  dimension(:,:),   allocatable :: output
    real(float),  dimension(:,:,:), allocatable :: bckgrnd_dens
    allocate(output(npa%counter,3))

    !! correction of fluxes for detector size
    npa%wght(:)=npa%wght(:)/(pi*npa%size(1)**2)
    !! correction of fluxes for gyro-angle restrictions
    npa%wght(:)=npa%wght(:)*(4.*npa%opening_angle/2./pi)

    filename=trim(adjustl(result_dir))//"/npa.bin"      
    open (66, file =filename,access='stream')
    write(66)inputs%shot_number
    write(66)real(inputs%time,float)
    write(66)inputs%nr_npa
    write(66)npa%counter
    output(:,:)=real(npa%ipos(:npa%counter,:),float)
    write(66)output
    output(:,:)=real(npa%fpos(:npa%counter,:),float)
    write(66)output
    output(:,:)=real(npa%v(:npa%counter,:)   ,float)
    write(66)output
    write(66)real(npa%wght(:npa%counter)  ,float)
    write(66)real(npa%kind(:npa%counter)  ,float)
    close (66)
    print*, 'NPA data written to: ',filename
    deallocate(output)
  end subroutine write_npa
  
  subroutine write_nbi_halo_spectra
    integer         :: i,j,k,ichan
    character(120)  :: filename
    real(double), dimension(:)  , allocatable :: lambda_arr
    !! ------------------------ calculate wavelength array ------------------ !!
    allocate(lambda_arr(diag%nlambda))
    do i=1,diag%nlambda 
       lambda_arr(i)=(i-0.5)*diag%dlambda*0.1d0 &
            +diag%lambdamin*0.1d0
    enddo 
    !! convert [Ph/(s*wavel_bin*cm^2*all_directions)] to [Ph/(s*nm*sr*m^2)]!
    result%spectra(:,:,:,:)=result%spectra(:,:,:,:)/(0.1d0*diag%dlambda) &
         /(4.d0*pi)*1.d4
    !! write to file
    filename=trim(adjustl(result_dir))//"/nbi_halo_spectra.bin"
    open (66, file =filename,access='stream')
    write(66)inputs%shot_number
    write(66)real(inputs%time,float)
    write(66)diag%nchan
    write(66)diag%nlambda
    !--- wavelength array: -----
    write(66)real(lambda_arr(:),float)
    !-------- spectra:----------
    write(66)real(result%spectra(:,:,:,nbif_type),float)
    write(66)real(result%spectra(:,:,:,nbit_type),float)
    write(66)real(result%spectra(:,:,:,nbih_type),float)
    write(66)real(result%spectra(:,:,:,halo_type),float)
    close (66)
    !! result arrays 
    deallocate(lambda_arr)
    print*, 'NBI and HALO spectra written to: ', filename
  end subroutine write_nbi_halo_spectra

 subroutine write_fida_spectra
    integer         :: i,j,k,ichan
    character(120)  :: filename
    real(double), dimension(:)  , allocatable :: lambda_arr
    !! ------------------------ calculate wavelength array ------------------ !!
    allocate(lambda_arr(diag%nlambda))
    do i=1,diag%nlambda 
       lambda_arr(i)=(i-0.5)*diag%dlambda*0.1d0 &
            +diag%lambdamin*0.1d0
    enddo 
    !! convert [Ph/(s*wavel_bin*cm^2*all_directions)] to [Ph/(s*nm*sr*m^2)]!
    result%spectra(:,:,:,afida_type)=result%spectra(:,:,:,afida_type) &
         /(0.1d0*diag%dlambda)/(4.d0*pi)*1.d4
    result%spectra(:,:,:,pfida_type)=result%spectra(:,:,:,pfida_type) &
         /(0.1d0*diag%dlambda)/(4.d0*pi)*1.d4
    !! write to file
    filename=trim(adjustl(result_dir))//"/fida_spectra.bin"  
    open (66, file =filename,access='stream')
    write(66)inputs%shot_number
    write(66)real( inputs%time,float)
    write(66)diag%nchan 
    write(66)diag%nlambda
    !--- wavelength array: -----
    write(66)real(lambda_arr(:),float)
    !-------- spectra:----------
    write(66)real(result%spectra(:,:,:,afida_type),float)
    write(66)real(result%spectra(:,:,:,pfida_type),float)
    close (66)
    !! result arrays 
    deallocate(lambda_arr)
    !print*, 'FIDA spectra written to ',filename
  end subroutine write_fida_spectra

  subroutine read_neutrals
    character(120)  :: filename
    real(float)     :: fdum
    real(float),dimension(:,:,:,:),allocatable   :: fdum_arr
    allocate(fdum_arr(grid%nx,grid%ny,grid%nz,nlevs))
    print*,'---- load neutrals  ----' 
    filename=trim(adjustl(result_dir))//"/neutrals.bin" 
    open (66, file =filename,access='stream')
    read(66)fdum
    read(66)fdum
    read(66)fdum
    read(66)fdum
    read(66)fdum
    read(66)fdum
    read(66)fdum_arr ;  result%neut_dens(:,:,:,:,nbif_type)=fdum_arr
    read(66)fdum_arr ;  result%neut_dens(:,:,:,:,nbih_type)=fdum_arr
    read(66)fdum_arr ;  result%neut_dens(:,:,:,:,nbit_type)=fdum_arr
    read(66)fdum_arr ;  result%neut_dens(:,:,:,:,halo_type)=fdum_arr
    deallocate(fdum_arr)
    close (66)
   end subroutine read_neutrals


 
  
  !*****************************************************************************
  !------------random number generator-----------------------------------------
  !*****************************************************************************
  function ran()
    !!uniform random number generator from NUMERICAL RECEPIES
    real(double)           :: ran
    !$OMP CRITICAL(random)
    ran_ix=ieor(ran_ix,ishft(ran_ix,13)) !Marsaglia shift sequence 
    ran_ix=ieor(ran_ix,ishft(ran_ix,-17))
    ran_ix=ieor(ran_ix,ishft(ran_ix,5))
    ran_k=ran_iy/ran_IQ !Park-Miller sequence by Schrage’s method,
    ran_iy=ran_IA*(ran_iy-ran_k*ran_IQ)-ran_IR*ran_k
    if(ran_iy.lt.0) ran_iy=ran_iy+ran_IM
    ran=ran_am*ior(iand(ran_IM,ieor(ran_ix,ran_iy)),1) !Combine the generators
    !$OMP END CRITICAL(random)
  end function ran  
  subroutine randn(randomn)
    !!Box Mueller Method to calculate normal distribution
    real(double),dimension(:),intent(inout):: randomn
    integer                                :: nran, i 
    real(double)                           :: x1,x2,w  
    randomn=0.d0 ;  nran=size(randomn) ; i=1
    do while (i.le.nran)
       w=1.d0
       do while (w.ge.1.d0) 
          x1=2.d0*ran()-1.d0
          x2=2.d0*ran()-1.d0
          w=x1*x1+x2*x2
       enddo
       w=sqrt((-2.d0*log(w))/w)
       randomn(i)=x1*w
       i=i+1
       if(i.gt.nran)exit
       randomn(i)=x2*w
       i=i+1
    enddo
  end subroutine randn 
  subroutine randu(randomu)
    real(double), dimension(:), intent(inout):: randomu
    integer                                  :: i
    randomu=0.d0
    do i=1,size(randomu) 
       randomu(i)=ran()
    enddo
  end subroutine randu  
  

  subroutine rho_interp(pos,rho)
    real(double),dimension(3), intent(in)  :: pos
    real(double),intent(out) :: rho
    integer(long):: ir,iz
    real(double) :: r,z,r0,z0,c00,c10,c01,c11,c0,c1
    rho=2.d0
    r=sqrt(pos(1)**2+pos(2)**2)
    z=pos(3)
    ir=floor((r-grid%rmin)/grid%dr) +1
    iz=floor((z-grid%zmin)/grid%dz) +1
    if((ir.lt.1.).or.(iz.lt.1).or.(ir.gt.grid%nr-1).or.(iz.gt.grid%nz-1))return
    r0=(ir-1)*grid%dr+grid%rmin
    z0=(iz-1)*grid%dz+grid%zmin
    c00 = grid%rho(ir  ,iz  )
    c10 = grid%rho(ir+1,iz  )
    c01 = grid%rho(ir  ,iz+1)
    c11 = grid%rho(ir+1,iz+1)
    !linear interpolation between C00 and C10 to find C0, C01 and C11 to find C1
    c0  = c00 + (c10 - c00) * (r - r0) /grid%dr
    c1  = c01 + (c11 - c01) * (r - r0) /grid%dr
    rho=  c0  + (c1  - c0)  * (z - z0) /grid%dz
  end subroutine rho_interp

  subroutine get_phi(pos,phi)
    real(double), dimension(3),intent(in) :: pos
    real(double),              intent(out):: phi
    phi=atan(pos(2)/pos(1))
    if(pos(1).gt.0.and.pos(2).lt.0)phi=phi+2.*pi
    if(pos(1).lt.0.)phi=phi+pi
  end subroutine get_phi

  subroutine bfield_interp(pos,b_abs,b_norm,a_norm,c_norm,brzt_out)
    real(double),dimension(3) , intent(in) :: pos
    real(double),intent(out)  :: b_abs
    real(double),dimension(3) ,intent(out):: b_norm
    real(double),dimension(3) ,intent(out), optional::a_norm,c_norm,brzt_out
    integer(long):: ir,iz
    real(double):: r,z,r0,z0,phi
    real(double),dimension(3)  ::c00,c10,c01,c11,c0,c1 
    real(double),dimension(3)  ::Brzt , Bxyz  
    b_abs=0.d0
    r=sqrt(pos(1)**2+pos(2)**2)
    z=pos(3)
    ir=floor((r-grid%rmin)/grid%dr) +1
    iz=floor((z-grid%zmin)/grid%dz) +1
    if((ir.lt.1.).or.(iz.lt.1).or.(ir.gt.grid%nr-1).or.(iz.gt.grid%nz-1))return

    r0=(ir-1)*grid%dr+grid%rmin
    z0=(iz-1)*grid%dz+grid%zmin
    !! b_abs
    c00 = grid%Brzt(ir  ,iz  ,:)
    c10 = grid%Brzt(ir+1,iz  ,:)
    c01 = grid%Brzt(ir  ,iz+1,:)
    c11 = grid%Brzt(ir+1,iz+1,:)
    c0  = c00 + (c10 - c00) * (r - r0) /grid%dr
    c1  = c01 + (c11 - c01) * (r - r0) /grid%dr
    Brzt=  c0  + (c1  - c0) * (z - z0) /grid%dz
    if(present(brzt_out))brzt_out=brzt
 

    !! convert to xyz coordinates
    call get_phi(pos,phi)
    bxyz(1) =- cos(pi*0.5d0-phi)*brzt(3) + cos(phi)*brzt(1)
    bxyz(2) =  sin(pi*0.5d0-phi)*brzt(3) + sin(phi)*Brzt(1) 
    bxyz(3) =  brzt(2)

    b_abs =sqrt(dot_product(bxyz,bxyz))
    b_norm= bxyz/b_abs

    if(present(a_norm))then
       if (abs(b_norm(3)).eq.1) then
          a_norm=(/1.d0,0.d0,0.d0/)
          c_norm=(/0.d0,1.d0,0.d0/)
       else 
          if (b_norm(3).eq.0.) then
             a_norm=(/0.d0,0.d0,1.d0/)
             c_norm=(/b_norm(2),-b_norm(1), 0.d0/) &
                  /sqrt(b_norm(1)**2+b_norm(2)**2)
          else
             a_norm=(/b_norm(2),-b_norm(1),0.d0/) &
                  /sqrt(b_norm(1)**2+b_norm(2)**2)
             c_norm=(/ a_norm(2) , &
                  -a_norm(1) , &
                  (a_norm(1)*b_norm(2)-a_norm(2)*b_norm(1))/b_norm(3)/)
             c_norm=c_norm/sqrt(dot_product(c_norm,c_norm))
          endif
       endif
    endif
  end subroutine bfield_interp



  subroutine prof_interp(pos,denp,ti,vtor,te,dene,deni,bckgrnd_dens)
    real(double),dimension(3), intent(in) :: pos
    real(double),              intent(out),optional:: ti,te,dene,denp,deni,vtor,bckgrnd_dens
    integer(long):: irho
    integer(long),dimension(3) :: ac
    real(double) :: rho,c0,c1,c2,arg,v
    call rho_interp(pos,rho)
    irho=floor(rho/prof%drho)+1  
    if(irho.gt.prof%nrho-1)irho=prof%nrho-1
    if(irho.lt.1)irho=1
    c2=(rho - prof%drho*(irho-1)) / prof%drho
    !! Denp:
    c0=prof%denp(irho)
    c1=prof%denp(irho+1)
    denp  = c0 + (c1 - c0) * c2
    !! Ti:
    if(present(ti))then
       c0=prof%ti(irho)
       c1=prof%ti(irho+1)
       ti  = c0 + (c1 - c0) * c2
    !! Vtor:
       c0=prof%vtor(irho)
       c1=prof%vtor(irho+1)
       vtor  = c0 +(c1 - c0) * c2 
       vtor  = vtor * sqrt(pos(1)**2+pos(2)**2) !! conversion from rad/s to cm/s
       if(present(te))then 
          !! Te:
          c0=prof%te(irho)
          c1=prof%te(irho+1)
          te  = c0 + (c1 - c0) * c2
          !! Dene:
          c0=prof%dene(irho)
          c1=prof%dene(irho+1)
          dene  = c0 + (c1 - c0) * c2
          !! Deni:
          c0=prof%deni(irho)
          c1=prof%deni(irho+1)
          deni  = c0 + (c1 - c0) * c2
          if(present(bckgrnd_dens))then
             !! Bckgrnd_dens:
             c0=prof%bckgrnd_dens(irho)
             c1=prof%bckgrnd_dens(irho+1)
             bckgrnd_dens  = c0 + (c1 - c0) * c2
          endif
       endif
    endif
  end subroutine prof_interp


  subroutine calc_vrot(pos,vtor,vrot)
    real(double), dimension(3),intent(in) :: pos
    real(double),              intent(in) :: vtor
    real(double), dimension(3),intent(out):: vrot
    real(double) :: phi,arg
    call get_phi(pos,phi)
    arg=pi*0.5d0-phi
    vrot(1) = - cos(arg)*vtor
    vrot(2) =   sin(arg)*vtor
    vrot(3) =   0.d0
  end subroutine calc_vrot

  subroutine denf_interp(pos,denf)
    !! bilinear interpolation of the effective rate coefficient tables !!
    real(double),dimension(3),              intent(in) :: pos
    real(double) ,intent(out) :: denf
    real(double)             :: c00, c10, c01, c11, c0, c1
    real(double) :: r,z,r0,z0
    integer(long):: ir,iz
    denf=0.d0
    if(fbm%dz.eq.0)return !! there is no fast-ion distribution stored
    r=sqrt(pos(1)**2+pos(2)**2)
    z=pos(3)
    ir=floor((r-fbm%rmin)/fbm%dr) +1
    iz=floor((z-fbm%zmin)/fbm%dz) +1
    if((ir.lt.1.).or.(iz.lt.1).or.(ir.gt.fbm%nr-1).or.(iz.gt.fbm%nz-1))return
    r0=(ir-1)*fbm%dr+fbm%rmin
    z0=(iz-1)*fbm%dz+fbm%zmin
    c00 = fbm%denf(ir  ,iz  )
    c10 = fbm%denf(ir+1,iz  )
    c01 = fbm%denf(ir  ,iz+1)
    c11 = fbm%denf(ir+1,iz+1)
    !linear interpolation between C00 and C10 to find C0, C01 and C11 to find C1
    c0  = c00 + (c10 - c00) * (r - r0) /fbm%dr
    c1  = c01 + (c11 - c01) * (r - r0) /fbm%dr
    denf=  c0  + (c1  - c0) * (z - z0) /fbm%dz
  end subroutine denf_interp

  subroutine fbm_interp(pos,fbm_out) 
    !! bilinear interpolation of the effective rate coefficient tables !!
    real(double), dimension(3),   intent(in) :: pos
    real(double), dimension(fbm%nenergy,fbm%npitch),intent(out) :: fbm_out
    real(double), dimension(fbm%nenergy,fbm%npitch):: c00, c10, c01, c11, c0, c1
    real(double) :: r,z,r0,z0
    integer(long):: ir,iz
    fbm_out=0.d0
    r=sqrt(pos(1)**2+pos(2)**2)
    z=pos(3)
    ir=floor((r-fbm%rmin)/fbm%dr) +1
    iz=floor((z-fbm%zmin)/fbm%dz) +1
    if((ir.lt.1.).or.(iz.lt.1).or.(ir.gt.fbm%nr-1).or.(iz.gt.fbm%nz-1))return
    r0=(ir-1)*fbm%dr+fbm%rmin
    z0=(iz-1)*fbm%dz+fbm%zmin
    c00 = fbm%fbm(:,:,ir  ,iz  )
    c10 = fbm%fbm(:,:,ir+1,iz  )
    c01 = fbm%fbm(:,:,ir  ,iz+1)
    c11 = fbm%fbm(:,:,ir+1,iz+1)
    !linear interpolation between C00 and C10 to find C0, C01 and C11 to find C1
    c0  = c00 + (c10 - c00) * (r - r0) /fbm%dr
    c1  = c01 + (c11 - c01) * (r - r0) /fbm%dr
    fbm_out=  c0  + (c1  - c0) * (z - z0) / fbm%dz
    !! normalize
    fbm_out=fbm_out/maxval(fbm_out)
  end subroutine fbm_interp

  
  subroutine get_ac(pos,ac)
    !! determine the actual cell at which a particle is located
    real(double),dimension(3), intent(in) :: pos
    integer(long),dimension(3),intent(out):: ac
    ac(1)=floor((pos(1)-grid%xmin)/grid%dx) +1
    if(ac(1).gt.grid%nx)ac(1)=grid%nx
    if(ac(1).lt.1)ac(1)=1
    ac(2)=floor((pos(2)-grid%ymin)/grid%dy) +1
    if(ac(2).gt.grid%ny)ac(2)=grid%ny
    if(ac(2).lt.1)ac(2)=1
    ac(3)=floor((pos(3)-grid%zmin)/grid%dz) +1
    if(ac(3).gt.grid%nz)ac(3)=grid%nz
    if(ac(3).lt.1)ac(3)=1
  end subroutine get_ac


  subroutine neut_rate(denn,vi,vn,rates)
    !GET neutralization rate from tables
    real(double), dimension(nlevs),intent(in) :: denn !!density of neutrals cm-3
    real(double), dimension(3),    intent(in) :: vi,vn!!of neutrals/ions (cm/s)
    real(double), dimension(nlevs),intent(out):: rates!! rates
    real(double), dimension(nlevs,nlevs)      :: neut !!rate coeff
    real(double)                :: eb, eb0           !! relative Energy
    real(double)                :: vrel          !! relative velocity
    integer                     :: ebi      !! indizes
    real(double), dimension(nlevs,nlevs) :: c0, c1
    !Eeff 
    vrel=sqrt(dot_product(vi-vn,vi-vn))
    eb=v_to_E*vrel**2  ! [kev/amu]
    ebi= floor(eb/tables%d_eb_neut)+1   
    if(ebi.ge.tables%nr_eb_neut) stop 'EB out of range of neuttable!'
    eb0=(ebi-1)*tables%d_eb_neut
    c0=tables%neut(:,:,ebi)
    c1=tables%neut(:,:,ebi+1)
    !linear interpolation between C0 and C1
    neut  = (c0 + (c1 - c0) * (eb - eb0) / tables%d_eb_neut)
    rates=matmul(neut(:,:),denn(:))*vrel
  end subroutine neut_rate

         
  subroutine neut_rate_beam_therm(denn,vi,pos,rates)
    !! GET neutralization rate for beam-therm reactions
    real(double),  dimension(nlevs),intent(in):: denn !!density of neutrals cm-3
    real(double), dimension(3),    intent(in) :: vi!!of neutrals/ions (cm/s)
    real(double),dimension(3),    intent(in)  :: pos
    real(double), dimension(nlevs),intent(out):: rates!! rates
    real(double)                :: vtor           !! toroidal rotation
    real(double), dimension(3)  :: vrot           !! Rotation velocity of plasma
    real(double)                :: vnet_square    !! netto velocity of neutrals 
    real(double)                :: ti             !! Ion temperature
    real(double)                :: denp           !! P density
    real(double)                :: eb             !! Energy of the fast neutral
    integer                     :: ebi, tii       !! bin postions in arrays
    real(double), dimension(3)  :: vnhalo
    integer                     :: in
    real(double), dimension(nlevs) :: rates2
    real(double),dimension(nlevs+1,nlevs) :: eff_cross_section  ! output 
    rates=0.d0
    call prof_interp(pos,denp,ti,vtor)
    call calc_vrot(pos,vtor,vrot)
    ! vnet_square=dot_product(vi-vrot,vi-vrot)          ![cm/s]
    vnet_square=dot_product(vi,vi)          ![cm/s]
    eb=v_to_E*inputs%ab*vnet_square ![kev]
    !! DEUTERIUM (Impact excitation/ionization,charge exchange)
    ebi= floor(eb/tables%d_eb_cx)+1  
    tii= floor(ti/tables%d_ti_cx)+1
    if((tii.ge.tables%nr_ti_cx).or.(ebi.ge.tables%nr_eb_cx))then
       !print*, 'Eb or Ti out of range of cxtable!'
       do in=1,int(nr_halo_neutrate)
          call mc_halo(pos, vnhalo(:))
          call neut_rate(denn,vi,vnhalo,rates2)
          rates=rates+rates2/nr_halo_neutrate
       enddo
       return
    endif
    call table_interp(tables%cx(:,:,:,:),tables%d_eb_cx,tables%d_ti_cx &
         ,ebi,tii,eb,ti,eff_cross_section)
    rates=matmul(eff_cross_section(1:nlevs,:),denn(:))
  end subroutine neut_rate_beam_therm


        

  !****************************************************************************
  !----------------mc_fastion------------------------------------------------
  !****************************************************************************
  subroutine mc_fastion(pos,nlaunch,vi_arr,nlaunch2)
    !!IN: ac,   OUT: vi
    !!mcbeam computes monte carlo velocity of fast ions from cell%fbm
    real(double), dimension(3), intent(in) :: pos
    real(double), intent(in)               :: nlaunch  
    real(double),dimension(:,:),intent(out):: vi_arr 
    real(double), intent(out)              :: nlaunch2  

    real(double), dimension(3)             :: a,b,c ! vectors relative to b
    real(double), dimension(fbm%nenergy,fbm%npitch) :: fbm2
    real(double), dimension(3)             :: vi
    real(double)                           :: eb ,ptch,b_abs
    integer(long)                          :: ienergy, ipitch ,ii, j
    real(double)                           :: vabs, phi, sinus
    real(double), dimension(3)             :: randomu3
    real(double), dimension(1)             :: randomu1
    call bfield_interp(pos,b_abs,b,a,c)
    call fbm_interp(pos,fbm2)
    !! -- use a rejection method to determine vi from distrbution function --!!
    nlaunch2=0
    do j=0,floor(nlaunch)
       vi=0.d0
       rejection_loop: do ii=1,10000
          call randu(randomu3)
          eb   = fbm%emin + fbm%eran * randomu3(1)
          ptch = fbm%pmin + fbm%pran * randomu3(2) 
          !! take point in FBM distribution closest to eb, ptch.
          ienergy=floor((eb-fbm%emin)  /fbm%dE)   +1
          ipitch =floor((ptch-fbm%pmin)/fbm%dP)+1
          if(fbm2(ienergy,ipitch).gt.randomu3(3))then
             call randu(randomu1)
             vabs          = sqrt(eb/(v_to_E*inputs%ab))
             phi           = 2.d0*pi*randomu1(1)
             sinus         = sqrt(1.d0-ptch**2)
             nlaunch2=nlaunch2+1
             vi_arr(int(nlaunch2),:)=vabs*(sinus*cos(phi)*a &
                  + ptch*b*inputs%btipsign &
                  + sinus*sin(phi)*c) 
             exit rejection_loop
          endif
       enddo rejection_loop
       !! print*, 'rejection method found no solution!'
    enddo

  end subroutine mc_fastion


  subroutine mc_fastion2(ac,pos,nlaunch,vi_arr,nlaunch2)
    !!IN: ac,   OUT: vi
    !!mcbeam computes monte carlo velocity of fast ions from cell%fbm
    real(double), dimension(3), intent(in) :: pos
    integer(long), dimension(3),intent(in) :: ac
    real(double), intent(in)               :: nlaunch  
    real(double),dimension(:,:),intent(out):: vi_arr 
    real(double), intent(out)              :: nlaunch2  

    real(double)                           :: b_abs
    real(double), dimension(3)             :: a,b,c ! vectors relative to b
    real(double), dimension(fbm%nenergy,fbm%npitch) :: fbm2,nr_arr
    integer(long)                          :: ienergy, ipitch ,nrr, ii

    real(double)                           :: drr, eb ,ptch
    real(double)                           :: vabs, phi, sinus
    real(double)                           :: phimin, phimax, phi0, a1, c1

    real(double), dimension(3)             :: randomu3
    real(double), dimension(1)             :: randomu1
    call bfield_interp(pos,b_abs,b,a,c)
    call fbm_interp(pos,fbm2)


    !gyro angle boundarys phimin phimax
    phimin = 0.
    phimax = 2.d0*pi
    !! -- use a rejection method to determine vi from distrbution function --!!
    vi_arr=0.d0
    nlaunch2=0.
    nr_arr=fbm2/sum(fbm2)*nlaunch
    do ienergy=1,fbm%nenergy
       do ipitch=1,fbm%npitch
          if(nr_arr(ienergy,ipitch).le.0)cycle
          nrr=floor(nr_arr(ienergy,ipitch))
          drr=nr_arr(ienergy,ipitch)-nrr
          call randu(randomu1)
          if(drr.gt.randomu1(1)) nrr=nrr+1
          do ii=1,nrr
             call randu(randomu3)
             nlaunch2=nlaunch2+1.
             eb=fbm%energy(ienergy)+fbm%dE*(randomu3(1)-0.5)
             ptch=fbm%pitch(ipitch)+fbm%dP*(randomu3(2)-0.5)
             vabs  = sqrt(eb/(v_to_E*inputs%ab))
             phi           = (phimax-phimin)*randomu3(3)+phimin
             sinus         = sqrt(1.d0-ptch**2)
             vi_arr(int(nlaunch2),:) = vabs * (sinus*cos(phi)*a + ptch*b*inputs%btipsign &
               + sinus*sin(phi)*c) 
          enddo
       enddo
    enddo
  end subroutine mc_fastion2

 
  !****************************************************************************
  !--------------------mc_halo------------------------------------------------
  !****************************************************************************
  subroutine mc_halo(pos,vhalo)
    real(double),dimension(3),intent(in)  :: pos
    real(double),dimension(3),intent(out) :: vhalo !! velocity [cm/s]
    real(double)              :: r,ti,denp,vtor
    real(double),dimension(3) :: vrot,randomn
    call randn(randomn)  
    r=sqrt(pos(1)**2+pos(2)**2)
    call prof_interp(pos,denp,ti,vtor)
    call calc_vrot(pos,vtor,vrot)
    vhalo(:)=vrot(:)+sqrt(ti*0.5/(v_to_E*inputs%ai))*randomn(:) !![cm/s]
  end subroutine mc_halo







  !****************************************************************************
  !----------------mc_nbi------------------------------------------------------
  !****************************************************************************
  subroutine rotate(uvw_vec,updown,xyz_vec)
    real(double), dimension(3), intent(in) :: uvw_vec !! vector in uvw coords
    integer             , intent(in)       :: updown  !! source has two plates
    real(double), dimension(3), intent(out):: xyz_vec !! vector in xyz coords
    real(double), dimension(3)             :: uvz_vec   
    !! rotate uvw vector in vertical dirction 
    if(updown.lt.0) uvz_vec(:)=matmul(nbi%Arot(:,:),uvw_vec(:))
    if(updown.ge.0) uvz_vec(:)=matmul(nbi%Brot(:,:),uvw_vec(:))
    !! rotate uvz_vec by phi_box onto xyz_coordinates
    xyz_vec=matmul(nbi%Crot(:,:),uvz_vec(:))
  end subroutine rotate
  subroutine mc_nbi(vnbi,efrac,rnbi)
    !!-- mc_nbi computes monte carlo velocity and initial start position of
    !!-- NBI neutrals on the FIDASIM grid
    integer             , intent(in)          :: efrac !! energy fraction
    real(double), dimension(3), intent(out)   :: vnbi  !! velocity [cm/s]
    real(double), dimension(3), intent(out), optional :: rnbi  !! postition
    integer                      :: jj, updown
    real(double)                 :: a_dx, a_dy, a_dz
    real(double), dimension(3)   :: uvw_pos    !! Start position on ion source
    real(double), dimension(3)   :: xyz_pos    !! Start position on ion source
    real(double), dimension(3)   :: uvw_ray    !! NBI veloicity in uvw coords
    real(double), dimension(2)   :: randomu    !! uniform random numbers
    real(double), dimension(2)   :: randomn    !! normal random numbers
    integer                      :: nstep
    !! ------------ Random start postion on ion source grid -------------- !!
    call randu(randomu)
    uvw_pos(1) =  0.d0
    uvw_pos(2) =  nbi%dv * 2.d0*(randomu(1)-0.5d0)
    uvw_pos(3) =  nbi%dw * 2.d0*(randomu(2)-0.5d0)
    if(uvw_pos(3).gt.0)then 
       updown=1
    else
       updown=-1
    endif
    call randn(randomn)
    uvw_ray(1)=-1.d0
    uvw_ray(2)=uvw_ray(1)*(uvw_pos(2)/nbi%focy &
         +tan(nbi%divay(efrac)*randomn(1)))
    uvw_ray(3)=uvw_ray(1)*(uvw_pos(3)/nbi%focz &
         +tan(nbi%divaz(efrac)*randomn(2)))
    call rotate(uvw_ray,updown,vnbi(:))
    vnbi(:)=vnbi(:)/sqrt(dot_product(vnbi(:),vnbi(:)))
    call rotate(uvw_pos,updown,xyz_pos)
    xyz_pos(:)=xyz_pos(:)+nbi%xyz_pos(:)
    !! ----------- Determine start postition on FIDASIM grid --------- !!
    if(present(rnbi)) then
       nstep=anint(2000./grid%dl)
       nbi_track: do jj=1,nstep
          xyz_pos(:) = xyz_pos(:) + grid%dl * vnbi(:)
          !! check if position is inside the grid
          if ( xyz_pos(1).gt.grid%xx(1) .and. &
               xyz_pos(1).lt.grid%xx(grid%nx)+grid%dx.and. & 
               xyz_pos(2).gt.grid%yy(1) .and. &
               xyz_pos(2).lt.grid%yy(grid%ny)+grid%dy.and. &
               xyz_pos(3).gt.grid%zz(1) .and. &
               xyz_pos(3).lt.grid%zz(grid%nz)+grid%dz) then
             exit nbi_track
          endif
       enddo nbi_track
       if (jj.ge.nstep) then
          !!print*, 'NBI marker outside of grid!'
          nbi_outside=nbi_outside+1
          rnbi(:)=(/-1,0,0/)
          return
       endif
       !! add random step tothe start_position
       call randu(randomu)
       rnbi(:)=xyz_pos(:)+vnbi(:)*grid%dl*randomu(1)
    endif
    !! ---- Determine velocity of neutrals corrected by efrac ---- !!
    vnbi(:) = vnbi(:)*nbi%vinj/sqrt(real(efrac))
    if(sqrt(dot_product(vnbi,vnbi)).le.0)then
       print*, vnbi
       print*, nbi%vinj
       print*, sqrt(real(efrac))
    endif
  end subroutine mc_nbi 


  !****************************************************************************
  !----------------mc_start-------------------------------------------------
  !****************************************************************************
  subroutine mc_start(ac,vi,ipos)
    !! determine random start position within a cell and correct for gyro
    !! orbit
    integer  , dimension(3)  , intent(in)          :: ac !ind of actual cell
    real(double)   , dimension(3)  , intent(in)    :: vi !velocity [cm/s]
    real(double)   , dimension(3)  , intent(out)   :: ipos !starting position
    real(double)   , dimension(3)    :: b_norm !Magnetic field vec
    real(double)   , dimension(3)    :: vxB      ! crossproduct
    real(double)   , dimension(3)    :: r_gyro! gyro-radius
    real(double)                     :: one_over_omega,b_abs! For gyro-radius
    real(double)   , dimension(3)    :: randomu  
    call randu(randomu)  
    ipos(1)=grid%xx(ac(1))+ grid%dx*randomu(1)
    ipos(2)=grid%yy(ac(2))+ grid%dy*randomu(2) 
    ipos(3)=grid%zz(ac(3))+ grid%dz*randomu(3)
    call bfield_interp(ipos,b_abs,b_norm)
    if(b_abs.gt.0)then
       one_over_omega=inputs%ab*mass_u/(b_abs*e0)
       vxB(1)= (vi(2) * b_norm(3) - vi(3) * b_norm(2))
       vxB(2)= (vi(3) * b_norm(1) - vi(1) * b_norm(3))
       vxB(3)= (vi(1) * b_norm(2) - vi(2) * b_norm(1))
       r_gyro(:)=vxB(:)*one_over_omega
       ipos(:)=ipos(:)-r_gyro(:) !! '-'because v x B is towards the gyrocenter 
    endif
  end subroutine mc_start

  !****************************************************************************
  !----------------------------- colrad  ------------------------------
  !****************************************************************************
  ! first subroutines for eigenvalue decomposition 
   subroutine RSWAP(a,b)
    real(double) :: a,b, t
    t=a; a=b; b=t
  end subroutine RSWAP
  subroutine balance(n,     &  !size of matrix         
       mat,   &  !input matrix
       scal,  &  !Scaling data
       low,   &  !first relevant row index
       high )   !last relevant row index                 
    integer, intent(in)  :: n
    real(double)         :: mat(0:n,0:n),scal(0:n)
    integer, intent(out) :: high, low
    integer, parameter   :: basis = 2
    real(double)         :: b2, r, c, f, g, s
    integer              :: m, k, i, j, iter
    !*====================================================================*
    !*  balance balances the matrix so that the rows with zero entries    *
    !*  off the diagonal are isolated and the remaining columns and rows  *
    !*  are resized to have one norm close to 1.                          *
    !*   Input parameters:                                                *
    !*      n        integer;  ( n > 0 )                                  *
    !*               Dimension of mat                                     *
    !*      mat      n x n input matrix                                   *
    !*                                                                    *
    !*   Output parameters:                                               *
    !*      mat      n x n scaled matrix                                  *
    !*      low      integer;                                             *
    !*      high     integer;                                             *
    !*               the rows 0 to low-1 and those from high to n-1       *
    !*               contain isolated eigenvalues (only nonzero entry on  *
    !*               the diagonal)                                        *
    !*      scal     vector of size n                                     *
    !*               the vector scal contains the isolated eigenvalues in *
    !*               the positions 0 to low-1 and high to n-1, its other  *
    !*               components contain the scaling factors for           *
    !*               transforming mat.                                    *
    !*====================================================================*
    scal=0.d0
    b2 = basis * basis
    m = 0
    k = n - 1
    iter=1
    do while(iter==1)
       iter = 0
       do j = k, 0, -1
          r = ZERO
          do i = 0, k
             if (i.ne.j)  r = r + DABS(mat(j,i))
          enddo
          if (r == ZERO) then
             scal(k) = j
             if (j.ne.k) then
                do i = 0, k 
                   call RSWAP(mat(i,j), mat(i,k))
                enddo
                do i = m, n-1 
                   call RSWAP(mat(j,i), mat(k,i))
                enddo
             endif
             k=k-1
             iter = 1
          endif
       enddo !j loop
    enddo !while iter=1
    iter=1
    do while (iter==1)
       iter = 0
       do j = m, k
          c = ZERO
          do i = m, k
             if (i.ne.j)  c = c + DABS(mat(i,j))
          enddo
          if (c == ZERO) then
             scal(m) = j
             if (j.ne.m) then
                do i = 0, k 
                   call RSWAP(mat(i,j), mat(i,m))
                enddo
                do i = m, n-1 
                   call RSWAP(mat(j,i), mat(m,i))
                enddo
             endif
             m = m + 1
             iter = 1
          endif
       enddo !j loop
    enddo !while iter=1
    low = m
    high = k
    do i = m, k 
       scal(i) = ONE
    enddo
    iter=1
    do while (iter==1)
       iter = 0
       do i = m, k
          c=ZERO; r=ZERO
          do j = m, k
             if (j.ne.i) then
                c = c + DABS(mat(j,i))
                r = r + DABS(mat(i,j))
             endif
          enddo
          g = r / basis
          f = ONE
          s = c + r
          do while (c < g)
             f = f * basis
             c = c * b2
          enddo
          g = r * basis
          do while (c >= g)
             f = f / basis
             c = c / b2
          enddo
          if ((c + r) / f < 0.95 * s) then
             g = ONE / f
             scal(i) = scal(i) * f
             iter = 1
             do j = m, n-1 
                mat(i,j) = mat(i,j) * g
             enddo
             do j = 0, k  
                mat(j,i) = mat(j,i) * f
             enddo
          endif
       enddo !i loop
    enddo !while iter=1
    return
  end subroutine balance
  subroutine balback(n,     &  !Dimension of matrix .........
       low,   &  !first nonzero row ...........
       high,  &  !last nonzero row ............
       scal,  &  !Scaling data ................
       eivec )   !Eigenvectors ................
    integer,intent(in)         :: high, low
    integer,intent(in)         :: n
    real(double), intent(in)   ::  scal(0:n)
    real(double), intent(inout):: eivec(0:n,0:n)
    real(double) :: s
    integer      :: i,j,k
    !*====================================================================*
    !*  balback reverses the balancing of balance for the eigenvactors.   *
    !*   Input parameters:                                                *
    !*   ----------------                                                 *
    !*      n        integer;  ( n > 0 )                                  *
    !*               Dimension of mat                                     *
    !*      low      integer;                                             *
    !*      high     integer;   see balance                               *
    !*      eivec    n x n matrix of eigenvectors, as computed in  qr2    *
    !*      scal     vector of size n;                                    *
    !*               Scaling data from  balance                           *
    !*   Output parameter:                                                *
    !*   ----------------                                                 *
    !*      eivec    n x n matrix;                                        *
    !*               Non-normalized eigenvectors of the original matrix   *
    !*====================================================================*
    do i = low, high
       s = scal(i)
       do j = 0, n-1  
	  eivec(i,j) = eivec(i,j) * s
       enddo
    enddo
    do i = low-1, 0, -1
       k = Int(scal(i))
       if (k.ne.i) then
          do j = 0, n-1
             call RSWAP(eivec(i,j), eivec(k,j))
          enddo
       endif
    enddo
    do i = high + 1, n-1
       k = Int(scal(i))
       if (k.ne.i) then
          do j = 0, n-1
             call RSWAP(eivec(i,j), eivec(k,j))
          enddo
       endif
    enddo
    return
  end subroutine balback
  subroutine elmhes(n,    &  !Dimension of matrix
       low,  &  !first nonzero row ...........
       high, &  !last nonzero row ............
       mat,  &  !input/output matrix .........
       perm )  !Permutation vector ..........
    integer,intent(in)   :: n
    integer,intent(in)   :: high, low
    real(double), intent(inout):: mat(0:n,0:n)
    integer,intent(out)  :: perm(0:n)
    integer              :: i, j, m
    real(double)         ::  x, y
    !*====================================================================*
    !*  elmhes transforms the matrix mat to upper Hessenberg form.        *
    !*   Input parameters:                                                *
    !*      n        integer;  ( n > 0 )                                  *
    !*               Dimension of mat                                     *
    !*      low      integer;                                             *
    !*      high     integer; see  balance                                *
    !*      mat      n x n matrix                                         *
    !*   Output parameter:                                                *
    !*      mat      n x n matrix;                                        *
    !*               upper Hessenberg matrix; additional information on   *
    !*               the transformation is stored in the lower triangle   *
    !*      perm     integer vector of size n;                            *
    !*               Permutation vector for elmtrans                      *
    !*====================================================================*
    do m = low + 1, high-1
       i = m
       x = ZERO
       do j = m, high
          if (DABS(mat(j,m-1)) > DABS (x)) then
             x = mat(j,m-1)
             i = j
          endif
       enddo
       perm(m) = i
       if (i.ne.m) then
          do j = m - 1, n-1 
             call RSWAP(mat(i,j), mat(m,j))
          enddo
          do j = 0, high 
             call RSWAP(mat(j,i), mat(j,m))
          enddo
       endif
       if (x.ne.ZERO) then
          do i = m + 1, high
             y = mat(i,m-1)
             if (y.ne.ZERO) then
                y = y / x
                mat(i,m-1) = y
                do j = m, n-1 
                   mat(i,j) = mat(i,j) - y * mat(m,j)
                enddo
                do j = 0, high 
                   mat(j,m) = mat(j,m) + y * mat(j,i)
                enddo
             endif
          enddo !i loop
       endif !x <> ZERO
    enddo !m loop
  end subroutine elmhes
  Subroutine elmtrans(n,    &  !Dimension of matrix .........
       low,  &  !first nonzero row ...........
       high, &  !last nonzero row ............
       mat,  &  !input matrix ................
       perm, &  !row permutations ............
       h  )     !Hessenberg matrix ...........
    integer,intent(in)         :: n
    integer,intent(in)         :: high, low
    real(double), intent(in)   :: mat(0:n,0:n)
    integer,intent(in)         :: perm(0:n)
    real(double),intent(out)   :: h(0:n,0:n)
    integer                    :: i, j, k
    !*====================================================================*
    !*  Elmtrans copies the Hessenberg matrix stored in mat to h.         *
    !*   Input parameters:                                                *
    !*      n        integer;  ( n > 0 )                                  *
    !*               Dimension of  mat and eivec                          *
    !*      low      integer;                                             *
    !*      high     integer; see  balance                                *
    !*      mat      n x n input matrix                                   *
    !*      perm     Integer vector of size n;                            *
    !*               Permutation data from  elmhes                        *
    !*   Output parameter:                                                *
    !*      h        n x n matrix;                                        *
    !*               Hessenberg matrix                                    *
    !*====================================================================*
    do i = 0, n-1
       do k = 0, n-1 
	  h(i,k) = ZERO
       enddo
       h(i,i) = ONE
    enddo
    do i = high - 1, low+1, -1
       j = perm(i)
       do k = i + 1, high 
	  h(k,i) = mat(k,i-1)
       enddo
       if (i.ne.j) then
          do k = i, high
             h(i,k) = h(j,k)
             h(j,k) = ZERO
          enddo
          h(j,i) = ONE
       endif
    enddo
  end subroutine elmtrans
  subroutine Comdiv(ar,     &       !Real part of numerator ..........
       ai,     &       !Imaginary part of numerator .....
       br,     &       !Real part of denominator ........
       bi,     &       !Imaginary part of denominator ...
       cr,     &       !Real part of quotient ...........
       ci,     &       !Imaginary part of quotient ......
       rc )            !return code .....................
    real(double) ::  ar,ai,br,bi,cr,ci
    integer      :: rc
    real(double) :: tmp
    !*====================================================================*
    !*  Complex division  c = a / b                                       *
    !*   Input parameters:                                                *
    !*   ================                                                 *
    !*      ar,ai    real, imaginary parts of numerator                   *
    !*      br,bi    real, imaginary parts of denominator                 *
    !*   Output parameters:                                               *
    !*   ==================                                               *
    !*      cr,ci     real , imaginary parts of the quotient              *
    !*====================================================================*
    if (br == ZERO.AND.bi == ZERO) then
       rc = 1
       return
    endif
    if (dabs(br) > dabs(bi)) then
       tmp = bi / br
       br  = tmp * bi + br
       cr  = (ar + tmp * ai) / br
       ci  = (ai - tmp * ar) / br
    else
       tmp = br / bi
       bi  = tmp * br + bi
       cr  = (tmp * ar + ai) / bi
       ci  = (tmp * ai - ar) / bi
    endif
    rc = 0
  end subroutine Comdiv !Comdiv
  function comabs(ar,ai)          !Real part ,Imaginary part ................. 
    real(double) :: ar,ai
    real(double) :: comabs
    !*====================================================================*
    !*   Input parameters:                                                *
    !*      ar,ai     Real, imaginary parts of  a                         *
    !*   Return value :                                                   *
    !*      Absolute value of a (real)                                    *
    !*====================================================================*
    if (ar == ZERO.and.ai == ZERO) then
       Comabs = ZERO
       return
    endif
    ar = DABS(ar)
    ai = DABS(ai)
    if (ai > ar) then                                  !Switch  ai and ar
       call RSWAP(ai, ar)
    endif
    if (ai == ZERO) then
       Comabs = ar
    else
       Comabs = ar * DSQRT(ONE + ai / ar * ai / ar)
    endif
  end function comabs
  subroutine  hqrvec(n,     & !Dimension of matrix .......
       low,   & !first nonzero row .........
       high,  & !last nonzero row ..........
       h,     & !upper Hessenberg matrix ...
       wr,    & !Real parts of evalues .....
       wi,    & !Imaginary parts of evalues 
       eivec, & !Eigenvectors ..............
       rc  )   !return code ...............
    integer,intent(in)         :: n
    integer,intent(in)         :: high, low
    real(double), intent(in)   :: wr(0:n),wi(0:n)
    real(double), intent(out)  :: eivec(0:n,0:n)
    real(double)  :: h(0:n,0:n)
    integer :: rc
    integer :: i, j, m, k, na, l
    integer :: code, en
    real(double)  :: p, q, r, s, t, w, x, y, z, ra, sa, vr, vi, norm, temp
    !*====================================================================*
    !*  hqrvec computes the eigenvectors for the eigenvalues found in hqr2*
    !*   Input parameters:                                                *
    !*   ================                                                 *
    !*      n        int n;  ( n > 0 )                                    *
    !*               Dimension of  mat and eivec, number of eigenvalues.  *
    !*      low      int low;                                             *
    !*      high     int high; see  balance                               *
    !*      h        n x n upper Hessenberg matrix                        *
    !*      wr       vector of size n;                                    *
    !*               Real parts of the n eigenvalues.                     *
    !*      wi       vector of size n;                                    *
    !*               Imaginary parts of the n eigenvalues.                *
    !*   Output parameter:                                                *
    !*   ================                                                 *
    !*      eivec    n x n matrix, whose columns are the eigenvectors     *
    !*====================================================================*
    r=ZERO; s=ZERO; z=ZERO; norm=ZERO
    do i = 0, n-1                               !find norm of h
       do j = i, n-1
          norm = norm + DABS(h(i,j))
       enddo
    enddo
    if (norm == ZERO) then
       rc = 1                                    !zero matrix
       return
    endif
    do en = n-1, 0, -1                          !transform back
       p = wr(en)
       q = wi(en)
       na = en - 1
       if (q == ZERO) then
          m = en
          h(en,en) = ONE
          do i = na, 0, -1
             w = h(i,i) - p
             r = h(i,en)
             do j = m, na 
                r = r + h(i,j) * h(j,en)
             enddo
             if (wi(i) < ZERO) then
                z = w
                s = r
             else
                m = i
                if (wi(i) == ZERO) then
                   if (w.ne.ZERO) then 
                      temp = w 
                   else 
                      temp=XMACH_EPS * norm
                   endif
                   h(i,en) = -r/temp            
                else
                   !Solve the linear system:
                   !| w   x |  | h[i][en]   |   | -r |
                   !|       |  |            | = |    |
                   !| y   z |  | h[i+1][en] |   | -s |
                   x = h(i,i+1)
                   y = h(i+1,i)
                   q = (wr(i) - p)**2 + wi(i)**2
                   h(i,en) = (x * s - z * r) / q
                   t = h(i,en)
                   if (DABS(x) > DABS(z)) then 
                      temp = (-r -w * t) / x 
                   else 
                      temp = (-s -y * t) / z
                   endif
                   h(i+1,en) = temp
                endif
             endif !wi[i] < 0
          enddo !i loop
       else if (q < ZERO) then
          m = na
          if (DABS(h(en,na)) > DABS(h(na,en))) then
             h(na,na) = - (h(en,en) - p) / h(en,na)
             h(na,en) = - q / h(en,na)
          else
             call Comdiv(-h(na,en),0.d0, h(na,na)-p, q, h(na,na), h(na,en),code)
          endif
          h(en,na) = ONE
          h(en,en) = ZERO
          do i = na - 1, 0, -1
             w = h(i,i) - p
             ra = h(i,en)
             sa = ZERO
             do j = m, na
                ra = ra + h(i,j) * h(j,na)
                sa = sa + h(i,j) * h(j,en)
             enddo
             if (wi(i) < ZERO) then
                z = w
                r = ra
                s = sa
             else
                m = i
                if (wi(i) == ZERO) then
                   call Comdiv(-ra, -sa, w, q, h(i,na), h(i,en),code)
                else
            !  solve complex linear system:
                   !| w+i*q     x | | h[i][na] + i*h[i][en]  |   | -ra+i*sa |
                   !|             | |                        | = |          |
            !|   y    z+i*q| | h[i+1][na]+i*h[i+1][en]|   | -r+i*s   |
                   x = h(i,i+1)
                   y = h(i+1,i)
                   vr = (wr(i) - p)**2 + wi(i)**2 - q*q
                   vi = TWO * q * (wr(i) - p)
                   if (vr == ZERO.AND.vi == ZERO) then
                      vr = XMACH_EPS * norm * (DABS(w) + DABS(q)  &
                           + DABS(x) + DABS(y) + DABS(z))
                   endif
                   
                   call Comdiv (x*r-z*ra+q*sa,x*s-z*sa-q*ra &
                        ,vr,vi,h(i,na),h(i,en),code)
                   if (DABS(x) > DABS(z) + DABS(q)) then
                      h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
                      h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
                   else
                      call Comdiv (-r - y * h(i,na), -s - y * h(i,en) &
                           , z, q, h(i+1,na), h(i+1,en),code)
                   endif
                endif !wi[i] = 0
             endif !wi[i] < 0
          enddo !i loop
       endif !else if q < 0
    enddo !en loop
    do i = 0, n-1                        !Eigenvectors for the evalues for
       if (i < low.or.i > high) then      !rows < low  and rows  > high
          do k = i + 1, n-1
             eivec(i,k) = h(i,k)
          enddo
       endif
    enddo
    j = n-1
    do while (j>=low)
       if(j<=high)then
	  m =j 
       else 
          j = high
       endif
       if (wi(j) < ZERO) then
          l=j-1
          do i = low, high
             y=ZERO; z=ZERO
             do k = low, m
                y = y + eivec(i,k) * h(k,l)
                z = z + eivec(i,k) * h(k,j)
             enddo
             eivec(i,l) = y
             eivec(i,j) = z
          enddo
       else
          if (wi(j) == ZERO) then
             do i = low, high
                z = ZERO
                do k = low, m
                   z = z + eivec(i,k) * h(k,j)
                enddo
                eivec(i,j) = z
             enddo
          endif
       endif
       j = j - 1
    enddo !j loop
    rc = 0
  end subroutine hqrvec
  subroutine hqr2(n,    &  !Dimension of matrix .........
       low,  &  !first nonzero row ...........
       high, &  !last nonzero row ............
       h,    &  !Hessenberg matrix ...........
       wr,   &  !Real parts of eigenvalues ...
       wi,   &  !Imaginary parts of evalues ..
       eivec,&  !Matrix of eigenvectors ......
       cnt,  &  !Iteration counter ...........
       rc   )    !return code .................              
    integer,intent(in)          :: n
    integer,intent(in)          :: high, low
    real(double) ,intent(out)   :: h(0:n,0:n)
    real(double), intent(out)   :: wr(0:n),wi(0:n)
    real(double), intent(out)   :: eivec(0:n,0:n)
    integer,intent(out)         :: rc
    integer,intent(out)         :: cnt(0:n)
    integer :: en
    integer :: i, j, na, iter, l, ll, m, k
    real(double)  :: p, q, r, s, t, w, x, y, z
    !**********************************************************************
    !* hqr2 computes the eigenvalues and (if vec = True) the eigenvectors *
    !* of an  n * n upper Hessenberg matrix.                              *
    !*   Input parameters:                                                *
    !*   ----------------                                                 *
    !*      n        integer;  ( n > 0 )                                  *
    !*               Dimension of  h and eivec,                           *
    !*               length of the real parts vector  wr and of the       *
    !*               imaginary parts vector  wi of the eigenvalues.       *
    !*      low      integer;                                             *
    !*      high     integer;  see balance                                *
    !*      h        n x n matrix;                                        *
    !*               upper Hessenberg matrix as output of Elmhes          *
    !*               (destroyed in the process).                          *
    !*   Output parameters:                                               *
    !*   -----------------                                                *
    !*      eivec    n x n matrix;  (only if vec = 1)                     *
    !*               Matrix, which for vec = 1 contains the               *
    !*               eigenvectors as follows:                             *
    !*               For real eigebvalues the corresponding column        *
    !*               contains the corresponding eigenvactor, while for    *
    !*               complex eigenvalues the corresponding column contains*
    !*               the real part of the eigenvactor with its imaginary  *
    !*               part is stored in the subsequent column of eivec.    *
    !*               The eigenvactor for the complex conjugate eigenvactor*
    !*               is given by the complex conjugate eigenvactor.       *
    !*      wr       vector of size n;                                    *
    !*               Real part of the n eigenvalues.                      *
    !*      wi       vector of size n;                                    *
    !*               Imaginary parts of the eigenvalues                   *
    !*      cnt      Integer vector of size n;                            *
    !*               vector of iterations used for each eigenvalue.       *
    !*               For a complex conjugate eigenvalue pair the second   *
    !*               entry is negative.                                   *
    !**********************************************************************
    p=ZERO; q=ZERO; r=ZERO 
    do i = 0, n-1
       if (i < low.or.i > high) then
          wr(i) = h(i,i)
          wi(i) = ZERO
          cnt(i) = 0
       endif
    enddo
    en = high
    t = ZERO
    do while (en >= low)
       iter = 0
       na = en - 1
       do while(1<2)
          ll=999                          
          do l = en, low+1, -1                      !search for small
             !subdiagonal element
             if(DABS(h(l,l-1))<=XMACH_EPS*(DABS(h(l-1,l-1))+DABS(h(l,l))))then
                ll=l;      !save current index
                goto 10    !exit l loop
             endif
          enddo
10        if(ll.ne.999)then 
             l=ll 
	  else 
             l=0          !restore l
          endif
          x = h(en,en)
          if (l == en) then                         !found one evalue
             wr(en) = x + t
             h(en,en) = x + t
             wi(en) = ZERO
             cnt(en) = iter
             en = en - 1
             goto 15      !exit from loop while(True)
          endif
          y = h(na,na)
          w = h(en,na) * h(na,en)
          if (l == na) then                         !found two evalues
             p = (y - x) * 0.5d0
             q = p * p + w
             z = DSQRT(DABS(q))
             x = x + t
             h(en,en) = x + t
             h(na,na) = y + t
             cnt(en) = -iter
             cnt(na) = iter
             if (q >= ZERO) then                     !real eigenvalues
                if (p<ZERO) then 
                   z=p-z 
                else 
                   z=p+z
                endif
                wr(na) = x + z
                wr(en) = x - w / z
                s = w - w / z
                wi(na) = ZERO
                wi(en) = ZERO
                x = h(en,na)
                r = DSQRT (x * x + z * z)
                p = x / r
                q = z / r
                do j = na, n-1
                   z = h(na,j)
                   h(na,j) = q * z + p * h(en,j)
                   h(en,j) = q * h(en,j) - p * z
                enddo
                do i = 0, en
                   z = h(i,na)
                   h(i,na) = q * z + p * h(i,en)
                   h(i,en) = q * h(i,en) - p * z
                enddo
                do i = low, high
                   z = eivec(i,na)
                   eivec(i,na) = q * z + p * eivec(i,en)
                   eivec(i,en) = q * eivec(i,en) - p * z
                enddo
             else                                  !pair of complex
                wr(na) = x + p
                wr(en) = x + p
                wi(na) =   z
                wi(en) = - z
             endif !if q>=ZERO
             en = en - 2
             goto 15                               !exit while(1<2)
          endif !if l = na
          if (iter >= MAXIT) then
             cnt(en) = MAXIT + 1
             rc = en
             write(*,*) ' stop at iter >= MAXIT.'
             return
          endif
          if (iter.ne.0.and.MOD(iter,10) == 0) then
             t = t + x
             do i = low, en 
                h(i,i) = h(i,i) - x
             enddo
             s = DABS(h(en,na)) + DABS(h(na,en-2))
             x = 0.75d0 * s; y = x
             w = -0.4375d0 * s * s
          endif
          iter = iter + 1
          do m = en - 2, l, -1
             z = h(m,m)
             r = x - z
             s = y - z
             p = ( r * s - w ) / h(m+1,m) + h(m,m+1)
             q = h(m + 1,m + 1) - z - r - s
             r = h(m + 2,m + 1)
             s = DABS(p) + DABS(q) + DABS (r)
             p = p / s
             q = q / s
             r = r / s
             if (m == l)  goto 12
             if (DABS(h(m,m-1)) * (DABS(q) + DABS(r)) <= XMACH_EPS * DABS(p) &
                  * (DABS(h(m-1,m-1)) + DABS(z) + DABS(h(m+1,m+1)))) then
                goto 12                !exit m loop
             endif
          enddo
12        do i = m + 2, en 
             h(i,i-2) = ZERO
          enddo
          do i = m + 3, en 
             h(i,i-3) = ZERO
          enddo
          do k = m, na
             if(k.ne.m)then!double QR step, for rows l to en and columns m to en
                p = h(k,k-1)
                q = h(k+1,k-1)
                if (k.ne.na) then 
                   r = h(k+2,k-1) 
                else 
                   r = ZERO
                endif
                x = DABS(p) + DABS(q) + DABS(r)
                if (x == ZERO) goto 30                  !next k
                p = p / x
                q = q / x
                r = r / x
             endif
             s = DSQRT(p * p + q * q + r * r)
             if (p < ZERO) s = -s
             if (k.ne.m) then
                h(k,k-1) = -s * x
             else if (l.ne.m) then
                h(k,k-1) = -h(k,k-1)
             endif
             p = p + s
             x = p / s
             y = q / s
             z = r / s
             q = q / p
             r = r / p
             do j = k, n-1                          !modify rows
                p = h(k,j) + q * h(k+1,j)
                if (k.ne.na) then
                   p = p + r * h(k+2,j)
                   h(k+2,j) = h(k+2,j) - p * z
                endif
                h(k+1,j) = h(k+1,j) - p * y
                h(k,j)   = h(k,j) - p * x
             enddo
             if (k+3 < en) then 
                j=k+3 
             else 
                j=en
             endif
             do i = 0, j                            !modify columns
                p = x * h(i,k) + y * h(i,k+1)
                if (k.ne.na) then
                   p = p + z * h(i,k+2)
                   h(i,k+2) = h(i,k+2) - p * r
                endif
                h(i,k+1) = h(i,k+1) - p * q
                h(i,k)   = h(i,k) - p
             enddo
             do i = low, high
                p = x * eivec(i,k) + y * eivec(i,k+1)
                if (k.ne.na) then
                   p = p + z * eivec(i,k+2)
                   eivec(i,k+2) = eivec(i,k+2) - p * r
                endif
                eivec(i,k+1) = eivec(i,k+1) - p * q
                eivec(i,k)   = eivec(i,k) - p
             enddo
30           continue
          enddo !k loop
       enddo !while(1<2)
15  continue
    enddo !while en >= low                         All evalues found
    !transform evectors back
    call hqrvec (n, low, high, h, wr, wi, eivec,rc)
  end subroutine hqr2
  
  subroutine eigen (matrix, eigvec, eigval)    
    real(double) ,intent(in),dimension(nlevs,nlevs)  :: matrix
    real(double) ,intent(out),dimension(nlevs,nlevs) :: eigvec
    real(double) ,intent(out),dimension(nlevs)       :: eigval   
    real(double)     :: mat(0:nlevs,0:nlevs)  
    real(double)     :: eivec(0:nlevs,0:nlevs)
    real(double)     :: valre(0:nlevs) !real parts of eigenvalues
    real(double)     :: valim(0:nlevs) !imaginary parts of eigenvalues
    integer          :: rc             !return code
    integer          :: cnt(0:nlevs)   !Iteration counter
    integer          :: high, low
    real(double)     :: d(0:nlevs), scale(0:nlevs)
    integer          :: perm(0:nlevs)
    integer          :: i,j,k ! counter
    real(double)     :: w,v,norm
    integer          :: n ! nlevels
    n=nlevs
    !**********************************************************************
    !* The subroutine eigen  determines all eigenvalues and (if desired)  *
    !* all eigenvectors of a real square  n * n  matrix via the QR method *
    !* in the version of Martin, Parlett, Peters, Reinsch and Wilkinson.  * 
    !*   Litterature:                                                     *
    !*   -----------                                                      *
    !*      1) Peters, Wilkinson: Eigenvectors of real and complex        *
    !*         matrices by LR and QR triangularisations,                  *
    !*         Num. Math. 16, p.184-204, (1970); [PETE70]; contribution   *
    !*         II/15, p. 372 - 395 in [WILK71].                           *
    !*      2) Martin, Wilkinson: Similarity reductions of a general      *
    !*         matrix to Hessenberg form, Num. Math. 12, p. 349-368,(1968)*
    !*         [MART 68]; contribution II,13, p. 339 - 358 in [WILK71].   *
    !*      3) Parlett, Reinsch: Balancing a matrix for calculations of   *
    !*         eigenvalues and eigenvectors, Num. Math. 13, p. 293-304,   *
    !*         (1969); [PARL69]; contribution II/11, p.315 - 326 in       *
    !*         [WILK71].                                                  *
    !*   Input parameters:                                                *
    !*   ----------------                                                 *
    !*      n        integer; ( n > 0 )                                   *
    !*               size of matrix, number of eigenvalues                *
    !*      mat      n x n matrix;                                        *
    !*               input matrix                                         *
    !*   Output parameters:                                               *
    !*   -----------------                                                *
    !*      eivec    n x n matrix;     (only if vec = 1)                  *
    !*               matrix, if  vec = 1  that holds the eigenvectors     *
    !*               thus :                                               *
    !*               If the jth eigenvalue of the matrix is real then the *
    !*               jth column is the corresponding real eigenvector;    *
    !*               if the jth eigenvalue is complex then the jth column *
    !*               of eivec contains the real part of the eigenvector   *
    !*               while its imaginary part is in column j+1.           *
    !*               (the j+1st eigenvector is the complex conjugate      *
    !*               vector.)                                             *
    !*      valre    vector of size n;                                    *
    !*               Real parts of the eigenvalues.                       *
    !*      valim    vector of size n;                                    *
    !*               Imaginary parts of the eigenvalues                   *
    !*      cnt      Integer vector of size n;                            *
    !*               vector containing the number of iterations for each  *
    !*               eigenvalue. (for a complex conjugate pair the second *
    !*               entry is negative).                                  *
    !**********************************************************************
    cnt=0 ; d=0.d0
    mat(0:n-1,0:n-1)=matrix(1:n,1:n)
    !balance mat for nearly
    call balance(n, mat, scale, low, high)      !equal row and column 
    !reduce mat to upper
    call elmhes(n, low, high, mat, perm)        !reduce mat to upper
    !Hessenberg form
    call elmtrans(n, low, high, mat, perm, eivec)
    !QR algorithm for eigenvalues and eigenvectors
    call hqr2(n, low, high, mat, valre, valim, eivec, cnt,rc)  
    !reverse balancing to determine eigenvectors
    call balback(n, low, high, scale, eivec) 
    if (rc.ne.0) stop 'problem in eigen!'
    eigval(1:n)=valre(0:n-1)
    eigvec(1:n,1:n)=eivec(0:n-1,0:n-1)
  end subroutine eigen

  function outerprod(a,b)
    real(double), dimension(:), intent(IN)   :: a,b
    real(double), dimension(size(a),size(b)) :: outerprod
    outerprod = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  end function outerprod
  subroutine swap(a,b)
    real(double), dimension(:), intent(INOUT) :: a,b
    real(double), dimension(size(a))          :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap
  subroutine ludcmp(a,indx,d)
    real(double), dimension(:,:),intent(INOUT):: a
    integer,dimension(:),  intent(OUT)        :: indx
    real(double),                intent(OUT)  :: d
    real(double), dimension(size(a,1))        :: vv
    integer,dimension(1)                      :: imaxloc
    integer :: j,n,imax
    n=size(indx)
    d=1.0
    vv=maxval(abs(a),dim=2)
    if(any(vv.eq.0.))stop 'singular matrix in ludcmp'
    vv=1.d0/vv
    do j=1,n
       imaxloc=maxloc(vv(j:n)*abs(a(j:n,j)))
       imax=(j-1)+imaxloc(1)
       if (j /= imax) then
          call swap(a(imax,:),a(j,:))
          d=-d
          vv(imax)=vv(j)
       endif
       indx(j)=imax
       if (a(j,j) == 0.0) a(j,j)=1.0d-20
       a(j+1:n,j)=a(j+1:n,j)/a(j,j)
       a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
    enddo
  end subroutine ludcmp
  subroutine lubksb(a,indx,b)
    real(double), dimension(:,:),intent(IN)   :: a
    integer,dimension(:),  intent(IN)         :: indx
    real(double), dimension(:),  intent(INOUT):: b
    integer       :: i,n,ii,ll
    real(double)  :: summ
    n=size(indx)
    ii=0
    do i=1,n
       ll=indx(i)
       summ=b(ll)
       b(ll)=b(i)
       if (ii /= 0) then
          summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
       else if (summ /= 0.0) then
          ii=i
       endif
       b(i)=summ
    enddo
    do i=n,1,-1
       b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
    enddo
  end subroutine lubksb
  subroutine matinv(a, b)
    !! - Matrix inversion with LU-decomposition
    !====================================================
    real(double), dimension(:,:), intent(IN)             :: a
    real(double), dimension(:,:), intent(OUT)            :: b
    real(double), dimension(size(a,dim=1),size(a,dim=2)) :: ah, y
    integer                                              :: i, N
    integer, dimension(size(a,dim=1))                    :: indx
    real(double)                                         :: d
    N = size(a,dim=1)
    if (N /= size(a,dim=2)) stop 'SUB matinv: ludcmp matrix must be square!'
    ah = a
    y  = 0.
    do i = 1, N
       y(i,i) = 1.d0
    enddo
    call ludcmp(ah,indx,d)
    do i = 1, N
       call lubksb(ah,indx,y(:,i))
    enddo
    b = y
  end subroutine matinv

  subroutine table_interp(coef_matrix,deb,dti,ebi,tii,eb,ti,interp_out)
    !! bilinear interpolation of the effective rate coefficient tables !!
    real(double),dimension(:,:,:,:), intent(in) :: coef_matrix ! (m,n,eb,ti)
    real(double),                intent(in)     :: deb, dti ! eb and ti grid
    integer,                    intent(in)      :: ebi, tii    ! lower indices
    real(double),                intent(in)     :: eb, ti      ! desired values
    real(double),dimension(nlevs+1,nlevs),intent(out) :: interp_out  ! output 
    real(double), dimension(nlevs+1,nlevs)      :: c00, c10, c01, c11, c0, c1
    real(double)    :: eb0, ti0
    eb0=(ebi-1)*deb
    ti0=(tii-1)*dti
    c00 = coef_matrix(:,:,ebi  ,tii  )
    c10 = coef_matrix(:,:,ebi+1,tii  )
    c01 = coef_matrix(:,:,ebi  ,tii+1)
    c11 = coef_matrix(:,:,ebi+1,tii+1)
    !linear interpolation between C00 and C10 to find C0, C01 and C11 to find C1
    c0  = c00 + (c10 - c00) * (eb - eb0) / deb
    c1  = c01 + (c11 - c01) * (eb - eb0) / deb
    interp_out=( c0  + (c1  - c0) * (ti - ti0) / dti)
  end subroutine table_interp
  

  subroutine colrad(pos,vn,dt,states,photons,neut_type,nlaunch,ipos)
    !colrad solves the collisional-radiative balance equations
    real(double),dimension(3),    intent(in)    :: pos
    real(double),dimension(3),    intent(in),optional  :: ipos   !! for npa 
    real(double) , dimension(:),  intent(in)    :: vn  !!velocitiy (cm/s)
    real(double)               ,  intent(in)    :: dt  !!time interval in cell
    integer                    ,  intent(in)    :: neut_type!!type of neutral
    real(double)               ,  intent(in)    :: nlaunch !! nr of markers
    real(double), dimension(:) ,  intent(inout) :: states  !!density of states
    real(double),                 intent(out)   :: photons !!emitted photons
    !! ---- to determine rate coefficients ---- !
    integer(long),dimension(3)                 :: ac
    real(double), dimension(nlevs+1,nlevs)     :: qp !! Proton rate coefficants
    real(double), dimension(nlevs+1,nlevs)     :: qi !! Impurity rate coefficant
    real(double), dimension(nlevs+1,nlevs)     :: qe !! Electron rate coefficant
    real(double), dimension(nlevs,nlevs)       :: matrix  !! Matrix
    real(double), dimension(nlevs,nlevs)       :: eigvec
    real(double), dimension(nlevs)             :: eigval

    real(double)                :: vtor           !! toroidal rotation
    real(double), dimension(3)  :: vrot           !! Rotation velocity of plasma
    real(double)                :: vnet_square    !! netto velocity of neutrals 
    real(double)                :: ti,te          !! Ion/electron temperature
    real(double)                :: denp,dene,deni !! P/impurity/electron density
    real(double)                :: eb             !! Energy of the fast neutral
    integer                     :: ebi, tii,tei   !! bin postions in arrays
    !! ---- Solution of differential equation  ---- ! 
    real(double),   dimension(nlevs,nlevs)  :: eigvec_inv
    real(double),   dimension(nlevs)        :: coef
    real(double),   dimension(nlevs)        :: exp_eigval_dt 
    real(double),   dimension(nlevs)        :: dens !! Density of neutrals 
    real(double)                            :: iflux !!Initial total flux
    real(double)                            :: dflux !! change of flux
    integer                                 :: n !! counter 
    real(double),   dimension(3)            :: b_norm  !! pitch of particle
    real(double)                            :: b_abs
    real(double)                            :: ptch  !! pitch of particle
    integer                                 :: ipitch !! index of pitch
    real(double)                            :: dray !! (for NPA)
    real(double), dimension(3)              :: ray     !! ray towards NPA
    integer                                 :: counter !!of velocity vectors
    integer                                 :: vel_vec_red_lev !! reduction level
    real(double), dimension(nvelocities/2)  :: velocity_dummy
    photons=0.d0
    iflux=sum(states)
    !! --------------- Check if inputs are valid for colrad -------------- !!
    if(iflux.lt.colrad_threshold .and. inputs%npa.eq.0)then
      !! print*, 'threshold!',ac
       return
    endif
    call prof_interp(pos,denp,ti,vtor,te,dene,deni)
    call calc_vrot(pos,vtor,vrot)
    call get_ac(pos,ac)
    !! IF the density is too low, stop simulation
    !! => particles in the SOL
    !! (stopped by return of photons=0.!)
    if(dene.le.0)then
       if(neut_type.le.3)then  !! Store density for NBI simulation!
          dens(:)=states*dt/nlaunch!![neutrals/(cm^3)]!!
          result%neut_dens(ac(1),ac(2),ac(3),:,neut_type)= & 
               result%neut_dens(ac(1),ac(2),ac(3),:,neut_type)+dens(:)
       endif
       if(inputs%npa.eq.1.and.present(ipos).and.neut_type.ge.5)then
          if(sqrt(pos(1)**2+pos(2)**2).lt.160.)then !! if on highfield side
             photons=1.
          else
             dray=sqrt(dot_product(ipos-diag%xyzhead(1,:) &
                  ,ipos-diag%xyzhead(1,:)))
             ray=ipos(:)+vn(:)/sqrt(dot_product(vn,vn))*dray
             !$OMP CRITICAL(col_rad_npa)
             npa%counter=npa%counter+1
             if(npa%counter.gt.inputs%nr_npa)stop'too many neutrals'
             npa%kind(npa%counter)=neut_type
             npa%v(npa%counter,:)=vn(:)
             npa%wght(npa%counter)=sum(states)/nlaunch*grid%dv/npa%npa_loop !![neutrals/s]
             npa%ipos(npa%counter,:)=ipos(:)
             npa%fpos(npa%counter,:)=ray(:)
             !$OMP END CRITICAL(col_rad_npa)
          endif
       endif
       return
    endif
    !! NORMAL START OF COLRAD
    !! ------------------ Get Matrix with rate coeffients ----------------- !!
    !! the rates  for  are computed for a thermal distribution!
    vnet_square=dot_product(vn-vrot,vn-vrot)          ![cm/s]
    eb=v_to_E*inputs%ab*vnet_square ![kev]

    !! DEUTERIUM (Impact excitation/ionization,charge exchange)
    ebi= floor(eb/tables%d_eb_qp)+1  
    if(ebi.ge.tables%nr_eb_qp)stop 'Eb out of range of qptable!'
    tii= floor(ti/tables%d_ti_qp)+1
    if(tii.ge.tables%nr_ti_qp)stop 'Ti out of range of qptable!'
    call table_interp(tables%qp(:,:,:,:),tables%d_eb_qp,tables%d_ti_qp &
         ,ebi,tii,eb,ti,qp)
    qp=qp*denp ![1/s]  

    !! IMPURITIES
    ebi= floor(eb/tables%d_eb_qi)+1   
    if(ebi.ge.tables%nr_eb_qi)stop 'Eb out of range of qitable!'
    tii= floor(ti/tables%d_ti_qi)+1
    if(tii.ge.tables%nr_ti_qi)stop 'Ti out of range of qitable!'
    call table_interp(tables%qi(:,:,:,:),tables%d_eb_qi,tables%d_ti_qi &
         ,ebi,tii,eb,ti,qi)
    qi=qi*deni ![1/s]  


    !! ELECTRONS
    ebi= floor(eb/tables%d_eb_qe)+1   
    if(ebi.ge.tables%nr_eb_qe)stop 'Eb out of range of qetable!'
    tei= floor(te/tables%d_te_qe)+1
    if(tei.ge.tables%nr_te_qe)stop 'Te out of range of qetable!'
    call table_interp(tables%qe(:,:,:,:),tables%d_eb_qe,tables%d_te_qe &
         ,ebi,tei,eb,te,qe)
    qe=qe*dene ![1/s] 
    !! - Write off-diagnonal elements (populating transitions) - !!
    matrix= tables%einstein(1:nlevs,1:nlevs)  &     
         +              qp(1:nlevs,1:nlevs)  &
         +              qi(1:nlevs,1:nlevs)  &
         +              qe(1:nlevs,1:nlevs)
    
 
    
    !! - Write diagonal elements (depopulating transitions) - !!
    do n=1,nlevs 	            
       matrix(n,n)=&
            - sum(tables%einstein(:,n)) &
            - sum(qp(:,n)) &
            - sum(qi(:,n)) &
            - sum(qe(:,n))
    enddo
    call eigen(matrix, eigvec, eigval)
    call matinv(eigvec, eigvec_inv)
    coef = matmul(eigvec_inv, states)!coeffs determined from states at t=0
    exp_eigval_dt = exp(eigval*dt)   ! to improve speed (used twice)
    states(:) = matmul(eigvec, coef * exp_eigval_dt)  ![neutrals/cm^3/s]!
    dens(:)   = matmul(eigvec,coef*(exp_eigval_dt-1.d0)/eigval)/nlaunch
    !$OMP CRITICAL(col_rad)
    result%neut_dens(ac(1),ac(2),ac(3),:,neut_type)= & 
         result%neut_dens(ac(1),ac(2),ac(3),:,neut_type)+dens(:)![neutrals/cm^3]
    result%velocity_counter(ac(1),ac(2),ac(3),neut_type)= & !!increase counter
         result%velocity_counter(ac(1),ac(2),ac(3),neut_type)+1
    counter=result%velocity_counter(ac(1),ac(2),ac(3),neut_type)
    vel_vec_red_lev=result%vel_vec_red_lev(ac(1),ac(2),ac(3),neut_type)
    if(mod(counter,vel_vec_red_lev).eq.0)then 
       counter=counter/vel_vec_red_lev
       if(counter.gt.nvelocities)then
          !! increase reduction level in storing velocity vectors
          result%vel_vec_red_lev(ac(1),ac(2),ac(3),neut_type)=vel_vec_red_lev*2
          !! move every second index to the first half of the array!
          do n=1,3
             velocity_dummy= &
                  result%velocity_vectors(ac(1),ac(2),ac(3),neut_type,every_second_index,n)
             result%velocity_vectors(ac(1),ac(2),ac(3),neut_type,first_half_indices,n)= &
                  velocity_dummy
          enddo
          !! don't save this one
       else 
          result%velocity_vectors(ac(1),ac(2),ac(3),neut_type,counter,:)=vn
       endif
    endif
    !$OMP END CRITICAL(col_rad)
    if(inputs%calc_birth.eq.1)then
       if(neut_type.le.3)then
          call bfield_interp(pos,b_abs,b_norm)
          ptch=dot_product(vn,b_norm)/sqrt(dot_product(vn,vn)) &
               *inputs%btipsign
          ipitch=int((ptch+1.)/(2./npitch_birth))+1
          dflux=(iflux-sum(states))*grid%dv/nlaunch !! [fast-ions/s]
          !$OMP CRITICAL(col_rad2)
          result%birth_dens(ac(1),ac(2),ac(3),neut_type,ipitch)= &
               result%birth_dens(ac(1),ac(2),ac(3),neut_type,ipitch) + dflux
          !$OMP END CRITICAL(col_rad2)
       endif
    endif
    !! -------------------- determine photon flux ------------------------ !!
    photons=dens(3)*tables%einstein(2,3) !! - [Ph/(s*cm^3)] - !!
  end subroutine colrad
  
  !***************************************************************************
  !-----------spectrum--------------------------------------------------------
  !***************************************************************************
  subroutine spectrum(vi,pos,photons,neut_type,wavout,intout)
    !!spectrum.pro computes the wavelengths of emitted photons
    real(double), dimension(:), intent(in) :: vi!!velocitiy of neutral [cm/s]
    integer             , intent(in)       :: neut_type!!type of neutral
    real(double)              , intent(in) :: photons !! photons from colrad
    real(double), dimension(3), intent(in) :: pos    !! mean position in cell
    real(double), dimension(n_stark), intent(out), optional :: intout!!intensity
    real(double), dimension(n_stark), intent(out), optional :: wavout !!wavlegth
    integer,dimension(3)       :: ac  !!actual cell
    real(double)               :: lambda_Doppler , cos_los_Efield, E
    real(double),dimension(n_stark) ::intens!!intensity vector
    real(double), dimension(n_stark)::wavel !!wavelength vector[A
    real(double), dimension(3) :: vp  !!unit vector of sight line
    real(double), dimension(3) :: vn  ! vi in m/s
    real(double), dimension(3) :: efield  !E-field (static + vxB)
    real(double), dimension(3) :: bfield,b_norm  !B-field
    real(double) :: b_abs
    integer                    :: i,ichan,bin,pisi !counter,wavelengths bins
    integer,parameter,dimension(n_stark)::stark_sign=+1*stark_sigma -1*stark_pi
               !sign in stark intensity formula:
               !- for Pi (linear), + for Sigma (circular)
    call get_ac(pos,ac)
    call bfield_interp(pos,b_abs,b_norm)
    loop_over_channels: do ichan=1,diag%nchan
       if(diag%los_wght(ac(1),ac(2),ac(3),ichan).le.0.)cycle loop_over_channels
       !! vector directing towards the optical head
       vp(1)=pos(1)-diag%xyzhead(ichan,1) 
       vp(2)=pos(2)-diag%xyzhead(ichan,2) 
       vp(3)=pos(3)-diag%xyzhead(ichan,3) 
       vp=vp/sqrt(dot_product(vp,vp))

       ! Calculate Doppler shift
       vn=vi*0.01d0 ! [m/s]
       lambda_Doppler = lambda0*(1.d0 + dot_product(vn,vp)/c0)
       !! Calculate Stark Splitting
       ! Calcualate E-field
       bfield(:) = b_norm*b_abs
       efield(:) = 0.d0
       efield(1) = efield(1) +  vn(2)*bfield(3) - vn(3)*bfield(2)
       efield(2) = efield(2) - (vn(1)*bfield(3) - vn(3)*bfield(1))
       efield(3) = efield(3) +  vn(1)*bfield(2) - vn(2)*bfield(1)
       E=sqrt(dot_product(efield,efield))
       !Stark Splitting
       wavel =  lambda_Doppler + E * stark_wavel ![A]
       !Intensities of stark components
       if (E .eq. 0.d0) then 
          cos_los_Efield = 0.d0 
       else 
          cos_los_Efield = dot_product(vp,efield) / E
       endif
       intens = stark_intens*(1.d0+ stark_sign* cos_los_Efield**2.d0)
       !! --- E.g. mirrors may change the pi to sigma intensity ratio  --- !!
       where (stark_sigma .eq. 1)
          intens = intens * diag%sigma_pi(ichan)
       endwhere
       !! --- normalize and multiply with photon density from colrad --- !!
       intens      = intens/sum(intens)*photons 
       if(present(wavout))then
          wavout=wavel
          intout=intens
          return
       endif
       !! ---------------------- Store spectra ---------------------- !!
       do i=1,n_stark                                                       
          bin=floor((wavel(i)-diag%lambdamin)/diag%dlambda)+1
          if (bin.lt.1)            bin = 1                        
          if (bin.gt.diag%nlambda) bin = diag%nlambda  
          pisi=1
          if(stark_sigma(i).eq.1)pisi=2
          !$OMP CRITICAL(spec_trum)  
          result%spectra(pisi,bin,ichan,neut_type)= &
               result%spectra(pisi,bin,ichan,neut_type) &
               +intens(i)*diag%los_wght(ac(1),ac(2),ac(3),ichan)
          !$OMP END CRITICAL(spec_trum)
       enddo
    enddo loop_over_channels
  end subroutine spectrum
  
  
  !*****************************************************************************
  !------------track------------------------------------------------------------
  !*****************************************************************************
  subroutine track(vin, posin, pos_arr,dt_arr,ntrack,icell)
    !!track computes the path of a neutral through a the FIDAcode grid 
    real(double), dimension(:)  , intent(in)   :: posin  ! initial position
    real(double), dimension(:)  , intent(in)   :: vin  ! velocitiy
    integer               , intent(out)        :: ntrack! number of cells
    real(double), dimension(:)  , intent(out)  :: dt_arr! time per cell
    integer,dimension(:,:),intent(out),optional:: icell ! cell indices
    real(double), dimension(:,:), intent(out)  :: pos_arr! mean position in cell
    integer                    :: cc    !!step number along the track
    integer,dimension(3)       :: p,l    !!indices of the cells
    real(double), dimension(3) :: dt_min !!time to cell boundary
    real(double)               :: dt     !!min time to cell boundary
    real(double), dimension(3) :: vn !!velocitiy that can be changed
    real(double), dimension(3) :: ipos !!position of ray  
    dt_arr=0.d0 ; pos_arr=0.d0 ; ntrack=0
    vn(:)=vin(:) ;  ipos(:)=posin(:)
    !! define actual cell
    call get_ac(ipos,p)
    !! Fudge zero velocity components to avoid overflow error
    where (vn(:).eq.0) vn(:) = 0.001   
    !! Start tracking routine
    if(present(icell))then
       icell=0
       icell(:,1)=p(:)
    endif
    cc=1 
    !!loop along track of neutral
    tracking: do while(cc.lt.(grid%ntrack))
       l(:)=p(:)
       where(vn(:).gt.0.d0) l(:)=p(:)+1 
       if ( l(1).gt.grid%nx.or.& 
            l(2).gt.grid%ny.or.&
            l(3).gt.grid%nz) exit tracking  
       !time needed to go to next cell
       dt_min(1)=(grid%xx(l(1))-ipos(1))/vn(1)
       dt_min(2)=(grid%yy(l(2))-ipos(2))/vn(2)
       dt_min(3)=(grid%zz(l(3))-ipos(3))/vn(3)
       minpos=minloc(dt_min) 
       dt=dt_min(minpos(1))
       pos_arr(:,cc) = ipos(:) + vn(:)*dt*0.5  !! mean postion in cell
       ipos(:)     = ipos(:) + vn(:)*dt
       if (vn(minpos(1)).gt.0.d0)  then 
          p(minpos(1))=p(minpos(1))+1
       else
          p(minpos(1))=p(minpos(1))-1
       endif
       if (any(p.le.0))exit tracking    
       dt_arr(cc)=dt
       if(present(icell))icell(:,cc+1)=p(:)
       cc=cc+1
    enddo tracking
    ntrack=cc-1
  end subroutine track
  


  
 !*****************************************************************************
 !------------track2-----------------------------------------------------------
 !*****************************************************************************
  subroutine track2(vin, posin, pos_arr,dt_arr,ntrack)
    !!track computes the path of a neutral through a the FIDAcode grid 
    real(double), dimension(3)  , intent(in)   :: posin  ! initial position
    real(double), dimension(3)  , intent(in)   :: vin    ! velocitiy
    real(double), dimension(:,:), intent(out)  :: pos_arr! positions along track
    real(double), dimension(:)  , intent(out)  :: dt_arr !time
    integer(long)               , intent(out)  :: ntrack     ! counter
    real(double)                               :: vabs   ! velocity
    real(double)                               :: dt
    real(double), dimension(3)                 :: tmin_arr ! time to grid bounda
    integer(long)                              :: i
    pos_arr=0.d0
    dt_arr =0.d0
    ntrack=0
    !! check if velocity is not zero
    vabs=sqrt(dot_product(vin,vin))
    if(vabs.eq.0)return
    dt=grid%dl/vabs
    pos_arr(:,1)=posin(:)
    tracking:do i=2,grid%ntrack
       dt_arr(i-1)=dt
       pos_arr(:,i)=pos_arr(:,i-1)+vin(:)*dt
       if((pos_arr(1,i).lt.grid%xmin).or. &
            (pos_arr(1,i).gt.grid%xmax))exit tracking
       if((pos_arr(2,i).lt.grid%ymin).or. &
            (pos_arr(2,i).gt.grid%ymax))exit tracking
       if((pos_arr(3,i).lt.grid%zmin).or. &
            (pos_arr(3,i).gt.grid%zmax))exit tracking
    enddo tracking
    if(i.eq.grid%ntrack)print*, 'grid%ntrack not large enough!'
    i=i-1
    !! the last position has a variable length that is determined by
    !! the grid boundary. Hence, find length to the grid boundary and 
    !! adjust the position and time
    if(vin(1).gt.0)then
       tmin_arr(1)=(grid%xmax-pos_arr(1,i))/vin(1)
    else
       tmin_arr(1)=(grid%xmin-pos_arr(1,i))/vin(1)
    endif
    if(vin(2).gt.0)then
       tmin_arr(2)=(grid%ymax-pos_arr(2,i))/vin(2)
    else
       tmin_arr(2)=(grid%ymin-pos_arr(2,i))/vin(2)
    endif
    if(vin(3).gt.0)then
       tmin_arr(3)=(grid%zmax-pos_arr(3,i))/vin(3)
    else
       tmin_arr(3)=(grid%zmin-pos_arr(3,i))/vin(3)
    endif
    dt_arr(i)=minval(tmin_arr)
    ntrack=i
    !! if there is only one position, check if this position is inside the 
    !! simultion grid! Otherwise, set ntrack to zero!
    if(ntrack.eq.1)then
       ntrack=0
       if((posin(1).lt.grid%xmin).or. &
            (posin(1).gt.grid%xmax))return
       if((posin(2).lt.grid%xmin).or. &
            (posin(2).gt.grid%xmax))return  
       if((posin(3).lt.grid%xmin).or. &
            (posin(3).gt.grid%xmax))return  
       ntrack=1
    endif
    !! Add 0.5*dl to every position!
    do i=1,ntrack
       pos_arr(:,i)=pos_arr(:,i)+0.5*vin(:)*dt_arr(i)
    enddo
  end subroutine track2


  
  function cross_product(a,b)
    real(double), dimension(3), intent(in)   :: a,b
    real(double), dimension(3) :: cross_product
    cross_product(1)=  a(2)*b(3)-a(3)*b(2)
    cross_product(2)=-(a(1)*b(3)-a(3)*b(1))
    cross_product(3)=  a(1)*b(2)-a(2)*b(1)
  end function cross_product
  
  subroutine orbit(ipos,iv,mass,charge,npos,pos_arr,dt)
    real(double), dimension(3)  , intent(in)   :: ipos ! initial position
    real(double), dimension(3)  , intent(in)   :: iv  ! velocitiy
    real(double)                , intent(in)   :: mass, charge
    integer(long)               , intent(in)   :: npos
    real(double), dimension(:,:), intent(out)  :: pos_arr
    real(double)                , intent(out)  :: dt

    real(double)               :: dt_sub, time
    real(double), dimension(3) :: v,v2
    integer(long)              :: i,j,nsub
    real(double)               :: b_abs,factor
    real(double), dimension(3) :: b_norm
    real(double), dimension(3) :: rk1,rk2,rk3,rk4
    real(double),   dimension(3,3)  :: eigvec
    real(double),   dimension(3,3)  :: B
    real(double),   dimension(3)    :: eigval
    real(double),   dimension(3,3)  :: eigvec_inv
    
    real(double),   dimension(3)        :: coef
    real(double),   dimension(3)        :: exp_eigval_dt 
    

    dt = grid%dl/sqrt(dot_product(iv,iv))
    nsub=grid%dl/0.01
    if(nsub.lt.2)nsub=2
    dt_sub=dt/nsub
    pos_arr(:,1)=ipos
    v(:)=iv(:)
    time = 0.0
    do i=1,npos-1
       call bfield_interp(pos_arr(:,i),b_abs,b_norm)
       if(b_abs.le.0)exit
       factor=charge/mass*b_abs
       pos_arr(:,i+1)=pos_arr(:,i)+v*0.5*dt

      ! b_norm=b_norm*factor

       B(1,:)=(/      0.d0, b_norm(3),-b_norm(2)/)
       B(2,:)=(/-b_norm(3),      0.d0, b_norm(1)/)
       B(3,:)=(/ b_norm(2),-b_norm(1),      0.d0/)
   
       if(b_abs.le.0)exit
       factor=charge/mass*b_abs
       !! a=charge* VxB/mass
       !! vnew=v+a*dt
       !! pos=pos+v*dt
       !! RUNGE KUTTA:
       do j=1,nsub
          rk1=cross_product(b_norm,v)*factor
          rk2=cross_product(b_norm,rk1(:)*dt_sub*0.5+v)*factor
          rk3=cross_product(b_norm,rk2(:)*dt_sub*0.5+v)*factor
          rk4=cross_product(b_norm,rk3(:)*dt_sub+v)*factor
          v(:) = v(:)+dt_sub*(rk1(:)+2.d0*(rk2(:)+rk3(:))+rk4(:))/6.d0
       enddo
       print*, v
       print*, '====='
       pos_arr(:,i+1)=pos_arr(:,i)+v*0.5*dt
       time=time+dt
    enddo
    print*, time
  end subroutine orbit






 !*****************************************************************************
 !----------- get_nlauch ------------------------------------------------------
 !*****************************************************************************
  subroutine get_nlaunch(nr_markers,papprox,papprox_tot,nlaunch)
    !! routine to define the number of MC particles started in one cell
    integer                       , intent(in)    :: nr_markers
    real(double), dimension(:,:,:), intent(in)    :: papprox
    real(double)                  , intent(in)    :: papprox_tot
    real(double), dimension(:,:,:), intent(out)   :: nlaunch  
    integer  :: i
    do i=1,1000
       nlaunch(:,:,:)=papprox(:,:,:)/papprox_tot*nr_markers*(1.+i*0.01)
       if(sum(nlaunch).gt.nr_markers)exit
    enddo
    nlaunch=ceiling(nlaunch)
  end subroutine get_nlaunch


  

  !*****************************************************************************
  !-----------ndmc (NBI)--------------------------------------------------------
  !*****************************************************************************
  subroutine ndmc
    integer                                :: indmc     !! counter for markers
    integer                                :: type      !! full half third En
    real(double)                           :: nlaunch   !! nr. of markers
    real(double)                           :: nneutrals !! # NBI particles 
    real(double), dimension(3)             :: vnbi      !! velocities(full..)
    real(double), dimension(3)             :: ipos      !! initial position
    !!Tracking routine output
    integer                                :: jj     !! counter for track
    integer                                :: ntrack !! number of track position
    real(double), dimension(3,grid%ntrack) :: pos_arr!! position along track
    real(double), dimension(grid%ntrack)   :: dt_arr !! time per track step
    !!collisional radiative model and spectrum calculation
    real(double), dimension(nlevs)         :: states
    real(double)                           :: photons,dt
    print*,'    # of markers: ',inputs%nr_ndmc
    !! ------------- calculate nr. of injected neutrals ---------------- !!
    !! # of injected neutrals = NBI power/energy_per_particle
    nneutrals=1.d6*nbi%pinj/ (1.d3*nbi%einj*e0 &
         *( nbi%species_mix(1)      &
         +  nbi%species_mix(2)/2.d0 &
            +  nbi%species_mix(3)/3.d0 ) )
    !! ------------------ loop over the markers ------------------------ !!
    nlaunch=real(inputs%nr_ndmc)
    !$OMP PARALLEL DO private(type,indmc,vnbi,ipos,pos_arr,dt_arr, &
    !$OMP& ntrack,states,jj,photons)
    energy_fractions: do type=1,3
       !! (type = 1: full energy, =2: half energy, =3: third energy
       loop_over_markers: do indmc=1,inputs%nr_ndmc
          call mc_nbi(vnbi(:),type,ipos(:))
          if(ipos(1).eq.-1)cycle loop_over_markers
          call track2(vnbi,ipos, pos_arr, dt_arr, ntrack)
          !! --------- solve collisional radiative model along track ---- !!
          states=0.d0
          states(1)=nneutrals*nbi%species_mix(type)/grid%dv 
          loop_along_track: do jj=1,ntrack
             call colrad(pos_arr(:,jj),vnbi,dt_arr(jj),states &
                  ,photons,type,nlaunch)
             if(photons.gt.0.d0) then
                if(inputs%calc_spec.eq.1) &
                     call spectrum(vnbi(:),pos_arr(:,jj),photons,type)
             endif
          enddo loop_along_track
       enddo loop_over_markers
    enddo energy_fractions
    !$OMP END PARALLEL DO
    if(nbi_outside.gt.0)then
       print*, 'Percent of markers outside the grid: ' &
            ,100.*nbi_outside/(3.*inputs%nr_ndmc)
       if(sum(result%neut_dens).eq.0)stop 'Beam does not intersect the grid!'
    endif
  end subroutine ndmc
  


  !*****************************************************************************
  !-------------- Direct charge exchange calculation---------------------------
  !*****************************************************************************
  subroutine dcx
    integer                                :: i,j,k   !! indices of cells
    real(double), dimension(3)             :: randomu   
    integer                                :: idcx    !! counter
    real(double), dimension(3)             :: ipos      !! start position
    real(double), dimension(3)             :: vhalo   !! velocity bulk plasma 
    real(double)                           :: rho
    !! Determination of the CX probability
    real(double), dimension(nlevs)         :: denn    !!  neutral dens (n=1-4)
    real(double) :: denf,denp
    real(double), dimension(nlevs)         :: prob    !!  Prob. for CX 
    integer                                :: counter !!of velocity vectors
    integer                                :: vel_vec_red_lev !! reduction level
    integer                                :: type    !! type of NBI neutral
    real(double), dimension(3)             :: vnbi    !! Velocity of NBIneutrals
    integer                                :: in      !! index neut rates
    real(double), dimension(nlevs)         :: rates   !! Rate coefficiants forCX
    !! Collisiional radiative model along track
    real(double), dimension(nlevs)         :: states  !! Density of n-states
    integer                                :: ntrack
    real(double), dimension(grid%ntrack)   :: dt_arr  !! time per track step
    real(double), dimension(3,grid%ntrack) :: pos_arr !! positions along track
    integer                                :: jj      !! counter along track
    real(double)                           :: photons !! photon flux
    real(double), dimension(grid%nx,grid%ny,grid%nz)::papprox !!approx.density
    real(double)                           :: papprox_tot
    real(double), dimension(grid%nx,grid%ny,grid%nz)::nlaunch
   
    papprox=0.d0
    papprox_tot=0.d0
    !! ------------- calculate papprox needed for guess of nlaunch --------!!
    do k=1,grid%Nz 
       do j=1,grid%Ny 
          do i=1,grid%Nx 
             ipos=(/grid%xxc(i),grid%yyc(j),grid%zzc(k)/)
             call prof_interp(ipos,denp)
             call denf_interp(ipos,denf)
             if(denp-denf.gt.0)then
                papprox(i,j,k)=    (sum(result%neut_dens(i,j,k,:,nbif_type))  &
                     +              sum(result%neut_dens(i,j,k,:,nbih_type))  &
                     +              sum(result%neut_dens(i,j,k,:,nbit_type))) &
                  *(denp-denf)
                call rho_interp(ipos,rho)
                if(rho.lt.1.1) papprox_tot=papprox_tot+papprox(i,j,k)  
             endif
          enddo
       enddo
    enddo
    call get_nlaunch(inputs%nr_dcx,papprox,papprox_tot,nlaunch)
    ! Loop through all of the cells
    print*,'    # of markers: ',int(sum(nlaunch))
    !$OMP PARALLEL DO private(i,j,k,idcx,randomu,ipos,vhalo,pos_arr,dt_arr, &
    !$OMP& ntrack,prob,type,counter,vel_vec_red_lev,vnbi, &
    !$OMP& denn,rates,denp,denf,states,jj,photons)
    loop_along_z: do k = 1, grid%Nz
       loop_along_y: do j = 1, grid%Ny
          loop_along_x: do i = 1, grid%Nx
             !! ------------ loop over the markers ---------------------- !!
             loop_over_dcx: do idcx=1,int(nlaunch(i,j,k))
                !! ---------------- calculate ri,vi and track -------------!!
                call randu(randomu)  
                ipos(1)=grid%xx(i)+ grid%dx*randomu(1)
                ipos(2)=grid%yy(j)+ grid%dy*randomu(2) 
                ipos(3)=grid%zz(k)+ grid%dz*randomu(3)
                call mc_halo(ipos,vhalo(:))
                if(sum(vhalo).eq.0)cycle loop_over_dcx 
                call track2(vhalo,ipos, pos_arr, dt_arr, ntrack)
                !! ---------------- calculate CX probability ------------- !!
                prob=0.d0
                energy_fractions: do type=1,3
                   call randu(randomu)  
                   !! (type = 1: full energy, =2: half energy, =3: third energy
                   !! get velocity vector of NBI neutrals ------ !!
                   counter=result%velocity_counter(i,j,k,type)
                   vel_vec_red_lev=result%vel_vec_red_lev(i,j,k,type)
                   counter=int(counter/vel_vec_red_lev)-1
                   jj=int(counter*randomu(1)+0.5)+1
                   vnbi=result%velocity_vectors(i,j,k,type,jj,:)
                   !! get density of NBI neutrals ------ !!
                   denn(:)=result%neut_dens(i,j,k,:,type)
                   !! calculate neutralization probability ------ !!
                   call neut_rate(denn,vhalo,vnbi,rates)
                   prob=prob + rates
                enddo energy_fractions
 
                if(sum(prob).le.0.)cycle loop_over_dcx
                !! --------- solve collisional radiative model along track-!!
                call prof_interp(ipos,denp)
                call denf_interp(ipos,denf)
                states=prob*(denp-denf)
                loop_along_track: do jj=1,ntrack
                   call colrad(pos_arr(:,jj),vhalo(:),dt_arr(jj),states &
                        ,photons,halo_type,nlaunch(i,j,k))
                   if(photons.le.0.d0)cycle loop_over_dcx 
                   if(inputs%calc_spec.eq.1)call spectrum(vhalo(:) &
                        ,pos_arr(:,jj),photons,halo_type)
                enddo loop_along_track
             enddo loop_over_dcx
          enddo loop_along_x
       enddo loop_along_y
    enddo loop_along_z
    !$OMP END PARALLEL DO
  end subroutine dcx
  
  !*****************************************************************************
  !-------------------------- halo -------------------------------------------
  !*****************************************************************************
  subroutine halo
    integer                                :: i,j,k !indices of cells 
    integer                                :: ihalo !! counter
    real(double), dimension(3)             :: randomu   
    real(double), dimension(3)             :: ipos    !! start position
    real(double), dimension(3)             :: vihalo!! velocity bulk plasma ion
    real(double)                           :: rho
    !! Determination of the CX probability
    real(double), dimension(nlevs)         :: denn    !! neutral dens (n=1-4)
    real(double)                           :: denf    !! fast-ion density
    real(double)                           :: denp    !! Proton density
    real(double), dimension(nlevs)         :: prob    !! Prob. for CX 
    integer                                :: counter !!of velocity vectors
    integer                                :: vel_vec_red_lev !! reduction level
    real(double), dimension(3)             :: vnhalo  !! v of halo neutral
    integer                                :: in      !! index over halo neutral
    !! Collisiional radiative model along track
    real(double), dimension(nlevs)         :: states  ! Density of n-states
    integer                                :: ntrack
    real(double), dimension(grid%ntrack)   :: dt_arr  !! time per track step
    real(double), dimension(3,grid%ntrack) :: pos_arr !! positions along track
    integer                                :: jj       !! counter along track
    real(double)                           :: photons  !! photon flux
    real(double), dimension(grid%nx,grid%ny,grid%nz)::papprox,nlaunch !! approx. density
    real(double)                           :: papprox_tot 
    !! Halo iteration
    integer                                :: hh !! counters
    real(double)                           :: dcx_dens, halo_iteration_dens
    integer  :: s1type  ! halo iteration
    integer  :: s2type  ! halo iteration
    real(double), dimension(nvelocities/2)  :: velocity_dummy
    integer                                 :: ncounters_s2,nn !! save velocoity vectors

    s1type=afida_type
    s2type=pfida_type
    dcx_dens=sum(result%neut_dens(:,:,:,:,halo_type))
    if(dcx_dens.eq.0)stop 'the denisty of DCX-neutrals is too small!'

    result%neut_dens(:,:,:,:,s1type)          = result%neut_dens(:,:,:,:,halo_type)
    result%velocity_vectors(:,:,:,s1type,:,:) = result%velocity_vectors(:,:,:,halo_type,:,:)
    result%velocity_counter(:,:,:,s1type)     = result%velocity_counter(:,:,:,halo_type)
    result%vel_vec_red_lev(:,:,:,s1type)      = result%vel_vec_red_lev(:,:,:,halo_type)
    iterations: do hh=1,20
       !! ------------- calculate papprox needed for guess of nlaunch --------!!
       papprox=0.d0
       papprox_tot=0.d0
       do k=1,grid%Nz 
          do j=1,grid%Ny 
             do i=1,grid%Nx   
                ipos=(/grid%xxc(i),grid%yyc(j),grid%zzc(k)/)
                call prof_interp(ipos,denp)
                papprox(i,j,k)=sum(result%neut_dens(i,j,k,:,s1type))*denp
                call rho_interp(ipos,rho)
                if(rho.lt.1.1) papprox_tot=papprox_tot+papprox(i,j,k) 
             enddo
          enddo
       enddo
       call get_nlaunch(inputs%nr_halo,papprox,papprox_tot,nlaunch)
       print*, '    # of markers: ' ,int(sum(nlaunch))
       
       !$OMP PARALLEL DO private(i,j,k,ihalo,randomu,ipos,denp,denf, &
       !$OMP& vihalo,pos_arr,dt_arr,ntrack,prob,counter,vel_vec_red_lev, &
       !$OMP& denn,vnhalo,states,jj,photons)
       loop_along_z: do k = 1, grid%Nz
          loop_along_y: do j = 1, grid%Ny
             loop_along_x: do i = 1, grid%Nx
                !! ------------- loop over the markers ---------------------- !!
                loop_over_halos: do ihalo=1,int(nlaunch(i,j,k))
                   !! ---------------- calculate ri,vhalo and track ----------!!
                   call randu(randomu)  
                   ipos(1)=grid%xx(i)+ grid%dx*randomu(1)
                   ipos(2)=grid%yy(j)+ grid%dy*randomu(2) 
                   ipos(3)=grid%zz(k)+ grid%dz*randomu(3)
                   call prof_interp(ipos,denp)
                   call denf_interp(ipos,denf)
                   if(denp.le.0)cycle loop_over_halos

                   call mc_halo(ipos, vihalo(:))
                   call track2(vihalo,ipos, pos_arr, dt_arr, ntrack)
                   !! ---------------- calculate CX probability --------------!!
                   prob=0.d0
                   !! get density of Halo neutrals of previous generation ------ !!
                   denn(:)=result%neut_dens(i,j,k,:,s1type)
                   !! there is no second iteration in the fast-ion simulation. 
                   !! So add something like the fi-contribution to denn
                   !! because fast-neutrals can also generate halos!
                   if(denf/denp.lt.0.99.or.denf/denp.gt.1.01)then 
                      denn(:)=denn(:)*(1.d0+denf/(denp-denf))
                   else
                      print*,'denf and denp almost the same (in halo routine)'
                   endif
                   !! get velocity vector of Halo neutrals of previous generation
                   counter=result%velocity_counter(i,j,k,s1type)
                   vel_vec_red_lev=result%vel_vec_red_lev(i,j,k,s1type)
                   counter=int(counter/vel_vec_red_lev)-1
                   jj=int(counter*randomu(1)+0.5)+1
                   vnhalo=result%velocity_vectors(i,j,k,s1type,jj,:)
                   !! calculate neutralization probability ------ !!
                   call neut_rate(denn,vihalo,vnhalo,prob)

                   if(sum(prob).le.0.) cycle loop_over_halos
                   !! --------- solve collisional radiative model along track-!!
                   states=prob*(denp-denf)
                   loop_along_track: do jj=1,ntrack
                      call colrad(pos_arr(:,jj),vihalo(:),dt_arr(jj),states &
                           ,photons,s2type,nlaunch(i,j,k))
                      if(photons.le.0.d0)cycle loop_over_halos 
                      if(inputs%calc_spec.eq.1) &
                           call spectrum(vihalo,pos_arr(:,jj),photons,halo_type)
                   enddo loop_along_track
                enddo loop_over_halos
             enddo loop_along_x
          enddo loop_along_y
       enddo loop_along_z
       !$OMP END PARALLEL DO
       !! ADD this generations density to the total halo density
       halo_iteration_dens=sum(result%neut_dens(:,:,:,:,s2type))
       result%neut_dens(:,:,:,:,halo_type)=result%neut_dens(:,:,:,:,halo_type) &
            + result%neut_dens(:,:,:,:,s2type)
       !! add this generations velocity vectors in halo_type 
       !! this is a bit more difficult as the vel_vec_red_level has to be taken into account
       do k = 1, grid%Nz
          do j = 1, grid%Ny
             do i = 1, grid%Nx
                vel_vec_red_lev=result%vel_vec_red_lev(i,j,k,halo_type)
                ncounters_s2=int(result%velocity_counter(i,j,k,s2type)/vel_vec_red_lev)
                do jj=1,ncounters_s2 !! add counter_s2 velocity vectors
                   result%velocity_counter(i,j,k,halo_type)= & !! increase halo_counter
                        result%velocity_counter(i,j,k,halo_type)+1
                   counter=result%velocity_counter(i,j,k,halo_type)
                   if(mod(counter,vel_vec_red_lev).eq.0)then 
                      if(counter.gt.nvelocities)then
                         !! increase reduction-level in storing velocity vectors
                         vel_vec_red_lev = vel_vec_red_lev*2
                         result%vel_vec_red_lev(i,j,k,halo_type)=vel_vec_red_lev
                         !! move every second index to the first half of the array!
                         do nn=1,3
                            velocity_dummy= &
                                 result%velocity_vectors(i,j,k,halo_type,every_second_index,nn)
                            result%velocity_vectors(i,j,k,halo_type,first_half_indices,nn)= &
                                 velocity_dummy
                         enddo
                         !! don't save this one
                      else 
                         result%velocity_vectors(i,j,k,halo_type,counter,:)= &
                              result%velocity_vectors(i,j,k,s2type,jj,:)
                      endif
                   endif
                enddo
             enddo
          enddo
       enddo
       !! reset storage arrays
       result%neut_dens(:,:,:,:,s1type)= result%neut_dens(:,:,:,:,s2type)
       result%neut_dens(:,:,:,:,s2type)= 0.
       result%velocity_vectors(:,:,:,s1type,:,:) = result%velocity_vectors(:,:,:,s2type,:,:)
       result%velocity_counter(:,:,:,s1type)     = result%velocity_counter(:,:,:,s2type)
       result%velocity_counter(:,:,:,s2type)     = 0
       result%vel_vec_red_lev(:,:,:,s1type)      = result%vel_vec_red_lev(:,:,:,s2type)
       result%vel_vec_red_lev(:,:,:,s2type)      = 1
       !! check if to exit the iteration 
       if(halo_iteration_dens/dcx_dens.gt.1)exit iterations
       inputs%nr_halo=inputs%nr_dcx*halo_iteration_dens/dcx_dens
       if(inputs%nr_halo.lt.inputs%nr_dcx*0.01)exit iterations
    enddo iterations 
    !! set the neutral density in s1type(fida_type) and s2type (dummy) to 0!
    result%neut_dens(:,:,:,:,s1type) = 0.d0
    result%neut_dens(:,:,:,:,s2type) = 0.d0
    result%velocity_counter(:,:,:,s1type)     = 0
    result%velocity_counter(:,:,:,s2type)     = 0
    result%vel_vec_red_lev(:,:,:,s1type)      = 1
    result%vel_vec_red_lev(:,:,:,s2type)      = 1
  end subroutine halo
  !*****************************************************************************
  !-----------FIDA simulation---------------------------------------------------
  !*****************************************************************************
  subroutine fida      
    integer                               :: i,j,k   !! indices  x,y,z  of cells
    integer                               :: iion
    real(double), dimension(3)            :: ipos    !! start position and index
    integer,      dimension(3)            :: iind    !! index of start cell
    real(double), dimension(3)            :: vi      !! velocity of fast ions
    real(double), dimension(:,:),allocatable  :: vi_arr  !! velocity array
    integer,dimension(3)                  :: ac      !! new actual cell
    real(double)                          :: rho
    !! Determination of the CX probability
    real(double), dimension(nlevs)        :: denn           !!  neutral dens (n=1-4)
    real(double)                          :: vtor           !! toroidal rotation
    real(double)                          :: ti,te          !! Ion/electron temperature
    real(double)                          :: denp,dene,deni !! P/impurity/electron densit
    real(double)                          :: denf    !! fast-ion density
    real(double), dimension(nlevs)        :: prob    !! Prob. for CX 
    !! Determine velocity vectors of neutrals
    integer                               :: counter !!of velocity vectors
    integer                               :: vel_vec_red_lev !! reduction level
    integer                               :: type    !! type of NBI neutral 
    real(double), dimension(3)            :: randomu   
    real(double), dimension(3)            :: vn      !!  velocity vector of neutrals  
    real(double), dimension(nlevs)        :: rates   !! Rate coefficiants for CX
    !! Collisiional radiative model along track
    real(double), dimension(nlevs)        :: states  ! Density of n-states
    integer                               :: ntrack
    real(double), dimension(grid%ntrack)  :: dt_arr  !! time per track step
    real(double), dimension(3,grid%ntrack):: pos_arr !! positions along track
    integer                               :: jj      !! counter along track
    real(double)                          :: photons !! photon flux 
    real(double), dimension(grid%nx,grid%ny,grid%nz)::papprox,nlaunch !! approx. density
    real(double)                          :: nlaunch2   
    real(double)                          :: vi_abs             !! (for NPA)
    real(double), dimension(3)            :: ray,ddet,hit_pos   !! ray towards NPA
    real(double)                          :: papprox_tot 
    integer                               :: inpa    
    real(double)                          :: alpha !! angle relative to detector LOS
    integer                               :: cnt,maxcnt !! counters to show progress of simualtion
    !! ------------- calculate papprox needed for guess of nlaunch --------!!
    papprox=0.d0
    papprox_tot=0.d0
    do k=1,grid%Nz 
       do j=1,grid%Ny 
          do i=1,grid%Nx
             ipos=(/grid%xxc(i),grid%yyc(j),grid%zzc(k)/)
             if (inputs%npa.eq.1) then
                ray=diag%xyzhead(1,:)-ipos
                alpha=acos(dot_product(npa%los,ray) &
                     /(npa%dlos*sqrt(dot_product(ray,ray))))
                if (alpha.gt.npa%opening_angle*3.) cycle
             endif
             call denf_interp(ipos,denf)
             papprox(i,j,k)=(sum(result%neut_dens(i,j,k,:,nbif_type))  &
                  +          sum(result%neut_dens(i,j,k,:,nbih_type))  &
                  +          sum(result%neut_dens(i,j,k,:,nbit_type))  &
                  +          sum(result%neut_dens(i,j,k,:,halo_type)))  &
                  *          denf
             call rho_interp(ipos,rho)
             if(rho.lt.1.1) papprox_tot=papprox_tot+papprox(i,j,k) 
          enddo
       enddo
    enddo
    call get_nlaunch(inputs%nr_fida,papprox,papprox_tot,nlaunch)
    print*,'    # of markers: ',int(sum(nlaunch))
    if(inputs%npa.eq.1)print*,'    # npa_loop:   ',int(npa%npa_loop)
    maxcnt=grid%Nz*grid%Ny*grid%Nz
    cnt=0
    !$OMP PARALLEL DO private(i,j,k,inpa,vi_arr,nlaunch2,iion,ac,ipos,vi, &
    !$OMP& vi_abs,alpha,iind,ray,hit_pos,ddet, &
    !$OMP& denp,ti,vtor,te,dene,deni, &
    !$OMP& pos_arr,dt_arr,ntrack,prob,type,counter,vel_vec_rel_lev,randomu,vn, &
    !$OMP& denn,rates,denf,states,jj,photons,cnt)
    loop_along_z: do k = 1, grid%Nz
       loop_along_y: do j = 1, grid%Ny
          loop_along_x: do i = 1, grid%Nx
             if(nlaunch(i,j,k).le.0)cycle
             ipos=(/grid%xxc(i),grid%yyc(j),grid%zzc(k)/)
             call denf_interp(ipos,denf)
             call prof_interp(ipos,denp,ti,vtor,te,dene,deni)
             if(allocated(vi_arr)) deallocate(vi_arr)
             allocate(vi_arr(int(nlaunch(i,j,k)+grid%ntrack),3))
             iind=(/i,j,k/)
             !! ------------- loop over the markers ---------------------- !!
             npa_loop: do inpa=1,int(npa%npa_loop)
                call mc_fastion2(iind,ipos,nlaunch(i,j,k),vi_arr(:,:),nlaunch2)
                loop_over_fast_ions: do iion=1,int(nlaunch2)
                   vi(:)=vi_arr(iion,:)
                   if(sum(vi).eq.0) cycle loop_over_fast_ions
                   !! -------- check if particle hits the NPA detector ---- !!
                   if(inputs%npa.eq.1)then  
                      vi_abs=sqrt(dot_product(vi,vi))
                      alpha=acos(dot_product(npa%los,vi(:))/(vi_abs*npa%dlos))
                      if (alpha.gt.npa%opening_angle) cycle loop_over_fast_ions
                   endif
                   call mc_start(iind(:), vi(:), ipos(:))
                   !! -------- check if track ends at the NPA detector ---- !!
                   if(inputs%npa.eq.1)then 
                      !! check if track ends at the NPA detector
                      ray=ipos-diag%xyzhead(1,:)
                      hit_pos(:)=ipos(:)+vi(:)/vi_abs*sqrt(dot_product(ray,ray))
                      ddet=hit_pos-diag%xyzhead(1,:)
                      if(sqrt(dot_product(ddet,ddet))&
                           .gt.npa%size(1)) cycle loop_over_fast_ions
                   endif
                   call track2(vi,ipos, pos_arr, dt_arr, ntrack)
                   if(ntrack.eq.0)cycle loop_over_fast_ions
                   !! ---------------- calculate CX probability --------------!!
                   call get_ac(ipos,ac) !! new cell maybe due to gyro orbit!
                   prob=0.d0
                   neutral_types: do type=1,4
                      !! = 1: full, =2: half, =3: one third, =4: halo
                      !! get velocity vectors of those neutrals ------ !!
                      call randu(randomu)  
                      counter=result%velocity_counter(ac(1),ac(2),ac(3),type)
                      vel_vec_red_lev=result%vel_vec_red_lev(ac(1),ac(2),ac(3),type)
                      counter=int(counter/vel_vec_red_lev)-1
                      jj=int(counter*randomu(1)+0.5)+1
                      vn=result%velocity_vectors(ac(1),ac(2),ac(3),type,jj,:)
                      !! get density of those neutrals ------ !!
                      denn(:)=result%neut_dens(ac(1),ac(2),ac(3),:,type)
                      if(type.eq.halo_type)then
                         !! there is no Impurity ion halo calculation. Therefore,  
                         !! assume that there were no impurities
                         if(denp.gt.0)denn(:)=denn(:)*dene/denp
                      endif
                      !! calculate neutralization probability ------ !!
                      call neut_rate(denn,vi,vn,rates)
                      prob=prob + rates
                   enddo neutral_types
                   if(sum(prob).le.0.)cycle loop_over_fast_ions
                   !! --------- solve collisional radiative model along track-!!
                   states=prob*denf
                   loop_along_track: do jj=1,ntrack     
                      call colrad(pos_arr(:,jj),vi(:),dt_arr(jj) &
                           ,states,photons,afida_type,nlaunch2,ipos)
                      if(photons.le.0.d0)cycle loop_over_fast_ions
                      if(inputs%calc_spec.eq.1) &
                           call spectrum(vi(:),pos_arr(:,jj),photons,afida_type)
                   enddo loop_along_track
                enddo loop_over_fast_ions
             enddo npa_loop

             !$OMP CRITICAL
             cnt=cnt+1
             if(mod(cnt,100).eq.0)then
                WRITE(*,'(f7.2,"%",a,$)') cnt/maxcnt*100,char(13)
             endif
             !$OMP END CRITICAL
          enddo loop_along_x
       enddo loop_along_y
    enddo loop_along_z
    !$OMP END PARALLEL DO 
  end subroutine fida


  subroutine passive_fida      
    integer                               :: i,j,k   !! indices  x,y,z  of cells
    integer                               :: iion
    real(double), dimension(3)            :: ipos    !! start position
    integer,      dimension(3)            :: iind    !! index of start cell
    real(double), dimension(3)            :: vi      !! velocity of fast ions
    real(double), dimension(:,:),allocatable  :: vi_arr  !! velocity array
    integer,dimension(3)                  :: ac      !! new actual cell
    real(double)                          :: rho
    !! Determination of the CX probability
    real(double), dimension(nlevs)        :: denn    !!  neutral dens (n=1-4)
    real(double)                          :: denf    !! fast-ion density
    real(double), dimension(nlevs)        :: prob    !! Prob. for CX 
    real(double), dimension(3)            :: vnhalo  !! v of halo neutral
    integer                               :: in      !! index of neut rates
    real(double), dimension(nlevs)        :: rates   !! Rate coefficiants for CX
    !! Collisiional radiative model along track
    real(double), dimension(nlevs)        :: states  ! Density of n-states
    integer                               :: ntrack
    real(double), dimension(grid%ntrack)  :: dt_arr  !! time per track step
    real(double), dimension(3,grid%ntrack):: pos_arr !! positions along track
    integer                               :: jj      !! counter along track
    real(double)                          :: photons !! photon flux 
    real(double), dimension(grid%nx,grid%ny,grid%nz)::papprox,nlaunch !! approx. density
    real(double)                          :: nlaunch2   
    real(double)                          :: vi_abs             !! (for NPA)
    real(double), dimension(3)            :: ray,ddet,hit_pos   !! ray towards NPA
    real(double)                          :: papprox_tot 
    integer                               :: inpa    
    real(double)                          :: alpha !! angle relative to detector LOS


    real(double)                          ::denp,ti,vtor,te,dene,deni,bckgrnd_dens

    print*, 'Passive fida simulation'
    

    inputs%nr_fida=inputs%nr_fida/2.
    !! ------------- calculate papprox needed for guess of nlaunch --------!!
    papprox=0.d0
    papprox_tot=0.d0
    do k=1,grid%Nz 
       do j=1,grid%Ny 
          do i=1,grid%Nx
             if (inputs%npa.eq.1) then
                ipos(:)=(/grid%xxc(i),grid%yyc(j),grid%zzc(k)/)
                ray=diag%xyzhead(1,:)-ipos
                alpha=acos(dot_product(npa%los,ray) &
                     /(npa%dlos*sqrt(dot_product(ray,ray))))
                if (alpha.gt.npa%opening_angle*3.) cycle
             endif
             ipos=(/grid%xxc(i),grid%yyc(j),grid%zzc(k)/)
             call denf_interp(ipos,denf)
             call prof_interp(ipos,denp,ti,vtor,te,dene,deni,bckgrnd_dens)
             papprox(i,j,k)= bckgrnd_dens*denf
             call rho_interp(ipos,rho)
             if(rho.lt.1.1) papprox_tot=papprox_tot+papprox(i,j,k) 
          enddo
       enddo
    enddo
    call get_nlaunch(inputs%nr_fida,papprox,papprox_tot,nlaunch)
    print*,'    # of markers: ',int(sum(nlaunch))
    if(inputs%npa.eq.1)print*,'    # npa_loop:   ',int(npa%npa_loop)
    !$OMP PARALLEL DO private(i,j,k,inpa,vi_arr,nlaunch2,iion,ac,ipos,vi, &
    !$OMP& vi_abs,alpha,iind,ray,hit_pos,ddet, &
    !$OMP& pos_arr,dt_arr,ntrack,prob,bckgrnd_dens,denn,rates,denf,vnhalo, &
    !$OMP& states,jj,photons)
    loop_along_z: do k = 1, grid%Nz
       loop_along_y: do j = 1, grid%Ny
          loop_along_x: do i = 1, grid%Nx
             if(nlaunch(i,j,k).le.0)cycle
             ac=(/i,j,k/)
             ipos=(/grid%xxc(i),grid%yyc(j),grid%zzc(k)/)
             !! ------------- loop over the markers ---------------------- !!
             npa_loop: do inpa=1,int(npa%npa_loop)
                if(allocated(vi_arr)) deallocate(vi_arr)
                allocate(vi_arr(int(nlaunch(i,j,k)+grid%ntrack),3))
                call mc_fastion2(ac,ipos,nlaunch(i,j,k),vi_arr(:,:),nlaunch2)
                call denf_interp(ipos,denf)
                loop_over_fast_ions: do iion=1,int(nlaunch2)
                   vi(:)=vi_arr(iion,:)
                   if(sum(vi).eq.0)then
                      print*,'vi is zero!'
                      cycle loop_over_fast_ions
                   endif
                   !! -------- check if particle flies into NPA detector ---- !!
                   if(inputs%npa.eq.1)then  
                      vi_abs=sqrt(dot_product(vi,vi))
                      alpha=acos(dot_product(npa%los,vi(:))/(vi_abs*npa%dlos))
                      if (alpha.gt.npa%opening_angle) cycle loop_over_fast_ions
                   endif
                   iind=(/i,j,k/)
                   call mc_start(iind(:), vi(:), ipos(:))
                   !! -------- check if track ends at the NPA detector ---- !!
                   if(inputs%npa.eq.1)then 
                      !! check if track ends at the NPA detector
                      ray=ipos-diag%xyzhead(1,:)
                      hit_pos(:)=ipos(:)+vi(:)/vi_abs*sqrt(dot_product(ray,ray))
                      ddet=hit_pos-diag%xyzhead(1,:)
                      if(sqrt(dot_product(ddet,ddet))&
                           .gt.npa%size(1)) cycle loop_over_fast_ions
                   endif
                   call track2(vi,ipos, pos_arr, dt_arr, ntrack)
                   if(ntrack.eq.0)cycle loop_over_fast_ions
                   !! ---------------- calculate CX probability --------------!!
                   prob=0.d0
                   denn(:)=0.
                   call prof_interp(ipos,denp,ti,vtor,te,dene,deni,bckgrnd_dens)
                   denn(1)=bckgrnd_dens
                   vnhalo=0.+vtor  ! Lab.frame! -> cold thermal ions move with vrot.
                   call neut_rate(denn,vi,vnhalo,rates)
                   prob=prob + rates

                   if(sum(prob).le.0.)cycle loop_over_fast_ions
                   !! --------- solve collisional radiative model along track-!!
                   states=prob*denf
                   loop_along_track: do jj=1,ntrack     
                      call colrad(pos_arr(:,jj),vi(:),dt_arr(jj) &
                           ,states,photons,pfida_type,nlaunch2,ipos)
                      if(photons.le.0.d0)cycle loop_over_fast_ions
                      if(inputs%calc_spec.eq.1) &
                           call spectrum(vi(:),pos_arr(:,jj),photons,pfida_type)
                   enddo loop_along_track
                enddo loop_over_fast_ions
             enddo npa_loop
          enddo loop_along_x
       enddo loop_along_y
    enddo loop_along_z
    !$OMP END PARALLEL DO 
  end subroutine passive_fida


 
  !*****************************************************************************
  !----------- Calculation of weight functions----------------------------------
  !*****************************************************************************
 subroutine gauss_funct(X, sigma,F)
    real(double), dimension(:),intent(in)  :: X!! x-axis
    real(double),              intent(in)  :: sigma!! width
    real(double), dimension(:),intent(out) :: F!! result
    real(double), dimension(3) :: A!! 
    real(double), dimension(size(X)) :: Z
    A(1)=0.
    A(2)=sigma
    Z=(X-A(1))/A(2)
    F = exp(-0.5d0*Z**2)/A(2)/sqrt(2.*pi)
  end subroutine gauss_funct  
  subroutine weight_function
    real(double), dimension(:), allocatable :: radius !! radial position
    real(double), dimension(:), allocatable :: phi_b_los,phi_a_bxlos !! angles
    real(double), dimension(:), allocatable :: Babs_arr  
    real(double)                   :: photons !! photon flux 
    real(double), dimension(nlevs) :: fdens,hdens,tdens,halodens
    real(double), dimension(3)     :: los_vec,bxl
    real(double), dimension(3)     :: a_norm,b_norm,c_norm
    real(double), dimension(3)     :: brzt,brzt_c
    real(double)                   :: b_abs
    real(double)                   :: ti_c,te_c,dene_c,denp_c,deni_c,vtor_c
    real(double), dimension(3)     :: vrot
    real(double)                   :: ti,te,dene,denp,deni,vtor,length
    real(double)                   :: rad
    integer                        :: nchan,nwav
    real(double),dimension(:)    ,allocatable  :: wav_arr,central_wavel
    integer                        :: ii,i,j,k,l  !! indices wavel,vx,vy,vz
    real(double),dimension(:,:,:,:,:),allocatable :: wfunct
    real(double),dimension(:)    ,allocatable :: ebarr,ptcharr,phiarr
    real(double)                   :: sinus,gyro_angle
    real(double),dimension(3)      :: vi,vi_norm
    real(double)                   :: vabs
    real(double),dimension(3)      :: efield
    real(double),dimension(n_stark):: intens !!intensity vector
    real(double),dimension(n_stark):: wavel  !!wavelength vector [A)
    real(double),dimension(3)      :: vn  ! vi in m/s
   !! Determination of the CX probability
    real(double),dimension(3)      :: vnbi_norm,vnbi_f,vnbi_h,vnbi_t !! Velocity of NBI neutrals 
    real(double),dimension(3)      :: vhalo  !! v of halo neutral
    integer                        :: in      !! index of neut rates
    real(double),dimension(nlevs)  :: rates   !! Rate coefficiants for CX
    real(double),dimension(nlevs)  :: states  ! Density of n-states
    !! COLRAD
    real(double)                   :: dt  !!time interval in cell
    !! ---- Solution of differential equation  ---- ! 
    integer,dimension(3)                  :: ac  !!actual cell
    real(double), dimension(3)            :: pos,ipos !! position of mean cell
    !! B field 
    integer(long):: ir,iz
    real(double):: r,z

    !! profiles 
    real(double):: rho
    integer(long):: irho

    integer                               :: cc 
    real(double),dimension(  grid%ntrack) :: wght   !! radiation wght per cell 
    real(double),dimension(  grid%ntrack) :: los_wght !! los wght 
    real(double),dimension(grid%nx,grid%ny,grid%nz,diag%nchan) :: los_weight !! los wght
    integer(long)                         :: ichan,iichan
    character(120)                        :: filename
    !! length through cloud of neutrals
    real(double), dimension(3,grid%ntrack):: pos_out
    real(double), dimension(3)            :: pos_edge
    integer                               :: jj
    integer                               :: ntrack  !! number of cells
    real(double), dimension(  grid%ntrack):: dt_arr  !! time per cell
    integer,dimension(3,grid%ntrack)      :: ic  !! index of cells
    real(double)                          :: half_max_wght
    !! instrument function
    integer                               :: iif ! counter
    integer, parameter                    :: ninstfu=11 ! elements in the instfument funcion
    real(double), dimension(ninstfu)      :: w_instfu,i_instfu ! wavelength and intenstiy array

    !! DEFINE wavelength array
    nwav=(inputs%wavel_end_wght-inputs%wavel_start_wght)/inputs%dwav_wght
    allocate(wav_arr(nwav+1))
    allocate(central_wavel(nwav))
    print*,nwav,' wavelengths to be simulated!'
    do i=1,nwav+1
       wav_arr(i)=real(i-1.)*inputs%dwav_wght+inputs%wavel_start_wght
    enddo
    do i=1,nwav
       central_wavel(i)=wav_arr(i)+0.5*inputs%dwav_wght
    enddo
    wav_arr=wav_arr*10. !![A]


    !! define pitch, energy and gyro angle arrays
    !! define energy - array
    print*, 'nr of energies, pitches and gyro angles' &
         ,inputs%nr_energy_wght,inputs%nr_pitch_wght,inputs%nr_gyro_wght
    print*, 'maximal energy: ', inputs%emax_wght
    allocate(ebarr(inputs%nr_energy_wght))  
    do i=1,inputs%nr_energy_wght
       ebarr(i)=real(i-0.5)*inputs%emax_wght/real(inputs%nr_energy_wght)
    enddo
    !! define pitch - array
    allocate(ptcharr(inputs%nr_pitch_wght))
    do i=1,inputs%nr_pitch_wght
       ptcharr(i)=real(i-0.5)*2./real(inputs%nr_pitch_wght)-1.
    enddo
    !! define gyro - array
    allocate(phiarr(inputs%nr_gyro_wght))
    do i=1,inputs%nr_gyro_wght
       phiarr(i)=real(i-1.)*2.d0*pi/real(inputs%nr_gyro_wght)
    enddo
    
    if(inputs%ichan_wght.gt.0) then
       nchan=1
    else
       nchan=diag%nchan
    endif
    
    allocate(radius(nchan),phi_b_los(nchan),phi_a_bxlos(nchan))
    allocate(babs_arr(nchan))
    phi_b_los=0.d0 ; phi_a_bxlos=0.d0 ; radius=0.d0
    !! define storage arrays
    if(inputs%wght_gyro_resolved.eq.1)then
       print*,'full weight functions with gyro motion'
       allocate(wfunct(nchan,nwav,inputs%nr_energy_wght &
            ,inputs%nr_pitch_wght,inputs%nr_gyro_wght))
    else
       print*,'average weight functions over gyro motion'
       allocate(wfunct(nchan,nwav,inputs%nr_energy_wght,inputs%nr_pitch_wght,1))
    endif
    wfunct       = 0.d0
    !!save the los-weights into an array
    !! because the structure is over-written
    los_weight(:,:,:,:)=diag%los_wght(:,:,:,:)
    loop_over_channels: do iichan=1,nchan
       
       if(inputs%ichan_wght.le.0) then
          ichan=iichan
       else 
          ichan=inputs%ichan_wght
       endif
       
       !! Define instfumental function
       do i=1,ninstfu
          w_instfu(i)=real(i-1.)/dble(ninstfu-1.)*6.d0*diag%instfu(ichan)-3.d0*diag%instfu(ichan)
       enddo
       call gauss_funct(w_instfu,diag%instfu(ichan),i_instfu)
       i_instfu=i_instfu/sum(i_instfu)
  



       print*,'channel:',ichan
       print*, diag%los_name(ichan)
       radius(iichan)=sqrt(diag%xyzlos(ichan,1)**2 &
            +diag%xyzlos(ichan,2)**2)
       print*,'Radius:',radius(iichan)
       !! Calcullate mean kinetic profiles...
       cc=0       ; los_wght=0.d0 ; wght=0.d0
       fdens=0.d0 ; hdens=0.d0    ; tdens=0.d0  ; halodens=0.d0  
       brzt=0.d0
       ti=0.d0    ; te=0.d0
       dene=0.d0  ; denp=0.d0     ; deni=0.d0
       vtor=0.d0  ; pos=0.d0 
       call read_field
       deallocate(prof%rho, prof%te, prof%ti,prof%dene &
            , prof%denp, prof%deni,prof%vtor, prof%zeff,prof%bckgrnd_dens)
       call read_profiles
       do k=1,grid%nz
          do j=1,grid%ny 
             do i=1,grid%nx 
                if(los_weight(i,j,k,ichan).gt.0.)then
                   cc=cc+1
                   if(cc.gt.grid%ntrack) stop 'ntrack value too low!!'
                   los_wght(cc)=los_weight(i,j,k,ichan)
                   !! integrate beam/halo density
                   fdens=fdens &
                        +result%neut_dens(i,j,k,:,nbif_type)*los_wght(cc) 
                   hdens=hdens &
                        +result%neut_dens(i,j,k,:,nbih_type)*los_wght(cc) 
                   tdens=tdens &
                        +result%neut_dens(i,j,k,:,nbit_type)*los_wght(cc)
                   halodens=halodens &
                        +result%neut_dens(i,j,k,:,halo_type)*los_wght(cc) 
                   !! determine mean values along LOS
                   wght(cc)=(result%neut_dens(i,j,k,3,nbif_type)   &
                        + result%neut_dens(i,j,k,3,nbih_type) &
                        + result%neut_dens(i,j,k,3,nbit_type)  &
                        + result%neut_dens(i,j,k,3,halo_type))*los_wght(cc)
                   ipos=(/grid%xxc(i),grid%yyc(j),grid%zzc(k)/)
                   call bfield_interp(ipos,b_abs,b_norm,a_norm,c_norm,brzt_c)
                   brzt    =brzt  + brzt_c * wght(cc)
                   call prof_interp(ipos,denp_c,ti_c,vtor_c &
                        ,te_c,dene_c,deni_c)
                   ti     =ti     +ti_c   * wght(cc)
                   te     =te     +te_c   * wght(cc)
                   dene   =dene   +dene_c * wght(cc)
                   denp   =denp   +denp_c * wght(cc)
                   deni   =deni   +deni_c * wght(cc)
                   vtor   =vtor   +vtor_c * wght(cc)/sqrt(ipos(1)**2+ipos(2)**2)
                   pos(:)=pos(:)+ipos*wght(cc)
                endif
             enddo
          enddo
       enddo

       !! there is no Impurity ion halo calculation. Therefore,  
       !! assume that there were no impurities
       if(denp.gt.0)halodens=halodens*dene/denp


       rad= sum(wght)
       pos= pos / rad
       !! Determine NBI vector
       vnbi_norm=pos(:)-nbi%xyz_pos(:)
       vnbi_norm=vnbi_norm/sqrt(dot_product(vnbi_norm,vnbi_norm))
       vnbi_f=vnbi_norm*nbi%vinj
       vnbi_h=vnbi_norm*nbi%vinj/sqrt(2.d0)
       vnbi_t=vnbi_norm*nbi%vinj/sqrt(3.d0)     
       !! store kinetic quantities in irho and irho+1
       call rho_interp(pos,rho)
       irho=floor(rho/prof%drho)+1  
       prof%denp(irho)  =denp/ rad
       prof%denp(irho+1)=denp/ rad
       prof%dene(irho)=dene/ rad
       prof%dene(irho+1)=dene/ rad
       prof%deni(irho)=deni/ rad
       prof%deni(irho+1)=deni/ rad
       prof%ti(irho)=ti/ rad
       prof%ti(irho+1)=ti/ rad
       prof%te(irho)=te/ rad
       prof%te(irho+1)=te/ rad
       prof%vtor(irho)=vtor/ rad
       prof%vtor(irho+1)=vtor/ rad
       !! store new bfield in structure
       r=sqrt(pos(1)**2+pos(2)**2)
       z=pos(3)
       ir=floor((r-grid%rmin)/grid%dr) +1
       iz=floor((z-grid%zmin)/grid%dz) +1
       grid%Brzt(ir  ,iz  ,:)=brzt/rad
       grid%Brzt(ir+1,iz  ,:)=brzt/rad
       grid%Brzt(ir  ,iz+1,:)=brzt/rad
       grid%Brzt(ir+1,iz+1,:)=brzt/rad
       call bfield_interp(pos,b_abs,b_norm,a_norm,c_norm)
       print*, '|B|: ',real(b_abs,float), ' T'
       print*,'pitch of NBI:',dot_product(vnbi_norm,b_norm)*inputs%btipsign
       babs_arr(iichan)=b_abs
       !! set los_wght to 1 only for one channel (needed for spectrum routine)
       call get_ac(pos,ac)
       diag%los_wght(ac(1),ac(2),ac(3),:)=0.
       diag%los_wght(ac(1),ac(2),ac(3),ichan)=1.
       !! Determine the angle between the B-field and the Line of Sight
       los_vec(1)=pos(1)-diag%xyzhead(ichan,1) 
       los_vec(2)=pos(2)-diag%xyzhead(ichan,2) 
       los_vec(3)=pos(3)-diag%xyzhead(ichan,3) 
       los_vec=los_vec/sqrt(dot_product(los_vec,los_vec)) 
            
       !!angle between a_norm and los_vec
       bxl(1) =  b_norm(2)*los_vec(3) - b_norm(3)*los_vec(2)
       bxl(2) =-(b_norm(1)*los_vec(3) - b_norm(3)*los_vec(1))
       bxl(3) =  b_norm(1)*los_vec(2) - b_norm(2)*los_vec(1)
       bxl=bxl/sqrt(dot_product(bxl,bxl)) 
       phi_a_bxlos(iichan)=acos(dot_product(a_norm,bxl))
       print*,'angle between a_norm and B x Los:', phi_a_bxlos(iichan)*180./pi


       phi_b_los(iichan)=180.-acos(dot_product(b_norm,los_vec))*180./pi
       print*,'Angle between B and LOS [deg]:', phi_b_los(iichan)
       !! START calculation of weight functions

       !! do the main simulation  !! 
       !$OMP PARALLEL DO private(i,vabs,j,sinus,k,gyro_angle,vi_norm, &
       !$OMP& dt_arr,ic,pos_out,ntrack,pos_edge, &
       !$OMP& wght,jj,length,dt,half_max_wght,  &
       !$OMP& states,vi,rates,in,vhalo,photons,wavel,intens,l,iif,ii)
 
       
       !! LOOP over the three velocity vector components 
       !! (pitch,gyro angle, energy)
       do j = 1, inputs%nr_pitch_wght !! pitch loop
          sinus = sqrt(1.d0-ptcharr(j)**2)
          do k = 1, inputs%nr_gyro_wght !! gyro angle
             !! cacluate velocity vector from energy,pitch, gyro angle
             gyro_angle=phiarr(k)-phi_a_bxlos(iichan)-.5*pi
             vi_norm=sinus*cos(gyro_angle)*a_norm &
                  +    ptcharr(j)*b_norm*inputs%btipsign &
                  +    sinus*sin(gyro_angle)*c_norm
             !! calcualte possible trajectory of fast-neutral that intersect 
             !! the measurment position with this velocity 
             !! The measurement postion is at maximum NBI density along LOS
             call track(-vi_norm,pos,pos_out,dt_arr,ntrack,ic)
             if(ntrack.gt.2)then
                pos_edge=pos_out(:,ntrack-2)
             else
                pos_edge=pos_out(:,ntrack)
             endif
             !! now determine the FWHM track length throught the grid!
             call track(vi_norm,pos_edge,pos_out,dt_arr,ntrack,ic)
             !! Calculate averge beam density seen by the fast-ion        
             wght=0.d0
             if(ntrack.gt.grid%ntrack) stop 'ntrack value too low!!'
             do jj=1,ntrack
                call rho_interp(pos_out(:,jj),rho)
                if(rho.lt.1)then
                   wght(jj)=sum( &
                        result%neut_dens(ic(1,jj),ic(2,jj),ic(3,jj),:,nbif_type)+ &
                        result%neut_dens(ic(1,jj),ic(2,jj),ic(3,jj),:,nbih_type)+ &
                        result%neut_dens(ic(1,jj),ic(2,jj),ic(3,jj),:,nbit_type)+ &
                        result%neut_dens(ic(1,jj),ic(2,jj),ic(3,jj),:,halo_type))
                endif
             enddo
             half_max_wght=0.5*maxval(wght)
             length=0.d0
             do jj=1,ntrack
                if (wght(jj).gt.half_max_wght)length=length+dt_arr(jj)
             enddo
             !! length=sum(wght(:)*dt_arr(:))/maxval(wght) ! (FWHM)
             !! determine time by length and velocity 
             !! calculate the average time until a fast-neutral is 
             !! detected after its neutralization
             do i = 1, inputs%nr_energy_wght !! energy loop
                vabs = sqrt(ebarr(i)/(v_to_E*inputs%ab))
                dt=length/vabs
                !! -------------- calculate CX probability -------!!
                ! CX with full energetic NBI neutrals
                states=0.d0
                vi(:) = vi_norm(:)*vabs
                call neut_rate(fdens,vi,vnbi_f,rates)
                states=states + rates
                ! CX with half energetic NBI neutrals
                call neut_rate(hdens,vi,vnbi_h,rates)
                states=states + rates
                ! CX with third energetic NBI neutrals
                call neut_rate(tdens,vi,vnbi_t,rates)
                states=states + rates
                ! CX with HALO neutrals
                call neut_rate_beam_therm(halodens,vi,ipos,rates)
                states=states + rates
                call colrad(pos,vi,dt,states,photons,pfida_type,1.d0)
                !! photons: [Ph*cm/s/fast-ion]-!!
                !! calcualte spectrum of this one fast-ion
                call spectrum(vi,pos,1.d0,pfida_type,wavel,intens)
                stark_components: do l=1,n_stark 
                   instrument_function:do iif=1,ninstfu
                      wavelength_ranges: do ii=1,nwav
                         if (wavel(l)+w_instfu(iif).ge.wav_arr(ii).and. &
                              wavel(l)+w_instfu(iif).lt.wav_arr(ii+1)) then
                            !$OMP CRITICAL(wf)
                            if(inputs%wght_gyro_resolved.eq.1)then
                               wfunct(iichan,ii,i,j,k) = wfunct(iichan,ii,i,j,k) &
                                    + intens(l)*photons*i_instfu(iif)
                            else
                               wfunct(iichan,ii,i,j,1) = wfunct(iichan,ii,i,j,1) &
                                    + intens(l)*photons*i_instfu(iif) &
                                    /real(inputs%nr_gyro_wght)
                            endif
                            !$OMP END CRITICAL(wf)
                         endif
                      enddo wavelength_ranges
                   enddo instrument_function
                enddo stark_components
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
       !!wfunct:[Ph*cm/s] !!
    enddo loop_over_channels

    print*,'end of weight function calcultaion'

 !! Open file for the outputs
    filename=trim(adjustl(result_dir))//"/weight_function.bin" 
    open (66, file =filename,access='stream')
    write(66)inputs%shot_number
    write(66)inputs%time
    write(66)inputs%ichan_wght 
    write(66)inputs%nr_energy_wght  !!Nr of energies
    write(66)ebarr
    write(66)ebarr(2)-ebarr(1)   
    write(66)inputs%nr_pitch_wght  !!Nr. of pitches
    write(66)ptcharr
    write(66)abs(ptcharr(2)-ptcharr(1))
    write(66)inputs%nr_gyro_wght  !!Nr. of gyro angles
    write(66)phiarr
    write(66)abs(phiarr(2)-phiarr(1))       
    write(66)nwav
    write(66)central_wavel
    write(66)inputs%dwav_wght
    write(66)inputs%wavel_start_wght
    write(66)inputs%wavel_end_wght
    write(66)nchan
    write(66)phi_b_los
    write(66)babs_arr
    write(66)radius
    write(66)real(wfunct(:,:,:,:,:),float)
    close (66)

    print*, 'weight function written to: ',filename
    deallocate(ebarr)  
    deallocate(ptcharr)
    deallocate(phiarr)
    deallocate(wav_arr)
    deallocate(central_wavel)
    deallocate(wfunct)
    deallocate(radius)
    deallocate(phi_b_los)
    deallocate(phi_a_bxlos)
    deallocate(babs_arr)
  end subroutine weight_function

  subroutine test_orbit
    real(double), dimension(3) :: ipos,iv
    real(double) :: mass, charge,dt
    real(double), dimension(:,:), allocatable :: pos_arr
    integer(long) :: npos

    npos=100
    allocate(pos_arr(3,npos))

    ipos=(/200.,0.,0./)
    mass=2.d0*mass_u
    charge=1.d0*e0
    print*, 'charge,mass:', charge,mass
    iv=  (/-1. , 0.3 , 0./) &
         *sqrt(2.d0*160.0*1.d3*e0/mass)*100.
    print*,'vi:' ,iv
    print*, 0.5*mass*dot_product(iv,iv)/100.**2/e0/1.e3
    call orbit(ipos,iv,mass,charge,npos,pos_arr,dt)
    !! write to file
    open (66, file ='trajectory.bin',access='stream')
    write(66)npos
    write(66)real(pos_arr,float)
    close(66)
    deallocate(pos_arr)
    stop 'end of orbit'

  end subroutine test_orbit
end module application
!*****************************************************************************
!-----------Main proagramm ---------------------------------------------------
!*****************************************************************************
program fidasim
  use application
  implicit none 
  integer, dimension(8)              :: time_arr,time_start,time_end !Time array
  integer                            :: i,j,k,n,los,seed
  integer                            :: hour,minu,sec
  real(double)                       :: random_init

 
  !! measure time
  call date_and_time (values=time_start)
  write(*,"(A,I2,A,I2.2,A,I2.2)") 'START: hour, minute, second: '  &
       ,time_start(5), ':', time_start(6), ':',time_start(7)
  !! get filename of input
  call getarg(1,result_dir)
  !result_dir='RESULTS/28746A10_fi_3'
  !! ----------------------------------------------------------
  !! ------ INITIALIZE THE RANDOM NUMBER GENERATOR  -----------
  !! ----------------------------------------------------------
  seed = 4 !! has to be negative
  ran_am=nearest(1.0,-1.0)/ran_IM
  ran_iy=ior(ieor(888889999,seed),1)
  ran_ix=ieor(777755555,seed)
  random_init=ran()

  !! ----------------------------------------------------------
  !! ------- READ INPUTS, PROFILES, LOS AND TABLES  -----------
  !! ----------------------------------------------------------
  call read_namelist
  call read_grid
  call read_nbi
  call read_field
  call read_profiles
  call read_diag
  call read_tables
  call read_fbm


  !! ----------------------------------------------------------
  !! --------------- TEST ORBIT SIMULATION ---------------
  !! ----------------------------------------------------------
  !call test_orbit
  !! ----------------------------------------------------------
  !! --------------- ALLOCATE THE RESULT ARRAYS ---------------
  !! ----------------------------------------------------------
  !! neutral density array!
  allocate(result%neut_dens(grid%Nx,grid%Ny,grid%Nz,nlevs,ntypes))
  result%neut_dens(:,:,:,:,:)=0.d0
  !! velocity vectors!
  !! nx,ny,nz,ntypes,ncounters,3 velocity dimensions
  allocate(result%velocity_vectors(grid%Nx,grid%Ny,grid%Nz,ntypes,nvelocities,3))
  result%velocity_vectors=0.d0
  allocate(result%velocity_counter(grid%Nx,grid%Ny,grid%Nz,ntypes))
  allocate(result%vel_vec_red_lev(grid%Nx,grid%Ny,grid%Nz,ntypes))
  result%vel_vec_red_lev=1.
  result%velocity_counter=0.
  do i=1,nvelocities/2
     every_second_index(i)=i*2
     first_half_indices(i)=i
  enddo


  !! birth profile
  if(inputs%calc_birth.eq.1)then
     allocate(result%birth_dens(grid%Nx,grid%Ny,grid%Nz,3,npitch_birth))
     result%birth_dens(:,:,:,:,:)=0.d0
  endif
  !! allocate the spectra array
  if(inputs%calc_spec.eq.1)then
     allocate(result%spectra(2,diag%nlambda,diag%nchan,ntypes))
     result%spectra(:,:,:,:)=0.d0
  endif
  if(inputs%npa.eq.1)then
     print*, 'this is a NPA simultation!'
     inputs%nr_npa=1000000
     allocate(npa%v(inputs%nr_npa,3)) 
     allocate(npa%ipos(inputs%nr_npa,3)) 
     allocate(npa%fpos(inputs%nr_npa,3)) 
     allocate(npa%wght(inputs%nr_npa))
     allocate(npa%kind(inputs%nr_npa))
     npa%counter=0    
  endif
  

  !! -----------------------------------------------------------------------
  !! --------------- CALCULATE/LOAD the BEAM and HALO DENSITY---------------
  !! -----------------------------------------------------------------------
  if(inputs%load_neutrals.eq.1)then
     call read_neutrals()
  else
     if(inputs%nr_ndmc.gt.0)then
        !! ----------- ndmc (neutral density monte carlo ------------------- !! 
        call date_and_time (values=time_arr)
        write(*,"(A,I2,A,I2.2,A,I2.2)") 'ndmc:   ' ,time_arr(5), ':' &
             , time_arr(6), ':',time_arr(7)
        call ndmc
        if(inputs%calc_birth.eq.1)then
           call write_birth_profile()
        endif
        !! do the HALO calcualtion only if enough markers are defined!
        if(inputs%nr_halo.gt.10)then
           !! -------------------------- DCX (Direct charge exchange) ------ !!
           call date_and_time (values=time_arr)
           write(*,"(A,I2,A,I2.2,A,I2.2)") 'dcx:    ' ,time_arr(5), ':' &
                , time_arr(6), ':',time_arr(7)
           call dcx
           !! ------------------------- HALO ------------------------------- !!
           call date_and_time (values=time_arr)
           write(*,"(A,I2,A,I2.2,A,I2.2)") 'halo:   ' ,time_arr(5), ':' &
                , time_arr(6), ':',time_arr(7)
           call halo
        endif
        !! ------------------ write output ---------------------------------- !!
        call write_neutrals()
        if(inputs%calc_spec.eq.1)then
           call write_nbi_halo_spectra()
        endif
     else
        call write_neutrals()
     endif
  endif


  !! -----------------------------------------------------------------------
  !! --------------- CALCULATE the FIDA RADIATION/ NPA FLUX ----------------
  !! -----------------------------------------------------------------------
  if(inputs%nr_fida.gt.0.and.inputs%simfa.ne.0)then    
     if(inputs%nr_ndmc.gt.0)then
        call date_and_time (values=time_arr)
        write(*,"(A,I2,A,I2.2,A,I2.2)") 'D-alpha main: ' ,time_arr(5), ':' &
             , time_arr(6), ':',time_arr(7)
        print*,'start fida'
        call fida
     endif
     if(sum(prof%bckgrnd_dens).gt.0)call passive_fida
     !! ------- Store Spectra and neutral densities in binary files ------ !!
     if(inputs%calc_spec.eq.1)then
        call write_fida_spectra()
     endif
     !! ---------------- Store NPA simulation ------------ !!
     if(inputs%npa.eq.1)then
        call write_npa()
     endif
  endif

  !! -------------------------------------------------------------------
  !! ----------- Calculation of weight functions -----------------------
  !! -------------------------------------------------------------------
  if(inputs%calc_wght.eq.1) then 
     colrad_threshold=0. !! to speed up simulation!
     call date_and_time (values=time_arr)
     write(*,"(A,I2,A,I2.2,A,I2.2)") 'weight function:    '  &
          ,time_arr(5), ':', time_arr(6), ':',time_arr(7)
     call weight_function()
  endif


  call date_and_time (values=time_arr)
  write(*,"(A,I2,A,I2.2,A,I2.2)") 'END: hour, minute, second: '  &
       ,time_arr(5), ':', time_arr(6), ':',time_arr(7)

  call date_and_time (values=time_end)
  hour = time_end(5) - time_start(5)
  minu = time_end(6) - time_start(6)
  sec  = time_end(7) - time_start(7)
  if (minu.lt.0.) then
    minu = minu +60
    hour = hour -1
  endif
  if (sec.lt.0.) then
    sec  = sec +60
    minu = minu -1
  endif
    
  write(*,"(A,I2,A,I2.2,A,I2.2)") 'duration:                  ' &
       ,hour, ':',minu, ':',sec

  !! -------------------------------------------------------------------
  !! --------------- Finally, deallocate allocated arrays --------------
  !! -------------------------------------------------------------------
  !! tables structure
  deallocate(tables%qp)
  deallocate(tables%neut)
  deallocate(tables%qi)
  deallocate(tables%qe)
  deallocate(tables%einstein) 
  !! distribution structure
  if(allocated(fbm%energy))then
     deallocate(fbm%energy)
     deallocate(fbm%pitch)
     deallocate(fbm%rgrid)
     deallocate(fbm%zgrid)
     deallocate(fbm%fbm)
     deallocate(fbm%denf)
  endif
  !! grid and cell structure
  deallocate(diag%los_wght)
  deallocate(grid%xx)
  deallocate(grid%yy)
  deallocate(grid%zz)
  deallocate(grid%rr)
  deallocate(grid%xxc)
  deallocate(grid%yyc)
  deallocate(grid%zzc)
  deallocate(grid%rrc)
  deallocate(grid%Brzt) 
  deallocate(grid%E) 
  deallocate(grid%rho) 
  !! Diagnostic inputs
  deallocate(prof%rho, prof%te, prof%ti,prof%dene &
       , prof%denp, prof%deni,prof%vtor, prof%zeff)
  deallocate(diag%xyzlos, diag%xyzhead, diag%headsize &
       , diag%opening_angle, diag%sigma_pi,diag%instfu &
       , diag%los_name)
  !! result arrays
  deallocate(result%neut_dens)
  deallocate(result%velocity_vectors)
  deallocate(result%velocity_counter)
  deallocate(result%vel_vec_red_lev)
  if(inputs%calc_spec.eq.1)deallocate(result%spectra)
  if(inputs%calc_birth.eq.1)deallocate(result%birth_dens) 

  !!deallocate npa arrays
  if(inputs%npa .eq. 1)then
     deallocate(npa%v) 
     deallocate(npa%ipos) 
     deallocate(npa%fpos) 
     deallocate(npa%wght)
     deallocate(npa%kind)
  endif
end program fidasim
 
