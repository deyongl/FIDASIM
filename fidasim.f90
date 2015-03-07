!!FIDASIM Version 0.4.0

!The main routine (fidasim) is at the end of the file!
module application
  use netcdf
  use fidasim_linalg
  implicit none
  !!                      Definition for the kind of the variables:
  integer , parameter   :: long      = kind(int(1))
  integer , parameter   :: float     = kind(1.e0)
  integer , parameter   :: double    = kind(1.d0)
  !!                      Indizes of the different components:
  character(120)        :: namelist_file
  character(120)        :: root_dir
  integer , parameter   :: nbif_type = 1 ! full energy NBI spectra/density
  integer , parameter   :: nbih_type = 2 ! half energy NBI spectra/density
  integer , parameter   :: nbit_type = 3 ! third energy NBI spectra/density
  integer , parameter   :: halo_type = 4 ! halo spectra/density
  integer , parameter   :: fida_type = 5 ! fida spectra/density
  integer , parameter   :: brems_type= 6 ! brems-strahlung
  integer , parameter   :: ntypes    = 6 ! number of different types of neutrals
  !! random number geneator save variables
  real(double)           :: ran_am
  integer                :: ran_ix=-1,ran_iy=-1,ran_k
  integer, parameter     :: ran_IA=16807,ran_IM=2147483647
  integer, parameter     :: ran_IQ=127773,ran_IR=2836
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

  !!Numerical Settings
  integer,parameter:: nlevs=6             !!nr of quantum states
  integer,parameter:: npitch_birth=100    !!nr pitches in birth profile
  real(double),parameter :: n_halo_neutrate=20. !! to average halo neut-rate
  real(double) :: colrad_threshold=1.d6 !! to speed up simulation!
  !! 1.d6 is a very low density rate per marker (typially ~1.d12)!
  !! it should not affect the calculation
  integer :: nbi_outside=0

  type beam_grid_type
    real(double), dimension(3) :: origin !! origin
    real(double), dimension(3) :: center !! Center of grid
    real(double)               :: alpha  !! rotation about z
    real(double)               :: beta   !! rotation about y/tilt
    real(double), dimension(3) :: dr     !! dx, dy, dz
    real(double), dimension(3) :: lwh    !! Grid length(x),width(y),height(z)
    real(double)               :: drmin  !! min(dx,dy,dz)
    real(double)               :: dv     !! volume of cells
    real(double)               :: volume !! Grid volume
    real(double)
    integer(long)              :: nx     !! Nr. of cells in x direction
    integer(long)              :: ny     !! Nr. of cells in y direction
    integer(long)              :: nz     !! Nr. of cells in z direction
    integer(long)              :: ntrack !! Maximum Nr. of cells for tracking
    integer(long)              :: ngrid  !! Nr. of cells
    real(double), dimension(:), allocatable :: xc !! X centers
    real(double), dimension(:), allocatable :: yc !! Y centers
    real(double), dimension(:), allocatable :: zc !! Z centers
  end type beam_grid_type

  type inter_grid_type
    integer(long) :: nr !! Number of radii
    integer(long) :: nw !! Number of W
    real(double)  :: dr !! Radial spacing
    real(double)  :: dw !! Vertical spacing
    real(double)  :: da !! Area
    real(double), dimension(:),   allocatable :: rr  !! Radius
    real(double), dimension(:),   allocatable :: ww  !! W/Z
    real(double), dimension(:,:), allocatable :: r2d !! 2D radius grid
    real(double), dimension(:,:), allocatable :: w2d !! 2D W grid
  end type inter_grid_type

  type profiles_type
    integer(long) :: nrho
    real(double)  :: drho
    real(double)  :: rhomax
    !! kinetic profiles
    real(double), dimension(:),  allocatable :: rho  !! Magnetic flux coordinate
    real(double), dimension(:),  allocatable :: te   !! Electron temperature
    real(double), dimension(:),  allocatable :: ti   !! Ion temperature
    real(double), dimension(:),  allocatable :: dene !! electron density
    real(double), dimension(:),  allocatable :: denp !! D-ions density
    real(double), dimension(:),  allocatable :: deni !! Carbon densitych
    real(double), dimension(:),  allocatable :: zeff !! zeff
    real(double), dimension(:),  allocatable :: omega!! Torodial rotation
  end type profiles_type

  type fbm_type
    integer(long) :: nenergy !! Number of energies
    integer(long) :: npitch  !! Number of pitches
    real(double)  :: deb     !! Energy spacing
    real(double)  :: dpitch  !! Pitch spacing
    real(double), dimension(:),       allocatable :: energy  !! Energy array
    real(double), dimension(:),       allocatable :: pitch   !! Pitch array
    real(double), dimension(:,:),     allocatable :: denf    !! Fast ion density (r,w)
    real(double), dimension(:,:),     allocatable :: fbm_norm!! Fast ion density (r,w)
    real(double), dimension(:,:,:,:), allocatable :: fbm     !! FBM (E,p,r,w)
  end type fbm_type

  type equil_type
    real(double), dimension(:,:), allocatable :: br !! Radial magnetic field component (r,w)
    real(double), dimension(:,:), allocatable :: bt !! Torodial magnetic field component (r,w)
    real(double), dimension(:,:), allocatable :: bw !! Magnetic field in W/Z direction (r,w)
    real(double), dimension(:,:), allocatable :: er !! Radial electric field component (r,w)
    real(double), dimension(:,:), allocatable :: et !! Torodial electric field component (r,w)
    real(double), dimension(:,:), allocatable :: ew !! Electrid field in W/Z direction (r,w)
    real(double), dimension(:,:), allocatable :: rho2d !! Magnetic flux coordinates (r,w)
  end type equil_type

  type nbi_type
    integer(long)  :: number                !! number of the NBI
    real(double)   :: dv                    !! half width in y direction
    real(double)   :: dw                    !! half width in z direction
    real(double)   :: focy                  !! focal lenght in y direction
    real(double)   :: focz                  !! focal lenght in z direction
    real(double), dimension(3)   :: divay   !! divergence in y direction
    real(double), dimension(3)   :: divaz   !! divergence in z direction
    real(double), dimension(3)   :: species_mix
    real(double), dimension(3)   :: xyz_src !! position of source
    real(double), dimension(3)   :: xyz_pos !! position along beam sightline
    real(double)                 :: einj    !! NBI voltage  [kV]
    real(double)                 :: pinj    !! NBI power    [MW]
    real(double)                 :: vinj    !! NBI velocity [cm/s]
    real(double)                 :: alpha   !! Z rotation not same as beam_grid%alpha
    real(double)                 :: beta    !! Tilt rotation not same as beam_grid%beta
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
    real(double), dimension(:,:,:,:), allocatable :: qp
    real(double), dimension(:,:,:,:), allocatable :: qi
    real(double), dimension(:,:,:,:), allocatable :: qe
    real(double), dimension(:,:,:,:), allocatable :: cx
    real(double), dimension(:,:,:),   allocatable :: neut
    real(double), dimension(:,:),     allocatable :: einstein
  end type tables_type

  type chords_type
    real(double), dimension(:,:),     allocatable :: xyzlos
    real(double), dimension(:,:),     allocatable :: xyzlens
    real(double), dimension(:),       allocatable :: ra
    real(double), dimension(:),       allocatable :: rd
    real(double), dimension(:),       allocatable :: h
    real(double), dimension(:),       allocatable :: chan_id
    real(double), dimension(:),       allocatable :: sigma_pi
    real(double), dimension(:),       allocatable :: radius
    real(double), dimension(:,:,:,:), allocatable :: los_wght
    integer(long) :: nchan
    integer(long) :: nlambda
    real(double)  :: dlambda
    real(double)  :: lambdamin
    real(double)  :: lambdamax
  end type chords_type

  type npa_type
    real(double), dimension(:,:,:) ,allocatable  :: v    !! velocity array
    real(double), dimension(:,:,:) ,allocatable  :: ipos !! initial position array
    real(double), dimension(:,:,:) ,allocatable  :: fpos !! final position array
    real(double), dimension(:,:)   ,allocatable  :: wght !! weight of ith mc particle
    real(double), dimension(:,:)   ,allocatable  :: E    !! Energy of ith mc particle
    real(double), dimension(:,:)   ,allocatable  :: flux !! flux
    real(double), dimension(:)     ,allocatable  :: energy !! energy array
    integer(long),dimension(:)     ,allocatable  :: counter
    real(double)                                 :: npa_loop
    integer(long)                                :: nchan
  end type npa_type

  type result_type
    real(double), dimension(:,:,:,:,:), allocatable :: neut_dens! Density
    real(double), dimension(:,:,:)    , allocatable :: spectra
    real(double), dimension(:,:,:,:,:), allocatable :: birth_dens!Deposition prof
  end type result_type

  type inputs_type
    integer(long) :: shot_number
    real(double)  :: time
    character(120):: runid
    character(4)  :: diag
    character(120):: result_dir
    !! Monte Carlo Settings
    integer(long) :: n_fast
    integer(long) :: n_nbi
    integer(long) :: n_dcx
    integer(long) :: n_halo
    integer(long) :: n_npa
    !! general settings
    integer(long) :: calc_spec
    integer(long) :: load_neutrals
    integer(long) :: load_fbm
    integer(long) :: calc_brems       !! 0 to use IDL v.b.
    integer(long) :: calc_npa
    integer(long) :: calc_fida_wght
    integer(long) :: calc_npa_wght
    integer(long) :: calc_birth
    integer(long) :: interactive
    !! Plasma parameters
    integer(long) :: impurity_charge
    real(double)  :: btipsign
    real(double)  :: ai   !! atomic mass of plasma ions
    real(double)  :: ab   !! atomic mass of beam neutrals
    real(double)  :: rho_max   !! Maximum normalized flux coord
    !! Settings for weight function calculation
    integer(long) :: ne_wght
    integer(long) :: np_wght
    integer(long) :: nphi_wght
    integer(long) :: ichan_wght
    real(double)  :: emax_wght
    real(double)  :: dwav_wght
    real(double)  :: wavel_start_wght
    real(double)  :: wavel_end_wght
  end type inputs_type

  !! definition of the structures:
  type(beam_grid_type)  :: beam_grid
  type(inter_grid_type) :: inter_grid
  type(profiles_type)    :: profiles
  type(fbm_type)        :: fbm
  type(equil_type)      :: equil
  type(nbi_type)        :: nbi
  type(tables_type)     :: tables
  type(npa_type)        :: npa
  type(chords_type)     :: chords
  type(inputs_type)     :: inputs
  type(result_type)     :: result

contains

  !****************************************************************************
  !*******************************I/O Routines*********************************
  !****************************************************************************
  subroutine check(stat)
    use netcdf
    integer, intent ( in) :: stat

    if(stat /= nf90_noerr) then
      print *, trim(nf90_strerror(stat))
      stop "Stopped: Failed to Write/Read netCDF file"
    end if
  end subroutine check

  subroutine read_inputs
    character(120) :: runid,result_dir
    integer       :: calc_spec,calc_npa,calc_birth,calc_fida_wght
    integer       :: calc_npa_wght,load_neutrals,load_fbm,interactive
    integer(long) :: shot,n_fast,n_nbi,n_halo,nlambda,ne_wght,np_wght,nphi_wght,ichan_wght
    real(double)  :: time,lambdamin,lambdamax,emax_wght,dwav_wght,wavel_start_wght,wavel_end_wght

    NAMELIST /fidasim_inputs/ runid,result_dir,calc_spec,calc_npa, &
        calc_birth,calc_fida_wght,calc_npa_wght, load_neutrals,load_fbm,shot, &
        n_fast,n_nbi,n_halo,nlambda,ne_wght,np_wght,nphi_wght,ichan_wght, &
        time,lambdamin,lambdamax,emax_wght,dwav_wght,wavel_start_wght,wavel_end_wght,interactive

    call getenv("FIDASIM_DIR",root_dir)
    print*,'---- loading inputs ----'
    print*,namelist_file
    open(13,file=namelist_file)
    read(13,NML=fidasim_inputs)
    close(13)

    !!General Information
    inputs%shot_number=shot
    inputs%time=time
    inputs%runid=runid
    inputs%result_dir=result_dir

    !!Simulation Switches
    inputs%calc_spec=calc_spec
    inputs%calc_npa=calc_npa
    inputs%calc_birth=calc_birth
    inputs%calc_fida_wght=calc_fida_wght
    inputs%calc_npa_wght=calc_npa_wght
    inputs%load_neutrals=load_neutrals
    inputs%load_fbm=load_fbm
    inputs%interactive=interactive

    !!Monte Carlo Settings
    inputs%n_fast=n_fast
    inputs%n_nbi=n_nbi
    inputs%n_halo=n_halo
    inputs%n_dcx=n_halo

    !!Weight Function Settings
    inputs%ne_wght=ne_wght
    inputs%np_wght=np_wght
    inputs%nphi_wght=nphi_wght
    inputs%ichan_wght=ichan_wght
    inputs%emax_wght=emax_wght
    inputs%dwav_wght=dwav_wght
    inputs%wavel_start_wght=wavel_start_wght
    inputs%wavel_end_wght=wavel_end_wght

    !!Wavelength Grid Settings
    chords%nlambda=nlambda
    chords%lambdamin=lambdamin*10. !convert to angstroms
    chords%lambdamax=lambdamax*10. !convert to angstroms
    chords%dlambda=(chords%lambdamax-chords%lambdamin)/chords%nlambda

    print*,'Shot:   ',inputs%shot_number
    print*,'Time: ',int(inputs%time*1.d3),'ms'
    print*,'Runid:        ',trim(adjustl(inputs%runid))

  end subroutine read_inputs

  subroutine read_beam_grid
    use netcdf
    character(120)  :: filename
    integer         :: nx_varid,ny_varid,nz_varid
    integer         :: ncid,alpha_varid,beta_varid,o_varid,xc_varid,yc_varid,zc_varid

    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_inputs.cdf"
    print*,'---- loading beam grid ----'
    print*,filename

    !!OPEN netCDF file
    call check( nf90_open(filename, nf90_nowrite, ncid) )

    !!GET THE VARIDS
    call check( nf90_inq_varid(ncid, "nx", nx_varid) )
    call check( nf90_inq_varid(ncid, "ny", ny_varid) )
    call check( nf90_inq_varid(ncid, "nz", nz_varid) )
    call check( nf90_inq_varid(ncid, "xc", xc_varid) )
    call check( nf90_inq_varid(ncid, "yc", yc_varid) )
    call check( nf90_inq_varid(ncid, "zc", zc_varid) )
    call check( nf90_inq_varid(ncid, "alpha", alpha_varid) )
    call check( nf90_inq_varid(ncid, "beta", beta_varid) )
    call check( nf90_inq_varid(ncid, "origin", o_varid) )

    !!Read in Variables
    call check( nf90_get_var(ncid, nx_varid, beam_grid%nx) )
    call check( nf90_get_var(ncid, ny_varid, beam_grid%ny) )
    call check( nf90_get_var(ncid, nz_varid, beam_grid%nz) )
    call check( nf90_get_var(ncid, alpha_varid, beam_grid%alpha) )
    call check( nf90_get_var(ncid, beta_varid, beam_grid%beta) )
    call check( nf90_get_var(ncid, o_varid, beam_grid%origin(:)) )

    allocate(beam_grid%xc(beam_grid%nx)  &
         ,   beam_grid%yc(beam_grid%ny)  &
         ,   beam_grid%zc(beam_grid%nz))

    call check( nf90_get_var(ncid, xc_varid, beam_grid%xc(:)) )
    call check( nf90_get_var(ncid, yc_varid, beam_grid%yc(:)) )
    call check( nf90_get_var(ncid, zc_varid, beam_grid%zc(:)) )

    !!CLOSE netCDF FILE
    call check( nf90_close(ncid) )

    beam_grid%dr(1) = abs(beam_grid%xc(2)-beam_grid%xc(1))
    beam_grid%dr(2) = abs(beam_grid%yc(2)-beam_grid%yc(1))
    beam_grid%dr(3) = abs(beam_grid%zc(2)-beam_grid%zc(1))

    beam_grid%lwh(1) = abs(beam_grid%xc(beam_grid%nx) - beam_grid%xc(1)) + beam_grid%dr(1)
    beam_grid%lwh(2) = abs(beam_grid%yc(beam_grid%ny) - beam_grid%yc(1)) + beam_grid%dr(2)
    beam_grid%lwh(3) = abs(beam_grid%zc(beam_grid%nz) - beam_grid%zc(1)) + beam_grid%dr(3)

    beam_grid%volume = beam_grid%lwh(1)*beam_grid%lwh(3)*beam_grid%lwh(3)

    beam_grid%center(1) = (minval(beam_grid%xc) - 0.5*beam_grid%dr(1)) + 0.5*beam_grid%lwh(1)
    beam_grid%center(2) = (minval(beam_grid%yc) - 0.5*beam_grid%dr(2)) + 0.5*beam_grid%lwh(2)
    beam_grid%center(3) = (minval(beam_grid%zc) - 0.5*beam_grid%dr(3)) + 0.5*beam_grid%lwh(3)

    beam_grid%drmin  = minval(beam_grid%dr)
    beam_grid%dv     = beam_grid%dr(1)*beam_grid%dr(2)*beam_grid%dr(3)
    beam_grid%ntrack = beam_grid%nx+beam_grid%ny+beam_grid%nz
    beam_grid%ngrid  = beam_grid%nx*beam_grid%ny*beam_grid%nz

    print*,'Nx: ',beam_grid%nx
    print*,'Ny: ',beam_grid%ny
    print*,'Nz: ',beam_grid%nz
    print*,'dV: ',real(beam_grid%dv,float),'[cm^-3]'
  end subroutine read_beam_grid

  subroutine read_inter_grid
    use netcdf
    character(120)  :: filename
    integer         :: nr_varid,nw_varid,r2d_varid,w2d_varid,rr_varid,ww_varid
    integer         :: ncid

    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_inputs.cdf"
    print*,'---- loading interpolation grid ----'

    !!OPEN netCDF file
    call check( nf90_open(filename, nf90_nowrite, ncid) )

    !!GET THE VARIDS
    call check( nf90_inq_varid(ncid, "nr", nr_varid) )
    call check( nf90_inq_varid(ncid, "nw", nw_varid) )
    call check( nf90_inq_varid(ncid, "r2d", r2d_varid) )
    call check( nf90_inq_varid(ncid, "w2d", w2d_varid) )
    call check( nf90_inq_varid(ncid, "rr", rr_varid) )
    call check( nf90_inq_varid(ncid, "ww", ww_varid) )

    !!Read in Variables
    call check( nf90_get_var(ncid, nr_varid, inter_grid%nr) )
    call check( nf90_get_var(ncid, nw_varid, inter_grid%nw) )

    allocate(inter_grid%rr(inter_grid.nr),inter_grid%ww(inter_grid.nw))
    allocate(inter_grid%r2d(inter_grid.nr,inter_grid.nw))
    allocate(inter_grid%w2d(inter_grid.nr,inter_grid.nw))

    call check( nf90_get_var(ncid, rr_varid, inter_grid%rr(:)) )
    call check( nf90_get_var(ncid, ww_varid, inter_grid%ww(:)) )
    call check( nf90_get_var(ncid, r2d_varid, inter_grid%r2d(:,:)) )
    call check( nf90_get_var(ncid, w2d_varid, inter_grid%w2d(:,:)) )

    !!CLOSE netCDF FILE
    call check( nf90_close(ncid) )

    inter_grid%dr = abs(inter_grid%rr(2)-inter_grid%rr(1))
    inter_grid%dw = abs(inter_grid%ww(2)-inter_grid%ww(1))
    inter_grid%da = inter_grid%dr*inter_grid%dw

    print*,'Nr: ',inter_grid%nr
    print*,'Nw: ',inter_grid%nw
    print*,'dA: ',real(inter_grid%da,float),'[cm^-3]'
  end subroutine read_inter_grid

  subroutine read_beam
    use netcdf
    character(120)  :: filename
    integer         :: ncid,ab_varid,divy_varid,divz_varid,focy_varid,focz_varid,bn_varid
    integer         :: einj_varid,pinj_varid,sm_varid,bwra_varid,bwza_varid
    integer         :: xyzsrc_varid,xyzpos_varid
    real(double)    :: dis

    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_inputs.cdf"
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
    call check( nf90_inq_varid(ncid, "xyz_pos", xyzpos_varid))
    !!Read the variables
    call check( nf90_get_var(ncid, bn_varid, nbi%number) )
    call check( nf90_get_var(ncid, ab_varid, inputs%ab) )
    call check( nf90_get_var(ncid, bwra_varid, nbi%dv) )
    call check( nf90_get_var(ncid, bwza_varid, nbi%dw) )
    call check( nf90_get_var(ncid, divy_varid, nbi%divay(:)) )
    call check( nf90_get_var(ncid, divz_varid, nbi%divaz(:)) )
    call check( nf90_get_var(ncid, focy_varid, nbi%focy) )
    call check( nf90_get_var(ncid, focz_varid, nbi%focz) )
    call check( nf90_get_var(ncid, einj_varid, nbi%einj) )
    call check( nf90_get_var(ncid, pinj_varid, nbi%pinj) )
    call check( nf90_get_var(ncid, sm_varid, nbi%species_mix(:)) )
    call check( nf90_get_var(ncid, xyzsrc_varid, nbi%xyz_src(:)) )
    call check( nf90_get_var(ncid, xyzpos_varid, nbi%xyz_pos(:)) )

    !!CLOSE netCDF FILE
    call check( nf90_close(ncid) )

    nbi%vinj=sqrt(2.d0*nbi%einj*1.d3 &
            *e0/(inputs%ab*mass_u))*1.d2 !! [cm/s]

    dis = sqrt(sum(nbi%xyz_pos-nbi%xyz_src)**2.0))
    nbi%beta = asin((nbi%xyz_src(3)-nbi%xyz_pos(3))/dis)
    nbi%alpha = atan2(nbi%xyz_pos(2)-nbi%xyz_src(2),nbi%xyz_pos(1)-nbi%xyz_src(1))

    print*,'NBI #',nbi%number+1
    print*,'NBI power   :', real(nbi%pinj,float)
    print*,'NBI voltage :', real(nbi%einj,float)

  end subroutine read_beam

  subroutine read_chords
    use netcdf
    character(120)  :: filename
    integer(long)   :: i, j, k
    integer         :: ncid,xlens_varid,ylens_varid,zlens_varid,rad_varid
    integer         :: xlos_varid,ylos_varid,zlos_varid,ra_varid,rd_varid
    integer         :: sig_varid,h_varid,wght_varid,nchan_varid,chan_id_varid

    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_inputs.cdf"
    print*,'---- loading detector information ----'

    !!OPEN netCDF file
    call check( nf90_open(filename, nf90_nowrite, ncid) )

    !!GET THE VARIDS
    call check( nf90_inq_varid(ncid, "nchan", nchan_varid) )
    call check( nf90_inq_varid(ncid, "xlens", xlens_varid) )
    call check( nf90_inq_varid(ncid, "ylens", ylens_varid) )
    call check( nf90_inq_varid(ncid, "zlens", zlens_varid) )
    call check( nf90_inq_varid(ncid, "xlos", xlos_varid) )
    call check( nf90_inq_varid(ncid, "ylos", ylos_varid) )
    call check( nf90_inq_varid(ncid, "zlos", zlos_varid) )
    call check( nf90_inq_varid(ncid, "ra", ra_varid) )
    call check( nf90_inq_varid(ncid, "rd", rd_varid) )
    call check( nf90_inq_varid(ncid, "sigma_pi", sig_varid) )
    call check( nf90_inq_varid(ncid, "rlos", rad_varid) )
    call check( nf90_inq_varid(ncid, "h", h_varid) )
    call check( nf90_inq_varid(ncid, "chan_id", chan_id_varid) )
    call check( nf90_inq_varid(ncid, "los_wght", wght_varid) )

    !!READ IN THE NUMBER OF CHANNELS
    call check( nf90_get_var(ncid, nchan_varid, chords%nchan) )

    !!ALLOCATE SPACE
    allocate(chords%xyzlens(chords%nchan,3))
    allocate(chords%xyzlos(chords%nchan,3))
    allocate(chords%ra(chords%nchan))
    allocate(chords%rd(chords%nchan))
    allocate(chords%h(chords%nchan))
    allocate(chords%chan_id(chords%nchan))
    allocate(chords%sigma_pi(chords%nchan))
    allocate(chords%radius(chords%nchan))
    allocate(chords%los_wght(beam_grid%nx,beam_grid%ny,beam_grid%nz,chords%nchan))

    !!READ IN OTHER PARAMETERS
    call check( nf90_get_var(ncid, xlens_varid,   chords%xyzlens(:,1)) )
    call check( nf90_get_var(ncid, ylens_varid,   chords%xyzlens(:,2)) )
    call check( nf90_get_var(ncid, zlens_varid,   chords%xyzlens(:,3)) )
    call check( nf90_get_var(ncid, xlos_varid,    chords%xyzlos(:,1)) )
    call check( nf90_get_var(ncid, ylos_varid,    chords%xyzlos(:,2)) )
    call check( nf90_get_var(ncid, zlos_varid,    chords%xyzlos(:,3)) )
    call check( nf90_get_var(ncid, ra_varid,      chords%ra) )
    call check( nf90_get_var(ncid, rd_varid,      chords%rd) )
    call check( nf90_get_var(ncid, h_varid,       chords%h) )
    call check( nf90_get_var(ncid, chan_id_varid, chords%chan_id) )
    call check( nf90_get_var(ncid, sig_varid,     chords%sigma_pi) )
    call check( nf90_get_var(ncid, rad_varid,     chords%radius) )
    call check( nf90_get_var(ncid, wght_varid,    chords%los_wght(:,:,:,:)) )

    !!CLOSE netCDF FILE
    call check( nf90_close(ncid) )

    npa%nchan=0
    if (inputs%calc_npa.eq.0) then
       npa%npa_loop=1.
    else
       npa%npa_loop=10000.
    endif

    do i=1,chords%nchan
       if(chords%chan_id(i).eq.1) npa%nchan=npa%nchan+1
    enddo

    print*,'FIDA Channels: ',chords%nchan-npa%nchan
    print*,'NPA Channels:  ',npa%nchan

  end subroutine read_chords

  subroutine read_profiles
    use netcdf
    character(120):: filename
    integer       :: rhomax_var,te_var,ti_var,dene_var,deni_var,denp_var,nrho_var
    integer       :: zeff_var,rho_var,ncid,ai_var,btip_var,impc_var,omega_var

    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_inputs.cdf"
    print*,'---- loading profiles ----'

    !!OPEN netCDF file
    call check( nf90_open(filename, nf90_nowrite, ncid) )

    !!GET THE VARIDS
    call check( nf90_inq_varid(ncid, "ai", ai_var) )
    call check( nf90_inq_varid(ncid, "btipsign", btip_var) )
    call check( nf90_inq_varid(ncid, "impurity_charge", impc_var) )
    call check( nf90_inq_varid(ncid, "nrho", nrho_var) )
    call check( nf90_inq_varid(ncid, "rho_max", rhomax_var) )
    call check( nf90_inq_varid(ncid, "rho", rho_var) )
    call check( nf90_inq_varid(ncid, "te", te_var) )
    call check( nf90_inq_varid(ncid, "ti", ti_var) )
    call check( nf90_inq_varid(ncid, "dene", dene_var) )
    call check( nf90_inq_varid(ncid, "deni", deni_var) )
    call check( nf90_inq_varid(ncid, "denp", denp_var) )
    call check( nf90_inq_varid(ncid, "zeff", zeff_var) )
    call check( nf90_inq_varid(ncid, "omega", omega_var) )

    !!READ IN PARAMETERS
    call check( nf90_get_var(ncid, ai_var, inputs%ai) )
    call check( nf90_get_var(ncid, btip_var, inputs%btipsign) )
    call check( nf90_get_var(ncid, impc_var, inputs%impurity_charge) )
    call check( nf90_get_var(ncid, rhomax_var, profiles%rho_max) )
    call check( nf90_get_var(ncid, nrho_var, profiles%nrho) )

    allocate(profiles%rho(profiles%nrho))
    allocate(profiles%te(profiles%nrho))
    allocate(profiles%ti(profiles%nrho))
    allocate(profiles%dene(profiles%nrho))
    allocate(profiles%denp(profiles%nrho))
    allocate(profiles%deni(profiles%nrho))
    allocate(profiles%zeff(profiles%nrho))
    allocate(profiles%omega(profiles%nrho))

    call check( nf90_get_var(ncid, rho_var,  profiles%rho) )
    call check( nf90_get_var(ncid, te_var,   profiles%te) )
    call check( nf90_get_var(ncid, ti_var,   profiles%ti) )
    call check( nf90_get_var(ncid, dene_var, profiles%dene) )
    call check( nf90_get_var(ncid, deni_var, profiles%deni) )
    call check( nf90_get_var(ncid, denp_var, profiles%denp) )
    call check( nf90_get_var(ncid, zeff_var, profiles%zeff) )
    call check( nf90_get_var(ncid, omega_var, profiles%omega) )

    !!CLOSE netCDF FILE
    call check( nf90_close(ncid) )

    profiles%drho = abs(profiles%rho(2) - profiles%rho(1))

  end subroutine read_profiles

  subroutine read_equilibrium
    character(120)  :: filename
    integer       :: rho2d_var,br_var,bt_var,bw_var,er_var,et_var,ew_var

    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_inputs.cdf"
    print*,'---- loading equilibrium ----'

    !!OPEN netCDF file
    call check( nf90_open(filename, nf90_nowrite, ncid) )

    !!GET THE VARIDS
    call check( nf90_inq_varid(ncid, "rho2d", rho2d_var) )
    call check( nf90_inq_varid(ncid, "br", br_var) )
    call check( nf90_inq_varid(ncid, "bt", bt_var) )
    call check( nf90_inq_varid(ncid, "bw", bw_var) )
    call check( nf90_inq_varid(ncid, "er", er_var) )
    call check( nf90_inq_varid(ncid, "et", et_var) )
    call check( nf90_inq_varid(ncid, "ew", ew_var) )

    !!Allocate arrays
    allocate(equil%rho2d(inter_grid.nr,inter_grid.nw))
    allocate(equil%br(inter_grid.nr,inter_grid.nw))
    allocate(equil%bt(inter_grid.nr,inter_grid.nw))
    allocate(equil%bw(inter_grid.nr,inter_grid.nw))
    allocate(equil%er(inter_grid.nr,inter_grid.nw))
    allocate(equil%et(inter_grid.nr,inter_grid.nw))
    allocate(equil%ew(inter_grid.nr,inter_grid.nw))

    !!READ IN PARAMETERS
    call check( nf90_get_var(ncid, rho2d_var, equil%rho2d(:,:)) )
    call check( nf90_get_var(ncid, br_var, equil%br(:,:)) )
    call check( nf90_get_var(ncid, bt_var, equil%bt(:,:)) )
    call check( nf90_get_var(ncid, bw_var, equil%bw(:,:)) )
    call check( nf90_get_var(ncid, er_var, equil%er(:,:)) )
    call check( nf90_get_var(ncid, et_var, equil%et(:,:)) )
    call check( nf90_get_var(ncid, ew_var, equil%ew(:,:)) )

    !!CLOSE netCDF FILE
    call check( nf90_close(ncid) )

  end subroutine read_equilibrium

  subroutine read_tables
    character(120)  :: filename
    integer         :: n,m !! initial/final state
    integer(long)   :: nlev

   !-------------------ELECTRON EXCITATION/IONIZATION TABLE--------
    filename=trim(adjustl(root_dir))//"TABLES/qetable.bin"
    open(66,form='unformatted',file=filename,access='stream')
    read(66) tables%nr_te_qe
    read(66) tables%d_te_qe
    read(66) tables%nr_eb_qe
    read(66) tables%d_eb_qe
    read(66) nlev
    if(nlev.ne.nlevs) stop 'stop at "read atomics qetable"'
    allocate(tables%qe(nlevs+1,nlevs,tables%nr_eb_qe,tables%nr_te_qe))
    read(66) tables%qe(:,:,:,:)
    close(66)

    !-------------------Deuterium EXCITATION/IONIZATION/CX TABLE------
    filename=trim(adjustl(root_dir))//"TABLES/qptable.bin"
    open(66,form='unformatted',file=filename,access='stream')
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
    open(66,form='unformatted',file=filename,access='stream')
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
    open(66,form='unformatted',file=filename,access='stream')
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
    open(66,form='unformatted',file=filename,access='stream')
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

  subroutine read_fbm
    use netcdf
    character(120)  :: filename
    real(double) :: fdens
    integer :: i,j
    integer :: ncid,np_var,ne_var,fbm_var,e_var,p_var,denf_var

    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_inputs.cdf"
    print*,'---- loading fast ion distribution function ----'

    !!OPEN netCDF file
    call check( nf90_open(filename, nf90_nowrite, ncid) )

    !!GET VARIABLES VARID
    call check( nf90_inq_varid(ncid, "fbm_nenergy", ne_var) )
    call check( nf90_inq_varid(ncid, "fbm_npitch" , np_var) )
    call check( nf90_inq_varid(ncid, "fbm_energy" , e_var) )
    call check( nf90_inq_varid(ncid, "fbm_pitch"  , p_var) )
    call check( nf90_inq_varid(ncid, "fbm"        , fbm_var) )
    call check( nf90_inq_varid(ncid, "denf"       , denf_var) )

    !!GET VARIABLES
    call check( nf90_get_var(ncid, ne_var, fbm%nenergy) )
    call check( nf90_get_var(ncid, np_var, fbm%npitch) )

    !!ALLOCATE SPACE
    allocate(fbm%energy(fbm%nenergy))
    allocate(fbm%pitch(fbm%npitch))
    allocate(fbm%denf(inter_grid%nr,inter_grid%nw))
    allocate(fbm%fbm_norm(inter_grid%nr,inter_grid%nw))
    allocate(fbm%fbm(fbm%nenergy,fbm%npitch,inter_grid%nr,inter_grid%nw))

    !!GET REST OF VARIABLES
    call check( nf90_get_var(ncid, e_var,     fbm%energy(:)) )
    call check( nf90_get_var(ncid, p_var,     fbm%pitch(:)) )
    call check( nf90_get_var(ncid, fbm_var,   fbm%fbm(:,:,:,:) ))
    call check( nf90_get_var(ncid, denf_var,  fbm%denf(:,:)) )

    !!CLOSE netCDF FILE
    call check( nf90_close(ncid) )

    fbm%deb = abs(fbm%energy(2) - fbm%energy(1))
    fbm%dpitch = abs(fbm%pitch(2) - fbm%pitch(1))

    !!Normalize FBM
    do j=1,inter_grid.nw
      do i=1,inter_grid.nr
        fdens = sum(fbm%fbm(:,:,i,j)*fbm%deb*fbm%dpitch)
        if (fdens.gt.0.0) then
          fbm%fbm(:,:,i,j) = fbm%fbm(:,:,i,j)/fdens
          fbm%fbm_norm(i,j) = maxval(fbm%fbm(:,:,i,j))
          fbm%fbm(:,:,i,j) = fbm%fbm(:,:,i,j)/fbm%fbm_norm(i,j)
        endif
      enddo
    enddo
  end subroutine read_fbm

  subroutine write_birth_profile
    use netcdf
    integer           :: ncid,varid,dimids(5)
    integer           :: dimid1,x_dimid,y_dimid,z_dimid,e_dimid,p_dimid,nr_varid
    character(120)    :: filename
    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_birth.cdf"

    !Create netCDF file
    call check( nf90_create(filename, cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid) )

    !Define Dimensions
    call check( nf90_def_dim(ncid,"dim_001",1,dimid1) )
    call check( nf90_def_dim(ncid,"x",beam_grid%nx,x_dimid) )
    call check( nf90_def_dim(ncid,"y",beam_grid%ny,y_dimid) )
    call check( nf90_def_dim(ncid,"z",beam_grid%nz,z_dimid) )
    call check( nf90_def_dim(ncid,"energy",3,e_dimid) )
    call check( nf90_def_dim(ncid,"pitch",npitch_birth,p_dimid) )
    dimids = (/ x_dimid, y_dimid, z_dimid, e_dimid, p_dimid /)

    !Define variable
    call check( nf90_def_var(ncid,"n_nbi",NF90_INT,dimid1,nr_varid) )
    call check( nf90_def_var(ncid,"birth_dens",NF90_DOUBLE,dimids,varid) )

    !Add unit attributes
    call check( nf90_put_att(ncid,varid,"units","fast-ions/(s*dx*dy*dz*dE*dP)") )
    call check( nf90_enddef(ncid) )

    !Write data to file
    call check( nf90_put_var(ncid, nr_varid,inputs%n_nbi) )
    call check( nf90_put_var(ncid, varid, result%birth_dens) )

    !Close netCDF file
    call check( nf90_close(ncid) )

    print*, 'birth profile written to:      ',filename
  end subroutine write_birth_profile

  subroutine write_neutrals
    use netcdf
    integer         :: x_dimid,y_dimid,z_dimid,l_dimid,dimids(4)
    integer         :: ncid,full_varid,half_varid,third_varid,halo_varid
    integer         :: dimid1
    character(120)  :: filename
    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_neutrals.cdf"

    !Create netCDF file
    call check( nf90_create(filename, cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid) )

    !Define Dimensions
    call check( nf90_def_dim(ncid,"dim001",1,dimid1) )
    call check( nf90_def_dim(ncid,"x",beam_grid%nx,x_dimid) )
    call check( nf90_def_dim(ncid,"y",beam_grid%ny,y_dimid) )
    call check( nf90_def_dim(ncid,"z",beam_grid%nz,z_dimid) )
    call check( nf90_def_dim(ncid,"nlevs",nlevs,l_dimid) )
    dimids = (/ x_dimid, y_dimid, z_dimid, l_dimid /)

    !Define variables
    call check( nf90_def_var(ncid,"fdens",NF90_DOUBLE,dimids,full_varid) )
    call check( nf90_def_var(ncid,"hdens",NF90_DOUBLE,dimids,half_varid) )
    call check( nf90_def_var(ncid,"tdens",NF90_DOUBLE,dimids,third_varid) )
    call check( nf90_def_var(ncid,"halodens",NF90_DOUBLE,dimids,halo_varid) )

    !Add unit attributes
    call check( nf90_put_att(ncid,full_varid,"units","neutrals/(dx*dy*dz)") )
    call check( nf90_put_att(ncid,half_varid,"units","neutrals/(dx*dy*dz)") )
    call check( nf90_put_att(ncid,third_varid,"units","neutrals/(dx*dy*dz)") )
    call check( nf90_put_att(ncid,halo_varid,"units","neutrals/(dx*dy*dz)") )
    call check( nf90_enddef(ncid) )

    !Write to file
    call check( nf90_put_var(ncid, full_varid, result%neut_dens(:,:,:,:,nbif_type)) )
    call check( nf90_put_var(ncid, half_varid, result%neut_dens(:,:,:,:,nbih_type)) )
    call check( nf90_put_var(ncid, third_varid, result%neut_dens(:,:,:,:,nbit_type)) )
    call check( nf90_put_var(ncid, halo_varid, result%neut_dens(:,:,:,:,halo_type)) )

    !Close netCDF file
    call check( nf90_close(ncid) )

    print*, 'neutral density written to:      ',filename
  end subroutine write_neutrals

  subroutine write_npa
    use netcdf
    integer         :: nchan_dimid,e_dimid,c_dimid,ncid,dimids(3),dimid1,dimid3,maxcnt
    integer         :: ipos_varid,fpos_varid,v_varid,e_varid,f_varid,wght_varid,ei_varid,cnt_varid,nchan_varid
    character(120)  :: filename
    real(float), dimension(:,:,:),allocatable :: output
    real(float), dimension(:,:),allocatable :: output1
    real(float), dimension(:),allocatable :: output2
    maxcnt=maxval(npa%counter(:))
    allocate(output(maxcnt,3,npa%nchan))
    allocate(output1(maxcnt,npa%nchan))
    allocate(output2(npa%nchan))

    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_npa.cdf"

    !Create netCDF file
    call check( nf90_create(filename, cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid) )

    !Define Dimensions
    call check( nf90_def_dim(ncid,"dim001",1,dimid1) )
    call check( nf90_def_dim(ncid,"dim003",3,dimid3) )
    call check( nf90_def_dim(ncid,"energies",fbm%nenergy,e_dimid) )
    call check( nf90_def_dim(ncid,"nchan",npa%nchan,nchan_dimid) )
    call check( nf90_def_dim(ncid,"max_counts",maxcnt,c_dimid) )
    dimids = (/ c_dimid, dimid3, nchan_dimid /)

    !Define variables
    call check( nf90_def_var(ncid,"nchan",NF90_INT,dimid1,nchan_varid) )
    call check( nf90_def_var(ncid,"ipos",NF90_FLOAT,dimids,ipos_varid) )
    call check( nf90_def_var(ncid,"fpos",NF90_FLOAT,dimids,fpos_varid) )
    call check( nf90_def_var(ncid,"v",NF90_FLOAT,dimids,v_varid) )
    call check( nf90_def_var(ncid,"wght",NF90_FLOAT,(/ c_dimid,nchan_dimid /),wght_varid) )
    call check( nf90_def_var(ncid,"E",NF90_FLOAT,(/ c_dimid,nchan_dimid /),ei_varid) )
    call check( nf90_def_var(ncid,"flux",NF90_FLOAT,(/ e_dimid,nchan_dimid /),f_varid) )
    call check( nf90_def_var(ncid,"energy",NF90_FLOAT,e_dimid,e_varid) )
    call check( nf90_def_var(ncid,"counts",NF90_FLOAT,nchan_dimid,cnt_varid) )

    !Add unit attributes
    call check( nf90_put_att(ncid,ipos_varid,"units","cm") )
    call check( nf90_put_att(ncid,fpos_varid,"units","cm") )
    call check( nf90_put_att(ncid,v_varid,"units","cm/s") )
    call check( nf90_put_att(ncid,wght_varid,"units","particle # / marker") )
    call check( nf90_enddef(ncid) )

    !Write to file
    call check( nf90_put_var(ncid, nchan_varid, npa%nchan) )

    output(:,:,:)=real(npa%ipos(:maxcnt,:,:),float)
    call check( nf90_put_var(ncid, ipos_varid, output) )

    output(:,:,:)=real(npa%fpos(:maxcnt,:,:),float)
    call check( nf90_put_var(ncid, fpos_varid, output) )

    output(:,:,:)=real(npa%v(:maxcnt,:,:)   ,float)
    call check( nf90_put_var(ncid, v_varid, output) )

    output1(:,:)=real(npa%wght(:maxcnt,:) ,float)
    call check( nf90_put_var(ncid, wght_varid, output1) )

    output1(:,:)=real(npa%E(:maxcnt,:) ,float)
    call check( nf90_put_var(ncid, ei_varid, output1) )

    call check( nf90_put_var(ncid,f_varid,real(npa%flux(:,:),float)) )

    call check( nf90_put_var(ncid,e_varid,real(npa%energy(:),float)) )

    output2(:)=real(npa%counter(:),float)
    call check( nf90_put_var(ncid, cnt_varid,output2) )

    !Close netCDF file
    call check( nf90_close(ncid) )
    deallocate(output)
    deallocate(output1)
    deallocate(output2)
    print*, 'NPA data written to: ',filename
  end subroutine write_npa

  subroutine write_spectra
    use netcdf
    integer         :: i,nchan
    integer         :: ncid,nchan_varid,fida_varid,brems_varid,halo_varid
    integer         :: full_varid,half_varid,third_varid,lam_varid,rad_varid
    integer         :: chan_dimid,lam_dimid,dimid1,dimids(2)
    character(120)  :: filename
    real(double), dimension(:)  , allocatable :: lambda_arr,radius

    !! ------------------------ calculate wavelength array ------------------ !!
    allocate(lambda_arr(chords%nlambda))
    do i=1,chords%nlambda
       lambda_arr(i)=(i-0.5)*chords%dlambda*0.1d0 &
            +chords%lambdamin*0.1d0
    enddo

    nchan=0
    do i=1,chords%nchan
       if(chords%chan_id(i).eq.0) nchan=nchan+1
    enddo

    allocate(radius(nchan))
    do i=1,chords%nchan
       if(chords%chan_id(i).eq.0) then
           radius(i) = chords%radius(i)
       endif
    enddo
    !! convert [Ph/(s*wavel_bin*cm^2*all_directions)] to [Ph/(s*nm*sr*m^2)]!
    result%spectra(:,:,:)=result%spectra(:,:,:)/(0.1d0*chords%dlambda) &
         /(4.d0*pi)*1.d4

    !! write to file
    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_spectra.cdf"

    !Create netCDF file
    call check( nf90_create(filename, cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid) )

    !Define Dimensions
    call check( nf90_def_dim(ncid,"dim001",1,dimid1) )
    call check( nf90_def_dim(ncid,"lambda",chords%nlambda,lam_dimid) )
    call check( nf90_def_dim(ncid,"chan",nchan,chan_dimid) )
    dimids = (/ lam_dimid, chan_dimid /)

    !Define variables
    call check( nf90_def_var(ncid,"nchan",NF90_INT,dimid1,nchan_varid) )
    call check( nf90_def_var(ncid,"lambda",NF90_DOUBLE,lam_dimid,lam_varid) )
    call check( nf90_def_var(ncid,"full",NF90_DOUBLE,dimids,full_varid) )
    call check( nf90_def_var(ncid,"half",NF90_DOUBLE,dimids,half_varid) )
    call check( nf90_def_var(ncid,"third",NF90_DOUBLE,dimids,third_varid) )
    call check( nf90_def_var(ncid,"halo",NF90_DOUBLE,dimids,halo_varid) )
    call check( nf90_def_var(ncid,"brems",NF90_DOUBLE,dimids,brems_varid) )
    call check( nf90_def_var(ncid,"fida",NF90_DOUBLE,dimids,fida_varid) )
    call check( nf90_def_var(ncid,"radius",NF90_DOUBLE,chan_dimid,rad_varid) )

    !Add unit attributes
    call check( nf90_put_att(ncid,lam_varid,"units","nm") )
    call check( nf90_put_att(ncid,full_varid,"units","Ph/(s*nm*sr*m^2)") )
    call check( nf90_put_att(ncid,half_varid,"units","Ph/(s*nm*sr*m^2)") )
    call check( nf90_put_att(ncid,third_varid,"units","Ph/(s*nm*sr*m^2)") )
    call check( nf90_put_att(ncid,halo_varid,"units","Ph/(s*nm*sr*m^2)") )
    call check( nf90_put_att(ncid,brems_varid,"units","Ph/(s*nm*sr*m^2)") )
    call check( nf90_put_att(ncid,fida_varid,"units","Ph/(s*nm*sr*m^2)") )
    call check( nf90_put_att(ncid,rad_varid,"units","cm") )

    call check( nf90_enddef(ncid) )

    !Write to file
    call check( nf90_put_var(ncid, nchan_varid, nchan) )
    call check( nf90_put_var(ncid, lam_varid, lambda_arr) )
    call check( nf90_put_var(ncid, rad_varid, radius) )
    call check( nf90_put_var(ncid, full_varid, result%spectra(:,:,nbif_type)) )
    call check( nf90_put_var(ncid, half_varid, result%spectra(:,:,nbih_type)) )
    call check( nf90_put_var(ncid, third_varid, result%spectra(:,:,nbit_type)) )
    call check( nf90_put_var(ncid, halo_varid, result%spectra(:,:,halo_type)) )
    call check( nf90_put_var(ncid, brems_varid, result%spectra(:,:,brems_type)) )
    call check( nf90_put_var(ncid, fida_varid, result%spectra(:,:,fida_type)) )

    !Close netCDF file
    call check( nf90_close(ncid) )

    deallocate(lambda_arr)
    deallocate(radius)
    print*, 'Spectra written to: ', filename
  end subroutine write_spectra

  subroutine read_neutrals
    use netcdf
    character(120)  :: filename
    integer :: ncid,full_var,half_var,third_var,halo_var

    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_neutrals.cdf"
    print*,'---- loading neutrals ----',filename

       !!OPEN netCDF file
    call check( nf90_open(filename, nf90_nowrite, ncid) )

    !!GET VARIABLE
    call check( nf90_inq_varid(ncid, "fdens", full_var) )
    call check( nf90_inq_varid(ncid, "hdens", half_var) )
    call check( nf90_inq_varid(ncid, "tdens", third_var) )
    call check( nf90_inq_varid(ncid, "halodens", halo_var) )
    call check( nf90_get_var(ncid, full_var,result%neut_dens(:,:,:,:,nbif_type)) )
    call check( nf90_get_var(ncid, half_var,result%neut_dens(:,:,:,:,nbih_type)) )
    call check( nf90_get_var(ncid, third_var,result%neut_dens(:,:,:,:,nbit_type)) )
    call check( nf90_get_var(ncid, halo_var,result%neut_dens(:,:,:,:,halo_type)) )

    !!CLOSE netCDF FILE
    call check( nf90_close(ncid) )
   end subroutine read_neutrals


  !*****************************************************************************
  !*****************************Geometry Routines*******************************
  !*****************************************************************************
  subroutine Ru(u,theta,R)
    real(double), dimension(3), intent(in)    :: u
    real(double), intent(in)                  :: theta
    real(double), dimension(3,3), intent(out) :: R
    real(double), dimension(3,3)              :: I,ux,uxu

    ! Identity matrix
    I(1,1) = 1.0 ; I(1,2) = 0.0 ; I(1,3) = 0.0
    I(2,1) = 0.0 ; I(2,2) = 1.0 ; I(2,3) = 0.0
    I(3,1) = 0.0 ; I(3,2) = 0.0 ; I(3,3) = 1.0

    ux(1,1) = 0.0  ; ux(1,2) =-u(3) ; ux(1,3) = u(2)
    ux(2,1) = u(3) ; ux(2,2) = 0.0  ; ux(2,3) =-u(1)
    ux(3,1) =-u(2) ; ux(3,2) = u(1) ; ux(3,3) = 0.0

    uxu(1,1) = u(1)**2.0  ; uxu(1,2) = u(1)*u(2)  ; uxu(1,3) = u(1)*u(3)
    uxu(2,1) = u(1)*u(2)  ; uxu(2,2) = u(2)**2.0  ; uxu(2,3) = u(2)*u(3)
    uxu(3,1) = u(1)*u(3)  ; uxu(3,2) = u(2)*u(3)  ; uxu(3,3) = u(3)**2.0

    R(:,:) = cos(theta)*I(:,:) + sin(theta)*ux(:,:) + (1.0-cos(theta))*uxu(:,:)
  end subroutine Ru

  subroutine Rz(theta,R)
    real(double), intent(in)                  :: theta
    real(double), dimension(3,3), intent(out) :: R
    real(double), dimension(3)                :: khat
    khat = (/0.0,0.0,1.0/)
    call Ru(khat,theta,R)
  end subroutine Rz

  subroutine xyz_to_uvw(xyz,uvw,origin_in,alpha_in,beta_in)
    real(double), dimension(3), intent(in)             :: xyz
    real(double), dimension(3), intent(out)            :: uvw
    real(double), dimension(3), intent(in), optional   :: origin_in
    real(double), intent(in), optional                 :: alpha_in
    real(double), intent(in), optional                 :: beta_in
    real(double), dimension(3)                         :: origin,jhat,jhat_p
    real(double), dimension(3,3)                       :: R_u,R_z,R
    real(double)                                       :: alpha,beta

    if(present(origin_in)) then
      origin(:) = origin_in(:)
    else
      origin(:) = (/0.0,0.0,0.0/)
    endif

    if(present(alpha_in)) then
      alpha = alpha_in
    else
      alpha = beam_grid%alpha
    endif

    if(present(beta_in)) then
      beta = beta_in
    else
      beta = beam_grid%beta
    endif

    jhat = (/0.0,1.0,0.0/)

    call Rz(alpha,R_z)
    jhat_p = matmul(R_z,jhat)

    call Ru(jhat_p,beta,R_u)
    R = matmul(R_u(,R_z)

    uvw = matmul(R,xyz()
    uvw = uvw + origin

  end subroutine xyz_to_uvw

  subroutine uvw_to_xyz(uvw,xyz,origin_in,alpha_in,beta_in)
    real(double), dimension(3), intent(in)             :: uvw
    real(double), dimension(3), intent(out)            :: xyz
    real(double), dimension(3), intent(in), optional   :: origin_in
    real(double), intent(in), optional                 :: alpha_in
    real(double), intent(in), optional                 :: beta_in
    real(double), dimension(3)                         :: origin, jhat,jhat_p,uvw_p
    real(double), dimension(3,3)                       :: R,R_z,R_u
    real(double)                                       :: alpha,beta

    if(present(origin_in)) then
      origin = origin_in(:)
    else
      origin = beam_grid%origin
    endif

    if(present(alpha_in)) then
      alpha = alpha_in
    else
      alpha = beam_grid%alpha
    endif

    if(present(beta_in)) then
      beta = beta_in
    else
      beta = beam_grid%beta
    endif

    uvw_p = uvw - origin

    jhat = (/0.0,1.0,0.0/)

    call Rz(alpha,R_z))
    jhat_p = matmul(R_z,jhat)

    call Ru(jhat_p,beta,R_u)
    R = transpose(matmul(R_u,R_z))

    xyz = matmul(R,uvw_p)

  end subroutine uvw_to_xyz

  subroutine chord_coor(xyz_i,xyz_f,origin,pos_in,pos_out)
    real(double), dimension(3)                    :: pos_in,pos_out
    real(double), dimension(3)                    :: xyz_i,xyz_f,origin
    real(double)                                  :: phi,theta,xp,yp,zp

    phi=atan2((xyz_f(2)-xyz_i(2)),(xyz_f(1)-xyz_i(1)))
    theta= -1*atan2(sqrt((xyz_f(1)-xyz_i(1))**2.0 + (xyz_f(2)-xyz_i(2))**2.0),(xyz_f(3)-xyz_i(3)))

    xp=(pos_in(1)-origin(1))*cos(phi)+(pos_in(2)-origin(2))*sin(phi)
    yp=-(pos_in(1)-origin(1))*sin(phi)+(pos_in(2)-origin(2))*cos(phi)
    zp=pos_in(3)-origin(3)
    pos_out(1)=xp*cos(theta)+zp*sin(theta)
    pos_out(2)=yp
    pos_out(3)=-xp*sin(theta)+zp*cos(theta)
  end subroutine chord_coor

  subroutine inv_chord_coor(xyz_i,xyz_f,origin,pos_in,pos_out)
    real(double), dimension(3)                    :: pos_in,pos_out
    real(double), dimension(3)                    :: xyz_i,xyz_f,origin
    real(double)                                  :: phi,theta

    phi=atan2((xyz_f(2)-xyz_i(2)),(xyz_f(1)-xyz_i(1)))
    theta= -1*atan2(sqrt((xyz_f(1)-xyz_i(1))**2.0 + (xyz_f(2)-xyz_i(2))**2.0),(xyz_f(3)-xyz_i(3)))
    pos_out(1)=origin(1)+pos_in(1)*cos(phi)*cos(theta) - pos_in(3)*cos(phi)*sin(theta) - pos_in(2)*sin(phi)
    pos_out(2)=origin(2)+pos_in(2)*cos(phi) + pos_in(1)*cos(theta)*sin(phi) - pos_in(3)*sin(phi)*sin(theta)
    pos_out(3)=origin(3)+pos_in(3)*cos(theta) + pos_in(1)*sin(theta)

  end subroutine inv_chord_coor

  subroutine grid_intersect(r0,v0,length,r_enter,r_exit)
    real(double), dimension(3), intent(in)  :: r0
    real(double), dimension(3), intent(in)  :: v0
    real(double), dimension(3), intent(out) :: r_enter
    real(double), dimension(3), intent(out) :: r_exit
    real(double), intent(out)               :: length
    real(double), dimension(3,6)            :: ipnts
    integer,      dimension(6)              :: side_inter
    integer,      dimension(2)              :: ind
    integer,                                :: i,j
    real(double)                            :: dis1,dis2

    do i=1,6
      j = int(ceiling(i/2.0))
      if j.eq.1 ind = (/2,3/)
      if j.eq.2 ind = (/1,3/)
      if j.eq.3 ind = (/1,2/)
      if (abs(v0[j]).gt.0.0) then
        ipnts(:,i) = r0 + v0*( ( (beam_grid%center(j) + &
                     (mod(i,2) - 0.5)*beam_grid%lwh(j)) - r0(j))/v0(j) )
        if ((abs(ipnts(ind(1),i) - beam_grid%center(ind(1))).le.(0.5*beam_grid%lwh(ind(1)))).and. &
            (abs(ipnts(ind(2),i) - beam_grid%center(ind(2))).le.(0.5*beam_grid%lwh(ind(2))))) then
          side_inter(i) = 1
        endif
      endif
    enddo

    length = 0.0
    if (sum(side_inter).lt.2) then
      r_enter = r0
      r_exit  = r0
    else
      i=1
      do while (i.le.6)
        if(side_inter(i).eq.1) exit
        i=i+1
      enddo
      j=i
      do while (sqrt(sum((ipnts(:,i)-ipnts(:,j)**2.0))).le.0.001)
        do while (j.le.6)
          j=j+1
          if(side_inter(j).eq.1) exit
        enddo
      enddo
      dis1 = sqrt(sum((ipnts(:,i)-r0)**2.0))
      dis2 = sqrt(sum((ipnts(:,j)-r0)**2.0))
      if (dis1.lt.dis2) then
        r_enter = ipnts(:,i)
        r_exit  = ipnts(:,j)
      else
        r_enter = ipnts(:,j)
        r_exit  = ipnts(:,i)
      endif
      length = abs(dis1-dis2)
    endif

  end subroutine grid_intersect

  subroutine get_index(pos,ac)
    real(double),  dimension(3), intent(in)  :: pos
    integer(long), dimension(3), intent(out) :: ac
    real(double),  dimension(3)              :: mini
    integer(long), dimension(3)              :: maxind
    integer,                                 :: i

    maxind(1) = beam_grid%nx
    maxind(2) = beam_grid%ny
    maxind(3) = beam_grid%nz

    mini(1) = minval(beam_grid%xc) - 0.5*beam_grid%dr(1)
    mini(2) = minval(beam_grid%yc) - 0.5*beam_grid%dr(2)
    mini(3) = minval(beam_grid%zc) - 0.5*beam_grid%dr(3)

    do i=1,3
      ac(i) = ceiling((pos(i)-mini(i))/beam_grid%dr(i))
      if (ac(i).gt.maxind(i)) ac(i)=maxind(i)
      if (ac(i).lt.1) ac(i)=1
    enddo

  end subroutine get_index

  subroutine track3D(rin, vin, tcell, icell, pcell, ncell)
    !!track computes the path of a neutral through a the beam grid
    real(double), dimension(:)  , intent(in)   :: rin  ! initial position
    real(double), dimension(:)  , intent(in)   :: vin  ! velocitiy
    integer                     , intent(out)  :: ncell! number of cells
    real(double), dimension(:)  , intent(out)  :: tcell! time per cell
    integer,dimension(:,:), intent(out)        :: icell! cell indices
    real(double), dimension(:,:), intent(out)  :: pcell! mean position in cell
    integer                    :: cc,i   !!step number along the track
    integer,dimension(3)       :: ac     !!indices of the cells
    real(double), dimension(3) :: dt_arr !!time to cell boundary
    real(double), dimension(3) :: vn     !!velocitiy that can be changed
    real(double), dimension(3) :: ri     !!position of ray
    real(double), dimension(3) :: sgn    !!sign
    integer,      dimension(3) :: gdims  !!grid dims
    integer,      dimension(1) :: minpos

    tcell=0.d0 ; icell=0 ;  ; pcell=0.d0 ; ncell=0
    vn(:)=vin(:) ;  ri(:)=rin(:) ; sgn=0

    gdims(1) = beam_grid%nx
    gdims(2) = beam_grid%ny
    gdims(3) = beam_grid%nz

    !! define actual cell
    call get_index(ri,ac)

    do i=1,3
      if (vn(i).gt.0.0) sgn(i) = 1
      if (vn(i).lt.0.0) sgn(i) =-1
      if (vn(i).eq.0.0) vn(i)  = 1.0d-30
    enddo

    cc=0
    do i=1,beam_grid%ntrack
      dt_arr(1) = ( (beam_grid%xc(ac(1)) + 0.5*sgn(1)*beam_grid%dr(1)) - ri(1))/vn(1)
      dt_arr(2) = ( (beam_grid%yc(ac(2)) + 0.5*sgn(2)*beam_grid%dr(2)) - ri(2))/vn(2)
      dt_arr(3) = ( (beam_grid%zc(ac(3)) + 0.5*sgn(3)*beam_grid%dr(3)) - ri(3))/vn(3)
      minpos = minloc(dt_arr)
      pcell(:,cc) = ri + 0.5*dt_arr(minpos(1))*vn !! path midpoint in cell
      ri = ri + dt_arr(minpos(1))*vn
      tcell(cc) = dt_arr(minpos(1))
      icell(:,cc) = ac
      ac(minpos(1)) = ac(minpos(1)) + sgn(minpos(1))
      cc = cc+1
      if (ac(minpos(1)).gt.gdims(minpos(1))) exit
      if (ac(minpos(1)).lt.1) exit
    enddo
    ncell = cc
  end subroutine track3D

  subroutine hit_npa_detector(pos,vel,det)
    integer                                       :: i,det,cnt
    real(double), dimension(3)                    :: pos,vel,xyz_i,xyz_f
    real(double), dimension(3)                    :: rpos,rvel
    real(double)                                  :: ra,rd,xa,ya,xd,yd

    det=0
    cnt=0
    loop_over_chan: do i=1,chords%nchan
      !!only do npa channels
      if(chords%chan_id(i).ne.1) cycle loop_over_chan
      cnt=cnt+1
      !!transform to chord coordinates
      xyz_i(1)=chords%xyzlens(i,1)
      xyz_i(2)=chords%xyzlens(i,2)
      xyz_i(3)=chords%xyzlens(i,3)
      xyz_f(1)=chords%xyzlos(i,1)
      xyz_f(2)=chords%xyzlos(i,2)
      xyz_f(3)=chords%xyzlos(i,3)

      call chord_coor(xyz_i,xyz_f,xyz_i,pos,rpos)
      call chord_coor(xyz_i,xyz_f,xyz_i*0.0,vel,rvel)

      !!if a neutral particle is moving in the opposite direction it will not hit detector
      if(rvel(3).ge.0) cycle loop_over_chan
      xa=rpos(1)+rvel(1)*(-rpos(3)/rvel(3))
      ya=rpos(2)+rvel(2)*(-rpos(3)/rvel(3))
      ra=sqrt(xa**2.0 + ya**2.0)
      xd=rpos(1)+rvel(1)*((-chords%h(i)-rpos(3))/rvel(3))
      yd=rpos(2)+rvel(2)*((-chords%h(i)-rpos(3))/rvel(3))
      rd=sqrt(xd**2.0 + yd**2.0)

      !!if a neutral particle pass through both the aperture and detector then it count it
      if( rd.le.chords%rd(i).and.ra.le.chords%ra(i) ) then
        det=cnt
        exit loop_over_chan
      endif
    enddo loop_over_chan

  end subroutine hit_npa_detector

  !*****************************************************************************
  !*************************Profiles and Fields Routines************************
  !*****************************************************************************
  subroutine calc_perp_vectors(b,a,c)
    !!Returns normalized vectors that are perpendicular to b
    real(double), dimension(3),intent(in)  :: b
    real(double), dimension(3),intent(out) :: a,c
    real(double), dimension(3)             :: bnorm

    bnorm=b/sqrt(dot_product(b,b))

    if (abs(bnorm(3)).eq.1) then
      a=(/1.d0,0.d0,0.d0/)
      c=(/0.d0,1.d0,0.d0/)
    else
      if (bnorm(3).eq.0.) then
        a=(/0.d0,0.d0,1.d0/)
        c=(/bnorm(2),-bnorm(1), 0.d0/)/sqrt(bnorm(1)**2+bnorm(2)**2)
      else
        a=(/bnorm(2),-bnorm(1),0.d0/)/sqrt(bnorm(1)**2+bnorm(2)**2)
        c=(/ a(2) , -a(1) , (a(1)*bnorm(2)-a(2)*bnorm(1))/bnorm(3) /)
        c=c/sqrt(dot_product(c,c))
      endif
    endif
  end subroutine calc_perp_vectors

  subroutine get_profiles(xyz,ti,te,dene,deni,denp,vrot)

  end subroutine get_profiles

  !*****************************************************************************
  !**************************Atomic Physics Routines****************************
  !*****************************************************************************
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

  subroutine neut_rates(denn,vi,vn,rates)
    !GET neutralization rate from tables
    real(double),  dimension(nlevs),intent(in):: denn !!density of neutrals cm-3
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
  end subroutine neut_rates

  subroutine colrad(ac,vn,dt,states,photons,neut_type,nlaunch,ri,det)
    !colrad solves the collisional-radiative balance equations
    real(double) , dimension(:),  intent(in)    :: vn  !!velocitiy (cm/s)
    real(double)               ,  intent(in)    :: dt  !!time interval in cell
    integer                    ,  intent(in)    :: neut_type!!type of neutral
    real(double)               ,  intent(in)    :: nlaunch !! nr of markers
    integer, dimension(:)      ,  intent(in)    :: ac   !!actual cell
    real(double), dimension(:) ,  intent(inout) :: states  !!density of states
    real(double),                 intent(out)   :: photons !!emitted photons
    real(double) , dimension(:),  intent(in), optional  :: ri  !! start position
    integer,                      intent(in), optional  :: det  !! npa detector index
    !! ---- to determine rate coefficients ---- !
    real(double), dimension(nlevs+1,nlevs)     :: qp !! Proton rate coefficants
    real(double), dimension(nlevs+1,nlevs)     :: qi !! Impurity rate coefficant
    real(double), dimension(nlevs+1,nlevs)     :: qe !! Electron rate coefficant
    real(double),dimension(nlevs,nlevs)        :: matrix  !! Matrix
    real(double), dimension(3)  :: vrot           !! Rotation velocity of plasma
    real(double)                :: vnet_square    !! netto velocity of neutrals
    real(double)                :: ti,te          !! Ion/electron temperature
    real(double)                :: denp,dene,deni !! P/impurity/electron density
    real(double)                :: eb             !! Energy of the fast neutral
    real(double)                :: rho            !! Normalized flux coord of cell
    integer                     :: ebi, tii,tei   !! bin postions in arrays
    !! ---- Solution of differential equation  ---- !
    real(double),   dimension(nlevs,nlevs)  :: eigvec, eigvec_inv
    real(double),   dimension(nlevs)        :: eigval, coef
    real(double),   dimension(nlevs)        :: exp_eigval_dt
    real(double),   dimension(nlevs)        :: dens !! Density of neutrals
    real(double)                            :: iflux !!Initial total flux
    real(double)                            :: dflux !! change of flux
    integer                                 :: n !! counter
    real(double),   dimension(3)            :: b_norm  !! pitch of particle
    real(double)                            :: ptch  !! pitch of particle
    integer                                 :: ipitch !! index of pitch
    real(double)                            :: dray !! (for NPA)
    real(double), dimension(3)              :: ray     !! ray towards NPA
    real(double), dimension(1)              :: ienergy !! (for NPA)
    photons=0.d0
    iflux=sum(states)
    !! --------------- Check if inputs are valid for colrad -------------- !!
    if(iflux.lt.colrad_threshold .and. inputs%calc_npa.eq.0)then
      ! print*, 'threshold!',ac
       return
    endif
    denp=cell(ac(1),ac(2),ac(3))%plasma%denp
    dene=cell(ac(1),ac(2),ac(3))%plasma%dene
    deni=cell(ac(1),ac(2),ac(3))%plasma%deni
    ti=cell(ac(1),ac(2),ac(3))%plasma%ti
    te=cell(ac(1),ac(2),ac(3))%plasma%te
    vrot=cell(ac(1),ac(2),ac(3))%plasma%vrot
    rho=cell(ac(1),ac(2),ac(3))%rho

    !! If cells rho is greater than rho_max then stop simulation
    if(rho.gt.profiles%rho_max) then
       if(neut_type.le.3.and.neut_type.ne.0)then  !! Store density for NBI simulation!
          dens(:)=states*dt/nlaunch!![neutrals/(cm^3)]!!
          result%neut_dens(ac(1),ac(2),ac(3),:,neut_type)= &
              result%neut_dens(ac(1),ac(2),ac(3),:,neut_type)+dens(:)
       endif
       if(inputs%calc_npa.eq.1.and.neut_type.eq.fida_type)then !! NPA simulation !!
          dray=sqrt(dot_product(ri-chords%xyzlens(1,:) &
               ,ri-chords%xyzlens(1,:)))
          ray=ri(:)+vn(:)/sqrt(dot_product(vn,vn))*dray
          !$OMP CRITICAL(col_rad_npa)
          npa%counter(det)=npa%counter(det)+1
          if(npa%counter(det).gt.inputs%n_npa)stop'too many neutrals'
          npa%v(npa%counter(det),:,det)=vn(:)
          npa%wght(npa%counter(det),det)=sum(states)/nlaunch*beam_grid%dv/npa%npa_loop !![neutrals/s]
          npa%E(npa%counter(det),det) = inputs%ab*v_to_E*dot_product(vn,vn)
          npa%ipos(npa%counter(det),:,det)=ri(:)
          npa%fpos(npa%counter(det),:,det)=ray(:)
          ienergy=minloc(abs(npa%energy - npa%E(npa%counter(det),det)))
          npa%flux(ienergy(1),det)=npa%flux(ienergy(1),det)+ npa%wght(npa%counter(det),det)/fbm%deb
          !$OMP END CRITICAL(col_rad_npa)
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
    call eigen(nlevs, matrix, eigvec, eigval)
    call matinv(eigvec, eigvec_inv)
    coef = matmul(eigvec_inv, states)!coeffs determined from states at t=0
    exp_eigval_dt = exp(eigval*dt)   ! to improve speed (used twice)
    do n=1,nlevs
       if(eigval(n).eq.0.0) eigval(n)=eigval(n)+1 !protect against dividing by zero
    enddo
    states(:) = matmul(eigvec, coef * exp_eigval_dt)  ![neutrals/cm^3/s]!
    dens(:)   = matmul(eigvec,coef*(exp_eigval_dt-1.d0)/eigval)/nlaunch

    !! - Protect against negative values
    if ((minval(states).lt.0).or.(minval(dens).lt.0)) then
        do n=1,nlevs
            if(states(n).lt.0) states(n)=0.0
            if(dens(n).lt.0) dens(n)=0.0
        enddo
    endif

    if(neut_type.ne.0) then
      !$OMP CRITICAL(col_rad)
      result%neut_dens(ac(1),ac(2),ac(3),:,neut_type)= &
           result%neut_dens(ac(1),ac(2),ac(3),:,neut_type)+dens(:)![neutrals/cm^3]
      !$OMP END CRITICAL(col_rad)
    endif
    if(inputs%calc_birth.eq.1)then
       if(neut_type.le.3.and.neut_type.ne.0)then
          b_norm=cell(ac(1),ac(2),ac(3))%plasma%b_norm(:)
          ptch=dot_product(vn,b_norm)/sqrt(dot_product(vn,vn))
          ipitch=int((ptch+1.)/(2./npitch_birth))+1
          dflux=(iflux-sum(states))*beam_grid%dv/nlaunch !! [fast-ions/s]
          !$OMP CRITICAL(col_rad2)
          result%birth_dens(ac(1),ac(2),ac(3),neut_type,ipitch)= &
               result%birth_dens(ac(1),ac(2),ac(3),neut_type,ipitch) + dflux
          !$OMP END CRITICAL(col_rad2)
       endif
    endif
    !! -------------------- determine photon flux ------------------------ !!
    photons=dens(3)*tables%einstein(2,3) !! - [Ph/(s*cm^3)] - !!
  end subroutine colrad

  subroutine spectrum(vi,ac,pos,photons,neut_type,wavout,intout)
    !!spectrum.pro computes the wavelengths of emitted photons
    real(double), dimension(:), intent(in) :: vi!!velocitiy of neutral [cm/s]
    integer,dimension(:), intent(in)       :: ac  !!actual cell
    integer             , intent(in)       :: neut_type!!type of neutral
    real(double)              , intent(in) :: photons !! photons from colrad
    real(double), dimension(3), intent(in) :: pos    !! mean position in cell
    real(double), dimension(n_stark), intent(out), optional :: intout!!intensity
    real(double), dimension(n_stark), intent(out), optional :: wavout !!wavelength
    real(double)               :: lambda_Doppler , cos_los_Efield, E
    real(double),dimension(n_stark) ::intens!!intensity vector
    real(double), dimension(n_stark)::wavel !!wavelength vector[A
    real(double), dimension(3) :: vp  !!unit vector of sight line
    real(double), dimension(3) :: vn  ! vi in m/s
    real(double), dimension(3) :: efield  !E-field (static + vxB)
    real(double), dimension(3) :: bfield  !B-field
    integer                    :: i, ichan, bin,cnt  !counter, wavelengths bins
    integer,parameter,dimension(n_stark)::stark_sign= +1*stark_sigma -1*stark_pi
               !sign in stark intensity formula:
               !- for Pi (linear), + for Sigma (circular)
    cnt=0
    loop_over_channels: do ichan=1,chords%nchan
       if(chords%chan_id(ichan).ne.0)cycle loop_over_channels
       cnt=cnt+1
       if(cell(ac(1),ac(2),ac(3))%los_wght(ichan).le.0.)cycle loop_over_channels
       !! vector directing towards the optical head
       vp(1)=pos(1)-chords%xyzlens(ichan,1)
       vp(2)=pos(2)-chords%xyzlens(ichan,2)
       vp(3)=pos(3)-chords%xyzlens(ichan,3)
       vp=vp/sqrt(dot_product(vp,vp))
       ! Calculate Doppler shift
       vn=vi*0.01d0 ! [m/s]
       lambda_Doppler = lambda0*(1.d0 + dot_product(vn,vp)/c0)
       !! Calculate Stark Splitting
       ! Calcualate E-field
       bfield(:) = cell(ac(1),ac(2),ac(3))%plasma%b_norm(:) &
            * cell(ac(1),ac(2),ac(3))%plasma%b_abs
       efield(:) = cell(ac(1),ac(2),ac(3))%plasma%E
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
       !intens(4:6) = intens(4:6)*inputs%sigma_pi_ratio
       where (stark_sigma .eq. 1)
          intens = intens * chords%sigma_pi(ichan)
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
          bin=int(((wavel(i)-chords%lambdamin)/chords%dlambda)+.5)
          if (bin.lt.1)            bin = 1
          if (bin.gt.chords%nlambda) bin = chords%nlambda
          !$OMP CRITICAL(spec_trum)
          result%spectra(bin,cnt,neut_type)= &
               result%spectra(bin,cnt,neut_type) &
               +intens(i)*cell(ac(1),ac(2),ac(3))%los_wght(cnt)
          !$OMP END CRITICAL(spec_trum)
       enddo
    enddo loop_over_channels
  end subroutine spectrum

  !*****************************************************************************
  !*************************Random Number Generators****************************
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

  !*****************************************************************************
  !***************************Monte Carlo Routines******************************
  !*****************************************************************************
  subroutine get_nlaunch(nr_markers,papprox,papprox_tot,nlaunch)
    !! routine to define the number of MC particles started in one cell
    integer                       , intent(in)    :: nr_markers
    real(double), dimension(:,:,:), intent(in)    :: papprox
    real(double)                  , intent(in)    :: papprox_tot
    real(double), dimension(:,:,:), intent(out)   :: nlaunch
    integer  :: i,j,k,cc
    real(double), dimension(:), allocatable :: randomu
    do i=1,1000
       nlaunch(:,:,:)=papprox(:,:,:)/papprox_tot*nr_markers*(1.+i*0.01)
       if(sum(nlaunch).gt.nr_markers)exit
    enddo
    allocate(randomu(count(nlaunch.gt.0)))
    call randu(randomu)
    cc=1
    do k = 1, beam_grid%nz
       do j = 1, beam_grid%ny
          do i = 1, beam_grid%nx
             if(nlaunch(i,j,k).gt.0.)then
                if(mod(nlaunch(i,j,k),1.).gt.randomu(cc))then
                   nlaunch(i,j,k)=nlaunch(i,j,k)+1.
                endif
                cc=cc+1
             endif
          enddo
       enddo
    enddo
    nlaunch=floor(nlaunch)
    deallocate(randomu)
  end subroutine get_nlaunch

  subroutine mc_fastion(ac,ri,vi)
    !!IN: ac,   OUT: ri,vi
    !!mcbeam computes monte carlo velocity of fast ions from cell%fbm
    integer,      dimension(3), intent(in) :: ac !ind of actual cell
    real(double), dimension(3), intent(out):: vi !velocity [cm/s]
    real(double), dimension(3), intent(out):: ri !starting position
    real(double), dimension(3)             :: a,b,c ! vectors relative to b
    real(double), dimension(3)             :: b_norm,vxB,r_gyro,rp,pc
    real(double)                           :: eb ,ptch
    integer                                :: ienergy, ipitch ,ii
    real(double)                           :: vabs, phi, sinus
    real(double)                           :: one_over_omega,b_abs
    real(double), dimension(3)             :: randomu3
    real(double), dimension(4)             :: randomu4
    integer, dimension(1) :: minpos  !! dummy array to determine minloc

    call randu(randomu3)
    ri(1)=beam_grid%xx(ac(1))+ beam_grid%dr(1)*randomu3(1)
    ri(2)=beam_grid%yy(ac(2))+ beam_grid%dr(2)*randomu3(2)
    ri(3)=beam_grid%zz(ac(3))+ beam_grid%dr(3)*randomu3(3)

    !! Assumes that magnetic field does not change appreciably over gyroradius
    a=cell(ac(1),ac(2),ac(3))%plasma%a_norm(:)
    b=cell(ac(1),ac(2),ac(3))%plasma%b_norm(:)
    c=cell(ac(1),ac(2),ac(3))%plasma%c_norm(:)
    b_abs=cell(ac(1),ac(2),ac(3))%plasma%b_abs
    b_norm=cell(ac(1),ac(2),ac(3))%plasma%b_norm(:)
    one_over_omega=inputs%ab*mass_u/(b_abs*e0)

    !! Use rejection method to determine velocity vector
    vi=0.d0
    rejection_loop: do ii=1,10000
       call randu(randomu4)
       !! Pick a random energy, pitch, and gyro angle
       eb   = fbm%emin + fbm%eran * randomu4(1)
       ptch = fbm%pmin + fbm%pran * randomu4(2)
       phi  = 2.d0*pi*randomu4(3)

       !! Calculate gyroradius
       vabs  = sqrt(eb/(v_to_E*inputs%ab))
       sinus = sqrt(1.d0-ptch**2)
       vi(:) = vabs * (sinus*cos(phi)*a + ptch*b + sinus*sin(phi)*c)
       vxB(1)= (vi(2) * b_norm(3) - vi(3) * b_norm(2))
       vxB(2)= (vi(3) * b_norm(1) - vi(1) * b_norm(3))
       vxB(3)= (vi(1) * b_norm(2) - vi(2) * b_norm(1))
       r_gyro(:)=vxB(:)*one_over_omega

       !! Move a gyrorbit away and sample distribution there
       rp(:)=ri(:)+r_gyro(:)
       !! find new cell
       minpos=minloc(abs(rp(1)-beam_grid%xxc))
       pc(1)=minpos(1)
       minpos=minloc(abs(rp(2)-beam_grid%yyc))
       pc(2)=minpos(1)
       minpos=minloc(abs(rp(3)-beam_grid%zzc))
       pc(3)=minpos(1)
       if(allocated(cell(pc(1),pc(2),pc(3))%fbm)) then
         !! take point in FBM distribution closest to eb, ptch.
         minpos=minloc(abs(eb   - fbm%energy))
         ienergy= minpos(1)
         minpos=minloc(abs(ptch - fbm%pitch ))
         ipitch = minpos(1)
         if((cell(pc(1),pc(2),pc(3))%fbm(ienergy,ipitch)).gt.randomu4(4))then
            return
         endif
       endif
       vi=0.d0
    enddo rejection_loop

    print*, 'rejection method found no solution!'
    !print*, cell(ac(1),ac(2),ac(3))%fbm
  end subroutine mc_fastion

  subroutine mc_halo(ac,vhalo)
    integer, dimension(3) , intent(in)          :: ac    !! index of actual cell
    real(double),    dimension(3) , intent(out) :: vhalo !! velocity [cm/s]
    real(double),    dimension(3)               :: randomn
    call randn(randomn)
    vhalo(:)=   cell(ac(1),ac(2),ac(3))%plasma%vrot(:) &
         + sqrt(cell(ac(1),ac(2),ac(3))%plasma%ti      &
         * 0.5/(v_to_E*inputs%ai)) &
         * randomn(:) !![cm/s]
  end subroutine mc_halo

  subroutine mc_nbi(vnbi,efrac,rnbi)
    !!-- mc_nbi computes monte carlo velocity and initial start position of
    !!-- NBI neutrals on the FIDASIM grid
    !! In this routine xyz refers to the beam sightline coordinate system where
    !! xyz_src is the origin uvw will then refer to beam grid coordinates(xyz)
    integer, intent(in)                               :: efrac !! energy fraction
    real(double), dimension(3), intent(out)           :: vnbi  !! velocity [cm/s]
    real(double), dimension(3), intent(out), optional :: rnbi  !! postition
    real(double), dimension(3)                        :: r_exit
    real(double), dimension(3)   :: uvw_src    !! Start position on ion source
    real(double), dimension(3)   :: xyz_src    !! Start position on ion source
    real(double), dimension(3)   :: uvw_ray    !! NBI velocity in uvw coords
    real(double), dimension(3)   :: xyz_ray    !! NBI velocity in xyz coords
    real(double), dimension(2)   :: randomu    !! uniform random numbers
    real(double), dimension(2)   :: randomn    !! normal random numbers
    integer                      :: nstep
    real(double)                 :: length
    !! ------------ Random start position on ion source grid -------------- !!
    call randu(randomu)
    xyz_src(1) =  0.d0
    xyz_src(2) =  nbi%dv * 2.d0*(randomu(1)-0.5d0)
    xyz_src(3) =  nbi%dw * 2.d0*(randomu(2)-0.5d0)

    !! Create random velocity vector
    call randn(randomn)
    xyz_ray(1)= 1.d0
    xyz_ray(2)=(xyz_src(2)/nbi%focy + tan(nbi%divay(efrac)*randomn(1)))
    xyz_ray(3)=(xyz_src(3)/nbi%focz + tan(nbi%divaz(efrac)*randomn(2)))

    call xyz_to_uvw(xyz_src,uvw_src,origin_in=nbi%xyz_src,alpha_in=nbi%alpha,beta_in=nbi%beta)
    call xyz_to_uvw(xyz_ray,uvw_ray,alpha_in=nbi%alpha,beta_in=nbi%beta)
    vnbi(:)=uvw_ray(:)/sqrt(dot_product(uvw_ray,uvw_ray))

    !! ----------- Determine start postition on FIDASIM grid --------- !!
    if(present(rnbi)) then
      call grid_intersect(uvw_src,vnbi,length,rnbi,r_exit)
      if(length.le.0.0) rnbi=(/-1.0,0.0,0.0/)
    endif

    !! ---- Determine velocity of neutrals corrected by efrac ---- !!
    vnbi(:) = vnbi(:)*nbi%vinj/sqrt(real(efrac))
  end subroutine mc_nbi

  !*****************************************************************************
  !************************Primary Simulation Routines**************************
  !*****************************************************************************
  subroutine ndmc
    integer                                :: indmc     !! counter for markers
    integer                                :: type      !! full half third En
    real(double)                           :: nlaunch   !! nr. of markers
    real(double)                           :: nneutrals !! # NBI particles
    real(double), dimension(3)             :: vnbi      !! velocities(full..)
    real(double), dimension(3)             :: rnbi      !! initial position
    !!Tracking routine output
    integer                                :: jj     !! counter for track
    integer                                :: ncell  !! number of cells
    real(double), dimension(  beam_grid%ntrack) :: tcell  !! time per cell
    integer,dimension(3,beam_grid%ntrack)       :: icell  !! index of cells
    real(double), dimension(3,beam_grid%ntrack) :: pos    !! mean position in cell
    integer,dimension(3)                   :: ac     !! actual cell
    !!collisional radiative model and spectrum calculation
    real(double), dimension(nlevs)         :: states
    real(double)                           :: photons
    print*,'    # of markers: ',inputs%n_nbi
    !! ------------- calculate nr. of injected neutrals ---------------- !!
    !! # of injected neutrals = NBI power/energy_per_particle
    nneutrals=1.d6*nbi%pinj/ (1.d3*nbi%einj*e0 &
         *( nbi%species_mix(1)      &
         +  nbi%species_mix(2)/2.d0 &
            +  nbi%species_mix(3)/3.d0 ) )
    !! ------------------ loop over the markers ------------------------ !!
    nlaunch=real(inputs%n_nbi)
    !$OMP PARALLEL DO schedule(guided) private(indmc,vnbi,rnbi,tcell,icell,pos,ncell,states,ac,photons,type,jj)
    energy_fractions: do type=1,3
       !! (type = 1: full energy, =2: half energy, =3: third energy
       loop_over_markers: do indmc=1,inputs%n_nbi
          call mc_nbi(vnbi(:),type,rnbi(:))
          if(rnbi(1).eq.-1)cycle loop_over_markers
          call track3D(rnbi,vnbi,tcell,icell,pos,ncell)
          if(ncell.eq.0) cycle loop_over_markers
          !! --------- solve collisional radiative model along track ----- !!
          states=0.d0
          states(1)=nneutrals*nbi%species_mix(type)/beam_grid%dv
          loop_along_track: do jj=1,ncell
             ac=icell(:,jj)
             call colrad(ac,vnbi,tcell(jj),states,photons,type,nlaunch)
             if(photons.gt.0.d0) then
                if(inputs%calc_spec.eq.1) &
                     call spectrum(vnbi(:),ac,pos(:,jj),photons,type)
             endif
          enddo loop_along_track
       enddo loop_over_markers
    enddo energy_fractions
    !$OMP END PARALLEL DO
    if(nbi_outside.gt.0)then
       print*, 'Percent of markers outside the grid: ' &
            ,100.*nbi_outside/(3.*inputs%n_nbi)
       if(sum(result%neut_dens).eq.0)stop 'Beam does not intersect the grid!'
    endif
  end subroutine ndmc

  subroutine bremsstrahlung
    integer        :: i,j,k,ichan,cnt   !! indices of cells
    real(double)   :: ne,zeff,te,gaunt,ccnt
    real(double), dimension(:)  , allocatable :: lambda_arr,brems
    !! ------------------------ calculate wavelength array ------------------ !!
    print*,'calculate the bremsstrahung!'
    allocate(lambda_arr(chords%nlambda),brems(chords%nlambda))
    do i=1,chords%nlambda
       lambda_arr(i)=(i-0.5)*chords%dlambda+chords%lambdamin ! [A]
    enddo
    ccnt=0.0
    loop_along_z: do k = 1, beam_grid%nz
       loop_along_y: do j = 1,beam_grid%ny
          loop_along_x: do i = 1, beam_grid%nx
             cnt=0
             loop_over_channels: do ichan=1,chords%nchan
                if(chords%chan_id(ichan).ne.0) cycle loop_over_channels
                if(cell(i,j,k)%los_wght(ichan).le.0.)cycle loop_over_channels
                cnt=cnt+1
                ne=cell(i,j,k)%plasma%dene     ![cm^3]
                zeff=cell(i,j,k)%plasma%zeff   !
                te=cell(i,j,k)%plasma%te*1000. ! [eV]
                if(te .le. 0.)cycle loop_over_channels
                gaunt=5.542-(3.108-log(te/1000.))*(0.6905-0.1323/zeff)
                brems(:)=7.57d-9*gaunt*ne**2*zeff/(lambda_arr(:)*sqrt(te)) &
                     *exp(-h_planck*c0/(lambda_arr(:)*te))
                result%spectra(:,cnt,brems_type) =  &
                     result%spectra(:,cnt,brems_type)  &
                     +brems(:)*cell(i,j,k)%los_wght(ichan)*1.d-2 & !!integration
                     *chords%dlambda*(4.d0*pi)*1.d-4 !! [ph/m^2/s/bin]
             enddo loop_over_channels
             ccnt=ccnt+1
             if (inputs%interactive.eq.1)then
             	!$OMP CRITICAL
             	WRITE(*,'(f7.2,"%",a,$)') ccnt/real(beam_grid%ngrid)*100,char(13)
             	!$OMP END CRITICAL
             endif
          enddo loop_along_x
       enddo loop_along_y
    enddo loop_along_z
    deallocate(lambda_arr,brems)
  end subroutine bremsstrahlung

  subroutine dcx
    integer                                :: i,j,k   !! indices of cells
    real(double), dimension(3)             :: randomu
    integer                                :: idcx    !! counter
    real(double), dimension(3)             :: ri      !! start position
    real(double), dimension(3)             :: vhalo   !! velocity bulk plasma
    integer,dimension(3)                   :: ac      !! actual cell
    !! Determination of the CX probability
    real(double), dimension(nlevs)         :: denn    !!  neutral dens (n=1-4)
    real(double), dimension(nlevs)         :: prob    !!  Prob. for CX
    real(double), dimension(3)             :: vnbi_f,vnbi_h,vnbi_t !! Velocity of NBIneutrals
    real(double), dimension(nlevs)         :: rates   !! Rate coefficiants forCX
    !! Collisiional radiative model along track
    real(double), dimension(nlevs)         :: states  !! Density of n-states
    integer                                :: ncell
    real(double), dimension(  beam_grid%ntrack) :: tcell   !! time per cell
    integer,dimension(3,beam_grid%ntrack)       :: icell   !! index of cells
    real(double), dimension(3,beam_grid%ntrack) :: pos     !! mean position in cell
    integer                                :: jj      !! counter along track
    real(double)                           :: photons !! photon flux
    real(double), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz)::papprox !!approx.density
    real(double)                           :: papprox_tot,ccnt
    real(double), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz)::nlaunch
    papprox=0.d0
    papprox_tot=0.d0
    !! ------------- calculate papprox needed for guess of nlaunch --------!!
    do k=1,beam_grid%nz
        do j=1,beam_grid%ny
            do i=1,beam_grid%nx
                if((cell(i,j,k)%plasma%denp-cell(i,j,k)%plasma%denf).gt.0) then
                    papprox(i,j,k)=    (sum(result%neut_dens(i,j,k,:,nbif_type))  &
                         +              sum(result%neut_dens(i,j,k,:,nbih_type))  &
                         +              sum(result%neut_dens(i,j,k,:,nbit_type))) &
                         *(cell(i,j,k)%plasma%denp-cell(i,j,k)%plasma%denf)
                    if(cell(i,j,k)%rho.lt.1.1)papprox_tot=papprox_tot+papprox(i,j,k)
                endif
            enddo
        enddo
    enddo
    call get_nlaunch(inputs%n_dcx,papprox,papprox_tot,nlaunch)

    ! Loop through all of the cells
    print*,'    # of markers: ',int(sum(nlaunch))
    ccnt=0.0
    !$OMP PARALLEL DO schedule(guided) private(i,j,k,idcx,randomu,ac,vhalo,ri,photons,rates, &
    !$OMP& prob,jj,states,vnbi_f,vnbi_h,vnbi_t,tcell,icell,pos,ncell,denn)
    loop_along_z: do k = 1, beam_grid%nz
       loop_along_y: do j = 1, beam_grid%ny
          loop_along_x: do i = 1, beam_grid%nx
             !! ------------- loop over the markers ---------------------- !!
             loop_over_dcx: do idcx=1,int(nlaunch(i,j,k))
                ac=(/i,j,k/)
                !! ---------------- calculate ri,vi and track -------------!!
                call mc_halo(ac,vhalo(:))
                call randu(randomu)
                ri(1)=beam_grid%xx(ac(1))+ beam_grid%dr(1)*randomu(1)
                ri(2)=beam_grid%yy(ac(2))+ beam_grid%dr(2)*randomu(2)
                ri(3)=beam_grid%zz(ac(3))+ beam_grid%dr(3)*randomu(3)
                call track3D(ri(:), vhalo(:), tcell, icell,pos, ncell)
                if(ncell.eq.0) cycle loop_over_dcx
                !! ---------------- calculate CX probability ------------- !!
                ac=icell(:,1)
                prob=0.d0
                vnbi_f=ri(:)-nbi%xyz_pos(:)
                vnbi_f=vnbi_f/sqrt(dot_product(vnbi_f,vnbi_f))*nbi%vinj
                vnbi_h=vnbi_f/sqrt(2.d0)
                vnbi_t=vnbi_f/sqrt(3.d0)
                ! CX with full energetic NBI neutrals ------ !!
                denn(:)=result%neut_dens(ac(1),ac(2),ac(3),:,nbif_type)
                call neut_rates(denn,vhalo,vnbi_f,rates)
                prob=prob + rates
                ! CX with half energetic NBI neutrals ------ !!
                denn(:)=result%neut_dens(ac(1),ac(2),ac(3),:,nbih_type)
                call neut_rates(denn,vhalo,vnbi_h,rates)
                prob=prob + rates
                ! CX with third energetic NBI neutrals ------ !!
                denn(:)=result%neut_dens(ac(1),ac(2),ac(3),:,nbit_type)
                call neut_rates(denn,vhalo,vnbi_t,rates)
                prob=prob + rates
                if(sum(prob).le.0.)cycle loop_over_dcx
                !! --------- solve collisional radiative model along track-!!
                states=prob*(cell(i,j,k)%plasma%denp-cell(i,j,k)%plasma%denf)
                loop_along_track: do jj=1,ncell
                   ac=icell(:,jj)
                   call colrad(ac,vhalo(:),tcell(jj),states &
                        ,photons,halo_type,nlaunch(i,j,k))
                   if(photons.le.0.d0)cycle loop_over_dcx
                   if(inputs%calc_spec.eq.1)call spectrum(vhalo(:),ac(:) &
                        ,pos(:,jj),photons,halo_type)
                enddo loop_along_track
             enddo loop_over_dcx
             ccnt=ccnt+1
             if (inputs%interactive.eq.1)then
                !$OMP CRITICAL
                WRITE(6,'(f7.2,"%",a,$)') ccnt/real(beam_grid%ngrid)*100,char(13)
                !$OMP END CRITICAL
             endif
          enddo loop_along_x
       enddo loop_along_y
    enddo loop_along_z
    !$OMP END PARALLEL DO
  end subroutine dcx

  subroutine halo
    integer                                :: i,j,k !indices of cells
    integer                                :: ihalo !! counter
    real(double), dimension(3)             :: randomu
    real(double), dimension(3)             :: ri    !! start position
    real(double), dimension(3)             :: vihalo!! velocity bulk plasma ion
    integer,dimension(3)                   :: ac    !! actual cell
    !! Determination of the CX probability
    real(double), dimension(nlevs)         :: denn    !!  neutral dens (n=1-4)
    real(double), dimension(nlevs)         :: prob    !!  Prob. for CX
    real(double), dimension(nlevs)         :: rates   !! Rate coefficiants forC
    real(double), dimension(3)             :: vnhalo  !! v of halo neutral
    integer                                :: in      !! index over halo neutral
    !! Collisiional radiative model along track
    real(double), dimension(nlevs)         :: states  ! Density of n-states
    integer                                :: ncell
    real(double), dimension(  beam_grid%ntrack) :: tcell  !! time per cell
    integer,dimension(3,beam_grid%ntrack)       :: icell    !! index of cells
    real(double), dimension(3,beam_grid%ntrack) :: pos    !! mean position in cell
    integer                                :: jj       !! counter along track
    real(double)                           :: photons  !! photon flux
    real(double), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz)::papprox,nlaunch !! approx. density
    real(double)                           :: papprox_tot ,ccnt
    !! Halo iteration
    integer                                :: hh !! counters
    real(double)                           :: dcx_dens, halo_iteration_dens
    integer  :: s1type  ! halo iteration
    integer  :: s2type  ! halo iteration

    s1type=fida_type
    s2type=brems_type
    dcx_dens=sum(result%neut_dens(:,:,:,:,halo_type))
    if(dcx_dens.eq.0)stop 'the density of DCX-neutrals is too small!'

    result%neut_dens(:,:,:,:,s1type) = result%neut_dens(:,:,:,:,halo_type)
    iterations: do hh=1,20

       !! ------------- calculate papprox needed for guess of nlaunch --------!!
       papprox=0.d0
       papprox_tot=0.d0
       do k=1,beam_grid%nz
          do j=1,beam_grid%ny
             do i=1,beam_grid%nx
                papprox(i,j,k)=sum(result%neut_dens(i,j,k,:,s1type)) &
                     *(cell(i,j,k)%plasma%denp-cell(i,j,k)%plasma%denf)
                if(cell(i,j,k)%rho.lt.1.1)papprox_tot=papprox_tot+papprox(i,j,k)
             enddo
          enddo
       enddo
       call get_nlaunch(inputs%n_halo,papprox,papprox_tot,nlaunch)
       print*, '    # of markers: ' ,int(sum(nlaunch))
       ccnt=0.0
       !$OMP PARALLEL DO schedule(guided) collapse(3) private(i,j,k,ihalo,ac,vihalo,randomu,ri,tcell,icell, &
       !$OMP& pos,ncell,prob,denn,in,vnhalo,rates,states,jj,photons)
       loop_along_z: do k = 1, beam_grid%nz
          loop_along_y: do j = 1, beam_grid%ny
             loop_along_x: do i = 1, beam_grid%nx
                !! ------------- loop over the markers ---------------------- !!
                loop_over_halos: do ihalo=1,int(nlaunch(i,j,k))
                   ac=(/i,j,k/)
                   !! ---------------- calculate ri,vhalo and track ----------!!
                   call mc_halo( ac, vihalo(:))
                   call randu(randomu)
                   ri(1)=beam_grid%xx(ac(1))+ beam_grid%dr(1)*randomu(1)
                   ri(2)=beam_grid%yy(ac(2))+ beam_grid%dr(2)*randomu(2)
                   ri(3)=beam_grid%zz(ac(3))+ beam_grid%dr(3)*randomu(3)
                   call track3D(ri(:),vihalo(:),tcell,icell,pos,ncell)
                   if(ncell.eq.0)cycle loop_over_halos
                   !! ---------------- calculate CX probability --------------!!
                   ac=icell(:,1) !! new actual cell maybe due to gyro orbit!
                   prob=0.d0
                   !CX with HALO neutrals
                   denn(:)=result%neut_dens(ac(1),ac(2),ac(3),:,s1type)
                   do in=1,int(n_halo_neutrate)
                      call mc_halo( ac(:), vnhalo(:))
                      call neut_rates(denn,vihalo,vnhalo,rates)
                      prob=prob+rates/n_halo_neutrate
                   enddo
                   if(sum(prob).le.0.)cycle loop_over_halos
                   !! --------- solve collisional radiative model along track-!!
                   states=prob*(cell(i,j,k)%plasma%denp)!p-cell(i,j,k)%plasma%denf)
                   loop_along_track: do jj=1,ncell
                      ac=icell(:,jj)
                      call colrad(ac,vihalo(:),tcell(jj),states &
                           ,photons,s2type,nlaunch(i,j,k))
                      if(photons.le.0.d0)cycle loop_over_halos
                      if(inputs%calc_spec.eq.1) call spectrum(vihalo(:),ac(:),pos(:,jj),photons,halo_type)
                   enddo loop_along_track
                enddo loop_over_halos
                ccnt=ccnt+1
                if (inputs%interactive.eq.1)then
                    !$OMP CRITICAL
                    WRITE(*,'(f7.2,"%",a,$)') ccnt/real(beam_grid%ngrid)*100,char(13)
                    !$OMP END CRITICAL
                endif
             enddo loop_along_x
          enddo loop_along_y
       enddo loop_along_z
       !$OMP END PARALLEL DO
       halo_iteration_dens=sum(result%neut_dens(:,:,:,:,s2type))
       result%neut_dens(:,:,:,:,halo_type)=result%neut_dens(:,:,:,:,halo_type) &
            + result%neut_dens(:,:,:,:,s2type)
       result%neut_dens(:,:,:,:,s1type)= result%neut_dens(:,:,:,:,s2type)
       result%neut_dens(:,:,:,:,s2type)= 0.
       if(halo_iteration_dens/dcx_dens.gt.1)exit iterations
       inputs%n_halo=inputs%n_dcx*halo_iteration_dens/dcx_dens
       if(inputs%n_halo.lt.inputs%n_dcx*0.01)exit iterations
    enddo iterations
    !! set the neutral density in s1type(fida_type) and s2type (brems) to 0!
    result%neut_dens(:,:,:,:,s1type) = 0.d0
    result%neut_dens(:,:,:,:,s2type) = 0.d0
  end subroutine halo

  subroutine fida
    integer                               :: i,j,k  !! indices  x,y,z  of cells
    integer                               :: iion,det,ip
    real(double), dimension(3)            :: ri      !! start position
    real(double), dimension(3)            :: vi      !! velocity of fast ions
    integer,dimension(3)                  :: ac      !! new actual cell
    integer,dimension(3,beam_grid%ngrid)       :: pcell
    !! Determination of the CX probability
    real(double), dimension(nlevs)        :: denn    !!  neutral dens (n=1-4)
    real(double), dimension(nlevs)        :: prob    !! Prob. for CX
    real(double), dimension(3)            :: vnbi_f,vnbi_h,vnbi_t    !! Velocity of NBI neutrals
    real(double), dimension(3)            :: vnhalo  !! v of halo neutral
    integer                               :: in      !! index of neut rates
    real(double), dimension(nlevs)        :: rates   !! Rate coefficiants for CX
    !! Collisiional radiative model along track
    real(double), dimension(nlevs)        :: states  ! Density of n-states
    integer                               :: ncell
    real(double), dimension(  beam_grid%ntrack):: tcell   !! time per cell
    integer,dimension(3,beam_grid%ntrack)      :: icell   !! index of cells
    real(double), dimension(3,beam_grid%ntrack):: pos     !! mean position in cell
    integer                               :: jj      !! counter along track
    real(double)                          :: photons !! photon flux
    real(double), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz)::papprox,nlaunch !! approx. density
    real(double)                          :: vi_abs             !! (for NPA)
    real(double), dimension(3)            :: ray,ddet,hit_pos   !! ray towards NPA
    real(double)                          :: papprox_tot,maxcnt,cnt,los_tot
    integer                               :: inpa,pcnt

    !! ------------- calculate papprox needed for guess of nlaunch --------!!
    papprox=0.d0
    papprox_tot=0.d0
    maxcnt=real(beam_grid%nx)*real(beam_grid%ny)*real(beam_grid%nz)
    pcnt=1
    los_tot=0.0
    do k=1,beam_grid%nz
       do j=1,beam_grid%ny
          loop_over_x: do i=1,beam_grid%nx
             papprox(i,j,k)=(sum(result%neut_dens(i,j,k,:,nbif_type))  &
                  +          sum(result%neut_dens(i,j,k,:,nbih_type))  &
                  +          sum(result%neut_dens(i,j,k,:,nbit_type))  &
                  +          sum(result%neut_dens(i,j,k,:,halo_type))) &
                  *          cell(i,j,k)%plasma%denf

             !!This saves time for mc NPA calculation
             if(inputs%calc_npa.eq.1) then
               do ip = 1,chords%nchan
                 if(chords%chan_id(ip).eq.1) los_tot=los_tot+cell(i,j,k)%los_wght(ip)
               enddo
             else
               los_tot=1
             endif

             if(papprox(i,j,k).gt.0.and.(los_tot.gt.0)) then
               pcell(:,pcnt)=(/i,j,k/)
               pcnt=pcnt+1
             endif
             if(cell(i,j,k)%rho.lt.1.1)papprox_tot=papprox_tot+papprox(i,j,k)
             los_tot=0
          enddo loop_over_x
       enddo
    enddo
    pcnt=pcnt-1
    maxcnt=real(pcnt)
    call get_nlaunch(inputs%n_fast,papprox,papprox_tot,nlaunch)
    print*,'    # of markers: ',int(sum(nlaunch))
    cnt=0.0
    !$OMP PARALLEL DO schedule(guided) private(ip,i,j,k,iion,ac,vi,ri,det,ray,inpa,hit_pos,ddet,&
    !$OMP vi_abs,tcell,icell,pos,ncell,jj,prob,denn,rates,vnbi_f,vnbi_h,vnbi_t,in,vnhalo,states,photons)
    loop_over_cells: do ip = 1, int(pcnt)
      npa_loop: do inpa=1,int(npa%npa_loop)
        loop_over_fast_ions: do iion=1,int(nlaunch(pcell(1,ip),pcell(2,ip),pcell(3,ip)))
          i=pcell(1,ip)
          j=pcell(2,ip)
          k=pcell(3,ip)
          ac=(/i,j,k/)
          !! ---------------- calculate vi, ri and track --------- !!
          call mc_fastion(ac,ri(:), vi(:))
          if(sum(vi).eq.0)cycle loop_over_fast_ions
          !! -------- check if particle flies into NPA detector ---- !!
          if(inputs%calc_npa.eq.1)then
            call hit_npa_detector(ri(:),vi(:),det)
            if(det.eq.0) cycle loop_over_fast_ions
          endif
          call track3D(ri(:), vi(:), tcell, icell,pos, ncell)
          if(ncell.eq.0)cycle loop_over_fast_ions
          !! ---------------- calculate CX probability --------------!!
          prob=0.d0
          vnbi_f=ri(:)-nbi%xyz_pos(:)
          vnbi_f=vnbi_f/sqrt(dot_product(vnbi_f,vnbi_f))*nbi%vinj
          vnbi_h=vnbi_f/sqrt(2.d0)
          vnbi_t=vnbi_f/sqrt(3.d0)
          ! CX with full energetic NBI neutrals
          denn(:)=result%neut_dens(ac(1),ac(2),ac(3),:,nbif_type)
          call neut_rates(denn,vi,vnbi_f,rates)
          prob=prob + rates
          ! CX with half energetic NBI neutrals
          denn(:)=result%neut_dens(ac(1),ac(2),ac(3),:,nbih_type)
          call neut_rates(denn,vi,vnbi_h,rates)
          prob=prob + rates
          ! CX with third energetic NBI neutrals
          denn(:)=result%neut_dens(ac(1),ac(2),ac(3),:,nbit_type)
          call neut_rates(denn,vi,vnbi_t,rates)
          prob=prob + rates
          ! CX with HALO neutrals
          denn(:)=result%neut_dens(ac(1),ac(2),ac(3),:,halo_type)
          do in=1,int(n_halo_neutrate)
            call mc_halo( ac(:), vnhalo(:))
            call neut_rates(denn,vi,vnhalo,rates)
            prob=prob + rates/n_halo_neutrate
          enddo

          if(sum(prob).le.0.)cycle loop_over_fast_ions
          !! --------- solve collisional radiative model along track-!!
          states=prob*cell(i,j,k)%plasma%denf
          loop_along_track: do jj=1,ncell
            ac=icell(:,jj)
            call colrad(ac(:),vi(:),tcell(jj) &
                ,states,photons,fida_type,nlaunch(i,j,k),ri(:),det)
            if(inputs%calc_npa.eq.1) then
               if(cell(ac(1),ac(2),ac(3))%rho.gt.profiles%rho_max)cycle loop_over_fast_ions
            else
               if(photons.le.0.d0)cycle loop_over_fast_ions
            endif
            if(inputs%calc_spec.eq.1) call spectrum(vi(:),ac(:),pos(:,jj),photons,fida_type)
          enddo loop_along_track
        enddo loop_over_fast_ions
      enddo npa_loop
      cnt=cnt+1
      if (inputs%interactive.eq.1)then
          !$OMP CRITICAL
          WRITE(*,'(f7.2,"%",a,$)') cnt/maxcnt*100,char(13)
          !$OMP END CRITICAL
      endif
    enddo loop_over_cells
    !$OMP END PARALLEL DO
  end subroutine fida

  subroutine fida_weight_function
    use netcdf
    real(double)                   :: radius
    real(double)                   :: photons !! photon flux
    real(double), dimension(nlevs) :: fdens,hdens,tdens,halodens
    real(double), dimension(3)     :: evec,vrot,los_vec,evec_sav,vrot_sav
    real(double), dimension(3)     :: a_norm,b_norm,c_norm,b_norm_sav
    real(double)                   :: b_abs,theta,b_abs_sav
    real(double)                   :: ti,te,dene,denp,deni,length,rho,rho_sav
    real(double)                   :: ti_sav,te_sav,dene_sav,denp_sav,deni_sav
    real(double)                   :: rad,max_wght
    integer                        :: nwav,nchan
    real(double),dimension(:)    ,allocatable  :: wav_arr,central_wavel
    integer                        :: ii,i,j,k,l   !! indices wavel,vx,vy,vz
    real(double),dimension(:,:,:,:),allocatable :: wfunct
    real(double),dimension(:)    ,allocatable :: ebarr,ptcharr,phiarr,rad_arr,theta_arr
    real(double),dimension(:,:), allocatable :: jacobian,e_grid,p_grid,vpa_grid,vpe_grid
    real(double)                   :: sinus
    real(double),dimension(3)      :: vi,vi_norm
    real(double)                   :: vabs
    real(double),dimension(n_stark):: intens !!intensity vector
    real(double),dimension(n_stark):: wavel  !!wavelength vector [A)
   !! Determination of the CX probability
    real(double),dimension(3)      :: vnbi_f,vnbi_h,vnbi_t !! Velocity of NBI neutrals
    real(double),dimension(3)      :: vhalo  !! v of halo neutral
    integer                        :: in      !! index of neut rates
    real(double),dimension(nlevs)  :: rates   !! Rate coefficiants for CX
    real(double),dimension(nlevs)  :: states  ! Density of n-states
    !! COLRAD
    real(double)                   :: dt  !!time interval in cell
    !! ---- Solution of differential equation  ---- !
    integer,dimension(3)                  :: ac  !!actual cell
    real(double), dimension(3)            :: pos !! position of mean cell
    integer                               :: cc
    real(double),dimension(  beam_grid%ntrack) :: wght   !! radiation wght per cell
    real(double),dimension(  beam_grid%ntrack) :: los_wght !! los wght
    real(double),dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz,chords%nchan) :: los_weight !! los wght
    integer(long)                         :: ichan,ind
    character(120)                        :: filename
    !! length through cloud of neutrals
    real(double), dimension(3,beam_grid%ntrack):: pos_out
    real(double), dimension(3)            :: pos_edge
    integer                               :: ic,jc,kc,jj,cnt
    integer                               :: ncell  !! number of cells
    real(double), dimension(  beam_grid%ntrack):: tcell  !! time per cell
    integer,dimension(3,beam_grid%ntrack)      :: icell  !! index of cells
    real(double)                          :: wght2
    !!netCDF variables
    integer :: ncid,dimid1,dimids(4),nwav_dimid,nchan_dimid,ne_dimid,np_dimid,nphi_dimid
    integer :: wfunct_varid,e_varid,ptch_varid,rad_varid,theta_varid,wav_varid,ichan_varid,nchan_varid
    integer :: vpe_varid,vpa_varid,j_varid

    !! DEFINE wavelength array
    nwav=(inputs%wavel_end_wght-inputs%wavel_start_wght)/inputs%dwav_wght
    allocate(wav_arr(nwav+1))
    allocate(central_wavel(nwav+1))
    print*,nwav,' wavelengths to be simulated!'
    do i=1,nwav+1
       wav_arr(i)=real(i-1.)*inputs%dwav_wght+inputs%wavel_start_wght
    enddo
    central_wavel=wav_arr+0.5*inputs%dwav_wght
    wav_arr=wav_arr*10. !![A]


    !! define pitch, energy and gyro angle arrays
    !! define energy - array
    print*, 'nr of energies, pitches and gyro angles', inputs%ne_wght,inputs%np_wght,inputs%nphi_wght
    print*, 'maximal energy: ', inputs%emax_wght
    allocate(ebarr(inputs%ne_wght))
    do i=1,inputs%ne_wght
       ebarr(i)=real(i-0.5)*inputs%emax_wght/real(inputs%ne_wght)
    enddo
    !! define pitch - array
    allocate(ptcharr(inputs%np_wght))
    do i=1,inputs%np_wght
       ptcharr(i)=real(i-0.5)*2./real(inputs%np_wght)-1.
    enddo
    !! define gyro - array
    allocate(phiarr(inputs%nphi_wght))
    do i=1,inputs%nphi_wght
       phiarr(i)=real(i-0.5)*2.d0*pi/real(inputs%nphi_wght)
    enddo

    !! define energy grid
    allocate(e_grid(inputs%ne_wght,inputs%np_wght))
    do i=1,inputs%ne_wght
        e_grid(i,:) = ebarr(i)
    enddo

    !! define pitch grid
    allocate(p_grid(inputs%ne_wght,inputs%np_wght))
    do i=1,inputs%np_wght
        p_grid(:,i) = ptcharr(i)
    enddo

    !! define velocity space grid
    allocate(vpe_grid(inputs%ne_wght,inputs%np_wght)) !! V perpendicular
    allocate(vpa_grid(inputs%ne_wght,inputs%np_wght)) !! V parallel
    vpa_grid(:,:) = sqrt((((2.0d3)*e0)/(mass_u*inputs%ab))*e_grid(:,:))*p_grid(:,:)
    vpe_grid(:,:) = sqrt((((2.0d3)*e0)/(mass_u*inputs%ab))*e_grid(:,:)*(1.0-p_grid(:,:)**2.0))

    !! define jacobian to convert between E-p to velocity
    allocate(jacobian(inputs%ne_wght,inputs%np_wght))
    jacobian(:,:) = ((inputs%ab*mass_u)/(e0*1.0d3)) *vpe_grid(:,:)/sqrt(vpa_grid(:,:)**2.0 + vpe_grid(:,:)**2.0)

    nchan=0
    do i=1,chords%nchan
      if(chords%chan_id(i).eq.0) nchan=nchan+1
    enddo
    print*,'Number of Channels: ',nchan

    !! define storage arrays
    allocate(wfunct(nwav,inputs%ne_wght,inputs%np_wght,nchan))
    allocate(rad_arr(nchan))
    allocate(theta_arr(nchan))

    !!save the los-weights into an array
    !! because the structure is over-written
    do k=1,beam_grid%nz
       do j=1,beam_grid%ny
          do i=1,beam_grid%nx
             los_weight(i,j,k,:)=cell(i,j,k)%los_wght(:)
          enddo
       enddo
    enddo
    !! use only cell 111 for the calculation!
    ac=(/1,1,1/)
    denp_sav=cell(ac(1),ac(2),ac(3))%plasma%denp
    dene_sav=cell(ac(1),ac(2),ac(3))%plasma%dene
    deni_sav=cell(ac(1),ac(2),ac(3))%plasma%deni
    ti_sav=cell(ac(1),ac(2),ac(3))%plasma%ti
    te_sav=cell(ac(1),ac(2),ac(3))%plasma%te
    vrot_sav=cell(ac(1),ac(2),ac(3))%plasma%vrot
    b_abs_sav=cell(ac(1),ac(2),ac(3))%plasma%b_abs
    b_norm_sav=cell(ac(1),ac(2),ac(3))%plasma%b_norm
    evec_sav=cell(ac(1),ac(2),ac(3))%plasma%E
    rho_sav=cell(ac(1),ac(2),ac(3))%rho

    cnt=1
    loop_over_channels: do ichan=1,chords%nchan
       if(inputs%ichan_wght.gt.0) then
          if(ichan.ne.inputs%ichan_wght)cycle loop_over_channels
       endif
       if(chords%chan_id(ichan).gt.0)cycle loop_over_channels

       print*,'channel:',ichan

       radius=chords%radius(ichan)
       print*,'Radius:',radius

       !! Calcullate mean kinetic profiles...
       cc=0       ; max_wght=0.d0 ; los_wght=0.d0 ; wght=0.d0
       fdens=0.d0 ; hdens=0.d0    ; tdens=0.d0    ; halodens=0.d0
       b_abs=0.d0 ; evec=0.d0
       b_norm=0.d0
       ti=0.d0    ; te=0.d0
       dene=0.d0  ; denp=0.d0     ; deni=0.d0
       vrot=0.d0  ; pos=0.d0      ; rho=0.d0
       do k=1,beam_grid%nz
          do j=1,beam_grid%ny
             do i=1,beam_grid%nx
                if(los_weight(i,j,k,ichan).gt.0.)then
                   cc=cc+1
                   los_wght(cc)=los_weight(i,j,k,ichan)
                   !! determine mean values like the halo density along LOS
                   wght(cc)=(result%neut_dens(i,j,k,3,nbif_type)   &
                        + result%neut_dens(i,j,k,3,nbih_type) &
                        + result%neut_dens(i,j,k,3,nbit_type)  &
                        + result%neut_dens(i,j,k,3,halo_type))*los_wght(cc)
                   if (wght(cc).gt.max_wght)max_wght=wght(cc)
                   fdens=fdens &
                        +result%neut_dens(i,j,k,:,nbif_type)*los_wght(cc)
                   hdens=hdens &
                        +result%neut_dens(i,j,k,:,nbih_type)*los_wght(cc)
                   tdens=tdens &
                        +result%neut_dens(i,j,k,:,nbit_type)*los_wght(cc)
                   halodens=halodens &
                        +result%neut_dens(i,j,k,:,halo_type)*los_wght(cc)
                   b_abs    =b_abs+cell(i,j,k)%plasma%b_abs  * wght(cc)
                   b_norm(:)=b_norm(:)+cell(i,j,k)%plasma%b_norm(:) * wght(cc)
                   evec(:)=evec(:)+cell(i,j,k)%plasma%E(:)  * wght(cc)
                   ti     =ti     +cell(i,j,k)%plasma%ti    * wght(cc)
                   te     =te     +cell(i,j,k)%plasma%te    * wght(cc)
                   dene   =dene   +cell(i,j,k)%plasma%dene  * wght(cc)
                   denp   =denp   +cell(i,j,k)%plasma%denp  * wght(cc)
                   deni   =deni   +cell(i,j,k)%plasma%deni  * wght(cc)
                   rho    =rho    +cell(i,j,k)%rho          * wght(cc)
                   vrot(:)=vrot(:)+cell(i,j,k)%plasma%vrot(:)*wght(cc)
                   pos(:)=pos(:)+(/beam_grid%xxc(i),beam_grid%yyc(j),beam_grid%zzc(k)/) &
                        *wght(cc)
                endif
             enddo
          enddo
       enddo
       if(max_wght.eq.0.) then
          print*,'Skipping Channel: Neutral Density is Zero'
          cnt=cnt+1
          cycle loop_over_channels
       endif
       length=sum(los_wght(:)*wght(:))/max_wght ! (FWHM)
       print*,'intersection length of NBI and LOS: ', length
       rad=sum(wght)
       pos =pos/rad
       vnbi_f=pos(:)-nbi%xyz_pos(:)
       vnbi_f=vnbi_f/sqrt(dot_product(vnbi_f,vnbi_f))*nbi%vinj
       vnbi_h=vnbi_f/sqrt(2.d0)
       vnbi_t=vnbi_f/sqrt(3.d0)
       !! normalize quantities
       rho=rho / rad
       cell(ac(1),ac(2),ac(3))%rho=rho
       denp=denp / rad
       cell(ac(1),ac(2),ac(3))%plasma%denp=denp
       dene=dene / rad
       cell(ac(1),ac(2),ac(3))%plasma%dene=dene
       deni=deni / rad
       cell(ac(1),ac(2),ac(3))%plasma%deni=deni
       ti=ti     / rad
       cell(ac(1),ac(2),ac(3))%plasma%ti=ti
       te=te     / rad
       cell(ac(1),ac(2),ac(3))%plasma%te=te
       vrot=vrot / rad
       cell(ac(1),ac(2),ac(3))%plasma%vrot=vrot
       b_abs=b_abs    / rad
       print*, '|B|: ',real(b_abs,float), ' T'
       cell(ac(1),ac(2),ac(3))%plasma%b_abs=b_abs

       b_norm=b_norm/ rad
       cell(ac(1),ac(2),ac(3))%plasma%b_norm=b_norm
       call calc_perp_vectors(b_norm,a_norm,c_norm)

       evec=evec    / rad
       cell(ac(1),ac(2),ac(3))%plasma%E=evec
       !! set los_wght to 1 only for one channel (needed for spectrum routine)
       cell(ac(1),ac(2),ac(3))%los_wght(:)=0.
       cell(ac(1),ac(2),ac(3))%los_wght(ichan)=1.
       !! Determine the angle between the B-field and the Line of Sight
       los_vec(1)=pos(1)-chords%xyzlens(ichan,1)
       los_vec(2)=pos(2)-chords%xyzlens(ichan,2)
       los_vec(3)=pos(3)-chords%xyzlens(ichan,3)
       !! normalize los_vec and bvec and determine angle between B and LOS
       los_vec=los_vec/sqrt(dot_product(los_vec,los_vec))
       theta=180.-acos(dot_product(b_norm,los_vec))*180./pi
       print*,'Angle between B and LOS [deg]:', theta
       !! write angle and radius into the output file
       rad_arr(cnt)=radius
       theta_arr(cnt)=theta
       !! START calculation of weight functions
       print*, 'nwav: ' ,nwav
       print*,''
       !! do the main simulation  !!
       !$OMP PARALLEL DO schedule(guided) private(i,j,k,ind,vabs,sinus,vi,states,    &
       !$OMP& rates,in,vhalo,dt,photons,wavel,intens,l,ii, &
       !$OMP& tcell,icell,pos_out,ncell,pos_edge,cc,max_wght,   &
       !$OMP& los_wght,wght,jj,ic,jc,kc,wght2,length,vi_norm)
       !! LOOP over the three velocity vector components
       !! (energy,pitch,gyro angle)
       do i = 1, inputs%ne_wght !! energy loop
          vabs = sqrt(ebarr(i)/(v_to_E*inputs%ab))
          do j = 1, inputs%np_wght !! pitch loop
             sinus = sqrt(1.d0-ptcharr(j)**2)
             do k = 1, inputs%nphi_wght !! gyro angle
                !! cacluate velocity vector from energy,pitch, gyro angle
                vi_norm(:)=sinus*cos(phiarr(k))*a_norm+ptcharr(j) &
                     *b_norm+sinus*sin(phiarr(k))*c_norm
                !! calcualte possible trajectory of a fast-neutral that intersects
                !! the measurment position with this velocity
                !! The measurement postion is at maximum NBI density along LOS
                call track3D(pos,-vi_norm,tcell,icell,pos_out,ncell)
                pos_edge=pos_out(:,ncell-2)
                !! now determine the track length throught the grid!
                call track3D(pos_edge,vi_norm,tcell,icell,pos_out,ncell)

                !! Calculate averge beam density seen by the fast-ion
                cc=0 ; max_wght=0.d0 ; los_wght=0.d0 ;wght=0.d0
                loop_along_track: do jj=1,ncell
                   ic=icell(1,jj)
                   jc=icell(2,jj)
                   kc=icell(3,jj)
                   wght2=result%neut_dens(ic,jc,kc,3,nbif_type) &
                        + result%neut_dens(ic,jc,kc,3,nbih_type) &
                        + result%neut_dens(ic,jc,kc,3,nbit_type) &
                        + result%neut_dens(ic,jc,kc,3,halo_type)
                   if (wght2.gt.0)then
                      cc=cc+1
                      los_wght(cc)=    tcell(jj)
                      wght(cc)=wght2 * tcell(jj)
                      if (wght(cc).gt.max_wght)max_wght=wght(cc)
                   endif
                enddo loop_along_track
                length=sum(los_wght(:)*wght(:))/max_wght ! (FWHM)
                !! determine time by length and velocity
                !! calculate the average time until a fast-neutral is
                !! detected after its neutralization
                dt=length/vabs
                !! -------------- calculate CX probability -------!!
                ! CX with full energetic NBI neutrals
                states=0.d0
                vi(:) = vi_norm(:)*vabs
                call neut_rates(fdens,vi,vnbi_f,rates)
                states=states + rates
                ! CX with half energetic NBI neutrals
                call neut_rates(hdens,vi,vnbi_h,rates)
                states=states + rates
                ! CX with third energetic NBI neutrals
                call neut_rates(tdens,vi,vnbi_t,rates)
                states=states + rates
                ! CX with HALO neutrals
                do in=1,int(n_halo_neutrate)
                   call mc_halo( ac(:),vhalo(:))
                   call neut_rates(halodens,vi,vhalo,rates)
                   states=states + rates/n_halo_neutrate
                enddo
                call colrad(ac,vi,dt,states,photons,0,1.d0)
                if(isnan(photons))photons=0
                !! photons: [Ph*cm/s/fast-ion]-!!
                !! calcualte spectrum of this one fast-ion
                call spectrum(vi,ac,pos,1.d0,nbif_type,wavel,intens)
                stark_components: do l=1,n_stark
                   wavelength_ranges: do ii=1,nwav
                      if (wavel(l).ge.wav_arr(ii).and. &
                           wavel(l).lt.wav_arr(ii+1)) then
                           !calc weight functions w/o cross-sections:
                           !wfunct(ii,i,j,ind) = wfunct(ii,i,j,ind) &
                           !    + intens(l)/real(inputs%ne_wght)
                           !normal calculation:
                           wght2 = intens(l)*photons/real(inputs%nphi_wght)
                           if (wght2.gt.2.22d-16) then
                               wfunct(ii,i,j,cnt) = wfunct(ii,i,j,cnt) + wght2
                           endif
                      endif
                   enddo wavelength_ranges
                enddo stark_components
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
       !!wfunct:[Ph*cm/s] !!
       cnt=cnt+1
    enddo loop_over_channels
    !! Put back plasma values so it doesn't possibly poison the rest of the code
    do k=1,beam_grid%nz
       do j=1,beam_grid%ny
          do i=1,beam_grid%nx
             cell(i,j,k)%los_wght(:)=los_weight(i,j,k,:)
          enddo
       enddo
    enddo

    !! use only cell 111 for the calculation!
    ac=(/1,1,1/)
    cell(ac(1),ac(2),ac(3))%rho=rho_sav
    cell(ac(1),ac(2),ac(3))%plasma%denp=denp_sav
    cell(ac(1),ac(2),ac(3))%plasma%dene=dene_sav
    cell(ac(1),ac(2),ac(3))%plasma%deni=deni_sav
    cell(ac(1),ac(2),ac(3))%plasma%ti=ti_sav
    cell(ac(1),ac(2),ac(3))%plasma%te=te_sav
    cell(ac(1),ac(2),ac(3))%plasma%vrot=vrot_sav
    cell(ac(1),ac(2),ac(3))%plasma%b_abs=b_abs_sav
    cell(ac(1),ac(2),ac(3))%plasma%b_norm=b_norm_sav
    cell(ac(1),ac(2),ac(3))%plasma%E=evec_sav

    !! Open file for the outputs
    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_fida_weights.cdf"

    !Create netCDF file
    call check( nf90_create(filename, cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid) )

    !Define Dimensions
    call check( nf90_def_dim(ncid,"dim001",1,dimid1) )
    call check( nf90_def_dim(ncid,"nwav",nwav,nwav_dimid) )
    call check( nf90_def_dim(ncid,"nchan",nchan,nchan_dimid) )
    call check( nf90_def_dim(ncid,"ne_wght",inputs%ne_wght,ne_dimid) )
    call check( nf90_def_dim(ncid,"np_wght",inputs%np_wght,np_dimid) )
    call check( nf90_def_dim(ncid,"nphi_wght",inputs%nphi_wght,nphi_dimid) )
    dimids = (/ nwav_dimid, ne_dimid, np_dimid, nchan_dimid /)

    !Define variables
    call check( nf90_def_var(ncid,"nchan",NF90_INT,dimid1,nchan_varid) )
    call check( nf90_def_var(ncid,"ichan",NF90_INT,dimid1,ichan_varid) )
    call check( nf90_def_var(ncid,"lambda",NF90_DOUBLE,nwav_dimid,wav_varid) )
    call check( nf90_def_var(ncid,"energy",NF90_DOUBLE,ne_dimid,e_varid) )
    call check( nf90_def_var(ncid,"pitch",NF90_DOUBLE,np_dimid,ptch_varid) )
    call check( nf90_def_var(ncid,"radius",NF90_DOUBLE,nchan_dimid,rad_varid) )
    call check( nf90_def_var(ncid,"theta",NF90_DOUBLE,nchan_dimid,theta_varid) )
    call check( nf90_def_var(ncid,"wfunct",NF90_DOUBLE,dimids,wfunct_varid) )
    call check( nf90_def_var(ncid,"vpa",NF90_DOUBLE,(/ne_dimid,np_dimid/),vpa_varid) )
    call check( nf90_def_var(ncid,"vpe",NF90_DOUBLE,(/ne_dimid,np_dimid/),vpe_varid) )
    call check( nf90_def_var(ncid,"jacobian",NF90_DOUBLE,(/ne_dimid,np_dimid/),j_varid) )

    !Add unit attributes
    call check( nf90_put_att(ncid,wav_varid,"units","nm") )
    call check( nf90_put_att(ncid,rad_varid,"units","cm") )
    call check( nf90_put_att(ncid,theta_varid,"units","deg") )
    call check( nf90_put_att(ncid,e_varid,"units","keV") )
    call check( nf90_put_att(ncid,wfunct_varid,"units","(Ph*cm)/(s*dE*dP)") )
    call check( nf90_enddef(ncid) )

    !Write to file
    call check( nf90_put_var(ncid, nchan_varid, nchan) )
    call check( nf90_put_var(ncid, ichan_varid, inputs%ichan_wght) )
    call check( nf90_put_var(ncid, wav_varid, central_wavel(:nwav)) )
    call check( nf90_put_var(ncid, e_varid, ebarr) )
    call check( nf90_put_var(ncid, ptch_varid, ptcharr) )
    call check( nf90_put_var(ncid, rad_varid, rad_arr) )
    call check( nf90_put_var(ncid, theta_varid, theta_arr) )
    call check( nf90_put_var(ncid, wfunct_varid, wfunct) )
    call check( nf90_put_var(ncid, vpe_varid, vpe_grid) )
    call check( nf90_put_var(ncid, vpa_varid, vpa_grid) )
    call check( nf90_put_var(ncid, j_varid, jacobian) )

    !Close netCDF file
    call check( nf90_close(ncid) )

    print*, 'fida weight function written to: ',filename

    !!Deallocate arrays
    deallocate(ebarr)
    deallocate(ptcharr)
    deallocate(phiarr)
    deallocate(wav_arr)
    deallocate(central_wavel)
    deallocate(wfunct)
    deallocate(rad_arr)
    deallocate(theta_arr)
    deallocate(jacobian)
    deallocate(e_grid)
    deallocate(p_grid)
    deallocate(vpe_grid)
    deallocate(vpa_grid)
  end subroutine fida_weight_function

  subroutine npa_weight_function
    use netcdf
    real(double)                    :: radius,theta
    real(double)                    :: photons !! photon flux
    real(double), dimension(nlevs)  :: fdens,hdens,tdens,halodens
    real(double), dimension(3)      :: los_vec
    real(double)                    :: pcxa,one_over_omega
    integer(long)                   :: nchan,cnt,det
    integer(long)                   :: ii,jj,kk,i,j,ic,jc,kc   !!indices
    integer,dimension(1)            :: minpitch,ipitch,ienergy,ix,iy,iz
    real(float),  dimension(:,:,:,:,:), allocatable :: emissivity
    real(float),  dimension(:,:,:,:,:), allocatable :: attenuation
    real(double), dimension(:,:,:),     allocatable :: wfunct_tot
    real(double), dimension(:,:),       allocatable :: flux_tot
    real(double), dimension(:),         allocatable :: ebarr,ptcharr,rad_arr
    real(double), dimension(100)    :: rd_arr,phid_arr !!npa detector array
    real(double), dimension(3)      :: vi,vi_norm,b_norm,vxB
    real(double)                    :: xcen,ycen,rshad,rs
    real(double)                    :: vabs,denf,fbm_denf,wght,b_abs,dE,dP,ccnt

    !! Determination of the CX probability
    real(double),dimension(3)       :: vnbi_f,vnbi_h,vnbi_t !! Velocity of NBI neutrals
    real(double),dimension(3)       :: vhalo  !! v of halo neutral
    integer                         :: in      !! index of neut rates
    real(double),dimension(nlevs)   :: rates,pcx   !! Rate coefficiants for CX
    real(double),dimension(nlevs)   :: states,states_i  ! Density of n-states

    !! ---- Solution of differential equation  ---- !
    integer,dimension(3)                  :: ac  !!actual cell
    real(double), dimension(3)            :: pos,rpos,dpos,rdpos,r_gyro,mrdpos !! position of mean cell
    integer(long)                         :: ichan,ind
    character(120)                        :: filename
    real(double), dimension(3,beam_grid%ntrack):: pos_out
    integer                               :: ncell  !! number of cells
    real(double), dimension(  beam_grid%ntrack):: tcell  !! time per cell
    integer,dimension(3,beam_grid%ntrack)      :: icell  !! index of cells

    !!netCDF variables
    integer :: ncid,dimid1,dimid3(3),dimid5(5),nchan_dimid,ne_dimid,np_dimid,nx_dimid,ny_dimid,nz_dimid
    integer :: wfunct_varid,e_varid,ptch_varid,rad_varid,flux_varid,emiss_varid,nchan_varid,atten_varid

    !! define pitch, energy arrays
    !! define energy - array
    print*, 'Nr of energies and pitches: ', inputs%ne_wght,inputs%np_wght
    write(*, '(A,f14.3)') ' Maximum energy: ', inputs%emax_wght

    allocate(ebarr(inputs%ne_wght))
    do i=1,inputs%ne_wght
       ebarr(i)=real(i-0.5)*inputs%emax_wght/real(inputs%ne_wght)
    enddo
    dE=abs(ebarr(2)-ebarr(1))

    !! define pitch - array
    allocate(ptcharr(inputs%np_wght))
    do i=1,inputs%np_wght
       ptcharr(i)=real(i-0.5)*2./real(inputs%np_wght)-1.
    enddo
    dP=abs(ptcharr(2)-ptcharr(1))

    nchan=0
    do i=1,chords%nchan
      if(chords%chan_id(i).eq.1) nchan=nchan+1
    enddo
    print*,'Number of Channels: ',nchan
    print*,'----'
    !! define storage arrays
    allocate(emissivity(beam_grid%nx,beam_grid%ny,beam_grid%nz,inputs%ne_wght,nchan))
    allocate(attenuation(beam_grid%nx,beam_grid%ny,beam_grid%nz,inputs%ne_wght,nchan))
    allocate(wfunct_tot(inputs%ne_wght,inputs%np_wght,nchan))
    allocate(flux_tot(inputs%ne_wght,nchan))
    allocate(rad_arr(nchan))

    emissivity(:,:,:,:,:)=emissivity(:,:,:,:,:)*0.
    attenuation(:,:,:,:,:)=attenuation(:,:,:,:,:)*0.
    wfunct_tot(:,:,:)=wfunct_tot(:,:,:)*0.
    flux_tot(:,:)=flux_tot(:,:)*0.
    cnt=1
    loop_over_channels: do ichan=1,chords%nchan
       if(chords%chan_id(ichan).ne.1)cycle loop_over_channels

       print*,'Channel:',ichan

       radius=chords%radius(ichan)
       write(*,'(A,f10.3)') ' Radius: ',radius
       rad_arr(cnt)=radius

       do i=1,100
         rd_arr(i)=(i-.5)*chords%rd(ichan)/100.
         phid_arr(i)=i*2*pi/100.
       enddo

       ccnt=0.0
       !$OMP PARALLEL DO schedule(guided) collapse(3) private(ii,jj,kk, &
       !$OMP& ic,jc,kc,ix,iy,iz,in,det,ind,ac,pos,rpos,rdpos,dpos,r_gyro,wght, &
       !$OMP& vnbi_f,vnbi_h,vnbi_t,b_norm,theta,radius,minpitch,ipitch,ienergy,mrdpos,rshad,rs,xcen,ycen, &
       !$OMP& vabs,fdens,hdens,tdens,halodens,vi,pcx,pcxa,rates,vhalo,icell,tcell,ncell,pos_out,   &
       !$OMP& states,states_i,los_vec,vi_norm,photons,denf,one_over_omega,vxB,fbm_denf,b_abs)
       loop_along_z: do kk=1,beam_grid%nz
         loop_along_y: do jj=1,beam_grid%ny
           loop_along_x: do ii=1,beam_grid%nx
            fdens=result%neut_dens(ii,jj,kk,:,nbif_type)
            hdens=result%neut_dens(ii,jj,kk,:,nbih_type)
            tdens=result%neut_dens(ii,jj,kk,:,nbit_type)
            halodens=result%neut_dens(ii,jj,kk,:,halo_type)
            b_norm(:) = cell(ii,jj,kk)%plasma%b_norm(:)
            b_abs=cell(ii,jj,kk)%plasma%b_abs
            one_over_omega=inputs%ab*mass_u/(b_abs*e0)

            if(cell(ii,jj,kk)%los_wght(ichan).gt.0) then
             pos(:) = (/beam_grid%xxc(ii), beam_grid%yyc(jj), beam_grid%zzc(kk)/)
             call chord_coor(chords%xyzlens(ichan,:),chords%xyzlos(ichan,:),chords%xyzlens(ichan,:),pos,rpos)
             xcen = -rpos(1)-rpos(1)*(rpos(3)+chords%h(ichan))/(-rpos(3))
             ycen = -rpos(2)-rpos(2)*(rpos(3)+chords%h(ichan))/(-rpos(3))
             rshad=abs(chords%ra(ichan)*(rpos(3)+chords%h(ichan))/(-rpos(3)))

             !!Loop over detector area to find mean detectorposition
             wght=0
             loop_along_xd: do i=1,100
               loop_along_yd: do j=1,100
                 rdpos(1)=rd_arr(j)*cos(phid_arr(i))
                 rdpos(2)=rd_arr(j)*sin(phid_arr(i))
                 rdpos(3)=-chords%h(ichan)
                 rs=sqrt((rdpos(1)+xcen)**2 + (rdpos(2)+ycen)**2)
                 if(rs.gt.rshad) cycle loop_along_yd
                 mrdpos(:)=mrdpos(:)+rdpos(:)
                 wght=wght+1
               enddo loop_along_yd
             enddo loop_along_xd
             if(wght.le.0) cycle loop_along_x
             mrdpos(:)=mrdpos(:)/wght
             call inv_chord_coor(chords%xyzlens(ichan,:),chords%xyzlos(ichan,:),chords%xyzlens(ichan,:),mrdpos,dpos)

             !!Determine velocity vector
             los_vec(:) = dpos(:) - pos(:)
             radius=sqrt(dot_product(los_vec,los_vec))
             los_vec=los_vec/radius
             vi_norm(:) = los_vec

             !!Check if it hits a detector just to make sure
             call hit_npa_detector(pos,vi_norm,det)
             if (det.eq.0) cycle loop_along_x

             !!Determine path particle takes throught the grid
             call track3D(pos,vi_norm,tcell,icell,pos_out,ncell)

             !!Set velocity vectors of neutral beam species
             vnbi_f(:)=pos(:) - nbi%xyz_pos(:)
             vnbi_f=vnbi_f/sqrt(dot_product(vnbi_f,vnbi_f))*nbi%vinj
             vnbi_h=vnbi_f/sqrt(2.d0)
             vnbi_t=vnbi_f/sqrt(3.d0)

             !! Determine the angle between the B-field and the Line of Sight
             theta=180.-acos(dot_product(b_norm,-los_vec))*180./pi
             minpitch=minloc(abs(ptcharr-cos(theta*pi/180.)))

             loop_over_energy: do ic = 1, inputs%ne_wght !! energy loop
               vabs = sqrt(ebarr(ic)/(v_to_E*inputs%ab))
               vi(:) = vi_norm(:)*vabs

               !!Correct for gyro orbit
               vxB(1)= (vi(2) * b_norm(3) - vi(3) * b_norm(2))
               vxB(2)= (vi(3) * b_norm(1) - vi(1) * b_norm(3))
               vxB(3)= (vi(1) * b_norm(2) - vi(2) * b_norm(1))
               r_gyro(:)=pos(:)+vxB(:)*one_over_omega
               ix=minloc(abs(r_gyro(1)-beam_grid%xxc))
               iy=minloc(abs(r_gyro(2)-beam_grid%yyc))
               iz=minloc(abs(r_gyro(3)-beam_grid%zzc))
               fbm_denf=0
               denf=0.
               if (allocated(cell(ix(1),iy(1),iz(1))%fbm)) then
                 ienergy=minloc(abs(fbm%energy-ebarr(ic)))
                 ipitch=minloc(abs(fbm%pitch-cos(theta*pi/180.)))
                 fbm_denf=cell(ix(1),iy(1),iz(1))%fbm(ienergy(1),ipitch(1))*cell(ix(1),iy(1),iz(1))%fbm_norm(1)
                 denf=cell(ix(1),iy(1),iz(1))%plasma%denf
               endif
               if (isnan(fbm_denf)) cycle loop_over_energy
               !! -------------- calculate CX probability -------!!
               ! CX with full energetic NBI neutrals
               pcx=0.d0
               call neut_rates(fdens,vi,vnbi_f,rates)
               pcx=pcx + rates
               ! CX with half energetic NBI neutrals
               call neut_rates(hdens,vi,vnbi_h,rates)
               pcx=pcx + rates
               ! CX with third energetic NBI neutrals
               call neut_rates(tdens,vi,vnbi_t,rates)
               pcx=pcx + rates

               ! CX with HALO neutrals
               do in=1,int(n_halo_neutrate)
                 call mc_halo((/ii, jj, kk/),vhalo(:))
                 call neut_rates(halodens,vi,vhalo,rates)
                 pcx=pcx + rates/n_halo_neutrate
               enddo
               if(sum(pcx).le.0) cycle loop_over_energy

               !!Calculate attenuation
               !!by looping along rest of track
               states = pcx*10000000000000. !!needs to be large aribitrary number so colrad works
               states_i=states
               do kc=1,ncell
                 ac=icell(:,kc)
                 call colrad(ac(:),vi(:),tcell(kc)/vabs,states,photons,0,1.d0)
                 if (cell(ac(1),ac(2),ac(3))%rho.gt.profiles%rho_max) then
                     exit
                 endif
               enddo
               !$OMP CRITICAL(npa_wght)
               pcxa=sum(states)/sum(states_i) !!This is probability of a particle not attenuating into plasma
               attenuation(ii,jj,kk,ic,cnt) = pcxa

               wfunct_tot(ic,minpitch(1),cnt) = wfunct_tot(ic,minpitch(1),cnt) + &
                 sum(pcx)*pcxa*cell(ii,jj,kk)%los_wght(ichan)*beam_grid%dv

               if (allocated(cell(ix(1),iy(1),iz(1))%fbm)) then
                 flux_tot(ic,cnt) = flux_tot(ic,cnt) + &
                   2*beam_grid%dv*denf*fbm_denf*sum(pcx)*pcxa*cell(ii,jj,kk)%los_wght(ichan)
                 emissivity(ii,jj,kk,ic,cnt)=emissivity(ii,jj,kk,ic,cnt)+ &
                   2*denf*fbm_denf*sum(pcx)*pcxa*cell(ii,jj,kk)%los_wght(ichan)
                   !Factor of 2 above is to convert fbm to ions/(cm^3 dE (domega/4pi))
               endif
               !$OMP END CRITICAL(npa_wght)
             enddo loop_over_energy
            endif
            ccnt=ccnt+1
            if (inputs%interactive.eq.1)then
                !$OMP CRITICAL
                WRITE(*,'(f7.2,"%",a,$)') ccnt/real(beam_grid%ngrid)*100,char(13)
                !$OMP END CRITICAL
            endif
           enddo loop_along_x
         enddo loop_along_y
       enddo loop_along_z
       !$OMP END PARALLEL DO
      write(*,'(A,ES14.5)'),' Flux:   ',sum(flux_tot(:,cnt))*dE
      write(*,'(A,ES14.5)'),' Weight: ',sum(wfunct_tot(:,:,cnt))*dE*dP
      write(*,*) ''
      cnt=cnt+1
    enddo loop_over_channels

    !! Open file for the outputs
    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_npa_weights.cdf"

    !Create netCDF file
    call check( nf90_create(filename, cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncid) )

    !Define Dimensions
    call check( nf90_def_dim(ncid,"dim001",1,dimid1) )
    call check( nf90_def_dim(ncid,"nchan",nchan,nchan_dimid) )
    call check( nf90_def_dim(ncid,"ne_wght",inputs%ne_wght,ne_dimid) )
    call check( nf90_def_dim(ncid,"np_wght",inputs%np_wght,np_dimid) )
    call check( nf90_def_dim(ncid,"Nx",beam_grid%nx,nx_dimid) )
    call check( nf90_def_dim(ncid,"Ny",beam_grid%ny,ny_dimid) )
    call check( nf90_def_dim(ncid,"Nz",beam_grid%nz,nz_dimid) )
    dimid3 = (/ ne_dimid, np_dimid, nchan_dimid /)
    dimid5 = (/ nx_dimid, ny_dimid, nz_dimid, ne_dimid, nchan_dimid /)

    !Define variables
    call check( nf90_def_var(ncid,"nchan",NF90_INT,dimid1,nchan_varid) )
    call check( nf90_def_var(ncid,"energy",NF90_DOUBLE,ne_dimid,e_varid) )
    call check( nf90_def_var(ncid,"pitch",NF90_DOUBLE,np_dimid,ptch_varid) )
    call check( nf90_def_var(ncid,"radius",NF90_DOUBLE,nchan_dimid,rad_varid) )
    call check( nf90_def_var(ncid,"wfunct",NF90_DOUBLE,dimid3,wfunct_varid) )
    call check( nf90_def_var(ncid,"flux",NF90_DOUBLE,(/ ne_dimid,nchan_dimid /),flux_varid) )
    call check( nf90_def_var(ncid,"emissivity",NF90_FLOAT,dimid5,emiss_varid) )
    call check( nf90_def_var(ncid,"attenuation",NF90_FLOAT,dimid5,atten_varid) )

    !Add unit attributes
    call check( nf90_put_att(ncid,rad_varid,"units","cm") )
    call check( nf90_put_att(ncid,e_varid,"units","keV") )
    call check( nf90_put_att(ncid,flux_varid,"units","Neutrals/(s*dE)") )
    call check( nf90_put_att(ncid,emiss_varid,"units","Neutrals/(s*dE*dV)") )
    call check( nf90_enddef(ncid) )

    !Write to file
    call check( nf90_put_var(ncid, nchan_varid, nchan) )
    call check( nf90_put_var(ncid, e_varid, ebarr) )
    call check( nf90_put_var(ncid, ptch_varid, ptcharr) )
    call check( nf90_put_var(ncid, rad_varid, rad_arr) )
    call check( nf90_put_var(ncid, wfunct_varid, wfunct_tot) )
    call check( nf90_put_var(ncid, flux_varid, flux_tot) )
    call check( nf90_put_var(ncid, emiss_varid, emissivity) )
    call check( nf90_put_var(ncid, atten_varid, attenuation) )

    !Close netCDF file
    call check( nf90_close(ncid) )

    print*, 'npa weight function written to: ',filename

    !!Deallocate arrays
    deallocate(ebarr)
    deallocate(ptcharr)
    deallocate(wfunct_tot)
    deallocate(emissivity)
    deallocate(attenuation)
    deallocate(flux_tot)
    deallocate(rad_arr)
  end subroutine npa_weight_function
end module application

!*****************************************************************************
!*******************************Main Program**********************************
!*****************************************************************************
program fidasim
  use application
  implicit none
  integer, dimension(8)              :: time_arr,time_start,time_end !Time array
  integer                            :: i,j,k,seed
  integer                            :: hour,minu,sec
  real(double)                       :: random_init
  !! measure time
  call date_and_time (values=time_start)
  !! get filename of input
  call getarg(1,namelist_file)
  !! ----------------------------------------------------------
  !! ------ INITIALIZE THE RANDOM NUMBER GENERATOR  -----------
  !! ----------------------------------------------------------
  seed = 4 !! has to be negative
  ran_am=nearest(1.0,-1.0)/ran_IM
  ran_iy=ior(ieor(888889999,seed),1)
  ran_ix=ieor(777755555,seed)
  random_init=ran()

  !! ----------------------------------------------------------
  !! ------- READ INPUTS, PROFILES, LOS, TABLES, & FBM --------
  !! ----------------------------------------------------------
  call read_inputs
  call read_beam_grid
  call read_inter_grid
  call read_beam
  call read_chords
  call read_equilibrium
  call read_profiles
  call read_tables

  !! ----------------------------------------------------------
  !! --------------- ALLOCATE THE RESULT ARRAYS ---------------
  !! ----------------------------------------------------------
  !! neutral density array!
  allocate(result%neut_dens(beam_grid%nx,beam_grid%ny,beam_grid%nz,nlevs,ntypes))
  result%neut_dens(:,:,:,:,:)=0.d0
  !! birth profile
  if(inputs%calc_birth.eq.1)then
     allocate(result%birth_dens(beam_grid%nx,beam_grid%ny,beam_grid%nz,3,npitch_birth))
     result%birth_dens(:,:,:,:,:)=0.d0
  endif
  !! allocate the spectra array
  if(inputs%calc_spec.eq.1)then
     j=0
     do i=1,chords%nchan
       if(chords%chan_id(i).eq.0) j=j+1
     enddo
     allocate(result%spectra(chords%nlambda,j,ntypes))
     result%spectra(:,:,:)=0.d0
  endif
  if(inputs%calc_npa.eq.1)then
     inputs%n_npa=1000000
     allocate(npa%v(inputs%n_npa,3,npa%nchan))
     allocate(npa%ipos(inputs%n_npa,3,npa%nchan))
     allocate(npa%fpos(inputs%n_npa,3,npa%nchan))
     allocate(npa%wght(inputs%n_npa,npa%nchan))
     allocate(npa%E(inputs%n_npa,npa%nchan))
     allocate(npa%counter(npa%nchan))
     npa%counter(:)=0
     npa%v(:,:,:)=0
     npa%ipos(:,:,:)=0
     npa%fpos(:,:,:)=0
     npa%wght(:,:)=0
     npa%E(:,:)=0
  endif



  !! -----------------------------------------------------------------------
  !! --------------- CALCULATE/LOAD the BEAM and HALO DENSITY---------------
  !! -----------------------------------------------------------------------
  if(inputs%load_neutrals.eq.1) then
     call read_neutrals()
  else
     !! ----------- ndmc (neutral density monte carlo ---------------- !!
     call date_and_time (values=time_arr)
     write(*,"(A,I2,A,I2.2,A,I2.2)") 'ndmc:   ' ,time_arr(5), ':' &
          , time_arr(6), ':',time_arr(7)
     call ndmc
     if(inputs%calc_birth.eq.1)then
        call write_birth_profile()
     endif
     !! calculate level of bremsstrahlung
     if(inputs%calc_spec.eq.1) then
        if(inputs%calc_brems.eq.1) then
           call bremsstrahlung
        else
           call read_bremsstrahlung
        end if
     end if
     !! do the HALO calcualtion only if enough markers are defined!
     if(inputs%n_halo.gt.10)then
        !! -------------------------- DCX (Direct charge exchange) ---------- !!
        call date_and_time (values=time_arr)
        write(*,"(A,I2,A,I2.2,A,I2.2)") 'dcx:    ' ,time_arr(5), ':' &
             , time_arr(6), ':',time_arr(7)
        call dcx
        !! ------------------------- HALO ----------------------------------- !!
        call date_and_time (values=time_arr)
        write(*,"(A,I2,A,I2.2,A,I2.2)") 'halo:   ' ,time_arr(5), ':' &
             , time_arr(6), ':',time_arr(7)
        call halo
     endif
     !! ---------- write output        ----------------------------------- !!
     call write_neutrals()
  endif

  !! -----------------------------------------------------------------------
  !!-----------------------LOAD FAST ION DISTRIBUTION-----------------------
  !! -----------------------------------------------------------------------
  if(inputs%load_fbm.eq.1) call read_fbm

  !! -----------------------------------------------------------------------
  !! --------------- CALCULATE the FIDA RADIATION/ NPA FLUX ----------------
  !! -----------------------------------------------------------------------
  if(inputs%calc_npa.eq.1 .or. inputs%calc_spec.eq.1 .and. inputs%n_fast.gt.10)then
     call date_and_time (values=time_arr)
     write(*,"(A,I2,A,I2.2,A,I2.2)") 'D-alpha main: ' ,time_arr(5), ':' &
          , time_arr(6), ':',time_arr(7)
     if(inputs%calc_npa.eq.1) then
       print*,'Starting NPA Simulation'
     else
       print*,'Starting FIDA Simulation'
     endif

     if(inputs%calc_npa.eq.1) then
       allocate(npa%energy(fbm%nenergy))
       allocate(npa%flux(fbm%nenergy,npa%nchan))
       npa%energy(:)=fbm%energy
     endif
     call fida

  endif

  !! -------------------------------------------------------------------
  !! ------------------- Write Spectrum/NPA to file --------------------
  !! -------------------------------------------------------------------
  if(inputs%calc_npa.eq.1.and.inputs%n_fast.gt.10) call write_npa()
  if(inputs%calc_spec.eq.1.and.inputs%n_fast.gt.10) call write_spectra()
  if(inputs%calc_spec.eq.1)deallocate(result%spectra)

  !! -------------------------------------------------------------------
  !! ----------- Calculation of weight functions -----------------------
  !! -------------------------------------------------------------------
  if(inputs%calc_fida_wght.eq.1) then
     colrad_threshold=0. !! to speed up simulation!
     call date_and_time (values=time_arr)
     write(*,"(A,I2,A,I2.2,A,I2.2)") 'fida weight function:    '  &
          ,time_arr(5), ':', time_arr(6), ':',time_arr(7)
     call fida_weight_function()
  endif

  if(inputs%calc_npa_wght.eq.1) then
     call date_and_time (values=time_arr)
     write(*,"(A,I2,A,I2.2,A,I2.2)") 'npa weight function:    '  &
          ,time_arr(5), ':', time_arr(6), ':',time_arr(7)
     call npa_weight_function()
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
  endif
  !! grid and cell structure
  do k = 1, beam_grid%nz
     do j = 1, beam_grid%ny
        do i = 1, beam_grid%nx
           if(allocated(cell(i,j,k)%fbm))deallocate(cell(i,j,k)%fbm)
           deallocate(cell(i,j,k)%los_wght)
        enddo
     enddo
  enddo
  deallocate(chords%xyzlos)
  deallocate(chords%xyzlens)
  !! result arrays
  deallocate(result%neut_dens)
  if(inputs%calc_birth.eq.1)deallocate(result%birth_dens)


  deallocate(cell)
  deallocate(beam_grid%xx)
  deallocate(beam_grid%yy)
  deallocate(beam_grid%zz)
  !!deallocate npa arrays
  if(inputs%calc_npa .eq. 1)then
     deallocate(npa%v)
     deallocate(npa%ipos)
     deallocate(npa%fpos)
     deallocate(npa%wght)
     deallocate(npa%E)
     deallocate(npa%energy)
     deallocate(npa%flux)
     deallocate(npa%counter)
  endif
end program fidasim
