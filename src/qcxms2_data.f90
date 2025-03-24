module qcxms2_data
   use xtb_mctc_accuracy, only: wp
   implicit none
   public

   type :: timer
      integer :: times = 0
      integer(wp) :: rate
      integer(wp), allocatable :: t(:, :)
      character(len=80), allocatable :: names(:)
   contains
      procedure :: init => init_timer
      procedure :: clear => clear_timer
      procedure :: start => start_timer
      procedure :: stop => stop_timer
   end type timer
   ! this type contains all runtype settings
   type :: runtypedata
      !--- GENERAL data
      integer :: chrg = 1          !molecular charge
      integer :: nfrags = 1 ! number of generated fragments including input structure
      integer :: gsmnodes = 5 ! maybe with input adjustable
      integer :: nsamples = 100000 ! 10**5 is bestnumber of runs for Monte Carlo simulation for Intensity calculation
      real(wp) :: pthr = 1.0_wp ! intensities below 1 percent are discarded for subsequent fragmentations
      integer :: temp = 298     !  thermal temperature ! make dependent from energy???
      integer :: tsnds = 9
      integer :: threads = 4 ! is variable number of threads for a single calculation
      integer :: cores = 4 ! overall number of cores specified by input

      !--- various names and flags
      character(len=80) :: infile !inputfilename in xyz format
      ! level of theory
      character(len=80) :: geolevel = 'gfn2' ! gfn2 and gfn1 possible
      character(len=80) :: iplevel = 'gfn2'  ! gfn1, gfn2, ptbrpbe, r2scan3c, and wb97x3c possible
      character(len=80) :: ip2level = 'gfn2'  ! higher level for IP calculation
      character(len=80) :: tslevel = 'gfn2'  ! gfn1, gfn2, ptbrpbe, r2scan3c, and wb97x3c possible
      character(len=80) :: pathlevel = 'none'  ! gxtb possible

      ! external programs
      character(len=80) :: qmprog = 'orca' ! currently only orca possible
      character(len=3)  :: tsfinder = 'neb' !xtb (not recommended) and neb
      character(len=80) :: topocheck = 'molbar' ! inchi or molbar for topo check after optimization
      logical :: tsgeodesic = .true. ! use geodesic interpolation as guess for transition state search
      logical :: nebnormal = .false. ! use normal settings for neb search
      !some global parameters
      character(len=1024) :: path !TODO FIXME can be removed
      character(len=1024) :: index !TODO FIXME can be removed
      character(len=1024) :: startdir
      ! global run data
      ! character(len=80) :: fraglist(10000) ! list of all fragment indices just ridiculously high
      logical :: restartrun = .false. ! restart run
      ! seems not to be that ridiculous maybe has to be even higher

      ! --- various parameter
      character(len=10) :: mode = 'ei' ! EI or CID, currently only ei mode implemented
      !CID thermal, auto, forced
      character(len=80) ::  edist = 'poisson' !'gaussian' or 'poisson'
      real(wp) :: eimp0 = 70.0_wp ! maximum of energy distribution
      real(wp) :: eimpw = 0.1_wp ! width of distribution
      real(wp) :: ieeatm = 0.8_wp ! most important for IEE distribution
      real(wp) :: fixe = 0.0_wp ! for testing purposes only one energy no distribution
      real(wp) :: tf = 50 !  50 mukroseconds usually time in MS, p 42 Gross
      real(wp) :: ehomo
      real(wp) :: hmass = 1.0801_wp ! change this to higher masses to reduce artificical H-abstractions
      real(wp) :: scaleeinthdiss = 0.5_wp ! this decreases the internal energy only for -H or -H2 abstractions
      real(wp) :: scaleker = 1.0_wp ! for testing, scale kinetic energy release
      real(wp) :: erelaxtime = 5.0e-12_wp ! relaxation time for IEE distribution for H-diss scaling
      real(wp) :: sumreacscale = 1.0_wp ! for testing, scale kinetic energy release
      logical :: eatomscale = .true. ! scale IEE after fragmentation according to number of atoms of fragments
      logical :: chargeatomscale = .false. ! larger fragment gets charge
      integer :: exstates = 0 ! number of excited states to include via TDDFT
      !--- general logical data

      !CREST MSREACT settings
      logical :: msnoiso = .false. ! no rearrangments are simulated
      logical :: msfulliso = .true. ! rearrengments also for subsequent structures
      logical :: msiso = .false. ! only rearrangements are simulated
      integer :: msnbonds = 3
      integer :: msnshifts = 200
      integer :: msnshifts2 = 0
      logical :: msinchi = .false. ! sort out according to inchi code (needs obabel in path)
      logical :: msmolbar = .false. ! sort out according to molbar (needs molbar in path)
      logical :: env%mskeepdir = .false. ! keep directories for MSREACT

      real(wp) :: msfragdist = 0.0_wp ! distance of fragments in msreact
      logical :: bhess = .true. ! add thermocontribution to barrier   ! is automatically included in de, barrier, sumreac
      logical :: notemp = .true. ! only take ZPVE no thermal effects G(RRHO)
      logical :: hotip = .false. ! use hot IP for fragmentation#
      logical :: fermi = .false. ! apply fermi smearing in SCF calculations

      logical :: nots = .false. ! for quickmode take only de instead of Ea
      logical :: esim = .false. ! simulate only mcsimu energy
      logical :: plot = .false. ! plot some stuff for analyzing purposes
      logical :: oplot = .false. ! only plotting tool
      logical :: int_masses = .false. ! only plot masses as integers
      logical :: noiso = .false. ! do not plot isotope pattern
      logical :: calcKER = .false. ! compute kinetic energy release
      ! special modes
      logical :: cneintscale = .false. ! scale internal energy according to CN of active atoms
      logical :: dxtb = .false. ! use special dxtb parameter requires dxtb_param.txt in starting DIR

      logical :: sortoutcascade = .false. ! sort out structures with more than one maximum in reaction path
      logical :: ignoreip = .false. ! ignore IPs, both fragments get full intensity
      logical :: reoptts = .true. ! reoptimize TS after GSM search
      logical :: eyring = .true. ! use eyr for rate constant calculation
      logical :: eyzpve = .false. ! use eyr but only with ZPVE for rate constant calculation
      real(wp) :: mthr = 0.0_wp ! m/z threshold for plotting
      real(wp) :: intthr = 0.0_wp ! peak intensity threshold for matching score and plotting
      integer :: sthr = 150 ! RRHO cutoff in thermo contribution
      integer :: ithr = 100 ! imaginary mode cutoff in thermo contribution, very critical
      ! number of subsequent fragmentations, at max 3
      integer :: nfrag = 3
      ! special stuff
      logical :: noirc = .false. ! if true set all IRC to 100 in esim mode
      integer :: printlevel = 1 ! 1: normal, 2: verbose, 3: debug mode
      logical :: logo = .false. ! print logo
      logical :: removedirs = .true. ! removedirs after calculation
      ! for timings
      real(wp) :: tcrest, tip, tpath, tts, thess, tbhess
      ! for testing
      real(wp) :: tfscale = 0.0_wp ! 0 means no scaling
      integer :: nthermosteps = 200 !number of increments to compute thermal corrections for IEE distribution
      logical :: noplotisos = .false. ! do not plot isomers
      logical :: ircrun = .false. ! do only an IRC analysis
      logical :: picktsrun = .false. ! do only an IRC analysis

      ! CID DATA
      ! TODO
      logical :: prot = .false.   ! for positve ESI Mode
      logical :: deprot = .false.    ! for negative ESI Mode
      integer  :: cid_mode = 1 ! 1: auto 2: temprun (no collisions) 3: only collisions
      real(wp) :: cid_elab = 40 ! collision energy in laboratory fram in eV
      real(wp) :: cid_esi = 0.0_wp ! ionization energy in eV by default 0.0,
      real(wp) :: cid_esiatom = 0.0_wp ! ionization energy per atom in eV by default 0.0,
      real(wp) :: cid_esiw = 0.2_wp ! width of ionization energy distribution in eV
      real(wp) :: cid_collw = 0.5_wp ! width of collision energy distribution in eV
      integer  :: cid_maxcoll = 10 ! maximum number of collisions
      real(wp) :: cid_lchamb = 0.25 ! chamber length in m
      real(wp) :: cid_pgas = 0.132 !Pa 0.000132 ! 1mtor   !0.5 ! gas pressure in Pa ! 1 torr QCxMS1
      real(wp) :: cid_TGas = 300 ! gas temperature in K
      real(wp) :: cid_mgas = 39.94800_wp ! argon ! mass of gas in u !He, Ne, Ar, Kr, Xe, N2 available TODO include them
      real(wp) :: cid_rgas = 3.55266638_wp !  bohr vdw-radius of Ar ! TODO include them
      real(wp) :: cid_scool = 1.0 ! scale collissional cooling
   end type runtypedata

   ! some more constants
   real(wp), parameter :: bohr = 0.52917726_wp

contains

!-----------------------------------------------------------------------------------------------------
   subroutine init_timer(self, n)
      implicit none
      class(timer) :: self
      integer, intent(in)  :: n
      integer(wp) :: dummy
      self%times = n
      call system_clock(dummy, self%rate)
      allocate (self%t(n, 3), source=0_wp)
      allocate (self%names(n))
      self%names = ''
   end subroutine init_timer
!-----------------------------------------------------------------------------------------------------
   subroutine clear_timer(self)
      implicit none
      class(timer) :: self
      deallocate (self%t)
      deallocate (self%names)
   end subroutine clear_timer
!-----------------------------------------------------------------------------------------------------
   subroutine start_timer(self, n, inp)
      implicit none
      class(timer) :: self
      integer, intent(in)  :: n
      character(len=*) :: inp
      integer(wp) :: dummy
      self%names(n) = inp
      call system_clock(self%t(n, 1), dummy)
   end subroutine start_timer
!-----------------------------------------------------------------------------------------------------
   subroutine stop_timer(self, n)
      implicit none
      class(timer) :: self
      integer, intent(in)  :: n
      integer(wp) :: dummy
      call system_clock(self%t(n, 2), dummy)
      self%t(n, 3) = self%t(n, 3) + (self%t(n, 2) - self%t(n, 1))
   end subroutine stop_timer
!================================================================================!
! Most of these useful subroutines are taken from CREST (github.com/crest-lab/crest)
! (using the GPL-3.0 license like QCxMS2)
! and are slightly modified for QCxMS2
!================================================================================!

!c  Transform a time in seconds to a sting in the format hh:mm:ss
   subroutine eval_time(time, printout)
      implicit none
      real(wp) :: time
      character(len=*) :: printout
      integer(wp) :: days, hours, minutes, seconds

      hours = floor(time/3600.0d0)
      time = time - (real(hours)*3600.0d0)
      minutes = floor(time/60.0d0)
      time = time - (real(minutes)*60.0d0)
      seconds = floor(time)

      write (printout, '(i0,''h :'',i2,''m :'',i2,''s'')') &
      &     hours, minutes, seconds
   end subroutine eval_time

   subroutine eval_sub_timer(tim)
      implicit none
      type(timer) :: tim
      character(len=80) :: ftime
      integer(wp) ::  t1, ttot1
      real(wp) :: t2, ttot2
      integer :: i, j
      j = tim%times
      if (j .lt. 1) return
      do i = 1, j
         t1 = tim%t(i, 3)
         if (t1 .ne. 0) then
            t2 = real(t1)/real(tim%rate)
            call eval_time(t2, ftime)
            if (trim(tim%names(i)) .ne. '') then
               write (*, '(a24,'' wall time : '',a20)') trim(tim%names(i)), trim(ftime)
            end if
         else
            cycle
         end if
      end do
      return
   end subroutine eval_sub_timer

   subroutine eval_timer(tim)
      implicit none
      type(timer) :: tim
      character(len=80) :: ftime
      integer(wp) ::  t1, ttot1
      real(wp) :: t2, ttot2
      integer :: j
      write (*, *)
      call smallhead('Wall Time Summary')
      call eval_sub_timer(tim)
      j = tim%times
      ttot1 = sum(tim%t(1:j, 3))
      ttot2 = real(ttot1)/real(tim%rate)
      call eval_time(ttot2, ftime)
      write (*, '(''--------------------'')')
      write (*, '(''Overall wall time  : '',a)') trim(ftime)

      call tim%clear
   end subroutine eval_timer

   subroutine propquit(tim)
      implicit none
      type(timer) :: tim
      call eval_timer(tim)
      write (*, *)
      write (*, *) 'QCxMS2 terminated normally.'
      stop
   end subroutine propquit
   subroutine smallhead(str)
      implicit none
      character(len=*) :: str
      integer :: strlen
      integer :: i, j
      character(len=:), allocatable :: str2
      strlen = len_trim(str)
      str2 = repeat('-', strlen)
      write (*, '(1x,a)') trim(str2)
      write (*, '(1x,a)') trim(str)
      write (*, '(1x,a)') trim(str2)
      return
   end subroutine smallhead

end module qcxms2_data
