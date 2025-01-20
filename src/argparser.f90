!======================================================================================================!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!c ARGUMENT PARSER FOR QCxMS2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!======================================================================================================!
module argparser
   use qcxms2_data
   use mctc_env, only: error_type, get_argument!, fatal_error
   use mctc_io, only: structure_type, read_structure, write_structure, &
     & filetype, get_filetype, to_symbol
   use cid
   implicit none
   character(len=:), allocatable :: arg(:)

contains
   subroutine parseflags(env, arg, nargs, iseed)
      integer :: nargs, i, l, j
      type(runtypedata) :: env
      character(len=:), allocatable :: arg(:)
      character(len=:), allocatable :: argument

      real(wp) :: xx(10)
      integer :: iseed(1)

      ! set defaults here to GFN2
      env%geolevel = "gfn2"
      env%iplevel = 'gfn2'
      env%ip2level = ''
      env%tslevel = 'gfn2'
      env%tsfinder = 'neb'
      env%nebnormal = .false.
      env%pathlevel = 'none'
      env%tsnds = 9
      !programsettings
      env%qmprog = 'orca' !'tm' ! or orca, orca default for now

      env%tsgeodesic = .true.

!
      env%mode = 'ei' ! compute EI spectrum default !cidtemp cidauto cidforced
      env%chrg = 1

      env%nfrag = 6
      env%fixe = 0.0_wp
      env%reoptts = .false.
      env%topocheck = 'molbar'
      env%sortoutcascade = .false.
! general runtype data
      env%nots = .false.
      env%esim = .false. ! simulate only different energies
      env%calcKER = .true.
      env%bhess = .true.
      env%hmass = 1.0801 ! change this to higher masses to reduce artificical H-abstractions
      env%hotip = .false.
      env%eimp0 = 70.0_wp
      env%ieeatm = 0.8_wp
      env%eimpw = 0.1_wp
      env%nsamples = 100000
      env%pthr = 1.0_wp
      ! crestms settings
      env%msnoiso = .false.
      env%msiso = .false.
      env%msfulliso = .true.
      env%msnbonds = 3
      env%msinchi = .false. ! requires obabel in path and can be problematic for the NCI complexes
      env%msmolbar = .true. ! requires mobar environment but is best for sorting out duplicates
      env%msnshifts = 0 ! only set higher for planar molecules
      env%msnshifts2 = 0
      env%msfragdist = 2.5_wp
      env%notemp = .true. ! only take ZPVE no thermal effects G(RRHO) ! the right way to do it ...
      ! for plotting
      env%noiso = .false.
      env%oplot = .false.
      env%mthr = 0.0_wp
      env%intthr = 0.0_wp
      env%dxtb = .false.
      env%tfscale = 0.0_wp ! no scaling of time of flight of subsequent fragmentations ! TODO FIXME THIS HAS TO BE 0.0 for no scaling!!!!
      env%sumreacscale = 1.0_wp ! no scaling of sumreac
      env%eatomscale = .true. ! scale energy of subsequent fragmentations according to number of atoms
      env%tf = 50.0_wp ! time of flight in mass spectrometer in microseconds
      iseed(1) = 42
! special settings for esim
      env%noirc = .false.
! for testing
      env%cneintscale = .false.
      env%printlevel = 1 ! 1: normal, 2: verbose, 3: debug mode
      env%ignoreip = .false.
      env%fermi = .false. ! TODO make default
      env%removedirs = .false. ! for now
      env%edist = 'poisson' ! poisson or gauss
      env%chargeatomscale = .false. ! scale IEE of subsequent fragmentations accoding to number of atoms, larger fragment gets charge
      env%exstates = 0 ! number of excited states to include via TDDFT
      env%eyring = .true. ! use eyring isntead of rrkm
      env%eyzpve = .false. ! use eyring instead of rrkm but only with ZPVE ...
      env%sthr = 150 ! RRHO cutoff in cm-1
      env%ithr = 100 ! imagniary mode RRHO cutoff in cm-1
      env%nthermosteps = 200 ! number of increments to compute thermal corrections for IEE distribution
      env%noplotisos = .false. ! plot isomers in mass spectrum
      env%ircrun = .false.
      env%picktsrun = .false.
      env%erelaxtime = 5.0e-12_wp !
      env%scaleeinthdiss = 0.5_wp

      ! CID specific settings
      env%cid_mode = 1 ! 1: auto 2: temprun (no collisions) 3: only collisions
      env%cid_elab = 40 ! collision energy in laboratory fram in eV
      env%cid_esi = 0.0_wp ! ionization energy in eV by default 0.0,
      env%cid_esiatom = 0.0_wp ! ionization energy in eV by default 0.0,
      env%cid_esiw = 0.2_wp ! width of ionization energy distribution in eV
      env%cid_collw = 0.5_wp ! width of collision energy distribution in eV
      env%cid_maxcoll = 10 ! maximum number of collisions
      env%cid_lchamb = 0.25 ! chamber length in m
      env%cid_pgas = 0.132 ! gas pressure in Pa ! 1 torr QCxMS1 0.132 Pa
      env%cid_TGas = 300 ! gas temperature in K
      env%cid_mgas = 39.94800_wp ! argon ! mass of gas in u !He, Ne, Ar, Kr, Xe, N2 available TODO include them
      env%cid_rgas = 3.55266638_wp !  bohr vdw-radius of Ar ! TODO include them
      env%cid_scool = 1.0 ! no cooling of ions

      do i = 1, nargs
         if (any((/character(7)::'-h', '-H', '--h', '--H', '-help', '--help'/) == trim(arg(i)))) then
            call printhelp()
         end if
      end do

      do i = 1, nargs
         if (any((/character(12)::'-advanced', '--advanced'/) == trim(arg(i)))) then
            call printhelp_advanced()
         end if
      end do

      do i = 1, nargs
         argument = trim(arg(i))
         if (argument .ne. '') then
         if (len(argument) .gt. 1) then
         if (argument(1:2) == '--') then
            argument = argument(2:)
         end if
         end if
         select case (argument)
            ! needed??
         case ('-ei ')
            env%mode = 'ei' !perform prediction of EI spectrum
            !todo
         case ('-cid ')
            env%mode = 'cid' ! perform prediction of CID spectrum (not yet implemented) ! TODO modify mcsimu for this, add collision energies
         case ('-chrg ')  ! set charge of molecule (not yet implemented) ! TODO
            call readl(arg(i + 1), xx, j)
            env%chrg = xx(1)

         case ('-geolevel ')
            env%geolevel = arg(i + 1) ! gfn2 and gfn1 possible and wb97m-3c possible
         case ('-path ')
            env%tsfinder = arg(i + 1) ! xtb and gsm possible
         case ('-nebnormal ')
            env%nebnormal = .true. ! use normal settings instead of loose settings for NEB search
         case ('-tslevel ')
            env%tslevel = arg(i + 1) ! gfn1, gfn2, ptbrpbe, r2scan3c, pbeh3c, and wb97x3c possible
         case ('-iplevel ')
            env%iplevel = arg(i + 1) ! gfn1, gfn2, ptbrpbe, r2scan3c, pbeh3c, wb97x3c, and ptb possible
         case ('-ip2level ')
            env%ip2level = arg(i + 1) ! gfn1, gfn2, ptbrpbe, r2scan3c, pbeh3c, wb97x3c, and ptb possible
         case ('-qmprog ')
            env%qmprog = arg(i + 1) ! orca or tm (turbomole) possible
         case ('-msnoiso ')
            env%msnoiso = .true. ! sort out fragments with same mass as input molecule in crestms
         case ('-msiso ')
            env%msiso = .true. ! sort out non-dissociated fragments
         case ('-msfulliso ')
            env%msfulliso = .true. ! rearrangments also allowed in subseqent fragmentations
         case ('-msmolbar ')
            env%msmolbar = .true. ! sort out duplicates by molbar
         case ('-msinchi ')
            env%msinchi = .true. ! sort out duplicates by inchi
         case ('-msfragdist ')
            call readl(arg(i + 1), xx, j)! seperate fragments in crestms from each other
            env%msfragdist = xx(1)
         case ('-topocheck ')
            env%topocheck = arg(i + 1) ! topochech after optimization
         case ('-nobhess ')
            env%bhess = .false. ! do not add thermocontribution to activation energy
         case ('-usetemp ')
            env%notemp = .false. ! take G_mRRHO contribution instead of only ZPVE
         case ('-plotk ')
            env%plot = .true. ! plot k(E) curves
         case ('-nots ')
            env%nots = .true. ! take de instead of ea -> no path search
         case ('-esim ')
            env%esim = .true. ! simulate only different energies !TODO not yet implemented
         case ('-oplot ')
            env%oplot = .true. ! only plot peaks give exp.dat as csv file
         case ('-noisotope ')
            env%noiso = .true. ! plot no isotope pattern
         case ('-int_masses ')
            env%int_masses = .true. ! only plot masses as integers
         case ('-dxtbparam ')
            env%dxtb = .true. ! use special dxtb parameter requires dxtb_param.txt in starting DIR
         case ('-notsgeo ')
            env%tsgeodesic = .false. ! do not use geodesic interpolation as guess for gsm
         case ('-noirc ')
            env%noirc = .true. ! set all imaginary frequencies to 100 cm-1
         case ('-edist ')
            env%edist = arg(i + 1) ! poisson or gaussian possible
         case ('-eimp0 ') ! impact excess energy (IEE) in eV (default 70 eV)
            call readl(arg(i + 1), xx, j)
            env%eimp0 = xx(1)
         case ('-eimpw ')  ! width of IEE distribution
            call readl(arg(i + 1), xx, j)
            env%eimpw = xx(1)
         case ('-ieeatm ')! energy per atom in eV for IEE
            call readl(arg(i + 1), xx, j)
            env%ieeatm = xx(1)
         case ('-tf ')  !time of flight in MS give in mikroseconds, default is 50
            call readl(arg(i + 1), xx, j)
            env%tf = xx(1)
         case ('-T ') ! number of cores
            call readl(arg(i + 1), xx, j)
            env%cores = xx(1)
         case ('-nfrag ') !number of subsequent fragmentations
            call readl(arg(i + 1), xx, j)
            env%nfrag = xx(1)
         case ('-tsnodes ')
            call readl(arg(i + 1), xx, j)
            env%tsnds = xx(1)
         case ('-pthr ') !threshold for further fragmentation
            call readl(arg(i + 1), xx, j)
            env%pthr = xx(1)
         case ('-nsamples ') !number of simulated runs in Monte Carlo
            call readl(arg(i + 1), xx, j)
            env%nsamples = xx(1)
         case ('-mthr ') ! m/z threshold for plotting
            call readl(arg(i + 1), xx, j)
            env%mthr = xx(1)
         case ('-intthr ') ! peak intensity threshold for plotting
            call readl(arg(i + 1), xx, j)
            env%intthr = xx(1)
         case ('-msnbonds ') ! bonds for constraint opt. in msreact
            call readl(arg(i + 1), xx, j)
            env%msnbonds = xx(1)
         case ('-msnshifts ') ! bonds for constraint opt. in msreact
            call readl(arg(i + 1), xx, j)
            env%msnshifts = xx(1)
         case ('-msnshifts2 ') ! bonds for constraint opt. in msreact
            call readl(arg(i + 1), xx, j)
            env%msnshifts2 = xx(1)
         case ('-atmassh ')  ! increase to scale H-R frequencies
            call readl(arg(i + 1), xx, j)
            env%hmass = xx(1)
         case ('-scaleeinthdiss ')   ! this decreases the internal energy only for -H or -H2 abstractions
            call readl(arg(i + 1), xx, j)
            env%scaleeinthdiss = xx(1)
         case ('-scaleker ')   ! This scale the kinetic energy release upon fragmentation
            call readl(arg(i + 1), xx, j)
            env%scaleker = xx(1)
         case ('-tfscale ') ! scale time of flight for subsequent fragmentations
            call readl(arg(i + 1), xx, j)
            env%tfscale = xx(1)
         case ('-noeatomscale ') ! scale energy of subsequent fragmentations according to number of atoms
            env%eatomscale = .false.
         case ('-sumreacscale ') ! scale time of flight for subsequent fragmentations
            call readl(arg(i + 1), xx, j)
            env%sumreacscale = xx(1)
         case ('-erelaxtime ') ! scale time of flight for subsequent fragmentations
            call readl(arg(i + 1), xx, j)
            env%erelaxtime = xx(1)*1.0e-12_wp
         case ('-noKER ') ! no kinetic energy release
            env%calcKER = .false.
         case ('-fixe ')  ! simulate only one energy for testing purposes
            call readl(arg(i + 1), xx, j)
            env%fixe = xx(1)
         case ('-reoptts ') !reoptimize TS after path search ''
            env%reoptts = .true.
         case ('-hotip ') !do not optimize fragments at charge 0 before computing IPs ''
            env%hotip = .true.
         case ('-sortoutcascade ') ! sort out fragments with same mass as input molecule in crestms
            env%sortoutcascade = .true.
         case ('-fermi ') !apply fermi smearing
            env%fermi = .true.
         case ('-rrkm ') ! use simplified RRKM instead of Eyring
            env%eyring = .false.
         case ('-eyzpve ') ! use eyring but with DH (without entropy) instead of DG
            env%eyring = .true.
            env%eyzpve = .true.
         case ('-cneintscale ') ! scale internal energy of subsequent fragmentations according to number of atoms
            env%cneintscale = .true.
            !SEED FOR RANDOM NO. GENERATOR
         case ('iseed ')
            call readl(arg(i + 1), xx, j)
            iseed = int(xx(1))
            ! printou
         case ('-printlevel ') ! do not delete files after calculation, for diagnostic purposes ! TODO implement this
            call readl(arg(i + 1), xx, j)
            env%printlevel = int(xx(1))
         case ('-logo ') !check in hessian calculation directory for an IRC
            env%logo = .true.
         case ('-keepdir ') ! do not delete sorted out fragment dirs for diagnostic purposes
            env%removedirs = .false.
         case ('-ignoreip ') ! ignore IPs, both fragments get full intensity
            env%ignoreip = .true.
         case ('-chargeatomscale ') ! ignore IPs, large fragment gets charge
            env%chargeatomscale = .true.
         case ('-exstates ') !number of excited states to include via TDDFT
            call readl(arg(i + 1), xx, j)
            env%exstates = xx(1)
         case ('-sthr ') !RRHO cutoff in cm-1
            call readl(arg(i + 1), xx, j)
            env%sthr = xx(1)
         case ('-ithr ') !imaginary frequency RRHO cutoff in cm-1
            call readl(arg(i + 1), xx, j)
            env%ithr = xx(1)
         case ('-nthermosteps ') !number of increments to compute thermal corrections for IEE distribution
            call readl(arg(i + 1), xx, j)
            env%nthermosteps = xx(1)
         case ('-noplotisos ') !isomers are assumed to fragmegnt further and are not plotted
            env%noplotisos = .true.
         case ('-ircrun ') !check in hessian calculation directory for an IRC
            env%ircrun = .true.
         case ('-picktsrun ') !check in hessian calculation directory for an IRC
            env%picktsrun = .true.
            !CID specific settings
            ! runtype
            ! for ESI-MS TODO
         case ('-protonate ')
            env%prot = .true. !perform search of favoured protomers
         case ('-deprotonate ')
            env%deprot = .true. !perform search of favoured deprotonated structures
         case ('-cidauto ')
            env%mode = 'cid'
            env%cid_mode = 1
         case ('-cidtemprun ')
            env%mode = 'cid'
            env%cid_mode = 2
            env%cid_maxcoll = 0
         case ('-cidforced ')
            env%mode = 'cid'
            env%cid_mode = 3
            ! various parameters
         case ('-elab ') !thermal temperature in biased hessian calculation default is 5000
            call readl(arg(i + 1), xx, j)
            env%cid_elab = xx(1)
         case ('-esi ')
            call readl(arg(i + 1), xx, j)
            env%cid_esi = xx(1)
         case ('-esiatom ')
            call readl(arg(i + 1), xx, j)
            env%cid_esiatom = xx(1)
         case ('-esiw ')
            call readl(arg(i + 1), xx, j)
            env%cid_esiw = xx(1)   ! width of ionization energy distribution in eV
         case ('-collw ')
            call readl(arg(i + 1), xx, j)
            env%cid_collw = xx(1)    ! width of collision energy distribution in eV
         case ('-maxcoll ')
            call readl(arg(i + 1), xx, j)
            env%cid_maxcoll = xx(1)
         case ('-lchamb ')
            call readl(arg(i + 1), xx, j)
            env%cid_lchamb = xx(1)
         case ('-mgas ') ! in u
            call readl(arg(i + 1), xx, j)
            env%cid_mgas = xx(1)
         case ('-rgas ') ! in u
            call readl(arg(i + 1), xx, j)
            env%cid_rgas = xx(1)
         case ('-pgas ') ! in u
            call readl(arg(i + 1), xx, j)
            env%cid_pgas = xx(1)
         case ('-cidscool ') ! in u
            call readl(arg(i + 1), xx, j)
            env%cid_scool = xx(1)
         case default
            continue
         end select
         end if
      end do
      call inputxyz(env, trim(arg(1)))
   end subroutine parseflags
   ! TODO, maybe use fixtemp for testing purposes?????
   !   case ('-temp '       ) /= 0) then  !thermal temperature in biased hessian calculation default is 5000
   !  call readl(arg(i + 1),xx,j)
   ! env%temp=xx(1)
   !end if

   subroutine inputxyz(env, arg)
      implicit none
      type(runtypedata) :: env
      character(len=*) :: arg
      logical :: ex

! in some way this is double
      inquire (file=arg, exist=ex)
      if (.not. ex) then
         error stop 'No (valid) input file! exit.'
      end if

      if (ex .and. arg(1:1) .ne. '-') then
         env%infile = trim(arg)
         write (*, '(a,x,a)') "input file given:", trim(env%infile)
      end if

   end subroutine inputxyz

   subroutine printhelp()
      implicit none
      write (*, *)
      write (*, '(1x, ''usage        :'')')
      write (*, '(1x, ''qcxms2 [input] [options]'')')
      write (*, *)
      write (*, '(1x, ''The FIRST argument has to be a coordinate file in the'')')
      write (*, '(1x, ''Xmol (*.xyz, Ang.) format.'')')
      write (*, *)
      write (*, '(/,1x,''runtype options:'')')
      write (*, '(5x,''-ei: compute a electron impact (EI) spectrum (default)'')')
      write (*, '(5x,''-esim: simulate only different IEE distributions for given reaction network from previous calculation '')')
      write (*, '(5x,''-oplot: plot only peaks from previous QCxMS2 calculation '')') ! give exp.dat as csv file
      write (*, *)
      write (*, '(/,1x,''General and technical options:'')')
      write (*, '(5x,''-chrg [int]: give charge of input molecule'')')
      write (*, '(5x,''-edist: choose distribution for IEE ("poisson" (default) or "gaussian")'')')
      write (*, '(5x,''-eimp0: give impact energy of electrons in eV (default 70)'')')
      write (*, '(5x,''-eimpw: give width of energy distribution (default 0.1)'')')
      write (*, '(5x,''-ieeatm: give energy per atom in eV (default 0.8)'')')
      write (*, '(5x,''-nots: take reaction energy instead of barrier -> no path search (quickmode for fragments)'')')
      write (*, '(5x,''-T  : select number of overall cores, (default 4) '')')
      write (*, '(5x,''-nfrag [integer]: select number of subsequent fragmentations you want to simulate (default 6) '')')
      write (*, '(5x,''-pthr [real] intensity at which fragment is further fragmented in % (default  1%)'')')
      write (*, *)
      write (*, '(/,1x,''Options for QC Calculations:'')')
      write (*, '(5x,''-geolevel [method]      : method for geometry optimization and path search'')')
      write (*, '(5x,''-tslevel  [method] : select level for computing reaction barriers'')')
      write (*, '(5x,''-iplevel  [method] : select level for computing IPs for charge assignment'')')
   write (*, '(5x,''-ip2level  [method] : select level for computing IPs for charge assignment of critical cases with close IPs'')')
      write (*, '(8x,''available methods: ("gfn2","gfn1","r2scan3c","pbeh3c","wb97x3c","pbe0","gxtb")'')') !dxtb gxtb wb97m3c
      !write(*,'(5x,''-geolevel [method]      : method for geometry optimization and path search (("gfn2","gfn1","pm6","gff","pbe","b973c","r2scan3c","pbeh3c","wb97x3c","pbe0")'')') !dxtb gxtb wb97m3c
      !write(*,'(5x,''-tslevel  [method] : select level for computing reaction barriers ("pbe","b973c","r2scan3c","pbeh3c","wb97x3c","kpr2scan50d4","wb97xd4tz","pbe0","ccsdt","gfn2","gfn1","pm6","gff")'')')
      !write(*,'(5x,''-iplevel  [method] : select level for computing IPs for charge assignment ("pbe","b973c","r2scan3c","pbeh3c","wb97x3c","kpr2scan50d4","wb97xd4tz","pbe0","ccsdt","gfn2","gfn1","pm6","gff")'')')
      !write(*,'(5x,''-ip2level  [method] : select level for computing IPs for charge assignment of critical cases with close IPs ("pbe","b973c","r2scan3c","pbeh3c","wb97x3c","kpr2scan50d4","wb97xd4tz","pbe0","gfn2","gfn1","pm6","dxtb","gff")'')')
      write (*, '(5x,''-nebnormal: use normal settings instead of loose settings for NEB search'')')
      write (*, '(5x,''-notsgeo  : do not use geodesic interpolation as guess for restarting not converged NEB runs '')')
      write (*, '(5x,''-tsnodes [integer]: select number of nodes for path (default 9) '')')
      write (*, *)
      write (*, '(/,1x,''Options for the fragment generation with CREST msreact:'')')
      write (*, '(5x,''-msnoiso: sort out non-dissociated fragments in crest msreact (i.e., rearranged structures) '')')
      write (*, '(5x,''-msiso: sort out dissociated fragments in crest msreact (i.e., rearranged structures) '')')
      write (*, '(5x,''-msfulliso: rearrangements also for subsequent fragmentations in crest msreact (activated by default)'')')
      write (*, '(5x,''-msnoattrh: deactivate attractive potential between hydrogen and LMO centers)'')')
      write (*, '(5x,''-msnshifts [int]: perform n optimizations with randomly shifted atom postions (default 0) '')')
      write(*,'(5x,''-msnshifts2 [int]: perform n optimizations with randomly shifted atom postions and repulsive potential applied to bonds (default 0) '')')
     write (*, '(5x ''-msnbonds [int]: maximum number of bonds between atoms pairs for applying repulsive potential (default 3)'')')
write (*, '(5x,''-msmolbar: sort out topological duplicates by molbar codes (activated by default - requires  sourced "molbar")'')')
      write (*, '(5x,''-msinchi: sort out topological duplicates by inchi codes (requires  sourced "obabel")'')')
      write (*, '(5x ''-msnfrag [int]: number of fragments that are printed by msreact (random selection)'')')
      write (*, '(5x ''-msfragdist [real]: seperate fragments before TS search from each other (default 2.5 [Angstrom]) '')')
      write (*, *)
      write (*, '(/,1x,''Special options:'')')
      write (*, '(5x,''-cneintscale  : cale internal energy of subsequent fragmentations according to number of atoms '')')
      write (*, '(5x,''-noKER  : do not compute kinetic energy release (KER) '')')
      write (*, '(5x,''-usetemp  : take G_mRRHO contribution instead of only ZPVE ( ZPVE is default)  '')')
      write (*, '(5x,''-scaleeinthdiss [real] this decreases the internal energy only for -H or -H2 abstractions (default  0.5)'')')
      write (*, '(5x,''-nsamples [int] number of simulated runs in Monte Carlo (default 10E05)'')')
      write (*, '(5x,''-rrkm : compute rat constants via RRKM equation instead of eyring (default eyring)'')')
      write (*, '(5x,''-sthr [int] : RRHO cutoff for thermo contribution  (default is 150 cm-1)'')')
      write (*, '(5x,''-ithr [int] : imaginary RRHO cutoff for thermo contribution  (default is 100 cm-1)'')')
      write (*, '(5x,''-nthermosteps [int] : number of increments to compute thermal corrections for IEE distribution (default 200, take a multiple of 10 )'')')
      write(*,'(5x,''-topocheck [string] : check topology after optimization (default "molbar", "inchi" with "obabel" in path also possible, but not thoroughly tested. Set empty (i.e.," ") to deactivate) '')')
      write (*, '(/,1x,''Options for plotting of mass spectrum:'')')
      write (*, '(5x,''-noisotope  : only plot peak with highest isotope propability'')')
      write (*, '(5x,''-int_masses  : only plot masses as integers'')')
      write (*, '(5x,''-mthr [real] m/z at which fragment is plotted (default 0)'')')
      write (*, '(5x,''-intthr [real] peak intensity threshold at which peak is plotted in percent (default 0)'')')
      write (*, *)
      ! write (*, '(/,1x,''for more options use  [-advanced] flag'')')
      stop '   [-h/-help] displayed. exit.'
   end subroutine printhelp

   subroutine printhelp_advanced
      implicit none

      write (*, '(/,1x,''Advanced Options for testing:'')')
      write (*, '(5x,''-cid: compute a CID spectrum (not yet implemented)'')')
      write (*, '(5x,''-tf [real] time of flight in spectrometer in mukroseconds default is 50 '')')
      !DEL write (*, '(5x,''-dxtbparam  : use special dxtb parameter - requires dxtb_param.txt in starting DIR '')')
      !DEL
      !DEL write (*, '(5x,''-reoptts  :reoptimize TS after path search (unstable so deactivated by default)) '')')
      !TODO DEL write (*, '(5x,''-atmassh [real] increase mass of hydrogen to scale R-H frequencies (default  1.08 u.)'')')
      !TODO write (*, '(5x,''-qmprog [program]      : external code for DFT calculations ("orca currently default")'')')
      write (*, '(5x,''-plotk  : plot k(E) curves '')')
      write (*, '(5x,''-nobhess  : do not compute thermo correction for barrier with bhess '')') ! currently not implemented
      write (*, '(5x,''-path [method]     : select pathfinder methods (default is "neb", "gsm" or "xtb" also possible)'')')
      write (*, '(5x,''-fermi : apply  fermi smearing (deactivated by default)'')')
      write (*, '(5x,''-tfscale [real] : scale time of flight for subseauent fragmentations (default 1.0 means no scaling)'')')
      write (*, '(5x,''-scaleker [real] : scale kinetic energy release upon fragmentation (default 1.0 means no scaling)'')')
      write (*, '(5x,''-sumreacscale [real] : scale sumreac in mcsimu (default 1.0 means no scaling)'')')
      write (*, '(5x,''-erelaxtime [real] : assumend relaxations time for H-diss scaling, very large number means always scaling, small only in first frag scaling, give in picoseconds (default 5 picoseconds)'')')
      write(*,'(5x,''-noeatomscale : do not scale IEE of subsequent fragmentations accoding to atom number ratio of generated fragments and charged fragment gets complete energy)'')')
      write(*,'(5x,''-sortoutcascade : sort out strucutres which have more than one maximum in the reaction path (turned of by default, as we can loose important fragments by this)'')')
      write (*, '(5x,''-fixe [real] : just use one energy in eV (default 0 means distribution is taken instead)'')')
      write (*, '(5x,''-chargeatomscale : larger fragment gets charge, IPs are ignored)'')')
      write (*, '(5x,''-hotip : compute IPs on unoptimized structures'')')
      write (*, '(5x,''-ignoreip : ignoreips, both fragments get full intensity'')')
      write (*, '(5x,''-noirc : set all imaginary freqs to 100 (only in esim mode)'')')
      write (*, '(5x,''-printlevel [int]: printout settings - 1: default, 2: verbose, 3: debug mode (very verbose!!!)'')')
   write (*, '(5x,''-keepdir : do not delete directories of sorted out fragments (for diagnostic purposes deactivated for now))'')')
      write (*, '(5x,''-exstates [int] : number of excited states to include via TDDFT (default 0)'')')
      write (*, '(5x,''-noplotisos : do not plot isomers in mass spectrum (default false)'')')
      write (*, '(5x,''-ircrun : search in hessian calculation directory for an IRC'')')
      write (*, '(5x,''-picktsrun : search reaction path for maxima/potential TSs'')')
      write (*, *)
      write (*, '(/,1x,''CID mode specific options (!!!!!Warning: highly experimental!!!!):'')')
      write (*, '(5x,''-cidauto : auto mode for CID simulation (default) '')')
      write (*, '(5x,''-cidtemprun : CID simulation without collisions (only thermal heating in ESI simulation) '')')
      write (*, '(5x,''-cidforced : CID simulation with only collisions (no thermal heating in ESI simulation) '')')
      write (*, '(5x,''-elab [real] : collision energy in laboratory frame in eV (default 40 eV) '')')
      write (*, '(5x,''-esi [real] : internal energy scaling in esi in eV (default dependent on number of atoms) '')')
 write (*, '(5x,''-esiatom [real] : internal energy scaling in esi per atom in eV (default  0.1 (for temprun 0.4) eV/per atom) '')')
      write (*, '(5x,''-esiw [real] : width internal energy scaling in esi in eV (default 0.2) '')')
      write (*, '(5x,''-collw [real] : width of collision energy distribution (default 0.5) '')')
      !write (*, '(5x,''-maxcoll [int] : maximum number of collisions (default 10) '')') ! TODO not yet implemented
      write (*, '(5x,''-lchamb [real] : length of collision chamber in m (default 0.25 m) '')')
      write (*, '(5x,''-mgas [real] : mass of gas in u (default Ar: 39.948 u) '')')
      write (*, '(5x,''-rgas [real] : vdw radius of gas in bohr (default Ar: 3.55266638 bohr) '')')
      write (*, '(5x,''-pgas [real] : pressure of gas in Pa (default Ar: 0.132 Pa) '')')
      write (*, '(5x,''-cidscool [real] : scale collisional cooling of ions in CID mode (default 1.0) '')')
      write (*, *)
      stop '   [-advanced] displayed. exit.'

   end subroutine printhelp_advanced

   ! check here input settings and print infos about run
   ! TODO add more checks, especially for settings that make no sense
   subroutine check_settings(env)
      type(runtypedata) :: env

      write (*, '(a)') "Settings:"
      write (*, '(60(''*''))')
      if (env%chrg .eq. 1) write (*, '(x,''           + Positive Ion mode +'')')
      if (env%chrg .eq. -1) write (*, '(x,''           - Negative Ion mode -'')')
      if (env%chrg .gt. 1 .or. env%chrg .lt. -1) then
         write (*, *) "Charge of input molecule was set to: ", env%chrg
         write (*, *) "!!!Warning multiply charged ions not yet tested!!!"
      end if

      write (*, '(a,x,a)') "IEE distribution:", env%edist
      write (*, '(a,1f5.2)') "width of IEE distribution (eimpw) is: ", env%eimpw
      write (*, '(a,1f5.2,x,a)') "energy per atom (ieeatm) is: ", env%ieeatm, "eV"
      write (*, '(a,1f5.2x,a)') "eimp0 is: ", env%eimp0, "eV"

     if (env%exstates > 0 .and. (.not. (env%tslevel == 'wb97x3c' .or. env%tslevel == 'r2scan3c' .or. env%tslevel == 'pbeh3c'))) then
         write (*, *) "Warning excited states via TD-DFT only possible for DFT functionals"
         error stop
      end if

      if (env%mode == "ei") then
         write (*, *) "spectral mode: EI-MS" !
      elseif (env%mode == "cid") then
         if (env%cid_mode == 1) then
            write (*, *) "spectral mode: CID auto (!WARNING: highly experimental!)"
         elseif (env%cid_mode == 2) then
            write (*, *) "spectral mode: CID temprun (no collisions simulated)"
         elseif (env%cid_mode == 3) then
            write (*, *) "spectral mode: CID forced collisions (!WARNING: not yet implemented)"
            write (*, *) "Aborting..."
            stop
         end if
      end if

      if (env%chrg .lt. 0) write (*, *) "negative ion mode chosen, Hybrid DFT with diffuse basis sets recommended to use!!!"
      if (env%qmprog .ne. "orca") write (*, *) "Warning! currently only ORCA is supported as external QM program"

      if (env%geolevel == "gxtb" .and. env%tsfinder == 'neb') then
         write (*, *) "Warning: special xTB version for NEB search necessary"
      end if

      write (*, '(a)') "QM methods used:"
      write (*, '(2x,a,x,a)') "Level for geometry optimizations and path searches:", trim(env%geolevel)
      write (*, '(2x,a,x,a)') "Level for reaction energy and barrier calculations:", trim(env%tslevel)
      write (*, '(2x,a,x,a)') "Level for IP prescreening:", trim(env%iplevel)
       if (env%ip2level .ne. '' .and. trim(env%ip2level) .ne. trim(env%iplevel)) write (*,'(2x,a,x,a)') "Level for IP refinement:", trim(env%ip2level)

      write (*, *) "Runtime settings used:"
      write (*, '(x,a,x,i0)') "Number of cores used:", env%cores
      write (*, '(x,a,x,i0)') "Number of allowed subsequent fragmentations:", env%nfrag
      write (*, '(x,a,f4.2,a)') "Intensity threshold for further fragmentation: ", env%pthr, "%"

      !write(*,*) "Number of nodes for path search:", env%tsnds
      if (env%exstates .gt. 0) write (*, *) "Number of excited states to include via TDDFT:", env%exstates
      write (*, '(60(''*''))')

   end subroutine check_settings

   subroutine check_progs(env)
      use iomod
      implicit none
      type(runtypedata) :: env
      character(len=80), dimension(4) :: levels
      character(len=1024) :: xtbpath
      integer :: i, io
      logical :: ex

      ! CHECK if external programs are available
      write (*, *) "External programs used:"
      !always needed
      call check_prog('xtb', .true., io)
      call check_prog('crest', .true., io)
      call check_prog('orca', .true., io)
      ! for xtb with orca TODO use setenv for this
      call execute_command_line("echo $(which xtb) > tmp.xtbpath")
      call rdshort_string('tmp.xtbpath', xtbpath)
      call remove('tmp.xtbpath')
      io = setenv('XTBEXE', trim(xtbpath))
      !call execute_command_line('export XTBEXE=$(which xtb)')
      call execute_command_line('echo "xtb path set for orca: $XTBEXE"')

      ! check python program paths
      if (env%msmolbar .or. env%topocheck == "molbar") then
         call check_prog('molbar', .true., io)
      end if

      if (env%tsgeodesic) then
         call check_prog('geodesic_interpolate', .true., io)
      end if

      levels = [env%tslevel, env%iplevel, env%ip2level, env%geolevel]
      ! special external programs
      do i = 1, 4
         if (levels(i) == 'gxtb') then
            call check_prog('gxtb', .true., io)
            inquire (file='~/.gxtb', exist=ex)
            if (.not. ex) then
               write (*, *) "Warning: requires ~/.gxtb file"
               error stop
            end if
            inquire (file='~/.ceh', exist=ex)
            if (.not. ex) then
               write (*, *) "Warning: requires ~/.ceh  file"
               error stop
            end if
            inquire (file='~/.basisq', exist=ex)
            if (.not. ex) then
               write (*, *) "Warning: requires ~/.basisq file"
               error stop
            end if
            exit
         end if
      end do

      write (*, '(60(''''))')
      write (*, *)
   end subroutine check_progs

end module argparser
