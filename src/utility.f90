! this module contains some utility functions
! Some of these useful subroutines are taken from CREST (github.com/crest-lab/crest)
! (using the GPL-3.0 license like QCxMS2)
! and are slightly modified for QCxMS2
module utility
   use xtb_mctc_accuracy, only: wp
   use qcxms2_data
   use iomod
   implicit none
contains

!=====================================================================!
! a general version of an omp-parallel calculation loop
!
! On Input: ndirs  - number of dirs to compute
!           base      - base name of the dirs, e.g., trim(env%path)//"/p"
!           jobcall   - the system call to be executed in all dirs
!           niceprint - boolean for nicer printout
!=====================================================================!
!simple loop to perform a job in parallel in multiple directories
   subroutine omp_jobcall(ndirs, basedir, jobcall)

      implicit none
      integer, intent(in) :: ndirs !number of directories
      character(len=1024), intent(in) :: jobcall

      character(len=512) :: tmppath
      character(len=*) :: basedir ! for example trim(env%path)//p
      integer :: vz, i, k, io
      k = 0
      !$omp parallel &
!$omp shared( vz,jobcall,ndirs,k)
!$omp single  !! done by single thread

      do i = 1, ndirs
         vz = i
         !$omp task firstprivate( vz ) private( tmppath,io )
         !$omp critical
         write (tmppath, '(a,i0)') trim(basedir), vz
         !$omp end critical
         call execute_command_line('cd '//trim(tmppath)//' && '//trim(jobcall), exitstat=io)
         !$omp critical
         k = k + 1
         !prints how many are finished but not which are still running
         write (6, '(1x,i0)', advance='no') k
         flush (6)
         !$omp end critical
         !$omp end task
      end do
!$omp taskwait
!$omp end single
!$omp end parallel
   end subroutine omp_jobcall

!loop to perform jobs in parallel
!jobcall is always the same job
! dirs is list of dirs where job should be executed relative to starting directory
   subroutine omp_samejobcall(njobs, dirs, jobcall, print)
      implicit none
      integer, intent(in) :: njobs
      character(len=1024), intent(in) :: jobcall
      !character(len=1024), intent(in) :: dirs(njobs)
      character(len=80), intent(in) :: dirs(:)
      character(len=1024) :: tmpjob
      character(len=2024) :: tmppath
      character(len=1024) :: pwd
      logical, optional :: print
      logical :: lprint
      integer :: i, io, vz, k
      integer :: shift ! first element often zero in dirs array

      ! print progress
      if (present(print)) then
         lprint = print
      else
         lprint = .true.
      end if

      ! prevent error here
      if (njobs .le. 0) return

      if (dirs(1) .eq. '') then
         shift = 1
      else
         shift = 0
      end if

      call getcwd(pwd)

      k = 0
!$omp parallel &
!$omp shared(pwd, vz,jobcall, dirs,njobs,k,shift,lprint)
!$omp single  !! done by single thread
      do i = 1, njobs
         vz = i + shift
         !$omp task firstprivate( vz ) private(tmppath,io )
         !$omp critical
         write (tmppath, '(a)') trim(pwd)//"/"//trim(dirs(vz))
         !write (tmppath, '(a)') trim(dirs(vz))
         !$omp end critical
         call execute_command_line('cd '//trim(tmppath)//' && '//trim(jobcall), exitstat=io)
         !$omp critical
         k = k + 1
         !prints how many are finished but not which are still running
         if (lprint) write (6, '(1x,i0)', advance='no') k
         flush (6)
         !$omp end critical
         !$omp end task
      end do
!$omp taskwait
!$omp end single  !! done by single thread
!$omp end parallel
      write (*, *)
   end subroutine omp_samejobcall

! threads is flexibel omp is fixed
   subroutine setompthreads(env, njobs)
      implicit none
      type(runtypedata) :: env
      character(len=200) :: exportcall
      integer, intent(in) :: njobs
      integer :: ncalc
      integer :: threads ! number of threads for each job
      integer :: omp ! number of parallel jobs
      integer :: io

      integer :: envompthreads, envparnodes

      if (njobs .eq. 0) then ! prevent error here integer divide by zero
         ncalc = 1
      else
         ncalc = njobs
      end if

      env%threads = int(env%cores/ncalc)
      if (env%threads .eq. 0) env%threads = 1
      !odd number of threads not efficient on intel CPUs, we go e.g. from 3 to 4 then
      if (env%threads .gt. 1 .and. modulo(env%threads, 2) .ne. 0) then
         if (env%threads .lt. 4) then ! up to 4 threads should be efficient, above we should reduce number of cores
            env%threads = env%threads + 1
         else
            env%threads = env%threads - 1
         end if
      end if
      if (env%threads .eq. 2) env%threads = 4 ! TODO FIXME 2 cores currently in ORCA buggy?

      !modulo check if cores are efficiently used ??
      !if (env%threads .gt. 28) env%threads = 26 ! 26 cores should be enough for calculations 28 can lead to bug ...
      if (env%threads .gt. 28) env%threads = 27! 26 cores should be enough for calculations 28 can lead to bug ...
      !FOR MARVIN everything on one CPU for cpus per task
!      env%threads = 1

      omp = int(env%cores/env%threads) ! ORCA can handle also many cores
      call OMP_Set_Num_Threads(omp)

      ! for slurm to be safe compute everything which would be cpu per taks on one core
      !envompthreads = env%threads
      !      if (envompthreads .gt. 8) envompthreads = 8 ! 8 cores should be enough for calculations with xtb (molbar??)
      envompthreads = 1
      io = setenv('OMP_NUM_THREADS', envompthreads)
      io = setenv('MKL_NUM_THREADS', envompthreads)
      ! for TURBOMOLE
      !envparnodes = env%threads
      envparnodes = env%threads

      if (envparnodes .gt. 28) envparnodes = 28 ! 28 cores should be enough for calculations with TM
      io = setenv('PARNODES', envparnodes)

   end subroutine setompthreads

   !for GSM and molbar, we need to set OMP_NUM_THREADS and MKL_NUM_THREADS to 1
   subroutine setomptoone()
      implicit none
      integer :: io

      io = setenv('OMP_NUM_THREADS', 1)
      io = setenv('MKL_NUM_THREADS', 1)
      io = setenv('PARNODES', 1)

   end subroutine setomptoone

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine velectrons_amount(nat, iat, chrg, nvel)

      integer, intent(in) :: nat, iat(nat), chrg
      integer, intent(out) :: nvel
      integer :: i
      real(wp)  :: z(nat)

      ! Set z as valence electrons
      do i = 1, nat
         z(i) = iat(i) - get_core_e(iat(i))
      end do

      nvel = idint(sum(z))

      nvel = nvel - chrg

      ! nb = nel/2          !nb= half the electrons, for uneven e-

   end subroutine velectrons_amount

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine electrons_amount(nat, iat, chrg, nel, nb, z)

      integer, intent(in) :: nat, iat(nat), chrg
      integer, intent(out) :: nel, nb
      integer :: i
      real(wp)  :: z(nat)

      ! Set z as valence electrons
      do i = 1, nat
         z(i) = iat(i) - get_core_e(iat(i))
      end do

      nel = idint(sum(z))
      nel = nel - chrg

      nb = nel/2          !nb= half the electrons, for uneven e-

   end subroutine electrons_amount
   function get_core_e(iat) result(ncore)
      integer :: iat, ncore

      if (iat <= 2) then
         ncore = 0
      elseif (iat <= 10) then
         ncore = 2
      elseif (iat <= 18) then
         ncore = 10
      elseif (iat <= 29) then   !zn
         ncore = 18
      elseif (iat <= 36) then
         ncore = 28
      elseif (iat <= 47) then
         ncore = 36
      elseif (iat <= 54) then
         ncore = 46
      elseif (iat <= 71) then
         ncore = 54
      elseif (iat <= 79) then
         ncore = 68
      elseif (iat <= 86) then
         ncore = 78
      end if
   end function get_core_e

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! calculate # of valence electrons in mopac and dftb
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine valel(at, el)

      integer  :: at
      real(wp) :: el

      if (at .le. 2) then
         el = at
      elseif (at .le. 10) then
         el = at - 2
      elseif (at .le. 18) then
         el = at - 10
      elseif (at .le. 36) then
         el = at - 18
         if (at .gt. 28) el = at - 28
      end if

   end subroutine valel

   subroutine printspin(env, fname, spin, chrgin)
      use mctc_env, only: error_type, get_argument!, fatal_error
      use mctc_io, only: structure_type, read_structure, write_structure, &
         & filetype, get_filetype, to_symbol, to_number
      use qcxms2_data
      implicit none
      character(len=*)      :: fname
      integer, intent(out) :: spin
      type(runtypedata) :: env
      integer :: nat, chrg
      integer, optional :: chrgin
      integer, allocatable :: iat(:)
      type(structure_type) :: mol
      type(error_type), allocatable :: error

      call read_structure(mol, fname, error, filetype%xyz)
      if (allocated(error)) then
      if (error%stat .eq. 1) then ! error reading structure
         write (*, *) "ERROR: Could not find file", fname, " for spin determination"
         spin = -1
         return
      end if
      end if
      nat = mol%nat
      if (present(chrgin)) then
         chrg = chrgin
      else
         chrg = env%chrg
      end if
      allocate (iat(nat))
      iat = mol%num(mol%id)
      call getspin(nat, iat, chrg, spin)
   end subroutine printspin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! spin value of molecule
! returns -1 if no electrons are found (e.g. H+)
! assumes always low-spin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine getspin(nat, iat, chrg, isp)

      integer nat, i, j, chrg, isp, iat(nat)

      j = 0

      do i = 1, nat
         j = j + iat(i)
      end do

      j = j - abs(chrg)
      !number of electrons
      isp = 1 + mod(j, 2)

      if (j < 1) isp = -1

   end subroutine getspin

   !=============================================================!
   ! classical quicksort algorithm, sort LOW-to-HIGH ! taken from CREST
   !=============================================================!
   recursive subroutine quicksort(n, arr)
      implicit none
      integer :: n, arr(n), i, j, k, m
      integer :: pivot
      integer, allocatable :: R(:), L(:)
      integer :: rr, ll, rc, lc, pp

      if (n .le. 1) return

      pivot = arr(1)
      if (arr(2) .lt. arr(1)) pivot = arr(2)
      pp = 0
      do i = 1, n
         if (arr(i) .eq. pivot) pp = pp + 1
      end do

      ll = 0
      do i = 1, n
         if (arr(i) .le. pivot) then
            ll = ll + 1
         end if
      end do
      ll = ll - pp
      rr = n - ll - pp
      allocate (L(ll), R(rr))

      lc = 0
      rc = 0
      do j = 1, n
         if (arr(j) .lt. pivot) then
            lc = lc + 1
            L(lc) = arr(j)
         else if (arr(j) .gt. pivot) then
            rc = rc + 1
            R(rc) = arr(j)
         end if
      end do

      call quicksort(ll, L)
      call quicksort(rr, R)

      do i = 1, ll
         arr(i) = L(i)
      end do
      do k = 1, pp
         m = k + ll
         arr(m) = pivot
      end do
      do j = 1, rr
         m = j + ll + pp
         arr(m) = R(j)
      end do

      deallocate (R, L)
   end subroutine quicksort

   !=============================================================!
   ! other variant of quicksort algos for real
   !=============================================================!
   recursive subroutine qsort(a, first, last, ind)
      implicit none
      real(wp) :: a(:), x, t
      integer ind(*)
      integer first, last
      integer i, j, ii

      x = a((first + last)/2)
      i = first
      j = last
      do
         do while (a(i) < x)
            i = i + 1
         end do
         do while (x < a(j))
            j = j - 1
         end do
         if (i >= j) exit
         t = a(i); a(i) = a(j); a(j) = t
         ii = ind(i); ind(i) = ind(j); ind(j) = ii
         i = i + 1
         j = j - 1
      end do
      if (first < i - 1) call qsort(a, first, i - 1, ind)
      if (j + 1 < last) call qsort(a, j + 1, last, ind)
   end subroutine qsort

   subroutine maskinvert(nall, mask)
      implicit none
      integer :: nall
      integer :: mask(nall)
      integer, allocatable :: imask(:)
      integer :: i
      allocate (imask(nall))
      do i = 1, nall
         imask(mask(i)) = i
      end do
      mask(:) = imask(:)
      deallocate (imask)
      return
   end subroutine maskinvert

!simple formula to approximate temperature
!doi.org/10.1063/5.0090205 very simple formula
   subroutine settemp(env, fname)
      implicit none
      type(runtypedata) :: env
      character(len=*) :: fname
      character(len=80) :: query
      integer :: nat, nvib
      integer :: chrg, uhf
      real(wp) :: sumreac, E
      logical :: ldum

      call rdshort_int('.CHRG', chrg)
      call rdshort_int('.UHF', uhf)

      ! read sum of reaction energies
      write (query, '(a,1x,i0,1x,i0)') trim(env%tslevel)//" sumreac", chrg, uhf
      call grepval('qmdata', trim(query), ldum, sumreac)

      ! read number of atoms
      call rdshort_int(trim(fname), nat)

      ! to get thermal energy of molecule in microcanonical ensemble
      ! just get here average energy minus sum of reaction energies
      E = env%ieeatm*nat - sumreac
      nvib = (3*nat - 6) ! number of harmonic oscillators

      env%temp = E/(nvib*8.617333262e-05_wp)
      if (env%temp < 0) then
         env%temp = 0
      end if
      write (*, *) "Average temperature set to", env%temp, " K"
   end subroutine settemp

!simple formula to approximate temperature
!doi.org/10.1063/5.0090205 very simple formula
   function calctemp(nvib, energy) result(T)
      implicit none
      integer, intent(in) :: nvib
      real(wp), intent(in) :: energy ! ev
      real(wp) :: T ! K

      T = energy/(nvib*8.617333262e-05_wp)
      if (T < 0) then
         T = 0.01_wp
         ! T = 298  ! assume at least room temperature ??? TODO DEL
      end if
   end function calctemp

   subroutine checkatom(fname, isatom)
      implicit none
      character(len=*) :: fname
      integer :: nat
      logical, intent(out) :: isatom

      ! read number of atoms
      call rdshort_int(trim(fname), nat)
      if (nat .eq. 1) then
         isatom = .true.
      else
         isatom = .false.
      end if

   end subroutine checkatom
! getsumform for xyz file
   subroutine getsumform(fname, sumformula)
      use mctc_env, only: error_type, get_argument!, fatal_error
      use mctc_io, only: structure_type, read_structure, write_structure, &
         & filetype, get_filetype, to_symbol, to_number
      implicit none
      type(runtypedata) :: env
      integer :: nat
      integer, allocatable :: iat(:)
      character(len=*), intent(in) :: fname
      character(len=80), intent(out) :: sumformula
      type(structure_type) :: mol
      type(error_type), allocatable :: error

      call read_structure(mol, fname, error, filetype%xyz)
      nat = mol%nat
      allocate (iat(nat))
      iat = mol%num(mol%id)
      sumformula = trim(sumform(nat, iat))

   end subroutine getsumform
!===================================================!
! get sumformula as a string from the AT array
!===================================================!
   character(len=80) function sumform(nat, at)
      implicit none
      integer :: nat
      integer :: at(nat)
      integer :: sumat(94)
      integer :: i
      character(len=6) :: str
      sumform = ''
      sumat = 0
      do i = 1, nat
         sumat(at(i)) = sumat(at(i)) + 1

      end do
      do i = 1, 94
         if (sumat(i) .lt. 1) cycle
         write (str, '(a,i0)') trim(adjustl(i2e(i, 'nc'))), sumat(i)
         sumform = trim(sumform)//trim(str)
      end do
      return
   end function sumform
!============================================================!
! e2i is used to map the element (as a string) to integer
!============================================================!
   integer function e2i(cin)
      implicit none
      character(len=*), intent(in) :: cin
      character(len=:), allocatable :: c
      integer :: iout
      c = trim(convertlable(cin))
      select case (c)
      case ('H'); iout = 1
      case ('D'); iout = 1
      case ('T'); iout = 1
      case ('HE'); iout = 2
      case ('LI'); iout = 3
      case ('BE'); iout = 4
      case ('B'); iout = 5
      case ('C'); iout = 6
      case ('N'); iout = 7
      case ('O'); iout = 8
      case ('F'); iout = 9
      case ('NE'); iout = 10
      case ('NA'); iout = 11
      case ('MG'); iout = 12
      case ('AL'); iout = 13
      case ('SI'); iout = 14
      case ('P'); iout = 15
      case ('S'); iout = 16
      case ('CL'); iout = 17
      case ('AR'); iout = 18
      case ('K'); iout = 19
      case ('CA'); iout = 20
      case ('SC'); iout = 21
      case ('TI'); iout = 22
      case ('V'); iout = 23
      case ('CR'); iout = 24
      case ('MN'); iout = 25
      case ('FE'); iout = 26
      case ('CO'); iout = 27
      case ('NI'); iout = 28
      case ('CU'); iout = 29
      case ('ZN'); iout = 30
      case ('GA'); iout = 31
      case ('GE'); iout = 32
      case ('AS'); iout = 33
      case ('SE'); iout = 34
      case ('BR'); iout = 35
      case ('KR'); iout = 36
      case ('RB'); iout = 37
      case ('SR'); iout = 38
      case ('Y'); iout = 39
      case ('ZR'); iout = 40
      case ('NB'); iout = 41
      case ('MO'); iout = 42
      case ('TC'); iout = 43
      case ('RU'); iout = 44
      case ('RH'); iout = 45
      case ('PD'); iout = 46
      case ('AG'); iout = 47
      case ('CD'); iout = 48
      case ('IN'); iout = 49
      case ('SN'); iout = 50
      case ('SB'); iout = 51
      case ('TE'); iout = 52
      case ('I'); iout = 53
      case ('XE'); iout = 54
      case ('CS'); iout = 55
      case ('BA'); iout = 56
      case ('LA'); iout = 57
      case ('CE'); iout = 58
      case ('PR'); iout = 59
      case ('ND'); iout = 60
      case ('PM'); iout = 61
      case ('SM'); iout = 62
      case ('EU'); iout = 63
      case ('GD'); iout = 64
      case ('TB'); iout = 65
      case ('DY'); iout = 66
      case ('HO'); iout = 67
      case ('ER'); iout = 68
      case ('TM'); iout = 69
      case ('YB'); iout = 70
      case ('LU'); iout = 71
      case ('HF'); iout = 72
      case ('TA'); iout = 73
      case ('W'); iout = 74
      case ('RE'); iout = 75
      case ('OS'); iout = 76
      case ('IR'); iout = 77
      case ('PT'); iout = 78
      case ('AU'); iout = 79
      case ('HG'); iout = 80
      case ('TL'); iout = 81
      case ('PB'); iout = 82
      case ('BI'); iout = 83
      case ('PO'); iout = 84
      case ('AT'); iout = 85
      case ('RN'); iout = 86
      case default; iout = 0
      end select
      e2i = iout
   end function e2i
!============================================================!
! i2e is used to map the element (as a integer) to a string
!============================================================!
   character(len=2) function i2e(iin, oformat)
      implicit none
      integer, intent(in) :: iin
      character(len=:), allocatable :: c
      character(len=*), optional :: oformat
      select case (iin)
      case (1); c = 'H'
      case (2); c = 'HE'
      case (3); c = 'LI'
      case (4); c = 'BE'
      case (5); c = 'B'
      case (6); c = 'C'
      case (7); c = 'N'
      case (8); c = 'O'
      case (9); c = 'F'
      case (10); c = 'NE'
      case (11); c = 'NA'
      case (12); c = 'MG'
      case (13); c = 'AL'
      case (14); c = 'SI'
      case (15); c = 'P'
      case (16); c = 'S'
      case (17); c = 'CL'
      case (18); c = 'AR'
      case (19); c = 'K'
      case (20); c = 'CA'
      case (21); c = 'SC'
      case (22); c = 'TI'
      case (23); c = 'V'
      case (24); c = 'CR'
      case (25); c = 'MN'
      case (26); c = 'FE'
      case (27); c = 'CO'
      case (28); c = 'NI'
      case (29); c = 'CU'
      case (30); c = 'ZN'
      case (31); c = 'GA'
      case (32); c = 'GE'
      case (33); c = 'AS'
      case (34); c = 'SE'
      case (35); c = 'BR'
      case (36); c = 'KR'
      case (37); c = 'RB'
      case (38); c = 'SR'
      case (39); c = 'Y'
      case (40); c = 'ZR'
      case (41); c = 'NB'
      case (42); c = 'MO'
      case (43); c = 'TC'
      case (44); c = 'RU'
      case (45); c = 'RH'
      case (46); c = 'PD'
      case (47); c = 'AG'
      case (48); c = 'CD'
      case (49); c = 'IN'
      case (50); c = 'SN'
      case (51); c = 'SB'
      case (52); c = 'TE'
      case (53); c = 'I'
      case (54); c = 'XE'
      case (55); c = 'CS'
      case (56); c = 'BA'
      case (57); c = 'LA'
      case (58); c = 'CE'
      case (59); c = 'PR'
      case (60); c = 'ND'
      case (61); c = 'PM'
      case (62); c = 'SM'
      case (63); c = 'EU'
      case (64); c = 'GD'
      case (65); c = 'TB'
      case (66); c = 'DY'
      case (67); c = 'HO'
      case (68); c = 'ER'
      case (69); c = 'TM'
      case (70); c = 'YB'
      case (71); c = 'LU'
      case (72); c = 'HF'
      case (73); c = 'TA'
      case (74); c = 'W'
      case (75); c = 'RE'
      case (76); c = 'OS'
      case (77); c = 'IR'
      case (78); c = 'PT'
      case (79); c = 'AU'
      case (80); c = 'HG'
      case (81); c = 'TL'
      case (82); c = 'PB'
      case (83); c = 'BI'
      case (84); c = 'PO'
      case (85); c = 'AT'
      case (86); c = 'RN'
      case default; c = 'XX'
      end select
      i2e = trim(c)
      if (present(oformat)) then
         select case (oformat)
         case ('lc', 'lowercase')
            i2e = lowerCase(trim(c))
         case ('nc', 'nicecase')
            if (len_trim(c) .gt. 1) then
               c(2:2) = lowerCase(c(2:2))
               i2e = trim(c)
            end if
         case default
            continue
         end select
      end if
   end function i2e
!============================================================!
! convert a string into lowercase
!============================================================!
   function lowerCase(s)
      implicit none
      character(len=*), intent(in) :: s
      character(len=:), allocatable :: sout
      character(len=:), allocatable :: lowerCase
      integer :: ic, i
      character(26), Parameter :: high = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
      sout = s
      do i = 1, LEN_TRIM(s)
         ic = INDEX(high, s(i:i))
         if (ic > 0) sout(i:i) = low(ic:ic)
      end do
      call move_alloc(sout, lowerCase)
   end function lowerCase

   !============================================================!
! split element lable if some isotope indicator was given
! and convert to uppercase
!============================================================!
   function convertlable(s)
      implicit none
      character(len=*), intent(in) :: s
      character(len=:), allocatable :: sout
      character(len=:), allocatable :: convertlable
      integer :: ic, i
      character(14), parameter :: lab = '0123456789*_+-'
      character(26), parameter :: high = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
      sout = s
      do i = 1, len_trim(s)
         ic = index(lab, s(i:i))
         if (ic > 0) sout(i:i) = ' '
         ic = index(low, s(i:i))
         if (ic > 0) sout(i:i) = high(ic:ic)
      end do
      sout = trim(adjustl(sout))
      call move_alloc(sout, convertlable)
   end function convertlable

   ! quick routine for debugging
   subroutine printpwd(string)
      implicit none
      character(len=1024) :: pwd
      character(len=*), intent(in), optional :: string

      call getcwd(pwd)
      if (present(string)) then
         write (*, *) "current working directory is", string, " ", trim(pwd)
      else
         write (*, *) "current working directory is", trim(pwd)
      end if
   end subroutine printpwd

   ! change element symbols in xyz file from capital notation to normal notation (e.g. "CL" -> "Cl")
   subroutine rewrite_xyz_elements(fname)
      implicit none
      character(len=*) :: fname
      character(len=256) :: line, output_line
      character(1) :: second_letter
      integer :: i, j, k, ich1, ich2, io, nat

      call rdshort(fname, nat)

      open (newunit=ich1, file=trim(fname), status='old', iostat=io)
      open (newunit=ich2, file='tmp', status='new')
      read (ich1, *) line
      write (ich2, *) trim(line)
      read (ich1, *) line
      write (ich2, *) trim(line)

      do i = 1, nat
         if (io .ne. 0) then
            write (*, *) "ERROR: Could not find file", fname, " for element rewriting"
            close (ich1)
            close (ich2)
            return
         end if
         read (ich1, '(a)') line
         if (len(line) .gt. 2) then
            do j = 1, len(line)
               if (line(j:j) .ne. ' ') then
                  k = j + 1
                  line(k:k) = lowerCase(line(k:k))
                  exit
               end if
            end do
         end if

         write (ich2, '(a)') trim(line)
      end do
      close (ich1)
      close (ich2)

      call move('tmp', trim(fname))
   end subroutine rewrite_xyz_elements

   ! todo delete this abomination TODO FIXME other subroutine much better !!!!
   subroutine sortout0elements(npairs_in, npairs_out, fragdirs_in, fragdirs_out, removedir)
      implicit none
      integer :: i, j
      integer, intent(in) :: npairs_in
      integer, intent(out) :: npairs_out
      character(len=80), intent(in) :: fragdirs_in(npairs_in, 3)
      character(len=80), allocatable, intent(out)   :: fragdirs_out(:, :)
      logical, intent(in) :: removedir

      npairs_out = 0
      !count number of 0 elements
      do i = 1, npairs_in
         if (index(fragdirs_in(i, 1), 'p') .ne. 0) then
            npairs_out = npairs_out + 1
         end if
      end do
      ! sort out 0 elements
      if (npairs_out .eq. 0) then
         write (*, *) "no products left"
         allocate (fragdirs_out(1, 3))
      else
         write (*, *) "remaining number of products is: ", npairs_out
         if (allocated(fragdirs_out)) then
            write (*, *) "deallocate fragdirs_out"
            deallocate (fragdirs_out)
         end if
         allocate (fragdirs_out(npairs_out, 3))

      end if
      fragdirs_out = ''
      j = 0
      do i = 1, npairs_in
         if (index(fragdirs_in(i, 1), 'p') .ne. 0) then
            j = j + 1
            fragdirs_out(j, 1) = fragdirs_in(i, 1)
            if (index(fragdirs_in(i, 3), 'p') .ne. 0) then ! fragmentpair!
               fragdirs_out(j, 2) = fragdirs_in(i, 2)
               fragdirs_out(j, 3) = fragdirs_in(i, 3)
            end if
         else
            if (removedir) then
               write (*, *) "Delete fragment directory", trim(fragdirs_in(i, 1))
               call rmrf(trim(fragdirs_in(i, 1)))
               if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
                  call rmrf(trim(fragdirs_in(i, 2)))
                  call rmrf(trim(fragdirs_in(i, 3)))
               end if
            end if
         end if
      end do
   end subroutine sortout0elements

   ! remove element from list
   ! todo make more elegant like in crestms
!subroutine remove_fragdir(npairs,fragdirs,index,removedir)
!        implicit none
!        integer :: i, j
!        integer, intent(inout) :: npairs
!        character(len=80), intent(inout) :: fragdirs(npairs_in,3)
!        logical, intent(in) :: removedir
!  do i = index, npairs-1
!    fragdirs(i,1) = fragdirs(i+1,1) ! replace the entry with the next one
!    if (index(fragdirs(i+1,3),'p') .ne. 0) then ! fragmentpair!
!        fragdirs(i,2) = fragdirs(i+1,2)
!        fragdirs(i,3) = fragdirs(i+1,3)
!    else
!        fragdirs(i,2) = ''
!        fragdirs(i,3) = ''
!    end if
!end do
!npairs = npairs - 1 ! ignore last entry
!fragdirs(npairs,1) = '' ! delete last entry
!if (index(fragdirs(npairs,3),'p') .ne. 0) then
!    fragdirs(npairs,2) = '' ! delete last entry
!    fragdirs(npairs,3) = '' ! delete last entry
!end if
!if (removedir) then
!    write(*,*) "Delete fragment directory", trim(fragdirs(index,1))
!    call rmrf(trim(fragdirs(index,1)))
!    if (index(fragdirs(index,3),'p') .ne. 0) then
!        call rmrf(trim(fragdirs(index,2)))
!        call rmrf(trim(fragdirs(index,3)))
!    end if
!end if
!end subroutine remove_fragdir

! returns true if topo has changed ! deprecated, not needed anymor TODO DEL
   subroutine checktopo(env, fname, changed)
      implicit none
      type(runtypedata) :: env
      logical, intent(out) :: changed
      character(len=80) :: fname
      character(len=100) :: ftopo
      character(len=1024) :: identifiertopo0, identifiertopo
      character(len=1024) :: jobcall
      character(len=1024) :: jobcall2 ! check when we need this
      integer :: ich
      integer :: nat
      logical :: ex

      call rdshort_int(fname, nat)
      if (nat .eq. 1) then
         changed = .false.
         return ! molbar does not work for single atoms and is anyway not needed
      end if

      call rewrite_xyz_elements(trim(fname))
      write (ftopo, '(a)') "topo"
      call rdshort_string(trim(ftopo), identifiertopo0)

      if (env%topocheck .eq. "molbar") then
         call cuttopology(identifiertopo0)
      end if

      if (env%topocheck .eq. "molbar") then
         write (jobcall, '(a)') "molbar "//trim(fname)//" > "//trim(ftopo)//" 2> /dev/null"
      elseif (env%topocheck .eq. "inchi") then
         write (jobcall, '(a)') "obabel -i xyz "//trim(fname)//" -o inchi >  "//trim(ftopo)//" 2> /dev/null"
      end if
      call execute_command_line(trim(jobcall))
      inquire (file=trim(ftopo), exist=ex)
      if (.not. ex) then
         write (*, *) "topo could not be calculated"
         call printpwd()
         stop
      end if

      call rdshort_string(trim(ftopo), identifiertopo)

      if (env%topocheck .eq. "molbar") then
         call cuttopology(identifiertopo)
      end if

      if (identifiertopo0 .ne. '' .and. identifiertopo .ne. '') then
         if (identifiertopo0 .ne. identifiertopo) then
            changed = .true.
            !DEL
            write (*, *) "topo has changed"
         else
            changed = .false.
            ! write(*,*) "topo remains the same"
         end if
      else ! if topo failed just assume it didnt change
         changed = .false.
      end if
   end subroutine checktopo

! give Boltzmann averaged energie for a set of energies in eV
   subroutine boltzmannweights(env, energies, nenergies, bwenergie)
      use xtb_mctc_constants
      use xtb_mctc_convert
      implicit none
      type(runtypedata) :: env
      real(wp), intent(in) :: energies(:)
      integer, intent(in) :: nenergies
      real(wp), intent(out) :: bwenergie
      real(wp), allocatable :: pop(:)
      real(wp) :: esum, f, temp, q
      integer :: i

      temp = env%temp
      esum = 0.0_wp
      f = temp*kB*autoev
      do i = 1, nenergies
         q = exp(-energies(i)/f)
         esum = esum + q
      end do
      bwenergie = 0.0_wp
      allocate (pop(nenergies))
      pop = 0.0_wp
      do i = 1, nenergies
         q = exp(-energies(i)/f)
         pop(i) = q/esum
         bwenergie = bwenergie + pop(i)*energies(i)
      end do
   end subroutine boltzmannweights

   ! check important fragments for printing ! TODO maybe add refinement step for them
   subroutine check_impfrags(env, npairs, fragdirs, impfrags, imp_fragdirs)
      implicit none
      integer, intent(in) :: npairs
      integer, intent(out) :: impfrags
      character(len=80), allocatable, intent(out) :: imp_fragdirs(:, :)
      character(len=80), allocatable :: fragdirs(:, :)
      integer :: i
      type(runtypedata) :: env
      logical :: ex, imp ! important fragment
      real(wp) :: pfrag1, pfrag2
      real(wp) :: pthr ! threshold for important fragments
      pthr = env%pthr ! 1 percent, maybe we take here more, 5 percent??

      allocate (imp_fragdirs(npairs, 3))
      imp_fragdirs = ''
      impfrags = 0
      do i = 1, npairs
         imp = .false.
         if (index(fragdirs(i, 3), 'p') .ne. 0) then
            call rdshort_real(trim(fragdirs(i, 2))//'/pfrag', pfrag1)
            call rdshort_real(trim(fragdirs(i, 3))//'/pfrag', pfrag2)
            if (pfrag1 .gt. pthr .or. pfrag2 .gt. pthr) imp = .true.
            ! one of the fragments has to be over the threshold
         else
            call rdshort_real(trim(fragdirs(i, 1))//'/pfrag', pfrag1)
            if (pfrag1 .gt. pthr) imp = .true.
         end if
         if (imp) then
            impfrags = impfrags + 1
            write (imp_fragdirs(impfrags, 1), '(a)') trim(fragdirs(i, 1))
            if (index(fragdirs(i, 3), 'p') .ne. 0) then
               write (imp_fragdirs(impfrags, 2), '(a)') trim(fragdirs(i, 2))
               write (imp_fragdirs(impfrags, 3), '(a)') trim(fragdirs(i, 3))
            end if
         end if
      end do

      return
   end subroutine check_impfrags

   ! cut topology of molbar output so that only topography is left
   ! check also for errors and set empty string if error is found
   subroutine cuttopology(identifiertopo)
      character(len=*), intent(inout) :: identifiertopo
      integer :: i, dashcount

      !DEL print *, "full molbar is", trim(identifiertopo)
      dashcount = 0
      do i = 1, len(trim(identifiertopo))
         if (identifiertopo(i:i) == '|') dashcount = dashcount + 1
         if (dashcount .eq. 5) then ! topography part goes up to 5th dash
            identifiertopo = identifiertopo(1:i)
            exit
         end if
      end do

      ! molbar failes sometimes, so we have to check if it is empty
      if (index(identifiertopo, 'Error') .ne. 0) then
         identifiertopo = ''
      end if

      !DEL print *, "topology of molbar is", trim(identifiertopo)

   end subroutine cuttopology

   subroutine append_char(array, element)
      implicit none
      character(len=80), allocatable, intent(inout) :: array(:)
      character(len=*), intent(in) :: element
      character(len=80), allocatable :: temp_array(:)
      integer :: old_size, i

      if (allocated(array)) then
         old_size = size(array)
         allocate (temp_array(old_size + 1))
         do i = 1, old_size
            write (temp_array(i), '(a)') array(i)
         end do
         write (temp_array(old_size + 1), '(a)') trim(element)
         call move_alloc(temp_array, array)
      else
         allocate (array(1))
         write (array(1), '(a)') trim(element)
      end if
   end subroutine append_char

   ! TODO move in another module for nicer structure
   ! compute relative rate constant with eyring equation
   function calc_releyring(energy, dbarrier, nvib) result(k)
      implicit none
      integer :: nvib
      real(wp) :: T
      real(wp) :: dbarrier ! in eV difference in barriers/energy
      real(wp) :: energy ! in eV
      real(wp), parameter :: h = 6.62606957e-34_wp ! J s
      real(wp), parameter :: kB = 1.3806488e-23_wp ! J / K
      real(wp) :: k ! rate in s-1
      real(wp), parameter :: R = 8.61733325916595e-5_wp ! eV / K

      T = calctemp(nvib, energy)
      ! if (barrier .lt. 0.0_wp) write(*,*) "WARNING NEGATIVE BARRIER"
      k = EXP(-dbarrier/(R*T))
      if (k/(kB*T/h) .gt. 1.0_wp) then
         k = 1.0_wp ! limit k to kB*T/h , should be limit for low barriers ...
      end if
   end function calc_releyring

   ! compute rate constant with eyring equation
   function calc_eyring(energy, barrier, nvib) result(k)
      implicit none
      integer :: nvib
      real(wp) :: T
      real(wp) :: barrier ! in eV
      real(wp) :: energy ! in eV
      real(wp), parameter :: h = 6.62606957e-34_wp ! J s
      real(wp), parameter :: kB = 1.3806488e-23_wp ! J / K
      real(wp) :: k ! rate in s-1
      real(wp) :: R = 8.61733325916595e-5_wp ! eV / K

      T = calctemp(nvib, energy)
      ! if (barrier .lt. 0.0_wp) write(*,*) "WARNING NEGATIVE BARRIER"
      k = kB*T/h*EXP(-barrier/(R*T))
      if (k .gt. kB*T/h) then
         k = kB*T/h ! limit k to kB*T/h , should be limit for low barriers ...
      end if
      ! if (k .gt. 10e20) then
      !  write(*,*) "k is very high k, barrier, energy, nvib, T", k, barrier, energy, nvib, T
      ! write(*,*) "we are setting it to 10^20"
      ! k = 10e20_wp ! limit k to 10^20 TODO find a better solution for these cases
      !end if
   end function calc_eyring

   subroutine calc_Eavg(nsamples, eiee, piee, Eavg)
      implicit none
      real(wp), intent(in) :: eiee(nsamples)
      real(wp), intent(in) :: piee(nsamples)
      real(wp), intent(out) :: Eavg
      integer, intent(in) :: nsamples
      integer :: i

      Eavg = 0.0_wp
      do i = 1, nsamples
         Eavg = Eavg + piee(i)*eiee(i)
      end do

      Eavg = Eavg/sum(piee)

   end subroutine calc_Eavg
! compute rate constant with simplified rrkm equation
   function calc_rrkm(energy, barrier, nvib, freq) result(k)
      implicit none
      integer :: nvib
      real(wp) :: T
      real(wp) :: barrier ! in eV
      real(wp) :: energy ! in eV
      real(wp) :: freq ! in s-1
      real(wp) :: k ! rate in s-1

      !alternative equation, leads in some tests to less distribution between peaks
      ! so low barrier peaks are getting more intensity
      !kabs(j)= freq *  ((IEE -barriers(j)) / IEE)** (nvib-1)
      ! Eq. 2.21 p. 39 Gross
      ! R. G. Cooks, J. H. Beynon, R. M. Caprioli, G. R. Lester - Metastable ions-Elsevier (1974) p.226

      freq = abs(freq)*299792458.0_wp*100 ! need to convert here from cm-1 to s-1 *c * 100

      k = freq*EXP(-(nvib - 1)*barrier/energy)

      if (k .gt. freq) then
         k = freq ! theoretical limit for reaction
      end if

      if (barrier .gt. energy) then
         k = 0.0_wp
      end if
      ! can happen for very low energies
      if (isnan(k)) k = 0.0_wp

   end function calc_rrkm

   function getindex(incr, maxiee, energy) result(ind)
      implicit none
      real(wp), intent(in) :: maxiee
      real(wp), intent(in) :: energy
      integer, intent(in) :: incr
      integer :: ind ! index of energy in array of barriers
      real(wp) :: estep

      estep = maxiee/(incr)
      ind = int((energy/estep))

      ! for higher energies we just take the highest temperature
      if (ind .gt. incr) ind = incr
      ! for lower energies we just take the lowest temperature
      if (ind .le. 0) ind = 1

   end function getindex

   subroutine printiee(nsamples, eiee, piee, fname)
      implicit none
      real(wp), intent(in) :: eiee(nsamples)
      real(wp), intent(in) :: piee(nsamples)
      character(len=*), intent(in) :: fname
      integer, intent(in) :: nsamples
      integer :: i, ich

      open (newunit=ich, file=fname, status='replace')
      do i = 1, nsamples
         write (ich, *) eiee(i), piee(i)
      end do
      close (ich)
   end subroutine printiee

   ! A simple routine to check if a program is available and stops QCxMS2 if not (taken from CREST)
   !=========================================================================================!
   !
   subroutine check_prog(pname, verbose, iostat)
      implicit none
      character(len=*) :: pname
      logical, intent(in), optional :: verbose
      integer, intent(out), optional :: iostat
      character(len=:), allocatable :: checkcall
      character(len=:), allocatable :: pipe
      character(len=:), allocatable :: pathcall
      integer :: io

      pipe = ' >/dev/null 2>/dev/null'
      checkcall = 'command -v '//trim(pname)//pipe
      call execute_command_line(trim(checkcall), exitstat=io)

      if (present(verbose)) then
         if (io .ne. 0 .and. verbose) then
            write (*, '(4x,a,a,a)') 'binary: "', trim(pname), '" not found, aborting QCxMS2'
            if (trim(pname) == "geodesic_interpolate") then
               write (*, *) "You can turn off the geodesic interpolation with the option -notsgeo"
            end if
            error stop 'Program not found, aborting QCxMS2'
         elseif (io .eq. 0 .and. verbose) then
            checkcall = 'which '//trim(pname)
            call execute_command_line(trim(checkcall), exitstat=io)
         end if
      end if

      if (present(iostat)) then
         iostat = io
      end if

      return
   end subroutine check_prog

! cleanup files for restart calculation that could mess up the calculation
   subroutine cleanup_qcxms2_for_restart(env)
      implicit none
      type(runtypedata) :: env

      call remove('pfrag')
      call execute_command_line('rm */sumreac_* > /dev/null 2>&1')
      call execute_command_line('rm pfrag > /dev/null 2>&1')
      call execute_command_line('rm */pfrag > /dev/null 2>&1')
      call execute_command_line('rm kerav > /dev/null 2>&1')
      call execute_command_line('rm */kerav > /dev/null 2>&1')
      call execute_command_line('rm keravold > /dev/null 2>&1')
      call execute_command_line('rm */keravold > /dev/null 2>&1')
      if (env%mode == "cid") then
         call execute_command_line('rm sumdekin > /dev/null 2>&1')
         call execute_command_line('rm */sumdekin > /dev/null 2>&1')
         call execute_command_line('rm sumdeint > /dev/null 2>&1')
         call execute_command_line('rm */sumdeint > /dev/null 2>&1')
         call execute_command_line('rm x_trav > /dev/null 2>&1')
         call execute_command_line('rm */x_trav > /dev/null 2>&1')
      end if

   end subroutine cleanup_qcxms2_for_restart

   ! read structures and get number of structures and atoms
! works only for xyz ensemble file that contains only xyz files
   subroutine readstruc(inp, nat, nstruc)
      implicit none
      character(len=*), intent(in) :: inp
      integer :: io, ich, i
      integer, intent(out) :: nat, nstruc
      open (newunit=ich, file=inp)
      read (ich, *, iostat=io) nat
      i = 1
      do
         read (ich, *, iostat=io)
         if (io < 0) exit
         i = i + 1
      end do
      nstruc = (i + 1)/(nat + 2)
   end subroutine

end module utility
