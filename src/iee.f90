! most of the code is taken and modified from qcxms (https://github.com/qcxms/qcxms) under the  LGPL_3.0 license
module qcxms_iee
   use xtb_mctc_accuracy, only: wp
   use qcxms_boxmuller
   use qcxms2_data
   use iomod
   implicit none

contains
   subroutine getenergy(env, nsamples, eiee, piee, plotlog)
      implicit none
      integer, intent(in) :: nsamples
      real(wp), allocatable, intent(out) :: eiee(:), piee(:)
      logical, intent(in) :: plotlog
      type(runtypedata) :: env

      if (env%mode == "ei") then
         call getiee(env, nsamples, eiee, piee, plotlog)
      elseif (env%mode == "cid") then
         call scaleesi(env, nsamples, eiee, piee, plotlog)
      end if

   end subroutine getenergy
   subroutine getiee(env, nsamples, eiee, piee, plotlog)
      implicit none
      integer, intent(in) :: nsamples
      real(wp), allocatable, intent(out) :: eiee(:), piee(:)
      real(wp) :: edum, ehomo, eimp0, eimpw, ieeel, randx
      real(wp) :: iee_a, iee_b, ieeatm, exc, pmax
      real(wp) :: energy, prob
      logical, intent(in) :: plotlog
      integer :: i, ich
      integer :: nat
      integer :: edistri
      type(runtypedata) :: env

      allocate (eiee(nsamples), piee(nsamples))
      eiee = 0.0_wp
      piee = 0.0_wp

      ! get necessary infos for the calculation of IEE
      call getspecs(env, nat, ieeel, ehomo)
      ieeatm = env%ieeatm
      edistri = 1 ! 1 for poisson 0 for gaussian
      eimp0 = env%eimp0 ! 70 is default
      exc = (eimp0 - ehomo)  ! maximum energy
      write (*, *) "Generate energy distribution for IEE"

      ! generate normal distributed IEE
      if (trim(env%edist) == 'gaussian') then
         write (*, *) "Use Gaussian distribution to approximate IEE"
         edistri = 0
      else
         edistri = 1
         write (*, *) "Poisson distribution to approximate IEE"
      end if

      call getieeab(iee_a, iee_b, ieeel, edistri, exc, nat, ieeatm, pmax)

      eimpw = env%eimpw ! 0.1 is default

      do i = 1, nsamples
         ! 1. generate an e- that can ionize any MO and vary it by boxuller
         do
            ! vary by normal distributed boxmuller random number
            Edum = vary_energies(eimp0, eimpw)
            if (Edum >= ehomo) exit
         end do

         ! 2. subtract E(HOMO) to get IEE
         Edum = Edum - ehomo

         ! the random IEE in the possible interval
         do
            call random_number(randx)

            energy = randx*Edum
            if (edistri == 0) then
               call gauss0(iee_a, iee_b, ieeel, energy, prob)
            elseif (edistri == 1) then
               call poiss0(iee_a, iee_b, ieeel, energy, prob)
            end if
            prob = prob/pmax
            call random_number(randx)
            if (prob >= randx) exit
         end do
         eiee(i) = energy
         piee(i) = prob

         ! just for testing
         if (env%fixe .ne. 0) then
            eiee = env%fixe
         end if
      end do
      if (plotlog) then
         open (newunit=ich, file='eimp.dat', status='replace')
         do i = 1, nsamples
            write (ich, *) eiee(i), piee(i)
         end do
         close (ich)
      end if
   end subroutine getiee

   subroutine getiee2(env, ieeel, iee_a, iee_b, eiee, piee, plotlog)
      implicit none
      logical, intent(in) :: plotlog
      type(runtypedata) :: env
      real(wp), intent(in) :: ieeel
      real(wp), intent(out) :: eiee, piee
      real(wp) :: edum, eimp0, eimpw, randx, ieeatm
      real(wp) :: ieemax, pmax, E_avg, ehomo
      real(wp) :: exc ! exc   = (eimp0 - ehomo) * autoev
      real(wp) :: iee_a, iee_b
      integer :: edistri, ich

      eimp0 = env%eimp0 ! 70 is default
      eimpw = env%eimpw
      ieeatm = env%ieeatm
      ehomo = env%ehomo

      ! generate normal distributed IEE

      do
         ! vary by normal distributed boxmuller random number
         Edum = vary_energies(eimp0, eimpw)
         ! if ( Edum >= ehomo ) exit
         if (Edum >= ehomo) exit
      end do

      ! 2. subtract E(HOMO) to get IEE
      Edum = Edum - ehomo

      ! the random IEE in the possible interval

      call random_number(randx)

      eiee = randx*Edum
      !if(edistri == 0)then
      ! call gauss0(iee_a,iee_b,ieeel,exc,E_rand(i))
      !else
      !call poiss0(iee_a,iee_b,ieeel,E_rand(i),P(i))
      call poiss0(iee_a, iee_b, ieeel, eiee, piee)
      !call gauss0(iee_a,iee_b,ieeel,eiee(i),piee(i))

      if (plotlog) then
         open (ich, file='eimp.dat', access='append')
         write (ich, *) eiee, piee
         close (ich)
      end if
   end subroutine getiee2

   ! get necessary information of your input molecule for IEE generation
   subroutine getspecs(env, nat, ieeel, ehomo)
      use mctc_env, only: error_type, get_argument!, fatal_error
      use mctc_io, only: structure_type, read_structure, write_structure, &
         & filetype, get_filetype, to_symbol, to_number
      use utility
      implicit none
      type(runtypedata) :: env
      real(wp), intent(out) :: ieeel !number of valence electrons
      real(wp), intent(out) :: ehomo
      integer, intent(out) :: nat
      integer :: i, j, nvel, chrg
      integer, allocatable :: iat(:)
      type(structure_type) :: mol
      type(error_type), allocatable :: error

      call read_structure(mol, env%infile, error, filetype%xyz)
      if (allocated(error)) then
         if (error%stat .eq. 1) then ! error reading structure
            write (*, *) "ERROR: Could not read input structure"
            error stop
         end if
      end if
      nat = mol%nat

      chrg = env%chrg
      allocate (iat(nat))
      iat = mol%num(mol%id)
      call velectrons_amount(nat, iat, chrg, nvel)
      ! want number of electrons of uncharged species
      nvel = nvel + chrg
      ! number of valence electrons of input molecule
      ieeel = dble(nvel) !env%ieel not sure why this has to be a real number
      ! for now, we take IP for ehomo
      ehomo = env%ehomo

   end subroutine getspecs

! determine automatically the two parameters in the P(E)
! function for given energy range exc, the number of bonds
! nbnd and the input parameter ieeatm (av energy/per bond in eV,
! normally 0.6 eV/bond)

   subroutine getieeab(iee_a, iee_b, ieeel, edistri, exc, nat, ieeatm, pmax) !in previous version nbnd instead of nat ??

      integer  :: edistri, nat, k

      real(wp) :: iee_a, iee_b, ieeel, pmax, ieemax, exc, E_avg, st, ieeatm

      st = 0.005_wp
      iee_a = 0.0_wp
      iee_b = 0.0_wp
      k = 0

      do
         k = k + 1

         iee_a = iee_a + st
         iee_a = min(iee_a, 0.3)
         iee_b = iee_b + st*7

         call getmaxiee(iee_a, iee_b, ieeel, edistri, exc, ieemax, pmax, E_avg)

         if (k > 10000) stop 'internal error inside getieeab'
         if (E_avg/nat >= ieeatm) exit

      end do

   end subroutine getieeab

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! get maximum value (pmax) at which energy (ieemax) and
   ! average energy E_avg for given parameters in the P(E)
   ! distribution function
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! in fact, only getting pmax here is the important part
   subroutine getmaxiee(iee_a, iee_b, ieeel, edistri, exc, ieemax, pmax, E_avg)

      integer  :: edistri
      real(wp) :: iee_a, iee_b, ieeel, val, x, pmax, ieemax, exc, E_avg, m

      x = 0.001_wp
      pmax = -1.0_wp
      ieemax = 0.0_wp
      E_avg = 0.0_wp
      m = 0.0_wp

      do
         if (edistri == 0) call gauss0(iee_a, iee_b, ieeel, x, val)
         if (edistri == 1) call poiss0(iee_a, iee_b, ieeel, x, val)

         if (val > pmax) then
            pmax = val
            ieemax = x
         end if

         x = x + 0.01_wp
         E_avg = E_avg + val*x
         m = m + val

         if (x >= exc) exit
      end do

      E_avg = E_avg/m

   end subroutine getmaxiee

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Gaussian energy distribution, ieeel= # val el, iee_a and iee_b are parameters
   subroutine gauss0(iee_a, iee_b, ieeel, x, p)

      real(wp) :: iee_a, iee_b, ieeel, x, p

      p = exp(-iee_a*(x - ieeel*iee_b)**2/ieeel)

   end subroutine gauss0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Poisson energy distribution, ieeel= # val el, iee_a and iee_b are parameters
   subroutine poiss0(iee_a, iee_b, ieeel, x, t18)

      real(wp) :: iee_a, iee_b
      real(wp) :: k, z, ieeel, x, t18
      real(wp) :: t2, t8, t14, t17
      !eq. (2) in Angewandte paper SI doi.org/10.1002/anie.201300158
      z = iee_b
      k = 1.0_wp/iee_a
      t2 = k/ieeel ! is c
      t8 = dlog(z/k*ieeel/x)
      t14 = EXP(t2*x*(0.1D1 + t8) - 1.D0*z)
      t17 = (t2*x + 0.1D1)**(-0.5_wp) !! Different then in original publication c instead of a in denominator
      ! but this seems right, for large molecules "a"/iee_a instead of "c"/t2 would give too large P(E) for values over 70eV
      t18 = t14*t17
!     write(*,*) "t18 is", t18

   end subroutine poiss0

   ! taken from QCxMS and modified to continous scaling of
   ! ESI energy instead of discrete scaling with number of atoms
   ! TODO test this with different molecules for CID mode
   subroutine scaleesi(env, nsamples, eiee, piee, plotlog)
      use cid
      implicit none
      integer, intent(in) :: nsamples
      real(wp), allocatable, intent(out) :: eiee(:), piee(:)
      real(wp) :: edum, eimpw, randx, escale
      real(wp) :: energy, prob
      logical, intent(in) :: plotlog
      integer :: i, ich
      integer :: nat
      integer :: b, dep
      type(runtypedata) :: env

      allocate (eiee(nsamples), piee(nsamples))
      eiee = 0.0_wp
      piee = 0.0_wp

      if (env%cid_mode == 3) then ! only collisions
         eiee = 0.0_wp
         piee = 1.0_wp/nsamples
         return
      end if

      call rdshort_int(trim(env%infile), nat)

      !scale dependent on the number of atoms
      ! if (env%cid_esi .le. 0.0_wp) then ! if esi is not set set dependent on the number of atoms

      ! taken from QCxMS
      call random_number(randx)

      ! b    = nat / 10.0_wp
      ! dep = int(b)
      !  if ( b < 2.0_wp ) dep = 1 ! hard coded for small molecules

      !>> Make random energy depending on molecular size
      !if (dep == 1) escale = 1 + FLOOR(2*randx)     ! 1-2
      !if (dep == 2) escale = 2 + FLOOR(2*randx)     ! 2-3
      !if (dep == 3) escale = 3 + FLOOR(3*randx)     ! 3-5
      !if (dep == 4) escale = 3 + FLOOR(4*randx)     ! 3-6
      !if (dep >= 5) escale = 4 + FLOOR(5*randx)     ! 4-8
      !else
      !    escale = env%cid_esi
      ! end if

      ! escale = dep
      !  escale = escale + env%cid_esi

      if (env%cid_esiatom == 0) then ! i.e. no input take default
         if (env%cid_mode == 2) then
            env%cid_esiatom = 0.4_wp ! temprun mode
         else
            env%cid_esiatom = 0.1_wp
         end if
      end if

      escale = env%cid_esiatom*nat + env%cid_esi
      write (*, *) "scaling energy in ESI by", env%cid_esiatom, " eV per atom and shift by", env%cid_esi, " eV"
      write (*, *) "ESCALE in ESI is", escale

      ! vary energies by Box-Muller
      eimpw = env%cid_esiw ! width of distribution 0.2 is default
      do i = 1, nsamples
         edum = vary_energies(escale, eimpw)
         eiee(i) = edum
         piee(i) = 1.0_wp/nsamples
      end do

      if (plotlog) then
         open (newunit=ich, file='eimp.dat', status='replace')
         do i = 1, nsamples
            write (ich, *) eiee(i), piee(i)
         end do
         close (ich)
      end if
   end subroutine scaleesi
end module qcxms_iee
