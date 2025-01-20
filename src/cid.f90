
! MODULE for CID mode: highly experimental and not yet fully implemented!!!!
module cid
   use qcxms2_data
   use iomod
   use tsmod
   use utility
   use structools
   use isotope_pattern
   implicit none
contains
   ! ALGO:
!!
!     ! one mcsimu run is:
!     ! 1. compute t between collission for this fragment! for subsequent fragmentations,
!    !    ekin should be elav - sum_dEkin + kerav/(1+m_1/m_2) / 2 (or *m1/m2 ??? both fragments get) ! energy and momentum conservation-> Ekin1 = m2/m1 Ekin2 = m1/m2 Ekin1
!     ! 2. compute shortest tav, if tav > t between collisions, do collision
!     ! 3. > modify eiee
!!
!     ! 4. repeat 2. and 3. until tav < t between collisions, after 3, tcoll gets also updated???
!     ! 5. do mcsimu run
!     next fragmentation: new ekin0, new tcoll0 ! but when do we stop?, need to track v and t throughout the complete process... or we simplify
!       tcoll, only once for every fragment, also v once, and with this we can track, t and v throughout the complete process... not sure if this that much easier...
!!    ??? HOW TO TRACK TF???? tav will get faster than t, after collision, but we obviously dont go back in time ...
!    ! track time of collision, and then compute tav for next fragment, and then do collision, and then compute tav for next fragment, and so on
!      and via time and v and lchamb we now, when to stop
!     ! 1.:
!     ! 2.:
!     ! need collision counter, no collisions after max number of collisions
!     ! how to compute fastes tav? have to save barriers once, and then
!     !    can do this
!     !    just take lowest barrier, and then compute it for this

! simulate CID mode by modifiying the energy distribution
   subroutine simcid(env, tf, maxiee, eiee, piee, nsamples, alldgs, nincr, nfragl, KERavold,&
      & de0, Ea0, nat0, nvib0, npairs, ishdiss, scaleeinthdiss, isrearr, noisos)
      use mctc_env, only: error_type, get_argument!, fatal_error
      use mctc_io, only: structure_type, read_structure, write_structure, &
         & filetype, get_filetype, to_symbol, to_number
      use xtb_mctc_convert, only: amutokg
      implicit none
      character(len=80) :: fname
      real(wp), intent(in) :: maxiee
      real(wp), intent(in) :: tf
      type(runtypedata) :: env
      type(structure_type) :: mol
      type(error_type), allocatable :: error
      integer, intent(in) :: nfragl
      integer, intent(in) :: nsamples
      integer, intent(in) :: nincr ! number of increments for the calculation of the thermo contributions to the barriers
      integer, allocatable :: iat(:)
      real(wp), intent(inout) :: eiee(nsamples)
      real(wp), intent(in) :: piee(nsamples)
      real(wp), intent(in), allocatable :: alldgs(:, :)
      integer :: nat
      integer :: nvib
      integer :: ncoll ! number of performed collisions for this ion
      integer :: maxcoll ! average number of collisions for this ion
      integer :: ind ! energy index for energy increment
      integer :: i, j
      integer :: ifast ! index of the fastest reaction
      real(wp) :: mp ! mass of ion
      real(wp) :: dEint ! internal energy change upon collision
      real(wp) :: dEkin ! kinetic energy change upon collision
      real(wp) :: minbarrier ! minimal barrier for the fastest reaction
      real(wp) :: Ekin ! kinetic energy of the ion
      real(wp) :: sum_dEint ! sum of kinetic energy changes upon collisions
      real(wp) :: sum_dEkin ! sum of internal energy changes upon collisions
      real(wp) :: x_trav ! distance traveled by the ion
      real(wp) :: qmass
      real(wp) :: calc_coll
      real(wp) :: mfp, kmax, tmin, tcoll
      real(wp) :: t_trav
      real(wp) ::  eavg, k, E
      real(wp) ::  prem ! survival yield of ion
      ! FOR KER computation
      real(wp) :: KER, KER_Eex, KER_Ear, KERav
      real(wp), intent(in) :: KERavold, de0, Ea0 ! KER, reaction energy and barrier of prior fragmentations
      integer, intent(in) :: nat0, nvib0, npairs
      real(wp), intent(in) :: scaleeinthdiss
      logical, intent(in) :: ishdiss(npairs)
      logical, intent(in) :: noisos ! if only fragmentations are considered
      logical, intent(in) :: isrearr(npairs)
      real(wp), parameter :: evtoJ = 1.602176634e-19_wp ! eV to Joule
      character(len=80) :: fcoll ! for plotting the energy distribution

      logical :: ex

      ! TODO FIXME, i.e. add freq here and compute kmax for the fastest reaction
      if (.not. env%eyring) then
         write (*, *) 'Error: CID mode only works currently only with Eyring enabled'
         stop
      end if
      write (*, *) "CID mode", env%cid_mode
      if (env%cid_mode == 2) then
         write (*, *) " CID Temprun mode, no collisions are performed"
         return
      end if

      inquire (file="fragment.xyz", exist=ex)
      if (ex) then
         fname = 'fragment.xyz'
      end if
      inquire (file="isomer.xyz", exist=ex)
      if (ex) then
         fname = 'isomer.xyz'
      end if
      call read_structure(mol, fname, error, filetype%xyz)
      allocate (iat(mol%nat))
      iat = mol%num(mol%id)
      nat = mol%nat

      call getmolmass(nat, iat, mp) ! TODO rewrite this function need xyz also for collision setup
      if (mp .le. 0.0_wp .or. mp .ge. 1e+5_wp) then
         write (*, *) 'Error: Could not determine mass of ion'
         stop
      end if

      nvib = 3*nat - 6 ! TODO correct for linear molecules

      ! TODO
      ! if ion is fragment of ion, we need to know mass of the other ion
      ! we need to know the barriers for the reactions

      ! read in data from prior fragmentations

      ! TODO we need this EKIN info already for the collision setup to compute tcoll
      if (nfragl == 1 .and. .not. noisos) then ! after isomers, we already have sumdeint ...
         sum_dEkin = 0.0_wp
         x_trav = 0.0_wp
         t_trav = 0.0_wp
         Ekin = env%cid_elab ! first collision first fragmentatio
         sum_dEint = 0.0_wp
      else
         call rdshort('sumdekin', sum_dEkin) ! read in from Dir from fragment
         call rdshort('sumdeint', sum_dEint) ! energy from previous collisions
         call rdshort('x_trav', x_trav)    ! distance traveled by the ion
         write (*, *) "read x_trav is ", x_trav
         call rdshort('qmass', qmass) ! mass quotient of two fragments
         call rdshort('t_trav', t_trav) ! time in collision chamber just for testing
         ! first collision subsequent fragmentations
         ! the fragment ion has to be faster due to the KER from the fragmentation
         ! energy and momentum conservation gives this energy partioning of KER

         ! correct in principle but wrong results with it
         Ekin = (env%cid_elab + sum_dEkin + KERavold)*1/(1 + qmass) ! for isomers qmass = 0 is correct
         ! write(*,*) "weight energy by qmass", 1/ (1+ qmass)

         !Ekin =  env%cid_elab + sum_dEkin + KERavold * 1/ (1+ qmass)
         !
         if (Ekin .lt. 0.0_wp) then
            write (*, *) "Ekin is negative, something went wrong, set it to 0"
            Ekin = 0.0_wp
         end if
      end if
      write (*, *) "Ekin is ", Ekin

      call collision_setup(env, nat, iat, mol%xyz, mfp, calc_coll)
      ncoll = 0
      maxcoll = calc_coll

   !! very important to consider this from prior fragmentations
      if (.not. noisos) eiee = eiee + sum_dEint
      ! TODO modify this with BOXMULLER???
      do i = 1, maxcoll

         if (x_trav .ge. env%cid_lchamb) then
            write (*, *)
            write (*, *) "Ion has traveled the chamber length, no collissions are possible"
            ! write down the data for the next fragmentation
            call wrshort_real('sumdekin', sum_dEkin) ! write out sum_dEkin
            call wrshort_real('sumdeint', sum_dEint)
            call wrshort_real('x_trav', x_trav) ! write out x_trav
            call wrshort_real('t_trav', t_trav)
            return ! stop if ion has traveled the chamber length
         end if
         ! compute half life time for fastest reaction at Eavg
         call calc_Eavg(nsamples, eiee, piee, Eavg)

         if (nfragl .gt. 1 .and. nat0 .gt. 0 .and. env%calcKER) then
            ! we compute here the KER of the last reaction which should be missing here
            ! for both fragments
            !TODO if we scale the energy of H-dissociations we should do this here too??? wrong if not !!!! TODO
            if (nat0 .eq. 2) then ! diatomic molecule
               KER_Eex = (Eavg - Ea0)
               write (*, *) "this case shouldnt exist"
            else
               KER_Eex = (Eavg - Ea0)/(0.44_wp*nvib0) ! J. Chem. Phys. 48, 4093 (1968)
            end if
            ! nvib0 vibrational degrees of freedom of transition state of last reaction
            KER_Ear = (0.33_wp*(Ea0 - de0)) ! J. Chem. Phys. 61, 1305 (1974)
            ! no negative kinetic energies allowed
            if (KER_Eex .lt. 0.0_wp) KER_Eex = 0.0_wp
            if (KER_Ear .lt. 0.0_wp) KER_Ear = 0.0_wp
            ! add here average KER of prior fragmentations
            KER = KER_Eex + KER_Ear + KERavold
       !!TODO KERexpl KER =  KER_Eex + KER_Ear + KERold(i)
            if (KER .lt. 0.0_wp) KER = 0.0_wp ! can sometimes happen

            if (env%scaleker .ne. 1.0_wp) KER = KER*env%scaleker
            if (env%printlevel .eq. 3) write (*, *) "KER, KER_Eex, KER_EAR are:", KER, KER_Eex, KER_EAR
            !Kerexpl TODO  write(ich,*)    KER
         else
            KER = 0.0_wp

         end if
         Eavg = Eavg - KER ! TODO check if this is correct, should be correct, but check again
         Ekin = Ekin + KER*1/(1 + qmass)
         ! compute the fastest reaction at Eavg

         ind = getindex(nincr, maxiee, Eavg)
         !write(*,*) "ind cid is", ind, nincr, maxiee, Eavg
         kmax = 0.0_wp
         ifast = -1
         do j = 1, npairs
            if (noisos .and. isrearr(j)) cycle
            if (.not. noisos .and. .not. isrearr(j)) cycle
            if (ishdiss(j)) then
               E = Eavg*scaleeinthdiss
            else
               E = Eavg
            end if
            k = calc_eyring(E, alldgs(j, ind), nvib)
            if (k .gt. kmax) then
               kmax = k
               ifast = j
            end if
         end do
         ! tmin = log(2.0_wp) / kmax
         tmin = log(2.0_wp)/kmax ! 50 percent of the fastest reaction remain
         tcoll = mfp*sqrt(mp*amutokg/(2.0_wp*Ekin*evtoJ))
         write (*, *) "tcoll,mfp, Ekin*evtoJ is", tcoll, mfp, mp, Ekin*evtoJ
         write (*, *) "Time between collisions is ", tcoll, "s" ! usually mikroseconds
         write (*, *) "Shortest tav is ", tmin, "s for reaction ", ifast, alldgs(ifast, ind), " at energy ", Eavg
         prem = EXP(-kmax*tf) ! tf contains sum of all tavs is tf remaining
         write (*, *) "prem is ", prem
         ! if fastest reaction is slower than time between collisions, do collision
         if (tmin .gt. tcoll) then

            ncoll = ncoll + 1
            ! do collision by modifying eiee and computing dEkin
            call simcoll(env, eiee, nsamples, mp, Ekin, dEkin, dEint)
            write (*, *) "perform collision", i, "with ekin, dekin, deint", Ekin, dEkin, dEint
            sum_dEkin = sum_dEkin + dEkin
            sum_dEint = sum_dEint + dEint
            !if (.not. env%cid_nocooling)

            Ekin = Ekin + dEkin*env%cid_scool
            write (*, *) "kinetic energy loss by ", dEkin*env%cid_scool
            write (*, *) "kinetic energy is ", Ekin, ",scaling of cooling is", env%cid_scool
            !Ekin = Ekin
            x_trav = x_trav + mfp*i !
            write (*, *) "x_trav is ", x_trav
            t_trav = t_trav + tcoll
            if (.not. noisos) write (fcoll, '(a,i0)') 'ecolliso_', i
            if (noisos) write (fcoll, '(a,i0)') 'ecollfrag_', i
            call printiee(nsamples, eiee, piee, fcoll)
         else
            write (*, *) "fastest reaction is faster than time between collisions, perform fragmentation instead"
            exit ! stop if fastest reaction is faster than time between collisions
         end if

      end do

      !write(*,*) "Performed on average", ncoll, " collisions"
      write (*, *) "Performed ", ncoll, " collisions"

      ! write down the data for the next fragmentation
      call wrshort_real('sumdekin', sum_dEkin) ! write out sum_dEkin
      call wrshort_real('sumdeint', sum_dEint)
      call wrshort_real('x_trav', x_trav) ! write out x_trav
      call wrshort_real('t_trav', t_trav)

      ! TODO print out energy distribution
      ! open (newunit=ich, file='eimp2.dat', status='replace')
      !    do i = 1, nsamples
      !      write (ich, *) eiee(i), piee(i)
      !  end do
      !  close (ich)
   end subroutine simcid
! formulas taken from: J. Am. Soc. Mass Spectrom. 2013, 24, 7, 1064–1071
   ! DOI: https://pubs.acs.org/doi/full/10.1007/s13361-013-0635-8

   ! simulate collision
   ! modifying the energy distribution
   subroutine simcoll(env, eiee, nsamples, mp, Ekin, dEkin, dEint)
      use qcxms_boxmuller
      implicit none
      type(runtypedata) :: env
      real(wp), intent(inout) :: eiee(nsamples)
      integer, intent(in) :: nsamples
      real(wp), intent(out) :: dEkin ! kinetic energy change upon collision
      real(wp), intent(out) :: dEint ! internal energy change upon collision
      real(wp), intent(in) :: ekin
      real(wp) :: mg ! mass of gas and ion
      real(wp), intent(in) :: mp
      real(wp) :: beta, gamma
      real(wp) :: ecom
      real(wp) :: eta !collision inelasticity

      real(wp), parameter :: evtoJ = 1.602176634e-19_wp ! eV to Joule
      real(wp) :: edist ! width of variation of number of collisions
      real(wp) :: edum ! dummy variable
      integer :: i, ncoll

      eta = 0.43_wp ! TODO tune this  value
      beta = env%cid_mgas/(env%cid_mgas + mp)

      ! Ecom is maximal amount of energy that can be transferred to the molecule
      !Ekin is input parameter (ELAB) subsequent EKIns are computed  with Ekin = Elab - deEINT
      Ecom = beta*Ekin

      dEint = eta*Ecom

      edist = env%cid_collw ! width of distribution 0.5 is default
      ! vary collision energy by boxmuller Distribution
      do i = 1, nsamples
         if (edist .gt. 0.0) then
            edum = vary_energies(dEint, edist)
         else
            edum = dEint
         end if
         if (edum .lt. 0.0_wp) edum = 0.0_wp
         eiee(i) = eiee(i) + edum
      end do

      ! do i = 1, nsamples
      !    ! displace energy distribution a little bit towards higher values
      !    ncoll = vary_collisions(1.0_wp,edist) ! for each nsample we vary the number of collisions
      !    !write(*,*) "Number of collisions is i", ncoll
      !    eiee(i) = eiee(i) + dEint * ncoll
      !
      ! end do

      gamma = (2*mp + (eta*env%cid_mgas))/(mp + env%cid_mgas)

      if (gamma*beta .ge. 1.0_wp) then
         write (*, *) "Error: gamma * beta is larger than 1, this is not possible"
      end if

      dEkin = -gamma*Ecom

      ! we assume elab/Ekin is the same for all samples, randomness already included in the esi distribution and by vary collisions

   end subroutine simcoll


!! Taken and slightly adapted from QCxMS !
!!#########################################################################
!! Calculate the molecular radius, cross section and mean-free-path
!!#########################################################################
   subroutine collision_setup(env, nat, iat, xyz, mfpath, calc_coll)
      use xtb_mctc_convert, only: autoaa
      use xtb_mctc_constants, only: pi
      implicit none
      type(runtypedata) :: env
      real(wp), intent(in)  :: xyz(3, nat)
      integer, intent(in)  :: iat(nat)
      real(wp) :: r_mol, cross
      real(wp), intent(out) :: mfpath
      real(wp), intent(out) :: calc_coll ! average number of collisions for this ion
      integer  :: nat
      integer  :: i
      real(wp) :: cg(3), mol_rad, rtot
      real(wp) :: r_atom
      real(wp), parameter :: kB_J = 1.38064852E-23

      ! TODO move this to a separate module for nicer code structure ..
      ! Radius used in QCxMS (in au)

      real(wp), parameter :: Rad(118) = aatoau *  [ &
      & 0.32_wp,0.37_wp, & ! H,He
      & 1.30_wp,0.99_wp,0.84_wp,0.75_wp,0.71_wp,0.64_wp,0.60_wp,0.62_wp, & ! Li-Ne
      & 1.60_wp,1.40_wp,1.24_wp,1.14_wp,1.09_wp,1.04_wp,1.00_wp,1.01_wp, & ! Na-Ar
      & 2.00_wp,1.74_wp, & ! K,Ca
      & 1.59_wp,1.48_wp,1.44_wp,1.30_wp,1.29_wp, & ! Sc-
      &  1.24_wp,1.18_wp,1.17_wp,1.22_wp,1.20_wp, & ! -Zn
      &  1.23_wp,1.20_wp,1.20_wp,1.18_wp,1.17_wp,1.16_wp, & ! Ga-Kr
      & 2.15_wp,1.90_wp, & ! Rb,Sr
      &  1.76_wp,1.64_wp,1.56_wp,1.46_wp,1.38_wp, & ! Y-
      &  1.36_wp,1.34_wp,1.30_wp,1.36_wp,1.40_wp, & ! -Cd
      &   1.42_wp,1.40_wp,1.40_wp,1.37_wp,1.36_wp,1.36_wp, & ! In-Xe
      & 2.38_wp,2.06_wp, & ! Cs,Ba
      & 1.94_wp,1.84_wp,1.90_wp,1.88_wp,1.86_wp,1.85_wp,1.83_wp, & ! La-Eu
      & 1.82_wp,1.81_wp,1.80_wp,1.79_wp,1.77_wp,1.77_wp,1.78_wp, & ! Gd-Yb
      &  1.74_wp,1.64_wp,1.58_wp,1.50_wp,1.41_wp, & ! Lu-
      &  1.36_wp,1.32_wp,1.30_wp,1.30_wp,1.32_wp, & ! -Hg
      &  1.44_wp,1.45_wp,1.50_wp,1.42_wp,1.48_wp,1.46_wp, & ! Tl-Rn
      & 2.42_wp,2.11_wp, & ! Fr,Ra
      &  2.01_wp,1.90_wp,1.84_wp,1.83_wp,1.80_wp,1.80_wp,& ! Ac-Pu
      ! from covalent 2009 covalent radii, such that it is complete up to 118
      &  1.49_wp, & ! Am
      &   1.49_wp,1.51_wp,1.51_wp,1.48_wp,1.50_wp,1.56_wp,1.58_wp, & ! Cm-No
      &   1.45_wp,1.41_wp,1.34_wp,1.29_wp,1.27_wp, & ! Lr-
      &  1.21_wp,1.16_wp,1.15_wp,1.09_wp,1.22_wp, & ! -Cn
      &   1.36_wp,1.43_wp,1.46_wp,1.58_wp,1.48_wp,1.57_wp ] ! Nh-Og   


      mol_rad = 0
      rtot = 0

      call center_of_geometry(nat, xyz, cg)

      ! Determine the radii by checking the largest distance from the center of mass
      do i = 1, nat
         mol_rad = sqrt((xyz(1, i) - cg(1))**2 + (xyz(2, i) - cg(2))**2 + (xyz(3, i) - cg(3))**2) + Rad(iat(i))
         if (mol_rad > rtot) rtot = mol_rad
      end do

      ! Calculate radii in meters
      r_atom = (env%cid_rgas*autoaa)*1E-10 ! m
      r_mol = (rtot*autoaa)*1E-10 ! m

      !write(*,*) 'R TOT', rtot
      !write(*,*) 'Radius Molecule', r_mol

      ! calculate collisional cross section and mean free path
      cross = pi*((r_mol + r_atom)**2)

      mfpath = (kB_J*env%cid_Tgas)/(cross*env%cid_Pgas)
      write (*, *) kB_J, env%cid_Tgas, cross, env%cid_Pgas
      write (*, *) "mean free path is ", mfpath, "meters"

      calc_coll = env%cid_lchamb/mfpath

      ! write(*,*) "calc_coll is", calc_coll

   end subroutine collision_setup

   ! TODO Maybe we want to add this routine later 
!   subroutine protonate(env)
!      implicit none
!      character(len=1024) :: jobcall
!      integer :: io, i
!      type(runtypedata) :: env
!      type(timer):: tim
!
!      !call crest protonate mode
!      write (jobcall, '(a)') 'crest '//trim(env%infile)//' --protonate --ewin 50 > protonate.out 2>/dev/null'
!      !DEL write(*,*) jobcall
!      !call tim%start(1,'Protonation')
!      call execute_command_line(trim(jobcall), exitstat=io)
!
!      !call eval_timer(tim)
!      !create directories for protonated structures
!      write (jobcall, '(a)') 'crest -splitfile protonated.xyz 1>/dev/null 2>/dev/null'
!      call execute_command_line(trim(jobcall), exitstat=io)
!
!      call findprotts(env)
!   end subroutine protonate
!   subroutine findprotts(env)
!      character(len=1024) :: jobcall
!      character(len=24) :: DIR
!      integer :: io, i, nat, nprot
!      logical :: ex
!      type(runtypedata) :: env
!      call readstruc('protonated.xyz', nat, nprot)
!      call chdir('SPLIT')
!      do i = 2, nprot
!         write (DIR, "(A,I1)") 'STRUC', i
!         call copy('STRUC1/struc.xyz', trim(DIR)//'/start.xyz')
!         call chdir(DIR)
!         call move('struc.xyz', 'end.xyz')
!         !call tssearch(env)
!         call chdir("..")
!      end do
!   end subroutine findprotts
!
!
end module cid
