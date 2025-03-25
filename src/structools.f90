module structools
   use xtb_mctc_accuracy, only: wp
   use iomod
   use qcxms2_data
   use utility
   implicit none
contains

   ! subroutine to detect topo changes
   ! based on wbo calculation of xtb
   ! for now only to detect changes in R-H bonds for scaling
   ! only for rearrangments important for now
   ! TODO maybe extend later so that we can perform hybrid NEB transition state search
   ! or use active atoms option in TS optimization

   subroutine readrhbond(fname, rhbond)
      use mctc_env, only: error_type, get_argument!, fatal_error
      use mctc_io, only: structure_type, read_structure, write_structure, &
         & filetype, get_filetype, to_symbol, to_number
      implicit none
      ! saves bond partner with highest wbo for all H, 0 if atom index no H
      integer, allocatable, intent(out) :: rhbond(:)
      type(runtypedata) :: env
      integer :: nat
      integer :: nbonds
      integer, allocatable :: iat(:)
      character(len=*), intent(in) :: fname
      type(structure_type) :: mol
      type(error_type), allocatable :: error
      integer :: i, j, k
      integer :: ich, io
      real(wp) :: wbo
      character(len=80) :: tmp

      call read_structure(mol, fname, error, filetype%xyz)
      nat = mol%nat
      allocate (iat(nat))
      iat = mol%num(mol%id)

      open (newunit=ich, file='wbo', iostat=io)
      nbonds = 0
      do
         read (ich, '(a)', iostat=io) tmp
         if (io < 0) exit
         nbonds = nbonds + 1
      end do
      close (ich)
      allocate (rhbond(nat))

      rhbond = 0
      open (newunit=ich, file='wbo', iostat=io)
      do i = 1, nbonds  ! TODO CHECK THIS; always the case? even for heavy elements? no its not
         if (io < 0) exit
         read (ich, *) j, k, wbo
         ! TODO tune this threshold
         if (wbo .gt. 0.3_wp) then  ! strange 2-pent of gxtb geometry has 0.47 wbo for C-H bond
         if (iat(j) .eq. 1) then
            rhbond(j) = k
         elseif (iat(k) .eq. 1) then
            rhbond(k) = j
         end if
         end if
      end do
      close (ich)

  !!do i = 1, nat
  !!   write(*,*) "hbond partner is: ",i,rhbond(i)
  !!end do

   end subroutine readrhbond

!input two vectors with rhbond information of start fragment and isomer
   subroutine compare_rhbond(nat, rhbond0, rhbond, ishdiss)

      implicit none
      integer, intent(in) :: nat
      integer, intent(in) :: rhbond0(nat), rhbond(nat)
      logical, intent(out) :: ishdiss
      integer :: i
      ishdiss = .false.

      do i = 1, nat
         if (rhbond(i) .ne. rhbond0(i)) then
            !write(*,*) "R-H bond has changed for nat", i
            ishdiss = .true.
         end if
      end do

   end subroutine compare_rhbond

! only if H or H2 is dissociated
   subroutine findhdiss(env, fname, nat, npairs, fragdirs, ishdiss, scaleeinthdiss)
      implicit none
      type(runtypedata) :: env
      character(len=80), intent(in) :: fragdirs(npairs, 3)
      integer, intent(in) :: npairs
      integer, intent(in) :: nat
      logical, allocatable, intent(out) :: ishdiss(:)
      integer :: i, j
      logical :: ex
      character(len=256) :: jobcall
      character(len=1024) :: thisdir
      character(len=80) :: fname
      character(len=80) :: sumform
      real(wp) :: tav0 !half life of reactions up to this point in ps
      real(wp) :: time_relax !TODO tune this parameter, maybe scale with size of vibrational DOF
      real(wp), intent(out) :: scaleeinthdiss
      integer, allocatable :: rhbond0(:), rhbond(:)

      !  scaleeinthdiss = env%scaleeinthdiss*1/(1+exp(-tav0*4.0_wp/time_relax)) ! alternative logistic scaling
      ! scaleeinthdiss = env%scaleeinthdiss*exp(tav0/time_relax*LOG(2)) ! exponential scaling
      ! tune this parameter estimated time of relaxation, dependent on number of vibrational degrees of freedom of system
      !if (env%scaleeinthdiss .ne. 1.0_wp)
      !scaleeinthdiss = 1.0_wp - env%scaleeinthdiss*EXP(-(tav0/time_relax)) ! exponential decay of scaling factor ! TODO formula has to be improved
      !scaleeinthdiss = env%scaleeinthdiss
      call getcwd(thisdir)
      !write(*,*) "scaleeinthdiss is", scaleeinthdiss, env%scaleeinthdiss, tav0, time_relax
      tav0 = 0.0_wp
      if (env%printlevel .eq. 3) write (*, *) "tav0 is", tav0, " s"
      call rdshort_real('tav', tav0)

      if (tav0 .lt. 0.0_wp) then
         write (*, *) "tav is negative, something went wrong, set it to 0"
         tav0 = 0.0_wp
      end if
      if (tav0 .gt. 50.0e-03_wp) then
         write (*, *) "tav is much to large, something went wrong, set it to 0"
         tav0 = 0.0_wp
      end if

      time_relax = env%erelaxtime
      scaleeinthdiss = 1.0_wp
      !if (tav0 .lt. time_relax) then ! a few picoseconds
      !   scaleeinthdiss = env%scaleeinthdiss*(1+(tav0/time_relax)) ! linear scaling
      !end if
      scaleeinthdiss = env%scaleeinthdiss
      ! write(*,*) "Internal energy of H-rearrangment scaled with ", scaleeinthdiss
      ! scale frequency for H-abstraction or H-rearrangment to account for inishomogenous IEE distribution in molecule

      allocate (ishdiss(npairs))
      ishdiss = .false.

      ! TODO rewrite this
      if (env%cneintscale) return

      do i = 1, npairs
         !check if  fragment  is isomer
         if (index(fragdirs(i, 3), 'p') .ne. 0) then
            do j = 2, 3
               call getsumform(trim(env%path)//"/"//trim(fragdirs(i, j))//"/fragment.xyz", sumform)
               if (sumform == "H1" .or. sumform == "H2") then
                  call touch(trim(env%path)//"/"//trim(fragdirs(i, j))//"/ishdiss")
                  call touch(trim(env%path)//"/"//trim(fragdirs(i, 1))//"/ishdiss")
                  ishdiss(i) = .true.
               end if
            end do
            if (ishdiss(i)) then
               write (*, '(a,a,f4.2)') "Internal energy of H-diss to fragpair ", &
               & trim(fragdirs(i, 1))//" scaled with ", scaleeinthdiss
            end if
         end if
      end do

   end subroutine findhdiss

! based on wbo calculation of xtb also for rearrangments
   subroutine findhdisswbo(env, fname, nat, npairs, fragdirs, ishdiss, scaleeinthdiss)
      implicit none
      type(runtypedata) :: env
      character(len=80), intent(in) :: fragdirs(npairs, 3)
      integer, intent(in) :: npairs
      integer, intent(in) :: nat
      logical, allocatable, intent(out) :: ishdiss(:)
      integer :: i, j
      logical :: ex
      character(len=256) :: jobcall
      character(len=1024) :: thisdir
      character(len=80) :: fname
      real(wp) :: tav0 !half life of reactions up to this point in ps
      real(wp) :: time_relax !TODO tune this parameter, maybe scale with size of vibrational DOF
      real(wp), intent(out) :: scaleeinthdiss
      integer, allocatable :: rhbond0(:), rhbond(:)

      !  scaleeinthdiss = env%scaleeinthdiss*1/(1+exp(-tav0*4.0_wp/time_relax)) ! alternative logistic scaling
      ! scaleeinthdiss = env%scaleeinthdiss*exp(tav0/time_relax*LOG(2)) ! exponential scaling
      ! tune this parameter estimated time of relaxation, dependent on number of vibrational degrees of freedom of system
      !if (env%scaleeinthdiss .ne. 1.0_wp)
      !scaleeinthdiss = 1.0_wp - env%scaleeinthdiss*EXP(-(tav0/time_relax)) ! exponential decay of scaling factor ! TODO formula is bullshit
      !scaleeinthdiss = env%scaleeinthdiss
      call getcwd(thisdir)
      !write(*,*) "scaleeinthdiss is", scaleeinthdiss, env%scaleeinthdiss, tav0, time_relax
      tav0 = 0.0_wp
      call rdshort_real('tav', tav0)

      write (*, *) "tav0 is", tav0
      if (env%printlevel .eq. 3) then
         if (tav0 .lt. 0.0_wp) then
            write (*, *) "tav is negative, something went wrong, set it to 0 s"
            tav0 = 0.0_wp
         end if
         if (tav0 .gt. 50.0e-03_wp) then
            write (*, *) "tav is much to large, something went wrong, set it to 0 s"
            tav0 = 0.0_wp
         end if
      end if

      time_relax = env%erelaxtime
      scaleeinthdiss = 1.0_wp
      if (tav0 .lt. time_relax) then ! a few picoseconds
         scaleeinthdiss = env%scaleeinthdiss*(1 + (tav0/time_relax)) ! linear scaling
      end if
      write (*, '(a,f4.2)') "Internal energy of H-rearrangment scaled with ", scaleeinthdiss
      ! scale frequency for H-abstraction or H-rearrangment to account for inishomogenous IEE distribution in molecule

      allocate (ishdiss(npairs))
      ishdiss = .false.
      inquire (file='wbo', exist=ex)
      if (.not. ex) then
         write (jobcall, '(a,i0,a)') "xtb "//trim(fname)//" --wbo --sp --gfn2 --chrg ", env%chrg, " > wbo.out 2>/dev/null"
         write (*, *) "jobcall is ", trim(jobcall)
         call printpwd()
         call execute_command_line(trim(jobcall))
      end if
      call readrhbond(fname, rhbond0)
      do i = 1, npairs
         !check if  fragment  is isomer
         if (index(fragdirs(i, 3), 'p') .ne. 0) then
            call chdir(trim(env%path)//"/"//trim(fragdirs(i, 1)))
            inquire (file='wbo', exist=ex)
            if (.not. ex) then
               write (jobcall, '(a,i0,a)') "xtb pair.xyz --wbo --sp --gfn2 --chrg ", env%chrg, " > wbo.out 2>/dev/null"
               call execute_command_line(trim(jobcall))
               ! write(*,*) "jobcall is ", trim(jobcall)
               !call printpwd()
            end if
            call readrhbond('pair.xyz', rhbond)
            call compare_rhbond(nat, rhbond0, rhbond, ishdiss(i))
            call chdir(trim(thisdir))
            if (ishdiss(i)) then
            if (scaleeinthdiss .ne. 1.0_wp) then
               write (*, '(a,a,f4.2)') "Internal energy of H-rearrangmet to fragpair ",&
                & trim(env%path)//"/"//trim(fragdirs(i, 1))//" scaled with ", scaleeinthdiss
            end if
            do j = 2, 3
               call touch(trim(env%path)//"/"//trim(fragdirs(i, j))//"/ishdiss")
            end do
            end if
         else
            call chdir(trim(env%path)//"/"//trim(fragdirs(i, 1)))
            inquire (file='wbo', exist=ex)
            if (.not. ex) then
               write (jobcall, '(a,i0,a)') "xtb isomer.xyz --wbo --sp --gfn2 --chrg ", env%chrg, " > wbo.out 2>/dev/null"
               call execute_command_line(trim(jobcall))
               ! write(*,*) "jobcall is ", trim(jobcall)
               !call printpwd()
            end if
            call readrhbond('isomer.xyz', rhbond)
            call compare_rhbond(nat, rhbond0, rhbond, ishdiss(i))
            if (scaleeinthdiss .ne. 1.0_wp) then
               write (*, '(a,a,f4.2)') "Internal energy of H-rearrangmet to isomer ",&
                & trim(env%path)//"/"//trim(fragdirs(i, 1))//" scaled with ", scaleeinthdiss
            end if
            if (ishdiss(i)) call touch("ishdiss")
            call chdir(trim(thisdir))
         end if
      end do

   end subroutine findhdisswbo
! determine effective CN of bond breakage for each pair
   subroutine get_active_cn(env, fname, nat, npairs, fragdirs, active_cn)
      implicit none
      type(runtypedata) :: env
      character(len=80), intent(in) :: fragdirs(npairs, 3)
      integer, intent(in) :: npairs
      integer, intent(in) :: nat
      integer, allocatable :: actatoms(:) ! atom indices of active atoms
      real(wp), intent(out) :: active_cn(npairs) ! CN of active atoms for each pair
      integer :: i, j
      integer :: n_act
      logical :: ex
      character(len=256) :: jobcall
      character(len=1024) :: thisdir
      character(len=80) :: fname
      real(wp) :: tav0 !half life of reactions up to this point in ps
      real(wp) :: time_relax !TODO tune this parameter, maybe scale with size of vibrational DOF
      real(wp), allocatable :: cn0(:), cn(:), bond0(:, :), bond(:, :)
      integer, allocatable :: rhbond0(:), rhbond(:)

      call getcwd(thisdir)

      ! just try CN numbers for now no WBO because more general if other level than GFN2 is used for geometries
      !inquire (file='wbo', exist=ex)
      !if (.not. ex) then
      !   write(jobcall,'(a,i0,a)') "xtb "//trim(fname)//" --wbo --sp --gfn2 --chrg ",env%chrg," > wbo.out 2>/dev/null"
      !   write(*,*) "jobcall is ", trim(jobcall)
      !   call printpwd()
      !   call execute_command_line(trim(jobcall))
      !end if

      ! get CN numbers of reactant and bond order matrix
      call get_CN(fname, cn0, bond0, .false.)

      do i = 1, npairs
         !check if  fragment  is isomer
         if (index(fragdirs(i, 3), 'p') .ne. 0) then
            write (*, *) "determine active atoms for pair ", i
            fname = "pair.xyz"
            call chdir(trim(env%path)//"/"//trim(fragdirs(i, 1)))
            call get_CN(fname, cn, bond, .false.)

            ! compare bond order matrices and get active atoms
            call compare_bond(nat, bond0, bond, actatoms)

            deallocate (bond)
            call get_CN(fname, cn, bond, .true.)
            ! effective CN of bond breakage
            active_cn(i) = 0
            n_act = 0
            do j = 1, size(actatoms)
               if (actatoms(j) .ne. 0) then
                  n_act = n_act + 1
                  write (*, *) "active atom is ", actatoms(j)
                  write (*, *) "CN of active atom is ", cn0(actatoms(j))
                  active_cn(i) = active_cn(i) + cn0(actatoms(j))
               end if
            end do
            active_cn(i) = active_cn(i)/n_act
            write (*, *) "effective CN of bond breakage ", i, " is ", active_cn(i)

            call chdir(trim(thisdir))

            ! for isomers TODO
         else
            fname = "isomer.xyz"
            write (*, *) "determine active atoms for isomer ", i
            call chdir(trim(env%path)//"/"//trim(fragdirs(i, 1)))
            call get_CN(fname, cn, bond, .false.)

            ! compare bond order matrices and get active atoms
            call compare_bond(nat, bond0, bond, actatoms)
            deallocate (bond)
            call get_CN(fname, cn, bond, .true.)

            ! effective CN of bond breakage
            active_cn(i) = 0
            n_act = 0
            do j = 1, size(actatoms)
               if (actatoms(j) .ne. 0) then
                  n_act = n_act + 1
                  write (*, *) "active atom is ", actatoms(j)
                  write (*, *) "CN of active atom is ", cn0(actatoms(j))
                  active_cn(i) = active_cn(i) + cn0(actatoms(j))
               end if
            end do
            active_cn(i) = active_cn(i)/n_act
            write (*, *) "effective CN of bond breakage ", i, " is ", active_cn(i)

            call chdir(trim(thisdir))
         end if
      end do

   end subroutine get_active_cn

   ! execute this only in "ts" directory!!!
   subroutine get_active_cnstring(env, act_atom_string)
      implicit none
      type(runtypedata) :: env
      character(len=80), intent(out) :: act_atom_string
      integer, allocatable :: actatoms(:) ! atom indices of active atoms
      integer :: i, nat
      logical :: ex
      character(len=1024) :: thisdir
      character(len=80) :: fname
      real(wp), allocatable :: cn0(:), cn(:), bond0(:, :), bond(:, :)
      integer, allocatable :: rhbond0(:), rhbond(:)

      act_atom_string = ""
      call getcwd(thisdir)
      call chdir("..") ! go in ts search directory
      call rdshort_int('start.xyz', nat)
      fname = "start.xyz"
      call get_CN(fname, cn0, bond0, .false.)
      fname = "end.xyz"
      call get_CN(fname, cn, bond, .false.)
      call compare_bond(nat, bond0, bond, actatoms)
      do i = 1, size(actatoms)
         actatoms(i) = actatoms(i) - 1 ! ORCA indices start with 0 ...
         if (actatoms(i) .gt. -1) then
            write(act_atom_string, '(a,x,i0)') trim(act_atom_string), actatoms(i)
         end if
      end do
      write(*,*) "active atoms are ", trim(act_atom_string)
      call chdir(trim(thisdir)) ! go back in "ts" dir

   end subroutine get_active_cnstring


! compare two bond order matrices and get active atoms of reaction
   subroutine compare_bond(nat, bond0, bond, actatoms)
      implicit none
      integer, intent(in) :: nat
      real(wp), intent(in) :: bond0(nat, nat), bond(nat, nat)
      integer, allocatable, intent(out) :: actatoms(:)
      integer :: i, j
      real(wp) :: diff, diffthr

      ! threshold for bond order change !TODO tune this parameter
      diffthr = 0.5_wp

      allocate (actatoms(nat)) ! maximum number of active atoms is nat
      ! 0 means no active atom
      actatoms = 0

      do i = 1, nat
         ! if any bond order changes more than diffthr, atom is active
         do j = 1, nat
            diff = ABS(bond0(i, j) - bond(i, j))
            if (diff .gt. diffthr) then
               actatoms(i) = i

            end if
         end do
      end do

   end subroutine compare_bond

! check if structure dissociated upon optimization
   subroutine isdissociated(env, fname, dissociated)
      use mctc_env, only: error_type !, fatal_error
      use mctc_io, only: structure_type, read_structure, write_structure, &
            & filetype, get_filetype, to_symbol, to_number
      use xtb_mctc_convert, only: aatoau
      implicit none
      type(runtypedata) :: env
      character(len=*), intent(in) :: fname
      type(structure_type) :: mol
      type(error_type), allocatable :: error
      integer :: nat
      integer, allocatable :: iat(:)
      integer, allocatable :: fragi(:)
      real(wp), allocatable :: xyz(:, :)
      integer :: fragcount
      logical, intent(out) :: dissociated

      dissociated = .false.

      call read_structure(mol, fname, error, filetype%xyz)

      nat = mol%nat

      allocate (iat(nat))

      allocate (xyz(3, nat))
      iat = mol%num(mol%id)
      xyz = mol%xyz(:, :)

      call fragment_structure(nat, iat, xyz, 3.0_wp, 1, 0, fragi, fragcount) ! works better than mrec
      if (fragcount .gt. 1) then
         dissociated = .true.
      end if

   end subroutine isdissociated

   !*****************************************************************************************
   !* Taken and slightly modified from QCxMS mass spectra code https://github.com/qcxms/QCxMS
   !*  Subroutine for definition of two or more fragments
   !*  if at1 = 0 :  global separation in (nonbonded parts), beginning with atom at2
   !*  if at1 and at2 <> 0 : define fragments only if a at1-at2 bond (and no other) exists
   !*  if at1 and at2 = 0 : delete all fragment assignments
   !*  no bond if rcut times the cov.dist.
   !*  works better than mrec
   !*****************************************************************************************
   subroutine fragment_structure(nat, oz, xyz, rcut, at1, at2, frag, fragcount)
      use xtb_mctc_convert, only: aatoau
      implicit none
      integer, intent(in)  :: at1, at2, nat
      integer  :: i, j
      integer  :: attotal, currentfrag
      integer, intent(in)  :: oz(nat)
      integer, allocatable, intent(out) :: frag(:)
      integer, intent(out)  :: fragcount
      real(wp), intent(in) ::  xyz(3, nat)
      real(wp) :: rcov, r
      real(wp) :: rcut

      logical  :: finish
      logical, allocatable  :: connect(:, :)

      !&<
      ! Radius used in QCxMS (in au)
      real(wp), parameter :: Rad(118) = aatoau *  [ &
      & 0.32_wp,0.37_wp, & ! H,He
      & 1.30_wp,0.99_wp,0.84_wp,0.75_wp,0.71_wp,0.64_wp,0.60_wp,0.62_wp, & ! Li-Ne
      & 1.60_wp,1.40_wp,1.24_wp,1.14_wp,1.09_wp,1.04_wp,1.00_wp,1.01_wp, & ! Na-Ar
      & 2.00_wp,1.74_wp, & ! K,Ca
      &                 1.59_wp,1.48_wp,1.44_wp,1.30_wp,1.29_wp, & ! Sc-
      &                 1.24_wp,1.18_wp,1.17_wp,1.22_wp,1.20_wp, & ! -Zn
      &                 1.23_wp,1.20_wp,1.20_wp,1.18_wp,1.17_wp,1.16_wp, & ! Ga-Kr
      & 2.15_wp,1.90_wp, & ! Rb,Sr
      &                 1.76_wp,1.64_wp,1.56_wp,1.46_wp,1.38_wp, & ! Y-
      &                 1.36_wp,1.34_wp,1.30_wp,1.36_wp,1.40_wp, & ! -Cd
      &                 1.42_wp,1.40_wp,1.40_wp,1.37_wp,1.36_wp,1.36_wp, & ! In-Xe
      & 2.38_wp,2.06_wp, & ! Cs,Ba
      &         1.94_wp,1.84_wp,1.90_wp,1.88_wp,1.86_wp,1.85_wp,1.83_wp, & ! La-Eu
      &         1.82_wp,1.81_wp,1.80_wp,1.79_wp,1.77_wp,1.77_wp,1.78_wp, & ! Gd-Yb
      &                 1.74_wp,1.64_wp,1.58_wp,1.50_wp,1.41_wp, & ! Lu-
      &                 1.36_wp,1.32_wp,1.30_wp,1.30_wp,1.32_wp, & ! -Hg
      &                 1.44_wp,1.45_wp,1.50_wp,1.42_wp,1.48_wp,1.46_wp, & ! Tl-Rn
      & 2.42_wp,2.11_wp, & ! Fr,Ra
      &         2.01_wp,1.90_wp,1.84_wp,1.83_wp,1.80_wp,1.80_wp,& ! Ac-Pu
      &                                                         1.49_wp, & ! Am   ! below: from covalent 2009 covalent radii, such that it is complete up to 118
      &         1.49_wp,1.51_wp,1.51_wp,1.48_wp,1.50_wp,1.56_wp,1.58_wp, & ! Cm-No
      &                 1.45_wp,1.41_wp,1.34_wp,1.29_wp,1.27_wp, & ! Lr-
      &                 1.21_wp,1.16_wp,1.15_wp,1.09_wp,1.22_wp, & ! -Cn
      &                 1.36_wp,1.43_wp,1.46_wp,1.58_wp,1.48_wp,1.57_wp ] ! Nh-Og

      allocate (frag(nat))
      allocate (connect(nat,nat))
      connect(1:nat,1:nat) = .false.

      do i = 1,nat-1
         do j = i+1,nat
            r = sqrt((xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2 &
            & +(xyz(3,i)-xyz(3,j))**2)
            rcov = rcut*0.5_wp*(Rad(oz(i))+Rad(oz(j)))
            if (r .lt. rcov) then
               connect(i,j) = .true.
               connect(j,i) = .true.
            end if
         end do
      end do
      if ((at1 .eq. 0).and.(at2 .eq. 0)) then
         do i = 1,nat
            frag(i) = 1
         end do
         return
      else

         do i = 1,nat
            frag(i) = 0
         end do

         frag(at1) = 1
         attotal = 1

         if (at2 .ne. 0) then
            connect(at1,at2) = .false.
            connect(at2,at1) = .false.
         end if

         finish = .false.
         currentfrag = 0

         do while (attotal .ne. nat)

            currentfrag = currentfrag+1

            ! cycle through atoms and find connected ones

            do while (.not. (finish))
               finish = .true.
               do i = 1,nat
                  if (frag(i) .eq. currentfrag) then
                     do j = 1,nat
                        if (connect(i,j)) then
                           if (frag(j) .eq. 0) then
                              frag(j) = currentfrag
                              attotal = attotal+1
                              finish = .false.
                           elseif (frag(j) .eq. currentfrag) then
                              cycle
                           end if
                        end if
                     end do
                  end if
               end do
            end do

            ! find the first atom in the next fragment

            do i = 1,nat
               if (frag(i) .eq. 0) then
                  frag(i) = currentfrag+1
                  attotal = attotal+1
                  exit
               end if
            end do
            finish = .false.
         end do

      end if

      do i = 1,3 ! is enough, we only need to know we have 1, 2 or more than 2 fragments
         if (count(frag == i) .gt. 0) then
            fragcount = i
         end if
      end do

      deallocate (connect)
      return
   end subroutine fragment_structure

!> TODO subroutine to split fragments and optimize them separately
!> currently not used
!   subroutine separate_fragments(env, fname)
!      use mctc_env, only: error_type, get_argument!, fatal_error
!      use mctc_io, only: structure_type, read_structure, write_structure, &
!         & filetype, get_filetype, to_symbol, to_number
!      use qcxms2_data
!
!      implicit none
!      character(len=*)      :: fname
!      type(runtypedata) :: env
!      integer, allocatable :: iat(:)
!      type(structure_type) :: mol
!      type(error_type), allocatable :: error
!      integer, allocatable ::  fragi(:)
!      integer, allocatable ::  oind(:) ! save indices, order is first first fragment then second fragment atoms index is original index of pair
!      integer :: i
!      integer :: io
!      integer :: fragcount
!      integer :: incount
!      integer :: nat1, nat2
!      integer :: ich1, ich2
!      real(wp) :: rcut
!
!      rcut = 1.3_wp
!
!      call read_structure(mol, fname, error, filetype%xyz)
!      allocate (iat(mol%nat))
!      iat = mol%num(mol%id)
!      allocate (fragi(mol%nat))
!      allocate (oind(mol%nat))
!      call fragment_structure(mol%nat, iat, mol%xyz, rcut, 1, 0, fragi, fragcount)
!      io = makedir('frag1')
!      open (newunit=ich1, file="frag1/frag1.xyz")
!      io = makedir('frag2')
!      open (newunit=ich2, file="frag2/frag2.xyz")
!      nat1 = count(fragi == 1)
!      nat2 = count(fragi == 2)
!      write (ich1, *) nat1
!      write (ich1, *) ""
!      write (ich2, *) nat2
!      write (ich2, *) ""
!      incount = 0
!      do i = 1, mol%nat
!         if (fragi(i) == 1) then
!            write (ich1, *) mol%id(i), mol%xyz(1, i), mol%xyz(2, i), mol%xyz(3, i)
!            incount = incount + 1
!            oind(incount) = i
!         end if
!      end do
!      close (ich1)
!      do i = 1, mol%nat
!         if (fragi(i) == 2) then
!            write (ich2, *) mol%id(i), mol%xyz(1, i), mol%xyz(2, i), mol%xyz(3, i)
!            incount = incount + 1
!            oind(incount) = i
!         end if
!      end do
!      close (ich2)
!
!   end subroutine separate_fragments
!
!! rejoin fragments after seperate optimization
!   subroutine rejoin_fragments(nat, oind, at, xyz)
!      integer, intent(in) ::  oind(nat)
!      integer, intent(in) ::  at(nat)
!      real(wp), intent(in) ::  xyz(3, nat)
!      integer, intent(in) :: nat
!
!      integer :: ich
!      integer :: i, ind, j
!
!      open (newunit=ich, file="end.xyz")
!      write (ich, *) nat
!      write (ich, *) '' ! TODO maybe write energy also here ...
!      do i = 1, nat
!         do j = 1, nat
!            ! find position of index in oind
!            if (oind(j) == i) then
!               ind = j
!               exit
!            end if
!         end do
!         write (ich, *) at(ind), xyz(1, ind), xyz(2, ind), xyz(3, ind)
!      end do
!      !TODO continue here ...
!   end subroutine rejoin_fragments
   ! increase distance between fragments, better for transition state search and topology tools?
!   subroutine increase_fragdistance(xyz, nat, at, fragdist)
!
!      implicit none
!      integer, intent(in) :: nat
!      integer, intent(in) :: at(nat)
!      real(wp), intent(in) :: fragdist
!      real(wp), intent(inout) :: xyz(3, nat)
!      integer, allocatable :: at1(:), at2(:)
!      integer, allocatable :: fragi(:)
!      integer :: frag1, frag2, i, j, fragcount
!      real(wp) :: norm
!      real(wp), allocatable :: xyz1(:, :), xyz2(:, :) ! xyz of fragment 1 and 2
!      real(wp) :: cmass1(3), cmass2(3) ! center of mass of fragment 1 and 2
!      real(wp) :: rcut
!      rcut = 1.3_wp
!      frag1 = 0
!      frag2 = 0
!      write (*, '(a,f10.8)') "Increasing distance of fragments by ", fragdist, " angstroem"
!      allocate (fragi(nat))
!      call fragment_structure(nat, at, xyz(:, :), rcut, 1, 0, fragi, fragcount) ! better than mrec
!      if (fragcount .gt. 1) then
!      do j = 1, nat
!         if (fragi(j) == 1) then
!            frag1 = frag1 + 1
!         end if
!         if (fragi(j) == 2) then
!            frag2 = frag2 + 1
!         end if
!      end do
!      allocate (at1(frag1), at2(frag2))
!      allocate (xyz1(3, frag1), xyz2(3, frag2))
!      frag1 = 0
!      frag2 = 0
!      do j = 1, nat
!         if (fragi(j) == 1) then
!            frag1 = frag1 + 1
!            xyz1(:, frag1) = xyz(:, j)
!            at1(frag1) = at(j)
!         end if
!         if (fragi(j) == 2) then
!            frag2 = frag2 + 1
!            xyz2(:, frag2) = xyz(:, j)
!            at2(frag2) = at(j)
!         end if
!      end do
!      call CMAv(frag1, at1, xyz1, cmass1)
!      call CMAv(frag2, at2, xyz2, cmass2)
!      norm = sqrt((cmass2(1) - cmass1(1))**2 + (cmass2(2) - cmass1(2))**2 + (cmass2(3) - cmass1(3))**2)
!      do j = 1, nat
!         if (fragi(j) == 2) then
!            xyz(:, j) = xyz(:, j) + (cmass2 - cmass1)/norm*fragdist
!         end if
!      end do
!      deallocate (at1, at2, xyz1, xyz2)
!      end if
!   end subroutine increase_fragdistance

!============================================================!
! Most of these useful subroutines are taken from CREST (github.com/crest-lab/crest)
! (using the GPL-3.0 license like QCxMS2)
! and are slightly modified for QCxMS2
!============================================================!

!============================================================!
! subroutine wrxyz_file
! this is the typical quick write routine for TM coord files
! version for writing directly to a new fi
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!           at   - atom number as integer
!           xyz  - coordinates (in Angström)
!           comment - (OPTIONAL) comment line
!
! On Output: file written to "fname"
!============================================================!
   subroutine wrxyz_file(fname, nat, at, xyz, comment)
      implicit none
      character(len=*) :: fname
      integer :: nat
      integer :: at(nat)
      integer :: ich
      integer :: i
      real(wp) ::  xyz(3, nat)
      character(len=*), optional :: comment

      open (newunit=ich, file=fname, status='replace')
      write (ich, '(2x,i0)') nat
      if (present(comment)) then
         write (ich, '(a)') trim(comment)
      else
         write (ich, *)
      end if
      do i = 1, nat
         write (ich, '(1x,a2,1x,3f20.10)') i2e(at(i), 'nc'), xyz(1:3, i)
      end do
      close (ich)
      return
   end subroutine wrxyz_file

!============================================================!
! subroutine wrxyz_file_mask
! this is the typical quick write routine for TM coord files
! version for writing directly to a new file
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!           at   - atom number as integer
!           xyz  - coordinates (in Angström)
!           mask - a mask to determine to write which atoms
!           comment - (OPTIONAL) comment line
!
! On Output: file written to "fname"
!============================================================!
   subroutine wrxyz_file_mask(fname, nat, at, xyz, mask, comment)
      implicit none
      character(len=*) :: fname
      integer :: nat
      integer :: at(nat)
      integer :: i
      integer :: ich
      real(wp) ::  xyz(3, nat)
      logical :: mask(nat)
      integer :: maskednat
      character(len=*), optional :: comment
      open (newunit=ich, file=fname, status='replace')
      maskednat = count(mask(:))
      write (ich, '(2x,i0)') maskednat
      if (present(comment)) then
         write (ich, '(a)') trim(comment)
      else
         write (ich, *)
      end if
      do i = 1, nat
         if (mask(i)) then
            write (ich, '(1x,a2,1x,3f20.10)') i2e(at(i), 'nc'), xyz(1:3, i)
         end if
      end do
      close (ich)
      return
   end subroutine wrxyz_file_mask

   !============================================================!
! subroutine wrxyz
! quick subroutine to write a xyz file
! On Input: fname  - name of the ensemble file
!           nat    - number of atoms
!           at   - atom number as integer
!           xyz  - coordinates (in Angström)
!           comment - (OPTIONAL) comment line
!
! On Output: file written to "fname"
!============================================================!
   subroutine wrxyz(fname, nat, at, xyz, comment)
      implicit none
      character(len=*) :: fname
      integer :: nat
      integer :: at(nat)
      integer :: ich
      integer :: i
      real(wp) ::  xyz(3, nat)
      character(len=*), optional :: comment
      logical :: ex
      inquire (file=fname, exist=ex)
      if (ex) then
         open (newunit=ich, file=fname, status='old', position="append")
      else
         open (newunit=ich, file=fname, status='replace')
      end if
      write (ich, '(2x,i0)') nat
      if (present(comment)) then
         write (ich, '(a)') trim(comment)
      else
         write (ich, *)
      end if
      do i = 1, nat
         write (ich, '(1x,a2,1x,3f20.10)') i2e(at(i), 'nc'), xyz(1:3, i)
      end do
      close (ich)
      return
   end subroutine wrxyz

!========================================================================================!
!> subroutine CMAv
!> calculate a CMA coordinats and save them to vec
!>----------------------------------
   subroutine CMAv(nat, at, coord, vec)
      use isotope_pattern, only: ams
      implicit none
      integer, intent(in) :: nat
      integer, intent(in) :: at(nat)
      real(wp), intent(in) :: coord(3, nat)
      real(wp), intent(out) :: vec(3)
      integer :: i
      real(wp) :: sumw, sumwx, sumwy, sumwz, atmass
      sumw = 1.d-20
      sumwx = 0.d0
      sumwy = 0.d0
      sumwz = 0.d0
      do i = 1, nat
         atmass = ams(at(i))
         sumw = sumw + atmass
         sumwx = sumwx + atmass*coord(1, i)
         sumwy = sumwy + atmass*coord(2, i)
         sumwz = sumwz + atmass*coord(3, i)
      end do
      vec(1) = sumwx/sumw
      vec(2) = sumwy/sumw
      vec(3) = sumwz/sumw
      return
   end subroutine CMAv

!========================================================================================!
!> subroutine center_of_geometry
!> calculate center of molecular coordinates and save them to vec
!>----------------------------------

   subroutine center_of_geometry(nuc, xyz, cg)

      integer, intent(in) :: nuc
      integer :: j
      integer :: normmass

      real(wp), intent(in) :: xyz(3, nuc)
      real(wp), intent(out) :: cg(3)

      normmass = 0
      cg(:) = 0
      do j = 1, nuc

         cg(:) = cg(:) + 1*xyz(:, j)

         normmass = normmass + 1
      end do

      cg(1) = cg(1)/normmass
      cg(2) = cg(2)/normmass
      cg(3) = cg(3)/normmass

   end subroutine center_of_geometry

!============================================================!
! subroutine rdxmol
! read a struncture in the *.xyz (Xmol) style.
! The commentary (second) line is ignored
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!
! On Output: at   - atom number as integer
!            xyz  - coordinates (in Angström)
!            comment - (OPTIONAL) commentary line of the file
!============================================================!
   subroutine rdxmol(fname, nat, at, xyz, comment)
      implicit none
      character(len=*), intent(in) :: fname
      integer, intent(in) :: nat
      integer, intent(inout)  :: at(nat)
      real(wp), intent(inout) :: xyz(3, nat)
      character(len=*), optional :: comment
      character(len=6) :: sym
      integer :: ich, io, i
      integer :: dum
      character(len=256) :: atmp
      open (newunit=ich, file=fname)
      read (ich, *, iostat=io) dum
      if (nat .ne. dum) then
         error stop 'error while reading input coordinates'
      end if
      read (ich, '(a)') atmp !--commentary line
      if (present(comment)) comment = trim(adjustl(atmp))
      do i = 1, nat
         read (ich, '(a)', iostat=io) atmp
         if (io < 0) exit
         atmp = adjustl(atmp)
         call coordline(atmp, sym, xyz(1:3, i), io)
         if (io < 0) then
            write (*, *) 'error while reading coord line. EOF'
            exit
         end if
         at(i) = e2i(sym)
      end do
      close (ich)
      return
   end subroutine rdxmol

!============================================================!
! read a line of coordinates and determine by itself
! if the format is x,y,z,at or at,x,y,z
!============================================================!
   subroutine coordline(line, sym, xyz, io)
      implicit none
      character(len=*) :: line
      character(len=*) :: sym
      real(wp) :: xyz(3)
      integer, intent(out) :: io

      io = 0
      read (line, *, iostat=io) xyz(1:3), sym
      if (io .ne. 0) then
         read (line, *, iostat=io) sym, xyz(1:3)
         !if(io.ne.0)then
         !  error stop 'error while reading coord line'
         !endif
      end if

      return
   end subroutine coordline

   subroutine get_CN(fname, cn, bond, lrange)
      use mctc_io, only: structure_type, read_structure, write_structure, &
             & filetype, get_filetype, to_symbol, to_number
      use mctc_env, only: error_type, get_argument!, fatal_error
      use xtb_mctc_convert, only: aatoau
      use qcxms2_data, only: bohr
      implicit none
      character(len=80) :: fname

      real(wp), allocatable, intent(out) :: cn(:), bond(:, :)
      real(wp), allocatable :: xyz(:, :)
      real(wp) :: rcovi, rcovj, rij(3), r2, r, rco
      real(wp) :: damp, den
      real(wp) :: cn_thr
      logical, intent(in) :: lrange

      integer :: nat, ati, atj
      integer :: i, j
      integer, allocatable :: iat(:)
      type(structure_type) :: mol
      type(error_type), allocatable :: error
      !> Pauling electronegativities, used for the covalent coordination number.
      real(wp), parameter :: pauling_en(118) = [ &
      & 2.20_wp, 3.00_wp, & ! H,He
      & 0.98_wp, 1.57_wp, 2.04_wp, 2.55_wp, 3.04_wp, 3.44_wp, 3.98_wp, 4.50_wp, & ! Li-Ne
      & 0.93_wp, 1.31_wp, 1.61_wp, 1.90_wp, 2.19_wp, 2.58_wp, 3.16_wp, 3.50_wp, & ! Na-Ar
      & 0.82_wp, 1.00_wp, & ! K,Ca
      &                 1.36_wp, 1.54_wp, 1.63_wp, 1.66_wp, 1.55_wp, & ! Sc-
      &                 1.83_wp, 1.88_wp, 1.91_wp, 1.90_wp, 1.65_wp, & ! -Zn
      &                 1.81_wp, 2.01_wp, 2.18_wp, 2.55_wp, 2.96_wp, 3.00_wp, & ! Ga-Kr
      & 0.82_wp, 0.95_wp, & ! Rb,Sr
      &                 1.22_wp, 1.33_wp, 1.60_wp, 2.16_wp, 1.90_wp, & ! Y-
      &                 2.20_wp, 2.28_wp, 2.20_wp, 1.93_wp, 1.69_wp, & ! -Cd
      &                 1.78_wp, 1.96_wp, 2.05_wp, 2.10_wp, 2.66_wp, 2.60_wp, & ! In-Xe
      & 0.79_wp, 0.89_wp, & ! Cs,Ba
      &         1.10_wp, 1.12_wp, 1.13_wp, 1.14_wp, 1.15_wp, 1.17_wp, 1.18_wp, & ! La-Eu
      &         1.20_wp, 1.21_wp, 1.22_wp, 1.23_wp, 1.24_wp, 1.25_wp, 1.26_wp, & ! Gd-Yb
      &                 1.27_wp, 1.30_wp, 1.50_wp, 2.36_wp, 1.90_wp, & ! Lu-
      &                 2.20_wp, 2.20_wp, 2.28_wp, 2.54_wp, 2.00_wp, & ! -Hg
      &                 1.62_wp, 2.33_wp, 2.02_wp, 2.00_wp, 2.20_wp, 2.20_wp, & ! Tl-Rn  !!!below: only dummies!!!
      & 1.50_wp, 1.50_wp, & ! Fr,Ra
      &         1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, & ! Ac-Am
      &         1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, & ! Cm-No
      &                 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, & ! Rf-
      &                 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, & ! Rf-Cn
      &                 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp] ! Nh-Og

      !> Covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
      !> 188-197), values for metals decreased by 10%.
      !> As used in D3
      real(wp), parameter :: rcov(1:118) = [ &
       & 0.32_wp, 0.46_wp, & ! H,He
       & 1.20_wp, 0.94_wp, 0.77_wp, 0.75_wp, 0.71_wp, 0.63_wp, 0.64_wp, 0.67_wp, & ! Li-Ne
       & 1.40_wp, 1.25_wp, 1.13_wp, 1.04_wp, 1.10_wp, 1.02_wp, 0.99_wp, 0.96_wp, & ! Na-Ar
       & 1.76_wp, 1.54_wp, & ! K,Ca
       &                 1.33_wp, 1.22_wp, 1.21_wp, 1.10_wp, 1.07_wp, & ! Sc-
       &                 1.04_wp, 1.00_wp, 0.99_wp, 1.01_wp, 1.09_wp, & ! -Zn
       &                 1.12_wp, 1.09_wp, 1.15_wp, 1.10_wp, 1.14_wp, 1.17_wp, & ! Ga-Kr
       & 1.89_wp, 1.67_wp, & ! Rb,Sr
       &                 1.47_wp, 1.39_wp, 1.32_wp, 1.24_wp, 1.15_wp, & ! Y-
       &                 1.13_wp, 1.13_wp, 1.08_wp, 1.15_wp, 1.23_wp, & ! -Cd
       &                 1.28_wp, 1.26_wp, 1.26_wp, 1.23_wp, 1.32_wp, 1.31_wp, & ! In-Xe
       & 2.09_wp, 1.76_wp, & ! Cs,Ba
       &         1.62_wp, 1.47_wp, 1.58_wp, 1.57_wp, 1.56_wp, 1.55_wp, 1.51_wp, & ! La-Eu
       &         1.52_wp, 1.51_wp, 1.50_wp, 1.49_wp, 1.49_wp, 1.48_wp, 1.53_wp, & ! Gd-Yb
       &                 1.46_wp, 1.37_wp, 1.31_wp, 1.23_wp, 1.18_wp, & ! Lu-
       &                 1.16_wp, 1.11_wp, 1.12_wp, 1.13_wp, 1.32_wp, & ! -Hg
       &                 1.30_wp, 1.30_wp, 1.36_wp, 1.31_wp, 1.38_wp, 1.42_wp, & ! Tl-Rn
       & 2.01_wp, 1.81_wp, & ! Fr,Ra
       &         1.67_wp, 1.58_wp, 1.52_wp, 1.53_wp, 1.54_wp, 1.55_wp, 1.49_wp, & ! Ac-Am
       &         1.49_wp, 1.51_wp, 1.51_wp, 1.48_wp, 1.50_wp, 1.56_wp, 1.58_wp, & ! Cm-No
       &                 1.45_wp, 1.41_wp, 1.34_wp, 1.29_wp, 1.27_wp, & ! Lr-
       &                 1.21_wp, 1.16_wp, 1.15_wp, 1.09_wp, 1.22_wp, & ! -Cn
       &                 1.36_wp, 1.43_wp, 1.46_wp, 1.58_wp, 1.48_wp, 1.57_wp] & ! Nh-Og
       & *aatoau*4.0_wp/3.0_wp

      cn_thr = 625.0_wp  !> 25.0^2 as in tblite

      call read_structure(mol, fname, error, filetype%xyz)
      allocate (iat(mol%nat))
      iat = mol%num(mol%id)
      nat = mol%nat
      ! leave it in Bohr !!!
      xyz = mol%xyz

      allocate (cn(nat), source=0.0_wp)

      allocate (bond(nat, nat), source=0.0_wp)

      cn(:) = 0.0_wp
      !>--- actual calculation

      do i = 1, nat
         ati = iat(i)
         rcovi = RCOV(ati)

         do j = 1, i - 1
            rij(:) = xyz(:, i) - xyz(:, j)
            r2 = sum(rij**2)
            !>--- cycle cutoff
            if (r2 > cn_thr) cycle
            r = sqrt(r2)
            atj = iat(j)
            rcovj = RCOV(atj)
            rco = rcovi + rcovj

            !>--- select the correct CN version

            !den = PAULING_EN(atj)-PAULING_EN(ati)
            !damp = cn_damp_d4(rco,r,den)

            if (lrange) then
               damp = cn_damp_exp_lrange(rco, r)
            else
               damp = cn_damp_exp(rco, r)
            end if

            bond(j, i) = damp
            bond(i, j) = damp

            cn(i) = cn(i) + damp

            if (i /= j) then
               cn(j) = cn(j) + damp*1.0_wp
            end if

         end do
      end do

      do i = 1, nat
         write (*, *) "CN of atom ", i, " is ", cn(i)
      end do

   end subroutine get_CN

   function cn_damp_d4(rco, r, den) result(damp)
      !*********************************************
      !* EN-dependent erf-CN damping factor from D4
      !*********************************************
      implicit none
      real(wp) :: damp
      real(wp), intent(in) :: rco, r, den
      real(wp) :: tmp, tmperf
      real(wp), parameter :: k4 = 4.10451_wp
      real(wp), parameter :: k5 = 19.08857_wp
      real(wp), parameter :: k6 = 2*11.28174_wp**2
      real(wp), parameter :: kcn = 7.50_wp
      tmp = k4*exp(-((abs(den) + k5)**2)/k6)
      tmperf = 0.5_wp*(1.0_wp + erf(-kcn*(r - rco)/rco))
      damp = tmp*tmperf
      !    damp = tmp*cn_damp_erf(rco,r)
   end function cn_damp_d4

   function cn_damp_exp(rco, r) result(damp)
      !*****************************************
      !* classic exponential CN damping factor
      !*****************************************
      implicit none
      real(wp) :: damp
      real(wp), intent(in) :: rco, r
      real(wp), parameter :: kcn = 4.0_wp !16.0_wp
      real(wp), parameter :: k2 = 1.5_wp
      damp = 1.0_wp/(1.0_wp + exp(-kcn*(k2*rco/r - 1.0_wp)))
   end function cn_damp_exp

   function cn_damp_exp_lrange(rco, r) result(damp)
      !*****************************************
      !* classic exponential CN damping factor
      !*****************************************
      implicit none
      real(wp) :: damp
      real(wp), intent(in) :: rco, r
      real(wp), parameter :: kcn = 8.0_wp !16.0_wp
      damp = 1.0_wp/(1.0_wp + exp(-kcn*(rco/r - 1.0_wp)))
   end function cn_damp_exp_lrange

end module structools
