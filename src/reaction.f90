module reaction
   use iomod
   use qcxms2_data
   use qmmod
   use xtb_mctc_convert
   use utility
   use mcsimu
   implicit none
contains

!> This routine takes energies of all fragment pairs (if set as option including rrho)
! as input and sorts out those which are too high in energy
! determined by same exponential formula as taken for barriers in mcsium to sort out fragment pairs
! with high relative energy but of course with conservative thresholds
! By default 3.0 time the average IEE is taken as threshold at which
! the realive intensity is computed and compared to the pthr threshold, set very low here for now
!TODO optimize the thresholds to reduce computation time
   subroutine sortouthighe(env, sortlevel, etots_in, rrho, npairs_in, npairs_out, fragdirs_in, fragdirs_out, scaleIEE)
      implicit none
      character(len=80) :: fname2
      character(len=80) :: sortlevel
      character(len=80) :: query
      character(len=80) :: fragdirs_in(:, :)
      character(len=80), allocatable, intent(out)   :: fragdirs_out(:, :)
      integer, intent(in) :: npairs_in
      integer, intent(out) :: npairs_out
      type(runtypedata) :: env
      real(wp) :: e_start, rrho_start ! energy and rrho contribution of starting fragment
      real(wp), allocatable, intent(in) :: etots_in(:)
      real(wp), allocatable :: etots! total energy of fragments and fragment pairs
      real(wp) :: e_reac, sumreac, dzpve ! sum of all reaction energies up to this point
      real(wp) :: pthr, IEE, IEE2
      real(wp), intent(in) :: scaleIEE ! scaling factor of avg IEE to account for occuring higher energies  ! multiply IEE by two here to be safe
      real(wp) :: de_min, krel
      real(wp), allocatable :: de(:) ! reaction energy
      integer :: nat
      integer :: i, j, io
      integer :: chrg, uhf, spin
      integer :: nvib ! vibrational degrees of freedom
      logical :: ex, ldum
      logical, intent(in) :: rrho ! use rrho thermocorrection for sorting
      logical, allocatable :: ishdiss(:) ! is this a H-dissociation reaction
      real(wp) :: scaleintdiss ! scaling factor for H-dissociation reactions
      character(len=80) :: fname

      write (*, *) "Sorting out fragmentpairs which are too high in energy"

      allocate (de(npairs_in))
      ! TODO maybe also dependent on level of theory but for now use this parameters
      pthr = 0.01_wp  ! probability threshold to account for fragmentation ! TODO TUNE ME

      ! first determine average energy of run
      call rdshort_int(trim(env%startdir)//"/in.xyz", nat)
      IEE = env%ieeatm*nat

      nvib = 3*nat - 6 ! vibrational degrees of freedom

      ! for diatomic molecules
      if (nvib .eq. 0) then
         nvib = 1 !
      end if

      ! get number of atoms of current fragment
      inquire (file="fragment.xyz", exist=ex)
      if (ex) then
         fname = 'fragment.xyz'
      end if
      inquire (file="isomer.xyz", exist=ex)
      if (ex) then
         fname = 'isomer.xyz'
      end if

      inquire (file=".CHRG", exist=ex)
      if (ex) then
         call rdshort_int('.CHRG', chrg)
      else
         chrg = env%chrg
      end if
      inquire (file=".UHF", exist=ex)
      if (ex) then
         call rdshort_int('.UHF', uhf)
      else
         call printspin(env, fname, spin, chrg)
         ! uhf in .UHF formalism of xtb is number of unpaired electrons
         uhf = spin - 1
         !correct for H+
         if (uhf == -2) uhf = 2
      end if

      write (query, '(a,1x,i0,1x,i0)') trim(sortlevel)//" sp", chrg, uhf
      call grepval('qmdata', trim(query), ldum, e_start)
      if (.not. ldum) then
         write (*, *) "ERROR: no energy found for starting fragment"
         stop
      end if
      if (rrho) then
         write (query, '(a,1x,i0,1x,i0)') trim(env%geolevel)//" bhess", chrg, uhf
         call grepval('qmdata', trim(query), ldum, rrho_start)
         if (.not. ldum) then
            write (*, *) "ERROR: no rrho contribution found for starting fragment"
            stop
         end if
      else
         rrho_start = 0.0_wp
      end if
      e_start = e_start + rrho_start
      call rdshort_real("sumreac_"//trim(sortlevel), sumreac)
      if (env%printlevel .eq. 3) write (*, *) "estart, rrhostart, and sum of reaction is", e_start, rrho_start, sumreac

      ! find lowest reaction energy and compute approximated relative rates
      ! with regard to fastest rate with lowest reaction energy
      ! negative reaction energies are set to zero as formula is not well defined for these energies

      do i = 1, npairs_in
         !de(i) = (etots_in(i) - e_start)*autoev
         inquire (file=trim(fragdirs_in(i, 1))//"/de_"//trim(env%tslevel), exist=ex)
         if (ex) then
            call rdshort(trim(fragdirs_in(i, 1))//"/de_"//trim(env%tslevel), de(i))
         else
            de(i) = 10000.0_wp ! just ridiculous high energy to sort this out later, the calculaiton failed here
         end if
         !if (de(i) .lt. 0.0_wp) de(i) = 0.0_wp ! set to 0 if negative TODO CHECKME negative barriers????
      end do

      de_min = 5.0_wp

      ! take only fragmentations here, not isomerizations
      do i = 1, npairs_in
         ! consider only fragmentation reactions for minimum energy here
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            if (de(i) .lt. de_min) de_min = de(i)
         end if
      end do

      !de_min = minval(de)
      if (de_min .lt. 5.0_wp) de_min = 5.0_wp ! to be save here, we take at least 5.0 eV as threshold
      if (env%mode == " cid") then
         if (de_min .lt. 5.0_wp) de_min = 5.0_wp ! for CID we need a higher threshold as Dissociation energies are usually higher just to be safe here
      end if
      if (env%printlevel .eq. 3) write (*, *) "dmin is", de_min, "at index", minloc(de)
      ! substract consumed energy up to this point
      IEE = IEE - sumreac

      call rdshort_int(trim(fname), nat)
      nvib = 3*nat - 6 ! vibrational degrees of freedom

      call findhdiss(env, fname, nat, npairs_in, fragdirs_in, ishdiss, scaleintdiss)
      do i = 1, npairs_in
         if (ishdiss(i)) then
            IEE2 = IEE*scaleintdiss
         else
            IEE2 = IEE
         end if

         if (env%eyring) then
            krel = calc_releyring(IEE2, de(i) - de_min, nvib)
         else
            krel = EXP(-(nvib - 1)*(de(i) - de_min)/(IEE2*scaleIEE))
         end if

         if (krel .lt. pthr) then
            write (*, '(x,a,x,a,x,a,f5.1,a,f5.1,a)') "fragment pair", trim(fragdirs_in(i, 1)) &
            & , " is sorted out with an reaction energy of ", de(i)*evtoau*autokcal, " kcal/mol (", de(i), " eV)"
            fragdirs_in(i, 1) = '' ! make it empty or no p in string
         end if
      end do
      ! rewrite fragdir list so that empty entries are removed
      call sortout0elements(npairs_in, npairs_out, fragdirs_in, fragdirs_out, env%removedirs)
      write (*, *) "from ", npairs_in, " fragment pairs ", npairs_out, " are left after sorting step"
      write (*, *) "Remaining Pairs:", npairs_out
      deallocate (de)
   end subroutine sortouthighe

!> compute reaction energy of all fragment pairs
! reaction can be A -> B (isomerization) or A -> B + C (fragmentation)
! TODO implement include check for H atom for all job-loops
! to these calculations at first so that more cores are available
! for real calculations
   subroutine calcreactionenergy(env, npairs, fragdirs, e_pairs)
      implicit none
      character(len=1024) :: jobcall, cleanupcall
      character(len=80) :: fname
      character(len=80) :: job
      character(len=80) :: query, fout, pattern
      character(len=80), intent(in)  :: fragdirs(:, :)
      character(len=80), allocatable :: dirs(:), dirs0(:)
      integer, intent(in) :: npairs
      integer :: njobs, njobs0
      integer :: i, ich, ierr, j, io
      integer :: chrg, uhf, spin
      real(wp), allocatable, intent(out) :: e_pairs(:) ! energy of fragment pair (can be sum of two fragments)
      type(runtypedata) :: env
      real(wp) :: edum, e_start, e_f1, e_f2 ! energy of starting fragments and fragment pairs
      real(wp), allocatable :: e_rrho(:), etot(:) ! total energy of fragments and fragment pairs
      real(wp) :: rrho_start, rrho_f1, rrho_f2 ! rrho contribution of starting fragments and fragment pairs
      real(wp), allocatable :: drrho(:), de(:)
      real(wp) :: sumreac, sumreac0 ! sum of all reaction energies to come to this point
      logical :: ex, ldum, failed, there
      write (*, *) "Computing reaction energy of fragmentation"

      allocate (e_rrho(npairs), drrho(npairs), de(npairs))

      inquire (file="fragment.xyz", exist=ex)
      if (ex) then
         fname = 'fragment.xyz'
      end if
      inquire (file="isomer.xyz", exist=ex)
      if (ex) then
         fname = 'isomer.xyz'
      end if

      inquire (file=".CHRG", exist=ex)
      if (ex) then
         call rdshort_int('.CHRG', chrg)
         if (env%printlevel .eq. 3) write (*, *) "read .CHRG charge is", chrg
      else
         chrg = env%chrg
      end if
      inquire (file=".UHF", exist=ex)
      if (ex) then
         call rdshort_int('.UHF', uhf)
      else
         call printspin(env, fname, spin, chrg)
         ! uhf in .UHF formalism of xtb is number of unpaired electrons
         uhf = spin - 1
         !correct for H+
         if (uhf == -2) uhf = 2
      end if

      if (env%exstates .gt. 0) then
         job = "tddft"
      else
         job = "sp"
      end if
      write (query, '(a,1x,i0,1x,i0)') trim(env%tslevel)//" "//trim(job), chrg, uhf

      call grepval('qmdata', trim(query), ldum, e_start)
      if (.not. ldum) then
         write (*, *) "ERROR: no energy found for starting fragment"
         call printpwd()
         write (*, *) "query is", trim(query)
         stop
      end if
      if (env%bhess) then
         write (query, '(a,1x,i0,1x,i0)') trim(env%geolevel)//" bhess", chrg, uhf
         call grepval('qmdata', trim(query), ldum, rrho_start)
         if (.not. ldum) then
            write (*, *) "ERROR: no GRRHO found for starting fragment" ! TODO FIXME we could calculate it but probably something was terrible wrong there
            stop
         end if
      end if

      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      do i = 1, npairs
         if (index(fragdirs(i, 3), 'p') .ne. 0) then
            fname = 'fragment.xyz'
            do j = 2, 3
               call chdir(trim(fragdirs(i, j)))
               call prepqm(env, fname, env%tslevel, job, jobcall, fout, pattern, cleanupcall, there, .false.)
               if (.not. there) then
                  njobs = njobs + 1
                  call append_char(dirs, trim(fragdirs(i, j)))
               end if
               call chdir(trim(env%path))
            end do
         else
            fname = 'isomer.xyz'
            call chdir(trim(fragdirs(i, 1)))
            call prepqm(env, fname, env%tslevel, job, jobcall, fout, pattern, cleanupcall, there, .false.)
            if (.not. there) then
               njobs = njobs + 1
               call append_char(dirs, trim(fragdirs(i, 1)))
            end if
            call chdir(trim(env%path))
         end if
      end do

      ! now actually prepare jobs
      call setompthreads(env, njobs)
      deallocate (dirs)
      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      do i = 1, npairs
         if (index(fragdirs(i, 3), 'p') .ne. 0) then
            fname = 'fragment.xyz'
            do j = 2, 3
               call chdir(trim(fragdirs(i, j)))
               call prepqm(env, fname, env%tslevel, job, jobcall, fout, pattern, cleanupcall, there, .false.)
               if (.not. there) then
                  njobs = njobs + 1
                  call append_char(dirs, trim(fragdirs(i, j)))
               end if
               call chdir(trim(env%path))
            end do
         else
            fname = 'isomer.xyz'
            call chdir(trim(fragdirs(i, 1)))
            call prepqm(env, fname, env%tslevel, job, jobcall, fout, pattern, cleanupcall, there, .false.)
            if (.not. there) then
               njobs = njobs + 1
               call append_char(dirs, trim(fragdirs(i, 1)))
            end if
            call chdir(trim(env%path))
         end if
      end do

      write (*, *) "Starting ", njobs, " fragment pair calculations in parallel at ", trim(env%tslevel), " level"
      call omp_samejobcall(njobs, dirs, jobcall)

      dirs0 = dirs
      njobs0 = njobs

      ! restart failed calculations
      deallocate (dirs)
      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      do i = 1, npairs
         if (index(fragdirs(i, 3), 'p') .ne. 0) then
            fname = 'fragment.xyz'
            do j = 2, 3
               call chdir(trim(fragdirs(i, j)))
               call readoutqm(env, fname, env%tslevel, job, fout, pattern, edum, failed)
               if (failed) then
                  njobs = njobs + 1
                  call append_char(dirs, trim(fragdirs(i, j)))
                  call prepqm(env, fname, env%tslevel, job, jobcall, fout, pattern, cleanupcall, there, .true.)
               end if
               call chdir(trim(env%path))
            end do
         else
            fname = 'isomer.xyz'
            call chdir(trim(fragdirs(i, 1)))
            call readoutqm(env, fname, env%tslevel, job, fout, pattern, edum, failed)
            if (failed) then
               call prepqm(env, fname, env%tslevel, job, jobcall, fout, pattern, cleanupcall, there, .true.)
               njobs = njobs + 1
               call append_char(dirs, trim(fragdirs(i, 1)))
            end if
            call chdir(trim(env%path))
         end if
      end do

      ! now actually prepare jobs
      call setompthreads(env, njobs)
      deallocate (dirs)
      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      do i = 1, npairs
         if (index(fragdirs(i, 3), 'p') .ne. 0) then
            fname = 'fragment.xyz'
            do j = 2, 3
               call chdir(trim(fragdirs(i, j)))
               call readoutqm(env, fname, env%tslevel, job, fout, pattern, edum, failed)
               if (failed) then
                  njobs = njobs + 1
                  call append_char(dirs, trim(fragdirs(i, j)))
                  call prepqm(env, fname, env%tslevel, job, jobcall, fout, pattern, cleanupcall, there, .true.)
               end if
               call chdir(trim(env%path))
            end do
         else
            fname = 'isomer.xyz'
            call chdir(trim(fragdirs(i, 1)))
            call readoutqm(env, fname, env%tslevel, job, fout, pattern, edum, failed)
            if (failed) then
               call prepqm(env, fname, env%tslevel, job, jobcall, fout, pattern, cleanupcall, there, .true.)
               njobs = njobs + 1
               call append_char(dirs, trim(fragdirs(i, 1)))
            end if
            call chdir(trim(env%path))
         end if
      end do

      write (*, *) "Starting ", njobs, " restart product calculations in parallel at ", trim(env%tslevel), " level"
      call omp_samejobcall(njobs, dirs, jobcall)

      ! read out energies
      allocate (e_pairs(npairs))
      do i = 1, npairs
         if (index(fragdirs(i, 3), 'p') .ne. 0) then
            fname = 'fragment.xyz'
            call chdir(trim(fragdirs(i, 2)))
            call readoutqm(env, fname, env%tslevel, job, fout, pattern, e_f1, failed)
            if (failed) then
               e_pairs(i) = 10000.0_wp ! jsut ridiculous high energy to sort this out later
               call chdir(trim(env%path))
               cycle
            end if
            call chdir(trim(env%path))
            call chdir(trim(fragdirs(i, 3)))
            call readoutqm(env, fname, env%tslevel, job, fout, pattern, e_f2, failed)
            if (failed) then
               e_pairs(i) = 10000.0_wp ! jsut ridiculous high energy to sort this out later
               call chdir(trim(env%path))
               cycle
            end if
            call chdir(trim(env%path))
            e_pairs(i) = e_f1 + e_f2
            ! DEL write(*,*) "e_pairs is", e_f1, e_f2, e_pairs(i)
         else
            fname = 'isomer.xyz'
            call chdir(trim(fragdirs(i, 1)))
            call readoutqm(env, fname, env%tslevel, job, fout, pattern, e_pairs(i), failed)
            if (failed) then
               e_pairs(i) = 10000.0_wp ! jsut ridiculous high energy to sort this out later
               call chdir(trim(env%path))
               cycle
            end if
            call chdir(trim(env%path))
         end if
         de(i) = (e_pairs(i) - e_start)*autoev

         ! if bhess then its written later
         if (.not. env%bhess) call wrshort_real("de_"//trim(env%tslevel), de(i))
      end do

      ! now cleanup all calculations
      call omp_samejobcall(njobs0, dirs0, cleanupcall, .false.)

      ! include rrho contribution
      if (env%bhess) then
         ! restart failed calculations
         deallocate (dirs)
         njobs = 0
         allocate (dirs(1)) ! first element always empty
         dirs = ['']

         do i = 1, npairs
            if (index(fragdirs(i, 3), 'p') .ne. 0) then
               fname = 'fragment.xyz'
               do j = 2, 3
                  call chdir(trim(fragdirs(i, j)))
                  call prepqm(env, fname, env%geolevel, 'bhess', jobcall, fout, pattern, cleanupcall, there, .false.)
                  if (.not. there) then
                     njobs = njobs + 1
                     call append_char(dirs, trim(fragdirs(i, j)))
                  end if
                  call chdir(trim(env%path))
               end do
            else
               fname = 'isomer.xyz'
               call chdir(trim(fragdirs(i, 1)))
               call prepqm(env, fname, env%geolevel, 'bhess', jobcall, fout, pattern, cleanupcall, there, .false.)
               if (.not. there) then
                  njobs = njobs + 1
                  call append_char(dirs, trim(fragdirs(i, 1)))
               end if
               call chdir(trim(env%path))
            end if
         end do

         call setompthreads(env, njobs)
         ! now actually prepare jobs
         deallocate (dirs)
         njobs = 0
         allocate (dirs(1)) ! first element always empty
         dirs = ['']

         do i = 1, npairs
            if (index(fragdirs(i, 3), 'p') .ne. 0) then
               fname = 'fragment.xyz'
               do j = 2, 3
                  call chdir(trim(fragdirs(i, j)))
                  call prepqm(env, fname, env%geolevel, 'bhess', jobcall, fout, pattern, cleanupcall, there, .false.)
                  if (.not. there) then
                     njobs = njobs + 1
                     call append_char(dirs, trim(fragdirs(i, j)))
                  end if
                  call chdir(trim(env%path))
               end do
            else
               fname = 'isomer.xyz'
               call chdir(trim(fragdirs(i, 1)))
               call prepqm(env, fname, env%geolevel, 'bhess', jobcall, fout, pattern, cleanupcall, there, .false.)
               if (.not. there) then
                  njobs = njobs + 1
                  call append_char(dirs, trim(fragdirs(i, 1)))
               end if
               call chdir(trim(env%path))
            end if
         end do

         write (*, *) "Starting ", njobs, " SPH calculations in parallel at ", trim(env%geolevel), " level"
         call omp_samejobcall(njobs, dirs, jobcall)

         !restart calculations
         deallocate (dirs)
         njobs = 0
         allocate (dirs(1)) ! first element always empty
         dirs = ['']

         do i = 1, npairs
            if (index(fragdirs(i, 3), 'p') .ne. 0) then
               fname = 'fragment.xyz'
               do j = 2, 3
                  call chdir(trim(fragdirs(i, j)))
                  call readoutqm(env, fname, env%geolevel, 'bhess', fout, pattern, edum, failed)
                  if (failed) then
                     call prepqm(env, fname, env%geolevel, 'bhess', jobcall, fout, pattern, cleanupcall, there, .false.)
                     njobs = njobs + 1
                     call append_char(dirs, trim(fragdirs(i, j)))
                  end if
                  call chdir(trim(env%path))
               end do
            else
               fname = 'isomer.xyz'
               call chdir(trim(fragdirs(i, 1)))
               call readoutqm(env, fname, env%geolevel, 'bhess', fout, pattern, edum, failed)
               if (failed) then
                  call prepqm(env, fname, env%geolevel, 'bhess', jobcall, fout, pattern, cleanupcall, there, .false.)
                  njobs = njobs + 1
                  call append_char(dirs, trim(fragdirs(i, 1)))
               end if
               call chdir(trim(env%path))
            end if
         end do

         call setompthreads(env, njobs)

         ! now actually prepare jobs
         deallocate (dirs)
         njobs = 0
         allocate (dirs(1)) ! first element always empty
         dirs = ['']

         do i = 1, npairs
            if (index(fragdirs(i, 3), 'p') .ne. 0) then
               fname = 'fragment.xyz'
               do j = 2, 3
                  call chdir(trim(fragdirs(i, j)))
                  call readoutqm(env, fname, env%geolevel, 'bhess', fout, pattern, edum, failed)
                  if (failed) then
                     call prepqm(env, fname, env%geolevel, 'bhess', jobcall, fout, pattern, cleanupcall, there, .false.)
                     njobs = njobs + 1
                     call append_char(dirs, trim(fragdirs(i, j)))
                  end if
                  call chdir(trim(env%path))
               end do
            else
               fname = 'isomer.xyz'
               call chdir(trim(fragdirs(i, 1)))
               call readoutqm(env, fname, env%geolevel, 'bhess', fout, pattern, edum, failed)
               if (failed) then
                  call prepqm(env, fname, env%geolevel, 'bhess', jobcall, fout, pattern, cleanupcall, there, .false.)
                  njobs = njobs + 1
                  call append_char(dirs, trim(fragdirs(i, 1)))
               end if
               call chdir(trim(env%path))
            end if
         end do

         write (*, *) "Starting ", njobs, " restart SPH calculations in parallel at ", trim(env%geolevel), " level"
         call omp_samejobcall(njobs, dirs, jobcall)

         ! readout bhess
         do i = 1, npairs
            if (index(fragdirs(i, 3), 'p') .ne. 0) then
               fname = 'fragment.xyz'
               call chdir(trim(fragdirs(i, 2)))
               ! TODO FIXME if geolevel is DFT this is correct ?? It should be correct
               call readoutqm(env, fname, env%geolevel, 'bhess', fout, pattern, rrho_f1, failed)
               if (failed) then ! if still failed give up
                  e_pairs(i) = 10000.0_wp ! jsut ridiculous high energy to sort this out later
                  call chdir(trim(env%path))
                  cycle
               end if
               call chdir(trim(env%path))
               call chdir(trim(fragdirs(i, 3)))
               call readoutqm(env, fname, env%geolevel, 'bhess', fout, pattern, rrho_f2, failed)
               if (failed) then ! give up
                  e_pairs(i) = 10000.0_wp ! just ridiculous high energy to sort this failed calculation out later
                  call chdir(trim(env%path))
                  cycle
               end if
               call chdir(trim(env%path))
               e_rrho(i) = rrho_f1 + rrho_f2
               e_pairs(i) = e_pairs(i) + e_rrho(i)
            else
               fname = 'isomer.xyz'
               call chdir(trim(fragdirs(i, 1)))
               call readoutqm(env, fname, env%geolevel, 'bhess', fout, pattern, e_rrho(i), failed)
               if (failed) then
                  e_pairs(i) = 10000.0_wp ! jsut ridiculous high energy to sort this out later
                  call chdir(trim(env%path))
                  cycle
               end if
               call chdir(trim(env%path))
               e_pairs(i) = e_pairs(i) + e_rrho(i)
            end if
            drrho(i) = (e_rrho(i) - rrho_start)*autoev
            de(i) = de(i) + drrho(i)
            call wrshort_real(trim(fragdirs(i, 1))//"/de_"//trim(env%tslevel), de(i))
            call wrshort_real(trim(fragdirs(i, 1))//"/drrho_"//trim(env%geolevel), drrho(i))
         end do
         call omp_samejobcall(njobs0, dirs0, cleanupcall, .false.)
      end if

      call rdshort_real("sumreac_"//trim(env%tslevel), sumreac0)

      do i = 1, npairs
         sumreac = sumreac0 + de(i)
         call wrshort_real(trim(fragdirs(i, 1))//"/sumreac_"//trim(env%tslevel), sumreac)
         if (index(fragdirs(i, 3), 'p') .ne. 0) then ! just copy file also to fragments ! TODO FIXME can be done better
            call wrshort_real(trim(fragdirs(i, 2))//"/sumreac_"//trim(env%tslevel), sumreac)
            call wrshort_real(trim(fragdirs(i, 3))//"/sumreac_"//trim(env%tslevel), sumreac)
         end if
      end do
      deallocate (dirs)
      deallocate (e_rrho, drrho, de)

   end subroutine calcreactionenergy

   !> optimize all fragments at charge 0 and input charge for relaxed IPs
   !> fragments can change topology upon optimization, but we tolerate this here
   !> dissociated structures upon optimization are sorted out
   !TODO implement topology check for duplicates upon optimization
   subroutine optfragments(env, npairs_in, npairs_out, fragdirs_in, fragdirs_out)
      implicit none
      character(len=80) :: fname
      integer, intent(in) :: npairs_in
      integer, intent(out) :: npairs_out
      character(len=80)  :: fragdirs_in(npairs_in, 3)
      character(len=80), allocatable, intent(out)   :: fragdirs_out(:, :)
      character(len=80) :: fout, pattern ! fout and pattern for readout
      character(len=1024) :: jobcall, cleanupcall
      character(len=80), allocatable :: dirs(:), dirs0(:)
      character(len=80)   ::  subdir
      type(runtypedata) :: env
      real(wp) :: edum
      real(wp) :: ip
      integer :: i, j, k, chrg1, io
      integer :: njobs, njobs0
      integer :: nat
      logical :: failed, there, dissociated0, dissociated1
      logical, allocatable :: changed(:, :)

      npairs_out = npairs_in

      fname = 'fragment.xyz'

      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']
      fname = 'fragment.xyz'

      write (subdir, '(a,i0)') "charge", 0

      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            do j = 2, 3
               call chdir(trim(fragdirs_in(i, j)))
               io = makedir(trim(subdir))
               call copysub(fname, subdir)
               call rdshort_int(fname, nat)
               if (nat .gt. 1) then
                  call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .false., env%chrg)
                  if (.not. there) then
                     call append_char(dirs, trim(fragdirs_in(i, j)))
                     njobs = njobs + 1
                  end if
                  ! charge 0 optimization
                  call chdir(subdir)
                  call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .false., 0)
                  if (.not. there) then
                     call append_char(dirs, trim(fragdirs_in(i, j))//'/charge0')
                     njobs = njobs + 1
                  end if
               end if
               call chdir(trim(env%path))
            end do
         end if
      end do

      call setompthreads(env, njobs)

      ! actually prepare jobs
      deallocate (dirs)
      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']
      fname = 'fragment.xyz'

      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            do j = 2, 3
               call rdshort(trim(fragdirs_in(i, j))//'/ip_'//trim(env%iplevel), ip)
               call chdir(trim(fragdirs_in(i, j)))
               call rdshort_int(fname, nat)
               if (nat .gt. 1) then
                  call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .false., env%chrg)
                  if (.not. there) then
                     call append_char(dirs, trim(fragdirs_in(i, j)))
                     njobs = njobs + 1
                  end if
                  ! charge 0 optimization
                  call chdir(subdir)
                  call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .false., 0)
                  if (.not. there) then
                     call append_char(dirs, trim(fragdirs_in(i, j))//'/charge0')
                     njobs = njobs + 1
                  end if
               end if
               call chdir(trim(env%path))
            end do
         end if
      end do

      write (*, *) "starting ", njobs, trim(env%geolevel), " optimizations in parallel"
      call omp_samejobcall(njobs, dirs, jobcall)

      ! TODO save here dirs as dirs0 and njobs0 and use this for cleanupcall
      dirs0 = dirs
      njobs0 = njobs

      ! restart calculations
      deallocate (dirs)
      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']
      fname = 'fragment.xyz'

      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            do j = 2, 3
               call rdshort(trim(fragdirs_in(i, j))//'/ip_'//trim(env%iplevel), ip)
               call chdir(trim(fragdirs_in(i, j)))
               call rdshort_int(fname, nat)
               if (nat .gt. 1) then
                  call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
                  if (failed) then
                     call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .true., env%chrg)
                     call append_char(dirs, trim(fragdirs_in(i, j)))
                     njobs = njobs + 1
                  end if
                  ! charge 0 optimization
                  call chdir(subdir)
                  call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
                  if (failed) then
                     call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .true., 0)
                     call append_char(dirs, trim(fragdirs_in(i, j))//'/charge0')
                     njobs = njobs + 1
                  end if
               end if
               call chdir(trim(env%path))
            end do
         end if
      end do

      call setompthreads(env, njobs)

      ! actually prepare jobs
      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            do j = 2, 3
               call rdshort(trim(fragdirs_in(i, j))//'/ip_'//trim(env%iplevel), ip)
               if (env%restartrun .and. abs(ip) .gt. 1.0e-2 .and. chrg1 .eq. 0) cycle ! TODO CHECKME, this is to avoid again optimization at charge0
               ! TODO if we want to change the ip level, but not the geo level, we want to copy here the geometry from the charge0 dir
               call chdir(trim(fragdirs_in(i, j)))
               call rdshort_int(fname, nat)
               if (nat .gt. 1) then
                  call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
                  if (failed) then
                     call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .true., env%chrg)
                     call append_char(dirs, trim(fragdirs_in(i, j)))
                     njobs = njobs + 1
                  end if
                  ! charge 0 optimization
                  call chdir(subdir)
                  call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
                  if (failed) then
                     call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .true., 0)
                     call append_char(dirs, trim(fragdirs_in(i, j))//'/charge0')
                     njobs = njobs + 1
                  end if
               end if
               call chdir(trim(env%path))
            end do
         end if
      end do

      write (*, *) "starting ", njobs, trim(env%geolevel), "restarted optimizations in parallel"
      call omp_samejobcall(njobs, dirs, jobcall)

      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            do j = 2, 3
               call rdshort(trim(fragdirs_in(i, j))//'/ip_'//trim(env%iplevel), ip)
               if (env%restartrun .and. abs(ip) .gt. 1.0e-2 .and. chrg1 .eq. 0) cycle ! TODO CHECKME, this is to avoid again optimization at charge0
               ! TODO if we want to change the ip level, but not the geo level, we want to copy here the geoemtry from the charge0 dir
               call chdir(trim(fragdirs_in(i, j)))
               call rdshort_int(fname, nat)
               if (nat .gt. 1) then
                  call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)

                  if (failed) then
                     write (*, *) "ERROR: optimization of fragment ", trim(fragdirs_in(i, j)), " failed, remove fragment"
                     npairs_out = npairs_out - 1
                     fragdirs_in(i, 1) = ''
                     call chdir(trim(env%path))
                     exit
                  end if
                  call isdissociated(env, fname, dissociated1)
                  ! charge 0 optimization
                  call chdir(subdir)
                  call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
                  if (failed) then
                     write (*, *) "ERROR: optimization of fragment ", trim(fragdirs_in(i, j)), " failed, remove fragment"
                     npairs_out = npairs_out - 1
                     fragdirs_in(i, 1) = ''
                     call chdir(trim(env%path))
                     exit
                  end if
                  call isdissociated(env, fname, dissociated0)

                  if (dissociated1 .and. dissociated0) then
                     write (*, *) "fragment ", trim(fragdirs_in(i, j)), " dissociated upon both optimizations, remove fragment"
                     npairs_out = npairs_out - 1
                     fragdirs_in(i, 1) = ''
                     call chdir(trim(env%path))
                     exit
                  elseif (dissociated0 .and. .not. dissociated1) then
                     write (*, *) "fragment ", trim(fragdirs_in(i, j)) &
                     &  , " dissociated upon neutrally charged optimization, take charged geometry instead"
                     call copy("../"//fname, "./"//fname)
                  elseif (dissociated1 .and. .not. dissociated0) then
                     write (*, *) "fragment ", trim(fragdirs_in(i, j)) &
                     & , " dissociated upon charged optimization, take neutrally charged geometry instead"
                     call copy(fname, "../"//fname)
                  end if
               end if
               call chdir(trim(env%path))
            end do
         end if
      end do

      call omp_samejobcall(njobs0, dirs0, cleanupcall, .false.)
      call sortout0elements(npairs_in, npairs_out, fragdirs_in, fragdirs_out, env%removedirs)

   end subroutine optfragments

   !> optimize fragments at given charge
   !> fragments can change topology upon optimization, but we tolerate this here
   !> dissociated structures upon optimization are sorted out
   !TODO implement topology check for duplicates upon optimization
   subroutine optfragmentsold(env, npairs_in, npairs_out, fragdirs_in, fragdirs_out, chrg)
      implicit none
      character(len=80) :: fname
      integer, intent(in) :: npairs_in
      integer, intent(out) :: npairs_out
      character(len=80)  :: fragdirs_in(npairs_in, 3)
      character(len=80), allocatable, intent(out)   :: fragdirs_out(:, :)
      character(len=80) :: fout, pattern ! fout and pattern for readout
      character(len=1024) :: jobcall, cleanupcall
      character(len=80), allocatable :: dirs(:), dirs0(:)
      type(runtypedata) :: env
      real(wp) :: edum
      real(wp) :: ip
      integer :: i, j, k, chrg1, io
      integer :: njobs, njobs0
      integer :: nat
      integer, intent(in), optional :: chrg
      logical :: failed, there, dissociated
      logical, allocatable :: changed(:, :)

      npairs_out = npairs_in
      if (present(chrg)) then
         chrg1 = chrg
      else
         chrg1 = env%chrg
      end if

      fname = 'fragment.xyz'

      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']
      fname = 'fragment.xyz'

      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            do j = 2, 3
               call rdshort(trim(fragdirs_in(i, j))//'/ip_'//trim(env%iplevel), ip)
               if (env%restartrun .and. abs(ip) .gt. 1.0e-2 .and. chrg1 .eq. 0) cycle ! TODO CHECKME, this is to avoid again optimization at charge0
               ! TODO if we want to change the ip level, but not the geo level, we want to copy here the geoemtry from the charge0 dir
               call chdir(trim(fragdirs_in(i, j)))
               call rdshort_int(fname, nat)
               if (nat .gt. 1) then
                  call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .false., chrg1)
                  if (.not. there) then
                     call append_char(dirs, trim(fragdirs_in(i, j)))
                     njobs = njobs + 1
                  end if
               end if
               call chdir(trim(env%path))
            end do
         end if
      end do

      call setompthreads(env, njobs)

      ! actually prepare jobs
      deallocate (dirs)
      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']
      fname = 'fragment.xyz'

      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            do j = 2, 3
               call rdshort(trim(fragdirs_in(i, j))//'/ip_'//trim(env%iplevel), ip)
               if (env%restartrun .and. abs(ip) .gt. 1.0e-2 .and. chrg1 .eq. 0) cycle ! TODO CHECKME, this is to avoid again optimization at charge0
               ! TODO if we want to change the ip level, but not the geo level, we want to copy here the geoemtry from the charge0 dir
               call chdir(trim(fragdirs_in(i, j)))
               call rdshort_int(fname, nat)
               if (nat .gt. 1) then
                  call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .false., chrg1)
                  if (.not. there) then
                     call append_char(dirs, trim(fragdirs_in(i, j)))
                     njobs = njobs + 1
                  end if
               end if
               call chdir(trim(env%path))
            end do
         end if
      end do

      write (*, *) "starting ", njobs, trim(env%geolevel), " optimizations in parallel"
      call omp_samejobcall(njobs, dirs, jobcall)

      ! TODO save here dirs as dirs0 and njobs0 and use this for cleanupcall
      dirs0 = dirs
      njobs0 = njobs

      ! restart calculations
      deallocate (dirs)
      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']
      fname = 'fragment.xyz'

      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            do j = 2, 3
               call rdshort(trim(fragdirs_in(i, j))//'/ip_'//trim(env%iplevel), ip)
               if (env%restartrun .and. abs(ip) .gt. 1.0e-2 .and. chrg1 .eq. 0) cycle ! TODO CHECKME, this is to avoid again optimization at charge0
               ! TODO if we want to change the ip level, but not the geo level, we want to copy here the geoemtry from the charge0 dir
               call chdir(trim(fragdirs_in(i, j)))
               call rdshort_int(fname, nat)
               if (nat .gt. 1) then
                  call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
                  if (failed) then
                     call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .true., chrg1)
                     call append_char(dirs, trim(fragdirs_in(i, j)))
                     njobs = njobs + 1
                  end if
               end if
               call chdir(trim(env%path))
            end do
         end if
      end do

      call setompthreads(env, njobs)

      ! actually prepare jobs
      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            do j = 2, 3
               call rdshort(trim(fragdirs_in(i, j))//'/ip_'//trim(env%iplevel), ip)
               if (env%restartrun .and. abs(ip) .gt. 1.0e-2 .and. chrg1 .eq. 0) cycle ! TODO CHECKME, this is to avoid again optimization at charge0
               ! TODO if we want to change the ip level, but not the geo level, we want to copy here the geoemtry from the charge0 dir
               call chdir(trim(fragdirs_in(i, j)))
               call rdshort_int(fname, nat)
               if (nat .gt. 1) then
                  call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
                  if (failed) then
                     call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .true., chrg1)
                     call append_char(dirs, trim(fragdirs_in(i, j)))
                     njobs = njobs + 1
                  end if
               end if
               call chdir(trim(env%path))
            end do
         end if
      end do

      write (*, *) "starting ", njobs, trim(env%geolevel), "restarted optimizations in parallel"
      call omp_samejobcall(njobs, dirs, jobcall)

      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            do j = 2, 3
               call rdshort(trim(fragdirs_in(i, j))//'/ip_'//trim(env%iplevel), ip)
               if (env%restartrun .and. abs(ip) .gt. 1.0e-2 .and. chrg1 .eq. 0) cycle ! TODO CHECKME, this is to avoid again optimization at charge0
               ! TODO if we want to change the ip level, but not the geo level, we want to copy here the geoemtry from the charge0 dir
               call chdir(trim(fragdirs_in(i, j)))
               call rdshort_int(fname, nat)
               if (nat .gt. 1) then
                  call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
                  call isdissociated(env, fname, dissociated)
                  if (failed) write (*, *) "ERROR: optimization of fragment ", trim(fragdirs_in(i, j)), " failed, remove fragment"
                  if (dissociated) write (*, *) "fragment ", trim(fragdirs_in(i, j)) &
                  & , " dissociated upon optimization, remove fragment"
                  if (failed .or. dissociated) then
                     npairs_out = npairs_out - 1
                     fragdirs_in(i, 1) = ''
                     call chdir(trim(env%path))
                     exit
                  end if
               end if
               call chdir(trim(env%path))
            end do
         end if
      end do

      call omp_samejobcall(njobs0, dirs0, cleanupcall, .false.)

      ! TODO sort out duplicate fragments, with double loop
      ! i.e. if two fragment pairs have the same toplogy for both fragments, remove one of them

      ! failed are removed
      ! now we could sort out duplicated fragments
      ! topology of fragment can change upon optimization at chrg 0 ... ! what should we do here?
      ! reoptimize at env%chrg

      !  ! topology check necessary? not really, it is far more important for isomers as they can rearrange back to starting structure
      !  ! TODO FIXME: maybe we should add this anyway
      !  !if (env%topocheck .eq. "inchi" .or. env%topocheck .eq. "molbar" ) then
      !  !    call checkproducttopo(env,npairs_in,fragdirs_in,.false.)
      !  !end if
      call sortout0elements(npairs_in, npairs_out, fragdirs_in, fragdirs_out, env%removedirs)

   end subroutine optfragmentsold

   ! copy correct geometry in directory with correct qmdata file
   ! and optimize isomers, TODO maybe we should move this into different routine for better parallelization
   subroutine optproducts(env, npairs_in, npairs_out, fragdirs_in, fragdirs_out)
      implicit none
      character(len=80) :: dir, fname
      integer, intent(in) :: npairs_in
      integer, intent(out) :: npairs_out
      integer :: nfrags
      character(len=80)  :: fragdirs_in(npairs_in, 3)
      character(len=80), allocatable, intent(out)   :: fragdirs_out(:, :)
      character(len=80) :: fout, pattern ! fout and pattern for readout
      character(len=1024) :: jobcall, cleanupcall
      character(len=80), allocatable :: dirs(:), dirs0(:)
      real(wp) :: edum
      type(runtypedata) :: env
      integer :: i, j, k, io
      integer :: nat, chrg
      integer :: njobs, njobs0
      logical :: failed, there, dissociated

      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      write (*, *)
      write (*, *) "Optimizing products at assigned charge"
      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            fname = 'fragment.xyz'
            do j = 2, 3
               call rdshort_int(trim(fragdirs_in(i, j))//'/.CHRG', chrg)
               ! else the correct geometry should be in the directory
               if (chrg .eq. 0) then
                  call chdir(trim(fragdirs_in(i, j)))
                  call copy("charge0/fragment.xyz", "./fragment.xyz")
                  call copy("charge0/qmdata", "./qmdata")
                  call chdir(trim(env%path))
               end if
            end do
         else ! isomers
            fname = 'isomer.xyz'
            ! no check necessary, has to have chrg = 1 or env%chrg
            call chdir(trim(fragdirs_in(i, 1)))
            call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .false., env%chrg)
            if (.not. there) then
               call append_char(dirs, trim(fragdirs_in(i, 1)))
               njobs = njobs + 1
            end if
            !DEL write(*,*) "dir is", trim(dirs(j))
            call chdir(trim(env%path))
         end if
      end do

      call setompthreads(env, njobs)
      ! now actually prepare jobs with correct number of threads
      deallocate (dirs)
      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            cycle
         else ! isomers
            fname = 'isomer.xyz'
            ! no check necessary, has to have chrg = 1 or env%chrg
            call chdir(trim(fragdirs_in(i, 1)))
            call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .false., env%chrg)
            if (.not. there) then
               call append_char(dirs, trim(fragdirs_in(i, 1)))
               njobs = njobs + 1
            end if
            !DEL write(*,*) "dir is", trim(dirs(j))
            call chdir(trim(env%path))
         end if
      end do

      write (*, *) "starting ", njobs, trim(env%geolevel), " optimizations in parallel"
      call omp_samejobcall(njobs, dirs, jobcall)
      njobs0 = njobs
      dirs0 = dirs

      ! restart failed calculations
      deallocate (dirs)
      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      ! TODO this entire nasty loop with condition checks has to be reworked
      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            cycle
         else ! isomers
            fname = 'isomer.xyz'
            ! no check necessary, has to have chrg = 1 or env%chrg
            call chdir(trim(fragdirs_in(i, 1)))
            call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
            if (failed) then ! retry
               call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .true.)
               njobs = njobs + 1
               call append_char(dirs, trim(fragdirs_in(i, 1)))
            end if
            call chdir(trim(env%path))
         end if
      end do

      call setompthreads(env, njobs)

      ! now actually prepare jobs
      deallocate (dirs)
      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      ! TODO this entire nasty loop with condition checks has to be reworked
      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            cycle
         else ! isomers
            fname = 'isomer.xyz'
            ! no check necessary, has to have chrg = 1 or env%chrg
            call chdir(trim(fragdirs_in(i, 1)))
            call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
            if (failed) then ! retry
               call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .true.)
               njobs = njobs + 1
               call append_char(dirs, trim(fragdirs_in(i, 1)))
            end if
            call chdir(trim(env%path))
         end if
      end do

      write (*, *) "starting ", njobs, trim(env%geolevel), " restarted optimizations in parallel"
      call omp_samejobcall(njobs, dirs, jobcall)

      ! now sort out failed and dissociated calculations ! what about topo check with fragments?
      npairs_out = npairs_in

      ! TODO this entire nasty loop with condition checks has to be reworked
      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            cycle
         else ! isomers
            fname = 'isomer.xyz'
            ! no check necessary, has to have chrg = 1 or env%chrg
            call chdir(trim(fragdirs_in(i, 1)))
            call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
            call isdissociated(env, fname, dissociated)
            if (failed .or. dissociated) then ! failed (means no scf convergence) or dissociated sort out
               npairs_out = npairs_out - 1
               fragdirs_in(i, 1) = ''
               call chdir(trim(env%path))
               cycle
            end if
            call chdir(trim(env%path))
         end if
      end do

      ! lets check only isomers here for topology, if they rearrange back
      ! TODO also check for fragments, if they could potentially become to duplicates after reoptimization
      ! (not really realistic and not that problematic like first case)
      if (env%topocheck .eq. "inchi" .or. env%topocheck .eq. "molbar") then
         call checkproducttopo(env, npairs_in, fragdirs_in, .true.)
      end if

      ! now cleanup calculations
      call omp_samejobcall(njobs0, dirs0, cleanupcall, .false.)

      call sortout0elements(npairs_in, npairs_out, fragdirs_in, fragdirs_out, env%removedirs)

   end subroutine optproducts

   ! reoptimize fragments after charge assignment
! and also isomers for reaction energy step
! if hotip everything has to be reoptimized
! PROBLEM isomers can also change topology, (H-shifts)
! TODO topo check for optimizations
   subroutine optproductsold(env, npairs_in, npairs_out, fragdirs_in, fragdirs_out)
      implicit none
      character(len=80) :: dir, fname
      integer, intent(in) :: npairs_in
      integer, intent(out) :: npairs_out
      integer :: nfrags
      character(len=80)  :: fragdirs_in(npairs_in, 3)
      character(len=80), allocatable, intent(out)   :: fragdirs_out(:, :)
      character(len=80) :: fout, pattern ! fout and pattern for readout
      character(len=1024) :: jobcall, cleanupcall
      character(len=80), allocatable :: dirs(:), dirs0(:)
      real(wp) :: edum
      type(runtypedata) :: env
      integer :: i, j, k, io
      integer :: nat, chrg
      integer :: njobs, njobs0
      logical :: failed, there, dissociated

      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            fname = 'fragment.xyz'
            if (env%hotip) then ! both fragments have to be reoptimized
               do j = 2, 3
                  call rdshort_int(trim(fragdirs_in(i, j))//"/fragment.xyz", nat)
                  call rdshort_int(trim(fragdirs_in(i, j))//'/.CHRG', chrg)
                  if (nat .gt. 1) then
                     call chdir(trim(fragdirs_in(i, j)))
                     call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .false., chrg)
                     if (.not. there) then
                        call append_char(dirs, trim(fragdirs_in(i, j)))
                        njobs = njobs + 1
                     end if
                     call chdir(trim(env%path))
                  end if
               end do
            else
               do j = 2, 3
                  call rdshort_int(trim(fragdirs_in(i, j))//"/fragment.xyz", nat)
                  call rdshort_int(trim(fragdirs_in(i, j))//'/.CHRG', chrg)
                  if (nat .gt. 1 .and. chrg .ne. 0) then
                     call chdir(trim(fragdirs_in(i, j)))
                     call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .false., chrg)
                     if (.not. there) then
                        call append_char(dirs, trim(fragdirs_in(i, j)))
                        njobs = njobs + 1
                     end if
                     call chdir(trim(env%path))
                  end if
               end do
            end if
         else ! isomers
            fname = 'isomer.xyz'
            ! no check necessary, has to have chrg = 1 or env%chrg
            call chdir(trim(fragdirs_in(i, 1)))
            call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .false., env%chrg)
            if (.not. there) then
               call append_char(dirs, trim(fragdirs_in(i, 1)))
               njobs = njobs + 1
            end if
            !DEL write(*,*) "dir is", trim(dirs(j))
            call chdir(trim(env%path))
         end if
      end do

      call setompthreads(env, njobs)
      ! now actually prepare jobs with correct number of threads
      deallocate (dirs)
      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            fname = 'fragment.xyz'
            if (env%hotip) then ! both fragments have to be reoptimized
               do j = 2, 3
                  call rdshort_int(trim(fragdirs_in(i, j))//"/fragment.xyz", nat)
                  call rdshort_int(trim(fragdirs_in(i, j))//'/.CHRG', chrg)
                  if (nat .gt. 1) then
                     call chdir(trim(fragdirs_in(i, j)))
                     call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .false., chrg)
                     if (.not. there) then
                        call append_char(dirs, trim(fragdirs_in(i, j)))
                        njobs = njobs + 1
                     end if
                     call chdir(trim(env%path))
                  end if
               end do
            else
               do j = 2, 3
                  call rdshort_int(trim(fragdirs_in(i, j))//"/fragment.xyz", nat)
                  call rdshort_int(trim(fragdirs_in(i, j))//'/.CHRG', chrg)
                  if (nat .gt. 1 .and. chrg .ne. 0) then
                     call chdir(trim(fragdirs_in(i, j)))
                     call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .false., chrg)
                     if (.not. there) then
                        call append_char(dirs, trim(fragdirs_in(i, j)))
                        njobs = njobs + 1
                     end if
                     call chdir(trim(env%path))
                  end if
               end do
            end if
         else ! isomers
            fname = 'isomer.xyz'
            ! no check necessary, has to have chrg = 1 or env%chrg
            call chdir(trim(fragdirs_in(i, 1)))
            call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .false., env%chrg)
            if (.not. there) then
               call append_char(dirs, trim(fragdirs_in(i, 1)))
               njobs = njobs + 1
            end if
            !DEL write(*,*) "dir is", trim(dirs(j))
            call chdir(trim(env%path))
         end if
      end do

      write (*, *) "starting ", njobs, trim(env%geolevel), " optimizations in parallel"
      call omp_samejobcall(njobs, dirs, jobcall)
      njobs0 = njobs
      dirs0 = dirs

      ! restart failed calculations
      deallocate (dirs)
      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      ! TODO this entire nasty loop with condition checks has to be reworked
      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            fname = 'fragment.xyz'
            if (env%hotip) then ! both fragments have to be reoptimized
            do j = 2, 3
               call rdshort_int(trim(fragdirs_in(i, j))//"/fragment.xyz", nat)
               if (nat .gt. 1) then
                  call chdir(trim(fragdirs_in(i, j)))
                  call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
                  if (failed) then ! retry
                     call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .true.)
                     njobs = njobs + 1
                     call append_char(dirs, trim(fragdirs_in(i, j)))
                  end if
                  call chdir(trim(env%path))
               end if
            end do
            else
            do j = 2, 3
               call rdshort_int(trim(fragdirs_in(i, j))//"/fragment.xyz", nat)
               call rdshort_int(trim(fragdirs_in(i, j))//'/.CHRG', chrg)
               if (nat .gt. 1 .and. chrg .ne. 0) then
                  !if (nat .gt. 1 .and. chrg .ne. env%chrg) then
                  call chdir(trim(fragdirs_in(i, j)))
                  call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
                  if (failed) then ! retry
                     call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .true.)
                     njobs = njobs + 1
                     call append_char(dirs, trim(fragdirs_in(i, j)))
                  end if
                  call chdir(trim(env%path))
               end if
            end do
            end if
         else ! isomers
            fname = 'isomer.xyz'
            ! no check necessary, has to have chrg = 1 or env%chrg
            call chdir(trim(fragdirs_in(i, 1)))
            call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
            if (failed) then ! retry
               call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .true.)
               njobs = njobs + 1
               call append_char(dirs, trim(fragdirs_in(i, 1)))
            end if
            call chdir(trim(env%path))
         end if
      end do

      call setompthreads(env, njobs)

      ! now actually prepare jobs
      deallocate (dirs)
      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      ! TODO this entire nasty loop with condition checks has to be reworked
      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            fname = 'fragment.xyz'
            if (env%hotip) then ! both fragments have to be reoptimized
            do j = 2, 3
               call rdshort_int(trim(fragdirs_in(i, j))//"/fragment.xyz", nat)
               if (nat .gt. 1) then
                  call chdir(trim(fragdirs_in(i, j)))
                  call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
                  if (failed) then ! retry
                     call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .true.)
                     njobs = njobs + 1
                     call append_char(dirs, trim(fragdirs_in(i, j)))
                  end if
                  call chdir(trim(env%path))
               end if
            end do
            else
            do j = 2, 3
               call rdshort_int(trim(fragdirs_in(i, j))//"/fragment.xyz", nat)
               call rdshort_int(trim(fragdirs_in(i, j))//'/.CHRG', chrg)
               if (nat .gt. 1 .and. chrg .ne. 0) then
                  !if (nat .gt. 1 .and. chrg .ne. env%chrg) then
                  call chdir(trim(fragdirs_in(i, j)))
                  call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
                  if (failed) then ! retry
                     call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .true.)
                     njobs = njobs + 1
                     call append_char(dirs, trim(fragdirs_in(i, j)))
                  end if
                  call chdir(trim(env%path))
               end if
            end do
            end if
         else ! isomers
            fname = 'isomer.xyz'
            ! no check necessary, has to have chrg = 1 or env%chrg
            call chdir(trim(fragdirs_in(i, 1)))
            call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
            if (failed) then ! retry
               call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .true.)
               njobs = njobs + 1
               call append_char(dirs, trim(fragdirs_in(i, 1)))
            end if
            call chdir(trim(env%path))
         end if
      end do

      write (*, *) "starting ", njobs, trim(env%geolevel), " restarted optimizations in parallel"
      call omp_samejobcall(njobs, dirs, jobcall)

      ! now sort out failed and dissociated calculations ! what about topo check with fragments?
      npairs_out = npairs_in

      ! TODO this entire nasty loop with condition checks has to be reworked
      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            fname = 'fragment.xyz'
            if (env%hotip) then ! both fragments have to be reoptimized
            do j = 2, 3
               call rdshort_int(trim(fragdirs_in(i, j))//"/fragment.xyz", nat)
               if (nat .gt. 1) then
                  call chdir(trim(fragdirs_in(i, j)))
                  call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
                  call isdissociated(env, fname, dissociated)
                  if (failed .or. dissociated) then ! failed (means no scf convergence) or dissociated sort out
                     npairs_out = npairs_out - 1
                     do k = 1, 3
                        fragdirs_in(i, k) = '' ! here all 3 for topocheck routine later
                     end do
                     call chdir(trim(env%path))
                     exit
                  end if
               end if
            end do
            else
            do j = 2, 3
               call rdshort_int(trim(fragdirs_in(i, j))//"/fragment.xyz", nat)
               call rdshort_int(trim(fragdirs_in(i, j))//'/.CHRG', chrg)
               if (nat .gt. 1 .and. chrg .ne. 0) then
                  !if (nat .gt. 1 .and. chrg .ne. env%chrg) then
                  call chdir(trim(fragdirs_in(i, j)))
                  call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
                  call isdissociated(env, fname, dissociated)
                  if (failed .or. dissociated) then ! failed (means no scf convergence) or dissociated sort out
                     npairs_out = npairs_out - 1
                     do k = 1, 3
                        fragdirs_in(i, k) = '' ! here all 3 for topocheck routine later
                     end do
                     call chdir(trim(env%path))
                     exit
                  end if
                  call chdir(trim(env%path))
               end if
            end do
            end if
         else ! isomers
            fname = 'isomer.xyz'
            ! no check necessary, has to have chrg = 1 or env%chrg
            call chdir(trim(fragdirs_in(i, 1)))
            call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed)
            call isdissociated(env, fname, dissociated)
            if (failed .or. dissociated) then ! failed (means no scf convergence) or dissociated sort out
               npairs_out = npairs_out - 1
               fragdirs_in(i, 1) = ''
               call chdir(trim(env%path))
               cycle
            end if
            call chdir(trim(env%path))
         end if
      end do

      ! lets check only isomers here for topology, if they rearrange back
      ! TODO also check for fragments, if they could potentially become to duplicates after reoptimization
      ! (not really realistic and not that problematic like first case)
      if (env%topocheck .eq. "inchi" .or. env%topocheck .eq. "molbar") then
         call checkproducttopo(env, npairs_in, fragdirs_in, .true.)
      end if

      ! now cleanup calculations
      call omp_samejobcall(njobs0, dirs0, cleanupcall, .false.)

      call sortout0elements(npairs_in, npairs_out, fragdirs_in, fragdirs_out, env%removedirs)

   end subroutine optproductsold

!> compute topology isomers, compare them, and sort out duplicates
!> only do this for isomers, as it is more crucial to do this here
!> TODO implement also for fragments
   subroutine checkproducttopo(env, npairs, fragdirs, checkisos)
      implicit none
      integer, intent(in) :: npairs
      character(len=80), intent(inout)  :: fragdirs(npairs, 3)
      character(len=512), allocatable :: topocodes(:, :) ! barcodes of structures
      character(len=80) :: ftopo
      character(len=1024) :: jobcall, cleanupcall, jobcall2
      ! character(len=1024) :: jobcall2 ! TODO check when this is necessary, atoms with two symbols and second is capital, but why
      ! does this happen?
      character(len=80) :: dir, fname
      character(len=80), allocatable :: dirs(:)
      type(runtypedata) :: env
      integer :: i, j, io
      integer ::  njobs
      logical, optional :: checkisos ! only isomers are checked
      logical :: ex

      write (*, *) "Check product topology after optimization"
      ftopo = 'topo'
      ! first get topo of input file
      inquire (file='isomer.xyz', exist=ex)
      if (ex) then
         fname = 'isomer.xyz'
      else
         inquire (file='fragment.xyz', exist=ex)
         if (ex) then
            fname = 'fragment.xyz'
         else
            write (*, *) "ERROR: no fragment or isomer file found"
            stop
         end if
      end if
      ! Molbar has problems with capital letters in atom names
      call rewrite_xyz_elements(trim(fname))

      if (env%topocheck .eq. "molbar") then
         write (jobcall, '(a)') "molbar "//trim(fname)//" > "//trim(ftopo)//" 2> /dev/null"
      elseif (env%topocheck .eq. "inchi") then
         write (jobcall, '(a)') "obabel -i xyz "//trim(fname)//" -o inchi >  "//trim(ftopo)//" 2> /dev/null"
      end if
      call execute_command_line(trim(jobcall))

      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      ! first isomer loop
      fname = 'isomer.xyz'
      if (env%topocheck .eq. "molbar") write (jobcall, '(a)') "molbar "//trim(fname)//" > "//trim(ftopo)//" 2> /dev/null"
      if (env%topocheck .eq. "inchi") write (jobcall, '(a)') &
      & "obabel -i xyz "//trim(fname)//" -o inchi >  "//trim(ftopo)//" 2> /dev/null"
      do i = 1, npairs
         if (index(fragdirs(i, 3), 'p') .ne. 0) cycle
         call append_char(dirs, trim(fragdirs(i, 1)))
         njobs = njobs + 1
      end do

      ! TODO solve this
      if (.not. checkisos) then
         write (*, *) "product topology check of fragments not implemented yet"
      end if

      ! prepare topocalls
      call setompthreads(env, njobs) !for molbar this is fine, only for ORCA threads have to be correctly set for prepqm call
      call omp_samejobcall(njobs, dirs, jobcall)

      ! now read out topologies
      if (checkisos) then
         call compareisotopos(env, npairs, fragdirs)
      else
         write (*, *) "product topology check of fragments not implemented yet"
         error stop
         call compareproducttopos(env, npairs, fragdirs)
      end if

      write (cleanupcall, '(a)') "rm topo"
      if (env%printlevel .lt. 3) then
         call omp_samejobcall(njobs, dirs, cleanupcall, .false.)
      end if

   end subroutine checkproducttopo

   ! TODO not only previous fragments but also all grand parents and
   subroutine compareisotopos(env, npairs, fragdirs)
      implicit none
      integer, intent(in) :: npairs
      character(len=80), intent(inout)  :: fragdirs(npairs, 3)
      character(len=512), allocatable :: topocodes(:, :) ! barcodes of structures
      character(len=80) :: ftopo
      character(len=512) :: grampatopo ! topo of gradnparent strucutre, dont want fragment to react back to initial structure
      character(len=:), allocatable :: path, grampadir, startdir
      character(len=80) :: fname
      character(len=1024) :: jobcall, jobcall2
      type(runtypedata) :: env
      integer :: i, j, io
      integer ::  nfragl ! fragmentation level
      integer :: ilen
      logical, allocatable :: double(:), iso(:)
      logical :: ex

      ftopo = 'topo'
      ilen = len_trim(env%startdir)
      path = trim(env%path)
      startdir = path(ilen + 2:)

      ilen = len(startdir)
      nfragl = 1
      ! level of fragmentation sequence  p1f1 is 1 meaning level 2 of fragmentation
      do i = 1, ilen
         if (startdir(i:i) == 'p') then
            nfragl = nfragl + 1
         end if
      end do

      grampadir = ""
      ! read also out topo of start fragment prevent fragmentation back to "parent parent structure"
      if (nfragl .gt. 1) then
         do i = ilen, 1, -1
            if (startdir(i:i) == 'p') then
               grampadir = startdir(1:i - 1)
               exit
            end if
         end do
         call rdshort_string(trim(env%startdir)//"/"//grampadir//"/"//trim(ftopo), grampatopo)
         if (env%topocheck .eq. "molbar") call cuttopology(grampatopo)
      end if

      allocate (double(0:npairs), iso(0:npairs), topocodes(0:npairs, 2))
      write (topocodes, '(a)') ''

      call rdshort_string(trim(ftopo), topocodes(0, 1)) ! read out startfragment barcode
      if (env%topocheck .eq. "molbar") call cuttopology(topocodes(0, 1))
      if (topocodes(0, 1) == '') then
         ! first get topo of input file
         inquire (file='isomer.xyz', exist=ex)
         if (ex) then
            fname = 'isomer.xyz'
         else
            inquire (file='fragment.xyz', exist=ex)
            if (ex) then
               fname = 'fragment.xyz'
            else
               write (*, *) "ERROR: no fragment or isomer file found"
               stop
            end if
         end if
         write (*, *) " topo file of start fragment is empty, retry"
         if (env%topocheck .eq. "molbar") then
            ! TODO change this in to in code routine, as molbar has problems with capital letters
            ! do this as molbar has problems if atoms are written in capital letters (e.g.,CL instead of Cl)
            call rewrite_xyz_elements(trim(fname))

            write (jobcall, '(a)') "molbar "//trim(fname)//" > "//trim(ftopo)//" 2> /dev/null"
         elseif (env%topocheck .eq. "inchi") then
            write (jobcall, '(a)') "obabel -i xyz "//trim(fname)//" -o inchi >  "//trim(ftopo)//" 2> /dev/null"
         end if
         call execute_command_line(trim(jobcall))
         call rdshort_string(trim(ftopo), topocodes(0, 1)) ! read out startfragment barcode
         if (env%topocheck .eq. "molbar") call cuttopology(topocodes(0, 1))
      end if

      ! we compare it with all isomers
      do i = 1, npairs
         if (index(fragdirs(i, 3), 'p') .ne. 0) then
            iso(i) = .false.
         else
            call chdir(trim(fragdirs(i, 1)))
            call rdshort_string(trim(ftopo), topocodes(i, 1)) !
            if (env%topocheck .eq. "molbar") call cuttopology(topocodes(i, 1))
            if (topocodes(i, 1) == '') then
               write (*, *) " topo file of ", trim(fragdirs(i, 1)), " is empty"
            end if
            call chdir(trim(env%path))
            iso(i) = .true.
         end if
      end do
      ! identify duplicates
      double = .false.
      do i = 0, npairs
         if (double(i)) cycle
         do j = i + 1, npairs
            if (.not. double(j)) then
               if (iso(j)) then
                  if (trim(topocodes(i, 1)) .ne. '' .and. trim(topocodes(i, 1)) == trim(topocodes(j, 1))) then
                     double(j) = .true.
                     if (i .gt. 0) write (*, *) trim(fragdirs(j, 1)), " is a duplicate of ", trim(fragdirs(i, 1))
                     fragdirs(j, 1) = ''
                     write (*, *) j, " is a duplicate of ", i

                  elseif (trim(topocodes(j, 1)) .ne. '' .and. trim(topocodes(j, 1)) == trim(grampatopo)) then ! also dont allow backreaction
                     double(j) = .true.
                     fragdirs(j, 1) = ''
                     write (*, *) j, "is a duplicate of its pre-precursor"
                     write (*, *) trim(fragdirs(j, 1)), " is a duplicate of its pre-precursor"
                  end if
               else
                  cycle
               end if
            end if
         end do
      end do
      deallocate (double, iso, topocodes)
   end subroutine compareisotopos

   ! TODO not only previous fragments but also all grand parent ions
   subroutine compareproducttopos(env, npairs, fragdirs)
      implicit none
      integer, intent(in) :: npairs
      character(len=80), intent(inout)  :: fragdirs(npairs, 3)
      character(len=512), allocatable :: topocodes(:, :) ! barcodes of structures
      character(len=80) :: ftopo
      character(len=512) :: grampatopo ! topo of gradnparent strucutre, dont want fragment to react back to initial structure
      character(len=:), allocatable :: path, grampadir, startdir
      type(runtypedata) :: env
      integer :: i, j, io
      integer ::  nfragl ! fragmentation level
      integer :: ilen
      logical, allocatable :: double(:), iso(:)

      ftopo = 'topo'
      ilen = len_trim(env%startdir)
      path = trim(env%path)
      startdir = path(ilen + 2:)

      ilen = len(startdir)
      nfragl = 1
      ! level of fragmentation sequence  p1f1 is 1 meaning level 2 of fragmentation
      do i = 1, ilen
         if (startdir(i:i) == 'p') then
            nfragl = nfragl + 1
         end if
      end do

      grampadir = ""
      ! read also out topo of start fragment prevent fragmentation back to "parent parent structure"
      if (nfragl .gt. 1) then
         do i = ilen, 1, -1
            if (startdir(i:i) == 'p') then
               grampadir = startdir(1:i - 1)
               exit
            end if
         end do
         call rdshort_string(trim(env%startdir)//"/"//grampadir//"/"//trim(ftopo), grampatopo)
         if (env%topocheck .eq. "molbar") call cuttopology(grampatopo)
      end if

      allocate (double(0:npairs), iso(0:npairs), topocodes(0:npairs, 2))
      topocodes = ''

      call rdshort_string(trim(ftopo), topocodes(0, 1)) ! read out startfragment barcode
      if (env%topocheck .eq. "molbar") call cuttopology(topocodes(0, 1))
      if (topocodes(0, 1) == '') then
         write (*, *) " topo file of start fragment is empty"
      end if
      ! we compare it with all isomers
      do i = 1, npairs
         call chdir(trim(fragdirs(i, 1)))
         call rdshort_string(trim(ftopo), topocodes(i, 1)) !
         if (env%topocheck .eq. "molbar") call cuttopology(topocodes(i, 1))
         if (topocodes(i, 1) == '') then
            write (*, *) " topo file of ", trim(fragdirs(i, 1)), " is empty"
         end if
         call chdir(trim(env%path))
      end do
      ! identify duplicates
      double = .false.
      do i = 0, npairs
         if (double(i)) cycle
         do j = i + 1, npairs
            if (.not. double(j)) then
               if (topocodes(i, 1) .ne. '' .and. trim(topocodes(i, 1)) == trim(topocodes(j, 1))) then
                  double(j) = .true.
                  fragdirs(j, 1) = ''
                  write (*, *) j, " is a duplicate of ", i
               elseif (topocodes(j, 1) .ne. '' .and. trim(topocodes(j, 1)) == trim(grampatopo)) then ! also dont allow backreaction
                  double(j) = .true.
                  fragdirs(j, 1) = ''
                  write (*, *) j, "is a duplicate of its pre-precursor"
               end if
            end if
         end do
      end do
   end subroutine compareproducttopos

! compare topos also for fragments
   subroutine comparefragtopos(env, npairs, fragdirs, checkisos)
      implicit none
      integer, intent(in) :: npairs
      character(len=80), intent(inout)  :: fragdirs(npairs, 3)
      character(len=512), allocatable :: topocodes(:, :) ! barcodes of structures
      character(len=512) :: grampatopo ! topo of gradnparent strucutre, dont want fragment to react back to initial structure
      character(len=80) :: grampadir
      character(len=100) :: ftopo
      type(runtypedata) :: env
      integer :: i, j, io
      integer ::  nfragl ! fragmentation level
      integer :: ilen
      logical, allocatable :: double(:), iso(:)
      logical, optional :: checkisos ! if all is set, all isomers are checked, otherwise only fragments are checked

      nfragl = 1
      ilen = len_trim(fragdirs(1, 1))
      ! level of fragmentation sequence  p1f1 is 1 meaning level 2 of fragmentation
      do i = 1, ilen
         if (fragdirs(1, 1) (i:i) == 'p') then
            nfragl = nfragl + 1
         end if
      end do

      allocate (double(0:npairs), iso(0:npairs), topocodes(0:npairs, 2))
      ftopo = 'topo'
      call rdshort_string(trim(ftopo), topocodes(0, 1)) ! read out startfragment barcode
      call cuttopology(topocodes(0, 1))
      ! read also out topo of start fragment prevent fragmentation back to "parent parent structure"
      if (nfragl .gt. 1) then
         do i = ilen, 1
            if (fragdirs(1, 1) (i:i) == 'p') then
               grampadir = trim(fragdirs(1, 1) (1:i - 1))
               exit
            end if
         end do
         ilen = len_trim(grampadir)
         do i = ilen, 1
            if (grampadir(i:i) == 'p') then
               grampadir = trim(grampadir(1:i - 1))
               exit
            end if
         end do
         call rdshort_string(trim(env%startdir)//"/"//trim(grampadir)//trim(ftopo), grampatopo)
         call cuttopology(grampatopo)
      end if

      ! we compare it with all isomers

      topocodes = 'null'
      do i = 1, npairs
         if (index(fragdirs(i, 3), 'p') .ne. 0) then
            call chdir(trim(fragdirs(i, 2)))
            call rdshort_string(trim(ftopo), topocodes(i, 1))
            call cuttopology(topocodes(i, 1))
            if (topocodes(i, 1) == '') then
               write (*, *) " topo file of ", trim(fragdirs(i, 2)), " is empty"
            end if
            call chdir(trim(env%path))
            call chdir(trim(fragdirs(i, 3)))
            call rdshort_string(trim(ftopo), topocodes(i, 2))
            call cuttopology(topocodes(i, 2))
            if (topocodes(i, 2) == '') then
               write (*, *) " topo file of ", trim(fragdirs(i, 3)), " is empty"
            end if
            call chdir(trim(env%path))
            iso(i) = .false.
         else
            call chdir(trim(fragdirs(i, 1)))
            call rdshort_string(trim(ftopo), topocodes(i, 1)) !
            call cuttopology(topocodes(i, 1))
            if (topocodes(i, 1) == '') then
               write (*, *) " topo file of ", trim(fragdirs(i, 1)), " is empty"
            end if
            call chdir(trim(env%path))
            iso(i) = .true.
         end if
      end do
      ! identify duplicates
      double = .false.
      do i = 0, npairs
         if (double(i)) cycle
         do j = i + 1, npairs
            if (.not. double(j)) then
               if (iso(j)) then
                  if (checkisos) then
                  if (topocodes(i, 1) .ne. '' .and. trim(topocodes(i, 1)) == trim(topocodes(j, 1))) then
                     double(j) = .true.
                     fragdirs(j, 1) = ''
                     write (*, *) j, " is a duplicate of ", i
                  elseif (trim(topocodes(i, 1)) == trim(grampatopo)) then ! also dont allow backreaction
                     double(j) = .true.
                     fragdirs(j, 1) = ''
                     write (*, *) j, "is a duplicate of its pre-precursor"
                  end if
                  end if
               else
                  if (topocodes(i, 1) .ne. '') then
                  if (trim(topocodes(i, 1)) == trim(topocodes(j, 1)) &
                  & .and. trim(topocodes(i, 2)) == trim(topocodes(j, 2)) &
                  & .or. trim(topocodes(i, 1)) == trim(topocodes(j, 2)) &
                  & .and. trim(topocodes(i, 2)) == trim(topocodes(j, 1))) then
                     double(j) = .true.
                     fragdirs(j, 1) = ''
                     write (*, *) j, " is a duplicate of ", i
                  end if
                  end if
               end if
            end if
         end do
      end do
   end subroutine comparefragtopos

! get topo of fname with molbar or inchi
   subroutine calctopo(env, fname)
      implicit none
      type(runtypedata) :: env
      character(len=80) :: fname
      character(len=100) :: ftopo
      character(len=1024) :: jobcall
      integer :: nat

      call rdshort_int(fname, nat)
      if (nat .eq. 1) return ! molbar not needed for single atoms

      write (ftopo, '(a)') "topo"
      call rewrite_xyz_elements(trim(fname))
      if (env%topocheck .eq. "molbar") then
         write (jobcall, '(a)') "molbar "//trim(fname)//" > "//trim(ftopo)//" 2> /dev/null"
      elseif (env%topocheck .eq. "inchi") then
         write (jobcall, '(a)') "obabel -i xyz "//trim(fname)//" -o inchi >  "//trim(ftopo)//" 2> /dev/null"
      end if
      !DEL swrite(*,*) "jobcall is", trim(jobcall)
      call execute_command_line(trim(jobcall))
   end subroutine calctopo

!> optimize all msreact products with GFN-FF and check again their topology
!> we only want to optimize here isomers and not fragment pairs !!!!!!!
!> this can 1.): prevent back-rearrangments of H-shifts, and 2.): prevent duplicates
!> which are not detected by molbar (e.g. other hybridizations of nitrogen)
!> TODO we could also do a double loop check with fragments and not only pairs
   subroutine gfftopocheck(env, npairs_in, npairs_out, fragdirs_in, fragdirs_out)
      implicit none
      character(len=80) :: dir, fname
      integer, intent(in) :: npairs_in
      integer, intent(out) :: npairs_out
      integer :: nfrags
      character(len=80)  :: fragdirs_in(npairs_in, 3)
      character(len=80), allocatable, intent(out)   :: fragdirs_out(:, :)
      character(len=80)  :: ftopo
      character(len=80) :: fout, pattern ! fout and pattern for readout
      character(len=1024) :: jobcall, cleanupcall
      character(len=80), allocatable :: dirs(:)
      logical :: there
      real(wp) :: edum
      type(runtypedata) :: env
      integer :: i, j
      integer :: nat, chrg
      integer :: njobs

      njobs = npairs_in

      ftopo = 'topo'

      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      fname = 'isomer.xyz'
      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) cycle !! TODO CHECK THIS
         call chdir(trim(fragdirs_in(i, 1)))
         njobs = njobs + 1
         call prepqm(env, fname, 'gff', 'opt', jobcall, fout, pattern, cleanupcall, there, .false.)
         call append_char(dirs, trim(fragdirs_in(i, 1)))
         call chdir(trim(env%path))
      end do

      write (*, *) "starting ", njobs, ' GFN-FF ', "optimizations in parallel before topology check"
      call setompthreads(env, njobs) ! for xtb threads have not to be correctly set for prepqm call, only for ORCA
      call omp_samejobcall(njobs, dirs, jobcall)

      ! lets just assume gff never fails an optimization
      ! todo maybe include check here anyways

      !cleanup
      call omp_samejobcall(njobs, dirs, cleanupcall, .false.)

      deallocate (dirs)
      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']
      ! determine topology of isomers
      fname = 'isomer.xyz'
      if (env%topocheck .eq. "molbar") write (jobcall, '(a)') "molbar "//trim(fname)//" > "//trim(ftopo)//" 2> /dev/null"
      if (env%topocheck .eq. "inchi") write (jobcall, '(a)') &
      & "obabel -i xyz "//trim(fname)//" -o inchi >  "//trim(ftopo)//" 2> /dev/null"
      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) cycle
         call append_char(dirs, trim(fragdirs_in(i, 1)))
         njobs = njobs + 1
      end do

      ! prepare topocalls
      call setompthreads(env, njobs) !for molbar this is fine, only for ORCA threads have to be correctly set for prepqm call
      write (*, *) "starting ", njobs, ' topo ', "calculations in parallel:"
      call omp_samejobcall(njobs, dirs, jobcall)

      ! lets check only isomers here for topology, if they rearrange back
      ! TODO also check for fragments, if they could potentially become to duplicates after reoptimization
      ! (not really realistic and not that problematic like first case)
      if (env%topocheck .eq. "inchi" .or. env%topocheck .eq. "molbar") then
         !call checkproducttopo(env, npairs_in, fragdirs_in, .false.) ! TODO CHECKME why should we check for pairs here??
         call checkproducttopo(env, npairs_in, fragdirs_in, .true.)
      end if
      call sortout0elements(npairs_in, npairs_out, fragdirs_in, fragdirs_out, env%removedirs)

   end subroutine gfftopocheck
end module reaction
