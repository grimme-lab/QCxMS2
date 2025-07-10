!> this module handels the transition state search and barrier calculations
module tsmod
   use iomod
   use qcxms2_data
   use qmmod
   use xtb_mctc_convert
   use utility
   use rmsd_ls
   !use rmsd, only: get_rmsd
   implicit none
contains

   !>  TODO implement other path search methods
   !>  CURRENTLY only for NEB with geodesic and without geodesic
   subroutine tssearch(env, fname, npairs_in, npairs_out2, fragdirs_in, fragdirs_out2)
      implicit none
      type(timer):: tim
      type(runtypedata) :: env
      character(len=80) :: reac, prod
      character(len=1024) :: jobcall, copycall, cleanupcall
      character(len=80), allocatable :: dirs(:), dirs0(:)
      character(len=80) :: job, fout, pattern, fout2, pattern2 ! fout and pattern for readout of qm data
      character(len=80) :: fragdirs_in(:, :)
      character(len=80), allocatable :: fragdirs_out(:, :)
      character(len=80), allocatable, intent(out)   :: fragdirs_out2(:, :)
      character(len=80) :: query, fname
      integer, intent(in) :: npairs_in
      integer :: npairs_out
      integer, intent(out) :: npairs_out2
      integer :: i, j, njobs, njobs0
      integer, allocatable :: nmode(:), nmode2(:)
      integer :: chrg, uhf, spin
      integer :: io
      logical :: ex, ldum, found, failed, there, conv
      logical, allocatable :: tsopt(:)
      character(len=80) :: pwd
      character(len=80) :: sumform
      logical :: hbonddiss
      real(wp), allocatable :: e_ts(:), e_rrho(:)
      real(wp) ::  e_start, rrho_start, drrho
      real(wp) :: ircmode !internal reaction coordinate at transition state
      real(wp) :: barrier ! barrier height
      real(wp) :: barrier_tddft, e_start_tddft ! if we include tddt excited states
      real(wp), allocatable :: e_ts_tddft(:)
      real(wp), allocatable :: e_ts_2uhf(:), e_ts_tddft_2uhf(:)
      real(wp) :: edum
      integer :: uhf1, uhf2
      real(wp) :: t1, w1, t2, w2

      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      !prepare and run path search for all pairs
      write (*, *) "Search for transition states for", npairs_in, " pairs"
      call setompthreads(env, npairs_in)

      ! get charge and multiplicity from reactant and perform path search at this state
      inquire (file='.CHRG', exist=ex)
      if (.not. ex) then
         write (*, *) "Error: Charge file not found, just using input charge ", env%chrg
         chrg = env%chrg
      else
         call rdshort_int('.CHRG', chrg)
      end if

      ! determine correct spin state
      ! always compute path at multiplicity of reactant
      inquire (file='.UHF', exist=ex)
         if (.not. ex) then
            write (*, *) "Error: UHF file not found, just computing uhf from input charge ", env%chrg
            call printspin(env, fname, spin, chrg)
            uhf = spin - 1
            if (uhf == -2) uhf = 0 ! for H+
         else
            call rdshort_int('.UHF', uhf)
         end if
         do i = 1, npairs_in
            call wrshort_int(trim(env%path)//"/"//trim(fragdirs_in(i, 1))//'/.UHF', uhf)
         end do


      ! for homolytic bond dissociation we have to sum up the uhf values of the fragments
      if (env%pathmult) then
         do i = 1, npairs_in
           if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
               call rdshort_int(trim(env%path)//"/"//trim(fragdirs_in(i, 2))//"/.UHF", uhf1)
               call rdshort_int(trim(env%path)//"/"//trim(fragdirs_in(i, 3))//"/.UHF", uhf2)
               uhf = uhf1 + uhf2
               if (uhf .gt. 1) then
                  write(*,*) "computing path of ",trim(env%path)//"/"//trim(fragdirs_in(i, 1))," with UHF ",uhf
               end if
               call wrshort_int(trim(env%path)//"/"//trim(fragdirs_in(i, 1))//'/.UHF', uhf)
           end if
         end do 
      end if


      

      ! prepare start and end geometries
      do i = 1, npairs_in
         call copy(fname, trim(fragdirs_in(i, 1))//'/start.xyz')
         call chdir(trim(fragdirs_in(i, 1)))
         inquire (file='pair.xyz', exist=ex)
         if (ex) prod = 'pair.xyz'
         inquire (file='isomer.xyz', exist=ex)
         if (ex) prod = 'isomer.xyz'
         call copy(trim(prod), "end.xyz")
         call wrshort_int('.CHRG', chrg)
         
         call chdir(trim(env%path))
      end do

      !> prepare geodesic inteprolation guess for neb
      !> for NEB we
      !> 1. compute with default guess
      !> 2. if SCF not converged we restart with fermi smearing
      !> 3. if NEB not converged we restart with geodesic interpolation
      !> 4. if ORCA not converged we restart with fermi smearing and geodesic interpolation
      if (env%tsgeodesic) then
         do i = 1, npairs_in
            call chdir(trim(fragdirs_in(i, 1)))
            inquire (file='interpolated.xyz', exist=ex)
            if (.not. ex) then
               call prepgeodesic(env, jobcall)
               call append_char(dirs, trim(fragdirs_in(i, 1)))
               njobs = njobs + 1
            end if
            call chdir(trim(env%path))
         end do
         call omp_samejobcall(njobs, dirs, jobcall)
      end if

      ! determine number of jobs
      njobs = 0
      do i = 1, npairs_in
         call chdir(trim(fragdirs_in(i, 1)))
         if (trim(env%tsfinder) == 'neb') then
            inquire (file='neb_finished', exist=ex)
            if (.not. ex) njobs = njobs + 1
          elseif (trim(env%tsfinder) == 'crest_tsgen') then
            inquire (file='crest_tsgen_finished', exist=ex)
            if (.not. ex) njobs = njobs + 1
         elseif (trim(env%tsfinder) == 'gsm') then
            inquire (file='gsm_finished', exist=ex)
            if (.not. ex) njobs = njobs + 1   
         else
            write (*, *) "Warning: other path searches than NEB not supported, stopping QCxMS2"
            stop
         !!!  !xtbpath
         !!!elseif (trim(env%tsfinder) == 'xtb') then
         !!!   inquire (file='xtb_finished', exist=ex)
         !!!   if (.not. ex) njobs = njobs + 1
         end if
         call chdir(trim(env%path))
      end do

      ! prepare calculations
      call setompthreads(env, njobs)
      deallocate (dirs)
      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      ! prepare calculations and conduct them later in parallel
      do i = 1, npairs_in

         hbonddiss = .false.   ! special settings for R-H bond dissociation, currently not used
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
         do j = 2, 3
            call getsumform(trim(env%path)//"/"//trim(fragdirs_in(i, j))//"/fragment.xyz", sumform)
            if (sumform == "H1") hbonddiss = .true.
         end do
         end if
         call chdir(trim(fragdirs_in(i, 1)))

         !NEB with ORCA
         if (trim(env%tsfinder) == 'neb') then
            inquire (file='neb_finished', exist=ex)
            if (.not. ex) then
               call prepneb(env, jobcall, hbonddiss, .false., .false.)
               njobs = njobs + 1
               call append_char(dirs, trim(fragdirs_in(i, 1)))
            end if
         end if

         !xtbpath
         if (trim(env%tsfinder) == 'xtb') then
            inquire (file='xtb_finished', exist=ex)
            if (.not. ex) then
               call prepxtbpath(env, jobcall)
               njobs = njobs + 1
               call append_char(dirs, trim(fragdirs_in(i, 1)))
            end if
         end if
         !crest tsguess 
         if (trim(env%tsfinder) == 'crest_tsgen') then
            inquire (file='crest_tsgen_finished', exist=ex)
            if (.not. ex) then
               call prepcrest_tsgen(env, jobcall)
               njobs = njobs + 1
               call append_char(dirs, trim(fragdirs_in(i, 1)))
            end if
         end if
 
          !xtbpath
         if (trim(env%tsfinder) == 'gsm') then
            inquire (file='gsm_finished', exist=ex)
            if (.not. ex) then
               call prepgsm(env, jobcall,hbonddiss, .false., .false.)
               njobs = njobs + 1
               call append_char(dirs, trim(fragdirs_in(i, 1)))
            end if
         end if

         call chdir(trim(env%path))
      end do
      !start TS search in parallel

      call timing(t1, w1)

      write (*, *) "Starting ", njobs, " "//trim(env%tsfinder)//" transition state searches"

      call setomptoone()
      call omp_samejobcall(njobs, dirs, jobcall)

      call timing(t2, w2)
      env%tpath = env%tpath + (w2 - w1)
      write (*,'(a,f10.1,a)')  "Time spent for path search and barrier calculation: ", (w2 - w1), "s"

      if (env%tsfinder == 'crest_tsgen') then
         do i = 1, npairs_in
            call chdir(trim(fragdirs_in(i, 1))) 
            ! restore .UHF and .CHRG deleted by crest
            !call move("uhftemp", ".UHF")
            !call move("chrgtemp", ".CHRG")
            ! TODO implement cleanup here
            call chdir(trim(env%path))
         end do
      end if



      if (trim(env%tsfinder) == 'neb') then
         njobs = 0
         deallocate (dirs)
         allocate (dirs(1)) ! first element always empty
         dirs = ['']

         do i = 1, npairs_in
            call chdir(trim(fragdirs_in(i, 1)))
            call minigrep('orca.out', 'THE NEB OPTIMIZATION HAS CONVERGED', conv)
            if (.not. conv) njobs = njobs + 1
            call chdir(trim(env%path))
         end do

         if (env%tsgeodesic) then
            ! prepare restart calculations
            call setompthreads(env, njobs)
            deallocate (dirs)
            njobs = 0
            allocate (dirs(1)) ! first element always empty
            dirs = ['']
            do i = 1, npairs_in
               call chdir(trim(fragdirs_in(i, 1)))
               call minigrep('orca.out', 'THE NEB OPTIMIZATION HAS CONVERGED', conv)
               if (.not. conv) then
                  ! retry path search
                  hbonddiss = .false.
                  njobs = njobs + 1
                  call append_char(dirs, trim(fragdirs_in(i, 1)))
                  if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
                  do j = 2, 3
                     call getsumform(trim(env%path)//"/"//trim(fragdirs_in(i, j))//"/fragment.xyz", sumform)
                     if (sumform == "H1") hbonddiss = .true.
                  end do
                  end if
                  call prepneb(env, jobcall, hbonddiss, .false., .true.)
               end if
               call chdir(trim(env%path))
            end do
            call setomptoone()
            write (*, *) "Retry NEB searches for ", njobs, " pairs"
            call omp_samejobcall(njobs, dirs, jobcall)
         end if
      end if
      ! for restart of path search

      ! check how many transition state searches were succesfull
      ! combine with pickts
      njobs = 0
      do i = 1, npairs_in
         call chdir(trim(fragdirs_in(i, 1)))
         ! control path search and take highest point on path as transition state
         call pickts(env, found) ! gsm neb and xtb should work
         call checkpath(env, failed)
         call chdir(trim(env%path))
         if (failed) njobs = njobs + 1
      end do

      call setompthreads(env, njobs)
      !retry path search
      deallocate (dirs)
      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      do i = 1, npairs_in
         call chdir(trim(fragdirs_in(i, 1)))
         ! control path search and take highest point on path as transition state
         call pickts(env, found) ! gsm neb and xtb should work
         call checkpath(env, failed)
         if (failed) then
            ! retry path search
            hbonddiss = .false.
            njobs = njobs + 1
            call append_char(dirs, trim(fragdirs_in(i, 1)))
            if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            do j = 2, 3
               call getsumform(trim(env%path)//"/"//trim(fragdirs_in(i, j))//"/fragment.xyz", sumform)
               if (sumform == "H1") hbonddiss = .true.
            end do
            end if
            if (trim(env%tsfinder) == 'neb') then
               call prepneb(env, jobcall, hbonddiss, .true., .false.)
            elseif (trim(env%tsfinder) == 'xtb') then   ! todo tune this
               call prepxtbpath(env, jobcall)
            end if
         end if
         call chdir(trim(env%path))
      end do

      if (njobs .gt. 0) write (*, *) "Retry transition state searches for ", njobs, " pairs"

      ! Set OMP parallelization of xtb call in ORCA to 1
      call setomptoone()
      call omp_samejobcall(njobs, dirs, jobcall)

      ! and now check again if neb converged now with fermi smearing at 5000K
      if (trim(env%tsfinder) == 'neb') then
         njobs = 0
         deallocate (dirs)
         allocate (dirs(1)) ! first element always empty
         dirs = ['']

         do i = 1, npairs_in
            call chdir(trim(fragdirs_in(i, 1)))
            call minigrep('orca.out', 'THE NEB OPTIMIZATION HAS CONVERGED', conv)
            if (.not. conv) njobs = njobs + 1
            call chdir(trim(env%path))
         end do

         if (env%tsgeodesic) then
            ! prepare restart calculations
            call setompthreads(env, njobs)
            deallocate (dirs)
            njobs = 0
            allocate (dirs(1)) ! first element always empty
            dirs = ['']
            do i = 1, npairs_in
               call chdir(trim(fragdirs_in(i, 1)))
               call minigrep('orca.out', 'THE NEB OPTIMIZATION HAS CONVERGED', conv)

               if (.not. conv) then
                  hbonddiss = .false.
                  njobs = njobs + 1
                  call append_char(dirs, trim(fragdirs_in(i, 1)))
                  if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
                  do j = 2, 3
                     call getsumform(trim(env%path)//"/"//trim(fragdirs_in(i, j))//"/fragment.xyz", sumform)
                     if (sumform == "H1") hbonddiss = .true.
                  end do
                  end if
                  call prepneb(env, jobcall, hbonddiss, .true., .true.)
               end if
               call chdir(trim(env%path))
            end do

            call setomptoone()
            if (njobs .gt. 0) write (*, *) "Retry NEB again at 5000 K searches for ", njobs, " pairs"
            call omp_samejobcall(njobs, dirs, jobcall)
         end if
      end if

      ! check how many transition state searches were succesfull
      do i = 1, npairs_in
         call chdir(trim(fragdirs_in(i, 1)))
         ! control path search and take highest point on path as transition state
         call pickts(env, found) ! gsm neb and xtb should work
         if (.not. found) then
            write (*, *) "Search for transition state for pair ", trim(fragdirs_in(i, 1)), " failed, is sorted out"
            fragdirs_in(i, 1) = ''
         end if
         call chdir(trim(env%path))
      end do

      ! cleanup NEB runs
      
     
      write (cleanupcall, '(a)') 'rm orca_im*.gbw orca_im*.citations.tmp orca_im*.xtbrestart'
      do i = 1, npairs_in
         call chdir(trim(fragdirs_in(i, 1)))
         call cleanup_nebcalc()
         call execute_command_line(trim(cleanupcall))
         if (env%printlevel .lt. 3) then
            !call remove('orca.out')
            !call remove('orcaerr.out')
            call remove('geoin.xyz')
            call remove('.data')
         end if
         call chdir(trim(env%path))
      end do

      ! rewrite fragdir list so that empty entries are removed
      call sortout0elements(npairs_in, npairs_out, fragdirs_in, fragdirs_out, env%removedirs) ! todo maybe we do this just at the end ...

      allocate (tsopt(npairs_out))

      write (*, *) "Number of succesfull ts searches is ", npairs_out
      if (npairs_in .ne. npairs_out) then
         write (*, *) "Warning: Some TS searches failed check output, corresponding fragments are removed from list"
         if (env%sortoutcascade) then
            write (*, *) "Warning: Some TS searches failed or turned out as cascade reactions &
            &, check output, corresponding fragments are removed from list!"
         end if
      end if

      ! more logical
      pattern = 'THE OPTIMIZATION HAS CONVERGED'
      fout = 'geo.out'

      !if optimization failed, ts.xyz should be highest point on path
      !pattern = 'ORCA-job'
      !fout = 'ts.xyz'
      njobs = npairs_out
      !copy ts in ts DIR
      do i = 1, npairs_out
         call chdir(trim(fragdirs_out(i, 1)))
         ! TODO CHECKME, this is for restarting old version calculations, in which the last geometry of not converged
         ! TS optimizations was taken as TS guess
         inquire (file='ts/ts.xyz', exist=ex)
         if (ex) then
            call chdir("ts")
            ! for restart run we dont want to overwrite optimized ts.xyz
            call readoutqm(env, 'ts.xyz', env%geolevel, 'optts', fout, pattern, edum, failed)
            if (.not. failed) then
               tsopt(i) = .true.
               njobs = njobs - 1
            else
               tsopt(i) = .true.
               ! overwrite not converged TS guess with last geometry of path
               write (*, *) "TS Optimization not converged, we have to take highest point of path as ts"
               call printpwd()
               call copy('../ts.xyz', 'ts.xyz') ! overwrite not converged ts with highest point of traj
               njobs = njobs - 1
            end if
            call chdir("..")
         else
            io = makedir('ts')
            call copysub('ts.xyz', 'ts')

            ! copy charge and multiplicity
            call copysub('.UHF', 'ts')
            call copysub('.CHRG', 'ts')
            tsopt(i) = .false.
            if (.not. env%reoptts) then
                tsopt(i) = .true. ! just set it to true if we do not want to reoptimize (avoid to use orca)
                njobs = 0
            end if
         end if
         call chdir(trim(env%path))
      end do

      call setompthreads(env, njobs)
      deallocate (dirs)
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      do i = 1, npairs_out
         if (tsopt(i)) cycle
         call chdir(trim(fragdirs_out(i, 1)))
         call chdir("ts")
         io = makedir('hess')
         call copysub('ts.xyz', 'hess')
         call copysub("../.CHRG", 'hess')
         call copysub("../.UHF", 'hess')
         call chdir("hess")
         ! do always GFN2-xTB  hessian calculation here, but for DFT geo, we do it in orca so that we can read in the hessian later
         call prepqm(env, 'ts.xyz', env%geolevel, 'hess', jobcall, fout, pattern, cleanupcall, there, .false.)
         call append_char(dirs, trim(fragdirs_out(i, 1))//'/ts/hess')
         call chdir("..")
         call chdir("..")
         call chdir(trim(env%path))
      end do

      write (*, *) "Starting ", njobs, " "//trim(env%geolevel)//" First Hessian calculations"
      call omp_samejobcall(njobs, dirs, jobcall)

      dirs0 = dirs
      njobs0 = njobs

      ! check freq calculations
      njobs = 0

      do i = 1, npairs_out
         if (tsopt(i)) cycle
         call chdir(trim(fragdirs_out(i, 1)))
         call chdir("ts")
         call chdir("hess")
         call readoutqm(env, 'ts.xyz', env%geolevel, 'hess', fout, pattern, edum, failed)
         if (failed) njobs = njobs + 1
         call chdir("..")
         call chdir("..")
         call chdir(trim(env%path))
      end do

      ! restart failed freq calculations
      call setompthreads(env, njobs)
      njobs = 0
      deallocate (dirs)
      allocate (dirs(1)) ! first element always empty
      dirs = ['']
      do i = 1, npairs_out
         if (tsopt(i)) cycle
         call chdir(trim(fragdirs_out(i, 1)))
         call chdir("ts")
         call chdir("hess")
         call readoutqm(env, 'ts.xyz', env%geolevel, 'hess', fout, pattern, edum, failed)
         if (failed) then ! retry
            call prepqm(env, 'ts.xyz', env%geolevel, 'hess', jobcall, fout, pattern, cleanupcall, there, .true.)
            call append_char(dirs, trim(fragdirs_out(i, 1))//'/ts/hess')
            njobs = njobs + 1
         end if
         call chdir("..")
         call chdir("..")
         call chdir(trim(env%path))
      end do

      call omp_samejobcall(njobs, dirs, jobcall)

      njobs = 0
      ! readout frequencies
      allocate (nmode(npairs_out))
      do i = 1, npairs_out
         if (tsopt(i)) cycle
         call chdir(trim(fragdirs_out(i, 1)))
         call chdir("ts")
         call chdir("hess")
         call readoutqm(env, 'ts.xyz', env%geolevel, 'hess', fout, pattern, edum, failed)
         if (failed) then
            write (*, *) "TS Hessian calculation failed, we have to take end as ts ircmode set to 0"
            tsopt(i) = .true.
            ircmode = 0.0_wp
            call chdir("..")
            call chdir("..")
            call wrshort_real('ircmode_'//trim(env%geolevel), ircmode)
            call chdir(trim(env%path))
            cycle
         end if
         call findirc(env, nmode(i), ircmode)
         if (nmode(i) .eq. 0) then
            write (*, *) "No imaginary mode found, we have to take end as ts ircmode set to 0"
            ircmode = 0.0_wp
         else
            njobs = njobs + 1 ! we can optimize the ts guess
         end if
         call chdir("..")
         call chdir("..")
         call wrshort_real('ircmode_'//trim(env%geolevel), ircmode)
         call chdir(trim(env%path))
      end do

      
      ! cleanup
      call omp_samejobcall(njobs0, dirs0, cleanupcall, .false.)

      deallocate (dirs)
      allocate (dirs(1)) ! first element always empty
      dirs = ['']
      call setompthreads(env, njobs)
      ! now optimize TS guess structures

      do i = 1, npairs_out
         if (tsopt(i)) cycle
         if (nmode(i) .le. 0) cycle
         call chdir(fragdirs_out(i, 1))
         call chdir("ts")
         call copy("../.CHRG", ".CHRG")
         call copy("../.UHF", ".UHF")
         call wrshort_int('nmode', nmode(i)) ! TODO rewrite routine for this ...
         if (env%qmprog == "orca") then ! TODO has to be ORCA in all cases ...
            call copy("hess/orca.hess", "orca.hess")
         end if
         call prepqm(env, 'ts.xyz', env%geolevel, 'optts', jobcall, fout, pattern, cleanupcall, there, .false.) ! give nmode here
         call append_char(dirs, trim(fragdirs_out(i, 1))//'/ts')
         call chdir("..")
         call chdir(trim(env%path))
      end do
      write (*, *) "Starting ", njobs, " "//trim(env%geolevel)//" TS optimizations"
      call omp_samejobcall(njobs, dirs, jobcall)

      do i = 1, npairs_out
         if (tsopt(i)) cycle
         if (nmode(i) .le. 0) cycle
         call chdir(trim(fragdirs_out(i, 1)))
         call chdir("ts")
         call readoutqm(env, 'ts.xyz', env%geolevel, 'optts', fout, pattern, edum, failed)
         if (failed) then
               nmode(i) = 0 
               write (*, *) "TS Optimization not converged we have to take last point of path as ts and ircmodeis set to 0"
               call copy("../ts.xyz", "ts.xyz")
               ircmode = 0.0_wp
         end if
         call chdir("..")
         call wrshort_real('ircmode_'//trim(env%geolevel), ircmode)
         call chdir(trim(env%path))
      end do

      ! cleanup ts optimizations
      call omp_samejobcall(njobs, dirs, cleanupcall, .false.)

      !  check for number of second hessian jobs
      njobs = 0
      do i = 1, npairs_out
         if (tsopt(i)) cycle
         if (nmode(i) .le. 0) cycle
         njobs = njobs + 1
      end do
      ! then we control again with findirc mode if we have a transition state
      call setompthreads(env, njobs)
      deallocate (dirs)
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      do i = 1, npairs_out
         if (tsopt(i)) cycle
         if (nmode(i) .le. 0) cycle
         call chdir(trim(fragdirs_out(i, 1)))
         call chdir("ts")
         io = makedir('hess2')
         call copysub('ts.xyz', 'hess2')
         call chdir("hess2")
         call prepqm(env, 'ts.xyz', env%geolevel, 'hess', jobcall, fout, pattern, cleanupcall, there, .false.)
         ! it cannot be already there anyway
         call append_char(dirs, trim(fragdirs_out(i, 1))//'/ts/hess2')
         call chdir("..")
         call chdir("..")
         call chdir(trim(env%path))
      end do
      write (*, *) "Starting ", njobs, " "//trim(env%geolevel)//" Hessian calculations"

      call omp_samejobcall(njobs, dirs, jobcall)

      njobs0 = njobs
      deallocate (dirs0)
      dirs0 = dirs

      ! check freq calculations
      njobs = 0
      do i = 1, npairs_out
         if (tsopt(i)) cycle
         if (nmode(i) .le. 0) cycle
         call chdir(trim(fragdirs_out(i, 1)))
         call chdir("ts")
         call chdir("hess2")
         call readoutqm(env, 'ts.xyz', env%geolevel, 'hess', fout, pattern, edum, failed)
         if (failed) njobs = njobs + 1
         call chdir("..")
         call chdir("..")
         call chdir(trim(env%path))
      end do

      ! restart failed freq calculations
      call setompthreads(env, njobs)
      deallocate (dirs)
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      do i = 1, npairs_out
         if (tsopt(i)) cycle
         if (nmode(i) .le. 0) cycle
         call chdir(trim(fragdirs_out(i, 1)))
         call chdir("ts")
         call chdir("hess2")
         call readoutqm(env, 'ts.xyz', env%geolevel, 'hess', fout, pattern, edum, failed)
         if (failed) then
            call prepqm(env, 'ts.xyz', env%geolevel, 'hess', jobcall, fout, pattern, cleanupcall, there, .true.)
            call append_char(dirs, trim(fragdirs_out(i, 1))//'/ts/hess2')
         end if
         call chdir("..")
         call chdir("..")
         call chdir(trim(env%path))
      end do

      call omp_samejobcall(njobs, dirs, jobcall)

      ! readout frequencies

      allocate (nmode2(npairs_out))
      nmode2 = 0
      do i = 1, npairs_out
         if (tsopt(i)) cycle
         if (nmode(i) .le. 0) cycle
         call chdir(trim(fragdirs_out(i, 1)))
         call chdir("ts")
         call chdir("hess2")
         ! failed calculations are handled in this routine
         call findirc(env, nmode2(i), ircmode)
         if (nmode2(i) .eq. 0) then ! includes also the cases without IRC mode in the first findirc check
            write (*, *) "No imaginary mode found, &
            & we just take highest point on reaction path as transition state"
            call copy("../../ts.xyz", "../ts.xyz") ! copy tsguess from highest point into tsdir overwrite failed ts optimization
         end if
         call chdir("..")
         call chdir("..")
         ! overwrite ircmode with optimized ts if we have still one, otherwise leave old ircmode ...
         if (nmode2(i) .ne. 0) then
            call wrshort_real('ircmode_'//trim(env%geolevel), ircmode)
            call wrshort_int('nmode', nmode2(i)) ! overwrite nmode
         end if
         call chdir(trim(env%path))
      end do

      ! cleanup hessian calculations
      call omp_samejobcall(njobs0, dirs0, cleanupcall, .false.)

      ! NOW COMPUTE BARRIERS
      call setompthreads(env, npairs_out) ! now again with all npair_out
      deallocate (dirs)
      allocate (dirs(1)) ! first element always empty
      dirs = ['']
      njobs = 0

      fname = 'ts.xyz'
      !TDDFT?
      if (env%exstates .gt. 0) then
         job = 'tddft'
      else
         job = 'sp'
      end if

      do i = 1, npairs_out
         call chdir(trim(fragdirs_out(i, 1)))
         call chdir("ts")
         call prepqm(env, fname, env%tslevel, job, jobcall, fout, pattern, cleanupcall, there, .false.)
         if (.not. there) njobs = njobs + 1
         call chdir("..")
         call chdir(trim(env%path))
      end do

      call setompthreads(env, njobs)
      deallocate (dirs)
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      do i = 1, npairs_out
         call chdir(trim(fragdirs_out(i, 1)))
         call chdir("ts")
         call prepqm(env, fname, env%tslevel, job, jobcall, fout, pattern, cleanupcall, there, .false.)
         if (.not. there) then
            call append_char(dirs, trim(fragdirs_out(i, 1))//'/ts')
         end if
         call chdir("..")
         call chdir(trim(env%path))
      end do

      write (*, *) "Starting ", njobs, " "//trim(env%tslevel)//" Singlepoint calculations on TS"
      call omp_samejobcall(njobs, dirs, jobcall)

      deallocate (dirs0)
      dirs0 = dirs
      njobs0 = njobs
      njobs = 0
      ! check failed singlepoints

      do i = 1, npairs_out
         call chdir(trim(fragdirs_out(i, 1)))
         call chdir("ts")
         call readoutqm(env, fname, env%tslevel, job, fout, pattern, edum, failed)
         if (failed) njobs = njobs + 1
         call chdir("..")
         call chdir(trim(env%path))
      end do

      ! actually restart failed singlepoints
      call setompthreads(env, njobs)
      deallocate (dirs)
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      do i = 1, npairs_out
         call chdir(trim(fragdirs_out(i, 1)))
         call chdir("ts")
         call readoutqm(env, fname, env%tslevel, job, fout, pattern, edum, failed)
         if (failed) then
            call prepqm(env, fname, env%tslevel, job, jobcall, fout, pattern, cleanupcall, there, .true.)
            call append_char(dirs, trim(fragdirs_out(i, 1))//'/ts')
         end if
         call chdir("..")
         call chdir(trim(env%path))
      end do

      if (njobs .gt. 0) write (*, *) "Restarting ", njobs, " "//trim(env%tslevel)//" singlepoint calculations on TS"
      call omp_samejobcall(njobs, dirs, jobcall)


      ! readout energies
      allocate (e_ts(npairs_out), e_rrho(npairs_out))
      if (env%exstates .gt. 0) then
         allocate (e_ts_tddft(npairs_out))
      end if
      
      do i = 1, npairs_out
         call chdir(trim(fragdirs_out(i, 1)))
         call chdir("ts")
         if (env%exstates .gt. 0) then
            call readoutqm(env, fname, env%tslevel, job, fout, pattern, e_ts_tddft(i), failed)
            if (failed) then 
               write (*, *) "Singl-point calculation of TS failed, pair ", trim(fragdirs_out(i, 1)), " is sorted out"
               fragdirs_out(i, 1) = ''
            end if
            call readoutqm(env, fname, env%tslevel, 'sp', fout, pattern, e_ts(i), failed)
            if (failed) then 
               write (*, *) "Singl-point calculation of TS failed, pair ", trim(fragdirs_out(i, 1)), " is sorted out"
               fragdirs_out(i, 1) = ''
            end if
         else
            call readoutqm(env, fname, env%tslevel, job, fout, pattern, e_ts(i), failed)
            if (failed) then 
               write (*, *) "Singl-point calculation of TS failed, pair ", trim(fragdirs_out(i, 1)), " is sorted out"
               fragdirs_out(i, 1) = ''
            end if
         end if
         call chdir("..")
         call chdir(trim(env%path))
      end do

      ! cleanup
      call omp_samejobcall(njobs0, dirs0, cleanupcall, .false.)

      !TODO DEL this option is spin forbidden anyways
      ! Check for higher spin states (uhf +2)
      if (env%checkmult) then 
         call setompthreads(env, npairs_out) ! now again with all npair_out
         deallocate (dirs)
         allocate (dirs(1)) ! first element always empty
         dirs = ['']
         njobs = 0

         fname = 'ts.xyz'
         !TDDFT?
         if (env%exstates .gt. 0) then
            job = 'tddft'
         else
            job = 'sp'
         end if

         do i = 1, npairs_out
            call chdir(trim(fragdirs_out(i, 1)))
            call chdir("ts")
            call prepqm(env, fname, env%tslevel, job, jobcall, fout, pattern, cleanupcall, there, .false.,chrg,uhf+2)
            if (.not. there) njobs = njobs + 1
            call chdir("..")
            call chdir(trim(env%path))
         end do

         call setompthreads(env, njobs)
         deallocate (dirs)
         allocate (dirs(1)) ! first element always empty
         dirs = ['']

         do i = 1, npairs_out
            call chdir(trim(fragdirs_out(i, 1)))
            call chdir("ts")
            call prepqm(env, fname, env%tslevel, job, jobcall, fout, pattern, cleanupcall, there, .false.,chrg,uhf+2)
            if (.not. there) then
               call append_char(dirs, trim(fragdirs_out(i, 1))//'/ts')
            end if
            call chdir("..")
            call chdir(trim(env%path))
         end do

         write (*, *) "Starting ", njobs, " "//trim(env%tslevel)//" Singlepoint calculations on TS with UHF+2"
         call omp_samejobcall(njobs, dirs, jobcall)

         deallocate (dirs0)
         dirs0 = dirs
         njobs0 = njobs
         njobs = 0
         ! check failed singlepoints

         do i = 1, npairs_out
            call chdir(trim(fragdirs_out(i, 1)))
            call chdir("ts")
            call readoutqm(env, fname, env%tslevel, job, fout, pattern, edum, failed)
            if (failed) njobs = njobs + 1
            call chdir("..")
            call chdir(trim(env%path))
         end do

         ! actually restart failed singlepoints
         call setompthreads(env, njobs)
         deallocate (dirs)
         allocate (dirs(1)) ! first element always empty
         dirs = ['']

         do i = 1, npairs_out
            call chdir(trim(fragdirs_out(i, 1)))
            call chdir("ts")
            call readoutqm(env, fname, env%tslevel, job, fout, pattern, edum, failed)
            if (failed) then
               call prepqm(env, fname, env%tslevel, job, jobcall, fout, pattern, cleanupcall, there, .true.)
               call append_char(dirs, trim(fragdirs_out(i, 1))//'/ts')
            end if
            call chdir("..")
            call chdir(trim(env%path))
         end do

         if (njobs .gt. 0) write (*, *) "Restarting ", njobs, " "//trim(env%tslevel)//" singlepoint calculations on TS with UHF+2"
         call omp_samejobcall(njobs, dirs, jobcall)

         ! readout energies
         allocate (e_ts_2uhf(npairs_out))
         if (env%exstates .gt. 0) then
            allocate (e_ts_tddft_2uhf(npairs_out))
         end if

         do i = 1, npairs_out
            call chdir(trim(fragdirs_out(i, 1)))
            call chdir("ts")
            if (env%exstates .gt. 0) then
               call readoutqm(env, fname, env%tslevel, job, fout, pattern, e_ts_tddft_2uhf(i), failed)
               if (failed) e_ts_tddft_2uhf(i) = 1000000.0_wp ! just set to high value TODO FIXME
               call readoutqm(env, fname, env%tslevel, 'sp', fout, pattern, e_ts_2uhf(i), failed)
            else
               call readoutqm(env, fname, env%tslevel, job, fout, pattern, e_ts_2uhf(i), failed)
               if (failed) then 
                  write (*, *) "single point calculation calculation of TS failed,reaction ", trim(fragdirs_out(i, 1)), " is sorted out"
                  fragdirs_out(i, 1) = ''
               end if
               if (failed) e_ts_2uhf(i) = 1000000.0_wp ! just set to high value TODO FIXME
            end if
            call chdir("..")
            call chdir(trim(env%path))
         end do
          ! cleanup
         call omp_samejobcall(njobs0, dirs0, cleanupcall, .false.)

         do i = 1, npairs_out 
            if (e_ts_2uhf(i) .lt. e_ts(i)) then
               write (*, *) "TS energy for pair ", trim(fragdirs_out(i, 1)), " is lower for UHF+2 spin state, we take this one"
               e_ts(i) = e_ts_2uhf(i)
               if (env%exstates .gt. 0) e_ts_tddft(i) = e_ts_tddft_2uhf(i)
            else
               ! write old multiplicity in directory 
               call wrshort_int('.UHF',uhf)
            end if
         end do
      end if 

     

      ! determine number of jobs not necessary for bhess calculations
      deallocate (dirs)
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      njobs = 0
      if (env%bhess) then
         do i = 1, npairs_out
            ! compute bhess for if pair is there and we cannot use the hessian
            if ((index(fragdirs_out(i, 1), 'p') .ne. 0 .and. nmode(i) .ne. 1 .and. nmode2(i) .ne. 1) &
            & .or. (env%geolevel .ne. "gfn2" .and. env%geolevel .ne. "gfn2spinpol" .and. env%geolevel .ne. "gfn2_tblite"))  then
               call chdir(trim(fragdirs_out(i, 1)))
               call chdir("ts")
               io = makedir('bhess')
               call copysub('ts.xyz', 'bhess')
               call copysub('.CHRG', 'bhess')
               ! compute bhess at correct spin state ! TODO reoptimization necessary?
               call copysub('.UHF', 'bhess')
               call chdir("bhess")
               call prepqm(env, fname, env%geolevel, 'bhess', jobcall, fout, pattern, cleanupcall, there, .false.)
               if (.not. there) then
                  njobs = njobs + 1
                  call append_char(dirs, trim(fragdirs_out(i, 1))//'/ts/bhess')
               end if
               call chdir("..")
               call chdir("..")
               call chdir(trim(env%path))
            end if
         end do

         call setompthreads(env, njobs)
         write (*, *) "Starting ", njobs, " bhess calculations on TS structures"
         call omp_samejobcall(njobs, dirs, jobcall)

         deallocate (dirs0)
         dirs0 = dirs
         njobs0 = njobs

         ! restart failed bhess calculations

         deallocate (dirs)
         allocate (dirs(1)) ! first element always empty
         dirs = ['']

         njobs = 0
         do i = 1, npairs_out
            if ((index(fragdirs_out(i, 1), 'p') .ne. 0 .and. nmode(i) .ne. 1 .and. nmode2(i) .ne. 1) &
            & .or. (env%geolevel .ne. "gfn2" .and. env%geolevel .ne. "gfn2spinpol" .and. env%geolevel .ne. "gfn2_tblite"))  then
               call chdir(trim(fragdirs_out(i, 1)))
               call chdir("ts")
               io = makedir('bhess')
               call copysub('ts.xyz', 'bhess')
               call chdir("bhess")
               call readoutqm(env, fname, env%geolevel, 'bhess', fout, pattern, edum, failed)
               if (failed) then
                  call prepqm(env, fname, env%geolevel, 'bhess', jobcall, fout, pattern, cleanupcall, there, .true.)
                  njobs = njobs + 1
                  call append_char(dirs, trim(fragdirs_out(i, 1))//'/ts/bhess')
               end if
               call chdir("..")
               call chdir("..")
               call chdir(trim(env%path))
            end if
         end do

         call setompthreads(env, njobs)
         if (njobs .gt. 0) write (*, *) "Restarting ", njobs, " failed bhess calculations on TS structures"
         call omp_samejobcall(njobs, dirs, jobcall)


         fout2 = 'orca.out'
         pattern2 = 'Zero point energy                ...'
         ! readout RRHO
         do i = 1, npairs_out
            if (index(fragdirs_out(i, 1), 'p') .ne. 0) then
               call chdir(trim(fragdirs_out(i, 1)))             
               call chdir("ts")
                  if ((nmode(i) .ne. 1 .and. nmode2(i) .ne. 1) & 
                  & .or. (env%geolevel .ne. "gfn2" .and. env%geolevel .ne. "gfn2spinpol" .and. env%geolevel .ne. "gfn2_tblite")) then 
                  call chdir("bhess")
                  call readoutqm(env, fname, env%geolevel, 'bhess', fout, pattern, e_rrho(i), failed)
                  if (failed) then 
                     write (*, *) "bhess calculation failed, fragment ", trim(fragdirs_out(i, 1)), " is sorted out"
                     fragdirs_out(i, 1) = ''
                  end if
                 
               elseif (nmode2(i) .eq. 1) then
                  call chdir("hess2")
                  call readoutqm(env, fname, env%geolevel, 'hess', fout2, pattern2, e_rrho(i), failed)
                  if (failed) then 
                     write (*, *) "hess2 calculation failed, fragment ", trim(fragdirs_out(i, 1)), " is sorted out"
                     fragdirs_out(i, 1) = ''
                  end if
               elseif (nmode(i) .eq. 1) then
                  call chdir("hess")
                  call readoutqm(env, fname, env%geolevel, 'hess', fout2, pattern2, e_rrho(i), failed)
                  if (failed) then 
                     write (*, *) "hess2 calculation failed, fragment ", trim(fragdirs_out(i, 1)), " is sorted out"
                     fragdirs_out(i, 1) = ''
                  end if      
               end if
               call chdir("..")
               call chdir("..")
               call chdir(trim(env%path))    
            end if
         end do

         ! cleanup
         call omp_samejobcall(njobs0, dirs0, cleanupcall, .false.)

      end if

      call rdshort_int('.CHRG', chrg)
      call rdshort_int('.UHF', uhf)
      ! now compute barriers
      if (env%exstates .gt. 0) then
         write (query, '(a,1x,i0,1x,i0)') trim(env%tslevel)//" tddft", chrg, uhf
         call grepval('qmdata', trim(query), ldum, e_start_tddft)
      end if

      write (query, '(a,1x,i0,1x,i0)') trim(env%tslevel)//" sp", chrg, uhf
      call grepval('qmdata', trim(query), ldum, e_start)
      if (env%bhess) then
         write (query, '(a,1x,i0,1x,i0)') trim(env%geolevel)//" bhess", chrg, uhf
         call grepval('qmdata', trim(query), ldum, rrho_start)
      end if

      do i = 1, npairs_out
         if (index(fragdirs_out(i, 1), 'p') .ne. 0) then
            call chdir(trim(fragdirs_out(i, 1)))
            barrier = (e_ts(i) - e_start)*autoev
            if (env%bhess) then
               drrho = (e_rrho(i) - rrho_start)*autoev
               barrier = barrier + drrho
               !DEL? not really needed just for diagnostice
               call wrshort_real('earrho_'//trim(env%geolevel), drrho)
            end if

            if (env%exstates .gt. 0) then
               barrier_tddft = (e_ts_tddft(i) - e_start_tddft)*autoev
               if (env%bhess) then
                  barrier_tddft = barrier_tddft + drrho
               end if
               call wrshort_real('barrier_'//trim(env%tslevel), barrier_tddft)
               call wrshort_real('barriergs_'//trim(env%tslevel), barrier) ! GS barrier, just for testing
               write (*, *) "Contributions of higher states to barrier:", barrier - barrier_tddft
            else
               call wrshort_real('barrier_'//trim(env%tslevel), barrier)
            end if
            call chdir(trim(env%path))
         end if
      end do

      call sortout0elements( npairs_out,npairs_out2,fragdirs_out,fragdirs_out2,env%removedirs)
   end subroutine tssearch

! read first 9 modes from g98 file, messy routine TODO FIXME
   subroutine rdg98modes(nat, freqs, modes, freqcount)
      implicit none

      character(len=512) :: tmp, dum, str
      real(wp), allocatable :: modes(:, :, :) ! number of mode, number of atom, x, y or z coordinate
      real(wp), allocatable :: freqs(:)
      integer, intent(out) :: freqcount ! number of read imaginary modes
      integer :: remfreq, freqcount2
      integer :: io, io1, io2, ich
      integer :: nat
      integer ::  maxfreq
      integer ::  i, j
      integer :: dum1, dum2
      logical :: ex

      freqcount = 0
      freqcount2 = 0
      maxfreq = 3*nat - 6
      modes = 0.0d0

      open (newunit=ich, file='g98.out')
      do

         read (ich, '(a)', iostat=io) tmp
         if (io < 0) exit
         io1 = index(tmp, 'Frequencies --', .true.)
         if (io1 .gt. 0) then
            do i = 0, 2  ! first 9 modes have to be enough ! extend later to more ...
               io2 = io1 + 14
               if (maxfreq - freqcount .ge. 3) then
                  remfreq = 3
               elseif (maxfreq - freqcount .eq. 2) then
                  remfreq = 2
               elseif (maxfreq - freqcount .eq. 1) then
                  remfreq = 1
               else
                  close (ich)
                  return
               end if
               do j = 1, remfreq
                  tmp = adjustl(tmp(io2 + 1:))
                  read (tmp, *) freqs(3*i + j)
                  freqcount2 = freqcount2 + 3
               end do
               do j = 1, 6
                  read (ich, '(a)', iostat=io) tmp
                  if (io < 0) then
                     close (ich)
                     return
                  end if
               end do
               ! TODO I hope we can write this nicer some day ..
               if (maxfreq - freqcount .ge. 3) then
                  do j = 1, nat
                     read (ich, '(a)', iostat=io) tmp
                     read (tmp, *) dum1, dum2, &
                     & modes(3*i + 1, 1, j), modes(3*i + 1, 2, j), modes(3*i + 1, 3, j), &
                     & modes(3*i + 2, 1, j), modes(3*i + 2, 2, j), modes(3*i + 2, 3, j), &
                     & modes(3*i + 3, 1, j), modes(3*i + 3, 2, j), modes(3*i + 3, 3, j)
                  end do
                  freqcount = freqcount + 3
               elseif (maxfreq - freqcount .eq. 2) then
                  do j = 1, nat
                     read (ich, '(a)', iostat=io) tmp
                     read (tmp, *) dum1, dum2, &
                     & modes(3*i + 1, 1, j), modes(3*i + 1, 2, j), modes(3*i + 1, 3, j), &
                     & modes(3*i + 2, 1, j), modes(3*i + 2, 2, j), modes(3*i + 2, 3, j)
                  end do
                  freqcount = freqcount + 2
               elseif (maxfreq - freqcount .eq. 1) then
                  do j = 1, nat
                     read (ich, '(a)', iostat=io) tmp
                     read (tmp, *) dum1, dum2, &
                     & modes(3*i + 1, 1, j), modes(3*i + 1, 2, j), modes(3*i + 1, 3, j)
                  end do
                  freqcount = freqcount + 1
               else
                  close (ich)
                  return
               end if
               do j = 1, 2
                  read (ich, '(a)', iostat=io) tmp
                  if (io < 0) then
                     close (ich)
                     return
                  end if
               end do
               read (ich, '(a)', iostat=io) tmp
               if (io < 0) then
                  close (ich)
                  return
               end if
            end do
            exit
         end if
      end do
      close (ich)
      return
   end subroutine rdg98modes
! read all modes from g98 file
   subroutine rdg98allmodes(nat, modes)
      implicit none

      character(len=512) :: tmp, dum, str
      real(wp), allocatable :: modes(:, :, :) ! number of mode, number of atom, x, y or z coordinate
      integer :: io, ich
      integer :: nat
      integer ::  i, j
      integer :: s
      integer :: dum1, dum2

      modes = 0.0d0
      open (newunit=ich, file='g98.out')
      do
         read (ich, '(a)', iostat=io) tmp
         if (io < 0) exit
         io = index(tmp, 'Atom AN', .true.)
         if (io .ne. 0) then

            do j = 1, nat
               read (ich, '(a)', iostat=io) tmp
               read (tmp, *) dum1, dum2, &
               & modes(1, j, 1), modes(1, j, 2), modes(1, j, 3), &
               & modes(2, j, 1), modes(2, j, 2), modes(2, j, 3), &
               & modes(3, j, 1), modes(3, j, 2), modes(3, j, 3)
               ! do i = 1, 3
               !  tmp=adjustl(tmp(io+1:))
               ! write(*,*) "mode is"
               write (*, *) "mode is", dum1, dum2, &
                  & modes(1, j, 1), modes(1, j, 2), modes(1, j, 3), &
                  & modes(2, j, 1), modes(2, j, 2), modes(2, j, 3), &
                  & modes(3, j, 1), modes(3, j, 2), modes(3, j, 3)
               !end do
            end do
            exit
         end if
      end do
      close (ich)
      return
   end subroutine rdg98allmodes

   ! fallback for RRKM calculations if no IRC was found
   subroutine rdlowestfreqg98(freq)
      implicit none

      character(len=512) :: tmp, dum, str

      real(wp) :: freq

      integer :: io, io1, io2, ich

      integer ::  i, j
      integer :: dum1, dum2

      open (newunit=ich, file='g98.out')
      do

         read (ich, '(a)', iostat=io) tmp
         if (io < 0) exit
         io1 = index(tmp, 'Frequencies --', .true.)
         if (io1 .gt. 0) then
            io2 = io1 + 14
            tmp = adjustl(tmp(io2 + 1:))
            read (tmp, *) freq
            exit
         end if
      end do
      close (ich)
      return
   end subroutine rdlowestfreqg98

! preprare interpolated.xyz file as guess for transition state search
   subroutine prepgeodesic(env, jobcall)
      implicit none
      type(runtypedata) :: env
      character(len=1024), intent(out) :: jobcall
      call copy('start.xyz', 'geoin.xyz')
      call appendto('end.xyz', 'geoin.xyz')
      write (jobcall, '(a,i0,a)') 'geodesic_interpolate geoin.xyz --nimages ', env%tsnds, ' >geoout 2>geo2out'
   end subroutine prepgeodesic


!
  
! TODO move into orca module
! Routine to prepare neb calculation
   subroutine prepneb(env, jobcall, hbonddiss, restart, usegeo)
      implicit none
      integer :: nnds, i, io, ich
      integer :: chrg, uhf, mult
      integer :: etemp
      character(len=80) :: levelkeyword
      character(len=1024), intent(out) :: jobcall

      character(len=80) :: prod
      character(len=80) :: xtbstring
      real(wp) :: maxfp, rmsfp ! TODO critical parameter
      integer :: maxiter
      integer    :: maxtime ! time for timeout command for neb search (gets stuck sometimes)
      type(runtypedata) :: env
      logical, intent(in)  :: restart
      logical, intent(in) :: usegeo ! use geodesic as guess
      logical :: ex, re
      logical :: fermi
      logical, intent(in) :: hbonddiss ! only H-dissociation, very loose settings, currently not used

      re = restart
      if (re) then
         write (*, *) "calculation failed, retrying with different settings"
      end if

      if (env%fermi) then
         fermi = .true.
      else
         fermi = .false.
      end if
      if (re) then
         fermi = .true.
      end if

      if (env%eltemp_hybrid .gt. 0 ) then 
         fermi = .true.
      end if

      if (env%geolevel == 'gfn2' .or. env%geolevel == 'gfn1') then
         if (env%eltemp_gga .gt. 0) then
            fermi = .true.
         end if
      end if

      !maxiter = 100 ! complicated rearr. anyway strange !500 default
      maxiter = 500 ! default

      ! hbond dissociations always same course but take converged settings for less randomness
      ! if (hbonddiss) then
      !   maxfp = 0.01_wp
      !    rmsfp = 0.005_wp
      !    maxiter = 15
      ! end if

      call rdshort_int('.CHRG', chrg)
      call rdshort_int(".UHF", uhf)
      !is this right ?
      mult = uhf + 1
      select case (env%geolevel)
      case ('pbeh3c')
         levelkeyword = 'UKS PBEh-3c'
         etemp = 13400
      case ('b973c')
         levelkeyword = 'UKS B97-3c'
         etemp = 5000
      case ('r2scan3c')
         levelkeyword = 'UKS R2SCAN-3c'
         etemp = 5000
      case ('wb97x3c')
         levelkeyword = 'UKS wB97X-3c'
         etemp = 15000
         if (env%eltemp_hybrid .gt. 0) etemp = env%eltemp_hybrid
     !!!!!!!! currently not needed REMOVE???
      case ('gfn2')
         levelkeyword = 'XTB2'
         etemp = 300
         if (env%eltemp_gga .gt. 0) etemp = env%eltemp_gga
         if (re) etemp = 5000 !force convergence
      case ('gfn2spinpol')
         levelkeyword = 'XTB2'
         etemp = 300
         if (env%eltemp_gga .gt. 0) etemp = env%eltemp_gga
         if (re) etemp = 5000
      case ('gfn2_tblite')
         levelkeyword = 'XTB2'
         etemp = 300
         if (env%eltemp_gga .gt. 0) etemp = env%eltemp_gga
         if (re) etemp = 5000         
      case ('gfn1')
         levelkeyword = 'XTB1'
         etemp = 300
         if (env%eltemp_gga .gt. 0) etemp = env%eltemp_gga
         if (re) etemp = 5000 !force convergence
      case ('gxtb')
         levelkeyword = 'XTB'
         fermi = .false.
        
         if (env%eltemp_hybrid .gt. 0) etemp = env%eltemp_hybrid
         if (re) then
            etemp = 15000 ! force convergence
            fermi = .true.
         end if
      case default
         levelkeyword = trim(env%geolevel)
         etemp = 5000
      end select

      if (usegeo) call prepare_neb_restart(env) ! prepare restart_ALLXYZFile for NEB

      nnds = env%tsnds - 2
      ! nnds = env%tsnds !  only if preopt end true selected reduce it by two to make it comparable to geodesic, xtb and gsm
     
      open (newunit=ich, file='orca.inp')
      open (newunit=ich, file='orca.inp')

      write (ich, '(a)') '! NEB   ' ! "LOOSE-NEB" doesnt work ....
      
      if (env%solv) then
         if (env%geolevel == 'gfn2' .or. env%geolevel == 'gfn1' .or. env%geolevel == 'gfn2spinpol' .or. env%geolevel == 'gfn2_tblite') then
         write (ich, *) "! ALPB(water)" 
         end if
     end if

      
      write (ich, *) "! "//trim(levelkeyword)
      if (env%geolevel .ne. 'gfn1' .and. env%geolevel .ne. 'gfn2' .and. env%geolevel .ne. 'gfn2spinpol') then
         write (ich, *) "! LOOSESCF UKS" ! DFT calculations with UKS for correct dissociation and LOOSESCF for faster convergence
      end if
      write (ich, *) "%maxcore 8000" ! TODO make parameter or read in orca sample input file


      !use  tblite for proper uhf scf
      if (env%geolevel == 'gfn2_tblite') then
         xtbstring = 'XTBINPUTSTRING2 "--tblite "'
         write (ich, *) "%xtb"
         write (ich, '(a)') trim(xtbstring)
         write (ich, *) "end"
      end if


      !use  spinpol
      if (env%geolevel == 'gfn2spinpol') then
         xtbstring = 'XTBINPUTSTRING2 "--tblite --spinpol"'
         write (ich, *) "%xtb"
         write (ich, '(a)') trim(xtbstring)
         write (ich, *) "end"
      end if

      if (env%geolevel == 'gfn2' .or. env%geolevel == 'gfn1') then 
         if (fermi)  then 
            write(xtbstring,'(a,i0,a)') 'XTBINPUTSTRING2 "--etemp ',etemp,'"'
            write (ich, *) "%xtb"
            write (ich, '(a)') trim(xtbstring)
            write (ich, *) "end"
         end if
      end if

      if (env%geolevel == 'gxtb') then
         xtbstring = 'XTBINPUTSTRING2 "--driver ''gxtb -c xtbdriver.xyz -symthr 0.0  -b  ~/.basisq ''"'
         if (fermi)  write(xtbstring,'(a,i0,a)') 'XTBINPUTSTRING2 "--driver ''gxtb -c orca.xtbdriver.xyz -tel  -b  ~/.basisq  ',etemp,' -symthr 0.0''"'
         write (ich, *) "%xtb"
         write (ich, '(a)') trim(xtbstring)
         write (ich, *) "end"
         call touch('.GRAD')
      end if
      ! TODO problem of xtbdriver here,dirty first solution for now
      write (ich, *) "%pal"
      if (trim(env%geolevel) == "gxtb") then
         write (ich, '(a)') "nprocs 1"
      else
         write (ich, '(a,i0)') "nprocs ", env%threads
      end if
      write (ich, *) "end"

      if (fermi .or. env%geolevel == "gfn2") then
         write (ich, *) "%scf"
         write (ich, *) "SmearTemp ", etemp
         write (ich, *) "end"
      end if
      write (ich, '(a)') '%NEB NEB_END_XYZFILE "end.xyz" '
      write (ich, '(a,i0)') 'NImages ', nnds
      !write (ich, '(a)') 'Free_End true'       !  can help to "force" TS search to maximum and prevents ending on a maximum !
      !write (ich, '(a)') ' Free_End_Type Perp' ! problematic, can lead to rearrangments of product, deactivate ...
      if (usegeo) then
         inquire (file='nebguess.xyz', exist=ex)
         if (ex) write (ich, '(a)') 'Restart_ALLXYZFile  "nebguess.xyz"' ! in tests geodesic showed to make results little better but it takes more time
      end if
      write (ich, '(a,i0)') 'MaxIter ', maxiter
      if (.not. env%nebnormal) then
         maxfp = 0.01_wp ! 0.001_wpdefault
         rmsfp = 0.005_wp
         ! rmsfp = 0.0025_wp ! 0.0005_wp default
         write (ich, '(a,f9.7)') 'Tol_MaxFP_I ', maxfp! 0.001  ! 0.01
         write (ich, '(a,f9.7)') 'Tol_RMSFP_I ', rmsfp !0.0005 ! Lose 0.005
      end if
      write (ich, '(a)') 'end'
      !   write(ich,'(a)')   '%geom     '
      !  write(ich,'(a)')   'MaxIter 200' ! we do this to prevent too many cycles here but 50 is default ! TODO critical parameter
      !  write(ich,'(a)')   'end'
      write (ich, *) "*xyzfile ", chrg, " ", mult, " start.xyz"
      close (ich)

      ! call orca with full path name
      write (jobcall, '(a)') '$(which orca) orca.inp > orca.out 2> orcaerr.out && touch neb_finished '
   end subroutine prepneb

   !> prepare gsm.orca calculation
   subroutine prepgsm(env, jobcall,hbonddiss, restart, usegeo)
      implicit none
      integer :: nnds, i, io, ich, tsre
      character(len=1024), intent(out) :: jobcall
      character(len=1024) :: copycall
      character(len=1024) :: cefinecall1, cefinecall2, adgcall
      character(len=80) :: geolevel
      integer :: maxtime
      integer :: maxiter
      integer :: chrg, mult, uhf
      real(wp) :: convtol, addnodetol
      integer :: etemp
      type(runtypedata) :: env
      logical, intent(in), optional  :: restart
      logical, intent(in) :: usegeo
      logical, intent(in) :: hbonddiss ! only H-dissociation, very loose settings
      logical :: ex, re
      logical :: fermi
    

      if (present(restart)) then
         re = restart
      else
         re = .false.
      end if

      if (env%fermi) then
         fermi = .true.
      else
         fermi = .false.
      end if
      if (re) then 
         fermi = .true.
      end if 

      maxiter = 20 ! 15
      convtol = 0.001_wp ! default 0.0005
      addnodetol = 0.25_wp ! default 0.1

      !if (hbonddiss) then 
      !   maxiter = 10 ! 10
      !   convtol = 0.01_wp ! default 0.0005
     !    addnodetol = 0.5_wp ! default 0.1
     ! end if
      ! use geodesic as guess ?
      if (usegeo) then
         tsre = 1
         call move('interpolated.xyz', 'restart.xyz0000')
      else
         tsre = 0
      end if

      nnds = env%tsnds
      write (*, *) "nnds is", nnds
      ! TODO FIXME, maybe use template files here instead
      !erstmal mit gsm.orca und tm2orca.py in path
      !inpfileq und ograd wird geschrieben
      ! prepare inpfileq
      open (newunit=ich, file='inpfileq')
      write (ich, '(a)') '------------- QCHEM Scratch Info ------------------------ '
      write (ich, '(a)') '$QCSCRATCH/                               '                !  path for scratch dir. end with "/"'
      write (ich, '(a)') '    GSM_go1q                              '                         ! # name of run
      write (ich, '(a)') ' --------------------------------------------------------- '
      write (ich, '(a)') ' ------------ String Info --------------------------------'
      write (ich, '(a)') ' SM_TYPE                 GSM      '    ! SSM, FSM or GSM
      write (ich, '(a,i0)') ' RESTART               ', tsre            !1  read restart.xyz ! geodesic guess?
      write (ich, '(a,i0)') ' MAX_OPT_ITERS      ',maxiter      !15 should be sufficient, complicated paths anyway strange! 100 ! maximum iterations
      write (ich, '(a)') ' STEP_OPT_ITERS          30       '    ! for FSM/SSM
     ! write (ich, '(a,f9.7)') ' CONV_TOL                       ',convtol  !0.0005       ! perp grad
     ! write (ich, '(a,f9.7)') ' ADD_NODE_TOL                    ',addnodetol !0.25      '   ! 0.1 ! for GSM
      write (ich, '(a)') ' CONV_TOL                       0.0005'       ! perp grad
      write (ich, '(a)') ' ADD_NODE_TOL                    0.1      '   ! 0.1 ! for GSM
      write (ich, '(a)') ' SCALING                                1.0      '    ! for opt steps
      write (ich, '(a)') ' SSM_DQMAX               0.8      '    ! add step size
      write (ich, '(a)') ' GROWTH_DIRECTION        2        '    ! normal/react/prod: 0/1/2 ! most reactions endothermic, TS near end Hammond's postulate timing for rearr worse
      write (ich, '(a)') ' INT_THRESH              2.0      '    ! intermediate detection
      write (ich, '(a)') ' MIN_SPACING             5.0      '    ! node spacing SSM
      write (ich, '(a)') ' BOND_FRAGMENTS          1        '    ! make ICs for fragments
      write (ich, '(a)') ' INITIAL_OPT             5        '   ! opt steps first node
      write (ich, '(a)') ' FINAL_OPT               100     '    ! opt steps last SSM node
      write (ich, '(a)') ' PRODUCT_LIMIT           300.0    '    ! kcal/mol
      write (ich, '(a)') ' TS_FINAL_TYPE           0        '    ! any/delta bond: 0/1
      write (ich, '(a,i0)') ' NNODES        ', nnds           ! including endpoints
      write (ich, '(a)') ' ---------------------------------------------------------'
      close (ich)
      ! prepare ograd
      ! write xtbcall here
      if (env%geolevel == "gfn1" .or. env%geolevel == "gfn2" .or. env%geolevel == "gfn2spinpol" &
      & .or. env%geolevel == "gfn2_tblite" .or. env%geolevel == "pm6" .or. env%geolevel == "gxtb" .or. env%geolevel == "aimnet2") then
         etemp = 300.0_wp
         if (re) etemp = 5000.0_wp
         write (jobcall, '(a)') 'xtb $ofile.xyz -grad --'//trim(env%geolevel)
         if (env%fermi) write (jobcall, '(a,i4)') trim(jobcall)//' --etemp ', etemp
         if (env%dxtb) write (jobcall, '(a)') trim(jobcall)//' --vparam dxtb_param.txt'
         if (env%geolevel == 'gfn2spinpol') write (jobcall, '(a)') trim(jobcall)//' --tblite --spinpol'
         if (env%geolevel == 'gfn2_tblite') write (jobcall, '(a)') trim(jobcall)//' --tblite'
   
         write (jobcall, '(a)') trim(jobcall)//' > $ofile.xtbout'
         if (env%geolevel == 'gxtb') then
            write (jobcall, '(a)') ' gxtb -c $ofile.xyz -grad   -b  ~/.basisq > $ofile.xtbout'
            if (env%fermi) write (jobcall, '(a)') ' gxtb -c $ofile.xyz -grad -tel 15000   -b  ~/.basisq > $ofile.xtbout'
         end if
         if (env%geolevel == 'aimnet2') then
            call remove('coord') ! necessary the way the xtb driver works at the momen TODO fix this
            write (jobcall, '(a)') 'xtb $ofile.xyz -grad --driver aimnet2grad  > $ofile.xtbout'
         end if
        
         open (newunit=ich, file='ograd')
         write (ich, '(a)') '#!/bin/bash                                                                               '
         write (ich, '(a)') '# The path to ORCA should be added to .bashrc or exported in command line                 '
         write (ich, '(a)') 'export PATH=/usr/software/orca_web:$PATH                                                  '
         write (ich, '(a)') 'export PATH=/usr/software/orca_web/openmpi/bin:$PATH                                      '
         write (ich, '(a)') 'export LD_LIBRARY_PATH=/usr/bin:/usr/software/orca_web/openmpi/lib:$LD_LIBRARY_PATH'
         write (ich, '(a)') '                                                                                          '
         write (ich, '(a)') 'if [ -z $2 ]                                                                              '
         write (ich, '(a)') 'then                                                                                      '
         write (ich, '(a)') '  echo " need two arguments "                                                             '
         write (ich, '(a)') '  exit                                                                                    '
         write (ich, '(a)') 'fi                                                                                        '
         write (ich, '(a)') '                                                                                          '
         write (ich, '(a)') '#echo " in ograd: $1 $2 "                                                                 '
         write (ich, '(a)') '                                                                                          '
         write (ich, '(a)') 'ofile=orcain$1.in                                                                         '
         write (ich, '(a)') 'ofileout=orcain$1.out                                                                     '
         write (ich, '(a)') 'molfile=structure$1                                                                       '
         write (ich, '(a)') 'ncpu=$2                                                                                   '
         write (ich, '(a)') 'basename="${ofile%.*}"                                                                    '
         write (ich, '(a)') '#echo " ofile: $ofile ofileout: $ofileout molfile: $molfile ncpu: $ncpu"                  '
         write (ich, '(a)') '                                                                                          '
         write (ich, '(a)') '########## XTB/TM settings: #################                                             '
         write (ich, '(a)') 'cd scratch                                                                                '
         write (ich, '(a)') 'wc -l < $molfile > $ofile.xyz                                                             '
         write (ich, '(a)') 'echo "Dummy for XTB/TM calculation" >> $ofile.xyz                                         '
         write (ich, '(a)') 'cat $molfile >> $ofile.xyz                                                                '
         write (ich, '(a)') 'mkdir -p $basename'
         write (ich, '(a)') 'cd $basename'
         write (ich, '(a)') 'cp ../$ofile.xyz .'
         write (ich, '(a)') 'cp ../.CHRG .'! important !!!!
         write (ich, '(a)') 'cp ../.UHF .' ! important !!!!
         write (ich, '(a)') trim(jobcall)
         write (ich, '(a)') 'tm2orca.py $basename '
         write (ich, '(a)') 'mv $basename.out ../$basename.out'
         write (ich, '(a)') 'mv $basename.engrad ../$basename.engrad'
         write (ich, '(a)') 'cd ..'
         write (ich, '(a)') 'cd ..'
         close (ich)
   elseif (env%geolevel == "pbe" .or. env%geolevel == "r2scan3c" .or. env%geolevel == "b973c" .or. env%geolevel == "pbeh3c" .or. env%geolevel == "wb97x3c" .or. env%geolevel == "pbe0") then ! todo make check for dft here
      select case (env%geolevel)
      case ('pbeh3c')
         geolevel = 'PBEh-3c'
         etemp = 13400
      case ('r2scan3c')
         geolevel = 'R2SCAN-3c'
         etemp = 5000
      case ('b973c')
         geolevel = 'B97-3c'
         etemp = 5000
      case ('pbe')
         geolevel = 'PBE SV(P)'
         etemp = 5000      
      case ('wb97x3c')
         geolevel = 'wB97X-3c'
         etemp = 15000 ! TODO tune this value
      case ('pbe0')
         geolevel = 'PBE0 def2-TZVP'
		 !if (env%chrg .lt. 0) levelkeyword = 'PBE0 ma-def2-TZVP'
       if (env%chrg .lt. 0) geolevel = 'PBE0 def2-TZVPD SLOWCONV'
         etemp = 10000
      case default
         geolevel = trim(env%geolevel)
         etemp = 5000
      end select
       
            call rdshort_int('.CHRG', chrg)
            call rdshort_int(".UHF", uhf)
      !is this right ?
            mult = uhf + 1
         open (newunit=ich, file='ograd')
            write (ich, '(a)') '#!/bin/bash  '     
            write (ich, '(a)') 'ofile=scratch/orcain$1.in'
            write (ich, '(a)') 'ofileout=scratch/orcain$1.out'
            write (ich, '(a)') 'molfile=scratch/structure$1'   ! for ORCA 5.0 
            write (ich, '(a)') 'ncpu=$2'  
            write (ich, '(a)') 'basename="${ofile%.*}"'                        
            write (ich, '(a)') 'echo "! UKS '//trim(geolevel)//'" > $ofile' ! just use development version of orca for this'
            ! leads to abortion due to bug in ORCA write (ich, '(a)') 'echo "! RIJONX def2/J" >> $ofile'
            write (ich, '(a)') 'echo "! LOOSESCF" >> $ofile'
            write (ich, '(a)') 'echo "! ENGRAD" >> $ofile'
            write (ich, '(a)') 'echo "! nomoprint" >> $ofile'
            write (ich, '(a)') 'echo "%scf" >> $ofile'
            write (ich, '(a)') 'echo "maxiter 350" >> $ofile'
            if (fermi) write (ich, '(a,i0,a)') 'echo " SmearTemp ',etemp,'" >> $ofile'
            write (ich, '(a)') 'echo "end" >> $ofile'
            write (ich, '(a)') 'echo "%pal" >> $ofile'
            write (ich, '(a,i0,a)') 'echo " nprocs ',env%threads,'" >> $ofile'
            write (ich, '(a)') 'echo "end" >> $ofile'
            write (ich, '(a,i0,1x,i0,1x,a)') 'echo "*xyz ', chrg, mult,'">> $ofile'
            write (ich, '(a)') 'cat $molfile  >> $ofile'
            write (ich, '(a)') 'echo "*" >> $ofile'
            write (ich, '(a)') 'orcabin=$(which orca)'
            write (ich, '(a)') '$orcabin $ofile > $ofileout'
            write (ich, '(a)') 'tm2orca.py $basename '
            !write (ich, '(a,i0,a)') 'srun --nodes 1 --ntasks ',env%threads, ' --cpus-per-task 1 $orcabin $ofile > $ofileout' 

            close(ich)
   end if
      call execute_command_line('chmod u+x ograd')
      io = makedir('scratch')
      call execute_command_line('sed -i "2s/.*/ /" start.xyz')
      call copy('start.xyz', 'initial0000.xyz')
      call copy('.CHRG', 'scratch/.CHRG')
      call copy('.UHF', 'scratch/.UHF')
      !somehow its important to delete this line to prevent bugs with gsm.orca
      ! TODO FIXME necessary??
      call execute_command_line('sed -i "2s/.*/ /" end.xyz')
      call appendto('end.xyz', 'initial0000.xyz')

      call copy('initial0000.xyz', 'scratch/initial0000.xyz')
      if (env%dxtb) then
         write (copycall, '(a)') "cp "//trim(env%startdir)//"/dxtb_param.txt scratch/dxtb_param.txt"
         call execute_command_line(trim(copycall))
      end if
     

      !somehow its important to delete this line to prevent bugs with gsm.orca
      ! TODO FIXME necessary??
      !  call execute_command_line('sed -i "2s/.*/ /" scratch/initial0000.xyz')
      write (*, *) "Starting GSM Run with", nnds, "nodes"
      !how to deal with GSM which does not converge?
      !timeout command, make dependent on number of atoms?
      maxtime = 3600000 ! 1000 hours
      if (trim(env%geolevel) == 'gfn2' .or. trim(env%geolevel) == 'gfn1') maxtime = 3600
      if (trim(env%geolevel) == 'gxtb') maxtime = 10000 ! TODO make dependent on number of atoms
       write(jobcall,'(a,i0,a)') 'timeout ',maxtime,' gsm > gsm.out 2> gsmerror.out'
      write (jobcall, '(a)') trim(jobcall)//' && touch gsm_finished '
      if (env%printlevel .le. 1) write (jobcall, '(a)') trim(jobcall)//' && rm -r scratch '

      ! && cp scratch/tsq0000.xyz ts.xyz'

   end subroutine prepgsm

!> check if TS search worked technically
   subroutine checkpath(env, failed)
      implicit none
      type(runtypedata) :: env
      logical :: ex
      logical, intent(out) :: failed
      integer :: nlines
      character(len=80) :: fname

      if (env%tsfinder == "gsm") then
         fname = 'stringfile.xyz0000'
      elseif (env%tsfinder == "xtb") then
         fname = 'xtbpath.xyz'
      elseif (env%tsfinder == "neb") then
         fname = 'orca_MEP_trj.xyz'
      elseif (env%tsfinder == "crest_tsgen") then
         fname = 'crestopt.xyz'      
      end if

      inquire (file=trim(fname), exist=ex, size=nlines)
      if (ex .and. nlines .gt. 0) then ! gsm can fail and produces and empty stringfile
         failed = .false.
      else
         failed = .true.
      end if

   end subroutine checkpath
!> check if TS search gave reasonable path and recognize multiple reaction steps in one path
!> TOOD build in to spllit multiple reaction steps in different paths
!> for now they can be sorted out with "sortoutcascade" keyword
   subroutine pickts(env, found)
      implicit none
      integer :: ich, io, i, j
      integer :: natoms, nnds
      integer :: nmax ! number of maxima in path
      integer :: lastpoint
      integer :: nlines
      real(wp), allocatable :: e_rels(:) ! relative energies
      real(wp) :: diff1, diff2, barrier, barrier_loc, emin, emin_loc, emax, tsthr, diffthr, diff4
      character(len=80), allocatable :: xyz(:)
      type(runtypedata) :: env
      character(len=512) :: tmp
      character(len=80) :: sdum
      character(len=80) :: fname
      logical :: ex, ex2
      logical :: found

      found = .false. ! successful?

      nnds = env%tsnds
      allocate (e_rels(nnds))
      if (env%tsfinder == "gsm") then
         fname = 'stringfile.xyz0000'
      elseif (env%tsfinder == "xtb") then
         fname = 'xtbpath.xyz'
      elseif (env%tsfinder == "neb") then
         fname = 'orca_MEP_trj.xyz'
      elseif (env%tsfinder == "crest_tsgen") then
         inquire(file='crestopt.xyz', exist=ex)
         if (ex) then 
            call copy('crestopt.xyz', 'ts.xyz')
            found = .true.
            return
         else
            write (*, *) "WARNING: no crestopt.xyz file found, search was not succesfull"
            found = .false.
            call printpwd() 
         end if   
      else
         write (*, *) "TS PICKER NOT YET IMPLEMENTED FOR OTHER SEARCHERS"
         write (*, *) "Just take end as TS"
         call copy('end.xyz', 'ts.xyz')
         return
      end if
      !even if search failed, stringfile usually there
      ! read xyz from stringfile
      e_rels = 0.0_wp

      if (env%tsfinder == "gsm" .and. env%tsgeodesic) then
         ! lastpoint = nnds - 1 !end point artificially low with geodesic sometimes -> pseudo maximum at nnds -1 po
         lastpoint = nnds
      else
         lastpoint = nnds
      end if

      tsthr = 0.0_wp
      inquire (file=trim(fname), exist=ex, size=nlines)
      if (ex .and. nlines .gt. 0) then ! gsm can fail and produces and empty stringfile
         open (newunit=ich, file=trim(fname))
         ! if file is corrupted we dont want to crash the whole programm
         read (ich, '(a)', iostat=io) tmp
         if (io < 0) GOTO 100
         backspace (ich)
         read (ich, *) natoms
         allocate (xyz(natoms))
         if (env%tsfinder == "gsm") read (ich, *) e_rels(1)
         if (env%tsfinder == "xtb") read (ich, *) sdum, e_rels(1), sdum, sdum, sdum
         if (env%tsfinder == "neb") read (ich, *) sdum, sdum, sdum, sdum, sdum, e_rels(1)
         emax = e_rels(1)
         do j = 1, natoms
            ! if file is corrupted we dont want to crash the whole programm
            read (ich, '(a)', iostat=io) tmp
            if (io < 0) GOTO 100
            backspace (ich)
            read (ich, '(a)') xyz(j)
         end do
         do i = 2, lastpoint
            ! if file is corrupted we dont want to crash the whole programm
            read (ich, '(a)', iostat=io) tmp
            if (io < 0) GOTO 100
            backspace (ich)
            read (ich, *) natoms
            if (env%tsfinder == "gsm") read (ich, *) e_rels(i)
            if (env%tsfinder == "xtb") read (ich, *) sdum, e_rels(i), sdum, sdum, sdum
            if (env%tsfinder == "neb") read (ich, *) sdum, sdum, sdum, sdum, sdum, e_rels(i)
            if (env%tsfinder == "neb") e_rels(i) = (e_rels(i) - e_rels(1))*autokcal ! transform to relative kcal/mol
            !DEL    write(*,*) "erel is", e_rels(i)
            if (e_rels(i) - emax .gt. tsthr) then
               emax = e_rels(i)
               do j = 1, natoms
                  read (ich, '(a)') xyz(j)
               end do
            else
               do j = 1, natoms
                  read (ich, '(a)') tmp
               end do
            end if
         end do
100      close (ich)

         ! find maxima of reaction path
         emin = e_rels(1)
         emin_loc = e_rels(1)
         nmax = 0

         ! but  threshold values is needed take 1 kcal?
         tsthr = 1.0_wp ! here different threshold than for taking highest point

         diffthr = 0.0_wp ! can not set too high value, if sampling is very dense -> no maximum is detected

         if (env%tsfinder == "gsm" .and. env%tsgeodesic) then
            lastpoint = nnds - 1
            ! lastpoint = nnds - 2  !TODO investigate this in more detail
            !end point artificially low with geodesic sometimes -> pseudo maximum at nnds -1 po
         else
            lastpoint = nnds - 1
         end if
         do i = 2, lastpoint !
            if (e_rels(i) .lt. emin) emin = e_rels(i)
            if (e_rels(i) .lt. emin_loc) emin_loc = e_rels(i)
            diff1 = e_rels(i) - e_rels(i - 1)
            diff2 = e_rels(i) - e_rels(i + 1)
            barrier = e_rels(i) - emin ! difference to prior minimum
            barrier_loc = e_rels(i) - emin_loc    ! add also condition that second maximum has to be significant higher than prior minimum
            if ((diff1 .gt. diffthr) .and. (diff2 .gt. diffthr) .and. (barrier .gt. tsthr) .and. (barrier_loc .gt. tsthr)) then
               nmax = nmax + 1
               if (env%printlevel .eq. 3) write (*, *) "found maximum for node", i
               emin_loc = e_rels(i) ! set pre minimum to highest point, so that it newly detected for the next barrier
            end if
         end do
         if (nmax .gt. 1) then
            if (env%sortoutcascade) then
               write (*, *) "WARNING: more than one maximum in path remove path because its a reaction cascade"
               found = .false.
            else
               write (*, *) "WARNING: more than one maximum in path, but we take highest one"
               found = .true. ! succesfull ts
            end if
            call touch("cascade")
         elseif (nmax .eq. 1) then
            if (env%printlevel .eq. 3) write (*, *) "only one maximum, search was succesfull"
            found = .true. ! succesfull ts
         elseif (nmax .eq. 0) then
            !DEL write (*, *) "no TS have to take endpoint"
            found = .true. ! succesfull ts
            call touch("endeqts") ! ToDO fixme, what to do with these???
            ! write last point of path to ts.xy xyz(j) contains then coordinates of last point
         end if
         ! write ts.xyz
         open (newunit=ich, file='ts.xyz')
         write (ich, *) natoms
         write (ich, *) "EMAX is: ", emax
         if (emax .gt. 200.0_wp) then
            write (*, *) "WARNING: TS energy unreasonably high, the search probably failed" ! DO not restart just sort out ..
            found = .false.
            call printpwd()
            return
         end if
         do j = 1, natoms
            write (ich, '(a)') xyz(j)
         end do
         close (ich)
      else
         write (*, *) "WARNING: no stringfile found, search was not succesfull"
         found = .false.
         call printpwd()
      end if
      if (allocated(xyz)) deallocate (xyz)
   end subroutine pickts

!Check IRC
! IRC is checked by comparing the RMSD of the TS with the start and end point
! As this is for complicated reaction paths not always clear, a scan is performed
! along 10 points forwards and backwards along each mode
! and the RMSD of the TS is compared for each point to start and end structure
! if at least half of the points exhibit changes in RMSD ("score" over > 4*10) that suit to
! the expected behavior of a TS, the IRC is considered as valid
   subroutine findirc(env, nmode, ircmode)
      use mctc_env, only: error_type, get_argument!, fatal_error
      use mctc_io, only: structure_type, read_structure, write_structure, &
     & filetype, get_filetype, to_symbol, to_number
      use structools, only: wrxyz
      implicit none
      logical :: ldum
      real(wp), intent(out) :: ircmode
      integer :: nat
      real(wp), allocatable :: freqs(:)
      real(wp), allocatable :: modes(:, :, :)
      real(wp), allocatable :: start(:, :), end(:, :), ts(:, :) ! geometry of start end, and ts
      real(wp), allocatable :: ts_for(:, :), ts_back(:, :) ! geometry of ts for- and backwards propagated
      integer, allocatable :: iat(:)
      integer :: nimags
      type(structure_type) :: mol
      type(error_type), allocatable :: error
      real(wp) :: rmsd_ts_start, rmsd_ts_end
      real(wp) :: rmsd_ts_for_start, rmsd_ts_for_end, rmsd_ts_back_start, rmsd_ts_back_end
      real(wp), allocatable :: diff_for_start(:), diff_back_start(:), diff_for_end(:), diff_back_end(:)
      real(wp) :: rmsd_thr ! threshold for significant change in RMSD, tune this parameter
      real(wp), allocatable :: gradient(:, :)
      real(wp) :: trafo(3, 3)
      integer :: i, j, ndispl
      integer, intent(out) :: nmode ! number of mode
      integer :: score ! a "score" to evaluate if we have a good IRC
      type(runtypedata) :: env
      character(len=1024) :: jobcall
      logical :: ex

      call rdshort_int('ts.xyz', nat)
      allocate (modes(9, 3, nat))
      allocate (freqs(9))
      call move('orca.g98.out', 'g98.out')
      inquire (file='g98.out', exist=ex)
      if (.not. ex) then
         write (*, *) "ERROR: Could not find orca.g98.out, Hessian calculation failed"
         call printpwd
         nmode = 0
         ircmode = 0.0_wp
         return
      end if
      call rdg98modes(nat, freqs, modes, nimags)
      if (nimags .eq. 0) then
         write (*, *) "NO IMAGINARY FREQUENCY FOUND"
         call printpwd
         nmode = 0
         return
      end if
      ! somehow the comment line makes problems here....
      call read_structure(mol, '../../start.xyz', error, filetype%xyz)
      if (allocated(error)) then
      if (error%stat .eq. 1) then ! error reading structure
         write (*, *) "ERROR: Could not read start.xyz"
         nmode = 0
         return
      end if
      end if
      ! have to convert it into angstroem
      start = mol%xyz*bohr
      call read_structure(mol, '../../end.xyz', error, filetype%xyz)
      if (allocated(error)) then
      if (error%stat .eq. 1) then ! error reading structure
         write (*, *) "ERROR: Could not read start.xyz"
         nmode = 0
         return
      end if
      end if
      end = mol%xyz*bohr
      call read_structure(mol, 'ts.xyz', error, filetype%xyz) ! todo we have to delet for this the content of the second line in the xyz file???
      if (allocated(error)) then
      if (error%stat .eq. 1) then ! error reading structure
         write (*, *) "ERROR: Could not read start.xyz"
         nmode = 0
         return
      end if
      end if
      ts = mol%xyz*bohr
      if (size(ts) .ne. size(start) .or. size(ts) .ne. size(end)) then
         write (*, *) "ERROR: TS, start and end have different number of atoms"
         nmode = 0
         return
      end if

      write (*, *) "Identifying IRC mode"

      call calcrmsd(nat, ts, start, rmsd_ts_start)

      if (rmsd_ts_start .lt. 1e-08_wp) then
         write (*, *) "WARNING: TS is identical to start, something went wrong aborting IRC determination"
         nmode = 0
         return
      end if

      call calcrmsd(nat, ts, end, rmsd_ts_end)
      if (rmsd_ts_end .lt. 1e-08_wp) then
         write (*, *) "WARNING: TS is identical to start, something went wrong aborting IRC determination"
         nmode = 0
         return
      end if

      ! distort along mode forward in steps of 0.1
      ndispl = 10
      allocate (diff_for_start(ndispl), diff_for_end(ndispl), diff_back_start(ndispl), diff_back_end(ndispl))
      do j = 1, nimags
         if (freqs(j) .ge. -1.0_wp) exit ! only check negative frequencies
         do i = 1, ndispl
            ts_for = ts + modes(j, :, :)*0.1_wp*i
            call calcrmsd(nat, ts_for, start, rmsd_ts_for_start)
            diff_for_start(i) = rmsd_ts_for_start - rmsd_ts_start
            call calcrmsd(nat, ts_for, end, rmsd_ts_for_end)
            diff_for_end(i) = rmsd_ts_for_end - rmsd_ts_end
            ts_back = ts - modes(j, :, :)*0.1_wp*i
            call calcrmsd(nat, ts_back, start, rmsd_ts_back_start)
            diff_back_start(i) = rmsd_ts_back_start - rmsd_ts_start
            call calcrmsd(nat, ts_back, end, rmsd_ts_back_end)
            diff_back_end(i) = rmsd_ts_back_end - rmsd_ts_end
         end do

         ! JUST PRINTING
         if (env%ircrun) then
            do i = 1, ndispl
               write (*, *) "diff_for_start is ", diff_for_start(i)
            end do
            do i = 1, ndispl
               write (*, *) "diff_for_end is ", diff_for_end(i)
            end do
            do i = 1, ndispl
               write (*, *) "diff_back_start is ", diff_back_start(i)
            end do
            do i = 1, ndispl
               write (*, *) "diff_back_end is ", diff_back_end(i)
            end do
         end if

         ! as we dont now if mode forward is actually forward we need to
         ! discriminate here between 16 possible cases as each criterion can fail sometimes
         ! e.g. if rmsd to end is poor because end is very different due to relaxation at the end
         score = 0
         rmsd_thr = 0.005_wp
         do i = 1, ndispl
            ! forwards ! end can be strange (due to relaxation at the end)
            !but to start should always be significant, so taking this as starting should be fine
            if (diff_for_start(i) .gt. rmsd_thr) then
               if (diff_back_start(i) .lt. -rmsd_thr) then
                  if (diff_for_end(i) .lt. -rmsd_thr) then
                     if (diff_back_end(i) .gt. rmsd_thr) then
                        score = score + 4
                     else
                        score = score + 3
                     end if
                  elseif (diff_for_end(i) .gt. rmsd_thr) then
                     if (diff_back_end(i) .gt. rmsd_thr) then
                        score = score + 3
                     else
                        score = score + 2
                     end if
                  end if
               elseif (diff_back_start(i) .gt. rmsd_thr) then
                  if (diff_for_end(i) .lt. -rmsd_thr) then
                     if (diff_back_end(i) .gt. rmsd_thr) then
                        score = score + 3
                     else
                        score = score + 2
                     end if
                  elseif (diff_for_end(i) .gt. rmsd_thr) then
                     if (diff_back_end(i) .gt. rmsd_thr) then
                        score = score + 1
                     else
                        score = score + 0
                     end if
                  end if
               end if
            elseif (diff_for_start(i) .lt. -rmsd_thr) then
               if (diff_back_start(i) .gt. rmsd_thr) then
                  if (diff_for_end(i) .gt. rmsd_thr) then
                     if (diff_back_end(i) .lt. -rmsd_thr) then
                        score = score + 4
                     else
                        score = score + 3
                     end if
                  elseif (diff_for_end(i) .lt. -rmsd_thr) then
                     if (diff_back_end(i) .lt. -rmsd_thr) then
                        score = score + 3
                     else
                        score = score + 2
                     end if
                  end if
               elseif (diff_back_start(i) .lt. -rmsd_thr) then
                  if (diff_for_end(i) .gt. rmsd_thr) then
                     if (diff_back_end(i) .lt. -rmsd_thr) then
                        score = score + 3
                     else
                        score = score + 2
                     end if
                  elseif (diff_for_end(i) .lt. -rmsd_thr) then
                     if (diff_back_end(i) .lt. -rmsd_thr) then
                        score = score + 1
                     else
                        score = score + 0
                     end if
                  end if
               end if
            end if
            if (env%ircrun) write (*, *) "SCORE IS", i, score
         end do

         ! score ranges from 0 to 40, 40 is best
         !write(*,*) "SCORE IS", score
         ! TODO CHECKME critical parameter ! score of 20 /one half of max means, that at least the rmsd to the start
         ! gets  larger (or smaller) going forward and smaller (or larger) going backwards along the mode
         if (score .ge. ndispl*4/2) then
            ircmode = freqs(j)
            write (*, *) "IRC MODE found! mode ", j, " with ", ircmode
            nmode = j
            return
         else ! continue search
            write (*, *) "have to continue search"
         end if
      end do
      write (*, *) "no suitable IRC mode found within first 9 frequencies for"
      call printpwd
      nmode = 0

      ! DEPRECATED CODE, instead of score we used thresholds for RMSD
      ! to evaluate how good the IRC mode is
      ! to be IRC, adding of the mode to TS geometry should reduce the RMSD
      ! to product and reactant
      ! change of sign for both with respect to forward and backward
      ! and between both
      ! and both changes between forward and backward for start and end
      ! have to be significant, e.g. 0.05A°
      ! if (difffor_start*diffback_start .lt. 0 .and. difffor_end*diffback_end .lt. 0 &
      !& .and. difffor_start*difffor_end .lt. 0 &
      !& .and. abs(difffor_start - diffback_start) .gt. 0.01 &
      !& .and. abs(difffor_end - diffback_end) .gt. 0.01) then
      !   write (*, *) "ircmode found its the first mode: ", freqs(1)
      !   ircmode = freqs(1)

      !END IF
      ! check next mode
      !       else
!
      !          ! distort along mode forward scaled by 0.5, empirically determined
      !          !  write(*,*) "MODE IS", modes(1,:,:)
      !          ts = ts + modes(1, :, :)*0.5_wp
      !  write(*,*) "TSfor XYZ is:"
      !          !  write(*,*) ts
      !          call calcrmsd(nat, ts, start, rmsdfor_start)
      !          write (*, *) "rmsdfor_start is ", rmsdfor_start
      !          call calcrmsd(nat, ts, end, rmsdfor_end)
      !          write (*, *) "rmsdfor_end is ", rmsdfor_end
      !          ! distort along mode backwards scaled by 0.5, empirically determined
      !          ts = ts - modes(1, :, :)*0.5_wp*2
      !          !   write(*,*) "TSback XYZ is:"
      !          !  write(*,*) ts
      !          call calcrmsd(nat, ts, start, rmsdback_start)
      !          write (*, *) "rmsdback_start is ", rmsdback_start
      !          call calcrmsd(nat, ts, end, rmsdback_end)
      !          write (*, *) "rmsback_end is ", rmsdback_end
!
      !          difffor_start = rmsdfor_start - rmsd_start
      !          diffback_start = rmsdback_start - rmsd_start
      !          difffor_end = rmsdfor_end - rmsd_end
      !          diffback_end = rmsdback_end - rmsd_end
!
      !          if (difffor_start*diffback_start .lt. 0 .and. difffor_end*diffback_end .lt. 0 &
      !          & .and. difffor_start*difffor_end .lt. 0 &
      !          & .and. abs(difffor_start - diffback_start) .gt. 0.05 &
      !          & .and. abs(difffor_end - diffback_end) .gt. 0.05) then
      !             write (*, *) "ircmode found its the second mode: ", freqs(2)
      !             ircmode = freqs(2)
      !   !          ! check next mode
      !       else
      !          ! distort along mode forward scaled by 0.5, empirically determined
      !          !  write(*,*) "MODE IS", modes(1,:,:)
      !          ts = ts + modes(1, :, :)*0.5_wp
      !          !  write(*,*) "TSfor XYZ is:"
      !          !  write(*,*) ts
      !          write (*, *) "compute RMSD"
      !          call calcrmsd(nat, ts, start, rmsdfor_start)
      !          write (*, *) "rmsdfor_start is ", rmsdfor_start
      !          call calcrmsd(nat, ts, end, rmsdfor_end)
      !          !  write(*,*) "rmsdfor_end is ", rmsdfor_end
      !          ! distort along mode backwards scaled by 0.5, empirically determined
      !          ts = ts - modes(1, :, :)*0.5_wp*2
      !          !   write(*,*) "TSback XYZ is:"
      !          !  write(*,*) ts
      !          call calcrmsd(nat, ts, start, rmsdback_start)
      !          ! write(*,*) "rmsdback_start is ", rmsdback_start
      !          call calcrmsd(nat, ts, end, rmsdback_end)
      !          !  write(*,*) "rmsback_end is ", rmsdback_end
!
      !          difffor_start = rmsdfor_start - rmsd_start
      !          diffback_start = rmsdback_start - rmsd_start
      !          difffor_end = rmsdfor_end - rmsd_end
      !          diffback_end = rmsdback_end - rmsd_end
!
      !          if (difffor_start*diffback_start .lt. 0 .and. difffor_end*diffback_end .lt. 0 &
      !          & .and. difffor_start*difffor_end .lt. 0 &
      !          & .and. abs(difffor_start - diffback_start) .gt. 0.05 &
      !          & .and. abs(difffor_end - diffback_end) .gt. 0.05) then
      !             write (*, *) "ircmode found its the third mode: ", freqs(3)
      !             ircmode = freqs(3)
      !             ! check next mode
      !          else
      !             write (*, *) "Mode not yet found!"
      !             write (*, *) "Have to check more but take lowest for now"
      !             ircmode = freqs(1)
      !          end if
      !       end if
      !    end if
      ! end if
   end subroutine findirc
! parse orca trajectory if optimization did not converge
! write ts.xyz
   ! just take first point of trajectory, as optimization probably
   ! run away from TS
   subroutine parseorcatrajectory()
      implicit none
      integer :: i, j, io, ich, ichout
      integer :: nat
      integer :: nmols
      integer :: nlines
      character(512) :: line
      character(80) :: fname, fout

      fname = 'orca_trj.xyz'
      fout = 'ts.xyz'
      call rdshort_int(trim(fname), nat)
      open (unit=ich, file=trim(fname), status='old', action='read', iostat=io)
      open (newunit=ichout, file=trim(fout))
      nlines = 0

      do i = 1, nat + 2
         read (ich, '(a)', iostat=io) line
         if (io /= 0) then
            exit
            call printpwd()
            !  error stop "Error reading orca_trj.xyz"
         end if
         write (ichout, '(a)') line
      end do

      close (ich)
      close (ichout)

      ! DEPRECATED took last structure of trajectory here
      ! do
      !    read (ich, *, iostat=io) line
      !    if (io /= 0) exit  !
      !    nlines = nlines + 1
      ! end do
      ! rewind (ich)
!
      ! nmols = nlines/(nat + 2)
      ! do i = 1, nmols - 1
      !    do j = 1, nat + 2
      !       read (ich, '(a)') line
      !    end do
      ! end do
      ! do j = 1, nat + 2
      !    read (ich, '(a)') line
      !    write (ichout, '(a)') line
      ! end do

   end subroutine parseorcatrajectory

!WARNING can lead to full stop of program with no convergence in svdcmp
   subroutine calcrmsd(nat, xyz1, xyz2, root_msd)
      integer, intent(in) :: nat
      real(wp), intent(in) ::  xyz1(:, :), xyz2(:, :)
      real(wp), allocatable ::  nxyz1(:, :), nxyz2(:, :)
      real(wp), allocatable :: cg(:, :) !centers of geometry
      real(wp) :: diff_cg(3)
      real(wp) :: normmass
      !real(wp) :: rmsd_check(3, nat, 2, 50)
      real(wp), intent(out) :: root_msd
      !real(wp) :: gradient(3, nat)
      !real(wp) :: trafo(3, 3)
      integer :: i, j

 !!! > this is the RMSD part - which I stole from Dr. Jay-Dog
      ! >> saved here for easy re-implementation

      allocate (cg(3, 2))
      allocate (nxyz1(3, nat), nxyz2(3, nat))
      !  !> compute the center-of-geometry of both structures

      normmass = 0
      cg(:, :) = 0

      do j = 1, nat
         !    !>> get the current fragment structure
         !  rmsd_check(:,j,i,cnt) = xyzf(:,j,i)
         !>> get the current center-of-geometry
         cg(:, 1) = cg(:, 1) + 1*xyz1(:, j)
         !
         !cg(:,i,cnt) = cg(:,i,cnt) + 1 * rmsd_check(:,j,i,cnt)
         normmass = normmass + 1
      end do

      cg(:, 1) = cg(:, 1)/normmass

      normmass = 0
      do j = 1, nat
         !    !>> get the current fragment structure
         !  rmsd_check(:,j,i,cnt) = xyzf(:,j,i)
         !>> get the current center-of-geometry
         cg(:, 2) = cg(:, 2) + 1*xyz2(:, j)
         !
         !cg(:,i,cnt) = cg(:,i,cnt) + 1 * rmsd_check(:,j,i,cnt)
         normmass = normmass + 1
      end do

      cg(:, 2) = cg(:, 2)/normmass

      !  !> calculate the difference betwen the two c-of-g
      diff_cg(:) = cg(:, 1) - cg(:, 2)

      !  !> shift the c-of-g by the difference of c-of-g
      do j = 1, nat
         nxyz1(:, j) = xyz1(:, j) - diff_cg(:)
      end do

      do j = 1, nat
         nxyz2(:, j) = xyz2(:, j) - diff_cg(:)
      end do
      !  !> transform into right array size for rmsd routine
      !  do j = 1, natf(i)
      !    !>> the right compare-structure is taken
      !    nxyz1(:,j) = rmsd_check(:,j,i,1)
      !    nxyz2(:,j) = rmsd_check(:,j,i,cnt)
      !  enddo

      !  !> calculate root-mean-sqare-deviation of the two structures
      !write(*,*) "xyz1 is", nxyz1
      !write(*,*) "xyz2 is", nxyz2
      !call get_rmsd(nxyz1, nxyz2, root_msd)
      call get_rmsd_for_coord(nxyz1,nxyz2,root_msd)

      deallocate (nxyz1, nxyz2, cg)

   end subroutine calcrmsd

   ! currently not used
   subroutine prepxtbpath(env, jobcall, nnodes)
      implicit none
      integer, intent(in), optional ::  nnodes
      integer :: nnds, i, io, ich
      character(len=1024) :: dir, thisdir
      character(len=1024), intent(out) :: jobcall
      character(len=80) :: reac, prod
      type(runtypedata) :: env
      logical :: ex
      integer :: etemp

      etemp = 300.0_wp

      call copy(trim(reac), 'start.xyz')
      inquire (file='pair.xyz', exist=ex)
      if (ex) prod = 'pair.xyz'
      inquire (file='isomer.xyz', exist=ex)
      if (ex) prod = 'isomer.xyz'
      call copy(trim(prod), "end.xyz")

      if (present(nnodes)) then
         nnds = nnodes
      else
         nnds = 5
      end if

      ! write path.inp file
      open (newunit=ich, file='path.inp')
      write (ich, '(a)') '$path'
      write (ich, '(a)') '    nrun=1'
      write (ich, '(a,i0)') '    npoint=', nnds
      write (ich, '(a)') '    anopt=10'
      write (ich, '(a)') '    kpush=0.003'
      write (ich, '(a)') '    kpull=-0.015'
      write (ich, '(a)') '    ppull=0.05'
      write (ich, '(a)') '    alp=1.2'
      write (ich, '(a)') '$end'
      close (ich)
      write (jobcall, '(a,i0,a)') ' xtb start.xyz --path end.xyz --input path.inp --etemp ', etemp, '  > path.out 2>/dev/null'
      write (jobcall, '(a)') trim(jobcall)//' && touch xtb_finished'
      !DEL write(jobcall,'(a)') trim(jobcall)//' && cp xtbpath_ts.xyz ts.xyz'
      write (jobcall, '(a)') trim(jobcall)//' && rm xtbpath_?.xyz'
      ! write(*,*) " jobcall is", trim(jobcall)
   end subroutine prepxtbpath

! this routine adds an > after each structure for the neb input guess file
   subroutine prepare_neb_restart(env)
      implicit none
      type(runtypedata) :: env
      character(len=100) :: input_file = "interpolated.xyz"  ! Input file name
      character(len=100) :: output_file = "nebguess.xyz"  ! Output file name
      integer :: nat, nlines  ! Insert character every nth line
      integer :: line_num, status
      integer :: ich1, ich2
      integer :: io
      character(len=1000) :: line
      character :: delimiter = '>'  !

      call rdshort_int("start.xyz", nat)

      nlines = nat + 2
      ! Open input file
      open (newunit=ich1, file=input_file, status='old', action='read', iostat=io)
      if (io .ne. 0) then
         write (*, *) "Error interpolated.xyz file from geodesic not found"
         return
      end if
      ! Open output file
      open (newunit=ich2, file=output_file, status='replace', action='write', iostat=io)
      if (io .ne. 0) then
         write (*, *) "Error opening output file"
         return
      end if
      line_num = 0

      ! Read lines from input file, insert character every nth line, and write to output file
      do
         read (ich1, '(a)', iostat=io) line
         if (io .ne. 0) then
            exit ! End of file
         end if

         line_num = line_num + 1

         write (ich2, '(a)') trim(line)

         if (mod(line_num, nlines) == 0) then
            write (ich2, '(a)') delimiter
         end if
      end do

      ! Close files
      close (ich1)
      close (ich2)
   end subroutine prepare_neb_restart

   subroutine cleanup_nebcalc
      implicit none
      
      call remove("orca_MEP.allxyz")
      call remove("orca.interp")
      call remove("geoout")
      call remove("geo2out")
      call remove("orca.property.txt")
      call remove("orca_initial_path_trj.xyz")
      call remove("orca.final.interp")
      call remove("orca.bibtex")
      call remove("orca.citations.tmp")
      call remove("orca.NEB.log")
      call remove("orca_NEB-HEI_converged.xyz")
      call remove("orca_MEP_ALL_trj.xyz")
      call execute_command_line("rm orca_atom*")
      call execute_command_line("rm orca_im*")
   
   end subroutine cleanup_nebcalc

   subroutine prepcrest_tsgen(env,jobcall)
      implicit none
      character(len=1024), intent(out) :: jobcall
      type(runtypedata) :: env
      integer :: chrg, uhf


      call rdshort_int('.CHRG', chrg)
      call rdshort_int(".UHF", uhf)
      

      call write_cresttoml(env)
      write(jobcall,'(a,i0,a,i0,a)') 'xtb start.xyz --hess --chrg ',chrg,' --uhf ',uhf,' > hess1.out 2> /dev/null'
      write(jobcall,'(a)') trim(jobcall)//' && mv hessian hess1 && mv wbo wbo1'
      write(jobcall,'(a,i0,a,i0,a)') trim(jobcall)//' && xtb end.xyz --hess --chrg ',chrg,' --uhf ',uhf,' > hess2.out 2> /dev/null'
      write(jobcall,'(a)') trim(jobcall)//' && mv hessian hess2 && mv wbo wbo2'
      write(jobcall,'(a,a)') trim(jobcall)//' && crest_tsgen --input input.toml > crest_tsgen.out 2> crest_tsgen_error.out', &
      &  ' && touch crest_tsgen_finished'
      !write(*,*) "jobcall is", trim(jobcall)
   end subroutine prepcrest_tsgen

   subroutine write_cresttoml(env)
      implicit none 
      type(runtypedata) :: env
      integer :: io, ich
      integer :: chrg, uhf


      call rdshort_int('.CHRG', chrg)
      call rdshort_int(".UHF", uhf)
      ! save .UHF and .CHRG which will be deleted by crest
      call move(".UHF", "uhftemp")
      call move(".CHRG", "chrgtemp")
      

      open (newunit=ich, file='input.toml')
      write(ich,'(a)') '#This is a CREST input file'
      write(ich,'(a)') 'input = "start.xyz"'
      write(ich,'(a)') 'runtype = "ancopt"'
      write(ich,'(a,i0)') 'threads =', env%threads
      write(ich,'(a)') ''
      write(ich,'(a)') '[calculation]'
      write(ich,'(a)') 'id = 3'
      write(ich,'(a)') 'maxcycle = 219'
      write(ich,'(a)') ''
      write(ich,'(a)') '[[calculation.level]]'
      write(ich,'(a)') 'method = "ssff"'
      write(ich,'(a)') 'refgeo = "start.xyz"'
      write(ich,'(a)') 'refhessian = "hess1"'
      write(ich,'(a)') 'refwbo = "wbo1"'
      write(ich,'(a)') 'refff = "ff1"'
      write(ich,'(a)') 'useff = false'
      write(ich,'(a,i0)') 'chrg =', chrg
      write(ich,'(a,i0)') 'uhf =', uhf 
      write(ich,'(a)') '[[calculation.level]]'
      write(ich,'(a)') 'method = "ssff"'
      write(ich,'(a)') 'refgeo = "end.xyz"'
      write(ich,'(a)') 'refhessian = "hess2"'
      write(ich,'(a)') 'refwbo = "wbo2"'
      write(ich,'(a)') 'refff = "ff2"'
      write(ich,'(a)') 'useff = false'
      write(ich,'(a,i0)') 'chrg =', chrg
      write(ich,'(a,i0)') 'uhf =', uhf  
      write(ich,'(a)') ''
      write(ich,'(a)') '[[calculation.level]]'
      write(ich,'(a)') 'method = "tsff"'
      write(ich,'(a)') 'refff = "tsff"'
      write(ich,'(a)') 'id1 = 1'
      write(ich,'(a)') 'id2 = 2'
      write(ich,'(a,i0)') 'chrg =', chrg
      write(ich,'(a,i0)') 'uhf =', uhf
      close(ich)



   end subroutine write_cresttoml
end module tsmod
