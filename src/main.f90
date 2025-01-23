!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! QCxMS-2
! author: Johannes Gorges
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
program QCxMS2
   use mcsimu
   use qcxms2_data
   use xtb_mctc_accuracy, only: wp
   use tsmod
   use argparser
   use fragmentation
   use qmmod
   use mcsimu
   use charges
   use plot
   use iomod
   use qcxms_iee
   use cid

   use mctc_env, only: error_type, get_argument!, fatal_error
   use mctc_io, only: structure_type, read_structure, write_structure, &
      & filetype, get_filetype, to_symbol, to_number

   implicit none
   character(len=1024) :: cmd, cdum, pwd, pwd2, dir, tmppath1, tmppath, pairdir, copycall
   character(len=80) :: fname
   character(len=80), allocatable :: fraglist(:), fraglist1(:), fraglist2(:) ! 1 D here p? or p?f? can occur here
   character(len=80), allocatable :: allfrags(:)
   character(len=80) :: startdir ! directory where the first fragment is located input molecule is ''
   character(len=80) :: fout, pattern ! fout and pattern for readout
   character(len=80) :: fraglevel
   character(len=1024) :: cleanupcall, jobcall
   integer :: nargs, i, l, j, k, m, n, io, ich, threads
   integer :: istat, npairs, nsamples, nat
   integer :: ilen
   integer :: iunit
   integer :: nfrags1, nfrags2, nfrags3, nfrags4, ntot, nfrags
   integer :: nfrags5, nfrags6, nfrags7, nfrags8
   integer :: iseed(1)
   real(wp) :: ip, pf1, pf2, pfrag, dum, edum
   real(wp), allocatable :: eiee(:), piee(:)
   real(wp) :: t1, w1, t2, w2, rrho, sumreac, dethr ! energy threshold for subsequent fragmentation
   logical :: ex, plotlog, ldum, there
   type(runtypedata) :: env
   type(timer):: tim

   call qcxms2_header()

   call disclamer(iunit)
   write (*, '(/,1x,a)') 'Command line input:'
   call get_command(cmd)
   write (*, '(1x,a,a,/)') '> ', trim(cmd)

   call timing(t1, w1)
   call tim%init(20)

   call getcwd(pwd)
   env%startdir = pwd
   env%path = pwd ! env%path always starting path for new fragmentation

   !set defaults and parse commandline arguments
   nargs = iargc()

   l = len_trim(cmd)
   allocate (arg(nargs), source=repeat(' ', l))
   do i = 1, nargs
      call getarg(i, arg(i))
   end do
   call parseflags(env, arg, nargs, iseed)
   deallocate (arg)

   if (env%logo) then
      call qcxms2_logo
   end if
   
   call citation(env)

   !> inialize random numbers
   !>> if not exlicitly set, use true random

   if (iseed(1) == 0) then
      call random_seed()
   else
      !>> if explicitly set, use seed number from input
      call random_seed(put=iseed)
   end if

   ! check if programs and paths are properly set
   ! TODO check if settings are consistent, and meaningful
   call check_settings(env)
   ! check if all requested programs are in PATH
   call check_progs(env)


   ! just for testing purposes here:
   if (env%ircrun) then
      call findirc(env, i, dum)
      write (*, *) "IRC is", dum!
      STOP
   end if

   if (env%picktsrun) then
      call pickts(env, ex)
      STOP
   end if

   ! convert pthr from percent in number
   env%pthr = env%pthr/100.0_wp*env%nsamples
   env%index = ''

   ! only plot
   if (env%oplot) GOTO 42

   if (env%esim) then
      env%cores = 1 ! just to be sure
      call setompthreads(env, 1)
      write (*, *) "Computing reaction network with Monte Carlo Simulation"
      plotlog = .true.
      nsamples = env%nsamples
      call rdshort_real("ip_"//trim(env%iplevel), ip)
      write (*, *) "read IP from ip_"//trim(env%iplevel)//": ", ip
      env%ehomo = ip ! Just use Koopmans theorem for now
      allocate (eiee(nsamples), piee(nsamples))
      call getenergy(env, nsamples, eiee, piee, plotlog)
      call mcesim(env, eiee, piee)
      GOTO 42
   end if

   !TODO cleanup the code here
   inquire (file='QCxMS2_finished', exist=ex)
   if (ex) then
      write (*, *) "Files of prior QCxMS2 calculation found: Only missing components are computed"
      call cleanup_qcxms2_for_restart(env)
      env%restartrun = .true.
   else
      env%restartrun = .false.
   end if

   fname = trim(env%infile)

   call setompthreads(env, 1)
   ! restart run
   inquire (file='ip_'//trim(env%iplevel), exist=ex)
   if (ex) then
      call rdshort_real('ip_'//trim(env%iplevel), ip)
   else
      write (*, *) "Optimizing starting structure at the ", trim(env%geolevel), " level"
      call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .false., 0)
      if (.not. there) call execute_command_line(trim(jobcall), exitstat=io)
      write (*, *) "Compute IP at level ", trim(env%iplevel)
      call getip1(env, fname, ip)
   end if

   write (*, *)
   write (*, '(a,f10.3,a)') "IP of input molecule is", ip, " eV"

   env%ehomo = ip ! Just use Koopmans theorem for now

   ! optimize and compute energies of starting fragment
   call calc_startfragment(env, fname, 1)

   call copy(fname, 'fragment.xyz')
   !generate energy distribution for IEE (impac excess energy)
   plotlog = .true.
   nsamples = env%nsamples

   call getenergy(env, nsamples, eiee, piee, plotlog) ! we use one IEE for complete run

   write (*, *)
   write (*, *) "--------------------------------"
   write (*, *) "| Starting first fragmentation |"
   write (*, *) "--------------------------------"
   write (*, *)
   call tim%start(1, 'first fragmentation')
   startdir = ''
   allocate (fraglist1(1))
   fraglist1 = ['']
   allfrags = ['']
   call calcfragments(env, fname, eiee, piee, 1, startdir, fraglist, nfrags)
   ntot = nfrags
   if (env%printlevel .eq. 3) write (*, *) "ntot and nfrags is:", ntot, nfrags
   nfrags1 = nfrags
   fraglist1 = [fraglist1, fraglist]
   allfrags = [allfrags, fraglist]

   call tim%stop(1)

   if (nfrags1 .eq. 0) then
      write (*, *) "No fragments were generated, something went wrong ..."
      GOTO 100
   end if
   ! GOTO statements are ugly but here its the easiest way
   if (env%nfrag .lt. 1) GOTO 100

   fraglist2 = ['']
   nfrags2 = 0

   ! MAIN FRAGMENTATION LOOP
   do k = 2, env%nfrag + 1
      write (*, *)
      write (*, *) "--------------------------------------------------"
      write (*, *) "|Starting fragmentation ", k, " for ", nfrags1, " fragments"
      write (*, *) "--------------------------------------------------"
      write (*, *)
      write (fraglevel, '(a,i0)') "Fragmentation level ", k
      call tim%start(k, trim(fraglevel))
      do i = 2, nfrags1 + 1 ! first element of array always empty
         inquire (file=trim(fraglist1(i))//"/pfrag", exist=ex)
         if (.not. ex) then
            write (*, *) "Fragment ", trim(fraglist1(i)), " is not fragmented further due to missing pfrag file" ! todo there is still a problem with this
            pfrag = 0.0_wp
         else
            call rdshort_real(trim(fraglist1(i))//"/pfrag", pfrag)
         end if
         if (pfrag .ge. env%pthr) then
            call chdir(trim(fraglist1(i)))
            call getcwd(env%path)
            inquire (file='isomer.xyz', exist=ex) !TODO make code nicer
            if (ex) then
               fname = 'isomer.xyz'
            else
               fname = 'fragment.xyz'
            end if
            startdir = fraglist1(i)
            write (*, *) "------------------------------------"
            write (*, *) "| Starting fragmentation of: ", trim(fraglist1(i))
            write (*, *) "------------------------------------"
            call calc_startfragment(env, fname, k)
            call calcfragments(env, fname, eiee, piee, k, startdir, fraglist, nfrags)
            call chdir(trim(env%startdir))
            ntot = ntot + nfrags
            write (*, *) "newly generated products: ", nfrags, " total number of products is now: ", ntot
            nfrags2 = nfrags2 + nfrags

            ! why do we need this?
            if (allocated(fraglist)) then
               fraglist2 = [fraglist2, fraglist]
               allfrags = [allfrags, fraglist]
            else
               !write(*,*) "fraglist not allocated"
            end if
         else
            write (*, '(a,a,a,f3.1,a)') "Fragment ", &
            & trim(fraglist1(i)), " is not fragmented further due to too low propability (", &
            & pfrag/env%nsamples*100, ")"
         end if
      end do
      call tim%stop(k)
      if (nfrags2 .eq. 0) then
         write (*, *) "No new fragments were generated"
         GOTO 100
      end if
      fraglist1 = fraglist2
      fraglist2 = ['']
      nfrags1 = nfrags2
      nfrags2 = 0
   end do

100 open (newunit=ich, file='allfrags.dat') ! write all fragments to file TODO maybe add info about fragmentation level
   write (ich, *) ntot
   do i = 2, ntot + 1 ! first element of array always empty
      write (ich, *) trim(allfrags(i))
   end do
   close (ich)

42 if (env%oplot .or. env%esim) then
      call chdir(env%startdir) !Just to be safe here
      call rdshort_int('allfrags.dat', ntot)
      if (ntot .eq. 0) then
         write (*, *) "No allfrags.dat file found, something went wrong ..."
         error stop
      end if
      allocate (allfrags(ntot + 1))
      open (newunit=ich, file='allfrags.dat')
      read (ich, *) ntot
      do i = 2, ntot + 1
         read (ich, *) allfrags(i)
      end do
      close (ich)
   end if

   call getpeaks(env, ntot, allfrags)
   call eval_timer(tim)

   !TODO  sometimes eval_timer crashes, this would be the fallback option
   call timing(t2, w2)
   ! write (*, '(/,'' wall time (seconds)'',F10.2  )') (w2 - w1)
   write (*, *)
   write (*, *) 'QCxMS2 terminated normally.'
   if (env%printlevel .eq. 3) then
      write (*, *) "Detailed Timings:"
      write (*,'(a,f10.1,a)') "time spent for Crest MSREACT: ", env%tcrest, "s"
      write (*,'(a,f10.1,a)') "time spent for IP calculations: ", env%tip, "s"
      write (*,'(a,f10.1,a)') "time spent for path calculations: ", env%tpath, "s"
   end if
   call touch('QCxMS2_finished')

end program QCxMS2
