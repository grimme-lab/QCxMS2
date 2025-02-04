! module to assign statistical charges to fragments
module charges
   use iomod
   use qcxms2_data
   use qmmod
   use utility
   use xtb_mctc_convert
   use xtb_mctc_constants
   implicit none

contains
! assign charges according to IPs from Delta SCF
! close IPs (defined by ipthr) are refined at higher level of theory
   subroutine assigncharges(env, npairs_in, npairs_out, fragdirs_in, fragdirs_out)
      implicit none
      type(runtypedata) :: env
      type(timer):: tim
      integer, intent(in) :: npairs_in
      character(len=80) :: dir
      character(len=80)   ::  fname
      character(len=80)  :: fragdirs_in(npairs_in, 3)
      character(len=80), allocatable, intent(out)   :: fragdirs_out(:, :)
      character(len=80), allocatable :: dirs(:), dirs0(:), dirs1(:)
      character(len=1024) :: jobcall, cleanupcall0, cleanupcall
      character(len=80) :: fout, pattern ! fout and pattern for readout
      character(len=80) :: ip3level ! level for repair of negative IPs, default is GFN1-xTB (robust) or GFN2-xTB if iplevel1 is GFN1-xTB
      integer :: i, njobs, njobs0, njobs1, j, io
      integer, intent(out) :: npairs_out
      logical :: ex, failed, there, ipok
      logical, allocatable :: closeips(:)
      real(wp) :: ip, ip1, ip2, qipb, diffip, ipthr
      real(wp) :: ips(2) ! TODO can modify this later to be more general for multiple charged fragments
      real(wp) :: e_charged, e_neutral, edum
      integer :: q1, q2

      ipthr = 2.0_wp ! TODO tune this values in eV threshold for difference in IPs to use higher level

      write (*, *) "Computing statistical charges via Delta SCF procedure"
      ! for now parallelization of charge1 and neut should be enough
      ! prepare calculations of charge1
      ! of course we only need to do this for for real fragmentpairs
      ! and not rearranged structures

      ! call calc_h_atoms(fragdirs_in,npairs_in,npairs_out,fragdirs_out)
      ! first find out how many real fragment pairs there are

      njobs = 0
      fname = 'fragment.xyz'
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then ! does fragment exist?
            do j = 2, 3
               call rdshort(trim(fragdirs_in(i, j))//'/ip_'//trim(env%iplevel), ip)
               if (env%restartrun .and. abs(ip) .gt. 1.0e-2) cycle
               call chdir(trim(fragdirs_in(i, j)))
               call prepip(env, fname, fragdirs_in(i, j), env%iplevel, dirs, jobcall, cleanupcall, fout, pattern, njobs)
               call chdir(trim(env%path))
            end do
         end if
      end do

      ! we first have to figure out the number of jobs
      call setompthreads(env, njobs)
      deallocate (dirs)
      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then ! does fragment exist?
            do j = 2, 3
               call rdshort(trim(fragdirs_in(i, j))//'/ip_'//trim(env%iplevel), ip)
               if (env%restartrun .and. abs(ip) .gt. 1.0e-2) cycle
               call chdir(trim(fragdirs_in(i, j)))
               call prepip(env, fname, fragdirs_in(i, j), env%iplevel, dirs, jobcall, cleanupcall0, fout, pattern, njobs)
               call chdir(trim(env%path))
            end do
         end if
      end do

      call setompthreads(env, njobs)
      write (*, *) "starting ", njobs, trim(env%iplevel), " calculations in parallel"
      call omp_samejobcall(njobs, dirs, jobcall)  

      dirs0 = dirs
      njobs0 = njobs

      ! determine number of jobs
      njobs = 0
      deallocate (dirs)
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            do j = 2, 3
               call rdshort(trim(fragdirs_in(i, j))//'/ip_'//trim(env%iplevel), ip)
               if (env%restartrun .and. abs(ip) .gt. 1.0e-2) cycle
               call chdir(trim(fragdirs_in(i, j)))
               call checkip(env, fname, fragdirs_in(i, j), env%iplevel, dirs, jobcall, cleanupcall0, fout, pattern, njobs)
               call chdir(trim(env%path))
            end do
         end if
      end do

      ! restart failed calculations
      call setompthreads(env, njobs)
      ! now actually prepare the jobs
      njobs = 0
      deallocate (dirs)
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            do j = 2, 3
               call rdshort(trim(fragdirs_in(i, j))//'/ip_'//trim(env%iplevel), ip)
               if (env%restartrun .and. abs(ip) .gt. 1.0e-2) cycle
               call chdir(trim(fragdirs_in(i, j)))
               call checkip(env, fname, fragdirs_in(i, j), env%iplevel, dirs, jobcall, cleanupcall0, fout, pattern, njobs)
               call chdir(trim(env%path))
            end do
         end if
      end do

      write (*, *) "starting ", njobs, trim(env%iplevel), "restart IP calculations in parallel"
      call omp_samejobcall(njobs, dirs, jobcall)

      ! now readout values

      ! refine close IPs at higher level of theory calculations
      njobs = 0
      deallocate (dirs)
      allocate (dirs(1)) ! first element always empty
      allocate (closeips(npairs_in))
      closeips = .false.
      dirs = ['']

      npairs_out = npairs_in

      do i = 1, npairs_in
         if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
            do j = 2, 3
               call rdshort(trim(fragdirs_in(i, j))//'/ip_'//trim(env%iplevel), ips(j - 1))
               if (env%restartrun .and. abs(ip) .gt. 1.0e-2) then
                  ipok = .true.
                  cycle
               end if
               call chdir(trim(fragdirs_in(i, j)))
               call readoutip(env, fname, env%iplevel, fout, pattern, ips(j - 1), ipok)
               call chdir(trim(env%path))
               if (.not. ipok) then ! todo maybe add a retry here
                  npairs_out = npairs_out - 1
                  fragdirs_in(i, 1) = ''
                  fragdirs_in(i, 2) = ''
                  fragdirs_in(i, 3) = ''
                  exit
               end if
            end do
            if (ipok) then
               call boltzmann(env, ips(1), ips(2), fragdirs_in(i, 1))
            end if
            if (env%ip2level .ne. '' .and. trim(env%ip2level) &
            & .ne. trim(env%iplevel) .and. ipok .and. abs(ips(1) - ips(2)) .lt. ipthr) then
               closeips(i) = .true.
               do j = 2, 3
                  call rdshort(trim(fragdirs_in(i, j))//'/ip_'//trim(env%ip2level), ip)
                  if (env%restartrun .and. abs(ip) .gt. 1.0e-2) cycle
                  call chdir(trim(fragdirs_in(i, j)))
                  call prepip(env, fname, fragdirs_in(i, j), env%ip2level, dirs, jobcall, cleanupcall, fout, pattern, njobs)
                  call chdir(trim(env%path))
               end do
            else
               closeips(i) = .false.
            end if
         end if
      end do

      ! refine closeips at higher level of theory
      if (env%ip2level .ne. '' .and. trim(env%ip2level) .ne. trim(env%iplevel) .and. njobs .gt. 0) then

         call setompthreads(env, njobs)
         njobs = 0
         deallocate (dirs)
         allocate (dirs(1)) ! first element always empty
         dirs = ['']

         do i = 1, npairs_in
            if (closeips(i)) then
               do j = 2, 3
                  call chdir(trim(fragdirs_in(i, j)))
                  call prepip(env, fname, fragdirs_in(i, j), env%ip2level, dirs, jobcall, cleanupcall, fout, pattern, njobs)
                  call chdir(trim(env%path))
               end do
            end if
         end do

         call setompthreads(env, njobs)
         write (*, *) "starting ", njobs, trim(env%ip2level), "refinement IP calculations in parallel"
         call omp_samejobcall(njobs, dirs, jobcall)

         njobs1 = njobs
         dirs1 = dirs
         ! restart close ips
         njobs = 0
         deallocate (dirs)
         allocate (dirs(1)) ! first element always empty
         dirs = ['']

         do i = 1, npairs_in
            if (.not. closeips(i)) cycle
            if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
               do j = 2, 3
                  call rdshort(trim(fragdirs_in(i, j))//'/ip_'//trim(env%ip2level), ip)
                  if (env%restartrun .and. abs(ip) .gt. 1.0e-2) cycle
                  call chdir(trim(fragdirs_in(i, j)))
                  call checkip(env, fname, fragdirs_in(i, j), env%ip2level, dirs, jobcall, cleanupcall, fout, pattern, njobs)
                  call chdir(trim(env%path))
               end do
            end if
         end do

         call setompthreads(env, njobs)
         njobs = 0
         deallocate (dirs)
         allocate (dirs(1)) ! first element always empty
         dirs = ['']

         do i = 1, npairs_in
            if (.not. closeips(i)) cycle
            if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
               do j = 2, 3
                  call rdshort(trim(fragdirs_in(i, j))//'/ip_'//trim(env%ip2level), ip)
                  if (env%restartrun .and. abs(ip) .gt. 1.0e-2) cycle
                  call chdir(trim(fragdirs_in(i, j)))
                  call checkip(env, fname, fragdirs_in(i, j), env%ip2level, dirs, jobcall, cleanupcall, fout, pattern, njobs)
                  call chdir(trim(env%path))
               end do
            end if
         end do

         write (*, *) "starting ", njobs, trim(env%ip2level), " restart IP2 calculations in parallel"
         call omp_samejobcall(njobs, dirs, jobcall)

         ! readout closeips
         do i = 1, npairs_in
            if (.not. closeips(i)) cycle
            if (index(fragdirs_in(i, 3), 'p') .ne. 0) then
               do j = 2, 3
                  call rdshort(trim(fragdirs_in(i, j))//'/ip_'//trim(env%ip2level), ip)
                  if (env%restartrun .and. abs(ip) .gt. 1.0e-2) cycle
                  call chdir(trim(fragdirs_in(i, j)))
                  call readoutip(env, fname, env%ip2level, fout, pattern, ips(j - 1), ipok)
                  call chdir(trim(env%path))
                  if (.not. ipok) then ! fall back on ip1level, if this failed it would have been sorted out before
                     call rdshort('ip_'//trim(env%iplevel), ips(j - 1))
                  end if
               end do
               call boltzmann(env, ips(1), ips(2), fragdirs_in(i, 1))
            end if
         end do
         ! cleanup ip2 level Calculations
         call omp_samejobcall(njobs1, dirs1, cleanupcall, .false.)
      end if

      ! cleanup ip1 level calculations
      call omp_samejobcall(njobs0, dirs0, cleanupcall0, .false.)

      ! now we sort out failed IP calculations
      call sortout0elements(npairs_in, npairs_out, fragdirs_in, fragdirs_out, env%removedirs)

      !now write .CHRG files
      do i = 1, npairs_out
         if (index(fragdirs_out(i, 3), 'p') .ne. 0) then ! only for fragments important
            call chdir(fragdirs_out(i, 2))
            call rdshort_real('statchrg', qipb)
            if (qipb .ge. 0.5_wp) then !if exact 0.5 then just give first fragment charge
               q1 = env%chrg
               q2 = 0
            else
               q1 = 0
               q2 = env%chrg
            end if
            call wrshort_int('.CHRG', q1)
            call chdir(trim(env%path))
            call chdir(fragdirs_out(i, 3))
            call wrshort_int('.CHRG', q2)
            call chdir(trim(env%path))
         end if
      end do
      deallocate (dirs, closeips)

   end subroutine assigncharges

! get ionisation potential via Delta SCF for first one we parallellize this explicitly

   subroutine getip1(env, fname, ip)
      implicit none
      type(runtypedata) :: env
      integer :: io
      real(wp), intent(out) :: ip
      character(len=80), intent(in)   ::  fname
      character(len=80), allocatable ::  dirs(:)
      character(len=80) :: dir
      character(len=1024) :: jobcall, cleanupcall
      integer :: njobs
      character(len=80) :: fout, pattern ! fout and pattern for readout
      logical :: ipok

      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']
      dir = '.'
      call prepip1(env, fname, dir, env%iplevel, dirs, jobcall, cleanupcall, fout, pattern, njobs)

      ! we first have to figure out the number of jobs
      call setompthreads(env, njobs)
      deallocate (dirs)
      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']

      call prepip1(env, fname, dir, env%iplevel, dirs, jobcall, cleanupcall, fout, pattern, njobs)
      call omp_samejobcall(njobs, dirs, jobcall)
      deallocate (dirs)

      njobs = 0
      allocate (dirs(1)) ! first element always empty
      dirs = ['']
      call checkip(env, fname, dir, env%iplevel, dirs, jobcall, cleanupcall, fout, pattern, njobs)
      call setompthreads(env, njobs)
      call omp_samejobcall(njobs, dirs, jobcall)
      call readoutip(env, fname, env%iplevel, fout, pattern, ip, ipok)
      if (.not. ipok) then
         stop "calculation of charged molecule for first IP failed, consider to change iplevel or check the calculation"
      end if

   end subroutine getip1

! get ionisation potential via Delta SCF for first IP for later we already have charge0 dirs with fragments
   subroutine prepip1(env, fname, dir, level, dirs, jobcall, cleanupcall, fout, pattern, njobs)
      implicit none
      type(runtypedata) :: env
      character(len=*), intent(in) :: dir, level
      character(len=1024) :: thisdir
      character(len=80)   ::  subdir
      integer :: io, chrg
      integer, intent(inout) :: njobs
      real(wp) :: e_charged, e_neutral
      character(len=1024), intent(out) :: jobcall, cleanupcall
      character(len=80), allocatable, intent(inout) :: dirs(:)
      logical :: there
      character(len=*), intent(in)   ::  fname
      character(len=80) :: fout, pattern ! fout and pattern for readout

      call getcwd(thisdir)
      ! make SP with charged molecule
      call prepqm(env, fname, level, 'sp', jobcall, fout, pattern, cleanupcall, there, .false., env%chrg)
      if (.not. there) then
         njobs = njobs + 1
         call append_char(dirs, trim(dir))
      end if
      ! make SP with neutral molecule
      write (subdir, '(a,i0)') "charge", 0
      io = makedir(trim(subdir))
      call copysub(fname, subdir)
      call chdir(subdir)
      ! give here charge=0
      call prepqm(env, fname, level, 'sp', jobcall, fout, pattern, cleanupcall, there, .false., 0)
      if (.not. there) then
         call append_char(dirs, trim(dir)//'/charge0')
         njobs = njobs + 1
      end if
      call chdir(thisdir)

   end subroutine prepip1

   ! get ionisation potential via Delta SCF
   subroutine prepip(env, fname, dir, level, dirs, jobcall, cleanupcall, fout, pattern, njobs)
      implicit none
      type(runtypedata) :: env
      character(len=*), intent(in) :: dir, level
      character(len=1024) :: thisdir
      character(len=80)   ::  subdir
      integer :: io, chrg
      integer, intent(inout) :: njobs
      real(wp) :: e_charged, e_neutral
      character(len=1024), intent(out) :: jobcall, cleanupcall
      character(len=80), allocatable, intent(inout) :: dirs(:)
      logical :: there
      character(len=*), intent(in)   ::  fname
      character(len=80) :: fout, pattern ! fout and pattern for readout

      call getcwd(thisdir)
      ! make SP with charged molecule
      call prepqm(env, fname, level, 'sp', jobcall, fout, pattern, cleanupcall, there, .false., env%chrg)
      if (.not. there) then
         njobs = njobs + 1
         call append_char(dirs, trim(dir))
      end if
      ! make SP with neutral molecule
      write (subdir, '(a,i0)') "charge", 0
      call chdir(subdir)
      ! give here charge=0
      call prepqm(env, fname, level, 'sp', jobcall, fout, pattern, cleanupcall, there, .false., 0)
      if (.not. there) then
         call append_char(dirs, trim(dir)//'/charge0')
         njobs = njobs + 1
      end if
      call chdir(thisdir)

   end subroutine prepip

   ! get ionisation potential via Delta SCF
   subroutine checkip(env, fname, dir, level, dirs, jobcall, cleanupcall, fout, pattern, njobs)
      implicit none
      type(runtypedata) :: env
      character(len=80), intent(in) :: dir, level
      character(len=80) :: subdir
      character(len=1024) :: thisdir
      integer :: io, chrg
      integer, intent(inout) :: njobs
      real(wp) :: e_charged, e_neutral, edum
      character(len=1024), intent(out) :: jobcall, cleanupcall
      character(len=80), allocatable, intent(inout) :: dirs(:)
      logical :: failed, there
      character(len=80), intent(in)   ::  fname
      character(len=80) :: fout, pattern ! fout and pattern for readout

      call getcwd(thisdir)
      ! make SP with charged molecule
      call readoutqm(env, fname, env%iplevel, 'sp', fout, pattern, edum, failed)
      if (failed) then
         call prepqm(env, fname, env%iplevel, 'sp', jobcall, fout, pattern, cleanupcall, there, .true.)
         call append_char(dirs, trim(dir))
         njobs = njobs + 1
      end if
      ! make SP with neutral molecule
      write (subdir, '(a,i0)') "charge", 0
      call chdir(subdir)
      ! give here charge=0
      call readoutqm(env, fname, env%iplevel, 'sp', fout, pattern, edum, failed)
      if (failed) then
         call prepqm(env, fname, env%iplevel, 'sp', jobcall, fout, pattern, cleanupcall, there, .true.)
         call append_char(dirs, trim(dir)//'/charge0')
         njobs = njobs + 1
      end if
      call chdir(thisdir)

   end subroutine checkip

   ! readout IPs
   subroutine readoutip(env, fname, level, fout, pattern, ip, ipok)
      implicit none
      type(runtypedata) :: env
      character(len=80), intent(in) :: level
      character(len=80) :: subdir
      character(len=1024) :: thisdir
      integer :: io, chrg
      real(wp) :: e_charged, e_neutral
      character(len=80), intent(in)   ::  fname
      character(len=80), intent(in) :: fout, pattern ! fout and pattern for readout
      logical, intent(out) :: ipok ! if ip is ok
      real(wp), intent(out) :: ip
      logical :: failed

      call getcwd(thisdir)
      !
      call readoutqm(env, fname, level, 'sp', fout, pattern, e_charged, failed)
      if (failed) then
         ipok = .false.
         return
      end if
      ! make SP with neutral molecule
      write (subdir, '(a,i0)') "charge", 0
      call chdir(subdir)
      ! give here charge=0
      call readoutqm(env, fname, level, 'sp', fout, pattern, e_neutral, failed)
      if (failed) then
         ipok = .false.
         return
      end if
      call chdir(thisdir)
      ip = (e_charged - e_neutral)*autoev ! total energies in au, all relative energies in eV
      if (ip .lt. 0 .and. env%chrg .ge. 0) then
         write (*, *) "warning negative IP"
         ipok = .false.
         return
      else
         ipok = .true.
      end if

      call writeip(level, e_neutral, e_charged, ip)
   end subroutine readoutip

! calculate and write ionisation potential via Delta SCF
   subroutine writeip(level, e_neutral, e_charged, ip)
      implicit none
      character(len=80) :: subdir, thisdir
      character(len=80) :: level
      integer :: io, chrg
      real(wp) :: e_charged, e_neutral
      real(wp), intent(out) :: ip

      ip = (e_charged - e_neutral)*autoev ! total energies in au, all relative energies in eV
      call wrshort_real('ip_'//trim(level), ip)
   end subroutine writeip

! determine statistical charge via Boltzmann Distribution at  for two fragments
! input two IPs in eV statchrg file is written
! TODO make nicer with general routine for boltzmann distribution in utility
   subroutine boltzmann(env, ip1, ip2, dir)
      implicit none
      type(runtypedata) :: env
      real(wp), intent(in)  :: ip1, ip2
      real(wp) :: esum, statchrg1, statchrg2, f, temp, q1, q2
      character(len=80), intent(in) :: dir
      character(len=512) :: tmpdir
      character(len=2024) :: thisdir
      integer :: i, ich

      call getcwd(thisdir)
      temp = env%temp
      f = temp*kB*autoev
      q1 = exp(-ip1/f)
      q2 = exp(-ip2/f)
      esum = q1 + q2
      statchrg1 = q1/esum
      statchrg2 = q2/esum

      write (tmpdir, '(a)') trim(dir)//'f1'
      call chdir(trim(tmpdir))
      call wrshort_real('statchrg', statchrg1)
      call chdir(trim(thisdir))
      write (tmpdir, '(a)') trim(dir)//'f2'
      call chdir(trim(tmpdir))
      call wrshort_real('statchrg', statchrg2)
      call chdir(trim(thisdir))

   end subroutine boltzmann

end module charges
