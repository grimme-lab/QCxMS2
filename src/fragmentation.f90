module fragmentation
   use qcxms2_data
   use iomod
   use tsmod
   use mcsimu
   use charges
   use plot
   use reaction
   use xtb_mctc_convert
   use structools

   implicit none
contains

   subroutine calc_startfragment(env, fname, nfragl)
      implicit none
      character(len=1024) :: jobcall, cleanupcall
      character(len=*) :: fname
      character(len=80) :: fout, pattern ! fout and pattern for readout
      character(len=80) :: job
      type(runtypedata) :: env
      integer :: io, chrg, uhf, statchrg
      integer :: nfragl ! fragmentation level
      real(wp) :: edum
      logical :: ex, ldum, there, failed1, failed2

      call setompthreads(env, 1)

      ! set thermal temperature for one starting fragment and subsequent fragments to the same value
      ! important for boltzmann weighting of charge assignment and if applicable for thermo contribution
      call settemp(env, fname)

      ! TODO FIXME implement charge assignment here
      ! Assumption is currently that only visible fragments (with input charge env%chrg) are fragmented further
      ! (e.g. a neutral fragment will not fragment to a charged fragment further)
      call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .false., env%chrg)
      if (.not. there) call execute_command_line(trim(jobcall), exitstat=io)
      call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed1)
      if (failed1) then
         call prepqm(env, fname, env%geolevel, 'opt', jobcall, fout, pattern, cleanupcall, there, .true., env%chrg)
         call execute_command_line(trim(jobcall), exitstat=io)
         call readoutqm(env, fname, env%geolevel, 'opt', fout, pattern, edum, failed2)
      end if
      if (failed1 .and. failed2) then
         if (nfragl .eq. 1) stop "first geometry opt failed, check input"
         ! just skip fragment
         write (*, *) "geometry optimization failed for this fragment, continuing with next fragment"
         return
      end if

      call execute_command_line(trim(cleanupcall), exitstat=io)

      !get first single point
      if (env%exstates .gt. 0) then
         job = 'tddft'
      else
         job = 'sp'
      end if
      call prepqm(env, fname, env%tslevel, job, jobcall, fout, pattern, cleanupcall, there, .false., env%chrg)

      if (.not. there) call execute_command_line(trim(jobcall), exitstat=io)
      call readoutqm(env, fname, env%tslevel, job, fout, pattern, edum, failed1)
      if (failed1) then
         call prepqm(env, fname, env%tslevel, job, jobcall, fout, pattern, cleanupcall, there, .true., env%chrg)
         call execute_command_line(trim(jobcall), exitstat=io)
         call readoutqm(env, fname, env%tslevel, job, fout, pattern, edum, failed2)
      end if
      if (failed1 .and. failed2) then
         if (nfragl .eq. 1) stop "first single point calculation failed, aborting run, something is wrong with the input geometry" ! this is ok, as we already calculated this fragment if we fragmentate it further, so this
         write (*, *) "single point calculation failed for this fragment, continuing with next fragment"
         return
      end if

      call execute_command_line(trim(cleanupcall), exitstat=io)

      if (env%bhess) then
         call prepqm(env, fname, env%geolevel, 'bhess', jobcall, fout, pattern, cleanupcall, there, .false., env%chrg)
         if (.not. there) call execute_command_line(trim(jobcall), exitstat=io)
         call readoutqm(env, fname, env%geolevel, 'bhess', fout, pattern, edum, failed1)
         if (failed1) then
            call prepqm(env, fname, env%geolevel, 'bhess', jobcall, fout, pattern, cleanupcall, there, .true., env%chrg)
            call execute_command_line(trim(jobcall), exitstat=io)
            call readoutqm(env, fname, env%geolevel, 'bhess', fout, pattern, edum, failed2)
         end if
         if (failed1 .and. failed2) stop "first SPH calucation failed, aborting run, something is wrong with the input geometry"
         call execute_command_line(trim(cleanupcall), exitstat=io)
      end if
   
   end subroutine calc_startfragment
! new fragmentation
   subroutine calcfragments(env, fname, eiee, piee, nfragl, startdir, fraglist, nfrags)
      implicit none
      character(len=1024) :: Jobcall
      character(len=80) :: fname
      character(len=80)   :: dir, dir2, dir1
      character(len=80) :: sortlevel
      character(len=80) :: startdir ! directory where the first fragment is located input molecule is ''
      character(len=80), allocatable   :: fragdirs(:, :), fragdirs_in(:, :)
      character(len=80), allocatable   :: imp_fragdirs(:, :) !
      integer :: impfrags
      character(len=80), allocatable, intent(out) :: fraglist(:)
      integer :: npairs, i, j, jj, io, npairs0, nfrags
      integer :: nfragl ! fragmentation level
      real(wp), intent(in) :: eiee(:), piee(:)
      real(wp), allocatable :: etots(:), e_pairs(:)
      type(runtypedata) :: env
      character(len=80) :: pwd
      logical :: ex, ex2
      logical :: is_atom
      real(wp) :: t1, w1, t2, w2
      logical :: restart

      ! if fragment is only one atome, it cannot be fragmented
      call checkatom(fname, is_atom)
      if (is_atom) then
         if (env%printlevel .eq. 3) then
            write (*, *) "fragment ", trim(startdir), " is only one atom, cannot be further fragmented"
         end if
         nfrags = 0
         return
      end if

      ! check if fragments already generated
      inquire (file='msreact.out', exist=ex)
      ! check if this fragmentation went through already
      if (ex) then
         call rdshort_int('npairs', npairs) ! CARE, we cannot delete directories of fragments here if we take npairs instead of npairs2
         write (*, *) "already generated fragment pairs", npairs
         ! lets just keep msreact.out for now, because we need the reaction energies
         call readfraglist(env, fname, npairs, nfragl, nfrags, fragdirs, etots)
         ! write(*,*) "fragdirs2 are ", fragdirs(:,1)
         ! because we moved all fragments back to the main directory
         ! lets just move the fragments back to fragdir
         if (nfragl .gt. 1) then
            do i = 1, npairs
               if (index(fragdirs(i, 3), 'p') .ne. 0) then
                  jj = 3 ! if fragmentpair
               else
                  jj = 1 ! if isomer
               end if
               do j = 1, jj
                  write (dir2, '(a)') trim(startdir)//trim(fragdirs(i, j))
                  write (jobcall, '(a)') 'mv  '//trim(env%startdir)//'/'//trim(dir2)
                  write (jobcall, '(a)') trim(jobcall)//' '//trim(env%startdir)//'/'//trim(startdir)//'/'//trim(fragdirs(i, j))
                  call execute_command_line(trim(jobcall), exitstat=io)
               end do
            end do
         end if
         restart = .true.
      else
         call timing(t1, w1)
         call genfragments(env, fname, npairs, nfragl)
         call timing(t2, w2)
         env%tcrest = env%tcrest + (w2 - w1)
         write (*,'(a,f10.1,a)') "Time spent for MSREACT calculation: ", (w2 - w1), "s"
         call readfraglist(env, fname, npairs, nfragl, nfrags, fragdirs, etots)
         restart = .false.
      end if

      ! topocheck bugy and not that helpful, dont want gff geometries
      ! optimize products first with GFF to prevent topology changes upon later  optimization
      if (.not. restart) then !> of course only, if structures are not already optimized
         npairs0 = npairs
         fragdirs_in = fragdirs
         deallocate (fragdirs)
         call gfftopocheck(env, npairs0, npairs, fragdirs_in, fragdirs)
         deallocate (fragdirs_in)
      end if

      ! optimize all fragments at charge0 for vertical IPs, except if hotip is set
      if (.not. env%hotip) then
         npairs0 = npairs
         fragdirs_in = fragdirs
         deallocate (fragdirs)
         ! This was used vor vertical IPs, but now we use the relaxed IPs
         !call optfragmentsold(env, npairs0, npairs, fragdirs_in, fragdirs, 0)
         ! currently relaxed IPs are set (meaning optimization at respective charges)
         call optfragments(env, npairs0, npairs, fragdirs_in, fragdirs)
         deallocate (fragdirs_in)
      end if

      !> assign charges to fragments
      npairs0 = npairs
      fragdirs_in = fragdirs
      deallocate (fragdirs)
      call timing(t1, w1)
      call assigncharges(env, npairs0, npairs, fragdirs_in, fragdirs) ! get charge assign    also writes right .CHRG files also for isomers
      deallocate (fragdirs_in)
      call timing(t2, w2)
      env%tip = env%tip + (w2 - w1)

      write (*,'(a,f10.1,a)') "Time spent for IP calculation: ", (w2 - w1), "s"

      ! but now we have to optimize the fragments again, because the charges are different
      ! optimize all fragments at assigned charge

      ! check here for convergence maybe? and topology change
      npairs0 = npairs
      fragdirs_in = fragdirs
      deallocate (fragdirs)
      !call optproductsold(env, npairs0, npairs, fragdirs_in, fragdirs)
      call optproducts(env, npairs0, npairs, fragdirs_in, fragdirs)
      deallocate (fragdirs_in)

      ! now compute reaction energy for all fragments
      ! we set here e_pairs for failed calculations just ridiculously high so that it is sorted out in sortouthighe
      call calcreactionenergy(env, npairs, fragdirs, e_pairs)

      npairs0 = npairs
      ! fragdirs is overwritten in sortout routine
      fragdirs_in = fragdirs
      deallocate (fragdirs)
      ! sort out high energy fragments ccording to reaction energy
      call sortouthighe(env, env%tslevel, e_pairs, env%bhess, npairs0, npairs, fragdirs_in, fragdirs, 3.0_wp) ! here higher scaling 3.0?? ! TODO critical parameter
      deallocate (fragdirs_in)

      if (npairs .eq. 0) then
         write (*, *) "No fragments left"
         nfrags = 0
         return
      end if

      ! We cannot optimize end structures as fragments drift too far apart ...
      ! TODO maybe put at the beginning of new fragmentation as topochanges could produce duplicates here
      ! call optend(env,npairs,fragdirs)

      ! get barrier
      if (.not. env%nots) then
         npairs0 = npairs
         fragdirs_in = fragdirs
         deallocate (fragdirs)
         call tssearch(env, fname, npairs0, npairs, fragdirs_in, fragdirs)
         deallocate (fragdirs_in)
      end if

      if (npairs .eq. 0) then
         write (*, *) "No fragments left"
         nfrags = 0
         return
      end if

      write (*, *) "npairs are", npairs
      ! update number of pairs
      call wrshort_int('npairs2', npairs) ! TODO FIXME make nicer

      ! save current path for later routines
      call getcwd(env%path)
      if (npairs .gt. 0) then
         call montecarlo(env, npairs, fragdirs, eiee, piee, nfragl)
         ! TODO maybe refinement of highes barriers here
         !call check_impfrags(env, npairs, fragdirs, impfrags, imp_fragdirs)
         ! write (*, "('impfrags are: ', *(a, ', '))") (trim(imp_fragdirs(i, 1)), i = 1, impfrags)
         call collectfrags(env, npairs, fragdirs, "fragments", startdir, fraglist, nfrags)
      else
         write (*, *) "No fragments were generated"
         nfrags = 0
         if (env%printlevel .eq. 3) call printpwd
      end if

      if (nfragl .gt. 1) then
         ! collect all fragments in allfragments file
         call appendto('fragments', '../allfragments')
         ! move everything back in main DIR
         do i = 1, npairs
            if (index(fragdirs(i, 3), 'p') .ne. 0) then
               jj = 3 ! if fragmentpair
            else
               jj = 1 ! if isomer
            end if
            do j = 1, jj
               write (dir2, '(a)') trim(startdir)//trim(fragdirs(i, j))
               write (jobcall, '(a)') 'mv '//trim(fragdirs(i, j))//' '//trim(env%startdir)//'/'//trim(dir2)  !//' 2>/dev/null'
               call execute_command_line(trim(jobcall), exitstat=io)
            end do
         end do
      else
         call copy('fragments', 'allfragments')
      end if

   end subroutine calcfragments

! generate fragments with CREST --msreact mode
   subroutine genfragments(env, fname, npairs, nfragl)
      implicit none
      character(len=1024) :: jobcall, dir
      character(len=80) :: fname
      character(len=4) :: fraglevel ! gfn1 or gfn2
      integer :: io, i, ich
      integer :: nshifts, nat, nfrags
      integer, intent(out) :: npairs
      integer, intent(in) :: nfragl ! fragmentation level
      integer :: uhf
      integer :: nat0 ! number of atoms of input molecule
      real(wp) :: sumreac
      real(wp) :: dethr, ewin ! dethr is de threshold ewin is energy window for msreact
      type(runtypedata) :: env
      logical :: ex
      logical :: ldum

      call rdshort_int(trim(env%startdir)//'/fragment.xyz', nat0)

      dethr = env%ieeatm*nat0*3.0_wp*evtoau*autokcal ! take three times of average starting IEE in kcall

      call rdshort_real("sumreac_"//trim(env%tslevel), sumreac)

      if (sumreac .lt. 0.0_wp) then
         if (env%printlevel .eq. 3) write (*, *) "sumreac is negative, set it to 0"
         sumreac = 0.0_wp
      end if
      ewin = dethr - (sumreac*evtoau*autokcal) ! in kcal!!!

      !TODO FIXME maybe use sumreac for this later
      ! save original xyz file
      call copy(fname, 'infrag.xyz')
      call rdshort_int(fname, nat)

      !TODO maybe add atom shifts here? can produce to many fragments at the moment
      !   nshifts = 3 * natoms
      !  if (nshifts .lt. 200)  nshifts = 200
      nshifts = env%msnshifts

      ! crest ignores .UHF in newest version, so we have to set it here
      call rdshort_int('.UHF', uhf)
      ! save .UHF and .CHRG which will be deleted by crest
      call move(".UHF", "uhftemp")
      call move(".CHRG", "chrgtemp")
      write (*, '(a,f10.2,a)') "Generate new fragments with CREST within energy window of", ewin, " kcal/mol"
      if (trim(env%geolevel) /= "gfn2" .and. trim(env%geolevel) /= "gfn1") then
         write (*, *) "using gfn2 instead of"//trim(env%geolevel)//" for CREST msreact" ! crest ignores simply unknown command flags and uses gfn2 as default ^^
         fraglevel = 'gfn2'
      else
         fraglevel = trim(env%geolevel) ! only gfn1 here possible beside gfn2
      end if

      ! each xtb calculation is done with 1 thread, so we have to set it here
      call setomptoone
      write (jobcall, '(a,i0)') 'crest infrag.xyz --msreact --mslargeprint --msnbonds ', env%msnbonds
      write (jobcall, '(a,i0,a,i0)') trim(jobcall)//' --chrg ', env%chrg, ' --uhf ', uhf
      write (jobcall, '(a,i0,a)') trim(jobcall)//' --T ', env%cores, ' --'//fraglevel
      write (jobcall, '(a,f10.5)') trim(jobcall)//' --ewin ', ewin

      !to prevent too many fragments isomers can be skipped
      if (env%msnoiso .or. (nfragl .gt. 1 .and. .not. env%msfulliso)) then
         write (jobcall, '(a)') trim(jobcall)//' --msnoiso'
      end if
      if (env%msiso) then
         write (jobcall, '(a)') trim(jobcall)//' --msiso'
      end if
      write (jobcall, '(a,i0)') trim(jobcall)//' --msnshifts ', nshifts
      write (jobcall, '(a,i0)') trim(jobcall)//' --msnshifts2 ', env%msnshifts2
      write (jobcall, '(a)') trim(jobcall)//' --msattrh'
      ! msmolbar for sorting out, need molbar in path
      if (env%msmolbar) then
         write (jobcall, '(a)') trim(jobcall)//' --msmolbar'
      end if
      ! note really tested yet
      if (env%msinchi) then
         write (jobcall, '(a)') trim(jobcall)//' --msinchi '
      end if
      ! fragment distance for better transition state search
      if (env%msfragdist .gt. 0.0_wp) then
         call wrshort_real('crestms.inp', env%msfragdist)
         open (newunit=ich, file='crestms.inp')
         write (ich, '(a,1x,f8.6)') 'fragdist', env%msfragdist
         close (ich)
         write (jobcall, '(a)') trim(jobcall)//' --msinput crestms.inp '
      end if
      write (jobcall, '(a)') trim(jobcall)//' > msreact.out 2>cresterror.out '
      write (*, *) "crestcall is: ", trim(jobcall)
      call execute_command_line(trim(jobcall), exitstat=io)

      ! count generated fragment pairs
      call rdshort_int('npairs', npairs)
      write (*, *) npairs, " new fragment pairs 'pairs.xyz' and 'isomers.xyz' were generated and written to pXX directories"
      write (*, *) npairs, " fragment structures fragment.xyz were written to pXXfX directories"

      ! restore .UHF and .CHRG deleted by crest
      call move("uhftemp", ".UHF")
      call move("chrgtemp", ".CHRG")
      ! Cleanup
      call rmrf('MSDIR')

    
      if (env%printlevel .lt. 2) then
         call remove('crest_msreact_products.xyz')
         call remove('infrag.xyz')
         call remove('isomers.xyz')
      end if
   end subroutine genfragments

! read fragment list from msreact.out file from CREST -msreact mode
   subroutine readfraglist(env, fname, npairs, nfragl, nfrags, fragdirs, etots)
      integer, intent(inout) :: npairs ! number of fragment pairs or isomers
      integer, intent(in) :: nfragl ! fragmentation level
      integer, intent(out) :: nfrags ! number of fragments
      character(len=80) :: tmp, dum
      character(len=80) :: fragkind
      character(len=80), allocatable, intent(out) :: fragdirs(:, :) ! contains list of fragment dirs 1 pair, 2 + 3 fragments
      character(len=80), allocatable :: fragdirs_in(:, :)
      character(len=80) :: fname
      real(wp), allocatable :: etots(:)
      type(runtypedata) :: env
      integer :: ich, io, i, nfrag, j
      integer :: nat
      integer :: npairs_out ! number of read fragmentpairs and isomers
      integer :: nisomers ! number of isomers
      logical :: ex

      allocate (fragdirs(npairs, 3), etots(npairs))

      fragdirs = '' ! empty array as default
      nfrags = 0
      nisomers = 0

      call minigrep('msreact.out', 'CREST terminated normally.', ex)

      if (ex) then
         open (newunit=ich, file='msreact.out', status="old", action="read")
         do while (.true.)
            read (ich, '(a)', iostat=io) tmp
            if (io < 0) exit
            if (index(tmp, 'directory | fragment type | rel. energy [kcal/mol]') .ne. 0) then
               if (.not. env%msiso) read (ich, '(a)', iostat=io) tmp

               do i = 1, npairs
                  read (ich, *) fragdirs(i, 1), fragkind, etots(i)
                  if (fragkind .eq. "isomer") then
                     nisomers = nisomers + 1
                     if (env%msnoiso .or. (nfragl .gt. 1 .and. .not. env%msfulliso)) then
                        fragdirs(i, 1) = ''
                        cycle ! skip isomers
                     end if
                  else
                     do j = 1, 2
                        nfrags = nfrags + 1
                        read (ich, *) fragdirs(i, j + 1), dum, dum
                     end do
                  end if
               end do
               if (index(tmp, '===================================================') .ne. 0) exit
            end if

         end do
         close (ich)
      else
         write (*, *) "fragment generation failed for fragment ", trim(env%path), " check msreact.out"
         stop "fragment generation failed, check msreact.out"
      end if

      write (*, *) "Generated fragments were read from msreact.out"
      if (env%msnoiso .or. (nfragl .gt. 1 .and. .not. env%msfulliso)) then
         write (*, *) "Isomers were skipped"
         npairs_out = npairs - nisomers
         write (*, *) "npairs_out is", npairs_out, npairs, nisomers
         fragdirs_in = fragdirs
         deallocate (fragdirs)
         call sortout0elements(npairs, npairs_out, fragdirs_in, fragdirs, .false.)
         npairs = npairs_out
      end if

      write (*, *) "Number of generated fragments and isomers is", nfrags, nisomers
   end subroutine readfraglist

   subroutine cleanup_crestmsreact
      implicit none

     

      call remove("cregen.out.tmp")
      call remove("cre_members")
      call remove("crest_allproducts.xyz")
      call remove("crest_best.xyz")
      call remove("cresterror.out")
      call remove("crest_input_copy.xyz")
      call remove("crestms.inp")
      call remove("crest.restart")
      call remove("crest_unique_products.xyz")
      call remove("crest_unique_products.sorted")
      call remove("fragmentpairs.xyz")
      call remove("coord")
      call remove("scoord.1")
      call remove("coordprot.0")
      call remove("lmocent.coord")

      

   end subroutine cleanup_crestmsreact

end module fragmentation
