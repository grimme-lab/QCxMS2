!> module to handle the external quantum mechanical program codes
module qmmod
   use iomod
   use qcxms2_data
   use utility
   implicit none

contains
!> basic procedure:
!> prepare a qm calculation -> steps
!> 1. copy the input file to the level of theory directory
!> 2. prepare the jobcall
!> 3. execute the jobcall
!> 4. read the output file and write the energy to a file
!> 5. clean up

   subroutine prepqm(env, fname, level, job, jobcall, fout, pattern, cleanupcall, there, restart, chrg, uhf)
      implicit none
      type(runtypedata) :: env
      character(len=1024), intent(out) :: jobcall
      character(len=1024) :: qmdelcall
      character(len=80), intent(out) :: fout, pattern ! fout and pattern for readout
      character(len=1024), intent(out) :: cleanupcall
      character(len=*), intent(in) :: job ! job = tddft (only ORCA at the moment), sp, opt, hess, (ohess, currently not used) bhess (only xtb) optts (only orca)
      character(len=*), intent(in) :: level ! level of theory
      integer, intent(in), optional :: chrg, uhf
      logical, intent(in), optional  :: restart
      logical, intent(out) :: there ! value already there?
      logical :: re
      character(len=*), intent(in)      :: fname

      character(len=80) :: sumform
      character(len=80) :: query
      character(len=80) :: fqm
      logical :: ex, ldum
      real(wp) :: edum
      integer :: io

      integer :: chrg1
      integer :: spin, uhf1

      if (present(restart)) then
         re = restart
      else
         re = .false.
      end if
      if (re) then
         write (*, *) "calculation failed, retrying with different settings"
      end if

      ! if no input charge is given, trust .CHRG file
      ! if not there just use the default

      if (present(chrg)) then
         chrg1 = chrg
         call wrshort_int('.CHRG', chrg1)
         call remove('.UHF') ! always remove .UHF file if charge changes
      else
         inquire (file='.CHRG', exist=ex)
         if (ex) then
            call rdshort_int('.CHRG', chrg1)
         else
            chrg1 = env%chrg
            call wrshort_int('.CHRG', chrg1)
            call remove('.UHF') ! always remove .UHF file if charge changes
         end if
      end if

      if (present(uhf)) then
         uhf1 = uhf
      else
         inquire (file='.UHF', exist=ex)
         if (ex) then
            call rdshort_int('.UHF', uhf1)
         else
            uhf1 = 0
            call printspin(env, fname, spin, chrg1)
            ! uhf in .UHF formalism of xtb is number of unpaired electrons
            uhf1 = spin - 1
         end if
         ! correct for H+
         if (uhf1 == -2) uhf1 = 0
      end if
      call wrshort_int('.UHF', uhf1)

      fqm = 'qmdata'
      inquire (file=fqm, exist=ex)
      if (ex) then
         write (query, '(a,1x,a,1x,i0,1x,i0)') trim(level), trim(job), chrg1, uhf1
         call grepval(fqm, trim(query), ldum, edum)
      else
         ldum = .false.
      end if
      there = .false.
      if (ldum .and. abs(edum) .gt. 0.00001_wp) then
         there = .true.
      end if

      ! special cases
      ! bhess only in xtb
      if (job == 'bhess') then
         ! bhess always gfn2 or gfn2spinpol
         if (env%geolevel == "gfn2spinpol") then 
            call prepxtb(env, fname, 'gfn2spinpol', job, jobcall, fout, pattern, cleanupcall, re)
         else
            call prepxtb(env, fname, 'gfn2', job, jobcall, fout, pattern, cleanupcall, re)
         end if
         !no TS optimization in xtb currently available
      elseif ((job == 'optts' .and. level == 'gfn2') &
         & .or. (job == 'optts' .and. level == 'gfn2spinpol') &
         & .or. (job == 'optts' .and. level == 'gfn1') &
         & .or. (job == 'optts' .and. level == 'gxtb')) then
         call preporca(env, fname, level, job, jobcall, fout, pattern, cleanupcall, re)
      elseif (job == 'tddft') then
         call preporca(env, fname, level, job, jobcall, fout, pattern, cleanupcall, re)
      elseif (job == "hess") then ! ORCA freqs are strange ....
         ! we go here always for gfn2 frequencies in orca
         if (env%geolevel == "gfn2spinpol") then 
            call preporca(env, fname, 'gfn2spinpol', job, jobcall, fout, pattern, cleanupcall, re)
         else
            call preporca(env, fname, 'gfn2', job, jobcall, fout, pattern, cleanupcall, re)
         end if
      else

         ! check for program
         select case (level)
            ! xtb calculations TODO move gfn2 and gfn1 to tblite calculation
         case ('gfn2','gfn2spinpol', 'gfn1', 'pm6', 'dxtb', 'gff')
            call prepxtb(env, fname, level, job, jobcall, fout, pattern, cleanupcall, re)
            !      DFT, set orca or turbomole and read template file orca_temp.inp and cefinecall?
         case ('pbe', 'b973c', 'r2scan3c', 'pbeh3c' &
            & , 'wb97x3c', 'pbe0', 'ccsdt', 'kpr2scan50d4', 'pbe0tzvpd', 'wb97xd4tz', 'wb97xd4qz' &
            & , 'wb97xd4matzvp') ! could be extended to other functionals
            call preporca(env, fname, level, job, jobcall, fout, pattern, cleanupcall, re)
            ! special methods
         case ('gxtb')
            if (job .eq. "sp" .or. job .eq. "opt") then
               call prepgxtb(env, fname, level, job, jobcall, fout, pattern, cleanupcall,re)
            end if
            ! special QM methods have to be included here
         case default
            write (*, *) "Warning! Keyword ", level, " is not supported!!!"
            write (*, *) "Please check the manual for supported keywords for QM levels"
            STOP
         end select
      end if
      if (env%printlevel .eq. 3) cleanupcall = 'echo "keepdir" >/dev/null' ! do not remove files if lprint is true

      ! if we optimize a structure and it is not already there, we have to remove the wrong singlepoint energy
      fqm = 'qmdata'
      inquire (file=fqm, exist=ex)
      if (job == 'opt' .and. ex .and. .not. there) then
         write (query, '(a,1x,i0,1x,i0)') "sp", chrg1, uhf1
         if (env%exstates .gt. 0) then
            write (query, '(a,1x,i0,1x,i0)') "tddft", chrg1, uhf1
         end if

         write (qmdelcall, '(a)') 'sed -i "/'//trim(query)//'/d" qmdata' ! leads to bug, TODO rewrite sed command with fortran routine ...
         call execute_command_line(trim(qmdelcall), exitstat=io)
         ! Also the optimization entry is not up do date anymore
         write (qmdelcall, '(a)') 'sed -i "/opt/d" qmdata'
         call execute_command_line(trim(qmdelcall), exitstat=io)
      end if

   end subroutine prepqm

! readout qm calculation and check if calculation was succesfull
   subroutine readoutqm(env, fname, level, job, fout, pattern, value, failed)

      implicit none
      type(runtypedata) :: env
      character(len=80), intent(in) :: fout, pattern
      character(len=*), intent(in) :: level, job
      character(len=80) :: sumform
      character(len=80) :: fqm, query
      character(len=162) :: qmdatentry
      character(len=*)      :: fname ! xyz structure file
      real(wp), intent(out) :: value
      real(wp) :: bwenergie, edisp, etot, escf ! for tddft
      integer :: chrg, uhf, ich, ierr
      logical ::   ex, ldum
      logical, intent(out) :: failed
      failed = .false.

      fqm = 'qmdata'

      call rdshort_int('.CHRG', chrg)
      call rdshort_int('.UHF', uhf)
      if (job == 'optts') then
         call minigrep(trim(fout), trim(pattern), ldum)
         if (.not. ldum) then
            write (*, *) "Warning! ", trim(pattern), " not found in ", trim(fout), "calculations failed in:"
            call printpwd
            failed = .true.
         else
            failed = .false.
            call move('opt.xyz', trim(fname))
         end if
         return
      end if

      write (query, '(a,1x,a,1x,i0,1x,i0)') trim(level), trim(job), chrg, uhf

      call grepval(fqm, trim(query), ldum, value) ! check if value already there
      if (ldum .and. abs(value) .gt. 0.000001_wp) then
         failed = .false.
         return
      end if

      ! checkfor H-atoms as QM codes sometimes don't like to calculate molecules without electrons
      call getsumform(trim(fname), sumform)
      if (sumform == "H1" .and. chrg == 1 .and. level .ne. 'gfn1' .and. level .ne. 'gfn2' .and. level .ne. 'gfn2spinpol' .and. level .ne. 'gxtb') then
         failed = .false.
         value = 0.0_wp
      else
         if (level .ne. 'gxtb' .or. job .ne. 'sp') then
            call grepval(trim(fout), trim(pattern), ldum, value)
            if (.not. ldum) then
               write (*, *) "Warning! ", trim(pattern), " not found in ", trim(fout), "calculations failed in:"
               call printpwd
               failed = .true.
            end if
         end if
      end if

      if (job == 'tddft') then
         call grepval('orca.out', 'Dispersion correction', ldum, edisp)
         call grepval('orca.out', 'E(SCF)  =', ldum, escf)
         etot = escf + edisp ! in au
         call readoutstates(env, bwenergie, failed) ! bwenergie in au, contribution of higher states to total energy
         value = etot + bwenergie
      end if

      ! write to qmdata
      if (level == 'gxtb' .and. job .eq. 'sp') then
         call readoutgxtb(ldum, value)
         !DEL write(*,*) "energy is ", value
         if (.not. ldum) then
            write (*, *) "Warning no energy found in energy file from gxtb"
            call printpwd
            failed = .true.
         else
            failed = .false.
         end if
      end if

      write (qmdatentry, *) trim(query)//"  ", value ! two spaces for easier reading
      inquire (file=fqm, exist=ex, iostat=ierr)
      if (ex) then !  This is done so that always the newest value is first found by the grepval routine
         call copy(fqm, 'qmdatatemp')
         call remove('qmdata')
         open (newunit=ich, file=fqm, status="new", action="write")
         write (ich, *) trim(qmdatentry)
         call appendto('qmdatatemp', fqm)
         call remove('qmdatatemp')
      else
         open (newunit=ich, file=fqm, status="new", action="write")
         write (ich, *) trim(qmdatentry)
      end if
      close (ich)
      ! write for comparison sp energy to qmdata
      if (job == 'tddft') then
         write (query, '(a,1x,a,1x,i0,1x,i0)') trim(level), "sp", chrg, uhf
         write (qmdatentry, *) trim(query)//"  ", etot
         call copy(fqm, 'qmdatatemp')
         call remove('qmdata')
         open (newunit=ich, file=fqm, status="new", action="write")
         write (ich, *) trim(qmdatentry)
         call appendto('qmdatatemp', fqm)
         call remove('qmdatatemp')
         close (ich)
      end if

      if (.not. failed) then
         if (job == 'optts' .or. job == 'opt' .or. job == 'ohess') then
            call move('opt.xyz', trim(fname))
         end if
      end if

   end subroutine readoutqm

   subroutine readoutgxtb(ldum, value)
      implicit none
      real(wp) :: value
      logical :: ldum
      character(len=100) :: line
      integer :: i, ich, io
      real :: temp

      ldum = .false.
      value = 0.0_wp

      open (newunit=ich, file='energy', status='old', action='read', iostat=io)

      if (io .ne. 0) then
         write (*, *) 'File energy not found or cannot be opened.'
         return
      end if

      ! Read lines up to the second line
      do i = 1, 2
         read (ich, '(a)', iostat=io) line
         if (io .ne. 0) then
            write (*, *) 'Error reading line ', i, ' or file is too short.'
            close (ich)
            return
         end if
      end do

      ! Try to read the second value in the second line
      read (line, *, iostat=io) temp, value
      if (io == 0) then
         ldum = .true.
      else
         ldum = .false.
      end if

      close (ich)
   end subroutine readoutgxtb

!DEL !TODO DELETEME
! stupid routine
! 'qmdata' saves all quantum mechanical data
! runtype specific data is stored in 'rundata', this will be deleted and not needed for restarts
! write entry to qmdata file
   subroutine writetoqm(level, job, value, chrg_in, uhf_in)
      implicit none

      character(len=*), intent(in) :: level, job
      character(len=80) :: fname, query
      real(wp), intent(in) :: value
      integer, intent(in), optional :: chrg_in, uhf_in
      integer :: chrg, uhf
      integer :: ich, ierr
      logical ::  ldum, ex

      if (.not. present(chrg_in)) then
         inquire (file='.CHRG', exist=ex)
         if (ex) then
            call rdshort_int('.CHRG', chrg)
         else
            error stop
         end if
      else
         chrg = chrg_in
      end if
      if (.not. present(uhf_in)) then
         inquire (file='.UHF', exist=ex)
         if (ex) then
            call rdshort_int('.UHF', uhf)
         else
            error stop
         end if
      else
         uhf = uhf_in
      end if

      fname = 'qmdata'

      write (query, '(a,1x,a,1x,i0,1x,i0)') trim(level), trim(job), chrg, uhf

      !write to qmdata
      inquire (file=fname, exist=ex, iostat=ierr)
      if (ex) then
         open (newunit=ich, file=fname, status="old", position="append", action="write")
      else
         open (newunit=ich, file=fname, status="new", action="write")
      end if
      ! check if value already there
      !call grepval('qmdata',trim(query),ex,edum)
      !if (ex) then ! remove entry of qmdata

      write (ich, '(a,2x,f10.5)') trim(query), value
      close (ich)

   end subroutine writetoqm
! read data from qmfile
   subroutine grepqmdata(level, item, ldum, value, chrg_in, uhf_in, dir)
      implicit none
      character(len=*), intent(in) :: level
      character(len=*), intent(in) :: item
      character(len=*), intent(in), optional :: dir
      integer, intent(in), optional :: chrg_in, uhf_in
      integer :: chrg, uhf
      character(len=1024) :: dir1
      character(len=80) :: query
      logical :: ldum, ex
      real(wp), intent(out) :: value

      if (present(dir)) then
         dir1 = dir
      else
         dir1 = '.'
      end if

      if (.not. present(chrg_in)) then
         inquire (file='.CHRG', exist=ex)
         if (ex) then
            call rdshort_int('.CHRG', chrg)
         else
            error stop
         end if
      else
         chrg = chrg_in
      end if
      if (.not. present(uhf_in)) then
         inquire (file='.UHF', exist=ex)
         if (ex) then
            call rdshort_int('.UHF', uhf)
         else
            error stop
         end if
      else
         uhf = uhf_in
      end if

      write (query, '(a,1x,i0,1x,i0)') trim(level)//" "//trim(item), chrg, uhf
      ! write(*,*) "query is ",trim(query)
      call grepval(trim(dir1)//'/qmdata', trim(query), ldum, value)

   end subroutine grepqmdata

! write entry to qmdata file
   subroutine writetoqmsub(subdir, level, job, value)
      implicit none

      character(len=*), intent(in) :: level, job
      character(len=*), intent(in) :: subdir
      character(len=80) :: fname, query
      real(wp), intent(in) :: value
      integer :: chrg, uhf, ich, ierr
      logical ::  ldum, ex

      call rdshort_int('.CHRG', chrg)
      call rdshort_int('.UHF', uhf)
      fname = trim(subdir)//'/qmdata'

      write (query, '(a,1x,a,1x,i0,1x,i0)') trim(level), trim(job), chrg, uhf

      !write to qmdata
      inquire (file=fname, exist=ex, iostat=ierr)
      if (ex) then
         open (newunit=ich, file=fname, status="old", position="append", action="write")
      else
         open (newunit=ich, file=fname, status="new", action="write")
      end if
      write (ich, '(a,2x,f10.5)') trim(query), value
      close (ich)

   end subroutine writetoqmsub
! check for entry in qmdata file
! grepvalue should be enough?

!if geometry changed, this should be deleted??
!subroutine delqmentry(env,job,level,value)

!end subroutine delqmentry
!

! prepare an ORCA calculation
   subroutine preporca(env, fname, level, job, jobcall, fout, pattern, cleanupcall, restart)
      use structools
      implicit none
      type(runtypedata) :: env
      character(len=*), intent(in)      :: fname
      character(len=*), intent(in)      :: job
      character(len=*), intent(in)      :: level
      character(len=1024), intent(out) :: jobcall
      character(len=80), intent(out) :: fout, pattern ! fout and pattern for readout
      character(len=1024), intent(out) :: cleanupcall ! remove uneccessary files
      character(len=200) ::  jobkeyword
      character(len=80) :: levelkeyword
      character(len=80) :: sumform
      character(len=80) :: xtbstring
      character(len=80) :: act_atom_string

      integer :: chrg
      integer :: uhf
      integer :: mult
      integer :: etemp
      logical :: ex
      integer :: io, ich
      integer :: nstates ! for TD-DFT
      integer :: nmode ! for TS optimization
      logical, intent(in) :: restart
      logical :: fermi ! apply fermi smearing

      if (env%fermi) then
         fermi = .true.
      else
         fermi = .false.
      end if
      if (restart) then
         fermi = .true.
      end if

      if (restart) then
         call execute_command_line("rm orca_*", exitstat=io)
      end if

      nstates = env%exstates ! TODO make this variable or tune this

      ! set electronic temperature for SCF etemp = 5000 + ax/100*20000 ;ax is percentage of Fock exchange
      select case (level)
      case ('pbeh3c')
         levelkeyword = 'PBEh-3c'
         etemp = 13400
      case ('r2scan3c')
         levelkeyword = 'R2SCAN-3c'
         etemp = 5000
      case ('b973c')
         levelkeyword = 'B97-3c'
         etemp = 5000
      case ('wb97x3c')
         levelkeyword = 'wB97X-3c'
         etemp = 15000 ! TODO tune this value
      case ('kpr2scan50d4')
         levelkeyword = ' kPr2SCAN50 RIJCOSX def2-QZVPP/C def2/J D4'
         etemp = 15000 ! TODO tune this value
      case ('wb97xd4tz')
         levelkeyword = 'wB97X-D4 def2-TZVP'
         etemp = 15000 ! TODO tune this value
      case('wb97xd4matzvp') ! 
         levelkeyword = 'wB97X-D4 ma-def2-TZVP'
         etemp = 15000 ! TODO tune this value
      case ('wb97xd4qz')
         levelkeyword = 'wB97X-D4 def2-QZVP'
         etemp = 15000 ! TODO tune this value
      case ('pbe0')
         levelkeyword = 'PBE0 def2-TZVP'
         etemp = 10000
      case ('pbe0tzvpd')
         levelkeyword = 'PBE0 def2-TZVP'
         etemp = 10000
      case ('pw6b95d4tzvpd')
         levelkeyword = 'PW6B95-D4 def2-TZVP'
         etemp = 10600
      case ('pbe')
         levelkeyword = 'PBE SV(P)'
         etemp = 5000
      case ('gfn2')
         levelkeyword = 'XTB2'
         etemp = 300
         fermi = .false. ! activated by default anyway
         if (restart) then
            etemp = 5000
            fermi = .true.
         end if
      case ('gfn2spinpol')
         levelkeyword = 'XTB2'
         etemp = 300
         fermi = .false. ! activated by default anyway
         if (restart) then
            etemp = 5000
            fermi = .true.
         end if   
      case ('gxtb')
         levelkeyword = 'XTB2'
         fermi = .false.
         if (restart) then
            etemp = 15000
            fermi = .true.
         end if
      case ('ccsdt')
         levelkeyword = 'cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/C TightSCF PModel'
         fermi = .false.
      case ('gfn1')
         levelkeyword = 'XTB1'
         etemp = 300
         fermi = .false. ! activated by default anyway
         if (restart .and. nstates .eq. 0) then
            etemp = 5000
            fermi = .true.
         end if
        !!!!!!!!!
      case default
         levelkeyword = trim(level)
         etemp = 5000
      end select
!TODO transform input to orca keyword?
      select case (job)
      case ('sp')
         jobkeyword = 'SP'
      case ('opt')
         jobkeyword = 'LOOSEOPT'
      case ('hess')
         jobkeyword = 'FREQ'
         fermi = .false. ! ORCA hessians do not work with fermi smearing
         !jobkeyword = 'NUMFREQ' ! too expensive
      case ('ohess')
         jobkeyword = 'OPT FREQ' ! not sure if this works
         fermi = .false. ! ORCA hessians do not work with fermi smearing
      case ('optts')
         
         jobkeyword = 'OptTS LOOSEOPT' !
         if (env%tsoptgmf) jobkeyword = 'OptTS(GMF) LOOSEOPT' !
      case ('tddft')
         jobkeyword = 'SP' !
      end select

      call rdshort_int('.CHRG', chrg)
      call rdshort_int(".UHF", uhf)
!is this right ?
      mult = uhf + 1
      ! TODO check here for orca template file

      OPEN (newunit=ich, FILE=' orca.inp')
      if (level .ne. "ccsdt") write (ich, *) "! DEFGRID2 "
      ! end if

      if (restart) then ! TODO tune this settings, maybe combine with FERMI smearing????
         write (ich, *) "! SlowConv"
         !if (job == 'optts' .or. job == 'opt') write (ich, *) "! LOOSEOPT"
      end if

      !if (level .ne. 'wb97x3c')
      write (ich, *) "! "//trim(levelkeyword)

      if (env%solv) then 
         if ( level =="gfn1" .or. level == 'gfn2' .or. level == 'gfn2spinpol') then 
            write (ich, *) "! ALPB(water)" 
         else
          write (ich, *) "! CPCM(WATER)" 
         end if
      end if

      write (ich, *) "! "//trim(jobkeyword)
             !use  tblite
      if (level == 'gfn2spinpol') then
         xtbstring = 'XTBINPUTSTRING2 "--tblite --spinpol"'
         write (ich, *) "%xtb"
         write (ich, '(a)') trim(xtbstring)
         write (ich, *) "end"
      end if
      if (job == 'optts') then
         ! for now we do it dirty like that, dont want to rewrite the whole routine ..
         call rdshort_int('nmode', nmode)
         write (ich, *) "%geom"
         write (ich, '(a,1x,i0,a)') "TS_Mode {M ", nmode - 1, "}" !orca counts from 0
         write (ich, *) " end"
         write (ich, *) ' inhess read  inhessname "orca.hess"'
         write (ich, *) " maxiter 100"
         if (env%tsoptact) then
            call get_active_cnstring(env, act_atom_string) 
            write (ich, *)  "TS_Active_Atoms { "//trim(act_atom_string)//" }" ! atoms that are involved in TS, e.g. for proton
         write (ich, *)  "end"
         write(ich,*)    "TS_Active_Atoms_Factor 1.5" ! factor by which the cutoff for bonds is increased for
         end if
         write (ich, *) "end"
         if (env%geolevel == 'gxtb') then
            xtbstring = 'XTBINPUTSTRING2 "--driver ''gxtb -c orca.xtbdriver.xyz -symthr 0.0 ''"'
            if (fermi)  xtbstring = 'XTBINPUTSTRING2 "--driver ''gxtb -c orca.xtbdriver.xyz -tel 15000 -symthr 0.0''"'
            write (ich, *) "%xtb"
            write (ich, '(a)') trim(xtbstring)
            write (ich, *) "end"
            call touch('.GRAD')
         end if
      end if

      ! TODO problem of xtbdriver here,dirty first solution for now
      write (ich, *) "%pal"
      write (ich, '(a,i0)') "nprocs ", env%threads
      write (ich, *) "end"

      write (ich, *) "%maxcore 8000" ! TODO make parameter or read in orca sample input file

      if (fermi) then
         write (ich, *) "%scf SmearTemp", etemp, "end"
      end if
      write (ich, *) "*xyzfile ", chrg, " ", mult, " "//trim(fname)

      ! for correct dissociation behavior
      ! for correct dissociation behavior, especially for CID/closed-shell systems

      if (job == 'tddft') write (ich, *) "%TDDFT  NROOTS ", nstates, " END"
      if (level .ne. "gfn2" .and. level .ne. "gfn1" .and. level .ne. "gxtb" .and. level .ne. "ccsdt") then
         if (trim(fname) == "ts.xyz") write (ich, *) "! UKS"
      end if

      ! first try to converge TZVP and then with diffuse functions added
      if (level == 'pw6b95tzvpd' .or. level == 'pbe0tzvpd') then
         write (ich, *) "$new_job"
         write (ich, *) "!def2-TZVPD"
         write (ich, *) "*xyzfile 0 1 in.xyz"
         write (ich, *) "*xyzfile ", chrg, " ", mult, " "//trim(fname)
      end if

      close (ich)

      ! ORCA has to be called with full path name
      write (jobcall, '(a)') '$(which orca) orca.inp > orca.out 2>/dev/null'

      fout = 'orca.out'

      write (cleanupcall, '(a)') "rm orca.gbw orca.prop orca.inp orca.bibtex orca.property.txt orca.engrad orca.xtbrestart"
      select case (job)
      case ('sp')
         pattern = 'FINAL SINGLE POINT ENERGY'
      case ('opt')
         !pattern = 'THE OPTIMIZATION HAS CONVERGED' ! TODO have to rewrite this for optimizations to check if they converged, just ignore it for now
         pattern = 'FINAL SINGLE POINT ENERGY'
         write (jobcall, '(a)') trim(jobcall)//' && cp orca.xyz opt.xyz'
         write (cleanupcall, '(a)') trim(cleanupcall)//" orca.opt orca.gu.tmp"
      case ('hess')
         if (env%notemp) then
            pattern = "Zero point energy                ..."
         else
            pattern = "Total correction"
         end if
         ! convert to g98 format
         call copy(trim(fname), 'xtbin.xyz')
         !TODO new xtb version
         ! inquire (file='g98.out', exist=ex)
         !if (.not. ex) then write (jobcall, '(a)') trim(jobcall)//' && xtb xtbin.xyz --hess  > xtbhess.out 2>/dev/null'
         write (jobcall, '(a)') trim(jobcall)//' && xtb thermo xtbin.xyz --orca orca.hess  > thermo.out 2>/dev/null'
         ! write (jobcall, '(a)') trim(jobcall)//' && xtb xtbin.xyz --hess  > xtbhess.out 2>/dev/null'
         ! write(cleanupcall,'(a)') trim(cleanupcall)//" charges xtbhess.out xtbrestart xtbtopo.mol xtbhess.xyz wbo"
      case ('ohess')
         if (env%notemp) then
            pattern = "zero point energy"
         else
            pattern = "G(RRHO) contrib."
         end if
         write (jobcall, '(a)') trim(jobcall)//' && cp orca.xyz opt.xyz'
      case ('optts')
         pattern = 'THE OPTIMIZATION HAS CONVERGED' !TODO have to rewrite this for optimizations to check if they converged, just ignore it for now
         !pattern = 'FINAL SINGLE POINT ENERGY'
         write (jobcall, '(a)') trim(jobcall)//' && cp orca.xyz opt.xyz' ! orca.xyz is last geometry of optimization
         write (jobcall, '(a)') trim(jobcall)//' && cp orca.out geo.out' ! orca.xyz is last geometry of optimization
         write (cleanupcall, '(a)') trim(cleanupcall)//" orca.opt orca.gu.tmp"
      case ('tddft') ! just for checking
         pattern = 'FINAL SINGLE POINT ENERGY'
      end select
      write (cleanupcall, '(a)') trim(cleanupcall)//" >/dev/null 2>/dev/null"

   end subroutine preporca

!prepare a xtb calculation via systemcall of xtb !TODO implement tblite here
! can use xtb as driver for other programms here
   subroutine prepxtb(env, fname, level, job, jobcall, fout, pattern, cleanupcall, restart)
      implicit none
      type(runtypedata) :: env
      character(len=*), intent(in)      :: fname
      character(len=*), intent(in)      :: level
      character(len=*), intent(in)      :: job
      character(len=1024), intent(out) :: jobcall
      character(len=80), intent(out) :: fout, pattern ! fout and pattern for readout
      character(len=1024), intent(out) :: cleanupcall
      character(len=200) ::  jobkeyword
      character(len=80) :: levelkeyword
      integer :: chrg
      integer :: etemp
      real(wp) :: temp ! temperature used for RRHO calculation
      logical :: ex
      logical :: fermi ! apply fermi smearing
      logical, intent(in) :: restart
      integer :: io, ich

      call rdshort_int('.CHRG', chrg)

      etemp = 300.0_wp
      if (restart) etemp = 5000.0_wp ! to force convergence
! set temperature higher for bhess
      if (.not. env%notemp .and. job == 'bhess') then
         temp = env%temp !TODO DEL * env%thermoscale
         OPEN (newunit=ich, FILE='xtb.inp')
         write (ich, *) "$thermo"
         write (ich, *) "  temp=", temp
         write (ich, *) "  sthr=", env%sthr
         write (ich, *) "end"
         close (ich)
      end if

      !!!>for TESTING
      !if (level .eq. "dxtb") call copydxtbparam(env)
      !if (level .eq. "pm6") call wrpm6input()
      call copy(trim(fname), 'xtbin.xyz')
      write (jobkeyword, '(a)') 'xtb xtbin.xyz --'//trim(job)//' --'//trim(level)

      ! use tblite for correct uhf
      if (level == "gfn2spinpol") write (jobkeyword, '(a)') trim(jobkeyword)//' --tblite --spinpol'
      ! for gff optimizations, to prevent H rearranges back to initial position
      if (level == "gff" .and. job == 'opt') write (jobkeyword, '(a)') 'xtb xtbin.xyz --gff --opt vtight'
      !if (env%fermi) ! TODO check if this is really needed
      ! apply here high fermi smearing for xTB, this needs to be robust and converge everytime ?
      write (jobkeyword, '(a,i4)') trim(jobkeyword)//' --etemp ', etemp

      if (env%solv) write (jobkeyword, '(a)') trim(jobkeyword)//' --alpb water '

       !!!>for TESTING
      !if (env%dxtb) write (jobkeyword, '(a)') trim(jobkeyword)//' --vparam dxtb_param.txt'
      !if (level .eq. "pm6") write (jobkeyword, '(a)') trim(jobkeyword)//' --input pm6.inp --mopac'

      ! TODO fix xtbdriver for hessian calculations
      ! if (job == 'bhess' .or. job == 'ohess') write (jobkeyword, '(a)') trim(jobkeyword)//' --input xtb.inp'
      !  if (job == 'bhess' .and. env%geolevel == 'gxtb') then
      !      write (jobkeyword, '(a)') trim(jobkeyword)//' --driver "gxtb -c xtbdriver.xyz"'
      !     call touch('.GRAD')
      !  end if

     !!!>FOR DFTBHESS
      !  if (job == 'bhess' .and. env%geolevel == 'wb97x3c') then
      !     write(jobcall,'(a)') 'x2t '//trim(fname)//' > coord 2>/dev/null'
      !     call execute_command_line(trim(jobcall), exitstat=io)
      !     write (jobkeyword, '(a)') 'xtb coord --bhess --driver "cefine --func wb97x-3c --sym c1 && ridft && rdgrad"'
      !     !inquire (file='basis', exist=ex) cefine --func wb97x-3c --sym c1 --noopt
      !     !call touch('.GRAD')
      !  end if

      write (jobkeyword, '(a)') trim(jobkeyword)//' > xtb.out 2>/dev/null'
      write (jobcall, '(a)') trim(jobkeyword)
      if (job == 'ohess' .or. job == 'opt') write (jobcall, '(a)') trim(jobcall)//' && cp xtbopt.xyz opt.xyz'
      cleanupcall = 'rm xtbrestart  xtb.out wbo charges'
! for readout and cleanup
      fout = 'xtb.out'
      select case (job)
      case ('sp')
         pattern = '| TOTAL ENERGY'
      case ('opt')
         pattern = '| TOTAL ENERGY'
         write (cleanupcall, '(a)') trim(cleanupcall)//' xtbopt.log'
      case ('bhess', 'hess', 'ohess')
         if (env%notemp) then
            pattern = "zero point energy"
         else
            pattern = "G(RRHO) contrib."
         end if
         write (cleanupcall, '(a)') trim(cleanupcall)//' xtbopt.log xtbhess.xyz'
         !cleanupcall='rm xtbrestart  xtb.out xtbopt.log g98.out vibspectrum wbo hessian'
      end select
      if (level .eq. "gff") write (cleanupcall, '(a)') trim(cleanupcall)//' gfnff_charges gfnff_topo'
      write (cleanupcall, '(a)') trim(cleanupcall)//' xtbtopo.mol .xtboptok'
      write (cleanupcall, '(a)') trim(cleanupcall)//'  >/dev/null 2>/dev/null'
   end subroutine prepxtb

!  employ g-xTB, preliminary version
   subroutine prepgxtb(env, fname, level, job, jobcall, fout, pattern, cleanupcall,restart)
      implicit none
      type(runtypedata) :: env
      character(len=*), intent(in)      :: fname
      character(len=*), intent(in)      :: level
      character(len=*), intent(in)      :: job
      character(len=1024), intent(out) :: jobcall
      character(len=80), intent(out) :: fout, pattern ! fout and pattern for readout
      character(len=1024), intent(out) :: cleanupcall
      character(len=80) :: levelkeyword
      integer :: chrg
      integer :: etemp
      logical :: ex
      integer :: io
      logical :: fermi ! apply fermi smearing
      logical, intent(in) :: restart


      if (restart) fermi = .true.

      call copy(trim(fname), 'xtbin.xyz')

      call remove('gxtbrestart')
      call remove('ceh.charges')

      cleanupcall = 'rm gxtbrestart energy ceh.charges gxtb.out '
      select case (job)
      case ('sp')

         call remove('.GRAD')
         call remove('.HESS')

         if (env%fermi .or. fermi) then ! apply fermi smearing
            write (jobcall, '(a)') 'gxtb -c xtbin.xyz -tel 15000 > gxtb.out 2>errorfile' ! TODO has currently no effect
         else
            write (jobcall, '(a)') 'gxtb -c xtbin.xyz > gxtb.out 2>errorfile'
         end if
         fout = 'gxtb.out'
         !pattern='total                   :'
         pattern = 'total                         '
      case ('opt')
         call touch('.GRAD')
         call remove('coord') ! necessary the way the xtb driver works at the momen TODO fix this
         write (jobcall, '(a)') 'xtb xtbin.xyz --opt --driver "gxtb -c xtbdriver.xyz -symthr 0.0" > opt.out 2>/dev/null'
         if (env%fermi .or. fermi) write (jobcall, '(a)') 'xtb xtbin.xyz --opt --driver "gxtb -c xtbdriver.xyz -tel 15000 -symthr 0.0" > opt.out 2>/dev/null'
         write (jobcall, '(a)') trim(jobcall)//' && cp xtbopt.xyz opt.xyz'
         fout = 'opt.out'
         write (cleanupcall, '(a)') trim(cleanupcall)//' xtbdriver.xyz opt.out .GRAD'
         pattern = "| TOTAL ENERGY"
      case ('hess')
         call touch('.HESS')

         call execute_command_line(trim(jobcall), exitstat=io)
         write (jobcall, '(a)') 'gxtb -c xtbin.xyz > gxtb.out 2>/dev/null'
         if (env%fermi .or. fermi)  write (jobcall, '(a)') 'gxtb -c xtbin.xyz tel 15000 > gxtb.out 2>/dev/null'
         fout = 'gxtb.out'
         !pattern='total                    :'
         ! pattern="total                         -115.33940557"
         !pattern="total                 "
         write (cleanupcall, '(a)') trim(cleanupcall)//' .GRAD .HESS'
      case default
         write (*, *) trim(job)//" not supported for gxtb"
         stop
      end select
      if (env%printlevel .lt. 2) write (cleanupcall, '(a)') trim(cleanupcall)//' errorfile'
      write (cleanupcall, '(a)') trim(cleanupcall)//' >/dev/null 2>/dev/null'

   end subroutine prepgxtb

!TODO FIXME
! H-atoms make problems for QM codes, so we have to remove them
! can also save time here, since we dont have to compute them
! only for sp and opt, bhess is ignored by xtb and fine
! possible charges are 0,1 and later possibly -1 ??
   subroutine calc_h_atoms(env, fname, level, job, npairs_in, fragdirs_in, npairs_out, fragdirs_out, chrg_in, uhf_in)
      implicit none
      type(runtypedata) :: env
      character(len=*), intent(in) :: job ! job = sp, opt, hess, (ohess, currently not used) bhess (only xtb) optts (only orca or tm) in xtb nomenclature
      character(len=80), intent(in) :: level ! level of theory
      integer, intent(in), optional :: chrg_in, uhf_in
      character(len=*)      :: fname
!logical, allocatable, intent(out)   :: hatomlist(:,:)
      character(len=80) :: fragdirs_in(:, :)
      character(len=80), allocatable, intent(out)   :: fragdirs_out(:, :)
      integer :: nhatoms
      integer, intent(in) :: npairs_in
      integer, intent(out) :: npairs_out
      character(len=80) :: sumform
      logical :: ex
      logical :: ishatom
      real(wp) :: edum
      integer :: io
      integer :: i, j

      integer :: chrg
      integer :: spin, uhf

!hatomlist = .false.
      nhatoms = 0

!allocate(hatomlist(npairs,3))
      do i = 1, npairs_in
         !check if  fragment  is isomer ! only fragments can be H-atom
         !  if (index(fragdirs(i,3),'p') .ne. 0) then
         do j = 2, 3
            ! call getsumform(trim(fragdirs(i,j))//"/fragment.xyz",sumform)
            call getsumform(trim(fname), sumform)
            if (sumform == "H1") then
               !hatomlist(i,j) = .true.
               ishatom = .true.
               nhatoms = nhatoms + 1
               npairs_out = npairs_in - nhatoms
               ! etot is zero for H-atom for DFT but not for SQM methods!!!!
               if (chrg == 1) then
                  if (level .eq. 'r2scan3c' .or. level .eq. 'pbeh3c' &
                  & .or. level .eq. 'wb97x3c' .or. level .eq. 'wb97xd4tz' &
                  & .or. level .eq. 'wb97xd4qz' .or. level .eq. 'kpr2scan50d4') then ! TODO make this more general if dft ..
                     call writetoqm(level, 'sp', 0.0_wp, chrg, uhf)
                     ! TODO we could add here total energies for the H-atom for common DFAs
                     !elseif(chrg1 == 1)
                  end if
               end if
            end if
         end do
         ! end if
      end do

   end subroutine calc_h_atoms

! only for ORCA currently
   subroutine readoutstates(env, bwenergie, failed)
      use xtb_mctc_constants
      use xtb_mctc_convert
      implicit none
      type(runtypedata) :: env
      real(wp), allocatable :: energies(:)
      real(wp), intent(out)    :: bwenergie
      integer :: i, j, io, ich, ichout
      integer :: nstates
      character(512) :: line
      character(80) :: fname
      character(80) :: dumc ! just dummy for reading
      real(wp) :: dumf ! just dummy for reading
      logical :: ex
      logical :: failed

      !nstates = env%ntddftstates
      nstates = 10

      fname = 'orca.out'
      inquire (file=fname, exist=ex)
      if (ex) then
         open (newunit=ich, file='orca.out')
      else
         write (*, *) "orca.out not found, no states read"
         failed = .true.
         call printpwd()
         return
      end if

      allocate (energies(nstates + 1))
      ! TODO CARE about units and absolute and relative energies
      energies = 0.0_wp
      ! call grepval('orca.out','E(SCF)',ex,energies(1)) ! GS energy energies are relative ...

      i = 1
      do
         read (ich, '(a)', iostat=io) line
         if (io /= 0) exit  ! Exit the loop when reaching the end of the file

         if (index(line, 'STATE ') .ne. 0) then
            i = i + 1
            !DEL   write(*,*) trim(line)
            j = index(line, "eV") ! read energy in eV
            !DEL   write(*,*) j
            write (*, *) line(j - 9:j - 1)
            read (line(j - 9:j - 1), *) energies(i)
            !DEL write(*,*) "energy of state", i, "is", energies(i)
            if (energies(i) < 0.0_wp) then
               write (*, *) "Negative energy in orca.out, set to 0, false state was chosen"
               energies(i) = 0.0_wp
            end if
         end if
      end do
      close (ich)
      call boltzmannweights(env, energies, nstates, bwenergie) ! define bwenergie as energie contribution of higher states to GS energy
      ! if (bwenergie > 0.0_wp) then
      !     failed = .false.
      ! end if

      write (*, *) "Boltzmann weighted energy in eV is", bwenergie
      ! transform bwenergie back to au
      bwenergie = bwenergie*evtoau
   end subroutine readoutstates

!!!>for TESTING
! just a dirty hack to test dxtbparam by copying it from startdir in current DIR
!   subroutine copydxtbparam(env)
!      implicit none
!      type(runtypedata) :: env
!      character(len=1024) :: copycall
!      write (copycall, '(a)') "cp "//trim(env%startdir)//"/dxtb_param.txt dxtb_param.txt"
!      call execute_command_line(trim(copycall))
!   end subroutine copydxtbparam
!
!   subroutine wrpm6input()
!      implicit none
!      integer :: io
!
!      open (newunit=io, file='pm6.inp')
!      write (io, *) "$external"
!      write (io, *) "mopac input=1scf pm6-d3h4X DISP aux(42,PRECISION=12,MOS=-99999,COMP) scfcrt=1.d-7 xyz geo-ok grad"
!      write (io, *) "mopac file=mopac"
!      write (io, *) "$end"
!      close (io)
!
!   end subroutine wrpm6input

end module qmmod
