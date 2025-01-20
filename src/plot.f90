module plot
   use qcxms2_data
   use iomod
   use utility
   use isotope_pattern
   implicit none
contains
! collect fragments after each fragmentationn in list
   subroutine readfragments(env, npairs, indexlist, nfrags)
      implicit none
      character(len=80) :: dir
      type(runtypedata) :: env
      character(len=80), allocatable, intent(out) :: indexlist(:)
      integer, intent(in) :: npairs
      integer, intent(out) :: nfrags
      integer :: i, j
      logical :: ex
      real(wp), allocatable :: exp_mass(:), exp_int(:)

      if (npairs .ne. 0) then
         nfrags = 0
         do i = 1, npairs
            do j = 1, 2
               write (dir, '(a,i0,a,i0)') trim(env%index)//"p", i, "f", j
               inquire (file=trim(dir)//"/fragment.xyz", exist=ex)
               if (ex) then
                  nfrags = nfrags + 1
               end if
            end do
         end do

         allocate (indexlist(nfrags))
         write (*, *) "nfrags1", nfrags, trim(dir)
         nfrags = 0
         do i = 1, npairs
            do j = 1, 2
               write (dir, '(a,i0,a,i0)') trim(env%index)//"p", i, "f", j
               inquire (file=trim(dir)//"/fragment.xyz", exist=ex)
               if (ex) then
                  nfrags = nfrags + 1
                  indexlist(nfrags) = dir
                  ! write(*,*) "trimdir is:", trim(dir)
                  ! write(*,*) "indexlist is:", trim(indexlist(nfrags))
               end if
            end do
         end do
         write (*, *) "nfrags2", nfrags, trim(dir)

      end if
   end subroutine readfragments

!subroutine to collect fragments after each fragmentation in file "fragments"
! directory (p1f1,p1f2p4f2), mass, pfrag,
   subroutine collectfrags(env, npairs, fragdirs, fname, startdir, fraglist, nfrags)
      implicit none
      character(len=80) :: dir, dir2, sumform
      character(len=*) :: fname
      character(len=80) :: startdir, product
      character(len=80), intent(in)   :: fragdirs(:, :)
      character(len=80), allocatable, intent(out) :: fraglist(:)
      character(len=80) :: query1, query2, query3
      integer :: i, j, k
      integer :: ich
      integer :: nfrags
      integer :: chrg, uhf
      type(runtypedata) :: env
      integer, intent(in) :: npairs
      real(wp) :: pfrag, mass, sumreac, irc, barrier, de
      logical :: ex, ldum
      logical :: restart

      restart = .false.
      inquire (file=trim(fname), exist=ex) ! check if file exists
      if (ex) restart = .true. ! if file exists, set restart to true

! first prepare fraglist
      nfrags = 0
      do i = 1, npairs
         if (index(fragdirs(i, 3), 'p') .ne. 0) then
            nfrags = nfrags + 2
         else
            nfrags = nfrags + 1
         end if
      end do

      allocate (fraglist(nfrags))
      j = 1

      do i = 1, npairs
         if (index(fragdirs(i, 3), 'p') .ne. 0) then
            fraglist(j) = trim(startdir)//trim(fragdirs(i, 2))
            j = j + 1
            fraglist(j) = trim(startdir)//trim(fragdirs(i, 3))
            j = j + 1
         else
            fraglist(j) = trim(startdir)//trim(fragdirs(i, 1))
            j = j + 1
         end if
      end do

      write (*, *) "Writing fragments and isomers to file ", trim(fname)
      write (*, *) "Important fragments are: "
      write (*, *) "dir | mass | sumformula | rel. I / %"
      if (npairs .ne. 0) then
         open (newunit=ich, file=fname)
         write (ich, *) "Dir  fragment_type sumreac  de  barrier  irc"
         do i = 1, npairs
            call rdshort_real(trim(fragdirs(i, 1))//"/barrier_"//trim(env%tslevel), barrier)
            call rdshort_real(trim(fragdirs(i, 1))//"/sumreac_"//trim(env%tslevel), sumreac)
            call rdshort_real(trim(fragdirs(i, 1))//"/de_"//trim(env%tslevel), de)
            call rdshort_real(trim(fragdirs(i, 1))//"/ircmode_"//trim(env%geolevel), irc)
            if (index(fragdirs(i, 3), 'p') .ne. 0) then
               write (product, *) trim(startdir)//trim(fragdirs(i, 1))//" fragmentpair"
               write (ich, '(a,2x,f6.3,2x,f6.3,2x,f6.3,2x,f10.3)') trim(product), sumreac, de, barrier, irc
               do j = 2, 3
                  call rdshort_real(trim(fragdirs(i, j))//"/mass", mass)
                  call rdshort_real(trim(fragdirs(i, j))//"/pfrag", pfrag)
                  call getsumform(trim(fragdirs(i, j))//"/fragment.xyz", sumform)
                  write (product, *) trim(startdir)//trim(fragdirs(i, j))
                  write (ich, '(a,2x,f9.3,2x,a,2x,e10.4)') trim(product), mass, trim(sumform), pfrag
                  if (pfrag .gt. env%pthr) then
                     write (*, '(a,2x,f9.3,2x,a,2x,f10.1)') trim(product), mass, trim(sumform), pfrag/env%nsamples*100
                  end if
               end do
            else
               write (product, *) trim(startdir)//trim(fragdirs(i, 1))//" isomer"
               write (ich, '(a,2x,f6.3,2x,f6.3,2x,f6.3,2x,f10.3)') trim(product), sumreac, de, barrier, irc
               call rdshort_real("mass", mass) ! has to be mass of starting structure
               call rdshort_real(trim(fragdirs(i, 1))//"/pfrag", pfrag)
               call getsumform(trim(fragdirs(i, 1))//"/isomer.xyz", sumform)
               write (product, *) trim(startdir)//trim(fragdirs(i, 1))
               write (ich, '(a,2x,f9.3,2x,a,2x,e10.4)') trim(product), mass, trim(sumform), pfrag
               if (pfrag .gt. env%pthr) then
                  write (*, '(a30,2x,f9.3,2x,a,2x,f10.1)') trim(product), mass, trim(sumform), pfrag/env%nsamples*100
               end if
            end if
         end do
         close (ich)
      end if
   end subroutine

   ! normalization of all obtained fragment peaks
   ! to get final spectrum
   subroutine getpeaks(env, nfrags, allfrags)
      implicit none
      integer :: i, j, npeaks, npeaks1, maxatm
      integer :: index_mass, count_mass
      integer, intent(in) :: nfrags
      character(len=80), allocatable :: allfrags(:)
      character(len=80) :: fname
      integer :: ich, ich2
      type(runtypedata) :: env
      real(wp) :: pfrag, mass, sumreac, sumea, imax
      real(wp) :: list_masses(100000)   !TODO for large mols (over 100) too small?
      real(wp) :: intensity(100000)
      real(wp), allocatable ::  added_masses(:), added_ints(:)
      real(wp), allocatable :: peak_ints(:), peak_masses(:), peak_ints1(:), peak_masses1(:)
      real(wp), allocatable :: isotope_masses(:), exact_intensity(:)
      logical :: ex1, ex2, ex
      logical, allocatable :: double(:)
      integer, allocatable :: ind(:) ! keep track of origin after sorting with index
      !real(wp) :: score, jaccard ! STEIN and SCOTT and JACCARD index score
      integer, parameter :: nrnd = 50000
      real(wp), allocatable :: rnd(:, :)
      real(wp) :: r

! for each fragment structure
      real(wp), allocatable :: fragment_masses(:), fragment_intensities(:)

! to be sure to start in correct directory if something was messed up in the beginning
      call chdir(trim(env%startdir))
      write (*, *)
      write (*, *) "-------------------------------------"
      write (*, *) "|Reaction network generation finished|"
      write (*, *) "-------------------------------------"
      write (*, *)

      write (*, *) "Normalize Intensities for all fragments"
      write (*, *) "Number of fragments is: ", nfrags
! set whole array to zero to avoid later problems in normalization off allpeaks.dat
      intensity = 0.0_wp
      fname = 'fragment.xyz'
      call rdshort_int(trim(fname), maxatm)

      ! initialize the random number array (efficiency)
      allocate (rnd(50000, maxatm))
      do i = 1, 50000
         do j = 1, maxatm ! maximum number of atoms of a fragment is nat of starting fragment
            call random_number(r)
            rnd(i, j) = r
         end do
      end do

      index_mass = 0
      count_mass = 0

      !first peak
      !> compute isotope patterns via monte carlo and save intensites
      !  of all possible combinations
      call getmass(env, fname, maxatm, rnd, index_mass, exact_intensity, isotope_masses)
      !> sort the single fragment intensities over the entire list of frags
      if (index_mass > 0) then
         call check_entries(index_mass, isotope_masses, exact_intensity, &
                            list_masses, intensity, count_mass)
      end if

      deallocate (exact_intensity)
      deallocate (isotope_masses)

      do i = 1, nfrags
         inquire (file=trim(allfrags(i + 1))//"/fragment.xyz", exist=ex1)
         if (.not. ex1) then
            inquire (file=trim(allfrags(i + 1))//"/isomer.xyz", exist=ex2)
            if (ex2) then
               fname = 'isomer.xyz'
               ex = .true.
            end if
         else
            fname = 'fragment.xyz'
            ex = .true.
         end if

         if (ex) then
            call chdir(trim(allfrags(i + 1)))
            inquire (file='fragment.xyz', exist=ex)
            if (ex) then
               fname = 'fragment.xyz'
            else
               fname = 'isomer.xyz'
            end if
            call getmass(env, fname, maxatm, rnd, index_mass, exact_intensity, isotope_masses)
            if (index_mass > 0) then
               call check_entries(index_mass, isotope_masses, exact_intensity, &
                                  list_masses, intensity, count_mass)
            end if
            deallocate (exact_intensity)
            deallocate (isotope_masses)
            call chdir(env%startdir)
         else
            write (*, *) "WARNING: no fragment or isomer found in", trim(allfrags(i + 1))
         end if
      end do

      if (env%int_masses) then
      do i = 1, count_mass
         list_masses(i) = nint(list_masses(i))

      end do
      end if

! add peaks with same m/z values
! in principle we only need this, if we round to integers
      allocate (added_masses(count_mass), added_ints(count_mass), double(count_mass))
      double = .false.
      npeaks = 0
      do i = 1, count_mass
         if (double(i)) then
            cycle
         else
            ! sort out masses below mass threshold
            ! arent visible in spectrum and distort intensities
            if (list_masses(i) .ge. env%mthr) then
               npeaks = npeaks + 1
               added_masses(npeaks) = list_masses(i)
               added_ints(npeaks) = intensity(i)
            end if
         end if
         do j = i + 1, count_mass
            if (.not. double(j)) then
               if (abs((list_masses(i) - list_masses(j))) .lt. 0.000001_wp) then ! TODO FIXME set mass spec tolerance here??
                  ! already in isotope routine?
                  if (list_masses(i) .ge. env%mthr) then
                     added_ints(npeaks) = added_ints(npeaks) + intensity(j)
                  end if
                  double(j) = .true.
               else
                  double(j) = .false.
               end if
            end if
         end do
      end do

! sort masses
      allocate (peak_masses(npeaks), peak_ints(npeaks), ind(npeaks))
      do i = 1, npeaks
         peak_masses(i) = added_masses(i)
         ind(i) = i ! array is not allowed to be empty somehow for qsort
      end do

      call qsort(peak_masses, 1, npeaks, ind)
      do i = 1, npeaks
         peak_ints(i) = added_ints(ind(i))
      end do

! peaks after ithr sorting
      allocate (peak_masses1(npeaks), peak_ints1(npeaks))
! number of peaks after ithr sorting
      npeaks1 = 0
! sort out intensities below certain threshold
      if (env%intthr .gt. 0.0_wp) write (*, *) "Sorting out peaks below intensity threshold of", env%intthr, " %"
      env%intthr = env%intthr*100 ! to 10000 like in NIST
!normalize to 10000 like in NIST

      imax = maxval(added_ints)
      do i = 1, npeaks
         peak_ints(i) = peak_ints(i)/(1.0_wp*imax)*10000
         npeaks1 = npeaks1 + 1
         peak_ints1(npeaks1) = peak_ints(i)
         peak_masses1(npeaks1) = peak_masses(i)
      end do

!write peaks to peaks.dat
      write (*, *)
      write (*, *) "Computed spectrum is written to peaks.dat and peaks.csv"
      write (*, *) "Important peaks (above 10 %)  are: "
      write (*, *) "(m/z | intensity - normalized to 10000)"
      open (newunit=ich, file='peaks.dat')
      open (newunit=ich2, file='peaks.csv')
      do i = 1, npeaks1
         if (peak_ints(i) .gt. env%intthr) then ! sort out here !TODO not for matchng score though ...
            write (ich, *) peak_masses1(i), peak_ints1(i)
            write (ich2, *) peak_masses1(i), ",", peak_ints1(i)
         end if
         !important fragments
         if (peak_ints1(i)/100 .gt. 10) then
            write (*, '(f10.6,2x,f10.1)') peak_masses1(i), peak_ints1(i)
         end if
      end do
      close (ich)
      close (ich2)

! read experimental spectrum given as exp.dat as csv file and compute matchscore (jaccard and old Stein and Scott score)
!inquire(file='exp.dat',exist=ex)
!if (ex) then
!call matchscore(env,npeaks1, peak_masses1,peak_ints1,score, jaccard)
!
!end if

! TODO rewrite in cleaner way
! write all peaks without isotope pattern to "allpeaks.dat", so that we know which fragment contributes
      env%noiso = .true. ! all peaks without isotope pattern here, so that we know which fragment contributes to intensity of m/z
      allocate (fragment_masses(nfrags + 1), fragment_intensities(nfrags + 1))

      index_mass = 0
      count_mass = 0

      !first peak
      !> compute isotope patterns via monte carlo and save intensites
      !  of all possible combinations
      call getmass(env, env%infile, maxatm, rnd, index_mass, exact_intensity, isotope_masses)
      ! no check entries here, we want intensity of each fragment
      fragment_masses(1) = isotope_masses(1)
      fragment_intensities(1) = exact_intensity(1)
      deallocate (exact_intensity)
      deallocate (isotope_masses)

      do i = 1, nfrags

         inquire (file=trim(allfrags(i + 1))//"/fragment.xyz", exist=ex1)
         if (.not. ex1) then
            inquire (file=trim(allfrags(i + 1))//"/isomer.xyz", exist=ex2)
            if (ex2) then
               fname = 'isomer.xyz'
               ex = .true.
               if (env%noplotisos) cycle
            end if
         else
            fname = 'fragment.xyz'
            ex = .true.
         end if
         if (ex) then
            if (env%printlevel .eq. 3) write (*, *) "GOING TO FRAGMENT", I, trim(allfrags(i + 1))
            call chdir(trim(allfrags(i + 1)))
            inquire (file='fragment.xyz', exist=ex)
            if (ex) then
               fname = 'fragment.xyz'
            else
               fname = 'isomer.xyz'
            end if
            call getmass(env, fname, maxatm, rnd, index_mass, exact_intensity, isotope_masses)
            fragment_masses(i + 1) = isotope_masses(1)
            fragment_intensities(i + 1) = exact_intensity(1)
            deallocate (exact_intensity)
            deallocate (isotope_masses)
            call chdir("..")
         else
            write (*, *) "WARNING: no fragment or isomer found in", trim(allfrags(i + 1))
         end if
      end do

!write all peaks without isotope pattern to "allpeaks.dat", so that we know which fragment contributes
      imax = maxval(fragment_intensities)
      allfrags(1) = "input structure"
      write (*, *)
      write (*, *) "Writing all fragment intensities without isotope pattern to allpeaks.dat"
      write (*, *) "Important fragments (above 10 %)  are: "
      write (*, *) "(fragment | m/z | intensity - normalized to 10000)"
      open (newunit=ich, file='allpeaks.dat')
      write (ich, *) "fragment|", "fragment mass|", "relative intensity"
      do i = 1, nfrags + 1
         fragment_intensities(i) = fragment_intensities(i)/(1.0_wp*imax)*10000
         write (ich, '(x,a,2x,f10.6,2x,f10.1)') trim(allfrags(i)), fragment_masses(i), fragment_intensities(i)
!important fragments
         if (fragment_intensities(i)/100 .gt. 10) then
            write (*, '(x,a30,2x,f10.6,2x,f10.1)') trim(allfrags(i)), fragment_masses(i), fragment_intensities(i)
         end if
      end do
      close (ich)

   end subroutine getpeaks

   ! currently not used
   subroutine matchscore(env, npeaks, peak_masses, peak_ints, score, jaccard)
      use xtb_mctc_accuracy, only: wp
      implicit none
      type(runtypedata) :: env
      real(wp) :: peak_masses(:), peak_ints(:)
      real(wp), allocatable :: exp_mass0(:), exp_int0(:), exp_mass(:), exp_int(:)
      integer ::  exp_entries
      real(wp) :: score, jaccard
      integer, intent(in) :: npeaks
      integer :: nexp ! number of experimental peaks
      integer :: i, j
      integer :: pp !peak pair number pp
      integer :: numcalc ! number of calculated peaks

      integer :: cnt

      real(wp) :: tmax
      real(wp), allocatable :: w_exp(:) !weighted vectors
      real(wp), allocatable :: w_thr(:) !weighted vectors
      real(wp) :: sum1, sum2, sum3, sum4, dot, fr
      real(wp) :: norm
      real(wp), allocatable :: p(:, :)
      real(wp) :: added_masses(npeaks)
      real(wp) :: added_intensities(npeaks)

      integer :: ich, iocheck, k
      integer :: imax, mmax
      character(len=80) :: line
      real(wp) :: xx(100)
      integer :: nn
! only csv file possible
      open (newunit=ich, file='exp.dat', status='old')
      exp_entries = 0

      !>> count the entries
      do
         read (ich, '(a)', iostat=iocheck)!line
         if (iocheck < 0) exit

         exp_entries = exp_entries + 1
      end do
      rewind (ich)
      !>> allocate memory
      allocate (exp_mass0(exp_entries), &
        &       exp_int0(exp_entries), &
        & exp_mass(exp_entries), &
        &      exp_int(exp_entries))

      iocheck = 0

      !>> count the intensities and masses
      do i = 1, exp_entries
         read (ich, '(a)', iostat=iocheck) line
         if (iocheck < 0) exit

         do k = 1, 80
            if (line(k:k) == ',') then
               line(k:k) = ' '
            end if
         end do

         call readl(line, xx, nn)

         exp_mass0(i) = xx(1)
         exp_int0(i) = xx(2)
      end do
      close (ich)

      imax = maxval(exp_mass0)
      !imin = minval(exp_mass)

      mmax = maxval(exp_int0)

      ! norm intensities to 10000
      do i = 1, exp_entries
         exp_int0(i) = (exp_int0(i)/mmax)*10000.0_wp

      end do

      ! sort out masses and intensities below threshold
      ! if m/z highest peak is sorted out, we get here bullshit
      nexp = 0
      do i = 1, exp_entries
         if (exp_int0(i) .ge. env%intthr .and. exp_mass0(i) .ge. env%mthr) then
            nexp = nexp + 1
            exp_int(nexp) = exp_int0(i)
            exp_mass(nexp) = exp_mass0(i)

         end if
      end do

      pp = 0
      !> get weighting of experimental spectrum
      allocate (w_exp(nexp))
      norm = 0.0_wp
      do i = 1, nexp
!!DEL  write(*,*) exp_mass(i)**1.0_wp , exp_int(i)**1.0_wp
         w_exp(i) = exp_mass(i)**1.0_wp*exp_int(i)**1.0_wp
         !  w_exp(i) = exp_mass(i)**3.0_wp * exp_int(i)**0.6_wp
         norm = norm + w_exp(i)**2
      end do

      norm = sqrt(norm)

      do i = 1, nexp
         w_exp(i) = w_exp(i)/norm
      end do

      !> get weighting of calculated spectrum
      allocate (w_thr(npeaks))
      norm = 0.0_wp
      do j = 1, npeaks
         w_thr(j) = (peak_masses(j))**1.0_wp*(peak_ints(j))**1.0_wp
         ! w_thr(j) = (peak_masses(j))**3.0_wp * (peak_ints(j))**0.6_wp
         norm = norm + w_thr(j)**2
      end do

      norm = sqrt(norm)

      do j = 1, npeaks
         w_thr(j) = w_thr(j)/norm
      end do

      dot = 0.0_wp
      sum1 = 0.0_wp
      sum2 = 0.0_wp
      sum3 = 0.0_wp

      do i = 1, exp_entries
         do j = 1, npeaks
            if (exp_mass(i) == peak_masses(j)) then
               pp = pp + 1
               sum1 = sum1 + w_exp(i)*w_thr(j)
               sum2 = sum2 + (w_exp(i))**2
               sum3 = sum3 + (w_thr(j))**2
            end if
         end do
      end do

      dot = sum1**2/(sum2*sum3)
      !write(*,*) 'DOT', dot

      deallocate (w_exp, w_thr)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! masses and intensity scalement for the second term
      ! m**0
      ! int**1
      allocate (w_exp(exp_entries), w_thr(npeaks))

      norm = 0.0_wp

      do i = 1, exp_entries
         !norm = norm + w_exp(i)**2 ! + exp_int(i) !w(1,i)**2
         norm = norm + exp_int(i)**2 !w(1,i)**2
      end do

      norm = sqrt(norm)

      do i = 1, exp_entries
         w_exp(i) = exp_int(i)/norm
      end do

  !!!!!!!!!!!!!!!!!!!!
      norm = 0.0_wp
      do j = 1, npeaks
         norm = norm + peak_ints(j)**2 !w(2,i)**2
      end do

      norm = sqrt(norm)

      do j = 1, npeaks
         w_thr(j) = peak_ints(j)/norm
      end do

      !pp   = 0
      !p    = 0.0_wp
      sum4 = 0.0_wp
      fr = 0.0_wp
      cnt = 0

      allocate (p(2, pp))

      do i = 1, exp_entries
         do j = 1, npeaks
            if (exp_mass(i) == peak_masses(j)) then
               cnt = cnt + 1
               p(1, cnt) = w_exp(i)
               p(2, cnt) = w_thr(j)
            end if
         end do
      end do

      ! ardous loop
      call calcfr(pp, p, sum4)
      if (pp == npeaks) sum4 = sum4 + 1.0_wp
      !write(*,*) 'SUM4', sum4
      !write(*,*) 'numcalc', numcalc
      !write(*,*) 'numcalc', npeaks

      fr = sum4/float(pp)
      score = (npeaks*dot + pp*fr)/(npeaks + pp)
      write (*, *) "Fr is:", fr
      write (*, *) "Dot is:", dot
      write (*, *) "npeaks", npeaks
      write (*, *) "PP is:", pp

      call getjaccard(npeaks, peak_masses, exp_entries, exp_mass, jaccard)

      write (*, *) "Composite match score, see "
      write (*, *) "Stein, S. E.; Scott, D. R. J. Am. Soc. Mass Spectrom. 1994, 5, 859-866"
      write (*, *) "Score is", score
      write (*, *) "Jaccard index is", jaccard
      write (*, *) "Product of Score and Jaccard index is", jaccard*score
      return

   end subroutine matchscore

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computes the second expression in the mass spec match score
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine calcfr(pp, pair, sum4)
      use xtb_mctc_accuracy, only: wp
      implicit none

      integer :: i
      integer :: pp

      real(wp) :: pair(2, pp)
      real(wp) :: sum4

      sum4 = 0.0_wp

      do i = 2, pp
         if (abs((pair(1, i)*pair(2, i - 1))/(pair(1, i - 1)*pair(2, i))) < 1.0_wp) then

            sum4 = sum4 + (pair(1, i)*pair(2, i - 1))/(pair(1, i - 1)*pair(2, i))

         elseif (abs((pair(1, i)*pair(2, i - 1))/(pair(1, i - 1)*pair(2, i))) > 1.0d0) then

            sum4 = sum4 + 1/((pair(1, i)*pair(2, i - 1))/(pair(1, i - 1)*pair(2, i)))

         elseif (abs((pair(1, i)*pair(2, i - 1))/(pair(1, i - 1)*pair(2, i))) == 1.0d0) then
            sum4 = sum4 + 1.0_wp
         end if
      end do

   end subroutine calcfr

! get jaccard index matching score between experimental and calculated spectrum
   subroutine getjaccard(npeaks, peak_masses, exp_entries, exp_mass, jaccard)
      use xtb_mctc_accuracy, only: wp
      implicit none

      real(wp) :: peak_masses(:), exp_mass(:)
      integer :: npeaks, exp_entries
      real(wp), intent(out) :: jaccard
      integer ::  sum_intersec, sum_union
      integer :: i, j, l

      write (*, *) "npeaks is", npeaks
      write (*, *) "exppeaks is", exp_entries
      sum_intersec = 0
      l = 0
      do i = 1, npeaks
         do j = 1, exp_entries
         if (peak_masses(i) .eq. exp_mass(j)) then
            sum_intersec = sum_intersec + 1
         end if
         end do
      end do
      sum_union = 0
      i = 1
      j = 1
      do while (i .le. npeaks .and. j .le. exp_entries)
      if (peak_masses(i) .lt. exp_mass(j)) then
         sum_union = sum_union + 1
         i = i + 1

      elseif (exp_mass(j) .lt. peak_masses(i)) then
         sum_union = sum_union + 1
         j = j + 1

      else
         sum_union = sum_union + 1
         i = i + 1
         j = j + 1
      end if
      end do
      do while (i .le. npeaks)
         sum_union = sum_union + 1
         i = i + 1
      end do
      do while (j .le. exp_entries)
         sum_union = sum_union + 1
         j = j + 1
      end do
      write (*, *) "union is", sum_union
      write (*, *) "intersection is", sum_intersec
      jaccard = float(sum_intersec)/float(sum_union)
   end subroutine getjaccard

   subroutine check_entries(index_mass, isotope_masses, exact_intensity, &
                            list_masses, intensity, count_mass)
      integer :: loop, loop2, index_mass
      integer :: count_mass
      integer :: sum_index, i
      integer :: check
      real(wp) :: xmass
      real(wp) :: isotope_masses(index_mass)
      real(wp) :: exact_intensity(index_mass)
      real(wp) :: chrg

      real(wp) :: list_masses(100000)
      real(wp) :: intensity(100000)

      real(wp) :: mass_diff
      logical :: there = .true.

      chrg = 1.0_wp
      xmass = 0
      loop = 0
      loop2 = 0
      check = 0
      !count_mass = count_mass + index_mass
      if (count_mass == 0) then
         sum_index = index_mass
      else
         sum_index = count_mass + index_mass
      end if

      !> loop over the list from the isotope subroutine
      outer: do loop2 = 1, index_mass
         loop = 0

         !> loop over all entries of new (check) list
         inner: do
            loop = loop + 1
            mass_diff = abs(list_masses(loop) - isotope_masses(loop2))

            if (list_masses(loop) == 0.0_wp) then
               there = .false.
               !write(*,*) 'NULL'
               exit inner

               !elseif  ( list_masses(loop) == isotope_masses(loop2) ) then
            elseif (mass_diff < 1.0d0 - 10 .or. mass_diff == 0.0_wp) then
               there = .true.
               intensity(loop) = intensity(loop) + 1*exact_intensity(loop2) &
                                 *abs(chrg)
               if (loop2 < index_mass) cycle outer
               if (loop2 == index_mass) exit inner

               !>> false if not in list, store
               !elseif ( list_masses(loop) /= isotope_masses(loop2) ) then
            elseif (mass_diff > 1.0d0 - 10) then
               there = .false.
               if (loop == sum_index) exit inner
            end if
         end do inner

         if (.not. there) then
            count_mass = count_mass + 1
            list_masses(count_mass) = isotope_masses(loop2)
            intensity(count_mass) = intensity(count_mass) + 1*exact_intensity(loop2) &
                                    *abs(chrg) ! TODO FIXME check if this is correct ??
         end if

      end do outer
   end subroutine check_entries
end module
