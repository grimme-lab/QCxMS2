module mcsimu
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_convert
   use qcxms_boxmuller
   use qcxms_iee
   use qcxms2_data
   use iomod
   use utility
   use qmmod
   use cid

   implicit none

contains
   ! perform a montecarlo simulation to calculate all intensities from barriers for one fragmentation step
   subroutine montecarlo(env, npairs, fragdirs, eieein, piee, nfragl)
      implicit none
      integer, intent(in) :: npairs
      integer, intent(in) :: nfragl ! number of subsequent fragmentation 1st 2nd 3rd ...
      real(wp), intent(in) :: eieein(:), piee(:)
      real(wp), allocatable :: eiee(:)
      type(runtypedata) :: env
      character(len=80), intent(in)   :: fragdirs(:, :) ! array of all fragment directories
      type(timer):: tim
      integer :: i, j, k, nsamples !nsamples:: number of runs of Montecarlo
      integer :: nat, nvib! number of vibrational degrees of freedom for reaction should always be 3N(atoms)-6 of start fragment
      real(wp), allocatable ::  pfrag(:, :)
      real(wp), allocatable :: barriers(:), de(:), irc(:)! barriers and imaginary modes at TS
      real(wp), allocatable :: barriersre(:) ! barrier of backreaction
      real(wp), allocatable        :: kabsre(:)
      real(wp), allocatable        :: kabs(:)
      real(wp), allocatable        :: k_ey(:) ! absolute rate constant but from eyring equation
      real(wp) :: krel, ktotal, maxpfrag, sumpfrag, kmax !relative and absolute rate constant
      real(wp) :: freq
      real(wp) :: sumpfrag0
      real(wp) :: tf ! time of flight in MS
      real(wp) :: prem ! proportion of fragment which remains stable, determined by ksum and tf
      real(wp) :: checksum, norm, sump ! just checking if all intensities add up to one
      real(wp) :: sumreac, p0 ! sum of all reaction energies to get to initial fragment, current intensity of initial fragment
      real(wp) :: qf1, qf2 ! statistical charge of fragment
      real(wp) :: KER, KER_Eex, KER_Ear! kinetic energy release upon fragmentation, consists of part of excess energy Eex and part of reverse reaction barrier Ear
      real(wp) :: KERav, KERavold, sumpiee
      real(wp) :: de0, Ea0 ! reaction energy and barrier of last reaction for computing KER
      real(wp) :: IEE
      real(wp) :: qmass, mf1, mf2 !
      real(wp) :: pfragtemp ! store temporary probability of fragment
      real(wp) :: Eavg, kavg, dgavg
      real(wp) :: drrho, earrho
      real(wp), allocatable :: tav(:) ! averaged half time of reaction
      real(wp) :: tav0 ! averaged half time of reaction of last reaction
      integer :: maxprob ! integer of most probable energy
      real(wp) :: maxiee
      real(wp) :: time_relax ! deprecated paramet TODO DEL
      real(wp) :: scaleeinthdiss ! scaling factor of IEE for H-dissociation
      integer :: nincr
      real(wp), allocatable :: f1_rrhos(:), f2_rrhos(:), ts_rrhos(:), start_rrhos(:)
      real(wp), allocatable :: dgs(:), dgeas(:), alldgs(:, :), alldgsre(:, :) ! array of free energies of activation for different energies

      real(wp) ::  pfrag_prec ! store temporary probability of fragment after isomerization eq.
      real(wp) :: kmax_frag, kmax_iso

      real(wp), allocatable :: active_cn(:) ! number of active coordination number
      real(wp) :: max_active_cn
      real(wp), parameter :: iee_cn_scale = 1.0_wp ! scaling factor for IEE  in dependence of the active coordination number

      integer :: ind

      ! real(wp), allocatable ::  KERold(:) ! array of old KERenergies ! todo not necessary, too much overhead just take average keravold

      ! for eatomscale
      real(wp) :: ef1, ef2, efold ! factor to scale energy of fragment
      integer :: natf1, natf2 ! number of atoms in generated fragments
      integer :: ich
      integer :: nmode ! number of IRC mode
      integer :: nvib0, nat0 ! number of atoms and vibrational degrees of freedom from last reaction
      ! and energy of reverse reaction barrier
      logical :: ex, ldum
      logical, allocatable :: ishdiss(:) ! is reaction H-dissociation
      logical, allocatable :: isrearr(:) ! is reaction rearrangement
      integer :: ilen, l, ndir
      character(len=80) ::  fname
      character(len=1024) :: pairdir, fname2
      character(len=1024) :: thisdir
      character(len=1024) :: jobcall

      write (*, *)
      write (*, *) "----------------------------------------------------------------------"
      write (*, *) "|Initializing the Monte Carlo Simulation to calculate all intensities|"
      write (*, *) "----------------------------------------------------------------------"
      write (*, *)
      call getcwd(thisdir)
 if (env%eyring) write (*, '(a,i3,a)') "Eyring equation is used to compute rate constants, with mRRHO-cutoff: ", env%sthr, " cm-1"
      if (.not. env%eyring) write (*, *) "Simplified RRKM equation is used to compute rate constants"

      tav0 = 0.0_wp
      call rdshort_real('tav', tav0)

      if (tav0 .lt. 0.0_wp) then
         if (env%printlevel .eq. 3) write (*, *) "tav is negative, something went wrong, set it to 0"
         tav0 = 0.0_wp
      end if
      if (tav0 .gt. 50.0e-03_wp) then
         if (env%printlevel .eq. 3) write (*, *) "tav is much to large, something went wrong, set it to 0"
         tav0 = 0.0_wp
      end if
      write (*, '(a,x,E10.2,x,a)') "sum of half lifes of prior reactions is:", tav0, "s"

      !For plotting k(E) curves
      if (env%plot) then
      do l = 1, npairs
         write (fname2, '(a,i0,a)') "k_", nfragl, "_"//trim(fragdirs(l, 1))
         call touch(fname2)
      end do
      end if

      if (env%scaleeinthdiss .ne. 1.0_wp) write (*, '(a,f4.2)') "Scaling of IEE for H-dissociation is", env%scaleeinthdiss
      if (env%tfscale .ne. 0.0_wp) write (*, *) "Scaling of time of flight for subsequent fragmentations is", env%tfscale
      if (env%scaleker .ne. 1.0_wp) write (*, *) "Scaling of KER for subsequent fragmentations is", env%scaleker
      if (env%sumreacscale .ne. 1.0_wp) then
         write (*, *) "Scaling of sum of reaction energies for subsequent fragmentations is", env%sumreacscale
      end if
      tf = env%tf*1.0e-06_wp/(nfragl**env%tfscale) - tav0 ! 50 mukroseconds usually time in MS, p 42 Gross
      if (tf .gt. env%tf*1.0e-06_wp) then
    if (env%printlevel .eq. 3) write (*, *) "time of flight is larger than 50 microseconds, something went wrong, set it to default"
         tf = env%tf*1.0e-06_wp
      elseif (tf .le. 0.0_wp) then
         if (env%printlevel .eq. 3) write (*, *) "time of flight is negative, something went wrong, set it to 0, no time left"
         tf = 0.0_wp
      end if
      write (*, '(a,x,E10.2,x,a)') "time of flight is", tf, "s"

      nsamples = env%nsamples
      eiee = eieein

      allocate (pfrag(npairs + 1, 2)) ! also for initial structure

      ! get vibrational degrees of freedom
      inquire (file="fragment.xyz", exist=ex)
      if (ex) then
         fname = 'fragment.xyz'
      end if
      inquire (file="isomer.xyz", exist=ex)
      if (ex) then
         fname = 'isomer.xyz'
      end if
      call rdshort_int(trim(fname), nat)
      nvib = 3*nat - 6 ! vibrational degrees of freedom
      if (nvib .eq. 0) nvib = 1 ! for linear molecules !TODO not completly true
      if (nvib .lt. 0) then
         write (*, *) "nvib is negative, something went wrong"
         call printpwd
         error stop "aborting run"
      end if
      ! read barriers
      allocate (barriers(npairs), de(npairs), irc(npairs), kabs(npairs))
      allocate (kabsre(npairs))
      if (env%eyring) allocate (k_ey(npairs))
      allocate (isrearr(npairs))
      isrearr = .false.
      write (*, *) "Reading barriers and frequencies"
      write (*, *) "reaction | barrier/kcal/mol | barrier/eV | irc mode /cm -1"
      do i = 1, npairs
         inquire (file=trim(env%path)//"/"//trim(fragdirs(i, 1))//"/isomer.xyz", exist=ex)
         if (ex) isrearr(i) = .true.
         call rdshort_real(trim(env%path)//"/"//trim(fragdirs(i, 1))//'/de_'//trim(env%tslevel), de(i))
         inquire (file=trim(env%path)//"/"//trim(fragdirs(i, 1))//"/drrho_"//trim(env%geolevel), exist=ex)
         if (ex) then
            ! WARNING drrho_"//trim(env%geolevel) was in earlier dev versions drrho_"//trim(env%tslevel)
            call rdshort_real(trim(env%path)//"/"//trim(fragdirs(i, 1))//"/drrho_"//trim(env%geolevel), drrho)
         else
            write (*, *) "WARNING no drrho file found for ", trim(env%path)//"/"//trim(fragdirs(i, 1)), " setting to 0"
            drrho = 0.0_wp
         end if
         if (env%nots) then
            barriers(i) = de(i)
     write (*, *) "Take reaction energy instead of barrier of "//trim(env%path)//"/"//trim(fragdirs(i, 1))//"which is ", barriers(i)
            irc(i) = 100.0_wp ! just take reaction energy for calculations, freqs all set to 100
         else
            de(i) = de(i) - drrho
            call rdshort_real(trim(env%path)//"/"//trim(fragdirs(i, 1))//'/barrier_'//trim(env%tslevel), barriers(i))
            call rdshort_real(trim(env%path)//"/"//trim(fragdirs(i, 1))//'/ircmode_'//trim(env%geolevel), irc(i))
            if (.not. env%eyring .and. irc(i) .le. 0.0_wp) then
               write (*, *) "no IRC mode found read in just lowest frequency from:"
               write (*, *) trim(env%path)//"/"//trim(fragdirs(i, 1))//"/ts/hess"
               call chdir(trim(env%path)//"/"//trim(fragdirs(i, 1))//"/ts/hess")
               call rdlowestfreqg98(irc(i))
               call chdir(trim(thisdir))
            end if
            write (*, '(2x,a,1x,f9.1,1x,f10.2,1x,f10.1)') trim(fragdirs(i, 1)), barriers(i)*evtoau*autokcal, barriers(i), irc(i)
         end if
      end do

      allocate (tav(npairs))
      tav = 0.0_wp
      maxprob = maxloc(piee, dim=1) ! get index of maximum probability

      call findhdiss(env, fname, nat, npairs, fragdirs, ishdiss, scaleeinthdiss)

      ! account for inhomogeneous IEE distribution
      if (env%cneintscale .and. nfragl .eq. 1) then ! scale IEE with active coordination number
         allocate (active_cn(npairs))
         call get_active_cn(env, fname, nat, npairs, fragdirs, active_cn)

         max_active_cn = 0.0_wp

         do i = 1, npairs
            write (*, *) "Active coordination number of fragment ", i, " is ", active_cn(i)
         end do

         max_active_cn = maxval(active_cn)
         write (*, *) "Maximal active coordination number is ", max_active_cn

         do i = 1, npairs
            write (*, *) "Scale IEE for reaction ", i, " with ", (active_cn(i)/max_active_cn)**iee_cn_scale
         end do
      end if

      if (env%esim) call getcwd(env%path)
      !pairdir is pairdir  of fragment which gets fragmented here
      if (nfragl .gt. 1) then
         if (env%printlevel .eq. 3) write (*, *) "env%path is", trim(env%path)
         inquire (file='isomer.xyz', exist=ex) !TODO FIXME, again this uglyness
         if (ex) then
            pairdir = env%path
         else
            ilen = len(trim(env%path)) ! ??TODO FIXME env%path equals cwd ??
            pairdir = env%path(1:ilen - 2) ! fragment number only goes up to 2 so this is ok
         end if

         if (nfragl .gt. 1 .and. env%esim) env%path = ".."  ! TODO FIXME this is ugly changing it backe here

         call rdshort_real(trim(pairdir)//"/sumreac_"//trim(env%tslevel), sumreac)
         inquire (file=trim(pairdir)//"/pair.xyz", exist=ex)
         if (ex) then
            call rdshort_int(trim(pairdir)//"/pair.xyz", nat0) ! pair.xyz is correct, if isomer we dont compute KER as there is no KER and nat0 is 0
         else
            nat0 = 0
         end if
         call rdshort_real(trim(pairdir)//"/barrier_"//trim(env%tslevel), ea0)
         call rdshort_real(trim(pairdir)//"/de_"//trim(env%tslevel), de0)

         sumreac = sumreac*env%sumreacscale
         if (sumreac .lt. 0.0_wp) sumreac = 0.0_wp ! fragment cannot gain energy from reaction ? TODO FIXME

         if (env%printlevel .eq. 3) write (*, *) "sumreac (scaled) ea0, de0 are read from:"//trim(pairdir), sumreac, Ea0, de0
         if (nat0 .gt. 0) then ! if parent dir is a fragment pair and no isomer
            nvib0 = (3*nat0) - 7 ! vibrational degrees of freedom of transition state of last reaction -7 ??
            if (env%printlevel .gt. 1) write (*, *) "nvib0 is", nvib0
         else
            if (env%printlevel .gt. 1) write (*, *) "fragment stems from rearrangement, no KER is computed"
         end if
         ! allocate(KERold(nsamples)) ! TODO remove, deprecated
         KERavold = 0.0_wp
         ! read old KERav
         ! As we compute KER always for the last fragmentation, keravold is only available for nfragl > 2
         if (nfragl .gt. 2) then ! check here also if parent dir is a fragment pair and no isomer
            if (env%calcKER) then
               !call rdshort_real(trim(pairdir)//"/keravold",keravold)
               inquire (file=trim(pairdir)//"/keravold", exist=ex) ! doesnt exist if fragment stems from rearragment
               if (ex) then
                  ! open(newunit=ich, file = trim(pairdir)//"/keravold")
                  ! for exact ker and not KERaveraged
                  ! do i = 1, nsamples
                  !     read (ich,*) KERold(i)
                  ! end do
                  ! close(ich)

                  call rdshort_real(trim(pairdir)//"/keravold", keravold)
                  if (env%printlevel .gt. 1) write (*, *) "KER of previous reactions is ", keravold, " eV"
               else
                  if (env%printlevel .gt. 1) write (*, *) "no keravold file found, setting to 0"
                  KERavold = 0.0_wp
               end if
            end if
         end if
      else
         sumreac = 0.0_wp
      end if
      pfrag(:, :) = 0
      KERav = 0
      sumpiee = 0

      if (env%eatomscale) then
         inquire (file='eatomscale', exist=ex)
         if (ex) then
            call rdshort_real('eatomscale', efold)
            write (*,'(a,f4.2)') "energy is scaled according to proportion of atom with other fragment with", efold
         else
            efold = 1.0_wp
         end if
      end if

      ! generate array (nincr x npairs) of barriers for different energies
      ! we then get for a given energy the index of the nearest temperature
      ! at which we computed the barrier
      ! call xtbthermo starting fragment

      if (env%eyring) then
         nincr = env%nthermosteps ! number of increments for calculating barriers TODO check for convergence, 500 doable 2000 very expensive, but still minor deviations here ...
      elseif (env%eyzpve) then
         nincr = 1 !only ZPVE taken into account
      else
         nincr = 1 !only ZPVE taken into account
      end if

      ! problematic for CID ... TODO FIXME
      maxiee = maxval(eiee) ! is buggy, dont know why ..., seems to work now...

      ! not the best solution but should be ok for now ..
      if (env%mode == "cid") then
         maxiee = maxiee + env%cid_elab
      end if

      ! energy distribution can exceed for very large molecules eimp0, so we have to check this
      if (env%mode == "ei") then
         if (maxiee .lt. 0 .or. maxiee .gt. env%eimp0) then
            if (env%printlevel .eq. 3) write (*, *) "maxiee is negative or larger than eimp0, something went wrong, set it to eimp0"
            maxiee = env%eimp0
         end if
      end if

      ! For CID, due to collisions we can have more energy but check also here
      ! that we dont exceed it by too much
      if (env%mode == "cid") then
         if (maxiee .lt. 0 .or. maxiee .gt. env%cid_elab*1.5_wp) then
          if (env%printlevel .eq. 3) write (*, *) "maxiee is negative or larger than elab, something went wrong, set it to 1.5 elab"
            maxiee = env%cid_elab*1.5_wp
         end if
      end if

      allocate (ts_rrhos(nincr))
      allocate (start_rrhos(nincr))
      allocate (f1_rrhos(nincr))
      allocate (f2_rrhos(nincr))
      allocate (alldgs(npairs, nincr)) ! array of barriers for different energies
      allocate (dgs(nincr))
      allocate (dgeas(nincr))

      ! barrier of backward reaction
      allocate (alldgsre(npairs, nincr)) ! array of barriers for different energies
      allocate (barriersre(npairs))
      call xtbthermo(env, nincr, nvib, maxiee, start_rrhos, .true.)

      do j = 1, npairs
         call rdshort_real(trim(env%path)//"/"//trim(fragdirs(j, 1))//"/earrho_"//trim(env%geolevel), earrho)
         ! for DFT spectra
         if (env%bhess) then
            call chdir(trim(env%path)//"/"//trim(fragdirs(j, 1))//"/ts/bhess")
            call xtbthermo(env, nincr, nvib, maxiee, ts_rrhos, .true.) ! for bhess we have to use constant ithr
            call chdir(trim(thisdir))
         else ! currently not supported TODO FIXME
            ! we cannot take hess here as we have no TS!! TODO CHECKME
            call rdshort_int(trim(env%path)//"/"//trim(fragdirs(j, 1))//"/ts/nmode", nmode)

            if (nmode .eq. 0) then
               write (*, *) "no irc mode found for ", trim(env%path)//"/"//trim(fragdirs(j, 1))
               write (*, *) "we take bhess instead"
               call chdir(trim(env%path)//"/"//trim(fragdirs(j, 1))//"/ts/bhess")
               call xtbthermo(env, nincr, nvib, maxiee, ts_rrhos, .true.) ! for bhess we have to use constant ithr
               call chdir(trim(thisdir))
            else
               call chdir(trim(env%path)//"/"//trim(fragdirs(j, 1))//"/ts/hess") ! try alternatively hess
               call xtbthermo(env, nincr, nvib, maxiee, ts_rrhos, .false., irc(j)) ! for hess we invert all except for the imaginary mode
               call chdir(trim(thisdir))
            end if

         end if

         ! we have to subtract the ZPVE from the barrier
         barriers(j) = barriers(j) - earrho
         ! now add thermal corrections
         dgeas = barriers(j) + (ts_rrhos - start_rrhos)*autoev !

         if (index(fragdirs(j, 3), 'p') .ne. 0) then
            call chdir(trim(env%path)//"/"//trim(fragdirs(j, 2)))
            call xtbthermo(env, nincr, nvib, maxiee, f1_rrhos, .true.)
            call chdir(trim(thisdir))
            call chdir(trim(env%path)//"/"//trim(fragdirs(j, 3)))
            call xtbthermo(env, nincr, nvib, maxiee, f2_rrhos, .true.)
            call chdir(trim(thisdir))
         else
            call chdir(trim(env%path)//"/"//trim(fragdirs(j, 1)))
            call xtbthermo(env, nincr, nvib, maxiee, f1_rrhos, .true.)
            f2_rrhos = 0.0_wp
            call chdir(trim(thisdir))
         end if

         dgs = de(j) + (f1_rrhos + f2_rrhos - start_rrhos)*autoev

         ! Due to PES errors of the geolevel method sometimes
         ! DG is larger than DGEAS, we take the larger one than
         ! as this is the correct barrier
         do k = 1, nincr
            if (dgs(k) .gt. dgeas(k)) then
               alldgs(j, k) = dgs(k)
            else
               alldgs(j, k) = dgeas(k)
            end if
         end do

         alldgsre(j, :) = alldgs(j, :) - dgs(:) ! reverse reaction barrier, can be 0 if de is taken instead of barrier, which is also ok
         call chdir(trim(thisdir))
      end do

      do i = 1, nsamples
         ! cannot go below 0
         eiee(i) = eiee(i) - sumreac
         if (eiee(i) .lt. 0.0_wp) eiee(i) = 0.0_wp
      end do

      ! simulate collisions by adding collision energy on top of IEE
      if (env%mode == "cid") then
         call simcid(env, tf, maxiee, eiee, piee, nsamples, alldgs, nincr, &
         & nfragl, KERavold, de0, Ea0, nat0, nvib0, npairs, ishdiss, scaleeinthdiss, isrearr, .false.)
      end if

      deallocate (ts_rrhos, start_rrhos, f1_rrhos, f2_rrhos, dgs, dgeas)
      !Kerexpl TODO if (nfragl.gt. 1 .and. .not. isrearr  .and. env%calcKER) open(newunit=ich, file = "kerav")

      !!!!!!!!!!!!! FIRST ISOMERS!!!!!!!!!!
      write (*, *) "Compute Isomer equilibrium"
      do i = 1, nsamples
         ! check if starting fragment stems from a fragmenation, not a rearrangement and if KER should be considered, nat0 is 0 if parent frag is isomer, stems from rearrangement
         if (nfragl .gt. 1 .and. nat0 .gt. 0 .and. env%calcKER) then
            ! we compute here the KER of the last reaction which should be missing here
            ! for both fragments
            !TODO if we scale the energy of H-dissociations we should do this here too??? wrong if not !!!! TODO
            if (nat0 .eq. 2) then ! diatomic molecule
               KER_Eex = (eiee(i) - Ea0)
               write (*, *) "this case shouldnt exist"
            else
               KER_Eex = (eiee(i) - Ea0)/(0.44_wp*nvib0) ! J. Chem. Phys. 48, 4093 (1968)
            end if
            ! nvib0 vibrational degrees of freedom of transition state of last reaction
            KER_Ear = (0.33_wp*(Ea0 - de0)) ! J. Chem. Phys. 61, 1305 (1974)
            ! no negative kinetic energies allowed
            if (KER_Eex .lt. 0.0_wp) KER_Eex = 0.0_wp
            if (KER_Ear .lt. 0.0_wp) KER_Ear = 0.0_wp
            ! add here average KER of prior fragmentations
            KER = KER_Eex + KER_Ear + keravold
       !!TODO KERexpl KER =  KER_Eex + KER_Ear + KERold(i)
            if (KER .lt. 0.0_wp) KER = 0.0_wp ! can sometimes happen

            if (env%scaleker .ne. 1.0_wp) KER = KER*env%scaleker
            if (env%printlevel .eq. 3) write (*, *) "KER, KER_Eex, KER_EAR are:", KER, KER_Eex, KER_EAR
            !Kerexpl TODO  write(ich,*)    KER
         else
            KER = 0.0_wp

         end if

         !  write(*,*) "KER is:", KER
         ! subtract energy which was already consumed to get to this point
         ! maybe we can add here later collision energy for CID ??TODO??
         !  if (KER .gt. 0.0_wp) then
         !      eiee(i) = (eiee(i) - sumreac - KER)  ! *0.5_wp ! energy gets divided on both fragments? !
        !! probably larger fragment more energy but for now test this  !TODO FIXME
        !! maybe try to scale energy with number of atoms in fragment: e1=nat1/nat e2 =nat2/nat
         ! else

         ! cannot go below 0

         eiee(i) = eiee(i) - KER
         if (eiee(i) .lt. 0.0_wp) eiee(i) = 0.0_wp

         ! we have to scale the energy of the fragment here I guess ....
         if (env%eatomscale) eiee(i) = eiee(i)*efold

         !end if
         KERav = KERav + KER*piee(i)
         sumpiee = sumpiee + piee(i)

         ! energy can be negative here, so we have o check this
         ! megative energy means, no reaction can happen anymore
         ! this unphysical relict stems from the way
         ! the simulation is performed
         ! The "physical" way would be to keep track of the nsamples energies
         ! and their path, but this is much more complicated
         ! and (hopefully) not necessary
         if (eiee(i) .ge. 0) then

            ktotal = 0
            kabs = 0
            kabsre = 0
            ! first find ksum, i.e., fastest rate constant which is maybe different for different energies
            ! so that we can normalize to it
            ! we have to normalize everything to ktotal so that the sum of all intensities of crossed barriers is always 1*piee(i) for everycase
            do j = 1, npairs
               freq = irc(j)
               IEE = eiee(i)

               if (env%cneintscale) then
                  if (nfragl .eq. 1) IEE = IEE*active_cn(j)/max_active_cn**iee_cn_scale
               else
                  if (ishdiss(j)) then
                     IEE = IEE*scaleeinthdiss
                  end if
               end if

               ind = getindex(nincr, maxiee, IEE) ! get index of energy in array of barriers
               if (.not. env%eyring) ind = 1 ! all barriers the same
               if (env%eyzpve) ind = 1 ! only ZPVE taken into account for Eyring equation
               barriers(j) = alldgs(j, ind) ! get barrier for this energy
               barriersre(j) = alldgsre(j, ind) ! get reverse reaction barrier for this energy

               ! compute RRKM rate constant
               kabs(j) = calc_rrkm(IEE, barriers(j), nvib, freq)
               kabsre(j) = calc_rrkm(IEE, barriersre(j), nvib, freq)
               !compute eyring rate constant
               if (env%eyring) then
                  k_ey(j) = calc_eyring(IEE, barriers(j), nvib)
                  kabsre(j) = calc_eyring(IEE, barriersre(j), nvib) ! we assume the same temperature for forward and backward reaction here
                  if (env%printlevel .eq. 3) write (*, *) "kabs, k_e, ratio", kabs(j), k_ey(j), kabs(j)/k_ey(j)
                  kabs(j) = k_ey(j)
               end if
               if (isrearr(j)) ktotal = ktotal + kabs(j)
            end do

            if (env%plot) then
               do l = 1, npairs
                  if (.not. isrearr(l)) cycle
                  write (fname2, '(a,i0,a)') "k_", nfragl, "_"//trim(fragdirs(l, 1))
                  open (l, file=fname2, STATUS='OLD', POSITION='APPEND')
                  write (l, '(e10.4,2x,e10.4)') kabs(l), eiee(i)
                  close (l)
               end do
            end if

            kmax_frag = 0.0_wp
            kmax_iso = 0.0_wp
            do j = 1, npairs
               if (isrearr(j)) then
                  if (kabs(j) .gt. kmax_iso) kmax_iso = kabs(j)
               else
                  if (kabs(j) .gt. kmax_frag) kmax_frag = kabs(j)
               end if
            end do

            kmax = maxval(kabs)  ! kmax is the fastest reaction

            ! time tracking for isomer reactions
            if (kmax_iso .gt. 1.0e-8_wp .and. kmax_frag .gt. 1.0e-8_wp) then
               ! TODO inconsistent with fragmentation step later, at least for CID if coll happens, as rates change then ...
               ! TODO this is a bit inconsistent, but I dont know how to do it better unless solving the master equations for the complete reaction network
               !prem = EXP(-kmax_iso*tf) * (1 - kmax_frag/kmax_iso) ! consider only part which goes into isomer reactions
               !prem = EXP(-kmax_iso*tf) * (1- ((kmax_iso/kmax_frag)/ (kmax_iso/kmax_frag + kmax_frag/kmax_iso)))
               prem = 1 - ((kmax_iso/kmax_frag)/(kmax_iso/kmax_frag + kmax_frag/kmax_iso))*(1 - EXP(-kmax_iso*tf))
               ! should be ok, we normalize again after isomerization eq. by reading in new p0
               !if kmax_iso very small, it remains the same here ...
            else
               prem = 1.0_wp
            end if

            if (prem .gt. 1.0_wp) then
               write (*, *) "prem is too large", prem, kmax_iso, kmax_frag, tf
               prem = 1.0_wp
            end if

            if (env%printlevel .eq. 3) write (*, *) "Prem is at eiee with kmax", prem, IEE, kmax
            !now compute relative rates at eiee(i)
            do j = 1, npairs
               if (.not. isrearr(j)) cycle
               if (ktotal .gt. 0.0_wp) then
                  krel = kabs(j)/ktotal ! normalize to ktotal so that sum of all rates equals 1
               else
                  krel = 0.0_wp ! set to 0 if nor reaction is possible
               end if

               if (env%printlevel .eq. 3) write (*, *) "KREL IS", krel, j
               ! compute relative probabilities

               ! normalization to 1 so that kabsre /kabs + kabs/kabsre = 1
               if (kabs(j) .gt. 1.0e-8_wp .and. kabsre(j) .gt. 1.0e-8_wp) then ! normalization constant diverges if either kabs or kabsre is 0
                  pfrag(j, 1) = pfrag(j, 1) + 1.0_wp*piee(i)*krel*(1 - prem)*(kabs(j)/kabsre(j))/(kabs(j)/kabsre(j) &
                  & + kabsre(j)/kabs(j))! part that reacts to isomer
                  pfrag(npairs + 1, 1) = pfrag(npairs + 1, 1) + 1.0_wp*piee(i)*krel*(1 - prem) &
                                        & *(kabsre(j)/kabs(j))/(kabs(j)/kabsre(j) + kabsre(j)/kabs(j)) ! part that reacts back to the original fragment
               elseif (kabs(j) .gt. 1.0e-8_wp .and. kabsre(j) .le. 1.0e-8_wp) then ! everything reacts to isomer
                  pfrag(j, 1) = pfrag(j, 1) + 1.0_wp*piee(i)*krel*(1 - prem)
               end if
               if (pfrag(j, 1) .lt. 0.0_wp) write (*, *) "pfrag(j,1) is small", j, krel, pfrag(j, 1)
               if (env%printlevel .eq. 3) write (*, *) "pair, krel, pfrag is", j, krel, pfrag(j, 1)
               if (isnan(pfrag(j, 1))) write (*, *) "pfrag(j,1) is NaN", j, piee(i), kabs(j), ktotal, prem
            end do

            ! rest of fragment remains intact
            pfrag(npairs + 1, 1) = pfrag(npairs + 1, 1) + piee(i)*prem ! fragment remains intact
            if (isnan(pfrag(npairs + 1, 1))) write (*, *) "initial peak is NaN", eiee(i), piee(i), prem, kmax, tf
            if (env%printlevel .eq. 3) write (*, *) "remains intact", pfrag(npairs + 1, 1)
         else
            !  no energy left for fragment prem = 1
            pfrag(npairs + 1, 1) = pfrag(npairs + 1, 1) + piee(i) ! fragment remains intact
            if (isnan(pfrag(npairs + 1, 1))) write (*, *) "initial peak is NaN", eiee(i), piee(i), prem, kmax, tf
            if (env%printlevel .eq. 3) write (*, *) "no energy", pfrag(npairs + 1, 1)
         end if
      end do

      ! after KER modifcation take average energy and compute tav
      call calc_Eavg(nsamples, eiee, piee, Eavg)
      ind = getindex(nincr, maxiee, Eavg)
      do i = 1, npairs
         if (.not. isrearr(i)) cycle
         dgavg = alldgs(i, ind)
         if (env%eyring) then
            kavg = calc_eyring(Eavg, dgavg, nvib)
         else
            kavg = calc_rrkm(Eavg, dgavg, nvib, irc(i))
         end if
         !  write(*,*) "kavg is", kavg
         tav(i) = log(2.0_wp)/kavg
         ! write(*,*) "tav is", tav(i)
      end do

      !Kerexpl TODO if (nfragl.gt. 1 .and. .not. isrearr  .and. env%calcKER) close(ich)
      ! for checking  purposes we save the sum of all probabilities
      sumpfrag0 = sum(pfrag)
      sumpfrag = sum(pfrag) ! not needed here (only for fragments), but for consistency

      !write half-life of reactions

      do i = 1, npairs
         if (.not. isrearr(i)) cycle
         tav(i) = tav(i) + tav0 ! add up tavs from one cascade
         call wrshort_real(trim(env%path)//"/"//trim(fragdirs(i, 1))//'/tav', tav(i))
         ! just easier to find this way ...
         if (index(fragdirs(i, 3), 'p') .ne. 0) then
            call wrshort_real(trim(env%path)//"/"//trim(fragdirs(i, 2))//'/tav', tav(i))
            call wrshort_real(trim(env%path)//"/"//trim(fragdirs(i, 3))//'/tav', tav(i))
         end if
      end do

      ! normalize to initial start intensity(p0)*nsamples and write
      ! at the final end we normalize everything to 1e05

      !read current intensity of fragment
      !somethings off here ????
      if (nfragl .gt. 1) then
         inquire (file="pfrag", exist=ex)
         if (ex) then
            call rdshort_real("pfrag", p0)
         else
            write (*, *) "!!!!WARNING!!!! no pfrag file found, setting to 1 but something went wrong!!!!!!"
            p0 = 1.0_wp*nsamples
         end if
      else ! for first framentation always set to 1
         p0 = 1.0_wp*nsamples
      end if

      checksum = 0.0_wp

      ! norm/sum(pfrag)
      if (sumpfrag .lt. 0.1e-10_wp) then
         write (*, *) "sumpfrag is very low, something went wrong"
         pfrag = 0.0_wp
         pfrag(npairs + 1, 1) = p0/nsamples !! what does this mean??
      else
         pfrag = pfrag*p0/sumpfrag
      end if

      do i = 1, npairs
         if (.not. isrearr(i)) cycle
         !DEL    write(*,*) trim(fragdirs(i,1)), trim(fragdirs(i,2)), trim(fragdirs(i,3))
         if (index(fragdirs(i, 3), 'p') .ne. 0) then

            call wrshort_real(trim(env%path)//"/"//trim(fragdirs(i, 2))//"/pfrag", pfrag(i, 1))
            checksum = checksum + pfrag(i, 1)
            call wrshort_real(trim(env%path)//"/"//trim(fragdirs(i, 3))//"/pfrag", pfrag(i, 2))
            checksum = checksum + pfrag(i, 2)
            ! this is a rearrangement
         else
            call wrshort_real(trim(env%path)//"/"//trim(fragdirs(i, 1))//"/pfrag", pfrag(i, 1))
            checksum = checksum + pfrag(i, 1)
         end if
      end do

      ! dont forget initial peak
      ! pfrag(npairs+1,1) = pfrag(npairs+1,1) !* norm / sumpfrag
      write (*, '(a,f10.1,a)') "Intensity of initial peak was ", p0/nsamples*100.0_wp, "%"
      write (*, '(a,f10.1,a)') "Intensity of initial peak is now ", pfrag(npairs + 1, 1)/nsamples*100.0_wp, "%"
      ! write (*, *) "sumpfrag is", sumpfrag
      !write (*, *) "Normed Initial peak is now", pfrag(npairs + 1, 1)
      call wrshort_real("pfrag", pfrag(npairs + 1, 1))
      checksum = checksum + pfrag(npairs + 1, 1)

      if (abs(checksum - p0) > 0.0001_wp) then
         write (*, *) "checksum is not qual to initial peak intensity, something went wrong"
         write (*, *) "checksum is", checksum
         write (*, *) "initial intensity", p0
         error stop "aborting run"
      end if

      !!!!!!!!!!!!!!!! END OF ISOMER MODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !return

      ! allow again collisions after equilibration of system via Isomerizations
      if (env%mode == "cid") then
         call simcid(env, tf, maxiee, eiee, piee, nsamples, alldgs, nincr, nfragl, KERavold, &
         & de0, Ea0, nat0, nvib0, npairs, ishdiss, scaleeinthdiss, isrearr, .true.)
      end if

      pfrag_prec = pfrag(npairs + 1, 1)
      pfrag(:, :) = 0

      sumpiee = 0

      ! NOW ONLY FRAGMENTATIONS
      write (*, *) "Compute fragmentation branching ratios"
      do i = 1, nsamples

         ! We donte need to compute KER here again
         ! KERav = KERav + KER*piee(i)

         ! scale energy again with eatomscale
         if (env%eatomscale) eiee(i) = eiee(i)*efold

         sumpiee = sumpiee + piee(i)
         ! energy can be negative here, so we have o check this
         ! negative energy means, no reaction can happen anymore
         ! this unphysical relict stems from the way
         ! the simulation is performed
         ! The "physical" way would be to keep track of the nsamples energies
         ! and their path, but this is much more complicated
         ! and (hopefully) not necessary

         if (eiee(i) .ge. 0) then

            ktotal = 0
            kabs = 0
            kabsre = 0
            ! first find ksum, i.e., fastest rate constant which is maybe different for different energies
            ! so that we can normalize to it
            ! we have to normalize everything to ktotal that the sum of all intensities of crossed barriers is always 1*piee(i) for everycase
            !srt out too high barriers here

            do j = 1, npairs

               freq = irc(j)
               IEE = eiee(i)

               if (env%cneintscale) then
                  if (nfragl .eq. 1) IEE = IEE*active_cn(j)/max_active_cn**iee_cn_scale
               else
                  if (ishdiss(j)) then
                     IEE = IEE*scaleeinthdiss
                  end if
               end if

               ind = getindex(nincr, maxiee, IEE) ! get index of energy in array of barriers

               if (.not. env%eyring) ind = 1 ! all barriers the same
               if (env%eyzpve) ind = 1 ! only ZPVE taken into account
               barriers(j) = alldgs(j, ind) ! get barrier for this energy

               barriersre(j) = alldgsre(j, ind) ! get reverse reaction barrier for this energy
               kabs(j) = calc_rrkm(IEE, barriers(j), nvib, freq)
               kabsre(j) = calc_rrkm(IEE, barriersre(j), nvib, freq)
               if (env%eyring) then
                  k_ey(j) = calc_eyring(IEE, barriers(j), nvib)
                  kabsre(j) = calc_eyring(IEE, barriersre(j), nvib) ! TODO we assume the same temperature for forward and backward reaction here
                  if (env%printlevel .eq. 3) write (*, *) "kabs, k_e, ratio", kabs(j), k_ey(j), kabs(j)/k_ey(j)
                  kabs(j) = k_ey(j)
               end if

               ! only consider fragmentations here
               if (.not. isrearr(j)) ktotal = ktotal + kabs(j)
            end do

            if (env%plot) then
               do l = 1, npairs
                  if (isrearr(l)) cycle
                  write (fname2, '(a,i0,a)') "k_", nfragl, "_"//trim(fragdirs(l, 1))
                  open (l, file=fname2, STATUS='OLD', POSITION='APPEND')
                  write (l, '(e10.4,2x,e10.4)') kabs(l), eiee(i)
                  close (l)
               end do
            end if

            kmax_frag = 0.0_wp
            kmax_iso = 0.0_wp
            do j = 1, npairs
               if (isrearr(j)) then
                  if (kabs(j) .gt. kmax_iso) kmax_iso = kabs(j)
               else
                  if (kabs(j) .gt. kmax_frag) kmax_frag = kabs(j)
               end if
            end do

            kmax = maxval(kabs)  ! kmax is the fastest reaction not exact but ok here? sum of k is wrong here

            ! time tracking for fragmentations
            ! check if there is any isomerization reaction
            if (kmax_iso .gt. 0.00001_wp) then
               prem = EXP(-kmax_frag*tf) ! consider now only part which goes into fragment reactions
               !part from isomer eq. is already considered in isomer eq. part we dont to consider it again would be double counting
               !* ((kmax_frag/kmax_iso) + pfrag_prec/nsamples)
            else
               prem = EXP(-kmax_frag*tf)
            end if

            ! just for checking purposes
            if (prem .gt. 1.0_wp) then
               if (env%printlevel .eq. 3) write (*, *) "prem is too large", prem, kmax_iso, kmax_frag, tf
               prem = 1.0_wp
            end if
            !
            if (env%printlevel .eq. 3) write (*, *) "Prem is at eiee with kmax", prem, IEE, kmax

            do j = 1, npairs
               if (isrearr(j)) cycle
               if (ktotal .gt. 0.0_wp) then
                  krel = kabs(j)/ktotal ! normalize to ktotal so that sum of all rates equals 1
               else
                  krel = 0.0_wp ! set to 0 if nor reaction is possible
               end if
               if (env%printlevel .eq. 3) write (*, *) "KREL IS", krel, j

               ! compute relative probabilities
               pfrag(j, 1) = pfrag(j, 1) + 1.0_wp*piee(i)*krel*(1 - prem) ! now no back reaction to parent ion
               if (pfrag(j, 1) .lt. 0.0_wp) write (*, *) "pfrag(j,1) is small", j, krel, pfrag(j, 1)
               if (env%printlevel .eq. 3) write (*, *) "pair, krel, pfrag is", j, krel, pfrag(j, 1)
               if (isnan(pfrag(j, 1))) write (*, *) "pfrag(j,1) is NaN", j, piee(i), kabs(j), ktotal, prem
            end do

            ! rest of fragment remains intact

            pfrag(npairs + 1, 1) = pfrag(npairs + 1, 1) + piee(i)*prem ! fragment remains intact
            if (isnan(pfrag(npairs + 1, 1))) write (*, *) "initial peak is NaN", eiee(i), piee(i), prem, kmax, tf
            if (env%printlevel .eq. 3) write (*, *) "remains intact", pfrag(npairs + 1, 1)
         else
            !   no energy left for fragment prem = 1
            pfrag(npairs + 1, 1) = pfrag(npairs + 1, 1) + piee(i) ! fragment remains intact
            if (isnan(pfrag(npairs + 1, 1))) write (*, *) "initial peak is NaN", eiee(i), piee(i), prem, kmax, tf
            if (env%printlevel .eq. 3) write (*, *) "no energy", pfrag(npairs + 1, 1)
         end if
      end do

      ! after KER modifcation take average energy and compute tav
      call calc_Eavg(nsamples, eiee, piee, Eavg)
      ind = getindex(nincr, maxiee, Eavg)

      do i = 1, npairs
         if (isrearr(i)) cycle
         dgavg = alldgs(i, ind)
         if (env%eyring) then
            kavg = calc_eyring(Eavg, dgavg, nvib)
         else
            kavg = calc_rrkm(Eavg, dgavg, nvib, irc(i))
         end if
         tav(i) = log(2.0_wp)/kavg ! average half-life of reactions
      end do

      !Kerexpl TODO if (nfragl.gt. 1 .and. .not. isrearr  .and. env%calcKER) close(ich)
      ! for checking  purposes we save the sum of all probabilities
      sumpfrag0 = sum(pfrag)

      !Kerexpl TODO if (nfragl.gt. 1 .and. .not. isrearr  .and. env%calcKER) close(ich)
      ! for checking  purposes we save the sum of all probabilities
      sumpfrag0 = sum(pfrag)

      ! special mode for testing, only meaningful for alkanes
      if (env%chargeatomscale) then
         write (*, *) "Assgining charges of fragments according to number of atoms"
         do i = 1, npairs
            if (index(fragdirs(i, 3), 'p') .ne. 0) then

               call rdshort_int(trim(env%path)//"/"//trim(fragdirs(i, 2))//'/fragment.xyz', natf1)
               call rdshort_int(trim(env%path)//"/"//trim(fragdirs(i, 3))//'/fragment.xyz', natf2)
               if (natf1 .gt. natf2) then
                  qf1 = 1.0_wp
                  qf2 = 0.0_wp
               else
                  qf1 = 0.0_wp
                  qf2 = 1.0_wp
               end if
               pfragtemp = pfrag(i, 1)
               pfrag(i, 1) = pfragtemp*qf1
               pfrag(i, 2) = pfragtemp*qf2
            end if
         end do
      else
         ! weight intensity according to statistical chage of fragments
         do i = 1, npairs
            if (index(fragdirs(i, 3), 'p') .ne. 0) then
               if (env%ignoreip) then
                  write (*, *) "Warning: ignoring IPs and setting statistical charge to 1"
                  qf1 = 1.0_wp
                  qf2 = 1.0_wp
               else
                  call rdshort_real(trim(env%path)//"/"//trim(fragdirs(i, 2))//'/statchrg', qf1)
                  call rdshort_real(trim(env%path)//"/"//trim(fragdirs(i, 3))//'/statchrg', qf2)
                  if (abs(qf1 + qf2 - 1) .gt. 0.1e-6_wp) then
                     write (*, *) "statistical charges are wrong for", fragdirs(i, 2), qf1, fragdirs(i, 3), qf2
                  end if
               end if

               pfragtemp = pfrag(i, 1)
               pfrag(i, 1) = pfragtemp*qf1
               pfrag(i, 2) = pfragtemp*qf2
            end if
         end do
      end if

      sumpfrag = sum(pfrag)

      if (abs(sumpfrag - sumpfrag0) .gt. 0.1e-6_wp) then
         write (*, *) "Warning, statistical charge evaluation went wrong, aborting"
         write (*, *) "sum is now: ", sumpfrag, "sumpfrag was: ", sumpfrag0
         if (.not. env%ignoreip) error stop "aborting run"
      end if

      ! weight intensity according to number of atoms ration
      ! the idea is that the larger fragment gets more energy but not all energy!!!
      ! write factor here to file
      if (env%eatomscale) then
         do i = 1, npairs
            ! only for fragmentpairs relevant
            if (index(fragdirs(i, 3), 'p') .ne. 0) then
               call rdshort_int(trim(env%path)//"/"//trim(fragdirs(i, 2))//'/fragment.xyz', natf1)
               call rdshort_int(trim(env%path)//"/"//trim(fragdirs(i, 3))//'/fragment.xyz', natf2)
               ef1 = 1.0_wp*natf1/(natf1 + natf2)
               ef2 = 1.0_wp*natf2/(natf1 + natf2)
               inquire (file='eatomscale', exist=ex)
               if (ex) then
                  call rdshort_real('eatomscale', efold) ! for subsequent fragmentations we have to multiply with old factor
                  ef1 = ef1*efold
                  ef2 = ef2*efold
               end if
               call wrshort_real(trim(env%path)//"/"//trim(fragdirs(i, 2))//'/eatomscale', ef1)
               call wrshort_real(trim(env%path)//"/"//trim(fragdirs(i, 3))//'/eatomscale', ef2)
            end if
         end do
      end if

      !write half-life of reactions

      do i = 1, npairs
         if (isrearr(i)) cycle
         tav(i) = tav(i) + tav0 ! add up tavs from one cascade
         call wrshort_real(trim(env%path)//"/"//trim(fragdirs(i, 1))//'/tav', tav(i))
         ! just easier to find this way ...
         if (index(fragdirs(i, 3), 'p') .ne. 0) then
            call wrshort_real(trim(env%path)//"/"//trim(fragdirs(i, 2))//'/tav', tav(i))
            call wrshort_real(trim(env%path)//"/"//trim(fragdirs(i, 3))//'/tav', tav(i))
         end if
      end do

      ! we read now p0 again produced by isomerization eq.
      call rdshort_real("pfrag", p0)

      checksum = 0.0_wp

      ! norm/sum(pfrag)
      if (sumpfrag .lt. 0.1e-10_wp) then
         write (*, *) "sumpfrag is very low, something went wrong"
         pfrag = 0.0_wp
         pfrag(npairs + 1, 1) = p0/nsamples !! what does this mean??
      else
         pfrag = pfrag*p0/sumpfrag
      end if

      do i = 1, npairs
         if (isrearr(i)) cycle
         if (index(fragdirs(i, 3), 'p') .ne. 0) then

            call wrshort_real(trim(env%path)//"/"//trim(fragdirs(i, 2))//"/pfrag", pfrag(i, 1))
            checksum = checksum + pfrag(i, 1)
            call wrshort_real(trim(env%path)//"/"//trim(fragdirs(i, 3))//"/pfrag", pfrag(i, 2))
            checksum = checksum + pfrag(i, 2)
         end if
      end do

      ! dont forget initial peak
      ! pfrag(npairs+1,1) = pfrag(npairs+1,1) !* norm / sumpfrag
      write (*, '(a,f10.1,a)') "Intensity of initial peak was ", p0/nsamples*100.0_wp, "%"
      write (*, '(a,f10.1,a)') "Intensity of initial peak is now ", pfrag(npairs + 1, 1)/nsamples*100.0_wp, "%"
      ! write (*, *) "sumpfrag is", sumpfrag
      ! write (*, *) "Normed Initial peak is now", pfrag(npairs + 1, 1)
      call wrshort_real("pfrag", pfrag(npairs + 1, 1))
      checksum = checksum + pfrag(npairs + 1, 1)

      if (abs(checksum - p0) > 0.0001_wp) then
         write (*, *) "checksum is not qual to initial peak intensity, something went wrong"
         write (*, *) "checksum is", checksum
         write (*, *) "initial intensity", p0
         error stop "aborting run"
      end if

      ! END OF FRAGMENTATION MODE

      ! copy kerav data for all fragments to pairdirs for mcsimu later so its easier to find and read
      if (nfragl .gt. 1) then
         KERav = KERav/sumpiee
         KERav = KERav + keravold ! kerav is sum of all old keravs
          !! KER
         if (env%printlevel .gt. 1) write (*, *) "average kinetic energy release is ", KERav, "eV"
         call wrshort_real("kerav", KERav)
         do i = 1, npairs
            if (index(fragdirs(i, 3), 'p') .ne. 0) then
               ndir = 3
            else
               ndir = 1
            end if
            do j = 1, ndir
               call copy('kerav', trim(env%path)//"/"//trim(fragdirs(i, j))//'/keravold')
            end do
         end do
      end if

      if (env%mode == "cid") then
         do i = 1, npairs
            if (index(fragdirs(i, 3), 'p') .ne. 0) then
               ndir = 3
            else
               ndir = 1
            end if
            do j = 1, ndir
               call copy('sumdekin', trim(env%path)//"/"//trim(fragdirs(i, j))//'/sumdekin')
               call copy('sumdeint', trim(env%path)//"/"//trim(fragdirs(i, j))//'/sumdeint')
               call copy('x_trav', trim(env%path)//"/"//trim(fragdirs(i, j))//'/x_trav')
            end do
         end do
      end if

      if (env%mode == "cid") then
         do i = 1, npairs
            if (index(fragdirs(i, 3), 'p') .ne. 0) then
               call rdshort_real(trim(env%path)//"/"//trim(fragdirs(i, 2))//"/mass", mf1)
               call rdshort_real(trim(env%path)//"/"//trim(fragdirs(i, 3))//"/mass", mf2)
               qmass = mf1/mf2 ! mass ratio for KER distribution on fragments
               call wrshort_real(trim(env%path)//"/"//trim(fragdirs(i, 2))//"/qmass", qmass)
               qmass = mf2/mf1
               call wrshort_real(trim(env%path)//"/"//trim(fragdirs(i, 3))//"/qmass", qmass)
            end if
         end do
      end if
      write (*, *)
      write (*, *) "---------------------------------"
      write (*, *) "|Monte Carlo simulation finished|"
      write (*, *) "---------------------------------"
      write (*, *)
   end subroutine montecarlo

   ! Pre-compute thermal contributions to barriers with xtb thermo module
   ! across all energy increments
   subroutine xtbthermo(env, nincr, nvib, maxiee, rrhos, bhess, irc)
      implicit none
      character(len=2024) :: jobcall
      character(len=2024) :: temps ! input of temperatures for xtb thermo (rounded to integers)
      real(wp) :: T
      real(wp), intent(in) :: maxiee
      real(wp), optional :: irc
      real(wp) :: energy
      real(wp), intent(out) :: rrhos(nincr)
      integer, intent(in) :: nvib ! number of vibrational degrees of freedom of starting fragment (temperature for one fragmentation level the same for all reactions)
      integer :: temp
      real(wp) :: estep
      integer :: sthr ! threshold for mRRHO cutoff
      integer :: ithr ! imag cutoff for thermo calculation, because small imaginary frequencies often present on not perfectly optimized structures
      integer, intent(in) :: nincr ! energy increments
      character(len=40) :: fname
      real(wp) :: freqscale ! freqeuncy scaling factor for xtb thermo
      integer :: i
      type(runtypedata) :: env
      logical :: ex, ex2
      logical :: bhess ! if true we have to compute bhess
      logical :: tsgeo ! if true we have to compute -tsgeo ! just for testing, not in use right now

      tsgeo = .false. ! just for testing

      sthr = env%sthr
      if (present(irc)) then
         ithr = abs(irc) - 1 ! invert all imaginary frequencies below IRC
      else
         ithr = env%ithr
      end if
      freqscale = 1.0_wp ! for wb97x-3c 0.95, for xtb 1.0
      if (env%geolevel .eq. 'wb97x-3c' .and. .not. bhess) freqscale = 1.0_wp ! 0.95_wp

      inquire (file="isomer.xyz", exist=ex)
      if (ex) then
         fname = "isomer.xyz"
         ithr = 100  ! shouldnt have any imaginary frequencies
      else
         !inquire (file="pair.xyz", exist=ex) ! just for testing now
         !if (ex) then
         !   fname="pair.xyz"
         !   ithr = 100 ! shouldnt have any imaginary frequencies

         ! equals end.xyz if there is no TS
         inquire (file="ts.xyz", exist=ex)
         if (ex) then
            fname = "ts.xyz"

            ! tsgeo was just for testing
            !tsgeo = .true.

            ! end.ts should also have no imaginary frequencies
            inquire (file="endeqts", exist=ex2)
            if (ex2) then
               ithr = 100 ! shouldnt have any imaginary frequencies
            end if
         else
            inquire (file="fragment.xyz", exist=ex)
            if (ex) then
               fname = "fragment.xyz"
               ithr = 100 ! shouldnt have any imaginary frequencies
            end if
         end if
      end if

      ! divide energy range into nincr parts
      estep = maxiee/(nincr)
      energy = 0.0_wp

      !warning check this!!!
      ithr = env%ithr
      write (temps, '(a)') ""
      if (nincr .gt. 1) then
         do i = 1, nincr
            energy = estep*i
            T = calctemp(nvib, energy)
            temp = int(T) ! round to integer
            if (nincr .eq. 1) temp = 0.1_wp ! for RRKM we only take ZPVE into account
            write (temps, '(a,i0)') trim(temps)//",", temp
            if (mod(i, 10) .eq. 0) then
               ! TODO implement here hessian of different level of theory -nobhess option
               !  if (env%geolevel .ne. 'gfn2' .and. env%geolevel .ne. 'gfn1' .and. env%geolevel .ne. 'gxtb') then
               !     if (env%qmprog == 'tm') then
               !        write (jobcall, *) 'xtb thermo '//trim(fname)//' hessian --temp '//trim(temps) ! ithr set to 50 (default 20, because Bhess is strange ..)
               !        write (jobcall, *) trim(jobcall)//' --ithr ', ithr, ' --sthr ', sthr, ' --scale ', freqscale
               !        write (jobcall, *) trim(jobcall)//' > thermo.out 2>/dev/null'
               !     elseif (env%qmprog == 'orca') then
               !        write (jobcall, *) 'xtb thermo '//trim(fname)//' orca.hess --orca  --temp '//trim(temps)
               !        write (jobcall, *) trim(jobcall)//' --ithr ', ithr, ' --sthr ', sthr, ' --scale ', freqscale ! ithr set to 50 (default 20, because Bhess is strange ..)
               !        write (jobcall, *) trim(jobcall)//' > thermo.out 2>/dev/null'
               !     end if
               !  else
               !     write (jobcall,*) 'xtb thermo '//trim(fname)//' hessian --temp '//trim(temps)//
               !     write (jobcall, *) trim(jobcall)//' --ithr ',ithr ,' --sthr ', sthr, ' --scale ',freqscale,' > thermo.out 2>/dev/null'
               !  end if
               if (bhess) then
                  write (jobcall, '(a)') 'xtb thermo '//trim(fname)//' hessian --temp '//trim(temps)
                  write (jobcall, '(a,i0,a,i0,a)') trim(jobcall)//' --ithr ', ithr, ' --sthr ', sthr, ' --bhess hessian_sph '
                  write (jobcall, '(a)') trim(jobcall)//' > thermo.out 2>/dev/null'
               end if
               call execute_command_line(trim(jobcall))

               call readoutthermo(i, nincr, rrhos)
               write (temps, '(a)') ""
            end if
         end do! ithr set to 50 (default 20, because Bhess is strange ..)
      elseif (nincr .eq. 1) then
         ! if (env%geolevel .ne. 'gfn2' .and. env%geolevel .ne. 'gfn1' .and. env%geolevel .ne. 'gxtb') then
         !    if (env%qmprog == 'tm') then
         !    write (jobcall, *) 'xtb thermo '//trim(fname)//' hessian --temp 0.1 --ithr ',ithr ,' --sthr ', sthr, ' --scale ',freqscale,' > thermo.out 2>/dev/null' ! just compute ZPVE, 0 K doesnt work
         !    elseif (env%qmprog == 'orca') then
         !    write (jobcall, *) 'xtb thermo '//trim(fname)//' orca.hess --orca --temp 0.1 --ithr ',ithr ,' --sthr ', sthr, ' --scale ',freqscale,' > thermo.out 2>/dev/null' ! just compute ZPVE, 0 K doesnt work
         !    end if
         ! else
         !    ! just compute ZPVE, 0 K doesnt work
         ! write (jobcall, *) 'xtb thermo '//trim(fname)//' hessian --temp 0.1 --ithr ',ithr ,' --sthr ', sthr, ' --scale ',freqscale,' > thermo.out 2>/dev/null'
         ! end if ! just compute ZPVE, 0 K doesnt work
         if (bhess) then
            write (jobcall, '(a)') 'xtb thermo '//trim(fname)//' hessian --temp 0.1'
            write (jobcall, '(a,i0,a,i0,a)') trim(jobcall)//' --ithr ', ithr, ' --sthr ', sthr, ' --bhess hessian_sph '
            write (jobcall, '(a)') trim(jobcall)//' > thermo.out 2>/dev/null'
         end if
         call execute_command_line(trim(jobcall))
         call readoutthermo(nincr, nincr, rrhos)
      end if
   end subroutine xtbthermo

! subroutine to read out thermo data from xtb thermo output
   subroutine readoutthermo(endind, incr, rrhos)
      implicit none
      character(len=1024) :: tmp
      integer :: i, io, ich
      integer, intent(in) :: incr
      integer, intent(in) :: endind
      real(wp) :: rrhos(incr)
      real(wp) :: dumr

      open (newunit=ich, file='thermo.out')
      do
         read (ich, '(a)', iostat=io) tmp
         if (io < 0) exit
         if (index(tmp, 'T/K') .ne. 0) then
            read (ich, '(a)', iostat=io) tmp
            if (incr .eq. 1) then
               read (ich, *) dumr, dumr, dumr, dumr, rrhos(1)
               exit
            else
               ! do it in blocks of 10
               do i = endind - 9, endind
                  read (ich, *) dumr, dumr, dumr, dumr, rrhos(i)
                  !DEL write(*,*) "rrhos are", rrhos(i)
               end do
               exit
            end if
         end if
      end do
      close (ich)
   end subroutine readoutthermo

   subroutine mcesim(env, eiee, piee)
      implicit none
      real(wp):: eiee(:), piee(:)
      type(runtypedata) :: env
      integer :: ntot, npairs, nfrag, npairs2
      integer :: nfragl! fragmentation level
      integer :: i, j, k, l
      integer :: ilen
      integer :: io, ich
      real(wp) :: pfrag
      character(len=80), allocatable :: allfrags(:)
      character(len=80), allocatable   :: fragdirs(:, :) ! array of all fragment directories
      character(len=80) :: fragdir, fragkind
      character(len=1024) :: startfragdir
      character(len=2048) :: restartfrags ! list of fragments which should be used for restart
      character(len=512) :: atmp
      character(len=80) :: dum
      logical :: ex

      restartfrags = ""
      write (*, *) "starting energy simulation up to fragmentation level", env%nfrag

      call cleanup_qcxms2_for_restart(env)

      inquire (file="allfrags.dat", exist=ex)
      if (.not. ex) then
         error stop "allfrags.dat file from previous qcxms run not found, no energy simulation possible!"
      end if
      write (*, *) "reading file <allfrags.dat>"
      call rdshort_int('allfrags.dat', ntot)
      allocate (allfrags(ntot + 1))
      allfrags(1) = './' ! first allfrags is startdir
      open (newunit=ich, file='allfrags.dat')
      read (ich, *) ntot
      do i = 2, ntot + 1
         read (ich, *) allfrags(i)
         !  read(ich,'(a)') allfrags(i) for restart better but reads in " " space in front of dirs which leads to problems then ....
         ! TODO FIXME also
         !DEL write(*,*) "allfrags is", trim(allfrags(i))
      end do
      close (ich)
      do k = 1, ntot + 1
         !write (*, *) "allfrags2 is:", trim(allfrags(k))
         call chdir(trim(allfrags(k)))
         ! initial fragment has intensity 1
         ! determine level of fragmentation sequence
         nfragl = 1
         ilen = len_trim(allfrags(k))
         ! level of fragmentation sequence  p1f1 is 1 meaning level 2 of fragmentation
         do l = 1, ilen
            if (allfrags(k) (l:l) == 'p') then
               nfragl = nfragl + 1
            end if
         end do
         call getcwd(startfragdir)
         call getcwd(env%path)
         !DEL write(*,*) "current directory is", trim(env%path)
         inquire (file="fragments", exist=ex)
         if (.not. ex) then
            ! write (*, *) "no fragments file found, skip this directory"
            call rdshort_real('pfrag', pfrag)
            if (pfrag .gt. env%pthr) then !!! TODO FIXME add here also info about sumreac, tav, kerav, ....
               write (restartfrags, '(a,1x,a)') trim(restartfrags), trim(allfrags(k))
            end if
            call chdir(trim(env%startdir)) ! here was env%path written???
            cycle
         end if

         call rdshort_int(trim(startfragdir)//'/npairs2', npairs)
         allocate (fragdirs(npairs, 3))
         fragdirs = ""
         open (newunit=ich, file=trim(startfragdir)//'/fragments')
         read (ich, '(a)', iostat=io) atmp
         npairs2 = npairs
         l = 0
         do i = 1, npairs
            read (ich, *) fragdir, fragkind, dum, dum, dum, dum
            if (fragkind .eq. "isomer") then
               ! TODO maybe this would be cool ???

               read (ich, *) fragdir, dum, dum, dum
               if (env%msnoiso) then
                  npairs2 = npairs2 - 1
                  cycle ! skip isomers
               end if
               if (.not. env%msfulliso .and. nfragl .gt. 1) then
                  npairs2 = npairs2 - 1
                  cycle ! skip isomers
               end if
               inquire (file='isomer.xyz', exist=ex) !TODO FIXME, again this uglyness
               !  write(*,*) "fragdir is", trim(fragdir)
               !   write(fragdirs(i,1),'(a)') trim(env%startdir)//"/"//trim(fragdir) ! set proper paths to directories
               l = l + 1
               write (fragdirs(l, 1), '(a)') trim(fragdir)
               fragdirs(l, 2) = ''
               fragdirs(l, 3) = ''
            else
               !write(fragdirs(i,1),'(a)') trim(env%startdir)//"/"//trim(fragdir)
               l = l + 1
               write (fragdirs(l, 1), '(a)') trim(fragdir)
               do j = 1, 2
                  read (ich, *) fragdir, dum, dum, dum
                  ! write(*,*) "fragdir is", trim(fragdir)
                  !   write(fragdirs(i,j+1),'(a)') trim(env%startdir)//"/"//trim(fragdir)
                  write (fragdirs(l, j + 1), '(a)') trim(fragdir)
               end do
            end if
         end do
         if (k == 1) call wrshort_int(trim(startfragdir)//'/pfrag', env%nsamples)
         call rdshort_real(trim(startfragdir)//'/pfrag', pfrag)

         !DEL do i = 1, npairs2
         !    write (*, '(a,1x,a,1x,a)') "fragdirs is:", trim(fragdirs(i, 1)), trim(fragdirs(i, 2)), trim(fragdirs(i, 3))
         !end do
         if (npairs2 .gt. 0 .and. pfrag .ge. env%pthr .and. nfragl .le. (env%nfrag + 1)) then
            if (nfragl .gt. 1) env%path = ".."
            call montecarlo(env, npairs2, fragdirs, eiee, piee, nfragl)
         end if
         close (ich)
         deallocate (fragdirs)
         call chdir(trim(env%startdir))
       !!!DEL  write(*,*) "current directory is", trim(env%startdir)
        !!!DEL  call printpwd()
      end do

      if (len(trim(restartfrags)) .gt. 0) then
         write (*, *) "!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!"
         write (*, *) "Some fragments have significant probability and energy left and should be restarted:", trim(restartfrags)
         write (*, *) "!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!"
      end if
   end subroutine mcesim
end module mcsimu
