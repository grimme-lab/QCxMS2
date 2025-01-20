!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! most of the code is taken from plotms (https://github.com/qcxms/PlotMS) under the  LGPL_3.0 license
! compute the isotopic mass distribution for
! given chemical formula (ntot atoms with OZ iat_save())
! returns nsig masses in mass() with probability
! this is a quick and dirty MC algorithm and its run-time depends
! critically on the number of trials nrnd

! TK See: ATOMIC WEIGHTS OF THE ELEMENTS: REVIEW 2000 (IUPAC Technical Report)
! TK Pure Appl. Chem., Vol. 75, No. 6, pp. 683800, 2003.
module isotope_pattern
   use xtb_mctc_accuracy, only: wp
   use iomod
   use qcxms2_data
   implicit none
   public
   real(wp) ::  ams(107)  ! TODO FIXME ugly
   data ams/1.00790_wp, 4.00260_wp, 6.94000_wp, 9.01218_wp, &
   &  10.81000_wp, 12.01100_wp, 14.00670_wp, 15.99940_wp, 18.99840_wp, &
   &  20.17900_wp, 22.98977_wp, 24.30500_wp, 26.98154_wp, 28.08550_wp, &
   &  30.97376_wp, 32.06000_wp, 35.45300_wp, 39.94800_wp, 39.09830_wp, &
   &  40.08000_wp, 44.95590_wp, 47.90000_wp, 50.94150_wp, 51.99600_wp, &
   &  54.93800_wp, 55.84700_wp, 58.93320_wp, 58.71000_wp, 63.54600_wp, &
   &  65.38000_wp, 69.73500_wp, 72.59000_wp, 74.92160_wp, 78.96000_wp, &
   &  79.90400_wp, 83.80000_wp, 85.46780_wp, 87.62000_wp, 88.90590_wp, &
   &  91.22000_wp, 92.90640_wp, 95.94000_wp, 98.90620_wp, 101.0700_wp, &
   &  102.9055_wp, 106.4000_wp, 107.8680_wp, 112.4100_wp, 114.8200_wp, &
   &  118.6900_wp, 121.7500_wp, 127.6000_wp, 126.9045_wp, 131.3000_wp, &
   &  132.9054_wp, 137.3300_wp,  &
   &  138.91, 140.12, 140.91, 144.24, 147.00, 150.36, 151.97, 157.25, &
   &  158.93, 162.50, 164.93, 167.26, 168.93, 173.04, 174.97, &
   &  178.4900_wp, 180.9479_wp, &
   &  183.8500_wp, 186.2070_wp, 190.2000_wp, 192.2200_wp, 195.0900_wp, &
   &  196.9665_wp, 200.5900_wp, 204.3700_wp, 207.2000_wp, 208.9804_wp, &
   &  209., 210., 222., 21*0.000_wp/

contains

   subroutine getmass(env, fname, maxatm, rnd, index_mass, exact_intensity, isotope_masses)
      use mctc_env, only: error_type, get_argument!, fatal_error
      use mctc_io, only: structure_type, read_structure, write_structure, &
         & filetype, get_filetype, to_symbol, to_number
      use utility

      implicit none
      type(runtypedata) :: env
      integer :: nat, iti, z_chrg
      integer, intent(in) :: maxatm
      integer :: i, iso, n, loop, index_mass
      integer, parameter :: nrnd = 50000
      integer :: store_int(1000)
      integer :: niso(200), nruns
      real(wp) :: rnd(nrnd, maxatm)
      real(wp) :: list_masses(nrnd)
      integer, allocatable :: iat(:)
      real(wp) :: atmasses(10), molmasses(10), p1, p2, x, xmass
      real(wp) :: mipmax, iipmax ! mass and intensity of highest isotope peak
      integer :: indexipmax ! index of peak with highest intensity
      real(wp) :: current_mass, pfrag
      real(wp), allocatable :: exact_intensity(:)
      real(wp), allocatable :: isotope_masses(:)
      character(len=*) :: fname
      !real(wp) :: exact_intensity(1000)

      real(wp) :: prob(200, 10), massiso(200, 10)
      real(wp) :: r, sum_prob
      logical  :: there

      type(structure_type) :: mol
      type(error_type), allocatable :: error

      !  1 H  (Hydrogen)
      niso(1) = 2
      prob(1, 1) = 0.0115_wp
      prob(1, 2) = 99.9885_wp
      massiso(1, 1) = 2.014101_wp
      massiso(1, 2) = 1.007825_wp

      !  2 He  (Helium)
      niso(2) = 2
      prob(2, 1) = 0.000134_wp
      prob(2, 2) = 99.99866_wp
      massiso(2, 1) = 3.016029309_wp
      massiso(2, 2) = 4.002603249_wp

      ! 3 Li (Lithium)

      niso(3) = 2
      prob(3, 1) = 7.59_wp
      prob(3, 2) = 92.41_wp
      massiso(3, 1) = 6.0151223_wp
      massiso(3, 2) = 7.0160041_wp

      ! 4 Be (Beryllium)

      niso(4) = 1
      prob(4, 1) = 100.0_wp
      massiso(4, 1) = 9.0121822_wp

      !  5 B  (Boron)
      niso(5) = 2
      prob(5, 1) = 19.9_wp
      prob(5, 2) = 80.1_wp
      massiso(5, 1) = 10.0129370_wp
      massiso(5, 2) = 11.0093055_wp

      !  6 C  (Carbon)
      niso(6) = 2
      prob(6, 1) = 1.07_wp
      prob(6, 2) = 98.93_wp
      massiso(6, 1) = 13.003354_wp
      massiso(6, 2) = 12.00_wp

      !  7 N  (Nitrogen)
      niso(7) = 2
      prob(7, 1) = 0.368_wp
      prob(7, 2) = 99.632_wp
      massiso(7, 1) = 15.000108_wp
      massiso(7, 2) = 14.003074_wp

      !  8 O  (Oxygen)
      niso(8) = 3
      prob(8, 1) = 0.038_wp
      prob(8, 2) = 0.205_wp
      prob(8, 3) = 99.757_wp
      massiso(8, 1) = 16.999131_wp
      massiso(8, 2) = 17.999159_wp
      massiso(8, 3) = 15.994914_wp

      !  9 F  (Fluorine)
      niso(9) = 1
      prob(9, 1) = 100.0_wp
      massiso(9, 1) = 18.998403_wp

      ! 10 Ne (Neon)
      niso(10) = 3
      prob(10, 1) = 90.48_wp
      prob(10, 2) = 0.27_wp
      prob(10, 3) = 9.25_wp
      massiso(10, 1) = 19.992440176_wp
      massiso(10, 2) = 20.99384674_wp
      massiso(10, 3) = 21.99138550_wp
      ! 11 Na (Sodium)
      niso(11) = 1
      prob(11, 1) = 100.0_wp
      massiso(11, 1) = 22.98976966_wp

      ! 12 Mg (Magnesium)
      niso(12) = 3
      prob(12, 1) = 78.99_wp
      prob(12, 2) = 10.00_wp
      prob(12, 3) = 11.01_wp
      massiso(12, 1) = 23.98504187_wp
      massiso(12, 2) = 24.98583700_wp
      massiso(12, 3) = 25.98259300_wp

      !  13 Al (Aluminium)
      niso(13) = 1
      prob(13, 1) = 100.0_wp
      massiso(13, 1) = 26.981538_wp

      ! 14 Si (Silicon)
      niso(14) = 3
      prob(14, 1) = 92.223_wp
      prob(14, 2) = 4.685_wp
      prob(14, 3) = 3.092_wp
      massiso(14, 1) = 27.976926_wp
      massiso(14, 2) = 28.976494_wp
      massiso(14, 3) = 29.973770_wp

      ! 15 P  (Phosphorus)
      niso(15) = 1
      prob(15, 1) = 100.0_wp
      massiso(15, 1) = 30.973761_wp

      ! 16 S  (Sulfur)
      niso(16) = 4
      prob(16, 1) = 0.02_wp
      prob(16, 2) = 0.76_wp
      prob(16, 3) = 4.29_wp
      prob(16, 4) = 94.93_wp
      massiso(16, 1) = 35.967080_wp
      massiso(16, 2) = 32.971458_wp
      massiso(16, 3) = 33.967867_wp
      massiso(16, 4) = 31.972071_wp

      ! 17 Cl (Chlorine)
      niso(17) = 2
      prob(17, 1) = 75.76_wp
      prob(17, 2) = 24.24_wp
      massiso(17, 1) = 34.968852_wp
      massiso(17, 2) = 36.965902_wp

      ! 18 Ar (Argon)
      niso(18) = 3
      prob(18, 1) = 0.3365_wp
      prob(18, 2) = 0.0632_wp
      prob(18, 3) = 99.6003_wp
      massiso(18, 1) = 35.9675462_wp
      massiso(18, 2) = 37.9627322_wp
      massiso(18, 3) = 39.96238312_wp

      ! 19 K  (Potassium)
      niso(19) = 3
      prob(19, 1) = 93.2581_wp
      prob(19, 2) = 0.0117_wp
      prob(19, 3) = 6.7302_wp
      massiso(19, 1) = 38.9637069_wp
      massiso(19, 2) = 39.96399867_wp
      massiso(19, 3) = 40.96182597_wp

      ! 20 Ca (Calcium)
      niso(20) = 5
      prob(20, 1) = 96.94_wp
      prob(20, 2) = 0.65_wp
      prob(20, 3) = 0.14_wp
      prob(20, 4) = 2.09_wp
      prob(20, 5) = 0.18_wp
      massiso(20, 1) = 39.962591_wp
      massiso(20, 2) = 41.958618_wp
      massiso(20, 3) = 42.958766_wp
      massiso(20, 4) = 43.955482_wp
      massiso(20, 5) = 45.953688_wp

      ! 21 Sc (Scandium)
      niso(21) = 1
      prob(21, 1) = 100.0_wp
      massiso(21, 1) = 44.9559102_wp

      ! 22 Ti (Titanium)
      niso(22) = 5
      prob(22, 1) = 8.0_wp
      prob(22, 2) = 7.3_wp
      prob(22, 3) = 73.8_wp
      prob(22, 4) = 5.5_wp
      prob(22, 5) = 5.4_wp
      massiso(22, 1) = 45.952632_wp
      massiso(22, 2) = 46.951763_wp
      massiso(22, 3) = 47.947946_wp
      massiso(22, 4) = 48.947870_wp
      massiso(22, 5) = 49.944791_wp

      ! 23 V  (Vanadium)
      niso(23) = 2
      prob(23, 1) = 0.250_wp
      prob(23, 2) = 99.750_wp
      massiso(23, 1) = 49.9471627_wp
      massiso(23, 2) = 50.943963_wp

      ! 24 Cr (Chromium)
      niso(24) = 4
      prob(24, 1) = 4.35_wp
      prob(24, 2) = 83.79_wp
      prob(24, 3) = 9.50_wp
      prob(24, 4) = 2.37_wp
      massiso(24, 1) = 49.946044_wp
      massiso(24, 2) = 51.940507_wp
      massiso(24, 3) = 52.940649_wp
      massiso(24, 4) = 53.938880_wp

      ! 25 Mn (Manganese)
      niso(25) = 1
      prob(25, 1) = 100.0_wp
      massiso(25, 1) = 54.938045_wp

      !  26 Fe (Iron)
      niso(26) = 4
      prob(26, 1) = 5.845_wp
      prob(26, 2) = 91.754_wp
      prob(26, 3) = 2.119_wp
      prob(26, 4) = 0.2_wp
      massiso(26, 1) = 53.939_wp
      massiso(26, 2) = 55.934_wp
      massiso(26, 3) = 56.935_wp
      massiso(26, 4) = 57.933_wp

      ! 27 Co (Cobalt)
      niso(27) = 1
      prob(27, 1) = 100.000_wp
      massiso(27, 1) = 58.933195_wp

      !  28 Ni (Nickel)
      niso(28) = 5
      prob(28, 1) = 68.08_wp
      prob(28, 2) = 26.22_wp
      prob(28, 3) = 1.14_wp
      prob(28, 4) = 3.63_wp
      prob(28, 5) = 0.93_wp
      massiso(28, 1) = 57.935343_wp
      massiso(28, 2) = 59.930786_wp
      massiso(28, 3) = 60.931056_wp
      massiso(28, 4) = 61.928345_wp
      massiso(28, 5) = 63.927966_wp

      ! 29 Cu (Copper)
      niso(29) = 2
      prob(29, 1) = 69.15_wp
      prob(29, 2) = 30.85_wp
      massiso(29, 1) = 62.929597_wp
      massiso(29, 2) = 64.927789_wp

      ! 30 Zn (Zinc)
      niso(30) = 5
      prob(30, 1) = 48.6_wp
      prob(30, 2) = 27.9_wp
      prob(30, 3) = 4.1_wp
      prob(30, 4) = 18.8_wp
      prob(30, 5) = 0.6_wp
      massiso(30, 1) = 63.929142_wp
      massiso(30, 2) = 65.926033_wp
      massiso(30, 3) = 66.927127_wp
      massiso(30, 4) = 67.924884_wp
      massiso(30, 5) = 69.925319_wp

      ! 31 Ga (Gallium)

      niso(31) = 2
      prob(31, 1) = 60.108_wp
      prob(31, 2) = 39.892_wp
      massiso(31, 1) = 68.925581_wp
      massiso(31, 2) = 70.9247073_wp

      ! 32 Ge (Germanium)
      niso(32) = 5
      prob(32, 1) = 20.38_wp
      prob(32, 2) = 27.31_wp
      prob(32, 3) = 7.76_wp
      prob(32, 4) = 36.72_wp
      prob(32, 5) = 7.83_wp
      massiso(32, 1) = 69.924247_wp
      massiso(32, 2) = 71.922076_wp
      massiso(32, 3) = 72.923459_wp
      massiso(32, 4) = 73.921178_wp
      massiso(32, 5) = 75.921402_wp

      ! 33 As (Arsenic)
      niso(33) = 1
      prob(33, 1) = 100.0_wp
      massiso(33, 1) = 74.921596_wp

      ! 34 Se (Selenium)
      niso(34) = 6
      prob(34, 1) = 0.89_wp
      prob(34, 2) = 9.37_wp
      prob(34, 3) = 7.63_wp
      prob(34, 4) = 23.77_wp
      prob(34, 5) = 49.61_wp
      prob(34, 6) = 8.73_wp
      massiso(34, 1) = 73.922476_wp
      massiso(34, 2) = 75.919213_wp
      massiso(34, 3) = 76.919914_wp
      massiso(34, 4) = 77.917309_wp
      massiso(34, 5) = 79.916521_wp
      massiso(34, 6) = 81.916699_wp

      ! 35 Br (Bromine)
      niso(35) = 2
      prob(35, 1) = 50.69_wp
      prob(35, 2) = 49.31_wp
      massiso(35, 1) = 78.91833_wp
      massiso(35, 2) = 80.91629_wp

      ! 36 Kr (Krypton)
      niso(36) = 6
      prob(36, 1) = 0.355_wp
      prob(36, 2) = 2.286_wp
      prob(36, 3) = 11.593_wp
      prob(36, 4) = 11.500_wp
      prob(36, 5) = 56.987_wp
      prob(36, 6) = 17.279_wp
      massiso(36, 1) = 77.920388_wp
      massiso(36, 2) = 79.916379_wp
      massiso(36, 3) = 81.9134850_wp
      massiso(36, 4) = 82.914137_wp
      massiso(36, 5) = 83.911508_wp
      massiso(36, 6) = 85.910615_wp

      ! 37 Rb (Rubidium)
      niso(37) = 2
      prob(37, 1) = 72.17_wp
      prob(37, 2) = 27.83_wp
      massiso(37, 1) = 84.9117924_wp
      massiso(37, 2) = 86.9091858_wp

      ! 38 Sr (Strontium)
      niso(38) = 4
      prob(38, 1) = 0.56_wp
      prob(38, 2) = 9.86_wp
      prob(38, 3) = 7.00_wp
      prob(38, 4) = 82.58_wp
      massiso(38, 1) = 83.913426_wp
      massiso(38, 2) = 85.9092647_wp
      massiso(38, 3) = 86.908816_wp
      massiso(38, 4) = 87.9056167_wp

      ! 39 Y  (Yttrium)

      niso(39) = 1
      prob(39, 1) = 100.0_wp
      massiso(39, 1) = 88.9058485_wp

      ! 40 Zr (Zirconium)
      niso(40) = 5
      prob(40, 1) = 51.45_wp
      prob(40, 2) = 11.22_wp
      prob(40, 3) = 17.15_wp
      prob(40, 4) = 17.38_wp
      prob(40, 5) = 2.80_wp
      massiso(40, 1) = 89.904704_wp
      massiso(40, 2) = 90.905646_wp
      massiso(40, 3) = 91.905040_wp
      massiso(40, 4) = 93.906315_wp
      massiso(40, 5) = 95.908273_wp

      ! 41 Nb (Niobium)
      niso(41) = 1
      prob(41, 1) = 100.0_wp
      massiso(41, 1) = 92.9063762_wp

      ! 42 Mo (Molybdenum)
      niso(42) = 7
      prob(42, 1) = 14.77_wp
      prob(42, 2) = 9.23_wp
      prob(42, 3) = 15.90_wp
      prob(42, 4) = 16.68_wp
      prob(42, 5) = 9.56_wp
      prob(42, 6) = 24.19_wp
      prob(42, 7) = 9.67_wp
      massiso(42, 1) = 91.906810_wp
      massiso(42, 2) = 93.9050876_wp
      massiso(42, 3) = 94.9058406_wp
      massiso(42, 4) = 95.9046780_wp
      massiso(42, 5) = 96.9060201_wp
      massiso(42, 6) = 97.9054069_wp
      massiso(42, 7) = 99.907476_wp

      ! 43 Tc (Technetium) no stable isotopes

      ! 44 Ru (Ruthenium)
      niso(44) = 7
      prob(44, 1) = 5.554_wp
      prob(44, 2) = 1.873_wp
      prob(44, 3) = 12.761_wp
      prob(44, 4) = 12.607_wp
      prob(44, 5) = 17.062_wp
      prob(44, 6) = 31.551_wp
      prob(44, 7) = 18.623_wp
      massiso(44, 1) = 95.907598_wp
      massiso(44, 2) = 97.905287_wp
      massiso(44, 3) = 98.905939_wp
      massiso(44, 4) = 99.904219_wp
      massiso(44, 5) = 100.905582_wp
      massiso(44, 6) = 101.906323_wp
      massiso(44, 7) = 103.905433_wp

      ! 45 Rh (Rhodium)

      niso(45) = 1
      prob(45, 1) = 100.0_wp
      massiso(45, 1) = 102.905504_wp

      ! 46 Pd (Palladium)
      niso(46) = 6
      prob(46, 1) = 1.02_wp
      prob(46, 2) = 11.15_wp
      prob(46, 3) = 22.34_wp
      prob(46, 4) = 27.33_wp
      prob(46, 5) = 26.47_wp
      prob(46, 6) = 11.73_wp
      massiso(46, 1) = 101.905609_wp
      massiso(46, 2) = 103.904036_wp
      massiso(46, 3) = 104.905085_wp
      massiso(46, 4) = 105.903486_wp
      massiso(46, 5) = 107.903892_wp
      massiso(46, 6) = 109.905153_wp

      ! 47 Ag (Silver)

      niso(47) = 2
      prob(47, 1) = 51.839_wp
      prob(47, 2) = 48.161_wp
      massiso(47, 1) = 106.905093_wp
      massiso(47, 2) = 108.904756_wp

      ! 48 Cd (Cadmium)

      niso(48) = 8
      prob(48, 1) = 1.25_wp
      prob(48, 2) = 0.89_wp
      prob(48, 3) = 12.49_wp
      prob(48, 4) = 12.80_wp
      prob(48, 5) = 24.13_wp
      prob(48, 6) = 12.22_wp
      prob(48, 7) = 28.73_wp
      prob(48, 8) = 7.49_wp
      massiso(48, 1) = 105.906458_wp
      massiso(48, 2) = 107.904183_wp
      massiso(48, 3) = 109.903006_wp
      massiso(48, 4) = 110.904182_wp
      massiso(48, 5) = 111.9027577_wp
      massiso(48, 6) = 112.9044014_wp
      massiso(48, 7) = 113.9033586_wp
      massiso(48, 8) = 115.904756_wp

      ! 49 In (Indium)
      niso(49) = 2
      prob(49, 1) = 4.29_wp
      prob(49, 2) = 95.71_wp
      massiso(49, 1) = 112.904062_wp
      massiso(49, 2) = 114.903879_wp

      ! 50 Sn (Tin)
      niso(50) = 10
      prob(50, 1) = 0.97_wp
      prob(50, 2) = 0.66_wp
      prob(50, 3) = 0.34_wp
      prob(50, 4) = 14.54_wp
      prob(50, 5) = 7.68_wp
      prob(50, 6) = 24.22_wp
      prob(50, 7) = 8.59_wp
      prob(50, 8) = 32.58_wp
      prob(50, 9) = 4.63_wp
      prob(50, 10) = 5.79_wp
      massiso(50, 1) = 111.904848_wp
      massiso(50, 2) = 113.902779_wp
      massiso(50, 3) = 114.903342_wp
      massiso(50, 4) = 115.901741_wp
      massiso(50, 5) = 116.902952_wp
      massiso(50, 6) = 117.901603_wp
      massiso(50, 7) = 118.903308_wp
      massiso(50, 8) = 119.902194_wp
      massiso(50, 9) = 121.903439_wp
      massiso(50, 10) = 123.905274_wp

      ! 51 Sb (Antimony)
      niso(51) = 2
      prob(51, 1) = 57.21_wp
      prob(51, 2) = 42.79_wp
      massiso(51, 1) = 120.903816_wp
      massiso(51, 2) = 122.904214_wp

      ! 52 Te (Tellurium)
      niso(52) = 8
      prob(52, 1) = 0.09_wp
      prob(52, 2) = 2.55_wp
      prob(52, 3) = 0.89_wp
      prob(52, 4) = 4.74_wp
      prob(52, 5) = 7.07_wp
      prob(52, 6) = 18.84_wp
      prob(52, 7) = 31.75_wp
      prob(52, 8) = 34.09_wp
      massiso(52, 1) = 119.90402_wp
      massiso(52, 2) = 121.903044_wp
      massiso(52, 3) = 122.904270_wp
      massiso(52, 4) = 123.902817_wp
      massiso(52, 5) = 124.904431_wp
      massiso(52, 6) = 125.903312_wp
      massiso(52, 7) = 127.904463_wp
      massiso(52, 8) = 129.906224_wp

      ! 53  I (Iodine)
      niso(53) = 1
      prob(53, 1) = 100.00_wp
      massiso(53, 1) = 126.904473_wp

      ! 54 Xe (Xenon)
      niso(54) = 9
      prob(54, 1) = 0.0952_wp
      prob(54, 2) = 0.0890_wp
      prob(54, 3) = 1.9102_wp
      prob(54, 4) = 26.4006_wp
      prob(54, 5) = 4.0710_wp
      prob(54, 6) = 21.2324_wp
      prob(54, 7) = 26.9086_wp
      prob(54, 8) = 10.4357_wp
      prob(54, 9) = 8.8573_wp
      massiso(54, 1) = 123.9058954_wp
      massiso(54, 2) = 125.904268_wp
      massiso(54, 3) = 127.9035305_wp
      massiso(54, 4) = 128.9047799_wp
      massiso(54, 5) = 129.9035089_wp
      massiso(54, 6) = 130.9050828_wp
      massiso(54, 7) = 131.9041546_wp
      massiso(54, 8) = 133.9053945_wp
      massiso(54, 9) = 135.907220_wp

      ! 55 Cs (Caesium)
      niso(55) = 1
      prob(55, 1) = 100.0_wp
      massiso(55, 1) = 132.905447_wp

      ! 56 Ba (Barium)
      niso(56) = 7
      prob(56, 1) = 0.106_wp
      prob(56, 2) = 0.101_wp
      prob(56, 3) = 2.417_wp
      prob(56, 4) = 6.592_wp
      prob(56, 5) = 7.854_wp
      prob(56, 6) = 11.232_wp
      prob(56, 7) = 71.698_wp
      massiso(56, 1) = 129.906311_wp
      massiso(56, 2) = 131.905056_wp
      massiso(56, 3) = 133.904504_wp
      massiso(56, 4) = 134.905684_wp
      massiso(56, 5) = 135.904571_wp
      massiso(56, 6) = 136.905822_wp
      massiso(56, 7) = 137.905242_wp

      ! 57 La (Lanthanum)
      niso(57) = 2
      prob(57, 1) = 0.09_wp
      prob(57, 2) = 99.91_wp
      massiso(57, 1) = 137.907112_wp
      massiso(57, 2) = 138.909477_wp

      ! 58 Ce (Cerium)
      niso(58) = 4
      prob(58, 1) = 0.185_wp
      prob(58, 2) = 0.251_wp
      prob(58, 3) = 88.450_wp
      prob(58, 4) = 11.114_wp
      massiso(58, 1) = 135.907140_wp
      massiso(58, 2) = 137.905986_wp
      massiso(58, 3) = 139.905435_wp
      massiso(58, 4) = 141.909241_wp

      ! 59 Pr (Praseodymium)
      niso(59) = 1
      prob(59, 1) = 100.0_wp
      massiso(59, 1) = 140.907648_wp

      ! 60 Nd (Neodymium)

      niso(60) = 7
      prob(60, 1) = 27.2_wp
      prob(60, 2) = 12.2_wp
      prob(60, 3) = 23.8_wp
      prob(60, 4) = 8.3_wp
      prob(60, 5) = 17.2_wp
      prob(60, 6) = 5.7_wp
      prob(60, 7) = 5.6_wp
      massiso(60, 1) = 141.907719_wp
      massiso(60, 2) = 142.909810_wp
      massiso(60, 3) = 143.910083_wp
      massiso(60, 4) = 144.912569_wp
      massiso(60, 5) = 145.913113_wp
      massiso(60, 6) = 147.916889_wp
      massiso(60, 7) = 149.920887_wp

      ! 61 Pm (Promethium) no stable isotopes

      ! 62 Sm (Samarium)
      niso(62) = 7
      prob(62, 1) = 3.07_wp
      prob(62, 2) = 14.99_wp
      prob(62, 3) = 11.24_wp
      prob(62, 4) = 13.82_wp
      prob(62, 5) = 7.38_wp
      prob(62, 6) = 26.75_wp
      prob(62, 7) = 22.75_wp
      massiso(62, 1) = 143.911996_wp
      massiso(62, 2) = 146.914894_wp
      massiso(62, 3) = 147.914818_wp
      massiso(62, 4) = 148.917180_wp
      massiso(62, 5) = 149.917272_wp
      massiso(62, 6) = 151.919729_wp
      massiso(62, 7) = 153.922206_wp

      ! 63 Eu (Europium)
      niso(63) = 2
      prob(63, 1) = 47.81_wp
      prob(63, 2) = 52.19_wp
      massiso(63, 1) = 150.919846_wp
      massiso(63, 2) = 152.921227_wp

      ! 64 Gd (Gadolinium)
      niso(64) = 7
      prob(64, 1) = 0.20_wp
      prob(64, 2) = 2.18_wp
      prob(64, 3) = 14.80_wp
      prob(64, 4) = 20.47_wp
      prob(64, 5) = 15.65_wp
      prob(64, 6) = 24.84_wp
      prob(64, 7) = 21.86_wp
      massiso(64, 1) = 151.919789_wp
      massiso(64, 2) = 153.920862_wp
      massiso(64, 3) = 154.922619_wp
      massiso(64, 4) = 155.922120_wp
      massiso(64, 5) = 156.923957_wp
      massiso(64, 6) = 157.924101_wp
      massiso(64, 7) = 159.927051_wp

      ! 65 Tb (Terbium)
      niso(65) = 1
      prob(65, 1) = 100.0_wp
      massiso(65, 1) = 158.925343_wp

      ! 66 Dy (Dysprosium)

      niso(66) = 7
      prob(66, 1) = 0.056_wp
      prob(66, 2) = 0.095_wp
      prob(66, 3) = 2.34_wp
      prob(66, 4) = 18.889_wp
      prob(66, 5) = 25.475_wp
      prob(66, 6) = 24.896_wp
      prob(66, 7) = 28.260_wp
      massiso(66, 1) = 155.924278_wp
      massiso(66, 2) = 157.924405_wp
      massiso(66, 3) = 159.925194_wp
      massiso(66, 4) = 160.926930_wp
      massiso(66, 5) = 161.926795_wp
      massiso(66, 6) = 162.928728_wp
      massiso(66, 7) = 163.929171_wp

      ! 67 Ho (Holmium)
      niso(67) = 1
      prob(67, 1) = 100.0_wp
      massiso(67, 1) = 164.930319_wp

      ! 68 Er (Erbium)
      niso(68) = 6
      prob(68, 1) = 0.139_wp
      prob(68, 2) = 1.601_wp
      prob(68, 3) = 33.503_wp
      prob(68, 4) = 22.869_wp
      prob(68, 5) = 26.978_wp
      prob(68, 6) = 14.910_wp
      massiso(68, 1) = 161.928775_wp
      massiso(68, 2) = 163.929197_wp
      massiso(68, 3) = 165.930290_wp
      massiso(68, 4) = 166.932046_wp
      massiso(68, 5) = 167.932368_wp
      massiso(68, 6) = 169.935461_wp
      ! 69 Tm (Thulium)
      niso(69) = 1
      prob(69, 1) = 100.0_wp
      massiso(69, 1) = 168.934211_wp

      ! 70 Yb (Ytterbium)
      niso(70) = 7
      prob(70, 1) = 0.13_wp
      prob(70, 2) = 3.04_wp
      prob(70, 3) = 14.28_wp
      prob(70, 4) = 21.83_wp
      prob(70, 5) = 16.13_wp
      prob(70, 6) = 31.83_wp
      prob(70, 7) = 12.76_wp
      massiso(70, 1) = 167.933895_wp
      massiso(70, 2) = 169.934759_wp
      massiso(70, 3) = 170.936323_wp
      massiso(70, 4) = 171.936378_wp
      massiso(70, 5) = 172.938207_wp
      massiso(70, 6) = 173.938858_wp
      massiso(70, 7) = 175.942569_wp

      ! 71 Lu (Lutetium)
      niso(71) = 2
      prob(71, 1) = 97.41_wp
      prob(71, 2) = 2.59_wp
      massiso(71, 1) = 174.9407682_wp
      massiso(71, 2) = 175.9426827_wp

      ! 72 Hf (Hafnium)
      niso(72) = 6
      prob(72, 1) = 0.16_wp
      prob(72, 2) = 5.26_wp
      prob(72, 3) = 18.60_wp
      prob(72, 4) = 27.28_wp
      prob(72, 5) = 13.62_wp
      prob(72, 6) = 35.08_wp
      massiso(72, 1) = 173.940042_wp
      massiso(72, 2) = 175.941403_wp
      massiso(72, 3) = 176.943204_wp
      massiso(72, 4) = 177.9436981_wp
      massiso(72, 5) = 178.9458154_wp
      massiso(72, 6) = 179.9465488_wp

      ! 73 Ta (Tantalum)
      niso(73) = 2
      prob(73, 1) = 0.012_wp
      prob(73, 2) = 99.988_wp
      massiso(73, 1) = 179.947466_wp
      massiso(73, 2) = 180.947996_wp

      ! 74 W (Tungsten)
      niso(74) = 5
      prob(74, 1) = 0.12_wp
      prob(74, 2) = 26.50_wp
      prob(74, 3) = 14.31_wp
      prob(74, 4) = 30.64_wp
      prob(74, 5) = 28.43_wp
      massiso(74, 1) = 179.946704_wp
      massiso(74, 2) = 181.948204_wp
      massiso(74, 3) = 182.950223_wp
      massiso(74, 4) = 183.950931_wp
      massiso(74, 5) = 185.954364_wp

      ! 75 Re (Rhenium)
      niso(75) = 2
      prob(75, 1) = 37.40_wp
      prob(75, 2) = 62.60_wp
      massiso(75, 1) = 184.952955_wp
      massiso(75, 2) = 186.955750_wp

      ! 76 Os (Osmium)
      niso(76) = 7
      prob(76, 1) = 0.02_wp
      prob(76, 2) = 1.59_wp
      prob(76, 3) = 1.96_wp
      prob(76, 4) = 13.24_wp
      prob(76, 5) = 16.15_wp
      prob(76, 6) = 26.26_wp
      prob(76, 7) = 40.78_wp
      massiso(76, 1) = 183.952491_wp
      massiso(76, 2) = 185.953838_wp
      massiso(76, 3) = 186.9557476_wp
      massiso(76, 4) = 187.9558357_wp
      massiso(76, 5) = 188.958145_wp
      massiso(76, 6) = 189.958445_wp
      massiso(76, 7) = 191.961479_wp

      ! 77 Ir (Iridium)
      niso(77) = 2
      prob(77, 1) = 37.3_wp
      prob(77, 2) = 62.7_wp
      massiso(77, 1) = 190.960591_wp
      massiso(77, 2) = 192.96293_wp

      ! 78 Pt (Platinum)
      niso(78) = 5
      prob(78, 1) = 0.782_wp
      prob(78, 2) = 32.97_wp
      prob(78, 3) = 33.83_wp
      prob(78, 4) = 25.24_wp
      prob(78, 5) = 7.16_wp
      massiso(78, 1) = 191.961038_wp
      massiso(78, 2) = 193.962680_wp
      massiso(78, 3) = 194.964791_wp
      massiso(78, 4) = 195.964951_wp
      massiso(78, 5) = 197.967893_wp

      ! 79 Au (Gold)
      niso(79) = 1
      prob(79, 1) = 100.0_wp
      massiso(79, 1) = 196.966551_wp

      ! 80 Hg (Mercury)
      niso(80) = 7
      prob(80, 1) = 0.15_wp
      prob(80, 2) = 10.0_wp
      prob(80, 3) = 16.87_wp
      prob(80, 4) = 23.10_wp
      prob(80, 5) = 13.19_wp
      prob(80, 6) = 29.86_wp
      prob(80, 7) = 6.87_wp
      massiso(80, 1) = 195.965833_wp
      massiso(80, 2) = 197.966769_wp
      massiso(80, 3) = 198.968279_wp
      massiso(80, 4) = 199.968326_wp
      massiso(80, 5) = 200.970302_wp
      massiso(80, 6) = 201.970643_wp
      massiso(80, 7) = 203.973494_wp

      ! 81 Tl (Thallium)
      niso(81) = 2
      prob(81, 1) = 29.52_wp
      prob(81, 2) = 70.48_wp
      massiso(81, 1) = 202.972344_wp
      massiso(81, 2) = 204.974427_wp

      ! 82 Pb (Lead)
      niso(82) = 4
      prob(82, 1) = 1.3_wp
      prob(82, 2) = 24.1_wp
      prob(82, 3) = 22.1_wp
      prob(82, 4) = 52.4_wp
      massiso(82, 1) = 203.973043_wp
      massiso(82, 2) = 205.974465_wp
      massiso(82, 3) = 206.975897_wp
      massiso(82, 4) = 207.976652_wp

      ! 83 Bi (Bismuth)
      niso(83) = 1
      prob(83, 1) = 100.0_wp
      massiso(83, 1) = 208.980398_wp

      ! 84 Po (Polonium) no stable isotopes

      ! 85 At (Astatine) no stable isotopes

      ! 86 Rn (Radon) no stable isotopes

      prob = prob*0.01_wp
      list_masses = 0.0_wp
      index_mass = 0
      loop = 0
      store_int = 0.0_wp

!TODO FIXME add here checks later
      molmasses = 0

      call rdshort_real("pfrag", pfrag)

      call read_structure(mol, fname, error, filetype%xyz)
      nat = mol%nat
      allocate (iat(nat))
      iat = mol%num(mol%id)

! loop over random number runs
      do n = 1, nrnd
         xmass = 0
         do i = 1, nat
            iti = iat(i)
            r = rnd(n, i)
            p1 = 0.0_wp
            p2 = prob(iti, 1)
            do iso = 1, niso(iti)
               if (r >= p1 .and. r <= p2) then
                  x = massiso(iti, iso)
                  exit
               end if
               p1 = p2
               p2 = p2 + prob(iti, iso + 1)
            end do
            xmass = xmass + x
         end do
         z_chrg = 1.0_wp ! TODO ADAPT ME
         current_mass = xmass/float(abs(z_chrg))

         there = .true.
         !> only values larger than user input
         if (current_mass > env%mthr) then ! warning this can lead to strange results ! TODO FIXME make nicer

            !> loop over all entries in the list
            inner: do
               loop = loop + 1

               !>> write if there is no entry
               if (list_masses(loop) == 0.0_wp) then
                  there = .false.
                  exit inner
                  !>> true if already in list, end
               elseif (abs(list_masses(loop) - current_mass) < 1.0d-10) then
                  there = .true.
                  store_int(loop) = store_int(loop) + 1
                  exit inner
                  !>> false if not in list, store
               elseif (abs(list_masses(loop) - current_mass) > 1.0d-10) then
                  there = .false.
                  cycle inner
               end if
            end do inner
            loop = 0

            !> if it is not in the list, add it
            if (.not. there) then
               index_mass = index_mass + 1
               list_masses(index_mass) = current_mass
               store_int(index_mass) = store_int(index_mass) + 1
               !write(*,*) list_masses(index_mass)
            end if
         end if

      end do

      allocate (exact_intensity(index_mass))
      allocate (isotope_masses(index_mass))

      do loop = 1, index_mass
         isotope_masses(loop) = list_masses(loop)
         exact_intensity(loop) = float(store_int(loop))/sum(float(store_int))*pfrag
      end do

      ! if no isotopes we take here the peak of the isotope pattern with the highest intensity
      ! still not sure if this is right .... ! TODO FIXME
      if (env%noiso) then
         iipmax = maxval(exact_intensity)
         ! iipmax = sum(exact_intensity) ! alternative would be to take the sum of all intensities/integral over all peaks
         indexipmax = maxloc(exact_intensity, dim=1)
         mipmax = isotope_masses(indexipmax)
         deallocate (exact_intensity, isotope_masses)
         index_mass = 1
         allocate (exact_intensity(1), isotope_masses(1))
         exact_intensity(1) = iipmax
         isotope_masses(1) = mipmax
      end if
      ! because of mthr this can habben
      if (sum(store_int) .ne. 50000) then
         exact_intensity = 0.0_wp
         isotope_masses = 0.0_wp
      end if

   end subroutine getmass

   subroutine getmolmass(nat, iat, molmass)
      implicit none

      integer, intent(in) :: nat
      integer, intent(in) :: iat(nat)
      real(wp), intent(out) :: molmass
      integer :: i

      molmass = 0.0_wp
      do i = 1, nat
         molmass = molmass + ams(iat(i))
      end do
   end subroutine getmolmass
end module isotope_pattern
