CHARMM Element source/util/clcg.src $Revision: 1.6 $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               Combined Linear Congruential Generator (CLCG)        C
C  Adapted from Pierre L'Ecuyer & Terry H Andres' C version code.    C
C  References                                                        C
C [1] P. L'Ecuyer and T. H. Andres,                                  C
C     ``A Random Number Generator Based on the Combination           C
C     of Four LCGs'', Mathematics and Computers in Simulation,       C
C     44 (1997), 99--107.                                            C
C [2] http://www.iro.umontreal.ca/~lecuyer/                          C
C                                                                    C
C For further information, please contact                            C
C              Tamar Schlick                                         C  
C              schlick@nyu.edu                                       C
C Converted to FORTRAN by                                            C
C               Xiaoliang Qian  10/7/99                              C
C                                                                    C
C     Fixed the usage of this code by                                C
C               Milan Hodoscek  6/6/04                               C
C                                                                    C
C     NOTE: The code needs to be improved for parallel               C
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  This pseudorandom number generator combines 4 linear              C
C  congruential generators (LCGs) to get the long period             C
C  of about 2^121 for the resulting sequence. This sequence has      C
C  passed many statistical tests of `randomness'.                    C
C  Essentially, the four LCGs defined as                             C
C                                                                    C
C      X{j,n} = lcgmul{j}X{j,n-1} mod lcgmod{j}                (1)   C
C                            j=1,2,3,4                               C
C                            n=1,2,3, ...                            C
C                                                                    C
C  with lcgmod{1}= 2147483647, lcgmod{2}=2147483543,                 C
C       lcgmod{3}=2147483423, lcgmod{4}=2147483323                   C
C  and lcgmul{1}=45991, lcgmul{2}=207707, lcgmul{3}=138556,          C
C  lcgmul{4}=49689.                                                  C
C                                                                    C
C  The construct                                                     C
C      Z{n} = (Sum [(-1)^{j+1}*X{j,n}/lcgmod{j}]) mod 1   (2)        C
C                                                                    C
C for n=1,2,3, ... is then a uniformly distributed random sequence   C
C in (0,1). It can be proved that the LCG corresponding to the       C
C combined generator has modulus, multiplier, and period length of   C
C      21267641435849934371830464348413044909,                       C
C      5494569482908719143153333426731027229,                        C
C      (2^{31}-2)(2^{31}-106)(2^{31}-226)(2^{31}-326) ~ 2^{121},     C
C  respectively.                                                     C
C                                                                    C
C  The default initial seed is the vector {11111111, 22222222,       C
C  33333333, 44444444} and can be changed by calling SetIniSD        C
C  after calling CLCGInit.                                           C
C                                                                    C
C  This RNG can be used under parallel conditions to give            C 
C  independent random number sequence when each processor            C
C  calls with different stream number g (e.g., RANDOM(g)).           C
C                                                                    C
C  To use these RNG routines, the user should proceed as follows:    C
C                                                                    C
C  1. Call routine CLCGInit() to initialize all Maxgen (100) streams C
C     using four default initial seeds.                              C
C  2. [Optional] Call SetiniSD(sd) with desired seed array           C
C     (4 values) to override the default values specified in Init(). C
C  3. Call function RANDOM(k), where k is an integer from 1 to 100   C
C     specifying the stream number. For parallel codes, k can be set C
C     to a processor id number.                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
           subroutine CLCGInit (sd)
C------------------------------------------------------------------------C
C   Initialize the RNG with seed values in vector sd of dimension 4.     C
C   Each such initial vector generates one stream of random variates     C
C   combining the 4 LCG sequences with resulting sequence length         C
C   2^{121}.                                                             C
C   The lcgaw and lcgavw arrays of dimension 4 have default values       C
C       lcgaw{j} = lcgmul{j}^{2^31} mod lcgmod{j}                        C
C   and lcgavw{j} = lcgmul{j}^{2^41} mod lcgmod{j},                      C
C   for j=1, ..., 4 corresponding to the 4 LCGs.                         C
C Converted to FORTRAN by                                                C
C               Xiaoliang Qian  10/7/99                                  C

           implicit none
           include 'clcg.h'
C
           integer MulMod
C
           integer sd(4)
           integer v,w,j,i
           parameter (v = 31, w = 41)
C
           do 100 j = 1,4
             lcgaw(j) = lcgmul(j)
                do 85 i = 1,w
                 lcgaw(j) = MulMod (lcgaw(j),lcgaw(j),lcgmod(j))
85              continue
             lcgavw(j) = lcgaw(j)
                do 90 i = 1,v
                 lcgavw(j) = MulMod (lcgavw(j),lcgavw(j),lcgmod(j))
90              continue
100        continue   
           call SetiniSD (sd)
           end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
          subroutine SetiniSD (s)
C------------------------------------------------------------------------C
C  Set initial seed values for all 100 (= Maxgen, defined in clcg.fcm)   C
C  streams using the initial seeds from the first stream.                C
C Converted to FORTRAN by                                                C
C               Xiaoliang Qian  10/7/99                                  C

          implicit none
C
          integer g,s(4),j 
          integer MulMod

          include 'clcg.h'

          do 70 j = 1,4
             lcgIg(j,1) = s(j)
 70       continue
          call IniGen (1,IniSD)
          do 80 g = 2, Maxgen
             do 75  j = 1,4
                lcgIg(j,g) = MulMod (lcgavw(j),lcgIg(j,g-1),lcgmod(j))
 75          continue
             call IniGen (g,IniSD)
 80       continue
          end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL*8 function RANDOM (g)
C------------------------------------------------------------------------C
C   Return a double precision uniformly distributed random number in     C
C   (0,1) from the gth stream and reset the current seed Cg accordingly  C
C   (i.e., using one of the 100 initial seed vectors generated in the    C
C   SetiniSD routine).                                                   C
C Converted to FORTRAN by                                                C
C               Xiaoliang Qian  10/7/99                                  C
C
C
      implicit none
      include 'clcg.h'

      integer g,k,s,j
      real*8 u(4)
      integer  dv(4),mv(4)
      data dv/46693,10339,15499,43218/ 
      data mv/25884, 870,3979,24121/ 
      data u/4.65661287524579692d-10,-4.65661310075985993d-10,
     &     4.65661336096842131d-10,-4.65661357780891134d-10/
C
      real*8 zero, one
      parameter (zero=0d0, one=1d0)
C
      if (g .le. 1)  g = 1
      g= mod(g-1,Maxgen) + 1
      RANDOM = zero
      do 110 j = 1,4
         s = lcgCg(j,g)
         k = s/dv(j)
         s = lcgmul(j) * (s - k * dv(j)) - k * mv(j)
         if (s .lt. 0) s = s + lcgmod(j)
         lcgCg(j,g) = s
         RANDOM = RANDOM + u(j) * s
         if (RANDOM .lt. zero)  RANDOM = RANDOM + one
         if (RANDOM .ge. one)   RANDOM = RANDOM - one
 110  continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
          subroutine  SetSeed (g,s)
C------------------------------------------------------------------------C
C  This optional routine uses the input seed value s for stream g        C
C  instead of the default settings (routine SetiniSD).                   C
C Converted to FORTRAN by                                                C
C               Xiaoliang Qian  10/7/99                                  C


          integer g,s(4),j

          include 'clcg.h'
              if (g .le. 1)  g = 1
               g= mod(g-1,Maxgen) + 1

          do 50 j = 1,4          
           lcgIg(j,g) = s(j)
50         continue
            call IniGen (g,IniSD)               
          end

          subroutine  GetSeed (g,s)
C------------------------------------------------------------------------C
C  This optional routine returns current seed value s for stream g       C
C Converted to FORTRAN by                                                C
C               Xiaoliang Qian  10/7/99                                  C


          integer g,s(4),j

          include 'clcg.h'
          if (g .le. 1)  g = 1
          g= mod(g-1,Maxgen) + 1
          
          do 50 j = 1,4          
             s(j)= lcgCg(j,g)
 50       continue
          return
          end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          subroutine IniGen (g,type)
C------------------------------------------------------------------------C
C  This optional routine resets the gth stream so that the initial seed  C
C  is either the original initial seed (if type = IniSD ) or the last    C
C   seed (if type = 3).                                                  C
C Converted to FORTRAN by                                                C
C               Xiaoliang Qian  10/7/99                                  C

          implicit none

                                                             
          integer g,type,j,ig
          integer MulMod
C
          include 'clcg.h'
C
          ig=g
          if (ig .le. 1)  ig = 1
          ig= mod(ig-1,Maxgen) + 1
          do 60 j = 1,4
             if (type .eq. IniSD) then
                lcgLg(j,ig) = lcgIg(j,ig)
             else
                if (type .eq. NewSD) 
     *             lcgLg(j,ig) = MulMod (lcgaw(j),lcgLg(j,ig),lcgmod(j)) 
             endif
             lcgCg(j,ig) = lcgLg(j,ig)
 60       continue
          return
          end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


         integer function MulMod (s,t,M)
C--------------------------------------------------------------------C
C  Return s*t mod M. All numbers out of range are truncated          C
C  before the mod operation to avoid overflow. See Park & Miller,    C
C  Comm. ACM 31, 1192 (1988) for this multiplication procedure.      C
C Converted to FORTRAN by                                            C
C               Xiaoliang Qian  10/7/99                              C
 
         integer s,t,M,H
         parameter (H = 32768) 
         integer   S0,S1,q,qh,rh,k 

         if (s .lt. 0) s = s + M
         if (t .lt. 0) t = t + M
         if (s .lt. H) then
            S0 = s
            MulMod = 0
         else
            S1 = s / H
            S0 = s - H * S1
            qh = M / H
            rh = M - H * qh

            if (S1 .ge. H) then
               S1 = S1 - H
               k = t / qh
               MulMod = H * (t - k * qh) - k * rh

 10            if (MulMod .lt. 0) then
                  MulMod = MulMod + M      
                  goto 10
               endif

            else
               MulMod = 0
            endif

            if (S1 .ne. 0) then
               q = M / S1
               k = t / q
               MulMod = MulMod - k * (M - S1 * q)
               if (MulMod .gt. 0) MulMod = MulMod - M
               MulMod = MulMod + S1 * (t - k * q) 

 20            if (MulMod .lt. 0) then
                  MulMod = MulMod + M      
                  goto 20
               endif
            endif

            k = MulMod / qh
            MulMod = H * (MulMod - k * qh) - k * rh
                 
 30         if (MulMod .lt. 0) then
               MulMod = MulMod + M      
               goto 30
            endif
         endif
         
         if (S0 .ne. 0) then
            Q = M / S0
            k = t / q                         
            MulMod = MulMod - k * (M - S0 * q)
            if (MulMod .gt. 0) MulMod = MulMod - M
            MulMod = MulMod + S0 * (t - k * q) 
 40         if (MulMod .lt. 0)  then
               MulMod = MulMod + M      
               goto 40
            endif
         endif
         return
         end
