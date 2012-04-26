CHARMM Element source/fcm/clcg.fcm $Revision: 1.3 $
C------------------------------------------------------------------CC
C  Header file for RNG: CLCG (source/util/clcg.src) X. Qian 12/00   C
C  The variables are:                                               C
C      Maxgen     512 (the same as MAXNODE in parallel.fcm)         C
C                 maximum number of independent streams of random   C
C                 number sequences  (This value can be increased    C
C                 as needed.)                                       C
C      IniSD      1                                                 C
C      LstSD      2                                                 C
C      NewSD      3                                                 C
C                 These three option flags are used to initialize   C
C                 the RNG with different initial conditions.        C
C                 By default, initial seeds lcgIg{(i=1, ..., 4),g}  C
C                 and last seeds  lcgLg{(i=1, ..., 4),g}            C
C                 (for g=1...Maxgen) are set to the original        C
C                 seeds, and previous seeds, respectively.          C
C                 Calls to IniGen(g,type)                           C
C                 (where type is IniSD or LstSD or NewSD)           C
C                 can reset the seeds for stream g to the initial   C
C                 values (type = IniSD), previous values            C
C                 (type = LstSD) or  new values (type = NewSD).     C
C      lcgIg      Initial seed values, dimension 4 by Maxgen        C
C                 for the four LCGs.                                C
C      lcgLg      Last seed values, dimension 4 by Maxgen           C
C      lcgCg      Current seed values, dimension 4 by Maxgen        C
C                                                                   C
C     lcgmul(4)   The multipliers for the 4 LCGs                    C
C     lcgmod(4)   The moduli for the 4 LCGs                         C
C                 THE MULTIPLIER AND MODULI VALUES                  C
C                     MUST NOT BE CHANGED                           C
C                                                                   C
C    lcgaw(4)     lcgmul{j}^{2^w}      w=41, j=1, ..., 4.           C
C    lcgavw(4)    lcgmul{j}^{2^(v+w)}, v=31, j=1, ..., 4.           C
C                 These two arrays are used to generate initial     C
C                 seeds for the specified Maxgen number of the      C
C                 streams  with the  initial seeds given by user    C
C                 or from default values.                           C
C                                                                   C
            integer IniSD,LstSD,NewSD,Maxgen 
            parameter (Maxgen = 512, IniSD = 1, LstSD = 2, NewSD = 3)
            integer lcgmul(4),lcgmod(4),lcgaw(4),lcgavw(4)
            integer lcgIg(4,Maxgen),lcgLg(4,Maxgen),lcgCg(4,Maxgen)
            data lcgmul/45991,207707,138556,49689/
            data lcgmod/2147483647,2147483543,2147483423,2147483323/ 
            common /CLCGRND/lcgIg,lcgLg,lcgCg,lcgaw,lcgavw 
C
