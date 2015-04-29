*DECLARE DOSUMS1
*INSERT DOSUMS1.234
         if (isnan(SUM_RESULTS(K))) then
            write(6,*) 'DO_SUM error at ',K
            STOP 11
         endif
