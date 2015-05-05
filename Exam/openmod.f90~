MODULE openmod
IMPLICIT NONE
        INTEGER, PARAMETER      :: test = 1, test2 = 2, twodim = 3, mode = 3
        INTEGER                 :: res 
CONTAINS
        SUBROUTINE openfile(filename, filenr)
        CHARACTER(LEN = *), INTENT (IN)  :: filename
        INTEGER, INTENT (IN) :: filenr

                OPEN(UNIT = filenr, FILE=filename, IOSTAT = res)
                IF (res /= 0) THEN
                        PRINT *, 'Error in opening file,',filename,', status: ',res
                        STOP
                END IF
 !               PRINT *, filename, ' opened successfully!'
        END
END MODULE

