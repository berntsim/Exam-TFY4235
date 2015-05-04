PROGRAM main
use parameters
use fractals
use openmod
IMPLICIT NONE
        integer :: i
        !call makePoints()
        call openfile("fracTest.dat",test)
        call getVariables()
        call generateOneDimFractal()
        call KochIsland()
        print*, size(FractalArray)
        call sneakPoints()
!        do i = 1,Nfrac*4
!                write(test,*) REAL(IslandArray(i)), IMAG(IslandArray(i))
        do i =1,8*Nfrac
                write(test,*) REAL(ExtendedArray(i)), IMAG(extendedArray(i))
        end do
        call writeGrid()
!        call makeGrid()

END PROGRAM
