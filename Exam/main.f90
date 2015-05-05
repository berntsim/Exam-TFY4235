PROGRAM main
use parameters
use fractals
use openmod
use F95_LAPACK, only : la_syev
IMPLICIT NONE
        integer :: i
        !call makePoints()
        call openfile("fracTest.dat",test)
        call getVariables()
        call generateOneDimFractal()
        call KochIsland()
        call sneakPoints()
        call writeFractal()
!        call writeGrid()
!        print*, Nreq*2 + 2
        call makeGrid()
        call iterateGrid()
        call generateAmatrix()
        allocate(U(Ncounter))
        call solveEVP()
        call plotMode()
        call plot_mode()
        call plot_fractal()
        do i = 1,10
                print*, U(i)
        end do
        do i = 1,10
                print*, sqrt(U(i))/delta
        end do
END PROGRAM
