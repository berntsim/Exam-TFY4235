PROGRAM main
use parameters
use fractals
use openmod
use F95_LAPACK, only : la_syev
IMPLICIT NONE
        integer :: i
        call openfile("fracTest.dat",test)
        call getVariables()
        call generateOneDimFractal()
        call KochIsland()
        call sneakPoints()
        call writeFractal()
        call iterateGrid()

!=======================================================
!        call generateAmatrix()
        allocate(U(Ncounter)) !Do not comment out
!        call solveEVP()
!        call plotMode()
!        call plot_mode()
!        call plot_fractal()

!=======================================================
!this section needs to be uncommented to do task 5
!and should be commented for task 3

        call generateBmatrix()
        call solveEVPB()
        call plotMode()
        call plot_biharmonic()
!=======================================================

        do i = 1,10
                print*, U(i)
        end do
        print*, 
        do i = 1,10
                print*, U(i)/(delta**4)
        end do
        print*,
        print*, delta
END PROGRAM
