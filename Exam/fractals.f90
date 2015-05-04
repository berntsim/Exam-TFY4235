MODULE fractals
use parameters
use openmod
contains
SUBROUTINE makePoints()
        FractalArray(1) = CMPLX(0,0)
        FractalArray(Nfrac) = CMPLX(4,0)
        call iterateLine(FractalArray(1),FractalArray(Nfrac),FractalArray)  
END SUBROUTINE

SUBROUTINE iterateLine(ILposStart, ILposEnd, ILArray)
!Actually does the l-iterations. This function takes in start end end point of the line
!segment, to which one applies the quadratic Koch generator. 
        complex, intent(in)                     :: ILposStart, ILposEnd
        !start and end position of the line to iterate over
        complex, dimension(9), intent(inout)      :: ILArray
        !The complex array of the line segment
        complex, dimension(2)                     :: ILsub
        !The parts of the line which is moved in dimension
        logical                                 :: ILisX
        !Keeping track of change in x or y direction
        logical                                 :: ILisPosX
        !Checking if the change is positive in x first
        logical                                 :: ILisPosY
        !Checking if the change is positive in y first
        complex(wp)                             :: ILdiff
        
        !Setting up the initial requirements
        ILArray(1) = ILposStart
        ILArray(9) = ILposEnd
        ILisx = .false.
        ILisPosX = .false.
        ILisPosY = .false.
        ILdiff = ILposEnd- ILposStart

        !testing what kind of position we are at, and what type
        !of iteration is appropriate for this line segment.
        if (IMAG(ILdiff) > 0) then
                ILisX = .true.
        else if (IMAG(ILdiff) < 0) then
                ILisX = .true.
                ILisPosX = .true.
        else if (REAL(ILdiff) > 0) then
                ILisPosY = .true.
        end if
        
        !Initiating the points which doesn't change in the iteration 
        if (ILisX .and. ILisPosX) then
                ILArray(2) = CMPLX(0,-1) + ILposStart
                ILArray(8) = CMPLX(0,1) + ILposEnd
                ILArray(5) = CMPLX(0,-2) + ILposStart
        else if (ILisx) then
                ILArray(2) = CMPLX(0,1) + ILposStart
                ILArray(8) = CMPLX(0,-1) + ILposEnd
                ILArray(5) = CMPLX(0,2) + ILposStart
        else if (ILisPosY) then
                ILArray(2) = CMPLX(1,0) + ILposStart
                ILArray(8) = CMPLX(-1,0) + ILposEnd
                ILArray(5) = CMPLX(2,0) + ILposStart
        else
                ILArray(2) = CMPLX(-1,0) + ILposStart
                ILArray(8) = CMPLX(1,0) + ILposEnd
                ILArray(5) = CMPLX(-2,0) + ILposStart              
        end if
        
        !Making the first part of the Koch line, which we know the shape of
        call addBetween(ILisX, IlisPosX, ILisPosY, ILArray(2), ILarray(5), ILsub)
        ILArray(3) = ILsub(1)
        ILArray(4) = ILsub(2) 

        !We now do the second part, which is by definition in the opposite direction
        if (ILisPosX .or. ILisPosY) then
                ILisPosX = .false.
                ILisPosY = .false.
        else
                ILisPosX = .true.
                ILisPosY = .true.
        end if
        
        call addBetween(ILisX, ILisPosX, ILisPosY, ILArray(5), ILarray(8), ILsub)
        ILArray(6) = ILsub(1)
        ILArray(7) = ILsub(2) 


END SUBROUTINE


SUBROUTINE addBetween(ABisX, ABisPosX, ABisPosY, ABpos1, ABpos2, ABArray)
!This routine creates Koch curve in between the argument points
        complex, intent(in)                      :: ABpos1, ABpos2
        !The points defining the line to change
        logical, intent(in)                      :: ABisX, ABisPosX
        !The check to orient what direction the curve should be in
        logical, intent(in)                      :: ABisPosY
        complex, dimension(2), intent(inout)     :: ABArray
        !The array to which the line is updated
        
        !In the following test, the orientation and navigation
        !determines what values to change in order to obtain
        !the Koch Island fractal
        if (ABisX .and. ABisPosX) then
                ABArray(1) = CMPLX(1,0) + ABpos1 
                ABArray(2) = CMPLX(1.0) + ABpos2
        else if (ABisX) then
                ABArray(1) = CMPLX(-1,0) + ABpos1
                ABArray(2) = CMPLX(-1,0) + ABpos2
        else if (ABisPosY) then
                ABArray(1) = CMPLX(0,1) + ABpos1
                ABArray(2) = CMPLX(0,1) + ABpos2
        else
                ABArray(1) = CMPLX(0,-1) + ABpos1
                ABArray(2) = CMPLX(0,-1) + ABpos2
        end if
END SUBROUTINE

SUBROUTINE generateOneDimFractal()
!This is the main routine for the fractal module.
!In this routine, the fractal is generated for l_dim
!iterations, and returns a Koch fractal at level l_dim.
!The process is so that the points are shiftet 4 places
!to the right, and new points are introduced in between,
!which are Koch curves. The arrays are allocated and
!deallocated accordingly.        
        complex, dimension(:), allocatable      :: ODFtmp
        complex, dimension(9)                   :: ODFline
        integer                                 :: ODF_i, ODF_j, ODF_k
        !iterators over 
        integer                                 :: Nlevel

        !This doing the first curve generation, l1        
        allocate(FractalArray(9))
        FractalArray(1) = CMPLX(0,0)
        FractalArray(9) = CMPLX(4,0)
        call iterateLine(FractalArray(1),FractalArray(9),FractalArray)
        
        !We noe test to see if we just do one iteration, or several
        if (l_dim < 3) then
                return
        end if
        !If we only want to do one iteration, we return at this point        

        !If we want to do more iterations, we iterate over the amount
        !of levels we want, and create the line segments for each part.

        do l_level = 3,l_dim
                Nlevel = 8**(l_level-1)+1
                allocate(ODFtmp(Nlevel))
                !We now shift the points 4 places to the right
                do ODF_i = 1,8**(l_level-2)+1
                        ODFtmp(1+(ODF_i-1)*8) = 4*FractalArray(ODF_i)
                end do
                deallocate(FractalArray)
                !We create a new and larger array for inclusion of the
                !new points, in effect extending the array.
                allocate(FractalArray(Nlevel))
                FractalArray = ODFtmp
                !At this point, we create Koch curves on each of the
                !line segments
                do ODF_j = 1,8**(l_level-2)
                        call iterateLine(FractalArray(1+(ODF_j-1)*8),&
                             FractalArray(1+(ODF_j)*8),ODFline)
                        !Here we insert the new points into the 
                        !FractalArray
                        do ODF_k = 1,7
                                FractalArray(1+(ODF_j-1)*8+ODF_k) =&
                                ODFline(ODF_k + 1)
                        end do
                end do
                deallocate(ODFtmp)
        end do
END SUBROUTINE

SUBROUTINE KochIsland()
        integer :: Ki
        integer :: Ky_max
        
        !The point where the line segments needs to be connected
        Ky_max = 4**(l_dim-1)

        !Iterating over the FractalArray and placing the elements
        !into the islandArray. Note that we iterate backwards for 
        !the second and third FractalArrays in order to obtain a
        !seamless connection.
        do Ki = 1, Nfrac
                IslandArray(Ki) = FractalArray(Ki)
                IslandArray(Ki + 2*NFrac) = &
                        FractalArray(NFrac-Ki) + CMPLX(0,-Ky_max)
        end do
        FractalArray = FractalArray*CMPLX(0,-1)+CMPLX(Ky_max,0)
        do Ki = 1, Nfrac
                IslandArray(Ki + Nfrac) = FractalArray(Ki)
                IslandArray(Ki + 3*Nfrac) = FractalArray(Nfrac+1-Ki)&
                + CMPLX(-Ky_max,0) 
        end do

        !Shifting the fractal to the 1.st quadrant:
        IslandArray = IslandArray - CMPLX(GVxmin-1,GVymin-1)
END SUBROUTINE

SUBROUTINE getVariables()
!Getting the variables necessary for creating
!the grid with the fractal

        integer :: GVi
        !iterator

        Nreq = int(1._wp/delta)
        GVxmin = 0
        GVyreq = 0
        do GVi = 2,l_dim
                Nreq = Nreq + int(2/(delta*4**(GVi-1)))
                GVxmin = GVxmin - int(1/(delta*4**(GVi-1)))
                GVyreq = GVyreq + int(1/(delta*4**(GVi-1)))
        end do
        GVymin = GVyreq - Nreq
        GVxreq = Nreq +  GVxmin
        
        
!        print*, delta
!        print*, GVymin, GVyreq, GVxmin, GVxreq        
END SUBROUTINE

SUBROUTINE sneakPoints()
!The purpose of this subroutine is to extend the
!fractal array so that one has a precission greater
!than the smallest fractal block.
!Notice that we now get three extra points where
!we connect the fractal lines!
        integer :: SPi
        logical :: SPisX, SPisPosX, SPisY

        IslandArray = IslandArray*2
        !Shifting the corners to make room for more points
        !In reality, this is more of a "strecth" than shift
        SPisX = .false.
        SPisPosX = .false.
        SPisY = .false.
        do SPi = 1, 4*Nfrac
                if (REAL(IslandArray(SPi+1))-REAL(IslandArray(SPi))>0) then
                        SPisX = .true.
                        SPisPosX = .true.
                else if (REAL(IslandArray(SPi+1))&
                        -REAL(IslandArray(SPi))<0) then
                        SPisX = .true.
                else if (IMAG(IslandArray(SPi+1))&
                        -IMAG(IslandArray(SPi))>0) then
                        SPIsPosX = .true.
                        SPisY = .true.
                else if (IMAG(IslandArray(SPi+1))&
                        -IMAG(IslandArray(SPi))<0) then
                        SPisY = .true.
                end if

                ExtendedArray(SPi*2-1) = IslandArray(SPi)
                if (SPIsPosX .and. SPisX) then
                       ExtendedArray(SPi*2) = IslandArray(SPi)&
                                            + CMPLX(1,0)
                else if (SPisX) then
                       ExtendedArray(SPi*2) = IslandArray(SPi)&
                                            + CMPLX(-1,0)    
                else if (SPisPosX) then
                       ExtendedArray(SPi*2) = IslandArray(SPi)&
                                            + CMPLX(0,1)
                else if (SPisY) then
                       ExtendedArray(SPi*2) = IslandArray(SPi)&
                                            + CMPLX(0,-1)
                else
                       ExtendedArray(SPi*2) = IslandArray(SPi)     
                end if
                SPisX = .false.
                SPisY = .false.
                SPisPosX = .false.
        end do
        ExtendedArray(8*Nfrac) = IslandArray(4*Nfrac)
        
END SUBROUTINE

SUBROUTINE makeGrid()
!In this routine, the grid is created with the boundary,
!in such a way that one keeps track of what points are 
!inside, outside or on the boundary
        integer  :: MGi
        !iterator 
        integer :: MGx, MGy
        
        allocate(Grid(2*Nreq+2,2*Nreq+2))
        Grid = 0
        do MGi = 1,Nfrac*8
                MGx = int(REALPART(ExtendedArray(MGi)))
                MGy = int(IMAG(ExtendedArray(MGi)))
                Grid(MGx,MGy) = 1
        end do
!        do MGi = 1,Nreq
                
!        end do

END SUBROUTINE


SUBROUTINE writeGrid()
        integer :: WGi,WGj
        real(wp) :: WGx,WGy
        call openfile('fracBound.dat',test2)
        call makeGrid()
        WGx = 0
        WGy = delta
        do WGi = 1,2*Nreq
                WGx = WGx + delta
                do WGj = 1,2*Nreq
                        write(test2,*) WGx, WGy, Grid(WGi,WGj)
                        WGy = WGy + delta
                end do
                write(test2,*)
                WGy = delta
        end do
END SUBROUTINE

END MODULE




