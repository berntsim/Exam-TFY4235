MODULE fractals
use parameters
use openmod
use F95_LAPACK, only: LA_SYEVD
contains




SUBROUTINE iterateLine(ILposStart, ILposEnd, ILArray)

!==============================================================
!This subroutine does the l-iterations. It takes in start
!and end point of the linesegment, to which one applies the
!quadratic Koch generator. 
!==============================================================
        complex, intent(in)                     :: ILposStart, ILposEnd
        !start and end position of the line to iterate over
        complex, dimension(9), intent(inout)    :: ILArray
        !The complex array of the line segment
        complex, dimension(2)                   :: ILsub
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

!==================================================================
!This routine creates Koch curve in between the argument points
!==================================================================
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

!==========================================================
!This is the main routine for the fractal creation.
!In this routine, the fractal is generated for l_dim
!iterations, and returns a Koch fractal at level l_dim.
!The process is so that the points are shiftet 4 places
!to the right, and new points are introduced in between,
!which are Koch curves. The arrays are allocated and
!deallocated accordingly.
!==========================================================
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

!==============================================================
!This subrotuine connectes the line segment which have 
!been through the fractal generator, and makes the 
!"Koch Island", i.e. the square fractal.
!==============================================================
        integer :: Ki
        integer :: K_max
        
        !The point where the line segments needs to be connected
        K_max = 4**(l_dim-1)

        !Iterating over the FractalArray and placing the elements
        !into the islandArray. Note that we iterate backwards for 
        !the second and third FractalArrays in order to obtain a
        !seamless connection (keeping in mind how gnuplot inter-
        !prets the data).
        do Ki = 1, Nfrac
                IslandArray(Ki) = FractalArray(Ki)
                IslandArray(Ki + 2*NFrac) = &
                        FractalArray(NFrac-Ki) + CMPLX(0,-K_max)
        end do
        FractalArray = FractalArray*CMPLX(0,-1)+CMPLX(K_max,0)
        do Ki = 1, Nfrac
                IslandArray(Ki + Nfrac) = FractalArray(Ki)
                IslandArray(Ki + 3*Nfrac) = FractalArray(Nfrac+1-Ki)&
                + CMPLX(-K_max,0) 
        end do

        !Shifting the fractal to the 1.st quadrant:
        IslandArray = IslandArray - CMPLX(GVxmin-1,GVymin-1)
END SUBROUTINE

SUBROUTINE getVariables()

!=====================================================
!Getting the variables necessary for creating
!the grid with the fractal
!=====================================================
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

!=======================================================
!The purpose of this subroutine is to extend the
!fractal array so that one has a precission greater
!than the smallest fractal block.
!Notice that we now get three extra points where
!we connect the fractal lines!
!=======================================================
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
        delta = delta/2 
        !The delta must be updated since we now sneak in new
        !points
END SUBROUTINE

SUBROUTINE makeGrid()

!==========================================================
!This subroutine creates the grid  with the boundary,
!in such a way that one keeps track of what points are 
!inside, outside or on the boundary
!==========================================================
        integer  :: MGi, MGj
        !iterators 
        logical :: KeepTrack
        !keeps track on wether the boundary is crossed or not
        !in the previous step
        integer :: counter
        !Counts steps taken since the boundary is crossed. This
        !value will be 0 outside, 1 on the boundary and keep 
        !counting while inside.
        complex :: MGxy
        integer :: KL
        
        allocate(Grid(2*Nreq+2,2*Nreq+2))
        Grid = 0
        do MGi = 1,Nfrac*8
                MGx = int(REALPART(ExtendedArray(MGi)))
                MGy = int(IMAG(ExtendedArray(MGi)))
                Grid(MGx,MGy) = 1
        end do

        KeepTrack = .false.
        counter = 0
        do MGi = 1,2*Nreq+2
                do MGj = 1,2*Nreq+2
                        MGxy = CMPLX(MGi,MGj)
                        !We now test if we cross the boundary
                        !and update values to keep track.
                        !The ANY statement checks if MGxy exist
                        !in any element in ExtendedArray. In 
                        !other words, it checks if one is at
                        !the boundary or not.
                        if (ANY(ExtendedArray == MGxy) .and.&
                           (counter == 0)) then
                                KeepTrack = .true.
                                counter = 1
                        !We now check if we are not on the boundary, and
                        !if the counter is 1, i.e. that we croseed it at
                        !the previous point.
                        else if (.not. ANY(ExtendedArray .eq. MGxy) .and.&
                                (counter == 1)) then
                                !We now check the point to the left, since 
                                !we iterate to the right, and the left point
                                !will indicate if we are inside or not. If
                                !we are inside -> update grid point and 
                                !counter. If not -> we are outside, and
                                !we update values accordingly.
                                if (Grid(MGi-1,MGj) == 1) then
                                        Grid(MGi,MGj) = 1
                                        counter = counter + 1
                                else
                                        KeepTrack = .false.
                                        counter = 0
                                end if
                        !Here we check if we are not on the boundary, and 
                        !if we are inside. We then update the grid point 
                        !that is inside, and keep the counter running
                        else if (.not. ANY(ExtendedArray .eq. MGxy) .and.&
                                (KeepTrack)) then
                                Grid(MGi,MGj) = 1
                                counter = counter + 1
                        !The final test checks whether we are on the
                        !boundary, and whether the previous point 
                        !also was on the boundary.
                        else if (ANY(ExtendedArray .eq. MGxy) .and.&
                                 (counter /= 1) ) then
                                KeepTrack = .false.
                                counter = 0
                        end if
                end do
        end do

END SUBROUTINE

SUBROUTINE writeFractal()

!===========================================================
!This subroutine is used for generating plots of the
!fractal shape. This is in 2d, since we have not yet
!put it on the grid, and updated inside/outside values.
!===========================================================

        call openfile('3dfractal.dat',twodim)
        do i =1,8*Nfrac
                write(twodim,*) REAL(ExtendedArray(i))*delta,&
                      IMAG(extendedArray(i))*delta, 0
        end do 
END SUBROUTINE

SUBROUTINE plotFractal()

!============================================================
!This subroutine writes the 2d-fractal to a file.
!============================================================
        call openfile('2dfractal.dat',twodim)
        do i =1,8*Nfrac
                write(twodim,*) REAL(ExtendedArray(i))*delta,&
                      IMAG(extendedArray(i))*delta
        end do 
END SUBROUTINE


SUBROUTINE writeGrid()

!===========================================================
!This sibroutine writes the data from the grid
!to a file.
!===========================================================
        integer :: WGi,WGj
        real(wp) :: WGx,WGy

        call openfile('fracBound.dat',test2)
        WGx = 0
        WGy = delta
        do WGi = 1,2*Nreq+2
                WGx = WGx + delta
                do WGj = 1,2*Nreq+2
                        write(test2,*) WGx, WGy, Grid(WGi,WGj)
                        WGy = WGy + delta
                end do
                write(test2,*)
                WGy = delta
        end do
END SUBROUTINE

SUBROUTINE iterateGrid()

!========================================================
!This subroutine assigns a counting value to every
!point inside the drum, not including the boundary.
!========================================================
        integer :: IGi, IGj
        !iterators
        integer :: IGcounter
        !to number the indices
        complex :: IGij        

        call makeGrid()
        allocate(Indices(2*Nreq+2,2*Nreq+2))
        
        Indices = 0
        IGcounter = 1
        do IGi = 1, 2*Nreq+2
                do IGj = 1, 2*Nreq+2
                        IGij = CMPLX(IGi,IGj)
                        !Checking if we are inside the drum and
                        !not on the boundary
                        if (.not. ANY(ExtendedArray .eq. IGij) .and.&
                           (Grid(IGi,IGj) .eq. 1)) then
                                Indices(IGi,IGj) = IGcounter
                                IGcounter = IGcounter + 1 
                        else if (Grid(IGi,IGj) .eq. 1) then
                                Indices(IGi,IGj) = -1
                        end if
                end do
        end do
        Ncounter = IGcounter -1
END SUBROUTINE

SUBROUTINE generateAmatrix()

!=======================================================
!This subroutine generates the matrix A
!for the Dirichlet eigenvalue problem.
!=======================================================
        integer :: GAi, GAj, GAu, GAd, GAl, GAr, GAp
        !iterators and neighbourcheckers
        
        allocate(AMatrix(Ncounter,Ncounter))
        allocate(Ucmplx(Ncounter))
!        allocate(U(Ncounter))
 
        AMatrix = 0
        do GAi = 2,2*Nreq+2
                do GAj = 2,2*Nreq+2
                        if (Indices(GAi,GAj) <= 0) then
                                cycle
                        end if
                        GAp = Indices(GAi,GAj)
                        GAu = Indices(GAi,GAj+1)
                        GAd = Indices(GAi,GAj-1)
                        GAl = Indices(GAi-1,GAj)
                        GAr = Indices(GAi+1,GAj)
                      AMatrix(GAp,GAp) = 4
                        Ucmplx(GAp) = CMPLX(GAi,GAj)
                        if (GAu .gt. 0) then
                                AMatrix(GAu,GAp) = -1
                        end if
                        if (GAd .gt. 0) then
                                AMatrix(GAd,GAp) = -1
                        end if
                        if (GAl .gt. 0) then
                                AMatrix(GAl,GAp) = -1
                        end if
                        if (GAr .gt. 0) then
                                AMatrix(GAr,GAp) = -1
                        end if
              end do
        end do
END SUBROUTINE

SUBROUTINE solveEVP()

!======================================================
!This routine solves the eigenvalueproblem using
!the LAPACK routine LA_SYEV, and maps the eigenvector
!values back to the corresponding points on the grid.
!======================================================
        integer :: SEVPi
        !iterators
        call LA_SYEVD(Amatrix,U,'V','L')
        allocate(modeGrid(Ncounter,Ncounter))
        modegrid = 0
        do SEVPi = 1, Ncounter
                modeGrid(nint(real(Ucmplx(SEVPi))),nint(IMAG(Ucmplx(SEVPi)))) = &
                AMatrix(SEVPi,modenumber)                
        end do
END SUBROUTINE

SUBROUTINE plotMode()

!====================================================
!This subroutine writes the data to file
!====================================================
        integer :: PMi, PMj
        !Iterators
        real(wp) :: PMx, PMy
        
        PMx = delta
        PMy = PMx
        call openfile("mode.dat",mode)
        do PMi = 1, 2*Nreq + 2
                do PMj = 1,2*Nreq + 2
                        write(mode,*) PMx, PMy, modeGrid(PMi,PMj)
                        PMy = PMy + delta 
                end do
                write(mode,*)
                PMy = delta
                PMx = PMx + delta
        end do 
        
END SUBROUTINE

SUBROUTINE plot_fractal()

!=======================================================
!This subroutine uses gnuplot to plot the 2-dim fractal
!and saves it as a .png file
!=======================================================
        integer, parameter :: gnuplotter = 28
        open (unit=gnuplotter, file = "plotModes.gnu")
        write(gnuplotter,*) 'set terminal png size 600,500 enhanced font "Helvetica,12"'
        write(gnuplotter,*) 'set output "figures/fractal',l_dim-1,'.png"'
        write(gnuplotter,*) 'set xlabel "x-points (x/L)"'
        write(gnuplotter,*) 'set title "The plot of mode number ',l_dim-1,'"' 
        write(gnuplotter,*) 'set ylabel "y-points (y/L)"'
        write(gnuplotter,*) 'plot "2dfractal.dat" w l notitle' 
        Call SYSTEM('gnuplot -p "plotModes.gnu"') 
        Call SYSTEM('rm plotModes.gnu')
END SUBROUTINE



SUBROUTINE plot_mode()

!=======================================================
!This subroutine uses gnuplot to plot the eigenmodes, and
!saves it to a .png file
!=======================================================
        integer, parameter :: gnuplotter = 28
        open (unit=gnuplotter, file = "plotModes.gnu")
        write(gnuplotter,*) 'set terminal png size 600,500 enhanced font "Helvetica,12"'
        write(gnuplotter,*) 'set output "figures/mode',modeNumber,'.png"'
        write(gnuplotter,*) 'set xlabel "x-points (x/L)"'
        write(gnuplotter,*) 'set title "The plot of mode number ',modeNumber,'"' 
        write(gnuplotter,*) 'set ylabel "y-points (y/L)"'
        write(gnuplotter,*) 'splot "mode.dat" w l notitle,&
                             "3dfractal.dat" w l fc rgb "black" notitle' 
        Call SYSTEM('gnuplot -p "plotModes.gnu"') 
        Call SYSTEM('rm plotModes.gnu')
END SUBROUTINE

SUBROUTINE findCorners()

!======================================================
!This subroutine examines the neighbouring points of
!the point of interest, and checks if its at a corner 
!or on a line part. 
!======================================================
        integer :: FCi, FCj
        !Iterators
        allocate(notCorners(Ncounter))
        notCorners = 0
        do FCi = 1,2*Nreq+2
                do FCj = 1, 2*Nreq+2
                        if (Indices(FCi,FCj) .gt. 0) then
                                if (Indices(FCi+2,FCj).eq.0) then
                                        notCorners(Indices(FCi,FCj)) = &
                                        notCorners(INdices(FCi,FCJ)) + 1
                                end if
                                if (Indices(FCi-2,FCj).eq.0) then 
                                        notCorners(Indices(FCi,FCj)) = &
                                        notCorners(INdices(FCi,FCJ)) + 1
                                end if                                
                                if (Indices(FCi,FCj+2).eq.0) then
                                        notCorners(Indices(FCi,FCj)) = &
                                        notCorners(INdices(FCi,FCJ)) + 1
                                end if                                
                                if (INdices(FCi,FCj-2).eq.0) then
                                        notCorners(Indices(FCi,FCj)) = &
                                        notCorners(INdices(FCi,FCJ)) + 1
                                end if                                
                        end if
                end do
        end do
END SUBROUTINE

SUBROUTINE generateBmatrix()

!=======================================================
!This subroutine generates the matrix for the biharmonic
!eigenvalue problem with homogeneous Dirichlet boundary
!conditions.
!=======================================================
         integer :: GBi, GBj, GBu, GBd, GBl, GBr, GBp
        !iterators and neighbourcheckers
        integer :: GBu2, GBd2, GBl2, GBr2
        integer :: GBur, GBul, GBdr, GBdl
        
        allocate(Bmatrix(Ncounter,Ncounter))
        allocate(UBcmplx(Ncounter))
        call findCorners()
        BMatrix = 0
        do GBi = 2,2*Nreq+2
                do GBj = 2,2*Nreq+2
                        if (Indices(GBi,GBj) <= 0) then
                                cycle
                        end if
                        GBp = Indices(GBi,GBj)
                        GBu = Indices(GBi,GBj+1)
                        GBd = Indices(GBi,GBj-1)
                        GBl = Indices(GBi-1,GBj)
                        GBr = Indices(GBi+1,GBj)
                        GBu2 = Indices(GBi,GBj+2)
                        GBd2 = Indices(GBi,GBj-2)
                        GBl2 = Indices(GBi-2,GBj)
                        GBr2 = Indices(GBi+2,GBj)
                        GBur = Indices(GBi+1,GBj+1)
                        GBul = Indices(GBi-1,GBj+1)
                        GBdl = Indices(GBi-1,GBj-1)
                        GBdr = Indices(GBi+1,GBj-1)
                        BMatrix(GBp,GBp) = 20 + notCorners(GBp)
                        UBcmplx(GBp) = CMPLX(GBi,GBj)
                        if (GBu .gt. 0) then
                                BMatrix(GBu,GBp) = -8
                        end if
                        if (GBd .gt. 0) then
                                BMatrix(GBd,GBp) = -8
                        end if
                        if (GBl .gt. 0) then
                                BMatrix(GBl,GBp) = -8
                        end if
                        if (GBr .gt. 0) then
                                BMatrix(GBr,GBp) = -8
                        end if
                         if (GBu2 .gt. 0) then
                                BMatrix(GBu2,GBp) = 1
                        end if
                        if (GBd2 .gt. 0) then
                                BMatrix(GBd2,GBp) = 1
                        end if
                        if (GBl2 .gt. 0) then
                                BMatrix(GBl2,GBp) = 1
                        end if
                        if (GBr2 .gt. 0) then
                                BMatrix(GBr2,GBp) = 1
                        end if
                         if (GBur .gt. 0) then
                                BMatrix(GBur,GBp) = 2
                        end if
                        if (GBul .gt. 0) then
                                BMatrix(GBul,GBp) = 2
                        end if
                        if (GBdl .gt. 0) then
                                BMatrix(GBdl,GBp) = 2
                        end if
                        if (GBdr .gt. 0) then
                                BMatrix(GBdr,GBp) = 2
                        end if
              end do
        end do
END SUBROUTINE

SUBROUTINE solveEVPB()

!======================================================
!This routine solves the eigenvalueproblem using
!the LAPACK routine LA_SYEV, and maps the eigenvector
!values back to the corresponding points on the grid.
!======================================================
        integer :: SEVPBi
        !iterators
        call LA_SYEVD(Bmatrix,U,'V','L')
        allocate(modeGrid(Ncounter,Ncounter))
        modegrid = 0
        do SEVPBi = 1, Ncounter
                modeGrid(nint(real(UBcmplx(SEVPBi))),nint(IMAG(UBcmplx(SEVPBi)))) = &
                BMatrix(SEVPBi,modenumber)                
        end do
END SUBROUTINE

SUBROUTINE plot_biharmonic()

!=======================================================
!This subroutine uses gnuplot to plot the biharmonic
!eigenmodes, and saves it to a .png file
!=======================================================
        integer, parameter :: gnuplotter = 28
        open (unit=gnuplotter, file = "plotModes.gnu")
        write(gnuplotter,*) 'set terminal png size 600,500 enhanced font "Helvetica,12"'
        write(gnuplotter,*) 'set output "bi/bimode',modeNumber,'.png"'
        write(gnuplotter,*) 'set xlabel "x-points (x/L)"'
        write(gnuplotter,*) 'set title "The plot of biharmonic mode number ',modeNumber,'"' 
        write(gnuplotter,*) 'set ylabel "y-points (y/L)"'
        write(gnuplotter,*) 'splot "mode.dat" w l notitle,&
                             "3dfractal.dat" w l fc rgb "black" notitle' 
        Call SYSTEM('gnuplot -p "plotModes.gnu"') 
        Call SYSTEM('rm plotModes.gnu')
END SUBROUTINE


END MODULE

