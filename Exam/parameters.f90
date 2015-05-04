MODULE parameters
IMPLICIT NONE
integer, parameter                              :: wp=kind(0.d0)        
!Set working precission

integer, parameter                              :: l_max = 1            
!Setting max iterations for generating Koch island

real(kind=wp), parameter                        :: length = 1.0_wp      
!Setting lenth of the fractal sides

integer, parameter                              :: l_dim = 3

integer                                         :: l_level

real(wp)                                        :: delta = 1._wp/(4**(l_dim-1))
!The grid step size for the square lattice

integer                                         :: Nreq
!The number of grid points in matrix enclosing fractal

integer, parameter                              :: Nfrac = 8**(l_dim-1) +1
!setting the dimension of the fractal array

complex, dimension(:), allocatable              :: FractalArray         
!The array for the points defining the fractal
complex, dimension(4*Nfrac)                     :: IslandArray

complex, dimension(8*Nfrac)                     :: ExtendedArray
!The extended Island array with greater precission

complex, dimension(9)                           :: lineSegment
!Declearing a sub-array for the line segments that
!are subject to change 

integer, dimension(:,:), allocatable              :: Grid
!the grid with the fractal

integer                                           :: GVxreq, GVxmin
!for calculating the max value of x in the grid
integer                                           :: GVyreq, GVymin
!For calculating the max value of y in the grid





!integer, parameter                      :: test = 29
END MODULE

