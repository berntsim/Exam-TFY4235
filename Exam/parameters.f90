MODULE parameters
IMPLICIT NONE

!============= General parameters section =============================
integer, parameter                              :: wp=kind(0.d0)        
!Set working precission

integer, parameter                              :: l_max = 1            
!Setting max iterations for generating Koch island

real(kind=wp), parameter                        :: length = 1.0_wp      
!Setting lenth of the fractal sides

integer, parameter                              :: l_dim = 2 
!Determine what generation of fractal to generate

integer, parameter                              :: modenumber = 10 
!To determine what mode to plot

!=============Fractal generator section==================================
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

!============= Putting the fractal on the grid ==========================


integer, dimension(:,:), allocatable              :: Grid
!the grid with the fractal

integer                                           :: GVxreq, GVxmin
!for calculating the max value of x in the grid
integer                                           :: GVyreq, GVymin
!For calculating the max value of y in the grid

!============ Solving the eigenvalue problem ============================

integer, dimension (:,:), allocatable           :: Indices
!To keep track of indecies on the drum

!real(wp), dimension(:), allocatable             :: U
!The U to 

complex, dimension(:), allocatable              :: Ucmplx
!Using complex for handy mapping of grid points indices
!to actual position on the grid.

real(wp), dimension(:,:), allocatable           :: Amatrix
!The matrix for solving the dirichlet eigenvalue problem

integer                                         :: Ncounter
!Dimension for some of the matrices used in solving 
!the eigenmodes.

real(wp), dimension(:), allocatable             :: U
!The eigenvector will be stored in this value 

real(wp), dimension(:,:), allocatable           :: modeGrid
!The grid with the mode vibrations

!===============Biharmonic section=================================
integer, dimension(:), allocatable              :: notCorners
!To keep track of neighbours in the nabla cubed operator

real(wp), dimension(:,:), allocatable           :: Bmatrix
!The matrix for the biharmonic eigenvalue problem

complex, dimension(:), allocatable              :: UBcmplx
END MODULE

