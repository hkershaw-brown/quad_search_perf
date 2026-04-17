program test_quad_search

use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities
use types_mod, only : r8
use quad_utils_mod, only : quad_interp_handle, init_quad_interp, &
                           set_quad_coords, quad_lon_lat_locate, &
                           quad_lon_lat_evaluate, quad_interp_handle, &
                           GRID_QUAD_FULLY_IRREGULAR, &
                           QUAD_LOCATED_CELL_CENTERS

implicit none

real(r8) :: lon, lat
integer :: Nx, Ny
real(r8), allocatable :: TLON(:,:), TLAT(:,:)
logical, allocatable :: TMSK(:,:)
integer :: lon_corner_index(4), lat_corner_index(4)
integer :: lstatus
type(quad_interp_handle) :: h



call allocate(TLON(Nx,Ny), TLAT(Nx,Ny), TMSK(Nx,Ny))

call initialize_mpi_utilities('test_quad_search')

! set up grid
Nx = 10
Ny = 12

call allocate(TLON(Nx,Ny), TLAT(Nx,Ny), TMSK(Nx,Ny))


! @todo HK not global
! set up interp handle
call init_quad_interp(GRID_QUAD_FULLY_IRREGULAR, Nx, Ny, &
                      QUAD_LOCATED_CELL_CENTERS,         &
                      global         = .false.,          & 
                      spans_lon_zero = .true.,           &
                      pole_wrap      = .true.,           & 
                      interp_handle  = h)
call set_quad_coords(h, TLON, TLAT, TMSK)

call quad_lon_lat_locate(h, lon, lat, lon_corner_index, lat_corner_index, lstatus) 

! check it is the correct quad. 

call finalize_mpi_utilities()

end program test_quad_search
