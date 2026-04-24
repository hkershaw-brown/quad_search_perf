program test_quad_search_acc

use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities, &
                              my_task_id
use types_mod, only : r8
use quad_utils_mod, only : quad_interp_handle, init_quad_interp, &
                           set_quad_coords, quad_lon_lat_locate, &
                           GRID_QUAD_FULLY_IRREGULAR, &
                           QUAD_LOCATED_CELL_CENTERS
use netcdf_utilities_mod, only : nc_open_file_readonly, nc_get_variable_size, &
                                 nc_get_variable, nc_close_file
use random_seq_mod, only : init_random_seq, random_uniform, random_seq_type

implicit none

integer :: Nx, Ny, N(2)
real(r8), allocatable :: TLON(:,:), TLAT(:,:)
integer :: lon_corner_index(4), lat_corner_index(4)
integer :: lstatus
type(quad_interp_handle) :: h
integer :: i, num_reps
real(r8), allocatable :: lon(:), lat(:)
character(len=256) :: grid_file
character(len=256) :: out_file
character(len=256) :: arg

type(random_seq_type) :: ran_seq
integer :: ncid, iunit

call initialize_mpi_utilities('test_quad_search_acc')

!------------------------------------------------------------
! Parse command-line arguments:  [num_reps] [grid_file] [out_file]
!   default num_reps  = 1000
!   default grid_file = 'irregular_grid.nc'
!   default out_file  = 'quad_results.txt'
!------------------------------------------------------------
num_reps  = 1000
grid_file = 'irregular_grid.nc'
out_file  = 'quad_results.txt'

if (command_argument_count() >= 1) then
   call get_command_argument(1, arg)
   read(arg, *) num_reps
endif

if (command_argument_count() >= 2) then
   call get_command_argument(2, grid_file)
endif

if (command_argument_count() >= 3) then
   call get_command_argument(3, out_file)
endif

if (my_task_id() == 0) then
   print *, 'Running test_quad_search_acc with:'
   print *, 'num_reps  = ', num_reps
   print *, 'grid_file = ', trim(grid_file)
   print *, 'out_file  = ', trim(out_file)
endif

! read grid
ncid = nc_open_file_readonly(trim(grid_file))

call nc_get_variable_size(ncid, 'longitudes', N)
Nx = N(1)
Ny = N(2)

if (my_task_id() == 0) then
   print*, 'Nx = ', Nx
   print*, 'Ny = ', Ny
endif

allocate(TLON(Nx,Ny), TLAT(Nx,Ny))

call nc_get_variable(ncid, 'longitudes', TLON)
call nc_get_variable(ncid, 'latitudes',  TLAT)
call nc_close_file(ncid)

! generate num_reps random lon/lat pairs (same seed as test_quad_search)
allocate(lon(num_reps), lat(num_reps))

call init_random_seq(ran_seq, 12345)
do i = 1, num_reps
   lon(i) = random_uniform(ran_seq)*360.0_r8
   lat(i) = random_uniform(ran_seq)*180.0_r8 - 90.0_r8
end do

! set up interp handle
call init_quad_interp(GRID_QUAD_FULLY_IRREGULAR, Nx, Ny, &
                      QUAD_LOCATED_CELL_CENTERS,         &
                      global         = .true.,           &
                      spans_lon_zero = .true.,           &
                      pole_wrap      = .true.,           &
                      interp_handle  = h)
call set_quad_coords(h, TLON, TLAT)

! open output file and write a header
iunit = 42
open(unit=iunit, file=trim(out_file), form='formatted', action='write', status='replace')
write(iunit,'(A)') '# obs_lon  obs_lat  lstatus  lon_idx(1:4)  lat_idx(1:4)'

do i = 1, num_reps
   call quad_lon_lat_locate(h, lon(i), lat(i), lon_corner_index, lat_corner_index, lstatus)
   write(iunit, '(2F12.6, I4, 4I8, 4I8)') &
         lon(i), lat(i), lstatus, lon_corner_index, lat_corner_index
end do

close(iunit)

if (my_task_id() == 0) then
   print*, 'Results written to ', trim(out_file)
endif

call finalize_mpi_utilities()

end program test_quad_search_acc
