program test_quad_search

use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities, &
                              my_task_id
use types_mod, only : r8
use quad_utils_mod, only : quad_interp_handle, init_quad_interp, &
                           set_quad_coords, quad_lon_lat_locate, &
                           GRID_QUAD_FULLY_IRREGULAR, &
                           QUAD_LOCATED_CELL_CENTERS, &
                           dump_coarse_box_data
use netcdf_utilities_mod, only : nc_open_file_readonly, nc_get_variable_size, &
                                 nc_get_variable, nc_close_file
use random_seq_mod, only : init_random_seq, random_uniform, random_seq_type

use mpi

implicit none

integer :: Nx, Ny, N(2)
real(r8), allocatable :: TLON(:,:), TLAT(:,:)
integer :: lon_corner_index(4), lat_corner_index(4)
integer :: lstatus
type(quad_interp_handle) :: h
real(r8) :: start
integer :: i, num_reps
real(r8), allocatable :: lon(:), lat(:)
character(len=256) :: grid_file
character(len=256) :: dump_prefix
character(len=256) :: arg
logical :: print_results, regional, dump_boxes
logical :: is_global, spans_zero, pole

type(random_seq_type) :: ran_seq
integer :: ncid

call initialize_mpi_utilities('test_quad_search')

!------------------------------------------------------------
! Parse command-line arguments:  [num_reps] [grid_file] [flags]
!   default num_reps      = 1000
!   default grid_file     = 'irregular_grid.nc'
!   flags (any order after grid_file):
!     'p'              -- print per-observation quad results
!     'r'              -- regional grid (disables global/cyclic/pole flags)
!     'd[prefix]'      -- dump coarse box data; optional prefix (default 'coarse_boxes')
!------------------------------------------------------------
num_reps     = 1000
grid_file    = 'irregular_grid.nc'
print_results = .false.
regional      = .false.
dump_boxes    = .false.
dump_prefix   = 'coarse_boxes'

if (command_argument_count() >= 1) then
   call get_command_argument(1, arg)
   read(arg, *) num_reps
endif

if (command_argument_count() >= 2) then
   call get_command_argument(2, grid_file)
endif

do i = 3, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg) == 'p')  print_results = .true.
   if (trim(arg) == 'r')  regional      = .true.
   if (arg(1:1)   == 'd') then
      dump_boxes = .true.
      if (len_trim(arg) > 1) dump_prefix = arg(2:)
   endif
end do

if (my_task_id() == 0) then
   print *, 'Running test_quad_search with:'
   print *, 'num_reps  = ', num_reps
   print *, 'grid_file = ', trim(grid_file)
   print *, 'regional  = ', regional
   print *, 'dump_boxes= ', dump_boxes
   if (dump_boxes) print *, 'dump_prefix= ', trim(dump_prefix)
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
call nc_get_variable(ncid, 'latitudes', TLAT)

! set global/cyclic/pole flags from regional flag
is_global  = .not. regional
spans_zero = .not. regional
pole       = .not. regional

! set up interp handle
call init_quad_interp(GRID_QUAD_FULLY_IRREGULAR, Nx, Ny, &
                      QUAD_LOCATED_CELL_CENTERS,         &
                      global         = is_global,        &
                      spans_lon_zero = spans_zero,       &
                      pole_wrap      = pole,             &
                      interp_handle  = h)
call set_quad_coords(h, TLON, TLAT)

! optionally write coarse box data for plotting
if (dump_boxes .and. my_task_id() == 0) then
   call dump_coarse_box_data(h, trim(dump_prefix))
   print *, 'Wrote coarse box data with prefix: ', trim(dump_prefix)
endif

! generate num_reps random lon/lat pairs within the grid domain for regional,
! or over the full globe for global
allocate(lon(num_reps), lat(num_reps))

call init_random_seq(ran_seq, 12345)
if (regional) then
   do i = 1, num_reps
      lon(i) = random_uniform(ran_seq) * (maxval(TLON) - minval(TLON)) + minval(TLON)
      lat(i) = random_uniform(ran_seq) * (maxval(TLAT) - minval(TLAT)) + minval(TLAT)
   end do
else
   do i = 1, num_reps
      lon(i) = random_uniform(ran_seq)*360.0_r8
      lat(i) = random_uniform(ran_seq)*180.0_r8 - 90.0_r8
   end do
endif

start = mpi_wtime()
do i = 1, num_reps
   call quad_lon_lat_locate(h, lon(i), lat(i), lon_corner_index, lat_corner_index, lstatus)
end do
if (my_task_id() == 0) then
   print*, 'locate time = ', mpi_wtime() - start
endif

if (print_results .and. my_task_id() == 0) then
   print*, 'lon corner index = ', lon_corner_index
   print*, 'lat corner index = ', lat_corner_index
endif

call finalize_mpi_utilities()

end program test_quad_search
